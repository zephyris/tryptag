from __future__ import annotations

import fnmatch
import io
import json
import logging
import os
import pathlib
import random
import tempfile
import time
from typing import Literal

import requests
from tqdm import tqdm

from tryptag.cache import Cache, FileNotCachedError, FileTypes
from tryptag.datasource import DataSource

logger = logging.getLogger("tryptag.bia")

BIA_API_ROOT = "https://ftp.ebi.ac.uk/biostudies/fire/S-BIAD/"
DOWNLOAD_MAX_TRIES = 5
DOWNLOAD_TIMEOUT: tuple[float, float] = (30.0, 60.0)
DOWNLOAD_BLOCK_SIZE = 65536

FILE_TYPES = {
    "Cell ROIs": FileTypes.CELL_ROIS,
    "Corrected image": FileTypes.IMAGE,
    "Thresholded image": FileTypes.THRESHOLDED,
    "Localisation annotation": FileTypes.LOCALISATION,
    "Other": FileTypes.OTHER,
}


class BIANotAvailableError(Exception):
    def __init__(self):
        super().__init__("Connection problem to BioImage Archive, please try "
                         "again later.")


class IncompleteDownloadError(OSError):
    """Raised when a streamed download terminates before all expected bytes
    are received."""
    pass


def _is_permanent_http_error(exc: Exception) -> bool:
    """Return True for non-retryable 4xx HTTP errors (excluding transient
    403/408/429 responses)."""
    return (
        isinstance(exc, requests.exceptions.HTTPError)
        and exc.response is not None
        and 400 <= exc.response.status_code < 500
        and exc.response.status_code not in (403, 408, 429)
    )


class BiaFile:
    path: str
    size: int
    file_type: FileTypes = FileTypes.OTHER
    plate: str | None = None
    gene_id: str | None = None
    terminus: Literal["C", "N"] | None = None
    field_index: int | None = None

    def __init__(self, data: dict, bia: BioimageArchive):
        self.path = data["path"]
        self.url = bia._file_url(self.path)
        self.size = data["size"]
        if "attributes" in data:
            attributes = {
                e["name"]: e["value"] for e in data["attributes"]
                if "value" in e
            }
            if "FileType" in attributes:
                self.file_type = FILE_TYPES[attributes["FileType"]]
            if "Plate" in attributes:
                self.plate = f"{attributes['Plate']}_{attributes['Date']}"
            if "Gene ID" in attributes:
                self.gene_id = attributes["Gene ID"]
            if "Terminus" in attributes:
                self.terminus = attributes["Terminus"]
            if "Field" in attributes:
                self.field_index = int(attributes["Field"])

    def download(
        self,
        outfile: io.BufferedIOBase,
    ):
        """
        Download the file from BIA and write to the given file object.

        :param outfile: file-like, file object the downloaded file should be
            written to.
        """

        logger.debug(f"Downloading {self.path} from URL {self.url}.")
        with requests.get(
            self.url,
            stream=True,
            headers={
                "Accept-Encoding": "gzip",
            },
            timeout=DOWNLOAD_TIMEOUT,
        ) as r:
            r.raise_for_status()
            content_length = int(r.headers.get("content-length", 0))
            content_encoding = r.headers.get("content-encoding", "").lower()
            is_compressed = content_encoding not in ("", "identity")
            if self.size > 0:
                total_size = self.size
            elif not is_compressed:
                total_size = content_length
            else:
                total_size = 0

            def log_progress(chunks):
                if total_size == 0:
                    for chunk in chunks:
                        yield chunk
                else:
                    progress = 0
                    delta = 0
                    for chunk in chunks:
                        chunk_len = len(chunk)
                        delta += chunk_len
                        progress += chunk_len
                        delta_percentage = 100 * (delta / total_size)
                        if delta_percentage >= 10:
                            logger.debug(
                                f"{self.path}: {progress} / {total_size}")
                            delta = 0
                        yield chunk
                    logger.debug(f"{self.path}: {progress} / {total_size}")

            bytes_written = 0
            with tqdm(
                desc=f"Downloading {self.path}",
                total=total_size,
                unit="B",
                unit_scale=True,
                disable=logger.getEffectiveLevel() != logging.INFO,
            ) as progress_bar:
                for chunk in log_progress(r.iter_content(DOWNLOAD_BLOCK_SIZE)):
                    if not chunk:
                        continue
                    chunk_len = len(chunk)
                    bytes_written += chunk_len
                    progress_bar.update(chunk_len)
                    outfile.write(chunk)

        if self.size > 0 and bytes_written != self.size:
            raise IncompleteDownloadError(
                f"Incomplete download for {self.path}: received "
                f"{bytes_written} of {self.size} expected bytes."
            )
        if (
            self.size <= 0
            and not is_compressed
            and content_length > 0
            and bytes_written != content_length
        ):
            raise IncompleteDownloadError(
                f"Incomplete download for {self.path}: received "
                f"{bytes_written} of {content_length} expected bytes."
            )


class BioimageArchive(DataSource):
    """Class to read TrypTag data from Bioimage Archive

    This subclasses the `DataSource` abstract base class and implements the
    methods to access data via the EMBL Bioimage Archive.
    """
    def __init__(
        self,
        cache: Cache,
        accession: str = "S-BIAD1866",
    ):
        super().__init__(cache)
        self.accession = accession
        try:
            self._get_file_index()
        except FileNotCachedError:
            raise BIANotAvailableError()
        self.__post_init__()

    @property
    def _url_root(self):
        return (
            f"{BIA_API_ROOT}/{self.accession.replace('S-BIAD1', '')}"
            f"/{self.accession}/Files/"
        )

    def _file_url(self, path: str):
        return f"{self._url_root}/{path}"

    def fetch_root_file(self, filename: str):
        logger.debug(f"Fetching root file {filename}.")
        biafile = self._file_index[filename]
        filepath = self.cache.file_path(biafile.path)
        for attempt in range(1, DOWNLOAD_MAX_TRIES + 1):
            tmp_path: pathlib.Path | None = None
            try:
                with tempfile.NamedTemporaryFile(
                    mode="wb",
                    dir=filepath.parent,
                    prefix=f".{filepath.name}.",
                    suffix=".tmp",
                    delete=False,
                ) as tmpfile:
                    tmp_path = pathlib.Path(tmpfile.name)
                    biafile.download(tmpfile)
                    tmpfile.flush()
                os.replace(tmp_path, filepath)
                tmp_path = None
                return
            except (requests.exceptions.RequestException, OSError) as e:
                if _is_permanent_http_error(e):
                    logger.warning(
                        f"Non-retryable error downloading {filename}: "
                        f"{e.__class__.__name__} - {e}."
                    )
                    raise
                if attempt >= DOWNLOAD_MAX_TRIES:
                    logger.warning(
                        f"Problem downloading {filename}, all "
                        f"{DOWNLOAD_MAX_TRIES} attempts failed, giving up."
                    )
                    raise
                delay = float(2 ** (attempt - 1)) * random.uniform(0.5, 1.5)
                logger.warning(
                    f"Problem downloading {filename}, attempt {attempt} of "
                    f"{DOWNLOAD_MAX_TRIES} failed ({e.__class__.__name__} - "
                    f"{e}). Retrying in {delay:.1f}s."
                )
                time.sleep(delay)
            finally:
                if tmp_path is not None:
                    try:
                        tmp_path.unlink(missing_ok=True)
                    except OSError:
                        pass

    def fetch_plate_file(self, plate, filename):
        self.fetch_root_file(f"{plate}/{filename}")

    def _get_file_index(self):
        logger.debug("Loading file index.")

        # Slightly hackish, but gets the job done without duplicating
        # download code.
        self._file_index: dict[str, BiaFile] = {
            "Files.json": BiaFile({"path": "Files.json", "size": -1}, self)
        }
        with self.load_root_file("Files.json") as filesfile:
            data: list[dict] = json.load(filesfile)

        self._file_index = {}
        for entry in data:
            biafile = BiaFile(entry, self)
            if biafile.path in self._file_index:
                raise ValueError(f"Duplicate file in index {biafile.path}")
            self._file_index[biafile.path] = biafile

    def glob_plate_files(self, plate: str, pattern: str):
        """
        Find and return a list of files matching the given pattern in the
        given plate.

        :param plate: str, the name of the plate the files are located in
        :param pattern: str, a shell-like glob pattern
        :return: a list of file matching the pattern
        """
        matches = fnmatch.filter(
            self._file_index,
            f"{plate}/{pattern}",
        )
        return [
            str(pathlib.Path(p).relative_to(plate)) for p in matches
        ]
