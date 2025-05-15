from __future__ import annotations
import fnmatch
import io
import json
import logging
import pathlib
from typing import Literal

import requests
from tqdm import tqdm

from tryptag.cache import Cache, FileNotCachedError, FileTypes
from tryptag.datasource import DataSource

logger = logging.getLogger("tryptag.bia")

BIA_API_ROOT = "https://ftp.ebi.ac.uk/pub/databases/biostudies/S-BIAD"
DOWNLOAD_MAX_TRIES = 3

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
        outfile: io.BufferedIOBase
    ):
        """
        Download the file from BIA and write to the given file object.

        :param outfile: file-like, file object the downloaded file should be
            written to.
        """

        logger.debug(f"Downloading {self.path} from URL {self.url}.")
        r = requests.get(
            self.url,
            stream=True,
            headers={
                "Accept-Encoding": "gzip",
            }
        )
        r.raise_for_status()
        total_size = int(r.headers.get("content-length", 0))
        block_size = 1024

        def log_progress(chunks):
            if total_size == 0:
                for chunk in chunks:
                    yield chunk
            else:
                progress = 0
                delta = 0
                for chunk in chunks:
                    delta += block_size
                    progress += block_size
                    delta_percentage = 100 * (delta / total_size)
                    if delta_percentage >= 10:
                        logger.debug(
                            f"{self.path}: {progress} / {total_size}")
                        delta = 0
                    yield chunk
                logger.debug(f"{self.path}: {total_size} / {total_size}")

        # TODO: check downloaded file size!

        with tqdm(
            desc=f"Downloading {self.path}",
            total=total_size,
            unit="B",
            unit_scale=True,
            disable=logger.getEffectiveLevel() != logging.INFO,
        ) as progress_bar:
            for chunk in log_progress(r.iter_content(block_size)):
                progress_bar.update(len(chunk))
                outfile.write(chunk)


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
        outfile = filepath.open("wb")
        try:
            with outfile:
                biafile.download(outfile)
        except Exception:
            filepath.unlink(missing_ok=True)
            raise

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
