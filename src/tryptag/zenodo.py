import hashlib
import io
import json
import logging
import typing
import requests
from tqdm import tqdm

from tryptag.cache import Cache
from tryptag.datasource import DataSource

ZENODO_API_ROOT = "https://zenodo.org/api/records"
DOWNLOAD_MAX_TRIES = 3

logger = logging.getLogger("tryptag.datasource.zenodo")


class ZenodoFileChecksumError(Exception):
    pass


class ZenodoFile:
    """Class describing a file in a Zenodo record."""
    def __init__(self, data: dict):
        """
        Initialise a Zenodo record file object.

        :param data: dict, the information about the file in the Zenodo record
        """
        self.id: str = data["id"]
        self.name: str = data["key"]
        self.url: str = data["links"]["self"]
        self.checksum: str = data["checksum"].split(":")[-1]

    def fetch_lines(self):
        """
        Fetch the file from Zenodo via HTTP and return a line iterator
        """
        r = requests.get(self.url, stream=True)

        md5 = hashlib.md5()

        for line in r.iter_lines():
            md5.update(line + b"\n")
            yield line

        if md5.hexdigest() != self.checksum:
            raise ZenodoFileChecksumError

    def download(self, outfile: io.BufferedIOBase):
        """
        Download the file and write to the given file object.

        :param outfile: file-like, file object the downloaded file should be
            written to.
        """
        logger.debug(f"Downloading {self.name} from URL {self.url}.")
        for try_nr in range(DOWNLOAD_MAX_TRIES):
            r = requests.get(self.url, stream=True)
            total_size = int(r.headers.get("content-length", 0))
            block_size = 1024

            md5 = hashlib.md5()

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
                                f"{self.name}: {progress} / {total_size}")
                            delta = 0
                        yield chunk
                    logger.debug(f"{self.name}: {total_size} / {total_size}")

            with tqdm(
                desc=f"Downloading {self.name}",
                total=total_size,
                unit="B",
                unit_scale=True,
                disable=logger.getEffectiveLevel() != logging.INFO,
            ) as progress_bar:
                progress_bar.n
                for chunk in log_progress(r.iter_content(block_size)):
                    progress_bar.update(len(chunk))
                    md5.update(chunk)
                    outfile.write(chunk)

            if self.checksum == md5.hexdigest():
                break
            logger.warning(f"Downloading file {self.name} from {self.url} "
                           "failed {try_nr+1} times.")
        else:
            raise ValueError(f"Downloading file {self.name} from {self.url} "
                             "failed and we reached the max number of "
                             "retries.")


class ZenodoRecord:
    """Class holding a Zenodo record."""
    def __init__(self, data: dict):
        """
        Initialise the object from the given data.

        :param data: dict, the Zenodo record information
        """
        self.doi = data["doi"]
        file_list = [ZenodoFile(f) for f in data["files"]]
        self.files = {f.name: f for f in file_list}
        self._original_data = data

    @staticmethod
    def fetch(record_id: int):
        """
        Fetch the record with the given ID from Zenodo.

        :param record_id: int, ID of the record
        :return: the initialised `ZenodoRecord` object.
        """
        logger.debug("Fetching record {record_id}.")
        r = requests.get(
            f"{ZENODO_API_ROOT}/{record_id}/versions/latest"
        )
        return ZenodoRecord(r.json())

    @staticmethod
    def from_file(json_file: typing.IO):
        """
        Load the record from the given file object.

        :param json_file: file-like, file object holding the JSON information
        :return: the initialised `ZenodoRecord` object.
        """
        logger.debug(f"Opening record from file {json_file.name}.")
        with json_file:
            return ZenodoRecord(json.load(json_file))

    def to_file(self, outfile: typing.IO):
        """
        Write the `ZenodoRecord` to a JSON file.

        :param outfile: file-like, file object the record should be written to
        """
        logger.debug(f"Saving record to file {outfile.name}")
        with outfile:
            json.dump(self._original_data, outfile)


class Zenodo(DataSource):
    """Class to read TrypTag data from Zenodo.

    This subclasses the `DataSource` abstract base class and implements the
    methods to access data via the Zenodo API.
    """
    def __init__(
        self,
        cache: Cache,
        master_record_id: int = 6862289,
    ):
        """
        Initialise access to Zenodo via the given master record ID.

        :param cache: `Cache`, the file system cache object to use
        :param master_record_id: int, the ID of the project master Zenodo
            record. The default is 6862289, the procyclic form data set.
        """
        DataSource.__init__(self, cache)
        self.master_record = self._load_or_fetch_record(master_record_id)
        self.plate_index = self._get_plate_index()
        self.__post_init__()

    def _load_or_fetch_record(self, record_id: int):
        """
        Loads the data for the record with the given ID. Fetches from Zenodo
        if not found in the local cache.

        :param record_id: int, Zenodo ID of the record to be loaded
        :return: `ZenodoRecord` object
        """
        logger.debug(f"Loading record {record_id}.")
        with self.cache.lock_file_name(f"zenodo_record_{record_id}"):
            def genfile_cb():
                path = self.cache.file_path(str(record_id), "_zenodo")
                outfile = path.open("w")

                record = ZenodoRecord.fetch(record_id)
                record.to_file(outfile)

            return ZenodoRecord.from_file(
                self.cache.load_file(
                    str(record_id),
                    "_zenodo",
                    genfile_cb,
                )
            )

    def fetch_root_file(self, filename: str):
        """
        Low-level fetching of a root (non-plate) file.

        Fetches the file with the given file name from the data source and
        stores it at the path given by `self.cache.file_path(filename)`.

        :param filename: str, name of the file to be fetched
        """
        logger.debug(f"Fetching root file {filename}.")
        zfile = self.master_record.files[filename]
        with open(self.cache.file_path(filename), "wb") as outfile:
            zfile.download(outfile)

    def fetch_plate_file(self, plate: str, file: str):
        """
        Low-level fetching of a plate file.

        Fetches the file with the given file name from the data source and
        stores it at the path given by `self.cache.file_path(filename, plate)`.

        :param plate: str, name of the plate the file is included in
        :param filename: str, name of the file to be fetched
        """
        logger.debug(f"Fetching file {file} from plate {plate}.")
        self._fetch_plate(plate)

    def _fetch_plate(self, plate: str):
        """
        Fetch a whole plate's worth of data. This looks up the relevant Zenodo
        record and then downloads the ZIP file from this record.

        We cannot download single files directly because of the way the data is
        organised into ZIP files on Zenodo.

        :param plate: str, name of the plate
        """
        logger.debug(f"Fetching plate {plate}.")
        record_id = self.plate_index[plate]
        record = self._load_or_fetch_record(record_id)

        with self.cache.lock_file_name(f"zenodo_platezip_{plate}"):
            # Check if the plate already exists
            # (it might have been downloaded in a different thread / process)
            if not self.cache.file_path(plate).exists():
                tmpfile = self.cache.temporary_file(suffix=".zip", delete=True)
                with tmpfile as outfile:
                    record.files[f"{plate}_processed.zip"].download(outfile)
                    outfile.flush()
                    self.cache.extract_plate_zip(plate, outfile.name)

    def _get_plate_index(self):
        logger.debug("Loading plate index.")
        plates = {}
        with self.load_root_file("plate_doi_index.tsv") as plate_index:
            for line in plate_index:
                doi, plate_id = line.split()
                plates[plate_id] = int(doi.split(".")[-1])
        return plates

    def glob_plate_files(self, plate: str, pattern: str):
        """
        Find and return a list of files matching the given pattern in the
        given plate.

        :param plate: str, the name of the plate the files are located in
        :param pattern: str, a shell-like glob pattern
        :return: a list of file matching the pattern
        """
        if not self.cache.file_path(plate).is_dir():
            self._fetch_plate(plate)
        return [
            str(p.relative_to(plate)) for p in
            self.cache.glob_files(f"{plate}/{pattern}")
        ]

    def __hash__(self):
        return hash(("datasource", "zenodo", self.cache))
