import hashlib
import io
import json
import logging
import requests
from tqdm import tqdm

from tryptag.cache import Cache
from tryptag.datasource import DataSource

ZENODO_API_ROOT = "https://zenodo.org/api/records"
DOWNLOAD_MAX_TRIES = 3

logger = logging.getLogger(__name__)


class ZenodoFile:
    def __init__(self, data: dict):
        self.checksum: str = data["checksum"]
        self.id: str = data["id"]
        self.name: str = data["key"]
        self.url: str = data["links"]["self"]
        self.checksum: str = data["checksum"].split(":")[-1]

    def fetch_lines(self):
        r = requests.get(self.url, stream=True)

        md5 = hashlib.md5()

        for line in r.iter_lines():
            md5.update(line + b"\n")
            yield line

        print(self.checksum)
        print(md5.hexdigest())

    def download(self, outfile: io.BufferedIOBase):
        for try_nr in range(DOWNLOAD_MAX_TRIES):
            r = requests.get(self.url, stream=True)
            total_size = int(r.headers.get("content-length", 0))
            block_size = 1024

            md5 = hashlib.md5()

            with tqdm(
                desc=f"Downloading {self.name}",
                total=total_size,
                unit="B",
                unit_scale=True,
                disable=logger.getEffectiveLevel() > logging.INFO,
            ) as progress_bar:
                for chunk in r.iter_content(block_size):
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
    def __init__(self, data: dict):
        self.doi = data["doi"]
        file_list = [ZenodoFile(f) for f in data["files"]]
        self.files = {f.name: f for f in file_list}
        self._original_data = data

    @staticmethod
    def fetch(record_id: int):
        r = requests.get(
            f"{ZENODO_API_ROOT}/{record_id}/versions/latest"
        )
        return ZenodoRecord(r.json())

    @staticmethod
    def from_file(json_file: io.BufferedIOBase):
        with json_file:
            return ZenodoRecord(json.load(json_file))

    def to_file(self, outfile: io.BufferedIOBase):
        with outfile:
            json.dump(self._original_data, outfile)


class Zenodo(DataSource):
    def __init__(
        self,
        cache: Cache,
        master_record_id: int = 6862289,
    ):
        DataSource.__init__(self, cache)
        self.master_record = self._load_or_fetch_record(master_record_id)
        self.plate_index = self._get_plate_index()
        self.__post_init__()

    def _load_or_fetch_record(self, record_id: int):
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
        zfile = self.master_record.files[filename]
        with open(self.cache.file_path(filename), "wb") as outfile:
            zfile.download(outfile)

    def fetch_plate_file(self, plate, _):
        self._fetch_plate(plate)

    def _fetch_plate(self, plate: str):
        logger.debug("Fetching plate ")
        record_id = self.plate_index[plate]
        record = self._load_or_fetch_record(record_id)

        zip_lock = self.cache.lock_file_name(f"zenodo_platezip_{plate}")
        if not zip_lock.is_locked():
            # The file isn't being downloaded in another thread
            # TODO: There is a potential race condition here between the
            # is_locked check above and acquiring the lock below!
            with zip_lock:
                with self.cache.temporary_file(suffix=".zip", delete=True) as outfile:
                    record.files[f"{plate}_processed.zip"].download(outfile)
                    outfile.flush()
                    # TODO: Check MD5
                    self.cache.extract_plate_zip(plate, outfile.name)
        else:
            with zip_lock:
                pass

    def _get_plate_index(self):
        logger.debug("Loading plate index.")
        plates = {}
        with self.load_root_file("plate_doi_index.tsv") as plate_index:
            for line in plate_index:
                doi, plate_id = line.split()
                plates[plate_id] = int(doi.split(".")[-1])
        return plates

    def glob_plate_files(self, plate: str, pattern: str):
        if not self.cache.file_path(plate).is_dir():
            self._fetch_plate(plate)
        return [
            str(p.relative_to(plate)) for p in
            self.cache.glob_files(f"{plate}/{pattern}")
        ]

    def __hash__(self):
        return hash(("datasource", "zenodo", self.cache))
