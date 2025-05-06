import logging
import pathlib
import shutil
import sys
import tempfile
from typing import Callable
import zipfile

from filelock import FileLock
from tqdm import tqdm

if sys.version_info >= (3, 11):
    from enum import StrEnum
else:
    from backports.strenum import StrEnum

logger = logging.getLogger("tryptag.cache")


class FileTypes(StrEnum):
    CELL_ROIS = "_roisCells.txt"
    THRESHOLDED = "_thr.tif"
    IMAGE = ".tif"


FILE_PATTERN = FileTypes.CELL_ROIS


class FileNotCachedError(Exception):
    def __init__(self, filename: str, plate: str | None = None, *args):
        super().__init__(*args)
        self.filename = filename
        self.plate = plate


# TODO: Label cache from data source used
# TODO: Cache integrity checks
class Cache:
    def __init__(self, cache_directory: str):
        self.root: pathlib.Path = pathlib.Path(cache_directory)
        if not self.root.exists():
            self.root.mkdir(parents=True)

    def load_file(
        self,
        filename: str,
        plate: str | None = None,
        genfile_cb: Callable | None = None,
        return_file_object: bool = True,
    ):
        local_path = self.file_path(filename, plate)
        logger.debug(f"Trying to load file {local_path}.")
        if not local_path.is_file():
            if genfile_cb is not None:
                logger.debug(
                    f"File {str(local_path.relative_to(self.root))} is not "
                    "cached, trying to generate it."
                )
                genfile_cb()
            else:
                logger.warning(
                    f"File {str(local_path.relative_to(self.root))} is not "
                    "cached and we have no way to generate it."
                )
                raise FileNotCachedError(filename, plate)
        else:
            logger.debug(
                f"File {str(local_path.relative_to(self.root))} is cached, "
                "using local version."
            )

        if return_file_object:
            return local_path.open("r", encoding='utf-8-sig')
        return local_path

    def file_path(self, filename: str, plate: str | None = None):
        if plate is None:
            path = self.root / filename
        else:
            path = (self.root / plate / filename)
        path.parent.mkdir(exist_ok=True)
        return path

    def temporary_file(self, suffix="", mode="wb", delete=False):
        tempdir = self.root / "_tmp"
        if not tempdir.exists():
            tempdir.mkdir()
        return tempfile.NamedTemporaryFile(
            suffix=suffix,
            mode=mode,
            delete=delete,
            dir=tempdir,
        )

    def lock_file_name(self, item: str):
        path = self.root / "lock"
        path.mkdir(exist_ok=True)
        return FileLock(path / item)

    def extract_plate_zip(self, plate: str, zipfilepath: pathlib.Path):
        logger.debug(f"Extracting zip file {zipfilepath} for plate {plate}.")
        with zipfile.ZipFile(zipfilepath) as zip:
            fields = [
                pathlib.Path(fname.replace(FILE_PATTERN, ""))
                for fname in zip.namelist()
                if fname.endswith(FILE_PATTERN)
            ]
            fields = [
                fname for fname in fields
                if not (
                    fname.stem.startswith("Control") or
                    fname.stem.startswith("ontrol") or
                    fname.stem.startswith("control")
                )
            ]
            for field in tqdm(
                fields,
                desc=f"Extracting {plate}",
                disable=logger.getEffectiveLevel() != logging.INFO,
            ):
                for filetype in FileTypes:
                    zippath = pathlib.Path(str(field) + filetype)
                    zipinfo = zip.getinfo(str(zippath))

                    filename = str(zippath.relative_to(zippath.parent))

                    zipinfo.filename = str(self.file_path(filename, plate))
                    logger.debug(f"Extracting file {zippath} to "
                                 f"{zipinfo.filename}.")
                    zip.extract(zipinfo)

    def glob_files(self, pattern: str):
        return [p.relative_to(self.root) for p in self.root.glob(pattern)]

    def delete(self):
        shutil.rmtree(self.root)

    def __hash__(self):
        return hash(("Cache", self.root))
