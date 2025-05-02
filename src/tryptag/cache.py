from enum import Enum
import logging
import pathlib
import shutil
import tempfile
from typing import Callable
import zipfile

from filelock import FileLock
from tqdm import tqdm


class FileTypes(Enum):
    CELL_ROIS = 0
    THRESHOLDED = 1
    IMAGE = 2


FILE_TYPE_ENDINGS = {
    FileTypes.CELL_ROIS: "_roisCells.txt",
    FileTypes.THRESHOLDED: "_thr.tif",
    FileTypes.IMAGE: ".tif",
}
FILE_PATTERN = FILE_TYPE_ENDINGS[FileTypes.CELL_ROIS]

logger = logging.getLogger(__name__)


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
            ):
                for filetype in FILE_TYPE_ENDINGS.values():
                    zippath = pathlib.Path(str(field) + filetype)
                    zipinfo = zip.getinfo(str(zippath))

                    filename = str(zippath.relative_to(zippath.parent))

                    zipinfo.filename = str(self.file_path(filename, plate))
                    zip.extract(zipinfo)

    def glob_files(self, pattern: str):
        return [p.relative_to(self.root) for p in self.root.glob(pattern)]

    def delete(self):
        shutil.rmtree(self.root)

    def __hash__(self):
        return hash(("Cache", self.root))
