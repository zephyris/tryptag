from __future__ import annotations
import logging
import pathlib
import shutil
import sys
import tempfile
from typing import Callable, TYPE_CHECKING
import zipfile

from filelock import FileLock
from tqdm import tqdm

if TYPE_CHECKING:
    from tryptag.datasource import DataSource

if sys.version_info >= (3, 11):
    from enum import StrEnum
else:
    from backports.strenum import StrEnum

logger = logging.getLogger("tryptag.cache")


class FileTypes(StrEnum):
    CELL_ROIS = "_roisCells.txt"
    THRESHOLDED = "_thr.tif"
    IMAGE = ".tif"
    LOCALISATION = ".json"
    OTHER = ""


FILE_PATTERN = FileTypes.CELL_ROIS


class FileNotCachedError(Exception):
    def __init__(self, filename: str, plate: str | None = None, *args):
        super().__init__(*args)
        self.filename = filename
        self.plate = plate


class IncompatibleCacheError(Exception):
    pass


# TODO: Label cache from data source used
# TODO: Cache integrity checks
class Cache:
    """Local filesystem cache object. This handles low-level loading of files
    and management of the cache directory."""
    def __init__(self, cache_directory: str):
        """
        Initialise the cache object.

        :param cache_directory: str, path to the local cache. This will be
            created if it does not exist. Insufficient permissions will lead
            to an error.
        """
        self.root: pathlib.Path = pathlib.Path(cache_directory)
        if not self.root.exists():
            self.root.mkdir(parents=True)

    def verify_cache_datasource(self, datasource: DataSource):
        stamp = str(type(datasource).__name__)

        def genfile_cb():
            with self.file_path("cache_origin").open("w") as f:
                f.write(stamp)
        with self.load_file(
            "cache_origin",
            genfile_cb=genfile_cb,
        ) as f:
            cache_origin = f.read()

        if cache_origin != stamp:
            raise IncompatibleCacheError(
                f"The cache at '{self.root}' has been created using the "
                f"'{cache_origin}' data source but you are trying to use it "
                f"with the '{stamp}' data source. Please either delete the "
                "cache or specify a different cache directory."
            )

    def load_file(
        self,
        filename: str,
        plate: str | None = None,
        genfile_cb: Callable | None = None,
        return_file_object: bool = True,
    ):
        """
        Low-level method to load a file from the cache. It will fall back to a
        given callback function to potentially fetch the file if it isn't
        found in the cache.

        :param filename: str, path (within the cache) to the file to be loaded
        :param plate: str or None, name of the plate the file resides in (or
            None if it is a root file)
        :param genfile_cb: function that is called if the file doesn't exist
            locally. This should take no argument and return None.
        :param return_file_object: bool (default `True`), whether to return an
            opened file object. Returns the path to the local file if `False.
        """
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
        """
        Return the local path to the given file.

        :param filename: str, name of the file
        :param plate: str or None, name of the plate or None if the file is a
            root file
        :return: absolute local path to the file (as a pathlib.Path object)
        """
        if plate is None:
            path = self.root / filename
        else:
            path = (self.root / plate / filename)
        path.parent.mkdir(exist_ok=True)
        return path

    def is_cached(self, filename: str, plate: str | None = None):
        local_path = self.file_path(filename, plate)
        return local_path.is_file()

    def temporary_file(self, suffix="", mode="wb", delete=False):
        """
        Create a local temporary file object within the cache (e.g. to
        download and temporarily store a ZIP file in).

        :param suffix: str (default ""), suffix of the file name
        :param mode: str (default "wb"), mode under which the file should be
            opened
        :param delete: bool (default `False`), whether the file should be
            deleted when it is closed. If `False`, the user must delete it.
        :return: opened file object
        """
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
        """
        Return the lock object for an item.

        :param item: str, name of the item to be locked
        :return: FileLock object
        """
        path = self.root / "lock"
        path.mkdir(exist_ok=True)
        return FileLock(path / item)

    def extract_plate_zip(self, plate: str, zipfilepath: pathlib.Path):
        """
        Extract a downloaded plate zip file into the cache (this is Zenodo
        specific).

        :param plate: str, name of the plate
        :param zipfilepath: Path, path to the zip file to be extracted.
        """
        # TODO: We should move this to the Zenodo class.
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
                for filetype in [
                    FileTypes.CELL_ROIS,
                    FileTypes.IMAGE,
                    FileTypes.THRESHOLDED,
                ]:
                    zippath = pathlib.Path(str(field) + filetype)
                    zipinfo = zip.getinfo(str(zippath))

                    filename = str(zippath.relative_to(zippath.parent))

                    zipinfo.filename = str(
                        self.file_path(filename, plate).relative_to(self.root))
                    logger.debug(f"Extracting file {zippath} to "
                                 f"{zipinfo.filename}.")
                    zip.extract(zipinfo, path=self.root)

    def glob_files(self, pattern: str):
        """
        Find and return a list of files matching the given pattern in the
        cache.

        :param pattern: str, a shell-like glob pattern
        :return: a list of files matching the pattern
        """

        return [p.relative_to(self.root) for p in self.root.glob(pattern)]

    def delete(self):
        """
        Delete the cache. This is non-reversible.
        """
        shutil.rmtree(self.root)

    def __hash__(self):
        return hash(("Cache", self.root))
