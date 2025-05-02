from __future__ import annotations
import logging
from typing import Callable

from tqdm.auto import tqdm

from .cache import Cache
from .datasource import TERMINI, CellLine, DataSource
from .zenodo import Zenodo
from .images import FieldImage, CellImage

logger = logging.getLogger(__name__)


class _tqdmDownload(tqdm):
    def urllib_callback(self, transferred_blocks, block_size, total_size):
        self.total = total_size
        return self.update(transferred_blocks * block_size - self.n)


STANDARD_DATASOURCES = {
    "procyclic": lambda cache: Zenodo(cache, master_record_id=6862289),
    "bloodstream": lambda cache: Zenodo(cache, master_record_id=7258722)
}


class TrypTag:
    """TrypTag abstract base class, this should never be instantiated."""

    datasource: DataSource

    def __init__(
        self,
        verbose: bool = True,
        data_cache_path: str = "./_tryptag_cache",
        um_per_px: float = 6.5 / 63,
        dataset_name: str | None = "procyclic",
        datasource: DataSource | None = None,
        life_stages: list[str] | None = None,  # deprecated
    ):
        """Initialise TrypTag Base object.

        Keyword arguments:
        verbose -- print verbose output from accessing data (default `True`)
        data_cache_path -- the directory that will hold the downloaded TrypTag
            data (default `./_tryptag_cache`, requires up to 8Tb)
        um_per_px -- physical pixel size / corrected magnification
        dataset_name -- str, the name of a known dataset (currently
            "procyclic" or "bloodstream"), default "procyclic"
        data_source -- DataSource object to manually specify the data source.
            `life_stage` needs to be set to `None` for this to have an effect.
        life_stages -- deprecated.
        """

        # TODO: Make logging object-local instead of global
        if verbose:
            logger.setLevel(logging.DEBUG)

        if life_stages is not None:
            logger.warning(
                "The life_stages argument is deprecated. Please "
                "either use dataset_name or explicitly specify a data source."
            )
            if len(life_stages) > 1:
                logger.warning(
                    "We only support specifying one life stage. Only using "
                    "the first one as the data set name."
                )
            elif len(life_stages) == 0:
                raise ValueError(
                    "If specified, life_stages needs to have exactly one "
                    "entry."
                )
            elif dataset_name is not None:
                raise ValueError(
                    "life_stages and dataset_name cannot be specified at the "
                    "same time."
                )
            elif datasource is not None:
                raise ValueError(
                    "life_stages and data_source cannot be specified at the "
                    "same time."
                )
            dataset_name = life_stages[0]
        if dataset_name is not None and datasource is not None:
            raise ValueError(
                "life_stage and data_source cannot be specified at the same "
                "time."
            )
        elif dataset_name is not None:
            cache = Cache(data_cache_path)
            self.datasource = STANDARD_DATASOURCES[dataset_name](cache)
        elif datasource is not None:
            self.datasource = datasource

        # image properties
        self.um_per_px = um_per_px

        # global variables for caching last field of view loaded
        self._field_base_path_sk = None
        self._thresholds_sk = None
        self._channels_sk = None

    @property
    def gene_list(self):
        return self.datasource.gene_collection

    def worklist_all(self, life_stage: str | None = None) -> list[CellLine]:
        """
        All `gene_id` and `terminus` combinations with data, as a `CellLine`
        object containing `life_stage`, `gene_id` and `terminus`.

        :param life_stage: Which life cycle stage to load data from, default
        is first entry in `self.life_stages`.
        """
        if life_stage is not None:
            logger.warning("life_stage is deprecated and ignored.")

        return [
            gene_entry[terminus]
            for gene_id, gene_entry in self.gene_list.items()
            for terminus in TERMINI
        ]

    def worklist_parental(self, life_stage: str | None = None) -> list:
        """
        All dummy `gene_id` and `terminus` combinations which correspond to
        parental cell line samples, as a `CellLine` object containing
        `life_stage`, `gene_id` and `terminus`.
        """
        if life_stage is not None:
            logger.warning("life_stage is deprecated and ignored.")

        return [
            gene_entry[terminus]
            for gene_id, gene_entry in self.gene_list.items()
            for terminus in TERMINI
            if "wild-type" in gene_id
        ]

    def localisation_search(
            self,
            query_term: str,
            life_stage: str = None,
            match_subterms: bool = True,
            exclude_modifiers: list = ["weak", "<10%"],
            required_modifiers: list = None
    ) -> list:
        """
        Get a worklist of `gene_id` and `terminus` hits where any of the
        localisation annotations match the query.

        :param life_stage: Life cycle stage, default to first entry in
            `self.life_stages`.
        :param query_term: Search query annotation term from the localisation
            ontology.
        :param match_subterms: Whether to also match
            child/subterm/substructures of `query_term`, default `True`.
        :param exclude_modifiers: List of modifier terms none of which can be
            matched, default `"weak"` and `"<10%"`.
        :param required_modifiers: List of modifier terms all of which must be
            matched, default `None`.
        :return: List of `CellLine` objects of the hits, containing
            `life_stage`, `gene_id` and `terminus`.
        """
        # determine life stage
        if life_stage is None:
            life_stage = next(iter(self.datasources))

        # check all against query
        hits = []
        for cell_line in self.worklist_all(life_stage):
            if cell_line.localisation.match(
                    query_term,
                    require_modifiers=required_modifiers,
                    exclude_modifiers=exclude_modifiers,
                    recursive=match_subterms,
            ):
                hits.append(cell_line)
        return hits

    def gene_id_search(
            self,
            gene_id_list: list,
            life_stage: str = None
    ) -> list:
        """
        Use a list of gene ids to build a worklist of all tagged termini for
        that gene.

        :param gene_id_list: List of strings which are gene IDs to search for.
        :param query_term: Search query annotation term from the localisation
            ontology.
        """
        # build list
        hits = []
        for cell_line in self.worklist_all(life_stage):
            if cell_line.gene_id in gene_id_list:
                hits.append(cell_line)
        return hits

    def open_field(
            self,
            cell_line: CellLine,
            field_index: int = 0,
            custom_field_image: FieldImage | None = None
    ):
        """
        Returns field of view image data.

        :param cell line: `CellLine` object containing `life_stage`, `gene_id`
            and `terminus`.
        :param field_index: Index of the field of view. If not set, then `0`.
        :param custom_field_image: `FieldImage` object containing custom field
            images to use. Images can be skimage image or `None`. Entries of
            `None` will use tryptag default. If not set or `None`, then use
            all tryptag defaults.
        :return: FieldImage object, containing the image channels as
            attributes `phase`, `mng`, `dna`, and the thresholds `phase_mask`
            and `dna_mask`.
        """
        return FieldImage.from_field(
            cell_line.fields[field_index],
            custom_field_image=custom_field_image
        )

    # open a cell, cropped from a field of view
    # uses the phase and dna threshold images from tryptag
    # cell x, y coordinate in phase threshold from tryptag
    def open_cell(
        self,
        cell_line: CellLine,
        field_index: int = 0,
        cell_index: int = 0,
        width: int = 323,
        rotate: bool = False,
        fill_centre: tuple[int, int] | None = None,
        crop_centre: tuple[int, int] | None = None,
        angle: float | None = None,
        custom_field_image: FieldImage | None = None,
    ):
        """
        Opens a cell from a `gene_id`, `terminus`, `field_index` and
        `cell_index`. Use `open_cell_custom` if you would like to inject
        custom images, coordinates and/or angle.

        :param cell line: `CellLine` object containing `life_stage`, `gene_id`
            and `terminus`.
        :param field_index: Index of the field of view. If not set, then `0`.
        :param cell_index: Index of the cell in the field of view. If not set,
            then `0`.
        :param width: If positive, width of cropped cell image (may clip very
            large cells). If negative, padding for crop around the
            `phase_mask`. Default, `323`.
        :param fill_centre: `(x, y)` tuple of a pixel which is in the target
            cell object (pixel value 255) in pth image.
        :param crop_centre: `(x, y)` tuple around which to crop, otherwise
            crop around `fill_centre`.
        :param rotate: Whether or not to rotate the cell. Default `False`. Set
            to `False` if `width < 0` (padded crop mode).
        :param angle: Angle in degrees to rotate cell clockwise. If not set or
            `None`, tryptag default.
        :param custom_field_image: `FieldImage` object containing custom field
            images to use. Images can be skimage image or `None`. Entries of
            None will use tryptag default. If not set or `None`, then use all
            tryptag defaults.
        :return: CellImage object, containing the image channels as attributes
            `phase`, `mng`, `dna`, and the thresholds `phase_mask`, `dna_mask`
            and `phase_mask_othercells`.
        """
        field = cell_line.fields[field_index]
        cell = field.cells[cell_index]
        return CellImage(
            cell,
            rotated=rotate,
            width=width,
            fill_centre=fill_centre,
            crop_centre=crop_centre,
            angle=angle,
            custom_field_image=custom_field_image,
        )

    def open_cell_custom(self, *args, **kwargs):
        """
        DEPRECATED - use `open_cell` instead.
        """
        return self.open_cell(*args, **kwargs)

    def _list_analysis_worker(
        self,
        cell_line: CellLine,
        analysis_function: callable,
        threading_mode: bool
    ) -> list:
        """
        Worker for multiprocess/thread parallel analysis of a `work_list`.

        :param cell_line: `CellLine` object defining the `life_stage`,
            `gene_id` and `terminus` combination to analyse.
        :param analysis_function: Function name to use for analysis.
            `analysis_function` should take exactly two arguments, `tryptag`
            (`TrypTag` instance) and `cell_line` (`CellLine` object) in this
            order.
        :param threading_mode: True if running in threading mode. This
            triggers a deep copy of `self` to avoid thread safety issues.
        :return: List of dicts in the form `{"life_stage": life_stage,
            "gene_id": gene_id, "terminus": terminus, "result":
            analysis_function_return}`.
        """
        # Deep copy tryptag object if running in threading mode (to avoid
        # thread safety issues)
        # TODO: I doubt this is needed anymore, but need to test!
        if threading_mode:
            from copy import deepcopy
            self = deepcopy(self)

        result = {
            "cell_line": cell_line,
            "result": analysis_function(self, cell_line),
        }
        return result

    def analyse_list(
        self,
        work_list,
        analysis_function,
        workers=None,
        multiprocess_mode="process"
    ):
        """
        Simple handler for parallel analysis of a `work_list`.

        :param work_list: List of CellLine objects, defining each `life_stage,
            `gene_id` and `terminus` combination to analyse.
        :param analysis_function: Function name to use for analysis.
            `analysis_function` should take exactly two arguments, `tryptag`
            (`TrypTag` instance) and `cell_line` (`CellLine` object) in this
            order.
        :param workers: Number of threads/processes to spawn, default is
            number of CPUs.
        :param multiprocess_mode: `"process"` for parallel processes,
            `"thread"` for parallel threads or `None` for no parallel
            processing (directly calls `analysis_function`).
        :return: List of dicts in the form `{"life_stage": life_stage,
            "gene_id": gene_id, "terminus": terminus, "result":
            analysis_function_return}`. These may be in a different order to
            `work_list`.
        """
        import concurrent.futures
        import multiprocessing

        logger.info("Analysing worklist")

        # deduplicate work_list
        dedup_work_list = set(work_list)

        # get number of workers, default to number of cpus
        if workers is None:
            workers = multiprocessing.cpu_count()

        if multiprocess_mode is None:
            # run in a single thread, still use the ThreadPoolExecutor since
            # that's equivalent
            logger.info("Single process")
            Executor = concurrent.futures.ThreadPoolExecutor
            workers = 1
        elif multiprocess_mode == "process":
            # setup executor as a process pool
            logger.info("Parallel processes with", workers, "workers")
            Executor = concurrent.futures.ProcessPoolExecutor
        elif multiprocess_mode == "thread":
            # setup executor as a thread pool
            logger.info("Parallel threads with", workers, "workers")
            Executor = concurrent.futures.ThreadPoolExecutor
        else:
            raise ValueError(f"Unknown multiprocess_mode '{multiprocess_mode}")

        with Executor(workers) as executor:
            futures = [
                executor.submit(
                    self._list_analysis_worker,
                    cell_line=cell_line,
                    analysis_function=analysis_function,
                    threading_mode=multiprocess_mode == "thread"
                ) for cell_line in dedup_work_list]
            results = [
                future.result()
                for future in tqdm(
                    concurrent.futures.as_completed(futures),
                    total=len(futures),
                    smoothing=0
                )
            ]
        return results
