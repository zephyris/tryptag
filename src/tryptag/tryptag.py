from __future__ import annotations
import logging
from typing import Callable, Literal
import warnings

from tryptag.processing import WorkList

from .cache import Cache
from .datasource import CellLine, CellLineStatus, DataSource
from .zenodo import Zenodo
from .images import FieldImage, CellImage

logger = logging.getLogger("tryptag")

STANDARD_DATASOURCES = {
    "procyclic": lambda cache: Zenodo(cache, master_record_id=6862289),
    "bloodstream": lambda cache: Zenodo(cache, master_record_id=7258722)
}


class TrypTag:
    """TrypTag API class."""

    datasource: DataSource

    def __init__(
        self,
        verbose: bool = False,
        data_cache_path: str = "./_tryptag_cache",
        um_per_px: float = 6.5 / 63,
        dataset_name: str | None = "procyclic",
        datasource: DataSource | None = None,
        life_stages: list[str] | None = None,  # deprecated
    ):
        """Initialise TrypTag API. This class can be instantiated multiple
        times for different data sets (e.g. the procyclic and bloodstream
        form data sets).

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

        if verbose:
            logger.setLevel(logging.DEBUG)

        if life_stages is not None:
            warnings.warn(
                "The life_stages argument is deprecated. Please "
                "either use dataset_name or explicitly specify a data source."
            )
            if len(life_stages) > 1:
                warnings.warn(
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
            logger.debug(
                f"Initialising cache and standard data source {dataset_name}."
            )
            cache = Cache(data_cache_path)
            self.datasource = STANDARD_DATASOURCES[dataset_name](cache)
        elif datasource is not None:
            logger.debug(
                "Initialising from given data source (assuming cache has "
                "been initialised)."
            )
            self.datasource = datasource

        # image properties
        self.um_per_px = um_per_px

        # global variables for caching last field of view loaded
        self._field_base_path_sk = None
        self._thresholds_sk = None
        self._channels_sk = None

    @property
    def gene_list(self):
        """The genes in the data set. This is a mapping, so the list of genes
        can be extracted via
        ```python
        genes = list(tt.gene_list)
        ```
        and individual genes can be accessed by
        ```python
        gene = tt.gene_list["Tb927.9.8570"]
        ```
        Individual cell lines can then be accessed through the "N" and "C"
        attributes:
        ```python
        cell_line = gene.C
        """
        return self.datasource.gene_collection

    @property
    def life_stages(self):
        warnings.warn("life_stages is deprecated.")
        return [None]

    @property
    def cell_lines(self):
        """A list of all cell lines that have been considered (even ones that
        were not successfully generated!)."""
        return WorkList(
            self,
            [
                cell_line
                for gene in self.gene_list.values()
                for cell_line in [gene.C, gene.N]
            ]
        )

    def worklist_all(self, life_stage: str | None = None):
        """
        All `gene_id` and `terminus` combinations with data, as a `CellLine`
        object containing `life_stage`, `gene_id` and `terminus`. In contrast
        to `TrypTag.cell_lines`, this only returns cell lines that have
        been successfully generated.
        """
        if life_stage is not None:
            warnings.warn("life_stage is deprecated and ignored.")

        return self.cell_lines.filter(
            lambda cell_line: cell_line.status == CellLineStatus.GENERATED
        )

    def worklist_parental(self, life_stage: str | None = None) -> list:
        """
        All dummy `gene_id` and `terminus` combinations which correspond to
        parental cell line samples, as a `CellLine` object containing
        `life_stage`, `gene_id` and `terminus`.
        """
        if life_stage is not None:
            warnings.warn("life_stage is deprecated and ignored.")

        return self.worklist_all().filter(
            lambda cell_line: cell_line.is_parental
        )

    def localisation_search(
            self,
            query_term: str,
            life_stage: str | None = None,
            match_subterms: bool = True,
            exclude_modifiers: list[str] = ["weak", "<10%"],
            required_modifiers: list[str] | None = None
    ):
        """
        Get a worklist of `gene_id` and `terminus` hits where any of the
        localisation annotations match the query.

        :param query_term: Search query annotation term from the localisation
            ontology.
        :param match_subterms: Whether to also match
            child/subterm/substructures of `query_term`, default `True`.
        :param exclude_modifiers: List of modifier terms none of which can be
            matched, default `"weak"` and `"<10%"`.
        :param required_modifiers: List of modifier terms all of which must be
            matched, default `None`.
        :return: List of `CellLine` objects of the hits, containing
            `gene_id` and `terminus`.
        """
        return self.worklist_all().localisation_search(
            query_term,
            match_subterms=match_subterms,
            exclude_modifiers=exclude_modifiers,
            required_modifiers=required_modifiers,
        )

    @property
    def localisation_ontology(self):
        """The localisation term ontology."""
        ontology = {}
        entries = self.datasource.localisation_ontology.entries
        for name, entry in entries.items():
            oentry = {}
            if entry.synonyms is not None:
                oentry["synonyms"] = list(entry.synonyms)
            if entry.comment is not None:
                oentry["comment"] = entry.comment
            if entry.ident is not None:
                oentry["ident"] = entry.ident
            if entry.goterm is not None:
                oentry["go"] = entry.goterm
            current = entry.parent
            ancestors = []
            while current is not None:
                ancestors.append(current.name)
                current = current.parent
            oentry["parent"] = ["root"] + ancestors[::-1]
            if entry.children:
                oentry["children"] = [child.name for child in entry.children]
            if entry.examples:
                oentry["examples"] = list(entry.examples)
            ontology[name] = oentry
        return ontology

    def gene_id_search(
            self,
            gene_id_list: list,
            life_stage: str | None = None
    ) -> list:
        """
        Use a list of gene ids to build a worklist of all tagged termini for
        that gene.

        :param gene_id_list: List of strings which are gene IDs to search for.
        """
        self.worklist_all().gene_id_search(gene_id_list)

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
        if not cell_line.initialised:
            cell_line = self.gene_list[cell_line.gene_id][cell_line.terminus]
        logger.debug(f"Opening field for cell line {cell_line}")

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
        if not cell_line.initialised:
            cell_line = self.gene_list[cell_line.gene_id][cell_line.terminus]

        field = cell_line.fields[field_index]
        cell = field.cells[cell_index]
        logger.debug(f"Opening field for cell {cell_line} field {field} "
                     f"cell {cell}.")

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
        warnings.warn("open_cell_custom is deprecated. Use open_cell instead.")
        return self.open_cell(*args, **kwargs)

    def analyse_list(
        self,
        work_list: WorkList,
        analysis_function: Callable,
        workers: int | None = None,
        multiprocess_mode: Literal["process", "thread"] | None = "process"
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
        return work_list.process(
            analysis_function,
            workers,
            multiprocess_mode
        )

    def cell_list(self, cell_line: CellLine):
        """Return a list of cells from the given `CellLine`."""
        return [
            cell
            for field in cell_line.fields.values()
            for cell in field.cells.values()
        ]
