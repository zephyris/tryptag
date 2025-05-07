from __future__ import annotations
from collections.abc import Mapping
import csv
from dataclasses import dataclass
from functools import cached_property
import json
import logging
import pathlib
import sys
from typing import Literal
import warnings

from .cache import Cache, FileTypes, FILE_PATTERN

from annotations import (
    Ontology,
    OntologyEntry,
    OntologyAnnotationCollection,
)

if sys.version_info >= (3, 11):
    from enum import StrEnum
else:
    from backports.strenum import StrEnum

logger = logging.getLogger("tryptag.datasource")


class HeadersGene(StrEnum):
    """Column headers that document gene properties in localisations.tsv"""

    GENE_ID = "Gene ID"
    GENE_ALIASES = "Gene aliases"


class HeadersCellLine(StrEnum):
    """Column headers that document cell line properties in
    localisations.tsv"""
    STATUS = "status"
    PLATE_AND_WELL = "plate and well"
    PRIMER_F = "primer F"
    PRIMER_R = "primer R"
    LOCALISATION = "localisation"
    FAINTER_THAN_PARENTAL = "fainter than parental"
    CLASSIFIED_FAINT = "classified as faint"


class CellLineStatus(StrEnum):
    """Status codes for the cell line generated in the TrypTag project"""

    NOT_ATTEMPTED = "not attempted"
    """Cell lines that haven't been attempted."""

    ATTEMPTED = "tagging attempted"
    """Cell lines that have been attempted but weren't successfully
    generated."""

    GENERATED = "cell line generated"
    """Cell lines that were successfully generated and imaged."""


class Cell:
    """Class describing metadata for a single cell in a `CellLine`'s
    `Field`."""

    def __init__(
        self,
        field: Field,
        index: int,
        wand: tuple[int, int],
        centre: tuple[float, float],
        extent: tuple[float, float],
        angle: float
    ):
        """
        Initialise a `Cell`.

        :param field: `Field` in which the `Cell` resides.
        :param index: int, Index of the cell in the field.
        :param wand: tuple[int, int], xy position that's guaranteed to be in
            the cell (for flood filling).
        :param centre: tuple[float, float], xy position of the centre of the
            cell
        :param extent: tuple[float, float], width and height of the cell
        :param angle: float, angle the cell subtends with the x axis, taking
            into account the cell anterior-posterior orientation
        """
        self.field = field
        self.index = index
        self.wand = wand
        self.centre = centre
        self.extent = extent
        self.angle = angle

    @staticmethod
    def from_line(field: Field, line: str):
        """
        Initialise a `Cell` from a line in the cell ROIs file for the given
        `Field`.
        """
        parts = line.split()
        return Cell(
            field,
            int(parts[0]),
            (int(parts[1]), int(parts[2])),
            (float(parts[3]), float(parts[4])),
            (float(parts[5]), float(parts[6])),
            float(parts[7])
        )

    def __hash__(self):
        return hash((self.field, self.index))

    def __getitem__(self, key):
        warnings.warn(
            "Accessing cell properties via keys is deprecated."
        )
        if key == "field_index":
            if self.field is not None:
                return self.field.index
            else:
                return None
        elif key == "cell_index":
            return self.index

    def __str__(self):
        return f"Cell from {self.field}"


class FieldDoesNotExistError(Exception):
    def __init__(
        self,
        cell_line: CellLine,
        index: int,
    ):
        self.cell_line = cell_line
        self.index = index


class Field:
    """Class describing metadata for a `CellLine`'s field of view."""
    def __init__(
        self,
        cell_line: CellLine,
        index: int,
        datasource: DataSource,
    ):
        """
        Initialise a `Field` object.

        :param cell_line`: the cell line the field of view shows
        :param index: int, the index of the field of view
        :param datasource: DataSource, the data source the main `TrypTag`
            object is using (needed for data access).
        """
        self.cell_line = cell_line
        self.index = index
        self.datasource = datasource

        if not self.exists():
            raise FieldDoesNotExistError(self.cell_line, self.index)

        cell_file = datasource.load_plate_file(
            cell_line.plate,
            self.filename(FileTypes.CELL_ROIS)
        )
        with cell_file:
            # Skip header line
            next(cell_file)
            cell_list = [Cell.from_line(self, line) for line in cell_file]
            self.cells = {cell.index: cell for cell in cell_list}

    def filename(self, file_type: FileTypes):
        """Return the name of the given file type for this field of view."""
        return (
            f"{self.cell_line.filename_stem()}{self.index+1}{file_type}"
        )

    def exists(self) -> bool:
        """Check that this field exists. Please note that this will trigger a
        download if the files are not found in the cache.
        """
        path: pathlib.Path = self.datasource.load_plate_file(
            self.cell_line.plate,
            self.filename(FileTypes.CELL_ROIS),
            return_file_object=False
        )
        return path.is_file()

    def __hash__(self):
        return hash((self.cell_line, self.index))

    def __repr__(self):
        return f"{self.cell_line} index = {self.index}"


class CellLine:
    """Class describing a cell line within the TrypTag project."""
    gene: Gene = None
    _gene_id: str = ""
    terminus: Literal["N", "C"]
    status: CellLineStatus
    plate: str
    well: str
    forward_primer: str
    reverse_primer: str
    localisation: OntologyAnnotationCollection
    fainter_than_parental: bool
    classified_faint: bool
    datasource: DataSource
    life_stage = None

    def __init__(
        self,
        gene_id: str,
        terminus: str
    ):
        """
        Create a `CellLine` object.

        Note that this constructor only fills the `gene_id` and the
        `terminus`. It does not perform data lookup.

        :param gene_id: str, Gene ID of the cell line.
        :param terminus: str, terminus of the tag on the cell line (either "C"
            or "N")
        """
        self._initialised = False
        self._gene_id = gene_id
        terminus = terminus.upper()
        if terminus not in ["N", "C"]:
            raise ValueError("terminus needs to be either 'N' or 'C'")
        self.terminus = terminus

    @property
    def initialised(self):
        """Whether the `CellLine` object has been populated with data from the
        data source."""
        return self._initialised

    @property
    def gene_id(self):
        """The Gene ID of the cell line."""
        if self.gene is not None:
            return self.gene.id
        return self._gene_id

    @property
    def is_parental(self):
        """Whether the `CellLine` stems from a parental tagging attempt (i.e.
        not tagging a specific protein)."""
        return "wild-type" in self.gene_id

    @staticmethod
    def from_data(
        terminus: Literal["N", "C"],
        status: CellLineStatus,
        plate_and_well: str,
        forward_primer: str,
        reverse_primer: str,
        localisation: str,
        fainter_than_parental: bool,
        classified_faint: bool,
        datasource: DataSource,
        ontology: Ontology,
    ):
        """Initialise a `CellLine` object and populate with the given data.

        :param terminus: the tagging terminus of the cell line (either "C" or
            "N")
        :param status: the status of the cell line (has to be one of the
            entries in `CellLineStatus`)
        :param plate_and_well: str, the plate and well the cell line was
            generated in (separated by a space)
        :param forward_primer: str, the sequence of the forward primer
        :param reverse_primer: str, the sequence of the reverse primer
        :param localisation: str, the localisation annotation string in the
            form "term1[mod1,mod2],term2[mod3,mod4]"
        :param fainter_than_parental: bool, whether the fluorescence signal
            was fainter than the parental cell lines
        :param classifed_faint: bool, whether the fluorescence signal was
            classified as faint
        :param datasource: `DataSource, the data source the `TrypTag` instance
            is using
        :param ontology: `Ontology`, the localisation annotation ontology the
            annotation terms are derived from
        """
        self = CellLine("", terminus)
        self.terminus = terminus
        self.status = CellLineStatus(status)
        if self.status != CellLineStatus.NOT_ATTEMPTED:
            self.plate, self.well = plate_and_well.split(" ")
        else:
            self.plate, self.well = "", ""
        self.forward_primer = forward_primer
        self.reverse_primer = reverse_primer
        self.localisation = OntologyAnnotationCollection(
            localisation,
            ontology,
        )
        self.fainter_than_parental = fainter_than_parental
        self.classified_faint = classified_faint
        self.datasource = datasource

        self._initialised = True

        return self

    def __hash__(self):
        return hash(
            (self.datasource, self.gene.id, self.terminus, self.life_stage))
    
    def __eq__(self, other: CellLine):
        return (
            self.gene_id == other.gene_id and
            self.terminus == other.terminus
        )

    def __repr__(self):
        gene_id = self.gene.id if self.gene is not None else self.gene_id
        return (
            f"gene_id = {gene_id} "
            f"terminus = {self.terminus} "
        )

    def filename_stem(self):
        return f"{self.gene.id}_4_{self.terminus}_"

    @cached_property
    def fields(self) -> dict[int, Field]:
        """The `Field`s of view belonging to this cell line."""
        if self.status != CellLineStatus.GENERATED:
            return {}
        stem = self.filename_stem()
        files = self.datasource.glob_plate_files(
            self.plate,
            f"{stem}*{FILE_PATTERN}",
        )
        indices = [
            int(fname.replace(stem, "").replace(FILE_PATTERN, "")) - 1
            for fname in files
        ]
        return {
            index: Field(self, index, self.datasource)
            for index in sorted(indices)
        }

    def __getitem__(self, key):
        warnings.warn(
            "Accessing cell line properties via keys is deprecated."
        )
        if key == "loc":
            return [
                {
                    "term": localisation.term,
                    "modifiers": list(localisation.modifiers),
                }
                for localisation in self.localisation
            ]
        elif key == "plate":
            return self.plate
        elif key == "well":
            return self.well
        elif key == "primer_f":
            return self.forward_primer
        elif key == "primer_r":
            return self.reverse_primer
        raise KeyError(f"unknown key {key}")


@dataclass
class Gene(Mapping):
    """Class describing the metadata for a gene that was part of the TrypTag
    project."""
    id: str
    aliases: str
    N: CellLine
    C: CellLine

    def __getitem__(self, terminus: str):
        terminus = terminus.upper()
        if terminus == "N":
            return self.N
        elif terminus == "C":
            return self.C
        else:
            raise KeyError(f"Terminus {terminus} unknown.")

    def __iter__(self):
        return iter(["C", "N"])

    def __len__(self):
        return 2


class GeneCollection(Mapping):
    """Collection holding the `Gene` objects."""
    def __init__(self, genes: dict[str, Gene]):
        """
        Initialise a `GeneCollection` object.

        :param genes: dict mapping gene IDs with `Gene` objects.
        """
        self.genes = genes

    def __getitem__(self, geneid: str) -> Gene:
        if geneid == "procyclic" or geneid is None:
            warnings.warn(
                "Specifying a life stage is deprecated."
            )
            return self
        return self.genes[geneid]

    def __iter__(self):
        return iter(self.genes)

    def __len__(self):
        return len(self.genes)


class DataSource:
    """The abstract data source class allowing the `TrypTag` class to fetch
    data."""
    def __init__(self, cache: Cache):
        """
        Initialise the `DataSource`.

        Please note that this is the abstract base class - this constructor
        should not be called on its own. When subclassing, this constructor
        needs to be called before anything else in the subclass' consttructor.

        :param cache: the file system cache object we'll be using
        """
        self.cache = cache

    def __post_init__(self):
        """Post init function of the `DataSource` base class.

        This needs to be called at the end of a subclass's constructor.
        """
        self.localisation_ontology = self._load_localisation_ontology()
        self._gene_collection = GeneCollection(self._load_gene_list())

    def fetch_root_file(self, filename: str) -> None:
        """
        Low-level fetching of a root (non-plate) file.

        Fetches the file with the given file name from the data source and
        stores it at the path given by `self.cache.file_path(filename)`.

        This method needs to be implemented in a subclass.

        :param filename: str, name of the file to be fetched
        """
        raise NotImplementedError

    def load_root_file(
        self,
        filename: str,
        return_file_object: bool = True,
    ):
        """
        Load the given root (non-plate) file.

        Tries to load it from the cache first and will revert back to
        `fetch_root_file` if it cannot be found in the cache.

        :param filename: str, name of the file to be loaded
        :param return_file_object: bool (default `True`), whether to return an
            opened file object. Returns the path to the local file if `False.
        """
        logger.debug(f"Loading root file {filename}")
        with self.cache.lock_file_name(filename):
            return self.cache.load_file(
                filename,
                genfile_cb=lambda: self.fetch_root_file(filename),
                return_file_object=return_file_object,
            )

    def fetch_plate_file(self, plate: str, filename: str) -> None:
        """
        Low-level fetching of a plate file.

        Fetches the file with the given file name from the data source and
        stores it at the path given by `self.cache.file_path(filename, plate)`.

        This method needs to be implemented in a subclass.

        :param plate: str, name of the plate the file is included in
        :param filename: str, name of the file to be fetched
        """
        raise NotImplementedError

    def load_plate_file(
        self,
        plate: str,
        filename: str,
        return_file_object: bool = True
    ):
        """
        Load the given plate file.

        Tries to load it from the cache first and will revert back to
        `fetch_root_file` if it cannot be found in the cache.

        :param plate: str, name of the plate the file is included in
        :param filename: str, name of the file to be loaded
        :param return_file_object: bool (default `True`), whether to return an
            opened file object. Returns the path to the local file if `False.
        """
        logger.debug(f"Loading file {filename} from plate {plate}.")
        with self.cache.lock_file_name(f"{plate}_{filename}"):
            return self.cache.load_file(
                filename,
                plate,
                genfile_cb=lambda: self.fetch_plate_file(plate, filename),
                return_file_object=return_file_object,
            )

    def glob_plate_files(
        self,
        plate: str,
        pattern: str,
    ) -> list[str]:
        """
        Find and return a list of files matching the given pattern in the
        given plate.

        This method needs to be implemented in a subclass.

        :param plate: str, the name of the plate the files are located in
        :param pattern: str, a shell-like glob pattern
        :return: a list of files matching the pattern
        """
        raise NotImplementedError

    def _load_gene_list(self) -> dict[str, Gene]:
        logger.debug("Loading gene list.")

        genes: dict[str, Gene] = {}
        with self.load_root_file("localisations.tsv") as localisations:
            locreader = csv.reader(
                localisations,
                delimiter="\t",
            )
            headers = next(locreader)
            for raw_row in locreader:
                row = dict(zip(headers, raw_row))
                N = CellLine.from_data(  # type: ignore[misc]
                    "N",
                    *[row["N " + h] for h in HeadersCellLine],  # type: ignore[arg-type]
                    datasource=self,
                    ontology=self.localisation_ontology
                )
                C = CellLine.from_data(  # type: ignore[misc]
                    "C",
                    *[row["C " + h] for h in HeadersCellLine],  # type: ignore[arg-type]
                    datasource=self,
                    ontology=self.localisation_ontology
                )
                gene = Gene(*[row[h] for h in HeadersGene], N, C)  # type: ignore[call-arg,arg-type]
                C.gene = gene
                N.gene = gene
                genes[gene.id] = gene
        return genes

    def _load_localisation_ontology(self):
        logger.debug("Loading ontology.")
        ontology = Ontology()

        def _populate(
            data: list[dict],
            parent: OntologyEntry | None = None,
            root: bool = False,
        ):
            for raw_entry in data:
                entry = OntologyEntry(
                    raw_entry["name"],
                    raw_entry["synonyms"] if "synonyms" in raw_entry else None,
                    raw_entry["comment"],
                    raw_entry["ident"] if "ident" in raw_entry else None,
                    raw_entry["go"] if "go" in raw_entry else None,
                    raw_entry["examples"] if "examples" in raw_entry else None,
                )

                if entry.name in ontology.entries:
                    raise ValueError(
                        f"Ontology entry '{entry.name}' is defined twice!")
                ontology.entries[entry.name] = entry
                if root:
                    ontology.root_entries.append(entry)

                if parent is not None:
                    entry.set_parent(parent)
                    parent.add_child(entry)

                if "sublocalisation" in raw_entry:
                    _populate(raw_entry["sublocalisation"], parent=entry)

        with self.load_root_file("localisation_ontology.json") as f:
            data = json.load(f)["localisation"]
        _populate(data, root=True)
        return ontology

    @property
    def gene_collection(self):
        """The collection of `Gene` objects present in the data source."""
        return self._gene_collection
