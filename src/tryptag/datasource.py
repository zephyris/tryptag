from __future__ import annotations
from collections.abc import Mapping
import csv
from dataclasses import InitVar, dataclass
import enum
from functools import cached_property
import json
import logging
import pathlib
from typing import Literal

from .cache import Cache, FileTypes, FILE_TYPE_ENDINGS, FILE_PATTERN

from annotations import (
    AnnotationCollection,
    Ontology,
    OntologyEntry,
    OntologyAnnotationCollection,
)

logger = logging.getLogger(__name__)

HEADERS_GENE = [
    "Gene ID",
    "Gene aliases",
]
HEADERS_CELL_LINE = [
    "status",
    "plate and well",
    "primer F",
    "primer R",
    "localisation",
    "fainter than parental",
    "classified as faint",
]

TERMINI = ["C", "N"]


class CellLineStatus(enum.Enum):
    NOT_ATTEMPTED = 0
    ATTEMPTED = 1
    GENERATED = 2


CELL_LINE_STATUS_MAP = {
    "not attempted": CellLineStatus.NOT_ATTEMPTED,
    "tagging attempted": CellLineStatus.ATTEMPTED,
    "cell line generated": CellLineStatus.GENERATED,
}


class Cell:
    def __init__(
        self,
        field: Field,
        index: int,
        wand: tuple[int, int],
        centre: tuple[float, float],
        extent: tuple[float, float],
        angle: float
    ):
        self.field = field
        self.index = index
        self.wand = wand
        self.centre = centre
        self.extent = extent
        self.angle = angle

    @staticmethod
    def from_line(field: Field, line: str):
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


class FieldDoesNotExistError(Exception):
    def __init__(
        self,
        cell_line: CellLine,
        index: int,
    ):
        self.cell_line = cell_line
        self.index = index


class Field:
    def __init__(
        self,
        cell_line: CellLine,
        index: int,
        datasource: DataSource,
    ):
        self.cell_line = cell_line
        self.index = index
        self.datasource = datasource

        if not self.exists():
            raise FieldDoesNotExistError(self.cell_line, self.index)

        cell_file = datasource.load_plate_file(
            cell_line.plate,
            self.filename(FileTypes.CELL_ROIS)
        )
        # Skip header line
        next(cell_file)
        cell_list = [Cell.from_line(self, line) for line in cell_file]
        self.cells = {cell.index: cell for cell in cell_list}

    def filename(self, file_type: FileTypes):
        return (
            f"{self.cell_line.filename_stem()}{self.index+1}"
            f"{FILE_TYPE_ENDINGS[file_type]}"
        )

    def exists(self):
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


@dataclass
class CellLine:
    terminus: Literal["N"] | Literal["C"]
    status_raw: InitVar[str]
    plate_and_well_raw: InitVar[str]
    forward_primer: str
    reverse_primer: str
    localisation_raw: InitVar[str]
    fainter_than_parental: bool
    classified_faint: bool
    datasource: DataSource
    ontology: InitVar[Ontology]
    status: CellLineStatus | None = None
    plate: str | None = None
    well: str | None = None
    gene_id: str | None = None
    life_stage: str | None = None
    localisation: AnnotationCollection | None = None

    def __post_init__(
        self,
        status_raw: str,
        plate_and_well_raw: str,
        localisation_raw: str,
        ontology: Ontology,
    ):
        self.status = CELL_LINE_STATUS_MAP[status_raw]
        if self.status != CellLineStatus.NOT_ATTEMPTED:
            self.plate, self.well = plate_and_well_raw.split(" ")

        self.localisation = OntologyAnnotationCollection(
            localisation_raw,
            ontology
        )

    def __hash__(self):
        return hash(
            (self.datasource, self.gene_id, self.terminus, self.life_stage))

    def __repr__(self):
        return (
            f"gene_id = {self.gene_id} terminus = {self.terminus} "
            f"life_stage = {self.life_stage}"
        )

    def filename_stem(self):
        return f"{self.gene_id}_4_{self.terminus}_"

    @cached_property
    def fields(self):
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


@dataclass
class Gene(Mapping):
    id: str
    aliases: str
    N: CellLine
    C: CellLine

    def __getitem__(self, terminus: Literal["N"] | Literal["C"]):
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


class DataSource(Mapping):
    def __init__(self, cache: Cache):
        self.cache = cache

    def __post_init__(self):
        self.localisation_ontology = self._load_localisation_ontology()
        self.genes = self._load_gene_list()

    def fetch_root_file(self, filename: str):
        raise NotImplementedError

    def load_root_file(
        self,
        filename: str,
        return_file_object: bool = True,
    ):
        with self.cache.lock_file_name(filename):
            return self.cache.load_file(
                filename,
                genfile_cb=lambda: self.fetch_root_file(filename),
                return_file_object=return_file_object,
            )

    def fetch_plate_file(self, plate: str, filename: str):
        raise NotImplementedError

    def load_plate_file(
        self,
        plate: str,
        filename: str,
        return_file_object: bool = True
    ):
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
    ):
        raise NotImplementedError

    def _load_gene_list(self):
        logger.debug("Loading gene list.")

        genes: dict[str: Gene] = {}
        with self.load_root_file("localisations.tsv") as localisations:
            locreader = csv.reader(
                localisations,
                delimiter="\t",
            )
            headers = next(locreader)
            for row in locreader:
                row = dict(zip(headers, row))
                N = CellLine(
                    "N",
                    *[row["N " + h] for h in HEADERS_CELL_LINE],
                    datasource=self,
                    ontology=self.localisation_ontology
                )
                C = CellLine(
                    "C",
                    *[row["C " + h] for h in HEADERS_CELL_LINE],
                    datasource=self,
                    ontology=self.localisation_ontology
                )
                gene = Gene(*[row[h] for h in HEADERS_GENE], N, C)
                C.gene_id = gene.id
                N.gene_id = gene.id
                genes[gene.id] = gene
        return genes

    def _load_localisation_ontology(self):
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

    def __getitem__(self, geneid: str):
        return self.genes[geneid]

    def __iter__(self):
        return iter(self.genes)

    def __len__(self):
        return len(self.genes)
