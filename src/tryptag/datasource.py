from __future__ import annotations
from collections.abc import Mapping
import csv
from dataclasses import dataclass
import enum
from functools import cached_property
import json
import logging
import pathlib
from typing import Literal
import warnings

from .cache import Cache, FileTypes, FILE_TYPE_ENDINGS, FILE_PATTERN

from annotations import (
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
        with cell_file:
            # Skip header line
            next(cell_file)
            cell_list = [Cell.from_line(self, line) for line in cell_file]
            self.cells = {cell.index: cell for cell in cell_list}

    def filename(self, file_type: FileTypes):
        return (
            f"{self.cell_line.filename_stem()}{self.index+1}"
            f"{FILE_TYPE_ENDINGS[file_type]}"
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
        self._initialised = False
        self._gene_id = gene_id
        terminus = terminus.upper()
        if terminus not in ["N", "C"]:
            raise ValueError("terminus needs to be either 'N' or 'C'")
        self.terminus = terminus

    @property
    def initialised(self):
        return self._initialised
    
    @property
    def gene_id(self):
        if self.gene is not None:
            return self.gene.id
        return self._gene_id

    @staticmethod
    def from_data(
        terminus: Literal["N", "C"],
        status: str,
        plate_and_well: str,
        forward_primer: str,
        reverse_primer: str,
        localisation: str,
        fainter_than_parental: bool,
        classified_faint: bool,
        datasource: DataSource,
        ontology: Ontology,
    ):
        self = CellLine("", terminus)
        self.terminus = terminus
        self.status = CELL_LINE_STATUS_MAP[status]
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

    def __repr__(self):
        return (
            f"gene_id = {self.gene.id if self.gene is not None else self.gene_id} "
            f"terminus = {self.terminus} life_stage = {self.life_stage}"
        )

    def filename_stem(self):
        return f"{self.gene.id}_4_{self.terminus}_"

    @cached_property
    def fields(self) -> dict[int, Field]:
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
    def __init__(self, genes: dict[str, Gene]):
        self.genes = genes

    def __getitem__(self, geneid: str):
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
    def __init__(self, cache: Cache):
        self.cache = cache

    def __post_init__(self):
        self.localisation_ontology = self._load_localisation_ontology()
        self._gene_collection = GeneCollection(self._load_gene_list())

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
                    *[row["N " + h] for h in HEADERS_CELL_LINE],  # type: ignore[arg-type]
                    datasource=self,
                    ontology=self.localisation_ontology
                )
                C = CellLine.from_data(  # type: ignore[misc]
                    "C",
                    *[row["C " + h] for h in HEADERS_CELL_LINE],  # type: ignore[arg-type]
                    datasource=self,
                    ontology=self.localisation_ontology
                )
                gene = Gene(*[row[h] for h in HEADERS_GENE], N, C)  # type: ignore[call-arg,arg-type]
                C.gene = gene
                N.gene = gene
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
        return self._gene_collection
