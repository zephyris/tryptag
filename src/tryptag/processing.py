from __future__ import annotations
from collections.abc import Callable, Sequence
import concurrent.futures
import logging
import multiprocessing
from typing import TYPE_CHECKING, Literal

from tqdm import tqdm

from .datasource import CellLine

if TYPE_CHECKING:
    from .tryptag import TrypTag

logger = logging.getLogger("tryptag.processing")


_tryptag_object: TrypTag | None = None
def _list_analysis_worker_processes_init(tryptag: TrypTag):
    global _tryptag_object
    _tryptag_object = tryptag


def _list_analysis_worker_processes(
    cell_line: CellLine,
    analysis_function: Callable,
):
    result = {
        # Give a shallow copy of the CellLine object only!
        # If the analysis function populates the CellLine, Field, Cell
        # objects, etc., the multiprocessing pickles the whole object
        # hierarchy which leads to very slow processing and a giant
        # memory leak in the parent process.
        "cell_line": CellLine(cell_line.gene_id, cell_line.terminus),
        "result": analysis_function(
            _tryptag_object,
            _tryptag_object.gene_list[cell_line.gene_id][cell_line.terminus],
        ),
    }
    return result


class WorkList(Sequence[CellLine]):
    """
    Class holding a list of `CellLine`s.

    This class holds a list of CellLine objects and implements facilities for
    querying / filtering this list and for running image and other analyses.
    """
    def __init__(
        self,
        tryptag: TrypTag,
        cell_line_list: list[CellLine],
    ):
        """
        Initialise a `WorkList` instance.

        :param tryptag: A reference to the core `TrypTag` object (to run
            analysis code).
        :param cell_line_list: A list of `CellLine` objects to include.
        """
        self.tryptag = tryptag
        self._cell_lines = list(set(cell_line_list))

    def __len__(self):
        return len(self._cell_lines)

    def __getitem__(self, index: int):
        return self._cell_lines[index]

    def __repr__(self):
        return f"WorkList([{','.join([str(c) for c in self])}])"

    def union(self, other: WorkList):
        return WorkList(
            self.tryptag,
            set(self).union(set(other))
        )

    def intersection(self, other: WorkList):
        return WorkList(
            self.tryptag,
            set(self).intersection(set(other))
        )

    def difference(self, other: WorkList):
        return WorkList(
            self.tryptag,
            set(self).difference(set(other))
        )

    def filter(
        self,
        filter_function: Callable[[CellLine], bool]
    ):
        """
        Advanced `CellLine` filtering functionality.

        Runs `filter_function` on each `CellLine` object and returns a
        `WorkList` of `CellLines` for which it evaluated as `True`.
        """
        return WorkList(
            self.tryptag,
            [
                cell_line
                for cell_line in self
                if filter_function(cell_line)
            ]
        )

    def gene_id_search(
        self,
        gene_id_list: list[str],
    ):
        """
        Takes a list of Gene IDs and returns a new `WorkList` containing all
        `CellLine`s in this current `WorkList` whose `gene_id` was in the
        given list.
        """
        return self.filter(
            lambda cell_line: cell_line.gene_id in gene_id_list
        )

    def localisation_search(
            self,
            query_term: str,
            match_subterms: bool = True,
            exclude_modifiers: list[str] = ["weak", "<10%"],
            required_modifiers: list[str] | None = None
    ):
        """
        Filters for localisation terms and modifiers.

        This returns a new `WorkList` containing all `CellLine`s in this
        current `WorkList` filtered by localisation annotation. A `CellLine`
        is only included if at least one annotation term matches the given
        `query_term` while containing all modifiers in `required_modifiers`
        and none of the modifiers in `exclude_modifiers`. If `match_subterms`
        is given, the ontology parents of all annotation terms of each
        `CellLine` are matched recursively until either a match is found or the
        root of the ontology is reached (for example if a `CellLine` has the
        annotation "nucleolus" it will be included when we search for
        "nucleus" if `match_subterms` is `True`).

        :param query_term: Search query annotation term from the localisation
            ontology.
        :param match_subterms: Whether to also match
            child/subterm/substructures of `query_term`, default `True`.
        :param exclude_modifiers: List of modifier terms none of which can be
            matched, default `"weak"` and `"<10%"`.
        :param required_modifiers: List of modifier terms all of which must be
            matched, default `None`.
        :return: A new `WorkList` containing the matching `CellLine`s.
        """
        if required_modifiers is None:
            set_of_required_mods = None
        else:
            set_of_required_mods = set(required_modifiers)
        if exclude_modifiers is None:
            set_of_excluded_mods = None
        else:
            set_of_excluded_mods = set(exclude_modifiers)

        return self.filter(
            lambda cell_line: cell_line.localisation.match(
                query_term,
                require_modifiers=set_of_required_mods,
                exclude_modifiers=set_of_excluded_mods,
                recursive=match_subterms,
            )
        )

    def _list_analysis_worker_thread(
        self,
        cell_line: CellLine,
        analysis_function: Callable,
    ):
        result = {
            "cell_line": cell_line,
            "result": analysis_function(self.tryptag, cell_line),
        }
        return result

    def analysis(
        self,
        analysis_function,
        workers: int | None = None,
        multiprocess_mode: Literal["process", "thread"] | None = "process",
    ):
        """
        Handles (potentially multi-process or multi-threaded) analysis of
        `CellLine`s included in this `WorkList`.

        :param analysis_function: Function to use for analysis.
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
        logger.debug("Running analysis of work list with "
                     f"{len(self._cell_lines)} entries")

        # get number of workers, default to number of cpus
        if workers is None:
            workers = multiprocessing.cpu_count()

        if multiprocess_mode is None:
            # run in a single thread, still use the ThreadPoolExecutor since
            # that's equivalent
            logger.debug("Single process")
            executor = concurrent.futures.ThreadPoolExecutor(1)
        elif multiprocess_mode == "process":
            # setup executor as a process pool
            logger.debug(f"Parallel processes with {workers} workers")
            executor = concurrent.futures.ProcessPoolExecutor(
                workers,
                initializer=_list_analysis_worker_processes_init,
                initargs=(self.tryptag,)
            )
        elif multiprocess_mode == "thread":
            # setup executor as a thread pool
            logger.debug("Parallel threads with", workers, "workers")
            executor = concurrent.futures.ThreadPoolExecutor(workers)
        else:
            raise ValueError(f"Unknown multiprocess_mode '{multiprocess_mode}")

        with executor:
            if multiprocess_mode == "process":
                futures = [
                    executor.submit(
                        _list_analysis_worker_processes,
                        cell_line=CellLine(cell_line.gene_id, cell_line.terminus),
                        analysis_function=analysis_function,
                    ) for cell_line in self
                ]
            else:
                futures = [
                    executor.submit(
                        self._list_analysis_worker_thread,
                        cell_line=cell_line,
                        analysis_function=analysis_function,
                    ) for cell_line in self
                ]
            results = []
            total = len(futures)
            for nr, future in enumerate(tqdm(
                concurrent.futures.as_completed(futures),
                total=total,
                smoothing=0,
                disable=logger.getEffectiveLevel() != logging.INFO,
            )):
                logger.debug(f"{nr + 1} / {total} done.")
                results.append(future.result())
        return results
