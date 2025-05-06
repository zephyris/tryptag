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


class WorkList(Sequence):
    def __init__(
        self,
        tryptag: TrypTag,
        work_list: list[CellLine],
    ):
        self.tryptag = tryptag
        self._cell_lines = list(set(work_list))

    def __len__(self):
        return len(self._cell_lines)

    def __getitem__(self, index: int):
        return self._cell_lines[index]

    def _list_analysis_worker(
        self,
        cell_line: CellLine,
        analysis_function: Callable,
    ):
        result = {
            "cell_line": cell_line,
            "result": analysis_function(self.tryptag, cell_line),
        }
        return result

    def process(
        self,
        analysis_function,
        workers: int | None = None,
        multiprocess_mode: Literal["process", "thread"] | None = "process",
    ):
        logger.debug(f"Analysing worklist with {len(self._cell_lines)} "
                     "entries")

        # get number of workers, default to number of cpus
        if workers is None:
            workers = multiprocessing.cpu_count()

        if multiprocess_mode is None:
            # run in a single thread, still use the ThreadPoolExecutor since
            # that's equivalent
            logger.debug("Single process")
            Executor = concurrent.futures.ThreadPoolExecutor
            workers = 1
        elif multiprocess_mode == "process":
            # setup executor as a process pool
            logger.debug(f"Parallel processes with {workers} workers")
            Executor = concurrent.futures.ProcessPoolExecutor
        elif multiprocess_mode == "thread":
            # setup executor as a thread pool
            logger.debug("Parallel threads with", workers, "workers")
            Executor = concurrent.futures.ThreadPoolExecutor
        else:
            raise ValueError(f"Unknown multiprocess_mode '{multiprocess_mode}")

        with Executor(workers) as executor:
            futures = [
                executor.submit(
                    self._list_analysis_worker,
                    cell_line=cell_line,
                    analysis_function=analysis_function,
                ) for cell_line in self]
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
