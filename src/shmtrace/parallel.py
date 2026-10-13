"""Run independent per-library work on several cores.

Results never depend on the number of processes: every task is a pure function
of its arguments (random draws use fixed seeds), and results come back in the
order the tasks were given.

`Params.jobs` sets the number of processes: 1 (the default) runs everything in
this process, 0 uses one fewer than the machine's cores, up to 8. Workers are
started by "spawn", which re-imports the calling script, so a script that asks
for more than one process must guard its work with `if __name__ == "__main__":`.
"""
from __future__ import annotations

import atexit
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context

_POOLS: dict[int, ProcessPoolExecutor] = {}


def resolve(jobs: int | None) -> int:
    if jobs is None or jobs <= 0:
        return max(1, min(8, (os.cpu_count() or 2) - 1))
    return jobs


def _pool(n: int) -> ProcessPoolExecutor:
    # one pool per size, kept for the life of the process: starting workers
    # ("spawn", which is safe with the threaded numerical libraries) costs ~1 s
    if n not in _POOLS:
        _POOLS[n] = ProcessPoolExecutor(max_workers=n, mp_context=get_context("spawn"))
    return _POOLS[n]


@atexit.register
def _shutdown() -> None:
    for pool in _POOLS.values():
        pool.shutdown(wait=False, cancel_futures=True)
    _POOLS.clear()


def pmap(fn, arg_tuples: list[tuple], jobs: int | None = 1, done=None) -> list:
    """[fn(*args) for args in arg_tuples], on up to `jobs` processes. `done(k)` is called
    as each of the k tasks finishes (for progress reports)."""
    n = min(resolve(jobs), len(arg_tuples))
    if n <= 1:
        out = []
        for a in arg_tuples:
            out.append(fn(*a))
            if done:
                done(len(out))
        return out
    futures = [_pool(n).submit(fn, *a) for a in arg_tuples]
    if done:
        for k, _f in enumerate(as_completed(futures), 1):
            done(k)
    return [f.result() for f in futures]
