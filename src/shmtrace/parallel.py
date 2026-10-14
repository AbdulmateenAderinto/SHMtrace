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
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context

_POOLS: dict[int, ProcessPoolExecutor] = {}


def resolve(jobs: int | None) -> int:
    if jobs is None or jobs <= 0:
        return max(1, min(8, (os.cpu_count() or 2) - 1))
    return jobs


def _pool(n: int) -> ProcessPoolExecutor:
    # one pool per size, kept for the life of the process: starting workers
    # ("spawn", which is safe with the threaded numerical libraries) costs ~1 s.
    # Callers ask for resolve(jobs) workers, not one per task, so a run uses one pool.
    if n not in _POOLS:
        _POOLS[n] = ProcessPoolExecutor(max_workers=n, mp_context=get_context("spawn"),
                                        initializer=_watch_parent, initargs=(os.getpid(),))
        _exit_on_sigterm()
    return _POOLS[n]


def _watch_parent(parent: int) -> None:
    # In each worker: leave as soon as the parent is gone, however it ended (SIGKILL, a
    # force-quit app, or SIGTERM in a program whose pool lives on a background thread,
    # where _exit_on_sigterm cannot install its handler). An orphan is re-parented.
    import threading
    import time

    def watch():
        while os.getppid() == parent:
            time.sleep(1)
        os._exit(0)
    threading.Thread(target=watch, daemon=True).start()


def _exit_on_sigterm() -> None:
    # Idle workers outlive a parent killed by SIGTERM (as a desktop app is when the system
    # quits it), each holding its memory. Turn SIGTERM into a normal exit, so _shutdown runs,
    # unless the program has a handler of its own (uvicorn run in the main thread does).
    import signal
    import threading
    if threading.current_thread() is threading.main_thread() and \
            signal.getsignal(signal.SIGTERM) is signal.SIG_DFL:
        signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))


@atexit.register
def _shutdown() -> None:
    for pool in _POOLS.values():
        procs = list((getattr(pool, "_processes", None) or {}).values())
        pool.shutdown(wait=False, cancel_futures=True)
        for proc in procs:          # a worker busy with a task would otherwise finish it first
            if proc.is_alive():
                proc.terminate()
    _POOLS.clear()


def pmap(fn, arg_tuples: list[tuple], jobs: int | None = 1, done=None) -> list:
    """[fn(*args) for args in arg_tuples], on up to `jobs` processes. `done(k)` is called
    as each of the k tasks finishes (for progress reports)."""
    n = resolve(jobs)
    if n <= 1 or len(arg_tuples) <= 1:
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
