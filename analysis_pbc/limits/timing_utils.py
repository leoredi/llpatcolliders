from __future__ import annotations

from contextlib import contextmanager
from time import perf_counter


@contextmanager
def _time_block(timing: dict | None, key: str):
    if timing is None:
        yield
        return
    start = perf_counter()
    try:
        yield
    finally:
        timing[key] = timing.get(key, 0.0) + (perf_counter() - start)
