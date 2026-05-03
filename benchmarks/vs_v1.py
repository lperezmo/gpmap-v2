"""Benchmark gpmap-v2 construction and binary_packed access.

To reproduce the v1 vs v2 comparison from the README, run in separate environments:

    # v1 environment
    uv venv .venv-v1 && uv pip install --python .venv-v1/Scripts/python.exe gpmap==0.7.0
    .venv-v1/Scripts/python benchmarks/vs_v1.py

    # v2 environment
    uv run python benchmarks/vs_v1.py
"""
from __future__ import annotations

import itertools
import json
import timeit
from pathlib import Path

import numpy as np

SIZES = [8, 10, 12, 14, 16]
N_CALLS = {8: 10, 10: 5, 12: 3, 14: 1, 16: 1}
N_REPEATS = 5


def make_data(L: int) -> tuple[str, list[str], np.ndarray]:
    genotypes = ["".join(g) for g in itertools.product("AT", repeat=L)]
    phenotypes = np.random.default_rng(42).normal(size=len(genotypes))
    return "A" * L, genotypes, phenotypes


def best_ms(fn, n_calls: int) -> float:
    return min(timeit.repeat(fn, number=n_calls, repeat=N_REPEATS)) / n_calls * 1000


from gpmap import GenotypePhenotypeMap

try:
    from importlib.metadata import version as _pkg_ver
    _ver = _pkg_ver("gpmap-v2")
except Exception:
    try:
        _ver = _pkg_ver("gpmap")
    except Exception:
        _ver = "unknown"

results: dict = {"version": _ver, "results": {}}
print(f"gpmap {_ver}")

print("construct")
construct_ms: dict[str, float] = {}
for L in SIZES:
    wt, gts, pts = make_data(L)
    t = best_ms(lambda: GenotypePhenotypeMap(wildtype=wt, genotypes=gts, phenotypes=pts), N_CALLS[L])
    construct_ms[f"L{L}"] = round(t, 4)
    print(f"  L={L:2d} ({2**L:6d} genotypes): {t:.3f} ms")
results["results"]["construct_ms"] = construct_ms

print("binary / binary_packed access")
binary_ms: dict[str, float] = {}
for L in SIZES:
    wt, gts, pts = make_data(L)
    gpm = GenotypePhenotypeMap(wildtype=wt, genotypes=gts, phenotypes=pts)
    try:
        t = best_ms(lambda: gpm.binary_packed, N_CALLS[L])
        label = "binary_packed"
    except AttributeError:
        t = best_ms(lambda: list(gpm.binary), N_CALLS[L])
        label = "binary"
    binary_ms[f"L{L}"] = round(t, 4)
    print(f"  L={L:2d} ({2**L:6d} genotypes) [{label}]: {t:.4f} ms")
results["results"]["binary_ms"] = binary_ms

out = Path(__file__).parent / f"results_{_ver.replace('.', '_')}.json"
out.write_text(json.dumps(results, indent=2))
print(f"Saved {out}")
