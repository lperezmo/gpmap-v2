"""I/O round-trip tests: CSV, JSON, pickle, Excel."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from gpmap import (
    GenotypePhenotypeMap,
    read_csv,
    read_excel,
    read_json,
    read_pickle,
    to_csv,
    to_excel,
    to_json,
    to_pickle,
)


def _assert_gpm_equal(a: GenotypePhenotypeMap, b: GenotypePhenotypeMap) -> None:
    assert a.wildtype == b.wildtype
    assert a.genotypes.tolist() == b.genotypes.tolist()
    np.testing.assert_array_equal(a.phenotypes, b.phenotypes)
    # stdeviations may have NaN; compare via nan-aware equality.
    np.testing.assert_array_equal(
        np.asarray(a.stdeviations, dtype=np.float64),
        np.asarray(b.stdeviations, dtype=np.float64),
    )
    np.testing.assert_array_equal(a.n_replicates, b.n_replicates)
    assert a.mutations == b.mutations
    assert a.site_labels == b.site_labels


def test_json_roundtrip(tiny_binary_gpm: GenotypePhenotypeMap, tmp_path: Path) -> None:
    p = tmp_path / "gpm.json"
    to_json(tiny_binary_gpm, p)
    gpm2 = read_json(p)
    _assert_gpm_equal(tiny_binary_gpm, gpm2)


def test_csv_roundtrip(tiny_binary_gpm: GenotypePhenotypeMap, tmp_path: Path) -> None:
    p = tmp_path / "gpm.csv"
    to_csv(tiny_binary_gpm, p)
    assert (tmp_path / "gpm.csv.meta.json").exists()
    gpm2 = read_csv(p)
    _assert_gpm_equal(tiny_binary_gpm, gpm2)


def test_pickle_roundtrip(tiny_binary_gpm: GenotypePhenotypeMap, tmp_path: Path) -> None:
    p = tmp_path / "gpm.pkl"
    to_pickle(tiny_binary_gpm, p)
    gpm2 = read_pickle(p)
    _assert_gpm_equal(tiny_binary_gpm, gpm2)


def test_excel_roundtrip(tiny_binary_gpm: GenotypePhenotypeMap, tmp_path: Path) -> None:
    p = tmp_path / "gpm.xlsx"
    to_excel(tiny_binary_gpm, p)
    gpm2 = read_excel(p)
    _assert_gpm_equal(tiny_binary_gpm, gpm2)


def test_pickle_type_guard(tmp_path: Path) -> None:
    import pickle

    p = tmp_path / "bad.pkl"
    with open(p, "wb") as f:
        pickle.dump({"not": "a gpm"}, f)
    with pytest.raises(TypeError, match="GenotypePhenotypeMap"):
        read_pickle(p)


def test_legacy_json_without_schema_version(tmp_path: Path) -> None:
    import json

    p = tmp_path / "legacy.json"
    p.write_text(
        json.dumps(
            {
                "wildtype": "A",
                "mutations": {"0": ["A", "T"]},
                "site_labels": ["0"],
                "genotypes": ["A", "T"],
                "phenotypes": [0.0, 1.0],
                "stdeviations": [float("nan"), float("nan")],
                "n_replicates": [1, 1],
            }
        ),
        encoding="utf-8",
    )
    with pytest.warns(UserWarning, match="schema_version"):
        gpm = read_json(p)
    assert gpm.n == 2
