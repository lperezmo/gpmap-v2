"""Tests for gpmap.encoding."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from gpmap import (
    ENCODING_COLUMNS,
    genotypes_to_binary,
    genotypes_to_binary_packed,
    get_encoding_table,
)
from gpmap.exceptions import UnknownLetterError


def test_encoding_table_columns_and_dtypes() -> None:
    table = get_encoding_table(wildtype="AK", mutations={0: ["A", "C"], 1: ["K", "R"]})
    assert list(table.columns) == list(ENCODING_COLUMNS)
    assert table["site_index"].dtype == "Int64"
    assert table["mutation_index"].dtype == "Int64"
    assert table["binary_index_stop"].dtype == "Int64"


def test_encoding_table_wildtype_rows_have_na_mutation_index() -> None:
    table = get_encoding_table(wildtype="AK", mutations={0: ["A", "C"], 1: ["K", "R"]})
    wt_rows = table[table["mutation_letter"] == table["wildtype_letter"]]
    assert wt_rows["mutation_index"].isna().all()


def test_encoding_table_mutation_index_is_monotonically_increasing() -> None:
    table = get_encoding_table(
        wildtype="AAA", mutations={0: ["A", "T"], 1: ["A", "T"], 2: ["A", "T"]}
    )
    mi = table["mutation_index"].dropna().astype(int).tolist()
    assert mi == [1, 2, 3]


def test_encoding_table_legacy_alias_for_genotype_index() -> None:
    table = get_encoding_table(wildtype="AA", mutations={0: ["A", "T"], 1: ["A", "T"]})
    with pytest.warns(DeprecationWarning, match="renamed to 'site_index'"):
        col = table["genotype_index"]
    pd.testing.assert_series_equal(
        col.reset_index(drop=True),
        table["site_index"].reset_index(drop=True),
        check_names=False,
    )


def test_genotypes_to_binary_packed_matches_encoding() -> None:
    table = get_encoding_table(wildtype="AK", mutations={0: ["A", "C"], 1: ["K", "R"]})
    packed = genotypes_to_binary_packed(["AK", "CK", "AR", "CR"], table)
    assert packed.shape == (4, 2)
    assert packed.tolist() == [[0, 0], [1, 0], [0, 1], [1, 1]]


def test_genotypes_to_binary_strings() -> None:
    table = get_encoding_table(wildtype="AK", mutations={0: ["A", "C"], 1: ["K", "R"]})
    s = genotypes_to_binary(["AK", "CR"], table)
    assert s.tolist() == ["00", "11"]


def test_genotypes_to_binary_unknown_letter_raises() -> None:
    table = get_encoding_table(wildtype="A", mutations={0: ["A", "C"]})
    with pytest.raises(UnknownLetterError):
        genotypes_to_binary_packed(["X"], table)


def test_encoding_frozen_site_has_empty_binary() -> None:
    table = get_encoding_table(wildtype="AX", mutations={0: ["A", "C", "G"], 1: None})
    frozen = table[table["site_index"] == 1]
    assert len(frozen) == 1
    assert frozen.iloc[0]["binary_repr"] == ""
    assert frozen.iloc[0]["binary_index_start"] == frozen.iloc[0]["binary_index_stop"]


def test_encoding_larger_alphabet() -> None:
    table = get_encoding_table(wildtype="A", mutations={0: ["A", "C", "G", "T"]})
    non_wt = table[table["mutation_letter"] != table["wildtype_letter"]]
    assert len(non_wt) == 3
    reprs = sorted(non_wt["binary_repr"].tolist())
    assert reprs == ["001", "010", "100"]


def test_packed_sum_equals_mutations_from_wt() -> None:
    table = get_encoding_table(
        wildtype="AAA", mutations={0: ["A", "T"], 1: ["A", "T"], 2: ["A", "T"]}
    )
    packed = genotypes_to_binary_packed(["AAA", "TAT", "TTT"], table)
    assert packed.sum(axis=1).tolist() == [0, 2, 3]


def test_genotypes_to_binary_packed_parallel_consistency() -> None:
    # Exercise the rayon-parallel hot path with enough rows to cross chunks.
    n = 5000
    table = get_encoding_table(
        wildtype="AAAA", mutations={0: ["A", "T"], 1: ["A", "T"], 2: ["A", "T"], 3: ["A", "T"]}
    )
    rng = np.random.default_rng(0)
    rows = ["".join("A" if b == 0 else "T" for b in rng.integers(0, 2, 4)) for _ in range(n)]
    packed = genotypes_to_binary_packed(rows, table)
    # Hamming weight should equal count of 'T' in each row.
    expected = np.array([g.count("T") for g in rows])
    assert (packed.sum(axis=1) == expected).all()
