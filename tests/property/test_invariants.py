"""Hypothesis-based property tests for schema invariants."""

from __future__ import annotations

import string

import numpy as np
from gpmap import (
    GenotypePhenotypeMap,
    genotypes_to_binary_packed,
    get_encoding_table,
)
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st

ALPHA = st.sampled_from(list(string.ascii_uppercase))


@st.composite
def small_mutations(draw) -> tuple[str, dict[int, list[str]]]:
    L = draw(st.integers(min_value=1, max_value=4))
    per_site: list[list[str]] = []
    for _ in range(L):
        size = draw(st.integers(min_value=2, max_value=4))
        letters = draw(st.lists(ALPHA, min_size=size, max_size=size, unique=True))
        per_site.append(sorted(letters))
    wildtype = "".join(draw(st.sampled_from(alpha)) for alpha in per_site)
    mutations = {i: alpha for i, alpha in enumerate(per_site)}
    return wildtype, mutations


@given(spec=small_mutations())
@settings(
    max_examples=30,
    suppress_health_check=[HealthCheck.too_slow],
    deadline=None,
)
def test_encoding_table_row_count(spec: tuple[str, dict[int, list[str]]]) -> None:
    wildtype, mutations = spec
    table = get_encoding_table(wildtype=wildtype, mutations=mutations)
    # One row per (site, letter-in-alphabet) combination.
    expected_rows = sum(len(mutations[i]) for i in range(len(wildtype)))
    assert len(table) == expected_rows


@given(spec=small_mutations())
@settings(
    max_examples=30,
    suppress_health_check=[HealthCheck.too_slow],
    deadline=None,
)
def test_binary_packed_sum_is_hamming_to_wt(
    spec: tuple[str, dict[int, list[str]]],
) -> None:
    wildtype, mutations = spec
    # Enumerate a subset of genotypes (up to 20) from the full space.
    from gpmap.enumerate import enumerate_genotypes_str

    all_gs = enumerate_genotypes_str(wildtype, mutations)
    assume(len(all_gs) <= 1024)
    sample = all_gs[: min(len(all_gs), 20)]
    table = get_encoding_table(wildtype=wildtype, mutations=mutations)
    packed = genotypes_to_binary_packed(sample, table)
    expected = np.array([sum(1 for ci, cg in zip(g, wildtype) if ci != cg) for g in sample])
    # packed sum should equal the hamming distance in sequence space,
    # since the unary encoding is one bit per non-WT letter chosen.
    assert (packed.sum(axis=1) == expected).all()


@given(spec=small_mutations())
@settings(
    max_examples=20,
    suppress_health_check=[HealthCheck.too_slow],
    deadline=None,
)
def test_gpm_roundtrip_via_from_dataframe(
    spec: tuple[str, dict[int, list[str]]],
) -> None:
    wildtype, mutations = spec
    from gpmap.enumerate import enumerate_genotypes_str

    gs = enumerate_genotypes_str(wildtype, mutations)
    assume(len(gs) <= 256)
    phenotypes = np.arange(len(gs), dtype=np.float64)
    gpm = GenotypePhenotypeMap(
        wildtype=wildtype,
        genotypes=gs,
        phenotypes=phenotypes,
        mutations=mutations,
    )
    gpm2 = GenotypePhenotypeMap.from_dataframe(gpm.data, wildtype=wildtype, mutations=mutations)
    assert gpm2.genotypes.tolist() == gs
    np.testing.assert_array_equal(gpm2.phenotypes, phenotypes)
