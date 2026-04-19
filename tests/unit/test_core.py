"""Unit tests for gpmap.core.GenotypePhenotypeMap."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from gpmap import GenotypePhenotypeMap
from gpmap.exceptions import SchemaError


def test_constructor_basic(tiny_binary_gpm: GenotypePhenotypeMap) -> None:
    gpm = tiny_binary_gpm
    assert gpm.n == 8
    assert gpm.length == 3
    assert gpm.wildtype == "AAA"
    assert gpm.genotypes[0] == "AAA"
    assert gpm.phenotypes.dtype == np.float64
    assert gpm.stdeviations.dtype == np.float64
    assert gpm.n_replicates.dtype == np.int64


def test_default_stdeviations_are_nan() -> None:
    gpm = GenotypePhenotypeMap(wildtype="AA", genotypes=["AA", "AT"], phenotypes=[0.0, 1.0])
    assert np.isnan(gpm.stdeviations).all()


def test_default_n_replicates_is_one() -> None:
    gpm = GenotypePhenotypeMap(wildtype="AA", genotypes=["AA", "AT"], phenotypes=[0.0, 1.0])
    assert (gpm.n_replicates == 1).all()


def test_infer_mutations_from_genotypes() -> None:
    gpm = GenotypePhenotypeMap(
        wildtype="AAA", genotypes=["AAA", "TAA", "ATA", "AAG"], phenotypes=[0, 1, 2, 3]
    )
    assert set(gpm.mutations[0]) == {"A", "T"}
    assert set(gpm.mutations[1]) == {"A", "T"}
    assert set(gpm.mutations[2]) == {"A", "G"}


def test_binary_packed_shape(tiny_binary_gpm: GenotypePhenotypeMap) -> None:
    B = tiny_binary_gpm.binary_packed
    assert B.shape == (8, 3)
    assert B.dtype == np.uint8


def test_binary_matches_binary_packed(tiny_binary_gpm: GenotypePhenotypeMap) -> None:
    strings = tiny_binary_gpm.binary
    packed = tiny_binary_gpm.binary_packed
    for i, s in enumerate(strings):
        assert s == "".join("1" if b else "0" for b in packed[i])


def test_n_mutations_equals_hamming_to_wildtype(aa_gpm: GenotypePhenotypeMap) -> None:
    # AK -> 0, CK -> 1, AR -> 1, CR -> 2
    assert aa_gpm.n_mutations.tolist() == [0, 1, 1, 2]


def test_phenotype_coercion_warns_on_lossy_cast() -> None:
    with pytest.warns(UserWarning):
        GenotypePhenotypeMap(
            wildtype="A",
            genotypes=["A"],
            phenotypes=np.array([2**53 + 1], dtype=np.int64),
        )


def test_wildtype_not_in_mutations_raises() -> None:
    with pytest.raises(ValueError, match="not present in mutations"):
        GenotypePhenotypeMap(
            wildtype="A",
            genotypes=["A"],
            phenotypes=[0.0],
            mutations={0: ["B", "C"]},
        )


def test_genotype_length_mismatch_raises() -> None:
    with pytest.raises(ValueError, match="same length as wildtype"):
        GenotypePhenotypeMap(wildtype="AA", genotypes=["A"], phenotypes=[0.0])


def test_from_dataframe_roundtrip(tiny_binary_gpm: GenotypePhenotypeMap) -> None:
    df = tiny_binary_gpm.data.copy()
    gpm2 = GenotypePhenotypeMap.from_dataframe(
        df, wildtype=tiny_binary_gpm.wildtype, mutations=tiny_binary_gpm.mutations
    )
    assert gpm2.genotypes.tolist() == tiny_binary_gpm.genotypes.tolist()
    np.testing.assert_array_equal(gpm2.phenotypes, tiny_binary_gpm.phenotypes)


def test_from_dataframe_requires_wildtype() -> None:
    df = pd.DataFrame({"genotypes": ["AA", "AT"], "phenotypes": [0.0, 1.0]})
    with pytest.raises(ValueError, match="wildtype"):
        GenotypePhenotypeMap.from_dataframe(df)


def test_data_property_columns(tiny_binary_gpm: GenotypePhenotypeMap) -> None:
    df = tiny_binary_gpm.data
    assert list(df.columns) == [
        "genotypes",
        "phenotypes",
        "stdeviations",
        "n_replicates",
        "binary",
        "n_mutations",
    ]


def test_subclass_friendly() -> None:
    class MyGPM(GenotypePhenotypeMap):
        pass

    g = MyGPM(wildtype="A", genotypes=["A"], phenotypes=[0.0])
    assert isinstance(g, GenotypePhenotypeMap)
    assert g.n == 1


def test_validate_encoding_table_good(tiny_binary_gpm: GenotypePhenotypeMap) -> None:
    from gpmap import validate_encoding_table

    validate_encoding_table(tiny_binary_gpm.encoding_table)


def test_validate_encoding_table_bad() -> None:
    from gpmap import validate_encoding_table

    bad = pd.DataFrame({"site_index": [0, 1]})
    with pytest.raises(SchemaError):
        validate_encoding_table(bad)
