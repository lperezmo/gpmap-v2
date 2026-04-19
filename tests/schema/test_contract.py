"""Contract tests: lock the surface that epistasis-v2 imports.

If these tests break, epistasis-v2 breaks. Only change them via a coordinated
breaking change bumping gpmap-v2's major version.
"""

from __future__ import annotations

import inspect

import gpmap
import numpy as np
from gpmap import ENCODING_COLUMNS


def test_public_api_surface_stable() -> None:
    required = {
        "GenotypePhenotypeMap",
        "genotypes_to_binary",
        "genotypes_to_binary_packed",
        "get_encoding_table",
        "upper_transform",
        "lower_transform",
        "StandardDeviationMap",
        "StandardErrorMap",
        "SpaceTooLargeError",
        "SchemaError",
        "UnknownLetterError",
        "read_csv",
        "read_json",
        "read_pickle",
        "read_excel",
        "to_csv",
        "to_json",
        "to_pickle",
        "to_excel",
    }
    missing = required - set(gpmap.__all__)
    assert not missing, f"missing public names: {missing}"


def test_encoding_table_columns_locked() -> None:
    assert ENCODING_COLUMNS == (
        "site_index",
        "site_label",
        "wildtype_letter",
        "mutation_letter",
        "mutation_index",
        "binary_repr",
        "binary_index_start",
        "binary_index_stop",
    )


def test_encoding_table_epistasis_consumer_pattern() -> None:
    """epistasis-v2 does:
    encoding_table[["mutation_index","site_index"]].dropna().astype(int)
    then groupby("site_index"). Verify both work.
    """
    gpm = gpmap.GenotypePhenotypeMap(
        wildtype="AA",
        genotypes=["AA", "TA", "AT"],
        phenotypes=[0, 1, 2],
        mutations={0: ["A", "T"], 1: ["A", "T"]},
    )
    table = gpm.encoding_table
    subset = table[["mutation_index", "site_index"]].dropna().astype(int)
    groups = subset.groupby("site_index").groups
    assert set(groups.keys()) == {0, 1}


def test_encoding_table_label_mapper_pattern() -> None:
    """epistasis-v2 does:
    t = gpm.encoding_table.dropna()
    for _, row in t.iterrows(): row.wildtype_letter, row.site_label, row.mutation_letter, row.mutation_index
    """
    gpm = gpmap.GenotypePhenotypeMap(
        wildtype="AA",
        genotypes=["AA", "TA"],
        phenotypes=[0, 1],
        mutations={0: ["A", "T"], 1: ["A", "T"]},
    )
    t = gpm.encoding_table.dropna()
    for _, row in t.iterrows():
        # All four attributes must be readable.
        _ = row["wildtype_letter"]
        _ = row["site_label"]
        _ = row["mutation_letter"]
        _ = row["mutation_index"]


def test_gpm_required_attributes_present() -> None:
    gpm = gpmap.GenotypePhenotypeMap(wildtype="AA", genotypes=["AA"], phenotypes=[0.0])
    for attr in (
        "wildtype",
        "genotypes",
        "phenotypes",
        "stdeviations",
        "n_replicates",
        "encoding_table",
        "binary",
        "binary_packed",
        "data",
        "mutations",
        "site_labels",
    ):
        assert hasattr(gpm, attr), f"missing attribute: {attr}"


def test_from_dataframe_exists_and_is_classmethod() -> None:
    fn = gpmap.GenotypePhenotypeMap.from_dataframe
    assert callable(fn)
    # Verify it's usable as a classmethod (bound to the class).
    assert inspect.ismethod(gpmap.GenotypePhenotypeMap.from_dataframe)


def test_phenotype_dtype_is_float64() -> None:
    gpm = gpmap.GenotypePhenotypeMap(wildtype="A", genotypes=["A"], phenotypes=[1])
    assert gpm.phenotypes.dtype == np.float64


def test_stdeviations_dtype_is_float64() -> None:
    gpm = gpmap.GenotypePhenotypeMap(
        wildtype="A", genotypes=["A"], phenotypes=[0.0], stdeviations=[0.5]
    )
    assert gpm.stdeviations.dtype == np.float64


def test_binary_packed_is_uint8_2d() -> None:
    gpm = gpmap.GenotypePhenotypeMap(wildtype="A", genotypes=["A"], phenotypes=[0.0])
    # Frozen site, no bits. Still a 2D uint8 array.
    assert gpm.binary_packed.dtype == np.uint8
    assert gpm.binary_packed.ndim == 2


def test_binary_is_object_array_of_strings() -> None:
    gpm = gpmap.GenotypePhenotypeMap(
        wildtype="AA",
        genotypes=["AA", "TA"],
        phenotypes=[0.0, 1.0],
        mutations={0: ["A", "T"], 1: ["A", "T"]},
    )
    assert gpm.binary.dtype == object
    assert all(isinstance(s, str) for s in gpm.binary)


def test_gpm_is_subclassable() -> None:
    class MySim(gpmap.GenotypePhenotypeMap):
        pass

    g = MySim(wildtype="A", genotypes=["A"], phenotypes=[0.0])
    assert isinstance(g, gpmap.GenotypePhenotypeMap)


def test_genotypes_to_binary_signature() -> None:
    sig = inspect.signature(gpmap.genotypes_to_binary)
    # Must accept (genotypes, encoding_table). Order matters for back-compat.
    params = list(sig.parameters.keys())
    assert params[:2] == ["genotypes", "encoding_table"]


def test_upper_lower_transform_are_distinct_callables() -> None:
    assert gpmap.upper_transform is not gpmap.lower_transform
    y = [1.0, 2.0]
    u = [1.5, 2.5]
    lo = [0.8, 1.9]
    # Different inputs, different outputs.
    np.testing.assert_allclose(gpmap.upper_transform(y, u), [0.5, 0.5])
    np.testing.assert_allclose(gpmap.lower_transform(y, lo), [0.2, 0.1])
