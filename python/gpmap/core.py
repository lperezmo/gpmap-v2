"""GenotypePhenotypeMap: the core container. Contract lives in SCHEMA.md."""

from __future__ import annotations

import warnings
from collections.abc import Mapping
from typing import Any

import numpy as np
import pandas as pd

from .encoding import (
    ENCODING_COLUMNS,
    genotypes_to_binary_packed,
    get_encoding_table,
    n_bits_from_encoding_table,
)
from .errors import StandardDeviationMap, StandardErrorMap
from .exceptions import SchemaError


def _infer_mutations(wildtype: str, genotypes: list[str]) -> dict[int, list[str] | None]:
    """Infer per-site alphabets from observed genotypes; WT is always included."""
    out: dict[int, list[str] | None] = {}
    L = len(wildtype)
    for i in range(L):
        seen: set[str] = {wildtype[i]}
        for g in genotypes:
            seen.add(g[i])
        letters = sorted(seen)
        # Frozen site: only WT observed and alphabet size is 1
        out[i] = letters
    return out


def _coerce_phenotypes(values: Any) -> np.ndarray:
    arr = np.asarray(values)
    coerced = np.asarray(values, dtype=np.float64)
    if arr.dtype != coerced.dtype:
        try:
            # Round-trip back to the original dtype. If we lose precision
            # (e.g. int64 values beyond 2**53 cannot be stored exactly as float64),
            # the two will disagree.
            exact = np.array_equal(coerced.astype(arr.dtype, copy=False), arr)
            if not exact:
                warnings.warn(
                    "phenotypes were coerced to float64 and some values changed "
                    "(e.g. int overflow or NaN handling).",
                    UserWarning,
                    stacklevel=3,
                )
        except Exception:  # pragma: no cover
            pass
    return coerced


def _coerce_float(values: Any, default: float, n: int) -> np.ndarray:
    if values is None:
        return np.full(n, default, dtype=np.float64)
    arr = np.asarray(values, dtype=np.float64)
    if arr.shape != (n,):
        raise ValueError(f"expected shape ({n},), got {arr.shape}")
    return arr


def _coerce_int(values: Any, default: int, n: int) -> np.ndarray:
    if values is None:
        return np.full(n, default, dtype=np.int64)
    arr = np.asarray(values, dtype=np.int64)
    if arr.shape != (n,):
        raise ValueError(f"expected shape ({n},), got {arr.shape}")
    return arr


class GenotypePhenotypeMap:
    """Container for a genotype-phenotype map. See SCHEMA.md for the full contract."""

    def __init__(
        self,
        wildtype: str,
        genotypes: list[str] | np.ndarray,
        phenotypes: list[float] | np.ndarray,
        *,
        stdeviations: list[float] | np.ndarray | None = None,
        n_replicates: list[int] | np.ndarray | None = None,
        mutations: Mapping[int, list[str] | None] | None = None,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> None:
        if not isinstance(wildtype, str) or len(wildtype) == 0:
            raise ValueError("wildtype must be a non-empty string")
        geno_list = [str(g) for g in np.asarray(genotypes, dtype=object).tolist()]
        if any(len(g) != len(wildtype) for g in geno_list):
            raise ValueError("all genotypes must have the same length as wildtype")
        n = len(geno_list)
        self._wildtype: str = wildtype
        self._genotypes: np.ndarray = np.asarray(geno_list, dtype=object)
        self._phenotypes: np.ndarray = _coerce_phenotypes(phenotypes)
        if self._phenotypes.shape != (n,):
            raise ValueError(f"phenotypes shape {self._phenotypes.shape} != genotypes count ({n},)")
        self._stdeviations: np.ndarray = _coerce_float(stdeviations, float("nan"), n)
        self._n_replicates: np.ndarray = _coerce_int(n_replicates, 1, n)

        if mutations is None:
            mutations = _infer_mutations(wildtype, geno_list)
        self._mutations: dict[int, list[str] | None] = dict(mutations)
        self._site_labels: list[str] = (
            [str(i) for i in range(len(wildtype))]
            if site_labels is None
            else [str(x) for x in site_labels]
        )
        if len(self._site_labels) != len(wildtype):
            raise ValueError(f"site_labels length {len(self._site_labels)} != wildtype length")

        # Cached lazy artifacts.
        self._encoding_table: pd.DataFrame | None = None
        self._binary_packed: np.ndarray | None = None
        self._binary: np.ndarray | None = None
        self._data: pd.DataFrame | None = None

        # Side-view maps.
        self.stdeviation_map = StandardDeviationMap(self)
        self.standard_error_map = StandardErrorMap(self)

        # Build binary right away if requested. Mirrors v1 default.
        if include_binary:
            _ = self.binary_packed

    # ------------------------------------------------------------------ props
    @property
    def wildtype(self) -> str:
        return self._wildtype

    @property
    def genotypes(self) -> np.ndarray:
        return self._genotypes

    @property
    def phenotypes(self) -> np.ndarray:
        return self._phenotypes

    @property
    def stdeviations(self) -> np.ndarray:
        return self._stdeviations

    @property
    def n_replicates(self) -> np.ndarray:
        return self._n_replicates

    @property
    def mutations(self) -> dict[int, list[str] | None]:
        return self._mutations

    @property
    def site_labels(self) -> list[str]:
        return self._site_labels

    @property
    def length(self) -> int:
        return len(self._wildtype)

    @property
    def n(self) -> int:
        return len(self._genotypes)

    @property
    def encoding_table(self) -> pd.DataFrame:
        if self._encoding_table is None:
            self._encoding_table = get_encoding_table(
                self._wildtype, self._mutations, self._site_labels
            )
        return self._encoding_table

    @property
    def binary_packed(self) -> np.ndarray:
        """Shape (n, n_bits) np.uint8 array. Computed once and cached."""
        if self._binary_packed is None:
            self._binary_packed = genotypes_to_binary_packed(
                self._genotypes.tolist(), self.encoding_table
            )
        return self._binary_packed

    @property
    def binary(self) -> np.ndarray:
        """Object-dtype array of '0'/'1' strings, one per genotype. Lazy, cached."""
        if self._binary is None:
            packed = self.binary_packed
            out = np.empty(packed.shape[0], dtype=object)
            for i in range(packed.shape[0]):
                out[i] = "".join("1" if b else "0" for b in packed[i])
            self._binary = out
        return self._binary

    @property
    def n_mutations(self) -> np.ndarray:
        """Per-genotype mutation count (Hamming weight of binary_packed)."""
        return self.binary_packed.sum(axis=1).astype(np.int64)

    @property
    def data(self) -> pd.DataFrame:
        if self._data is None:
            df = pd.DataFrame(
                {
                    "genotypes": pd.array(self._genotypes.tolist(), dtype="string"),
                    "phenotypes": self._phenotypes,
                    "stdeviations": self._stdeviations,
                    "n_replicates": pd.array(self._n_replicates, dtype="Int64"),
                    "binary": pd.array(self.binary.tolist(), dtype="string"),
                    "n_mutations": pd.array(self.n_mutations, dtype="Int64"),
                }
            )
            self._data = df
        return self._data

    # --------------------------------------------------------- constructors
    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        *,
        wildtype: str | None = None,
        mutations: Mapping[int, list[str] | None] | None = None,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> GenotypePhenotypeMap:
        """Hydrate a GenotypePhenotypeMap from a pandas DataFrame.

        Required columns: ``genotypes``, ``phenotypes``. Optional:
        ``stdeviations``, ``n_replicates``. Other columns are ignored.

        ``wildtype`` may be omitted if the DataFrame attrs contain ``"wildtype"``
        (e.g. when hydrating the output of ``gpm.data`` round-tripped through
        ``pd.read_*``). Otherwise it must be provided.
        """
        if "genotypes" not in df.columns or "phenotypes" not in df.columns:
            raise SchemaError("DataFrame must have columns 'genotypes' and 'phenotypes'")
        wt = wildtype or df.attrs.get("wildtype")
        if not wt:
            # Heuristic last resort: a genotype that appears in the frame and
            # has zero Hamming distance to itself, i.e. first row's genotype.
            # Only use this if mutations is also provided so we can verify.
            if mutations is None:
                raise ValueError("wildtype must be provided (or stored in df.attrs['wildtype'])")
            raise ValueError("wildtype must be provided")
        return cls(
            wildtype=str(wt),
            genotypes=df["genotypes"].tolist(),
            phenotypes=df["phenotypes"].to_numpy(),
            stdeviations=(df["stdeviations"].to_numpy() if "stdeviations" in df.columns else None),
            n_replicates=(df["n_replicates"].to_numpy() if "n_replicates" in df.columns else None),
            mutations=mutations,
            site_labels=site_labels,
            include_binary=include_binary,
        )

    # ----------------------------------------------------- missing-genotypes
    def get_missing_genotypes(self) -> np.ndarray:
        """All genotypes in the Cartesian product of `mutations` not present in `genotypes`.

        Enforces the SpaceTooLargeError guard via enumerate.
        """
        from .enumerate import enumerate_genotypes_str

        full = enumerate_genotypes_str(self._wildtype, self._mutations)
        present = set(self._genotypes.tolist())
        missing = [g for g in full if g not in present]
        return np.asarray(missing, dtype=object)

    # ---------------------------------------------------------------- repr
    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"GenotypePhenotypeMap(wildtype={self._wildtype!r}, "
            f"n_genotypes={self.n}, length={self.length})"
        )


# Defensive schema sanity check for callers that construct their own tables.
def validate_encoding_table(table: pd.DataFrame) -> None:
    missing = [c for c in ENCODING_COLUMNS if c not in table.columns]
    if missing:
        raise SchemaError(
            f"encoding_table missing columns {missing}; expected {list(ENCODING_COLUMNS)}"
        )
    if "site_index" in table.columns and table["site_index"].isna().any():
        raise SchemaError("encoding_table.site_index contains NaN (not allowed)")
    if "binary_index_stop" in table.columns and table["binary_index_stop"].isna().any():
        raise SchemaError("encoding_table.binary_index_stop contains NaN")
    _ = n_bits_from_encoding_table(table)
