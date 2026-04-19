"""Encoding-table construction and genotype -> binary conversion.

The encoding table follows the schema locked in SCHEMA.md section 3.
Column order: site_index, site_label, wildtype_letter, mutation_letter,
mutation_index, binary_repr, binary_index_start, binary_index_stop.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd

from .exceptions import SchemaError, UnknownLetterError

try:
    from gpmap import _rust as _r
except ImportError:  # pragma: no cover
    _r = None  # Rust extension not yet built


ENCODING_COLUMNS: tuple[str, ...] = (
    "site_index",
    "site_label",
    "wildtype_letter",
    "mutation_letter",
    "mutation_index",
    "binary_repr",
    "binary_index_start",
    "binary_index_stop",
)


def _validate_mutations(wildtype: str, mutations: Mapping[int, list[str] | None]) -> None:
    if not isinstance(wildtype, str) or len(wildtype) == 0:
        raise ValueError("wildtype must be a non-empty string")
    expected = set(range(len(wildtype)))
    if set(mutations.keys()) != expected:
        raise ValueError(
            f"mutations keys must be exactly {{0..{len(wildtype) - 1}}}, "
            f"got {sorted(mutations.keys())}"
        )
    for i, alphabet in mutations.items():
        if alphabet is None:
            continue
        if not isinstance(alphabet, list) or len(alphabet) == 0:
            raise ValueError(f"mutations[{i}] must be a non-empty list or None")
        if any(not isinstance(c, str) or len(c) != 1 for c in alphabet):
            raise ValueError(f"mutations[{i}] must be a list of single-char strings")
        if wildtype[i] not in alphabet:
            raise ValueError(
                f"wildtype[{i}] = {wildtype[i]!r} not present in mutations[{i}] = {alphabet!r}"
            )


def _resolve_site_labels(wildtype: str, site_labels: list[str] | None) -> list[str]:
    if site_labels is None:
        return [str(i) for i in range(len(wildtype))]
    if len(site_labels) != len(wildtype):
        raise ValueError(
            f"site_labels length {len(site_labels)} != wildtype length {len(wildtype)}"
        )
    return [str(x) for x in site_labels]


def get_encoding_table(
    wildtype: str,
    mutations: Mapping[int, list[str] | None],
    site_labels: list[str] | None = None,
) -> pd.DataFrame:
    """Build the encoding lookup table described in SCHEMA.md.

    Unary-minus-one encoding: a site with alphabet size k gets k-1 binary bits.
    WT row has `mutation_letter == wildtype_letter`, `mutation_index` is <NA>,
    `binary_repr` is `"0" * (k - 1)`. Non-WT letters each own one bit.
    """
    _validate_mutations(wildtype, mutations)
    labels = _resolve_site_labels(wildtype, site_labels)

    rows: list[dict[str, object]] = []
    mutation_index_counter = 0
    binary_index_counter = 0

    for site_i in range(len(wildtype)):
        alphabet = mutations[site_i]
        wt_letter = wildtype[site_i]

        if alphabet is None:
            rows.append(
                {
                    "site_index": site_i,
                    "site_label": labels[site_i],
                    "wildtype_letter": wt_letter,
                    "mutation_letter": pd.NA,
                    "mutation_index": pd.NA,
                    "binary_repr": "",
                    "binary_index_start": binary_index_counter,
                    "binary_index_stop": binary_index_counter,
                }
            )
            continue

        n_bits_here = len(alphabet) - 1

        # WT row first.
        rows.append(
            {
                "site_index": site_i,
                "site_label": labels[site_i],
                "wildtype_letter": wt_letter,
                "mutation_letter": wt_letter,
                "mutation_index": pd.NA,
                "binary_repr": "0" * n_bits_here,
                "binary_index_start": binary_index_counter,
                "binary_index_stop": binary_index_counter + n_bits_here,
            }
        )

        # Mutant rows in alphabet order, skipping WT.
        non_wt = [c for c in alphabet if c != wt_letter]
        for j, letter in enumerate(non_wt):
            bits = ["0"] * n_bits_here
            bits[j] = "1"
            mutation_index_counter += 1
            rows.append(
                {
                    "site_index": site_i,
                    "site_label": labels[site_i],
                    "wildtype_letter": wt_letter,
                    "mutation_letter": letter,
                    "mutation_index": mutation_index_counter,
                    "binary_repr": "".join(bits),
                    "binary_index_start": binary_index_counter,
                    "binary_index_stop": binary_index_counter + n_bits_here,
                }
            )

        binary_index_counter += n_bits_here

    df = pd.DataFrame(rows, columns=list(ENCODING_COLUMNS))
    df["site_index"] = df["site_index"].astype("Int64")
    df["site_label"] = df["site_label"].astype("string")
    df["wildtype_letter"] = df["wildtype_letter"].astype("string")
    df["mutation_letter"] = df["mutation_letter"].astype("string")
    df["mutation_index"] = df["mutation_index"].astype("Int64")
    df["binary_repr"] = df["binary_repr"].astype("string")
    df["binary_index_start"] = df["binary_index_start"].astype("Int64")
    df["binary_index_stop"] = df["binary_index_stop"].astype("Int64")
    return _LegacyAliasFrame(df)


class _LegacyAliasFrame(pd.DataFrame):
    """DataFrame subclass exposing `genotype_index` as a deprecated alias for `site_index`.

    The alias is surfaced via __getitem__ so users can still do
    `table["genotype_index"]` or `table[["mutation_index","genotype_index"]]`
    during the epistasis-v2 transition, with a DeprecationWarning.
    """

    @property
    def _constructor(self):  # type: ignore[override]
        return _LegacyAliasFrame

    def __getitem__(self, key):  # type: ignore[override]
        import warnings

        if isinstance(key, str) and key == "genotype_index":
            warnings.warn(
                "encoding_table column 'genotype_index' was renamed to 'site_index' in gpmap-v2; "
                "use 'site_index' instead. Alias will be removed in a future minor release.",
                DeprecationWarning,
                stacklevel=2,
            )
            return super().__getitem__("site_index")
        if isinstance(key, list) and "genotype_index" in key:
            warnings.warn(
                "encoding_table column 'genotype_index' was renamed to 'site_index' in gpmap-v2; "
                "use 'site_index' instead. Alias will be removed in a future minor release.",
                DeprecationWarning,
                stacklevel=2,
            )
            key = ["site_index" if k == "genotype_index" else k for k in key]
        return super().__getitem__(key)


def _schema_check(encoding_table: pd.DataFrame) -> None:
    missing = [c for c in ENCODING_COLUMNS if c not in encoding_table.columns]
    if missing:
        raise SchemaError(
            f"encoding_table is missing required columns: {missing}; "
            f"expected {list(ENCODING_COLUMNS)}"
        )


def n_bits_from_encoding_table(encoding_table: pd.DataFrame) -> int:
    _schema_check(encoding_table)
    return int(encoding_table["binary_index_stop"].max())


def _build_rust_inputs(encoding_table: pd.DataFrame):
    """Flatten encoding_table into Rust-friendly structures.

    Returns
    -------
    n_sites, n_bits, letter_to_bits_flat, letter_to_bits_offsets, letter_index
    """
    _schema_check(encoding_table)
    df = encoding_table
    n_sites = int(df["site_index"].max()) + 1
    n_bits = int(df["binary_index_stop"].max())

    # letter_index[site_i] : dict[letter_str, local_row_index]
    letter_index: list[dict[str, int]] = []
    # rows: for site i, all letter rows in alphabet order including WT.
    # We flatten by walking sites in order.
    bits_rows: list[list[int]] = []
    offsets: list[int] = [0]

    for site_i in range(n_sites):
        site_rows = df[df["site_index"] == site_i]
        mapping: dict[str, int] = {}
        for local, (_, row) in enumerate(site_rows.iterrows()):
            letter = row["mutation_letter"]
            if pd.isna(letter):
                continue
            repr_str = row["binary_repr"]
            start = int(row["binary_index_start"])
            stop = int(row["binary_index_stop"])
            vec = [0] * n_bits
            for k, ch in enumerate(repr_str):
                vec[start + k] = 1 if ch == "1" else 0
            # Guard against malformed repr.
            if stop - start != len(repr_str):
                raise SchemaError(
                    f"binary_repr length {len(repr_str)} mismatched for site {site_i}: "
                    f"start={start} stop={stop}"
                )
            mapping[str(letter)] = len(bits_rows) - offsets[-1]
            # Actually, local index is the index within this site's row group; use enumerate.
            # Recompute correctly:
            # local index is the row's position within this site's slice (0-based).
            # But iterrows enumerate gives us that directly if we used `local`.
            mapping[str(letter)] = local
            bits_rows.append(vec)
        offsets.append(len(bits_rows))
        letter_index.append(mapping)

    # Empty site (frozen) can produce no rows; guard that offsets still advance.
    if len(letter_index) != n_sites:
        raise SchemaError("letter_index size mismatch after flattening")

    letter_to_bits_flat = (
        np.asarray(bits_rows, dtype=np.uint8)
        if bits_rows
        else np.zeros((0, n_bits), dtype=np.uint8)
    )
    letter_to_bits_offsets = np.asarray(offsets, dtype=np.int64)

    return n_sites, n_bits, letter_to_bits_flat, letter_to_bits_offsets, letter_index


def genotypes_to_binary_packed(
    genotypes: Iterable[str],
    encoding_table: pd.DataFrame,
) -> np.ndarray:
    """Convert iterable of genotype strings to packed uint8 2D array.

    Returns
    -------
    np.ndarray[uint8] of shape (n_genotypes, n_bits).
    """
    geno_list: list[str] = list(genotypes)
    n_sites, n_bits, flat, offsets, letter_index = _build_rust_inputs(encoding_table)

    if _r is not None:
        try:
            return _r.genotypes_to_binary_packed(
                geno_list,
                n_sites,
                n_bits,
                flat,
                offsets,
                letter_index,
            )
        except ValueError as exc:
            # Translate Rust-side "unknown letter" to UnknownLetterError.
            msg = str(exc)
            if "unknown letter" in msg:
                raise UnknownLetterError(msg) from exc
            raise

    # Pure-Python fallback (for developers running without maturin develop).
    out = np.zeros((len(geno_list), n_bits), dtype=np.uint8)
    for i, g in enumerate(geno_list):
        if len(g) != n_sites:
            raise ValueError(f"genotype {g!r} has length {len(g)}, expected {n_sites}")
        for site_i, letter in enumerate(g):
            mapping = letter_index[site_i]
            if letter not in mapping:
                raise UnknownLetterError(f"unknown letter {letter!r} at site {site_i}")
            local = mapping[letter]
            flat_row = int(offsets[site_i]) + local
            out[i] |= flat[flat_row]
    return out


def genotypes_to_binary(
    genotypes: Iterable[str],
    encoding_table: pd.DataFrame,
) -> np.ndarray:
    """Convert iterable of genotype strings to a 1-D np.ndarray[object] of binary strings.

    Back-compat surface: returns an array of strings like "010101..." to match
    epistasis-v2's consumer pattern. Prefer `genotypes_to_binary_packed` for new code.
    """
    packed = genotypes_to_binary_packed(genotypes, encoding_table)
    out = np.empty(packed.shape[0], dtype=object)
    for i in range(packed.shape[0]):
        out[i] = "".join("1" if b else "0" for b in packed[i])
    return out
