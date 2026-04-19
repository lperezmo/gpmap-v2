# gpmap-v2 SCHEMA

Contract between `gpmap-v2` (producer) and consumers (`epistasis-v2`, downstream packages).

Status: DRAFT, pending Luis sign-off. File lives at the root of `lperezmo/gpmap-v2` and is load-bearing for any consumer.

Versioning: this document is versioned with the package. Breaking changes to the schema bump the major version. Additive-only changes bump the minor version.

---

## 1. `GenotypePhenotypeMap` public surface

```python
class GenotypePhenotypeMap:
    wildtype: str                                   # length L, each char in alphabet[i]
    genotypes: np.ndarray                           # shape (n,), dtype object (str), length L each
    phenotypes: np.ndarray                          # shape (n,), dtype float64
    stdeviations: np.ndarray                        # shape (n,), dtype float64, default np.nan
    n_replicates: np.ndarray                        # shape (n,), dtype int64, default 1
    mutations: dict[int, list[str] | None]          # site index -> alphabet at that site (None = frozen site)
    site_labels: list[str]                          # length L, string labels
    encoding_table: pd.DataFrame                    # see section 3
    binary: np.ndarray                              # shape (n,), dtype object (str of '0'/'1'), lazy, cached
    binary_packed: np.ndarray                       # shape (n, n_bits), dtype uint8, lazy, cached
    data: pd.DataFrame                              # pandas view built on demand, see section 6

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, *, wildtype: str | None = None,
                       mutations: dict | None = None, site_labels: list[str] | None = None,
                       include_binary: bool = True) -> "GenotypePhenotypeMap": ...

    # I/O
    @classmethod
    def read_csv(cls, path, **kwargs) -> "GenotypePhenotypeMap": ...
    @classmethod
    def read_json(cls, path) -> "GenotypePhenotypeMap": ...
    @classmethod
    def read_pickle(cls, path) -> "GenotypePhenotypeMap": ...
    @classmethod
    def read_excel(cls, path, **kwargs) -> "GenotypePhenotypeMap": ...
    def to_csv(self, path, **kwargs) -> None: ...
    def to_json(self, path) -> None: ...
    def to_pickle(self, path) -> None: ...
    def to_excel(self, path, **kwargs) -> None: ...
    def to_dict(self) -> dict: ...
```

Rules:
- `genotypes[i]` has length `L = len(wildtype)` for all i.
- All arrays are length `n` or shape `(n, ...)` where `n = len(genotypes)`.
- The class is subclassable. `BaseSimulation(GenotypePhenotypeMap)` in simulators stays a supported pattern, but internal code prefers composition over inheritance when extending.

---

## 2. Dtype contract

| Attribute | dtype | Shape | Notes |
|---|---|---|---|
| `phenotypes` | `np.float64` | `(n,)` | coerced at construction; `UserWarning` once if coercion changed values |
| `stdeviations` | `np.float64` | `(n,)` | default `np.nan` when not provided |
| `n_replicates` | `np.int64` | `(n,)` | default 1 when not provided |
| `genotypes` (public) | `np.ndarray[object]` of `str` | `(n,)` | |
| `binary` (public) | `np.ndarray[object]` of `str` | `(n,)` | strings like `"0010..."`; lazy |
| `binary_packed` (public) | `np.uint8` | `(n, n_bits)` | 0 or 1 per cell; lazy |
| `_genotypes_int` (internal) | `np.uint8` | `(n, L)` | 0..alphabet_size(i)-1 per site; lazy |

`n_bits = sum(len(mutations[i]) - 1 for i in range(L) if mutations[i] is not None)`.

---

## 3. `encoding_table` schema

A pandas DataFrame, one row per (site, mutation letter) combination, including one row per site that carries the wildtype letter (`mutation_index` is `<NA>` on that row).

Columns, in this order:

| Column | dtype | Nullable | Meaning |
|---|---|---|---|
| `site_index` | `Int64` (pandas nullable) | no (never NaN) | 0..L-1, which position in the wildtype this row describes |
| `site_label` | `string` | no | user-visible label for the site (defaults to `str(site_index)`) |
| `wildtype_letter` | `string` | no | single character, the wildtype at `site_index` |
| `mutation_letter` | `string` | yes (NaN on frozen sites with no alphabet) | the letter this row describes; for wildtype rows it equals `wildtype_letter` |
| `mutation_index` | `Int64` | yes (NaN on wildtype rows and frozen sites) | global per-mutation index, starts at 1, monotonically increasing across sites |
| `binary_repr` | `string` | no (empty string on frozen sites) | unary-minus-one encoding, length = `alphabet_size - 1`, WT is `"00..0"` |
| `binary_index_start` | `Int64` | no | left edge of this site's slot inside the concatenated binary genotype |
| `binary_index_stop` | `Int64` | no | right edge (exclusive) |

Legacy alias: the column `genotype_index` (v1 name, actually a site index) is available as a read-only alias to `site_index` for one minor version. Writes go to `site_index` only. `DeprecationWarning` on read.

Consumer contract (`epistasis-v2` specifically reads):
- `encoding_table[["mutation_index", "site_index"]].dropna().astype(int)` then `.groupby("site_index")` (equivalent of v1's `groupby("genotype_index")`).
- `encoding_table.dropna()` then `.iterrows()` with access to `wildtype_letter`, `site_label`, `mutation_letter`, `mutation_index`.

Example (wildtype `"AX"`, `mutations={0: ["A","C","G"], 1: None}` where site 1 is frozen):

| site_index | site_label | wildtype_letter | mutation_letter | mutation_index | binary_repr | binary_index_start | binary_index_stop |
|---|---|---|---|---|---|---|---|
| 0 | "0" | "A" | "A" | `<NA>` | "00" | 0 | 2 |
| 0 | "0" | "A" | "C" | 1 | "10" | 0 | 2 |
| 0 | "0" | "A" | "G" | 2 | "01" | 0 | 2 |
| 1 | "1" | "X" | `<NA>` | `<NA>` | "" | 2 | 2 |

Construction cost guarantee: `O(sum alphabet_size_i)`, does not depend on `n_genotypes`. Must be memoized by `(wildtype, mutations, site_labels)` key so epistasis-v2's `add_gpm()` pays it at most once per unique map.

---

## 4. `genotypes_to_binary`

```python
def genotypes_to_binary(
    genotypes: Iterable[str],
    encoding_table: pd.DataFrame,
) -> np.ndarray:
    """Convert iterable of genotype strings to a 1-D np.ndarray[object] of binary strings
    of length n_bits, in the same order as the input.
    """
```

Contract:
- Input `genotypes[i]` is a length-L string over the union of alphabets; unknown letters raise `ValueError`.
- Output `out[i]` is a `str` of length `n_bits`, containing only `'0'` and `'1'`.
- Pure; no side effects; does not touch the `GenotypePhenotypeMap` object.
- Rust-backed (`gpmap._rust.genotypes_to_binary_packed`) for the inner loop, parallelized over rows with rayon. Python wrapper packages the result into the `np.ndarray[object]` of strings for back-compat.
- A sibling `genotypes_to_binary_packed(genotypes, encoding_table) -> np.ndarray[uint8]` of shape `(n, n_bits)` is exposed for consumers who want the packed form directly (epistasis-v2's `build_model_matrix` should prefer this).

Performance target:
- L=20, alphabet_size=4, n=1e6 genotypes: <200 ms on a 2024-class laptop (parallel). Pure-Python v1 equivalent is ~minutes.

---

## 5. Error bar transforms

```python
def upper_transform(bar_y: np.ndarray, upper: np.ndarray, *,
                    log_transform: bool = False, logbase: float = 10.0) -> np.ndarray: ...

def lower_transform(bar_y: np.ndarray, lower: np.ndarray, *,
                    log_transform: bool = False, logbase: float = 10.0) -> np.ndarray: ...
```

Semantics:
- Linear mode (`log_transform=False`): returns the absolute distance from `bar_y` to the bound, i.e. `upper_transform = abs(upper - bar_y)`, `lower_transform = abs(bar_y - lower)`. Matplotlib's `errorbar` expects positive offsets.
- Log mode (`log_transform=True`): returns the log-space offset, i.e. `upper_transform = log_b(upper) - log_b(bar_y)` and `lower_transform = log_b(bar_y) - log_b(lower)`, clipped at 0.
- Both are pure functions, no state, no container. No call back into pandas.

Consumer: `epistasis/pyplot/coefs.py:346-348`. New signatures can be adopted there; we will annotate the change in CHANGELOG when the functions land.

Note: v1 `errors.py` has identical implementations for `upper` and `lower` (copy-paste bug). v2 fixes this; `lower_transform` is truly the lower-bound distance.

---

## 6. `data` property (pandas view)

```python
@property
def data(self) -> pd.DataFrame:
    """Pandas DataFrame with columns: genotypes, phenotypes, stdeviations,
    n_replicates, binary, n_mutations. Built lazily, cached."""
```

Columns, in this order:

| Column | dtype |
|---|---|
| `genotypes` | `string` (object of str) |
| `phenotypes` | `float64` |
| `stdeviations` | `float64` |
| `n_replicates` | `Int64` |
| `binary` | `string` |
| `n_mutations` | `Int64` |

Invariant: `data` is derived state. Mutating `data` does not mutate the container. To change the container, construct a new `GenotypePhenotypeMap` from `from_dataframe(data)`.

---

## 7. Serialization formats

### 7.1 CSV

Columns: `genotypes`, `phenotypes`, `stdeviations`, `n_replicates`. Wildtype and mutations are stored in a sidecar `<basename>.meta.json` with fields `wildtype`, `mutations`, `site_labels`. Round-trip preserves all schema.

### 7.2 JSON

```json
{
  "wildtype": "AAA",
  "mutations": {"0": ["A","T"], "1": ["A","T"], "2": ["A","T"]},
  "site_labels": ["0","1","2"],
  "genotypes": ["AAA", "AAT", ...],
  "phenotypes": [0.1, 0.2, ...],
  "stdeviations": [0.05, 0.05, ...],
  "n_replicates": [1, 1, ...]
}
```

Version marker: top-level `"schema_version": "1"` is added in v2 and required. v1 JSON files are auto-upgraded on read with a one-time `UserWarning`.

### 7.3 Pickle

`pickle.dump(gpm, f)` / `GenotypePhenotypeMap.read_pickle(path)`. The pickled object is the `GenotypePhenotypeMap` instance. We do not change the class path between v2 releases within the same major.

### 7.4 Excel

Same schema as CSV (data sheet) with a second sheet `meta` holding `wildtype`, `mutations`, `site_labels`.

### 7.5 Parquet (optional extra, `gpmap-v2[parquet]`)

Preferred for large maps. Same column set as CSV, metadata sheet collapsed into `pyarrow` table-level metadata (`schema.metadata`).

---

## 8. Size guards

Every function that materializes the full genotype space has a size cap:

- `enumerate_genotypes(alphabet_sizes, max=2**28)`: Rust-backed; refuses if `prod(alphabet_sizes) > max`. Caller passes `allow_huge=True` to override.
- `mutations_to_genotypes(mutations, max=2**28)`: thin Python wrapper around the above.
- `get_missing_genotypes(gpm)`: same cap, applied to the underlying enumeration.

Raised exception: `gpmap.SpaceTooLargeError` (defined in `gpmap.exceptions`), subclass of `ValueError`.

---

## 9. Invariants guaranteed by the container

On successful construction, `GenotypePhenotypeMap` guarantees:

1. `len(wildtype) == L` and every `genotypes[i]` has length `L`.
2. Every character in `genotypes[i][j]` is in `mutations[j]` (if `mutations[j] is not None`), otherwise it equals `wildtype[j]`.
3. `phenotypes.shape == stdeviations.shape == n_replicates.shape == (len(genotypes),)`.
4. `encoding_table` has exactly `sum(alphabet_size_i or 1 for i in range(L))` rows.
5. `binary_packed[i, :].sum()` equals the Hamming distance (in the unary-encoding sense) from `genotypes[i]` to `wildtype`.
6. Round-trip: `from_dataframe(gpm.data).data.equals(gpm.data)` is `True` for every map.

Contract tests in `tests/schema/` enforce all six on every release.

---

## 10. What epistasis-v2 needs to know at a glance

A tiny migration checklist for the epistasis-v2 author:

- Replace `encoding_table[["mutation_index","genotype_index"]]` with `encoding_table[["mutation_index","site_index"]]`. (Or keep the old form during the transition; the alias is live for one minor version.)
- Replace `GenotypePhenotypeMap.read_dataframe(df)` with `GenotypePhenotypeMap.from_dataframe(df)`.
- Prefer `gpm.binary_packed` over `gpm.binary` when building model matrices. Packed form is `np.uint8` 2D, trivially fast to feed into the Rust `build_model_matrix`.
- `upper_transform`/`lower_transform` now honor the sign of the offset and take a `log_transform=False` kwarg. Calls `upper_transform(bar_y, upper)` continue to work.
- Expect `genotypes_to_binary` to be roughly 50-200x faster than v1 on L>=16 maps. If you were splitting calls to cap memory, you can probably stop.

---

## 11. Revision log

- 2026-04-19 - v0 draft by Claude (gpmap-v2 session). Pending Luis review. Pending epistasis-v2 session consumption check.
