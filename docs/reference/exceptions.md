---
title: "Exceptions"
description: "GPMapError, SpaceTooLargeError, SchemaError, UnknownLetterError."
---

# `gpmap.exceptions`

A small exception hierarchy. All gpmap-specific errors inherit from `GPMapError`; the value-error subclasses also inherit from `ValueError` so `except ValueError` keeps working.

## Hierarchy

```
Exception
└── GPMapError
    ├── SpaceTooLargeError (also ValueError)
    ├── SchemaError (also ValueError)
    └── UnknownLetterError (also ValueError)
```

## `GPMapError`

```python
class GPMapError(Exception):
    """Base class for gpmap errors."""
```

Catch-all for any gpmap-specific failure. Useful when you want to differentiate gpmap bugs from generic `ValueError` raised by NumPy or pandas.

## `SpaceTooLargeError`

```python
class SpaceTooLargeError(GPMapError, ValueError):
    """Raised when a caller asks to materialize a genotype space above the safety cap."""
```

Raised by `enumerate_genotypes_int`, `enumerate_genotypes_str`, and any function that materializes the full Cartesian product when it would exceed `max_genotypes`. Default cap is `2**28` (about 268 million genotypes). Pass `allow_huge=True` to the calling function to bypass.

## `SchemaError`

```python
class SchemaError(GPMapError, ValueError):
    """Raised when an encoding_table / DataFrame does not match the locked schema."""
```

Raised by `validate_encoding_table` and by internal flatten routines when an encoding table is missing required columns, has NaN in `site_index` or `binary_index_stop`, or carries `binary_repr` strings whose length disagrees with `binary_index_stop - binary_index_start`.

Also raised by `from_dataframe` when the DataFrame lacks the required `genotypes` or `phenotypes` columns.

## `UnknownLetterError`

```python
class UnknownLetterError(GPMapError, ValueError):
    """Raised when a genotype contains a letter not present in the per-site alphabet."""
```

Raised by `genotypes_to_binary_packed` (and indirectly by `genotypes_to_binary`) when a genotype contains a letter that the encoding table does not know about. The error message includes the site index and the offending letter.

## Catching pattern

```python
from gpmap import SpaceTooLargeError, SchemaError, UnknownLetterError

try:
    gpm.get_missing_genotypes()
except SpaceTooLargeError as e:
    log.warning("Space too large, skipping: %s", e)
except SchemaError as e:
    log.error("Bad encoding table: %s", e)
```

For library code that wants to swallow only gpmap-related failures, catch `GPMapError` and re-raise the rest.
