"""gpmap-v2 public API.

See SCHEMA.md at the repo root for the contract epistasis-v2 depends on.
"""

from __future__ import annotations

from ._version import __version__
from .core import GenotypePhenotypeMap, validate_encoding_table
from .encoding import (
    ENCODING_COLUMNS,
    genotypes_to_binary,
    genotypes_to_binary_packed,
    get_encoding_table,
)
from .enumerate import enumerate_genotypes_int, enumerate_genotypes_str
from .errors import (
    StandardDeviationMap,
    StandardErrorMap,
    lower_transform,
    upper_transform,
)
from .exceptions import (
    GPMapError,
    SchemaError,
    SpaceTooLargeError,
    UnknownLetterError,
)
from .io import (
    read_csv,
    read_excel,
    read_json,
    read_pickle,
    to_csv,
    to_excel,
    to_json,
    to_pickle,
)

__all__ = [
    "ENCODING_COLUMNS",
    "GPMapError",
    "GenotypePhenotypeMap",
    "SchemaError",
    "SpaceTooLargeError",
    "StandardDeviationMap",
    "StandardErrorMap",
    "UnknownLetterError",
    "__version__",
    "enumerate_genotypes_int",
    "enumerate_genotypes_str",
    "genotypes_to_binary",
    "genotypes_to_binary_packed",
    "get_encoding_table",
    "lower_transform",
    "read_csv",
    "read_excel",
    "read_json",
    "read_pickle",
    "to_csv",
    "to_excel",
    "to_json",
    "to_pickle",
    "upper_transform",
    "validate_encoding_table",
]
