"""Custom exceptions for gpmap."""

from __future__ import annotations


class GPMapError(Exception):
    """Base class for gpmap errors."""


class SpaceTooLargeError(GPMapError, ValueError):
    """Raised when a caller asks to materialize a genotype space above the safety cap."""


class SchemaError(GPMapError, ValueError):
    """Raised when an encoding_table / DataFrame does not match the locked schema."""


class UnknownLetterError(GPMapError, ValueError):
    """Raised when a genotype contains a letter not present in the per-site alphabet."""
