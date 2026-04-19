"""Shared pytest fixtures."""

from __future__ import annotations

import numpy as np
import pytest
from gpmap import GenotypePhenotypeMap


@pytest.fixture
def tiny_binary_gpm() -> GenotypePhenotypeMap:
    return GenotypePhenotypeMap(
        wildtype="AAA",
        genotypes=[
            "AAA",
            "AAT",
            "ATA",
            "TAA",
            "ATT",
            "TAT",
            "TTA",
            "TTT",
        ],
        phenotypes=[0.1, 0.2, 0.2, 0.6, 0.4, 0.6, 1.0, 1.1],
        stdeviations=[0.05] * 8,
    )


@pytest.fixture
def aa_gpm() -> GenotypePhenotypeMap:
    return GenotypePhenotypeMap(
        wildtype="AK",
        genotypes=["AK", "CK", "AR", "CR"],
        phenotypes=[0.0, 0.5, 0.25, 1.0],
        mutations={0: ["A", "C"], 1: ["K", "R"]},
    )


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng(0x5EED)
