"""Toy-landscape simulators built on GenotypePhenotypeMap."""

from __future__ import annotations

from .base import BaseSimulation, random_mutation_set
from .fuji import MountFujiSimulation
from .hoc import HouseOfCardsSimulation
from .mask import MaskedGPM, mask
from .multipeak_fuji import MultiPeakMountFujiSimulation
from .nk import NKSimulation
from .random import RandomPhenotypesSimulation

__all__ = [
    "BaseSimulation",
    "HouseOfCardsSimulation",
    "MaskedGPM",
    "MountFujiSimulation",
    "MultiPeakMountFujiSimulation",
    "NKSimulation",
    "RandomPhenotypesSimulation",
    "mask",
    "random_mutation_set",
]
