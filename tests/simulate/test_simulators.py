"""Smoke tests for every simulator."""

from __future__ import annotations

import numpy as np
import pytest
from gpmap.simulate import (
    BaseSimulation,
    HouseOfCardsSimulation,
    MountFujiSimulation,
    MultiPeakMountFujiSimulation,
    NKSimulation,
    RandomPhenotypesSimulation,
    mask,
    random_mutation_set,
)


def test_random_mutation_set_does_not_mutate_module_state() -> None:
    from gpmap.simulate.base import AMINO_ACIDS

    snapshot = tuple(AMINO_ACIDS)
    _ = random_mutation_set(length=4, alphabet_size=3, kind="AA")
    _ = random_mutation_set(length=4, alphabet_size=3, kind="AA")
    # Module-level tuple never changes.
    assert tuple(AMINO_ACIDS) == snapshot


def test_random_mutation_set_sizes() -> None:
    m = random_mutation_set(length=5, alphabet_size=4, kind="AA")
    assert len(m) == 5
    for alphabet in m.values():
        assert len(alphabet) == 4


def test_random_mutation_set_rejects_bad_size() -> None:
    with pytest.raises(ValueError):
        random_mutation_set(length=3, alphabet_size=25, kind="AA")


def test_random_phenotypes_simulator_shape() -> None:
    sim = RandomPhenotypesSimulation(
        wildtype="AA",
        mutations={0: ["A", "T"], 1: ["A", "T"]},
        low=0.0,
        high=1.0,
        rng=np.random.default_rng(0),
    )
    assert sim.n == 4
    assert sim.phenotypes.min() >= 0.0
    assert sim.phenotypes.max() <= 1.0


def test_nk_simulator_invariants() -> None:
    sim = NKSimulation(
        wildtype="AAAA",
        mutations={i: ["A", "T"] for i in range(4)},
        K=2,
        rng=np.random.default_rng(1),
    )
    assert sim.n == 16
    assert sim.phenotypes.shape == (16,)
    assert np.isfinite(sim.phenotypes).all()


def test_house_of_cards_smoke() -> None:
    sim = HouseOfCardsSimulation(
        wildtype="AAA",
        mutations={i: ["A", "T"] for i in range(3)},
        rng=np.random.default_rng(2),
    )
    assert sim.n == 8
    assert np.isfinite(sim.phenotypes).all()


def test_mount_fuji_phenotype_decreases_with_hamming() -> None:
    # Use field_strength > 0 and zero roughness; phenotype should be -|hamming|.
    sim = MountFujiSimulation(
        wildtype="AAAA",
        mutations={i: ["A", "T"] for i in range(4)},
        field_strength=1.0,
        roughness_width=0.0,
        rng=np.random.default_rng(3),
    )
    hammings = sim.n_mutations
    expected = -hammings.astype(np.float64)
    np.testing.assert_allclose(sim.phenotypes, expected)


def test_multipeak_mount_fuji_produces_finite_phenotypes() -> None:
    sim = MultiPeakMountFujiSimulation(
        wildtype="AAAA",
        mutations={i: ["A", "T"] for i in range(4)},
        peak_n=2,
        min_peak_distance=1,
        max_peak_distance=4,
        field_strength=1.0,
        rng=np.random.default_rng(4),
    )
    assert np.isfinite(sim.phenotypes).all()
    assert len(sim.peak_genotypes) == 2


def test_multipeak_retry_cap_triggers() -> None:
    with pytest.raises(RuntimeError, match="failed to find"):
        MultiPeakMountFujiSimulation(
            wildtype="A",
            mutations={0: ["A", "T"]},
            peak_n=5,  # Impossible with only 2 genotypes.
            max_proposal_retries=50,
            rng=np.random.default_rng(5),
        )


def test_mask_returns_named_tuple(tiny_binary_gpm) -> None:
    rng = np.random.default_rng(6)
    result = mask(tiny_binary_gpm, 0.5, rng=rng)
    assert result.fraction == 0.5
    assert result.gpm.n == 4


def test_mask_rejects_bad_fraction(tiny_binary_gpm) -> None:
    with pytest.raises(ValueError):
        mask(tiny_binary_gpm, 0.0)
    with pytest.raises(ValueError):
        mask(tiny_binary_gpm, 1.5)


def test_base_simulation_subclass_must_implement_build() -> None:
    class BadSim(BaseSimulation):
        pass

    with pytest.raises(NotImplementedError):
        BadSim(wildtype="A", mutations={0: ["A", "T"]})


def test_base_from_length() -> None:
    sim = RandomPhenotypesSimulation.from_length(
        length=4,
        alphabet_size=2,
        kind="BINARY",
        rng=None,
    )
    assert sim.length == 4
    assert sim.n == 16
