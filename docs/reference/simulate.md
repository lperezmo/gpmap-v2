---
title: "Simulators reference"
description: "API reference for NK, Mount Fuji, multi-peak Fuji, House of Cards, random phenotypes, and the masking helper."
---

# `gpmap.simulate`

Six toy-landscape generators plus a masking helper. See [Simulators guide](../guides/simulators.md) for usage examples.

## `BaseSimulation`

```python
class BaseSimulation(GenotypePhenotypeMap):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        site_labels: list[str] | None = None,
        include_binary: bool = True,
    ) -> None: ...

    def build(self) -> None:
        """Subclass hook: fill self._phenotypes[:] = ..."""

    @classmethod
    def from_length(
        cls,
        length: int,
        alphabet_size: int = 2,
        kind: str = "BINARY",
        *,
        rng: random.Random | None = None,
        **kwargs: object,
    ) -> "BaseSimulation": ...

    def set_stdeviations(self, sigma: float) -> None: ...
```

`BaseSimulation` enumerates the full Cartesian product of `mutations`, calls `super().__init__(...)` to set up the container, and then calls `self.build()` to fill `phenotypes`. Subclasses override `build()`.

`from_length(length, alphabet_size, kind)` builds the simulator with a randomly chosen per-site alphabet drawn from `"AA"`, `"DNA"`, `"RNA"`, or `"BINARY"` source. Wildtype is the first letter of each site's chosen alphabet.

`set_stdeviations(sigma)` assigns a uniform standard deviation across all genotypes.

## `NKSimulation`

```python
class NKSimulation(BaseSimulation):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        K: int = 1,
        rng: np.random.Generator | None = None,
        ...
    ) -> None: ...
```

Kauffman NK model: each bit's fitness contribution depends on a wrap-around window of K nearest bits. Phenotypes are normalized to roughly `[0, 1]` via averaging.

## `MountFujiSimulation`

```python
class MountFujiSimulation(BaseSimulation):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        field_strength: float = 1.0,
        roughness_width: float = 0.0,
        roughness_dist: str = "normal",
        rng: np.random.Generator | None = None,
        ...
    ) -> None: ...
```

Single-peak additive landscape. The phenotype at genotype g is `-field_strength * hamming(g, wildtype) + noise`. WT is the peak at `phenotype = noise`.

| Parameter | Type | Notes |
|---|---|---|
| `field_strength` | `float` | Slope of the additive component |
| `roughness_width` | `float` | Std (normal) or half-width (uniform) of the noise |
| `roughness_dist` | `"normal" \| "uniform"` | Noise distribution |

## `MultiPeakMountFujiSimulation`

```python
class MultiPeakMountFujiSimulation(BaseSimulation):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        peak_n: int = 2,
        min_peak_distance: int = 1,
        max_peak_distance: int | None = None,
        field_strength: float = 1.0,
        roughness_width: float = 0.0,
        roughness_dist: str = "normal",
        max_proposal_retries: int = 10_000,
        rng: np.random.Generator | None = None,
        ...
    ) -> None: ...

    @property
    def peak_genotypes(self) -> list[str]: ...
```

Max-of-Fujis landscape. Picks `peak_n` genotypes as peaks subject to pairwise Hamming-distance constraints, builds a single-peak Fuji around each, takes the per-genotype max. Raises `RuntimeError` after `max_proposal_retries` if the constraints are infeasible.

## `HouseOfCardsSimulation`

```python
class HouseOfCardsSimulation(NKSimulation):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        rng: np.random.Generator | None = None,
        ...
    ) -> None: ...
```

Special-case NK with `K = n_bits - 1`. Every site's fitness contribution depends on the entire genotype, so neighboring genotypes are uncorrelated. The roughest possible NK landscape.

## `RandomPhenotypesSimulation`

```python
class RandomPhenotypesSimulation(BaseSimulation):
    def __init__(
        self,
        wildtype: str,
        mutations: Mapping[int, list[str] | None],
        *,
        low: float = 0.0,
        high: float = 1.0,
        rng: np.random.Generator | None = None,
        ...
    ) -> None: ...
```

Phenotypes drawn uniformly from `[low, high]`. Simplest possible null model.

## `random_mutation_set`

```python
def random_mutation_set(
    length: int,
    alphabet_size: int = 2,
    kind: str = "AA",
    *,
    rng: random.Random | None = None,
) -> dict[int, list[str]]
```

Build a random per-site alphabet from one of the canonical letter pools.

| `kind` | Source pool |
|---|---|
| `"AA"` | 20-letter amino-acid alphabet |
| `"DNA"` | `A`, `C`, `G`, `T` |
| `"RNA"` | `A`, `C`, `G`, `U` |
| `"BINARY"` | `0`, `1` |

The source pool is sliced and shuffled per call; v1's in-place mutation of the module-level constant is fixed.

## `mask` and `MaskedGPM`

```python
class MaskedGPM(NamedTuple):
    fraction: float
    gpm: GenotypePhenotypeMap


def mask(
    gpm: GenotypePhenotypeMap,
    fraction: float,
    *,
    rng: np.random.Generator | None = None,
) -> MaskedGPM
```

Subsample `gpm` by keeping a uniform random fraction of its genotypes. Returns a named tuple so callers can index by field instead of by position.

`fraction` must be in `(0, 1]`. The actual fraction kept (`MaskedGPM.fraction`) is `keep_count / n` where `keep_count = max(1, round(fraction * n))`.
