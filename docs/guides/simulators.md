---
title: "Simulators"
description: "Toy fitness landscapes: NK, Mount Fuji, multi-peak Fuji, House of Cards, random phenotypes, masked subsamples."
---

# Simulators

`gpmap.simulate` ships a small zoo of toy-landscape generators that subclass `GenotypePhenotypeMap`. Each generator enumerates the full Cartesian product of per-site alphabets and fills in `phenotypes` according to its model. Use them when you need a deterministic landscape for testing fits, benchmarking, or pedagogy.

## Shared interface

Every simulator takes at minimum a wildtype and a mutations dict. Optional kwargs:

- `rng`: a `np.random.Generator`. Defaults to `np.random.default_rng()`. Pass a seeded generator for reproducibility.
- `site_labels`: optional per-site labels.
- `include_binary`: whether to build `binary_packed` eagerly (default `True`).

All simulators are `GenotypePhenotypeMap` subclasses, so the full container API (`gpm.data`, `gpm.binary_packed`, `to_json`, etc.) is available on every instance.

## NKSimulation

The Kauffman NK model: each site's fitness contribution depends on its K nearest neighbors (wrap-around). Total phenotype is the average per-site contribution.

```python
import numpy as np
from gpmap.simulate import NKSimulation

sim = NKSimulation(
    wildtype="AAAA",
    mutations={i: ["A", "T"] for i in range(4)},
    K=2,
    rng=np.random.default_rng(0),
)
sim.phenotypes.shape  # (16,)
```

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `K` | `int` | `1` | Neighborhood size, in bits |
| `rng` | `np.random.Generator \| None` | `None` | Seed source for the per-neighborhood random table |

Rough phenotype range: `[0, 1]`. Higher K means rougher landscapes with more local optima.

![NK landscapes at K=0 (smooth, additive) versus K=3 (rugged, many local optima)](../assets/nk-panel-light.png#only-light)
![NK landscapes at K=0 (smooth, additive) versus K=3 (rugged, many local optima)](../assets/nk-panel-dark.png#only-dark)

## MountFujiSimulation

A single-peak additive landscape with optional Gaussian or uniform roughness. The peak is at the wildtype.

```python
from gpmap.simulate import MountFujiSimulation

sim = MountFujiSimulation(
    wildtype="AAAA",
    mutations={i: ["A", "T"] for i in range(4)},
    field_strength=1.0,
    roughness_width=0.1,
    roughness_dist="normal",
)
```

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `field_strength` | `float` | `1.0` | Slope of the additive component (units: phenotype per Hamming step from WT) |
| `roughness_width` | `float` | `0.0` | Standard deviation (`"normal"`) or half-width (`"uniform"`) of the noise term |
| `roughness_dist` | `"normal" \| "uniform"` | `"normal"` | Noise distribution |

Phenotype is `-field_strength * hamming_distance + noise`, so the wildtype is the global maximum at `0 + noise`.

## MultiPeakMountFujiSimulation

A max-of-Fujis landscape: pick `peak_n` peaks at random (subject to a minimum pairwise Hamming distance), build a single-peak Fuji around each, and the genotype's phenotype is the max over those.

```python
from gpmap.simulate import MultiPeakMountFujiSimulation

sim = MultiPeakMountFujiSimulation(
    wildtype="AAAA",
    mutations={i: ["A", "T"] for i in range(4)},
    peak_n=3,
    min_peak_distance=2,
)

sim.peak_genotypes  # list[str] of length 3
```

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `peak_n` | `int` | `2` | Number of peaks to plant |
| `min_peak_distance` | `int` | `1` | Minimum Hamming distance between any two peaks |
| `max_peak_distance` | `int \| None` | `None` (= `len(wildtype)`) | Maximum allowed pairwise distance |
| `max_proposal_retries` | `int` | `10_000` | Hard cap on peak-search retries before raising |

!!! note "Retry cap"

    v1 had an unguarded while loop that could spin forever on infeasible constraints. v2 caps the search and raises `RuntimeError` with a clear message if it cannot place all peaks.

## HouseOfCardsSimulation

NK with `K = n_bits - 1`: every site sees the full genotype, so neighboring genotypes are uncorrelated. The roughest possible NK landscape.

```python
from gpmap.simulate import HouseOfCardsSimulation

sim = HouseOfCardsSimulation(
    wildtype="AAAA",
    mutations={i: ["A", "T"] for i in range(4)},
)
```

No tunable knobs beyond `rng`; the K is fixed by the geometry.

![Phenotype versus Hamming distance for a smooth Mount Fuji landscape against a fully random House of Cards landscape](../assets/landscape-compare-light.png#only-light)
![Phenotype versus Hamming distance for a smooth Mount Fuji landscape against a fully random House of Cards landscape](../assets/landscape-compare-dark.png#only-dark)

The same contrast on the graph itself: a single-peak Fuji funnels toward one optimum, while House of Cards scatters local peaks (red rings) across the map.

![Mount Fuji versus House of Cards drawn as genotype graphs, local peaks ringed in red](../assets/landscape-graphs-light.png#only-light)
![Mount Fuji versus House of Cards drawn as genotype graphs, local peaks ringed in red](../assets/landscape-graphs-dark.png#only-dark)

## RandomPhenotypesSimulation

Phenotypes are drawn uniformly from `[low, high]`. The simplest possible landscape, useful as a null model.

```python
from gpmap.simulate import RandomPhenotypesSimulation

sim = RandomPhenotypesSimulation(
    wildtype="AAAA",
    mutations={i: ["A", "T"] for i in range(4)},
    low=0.0,
    high=1.0,
)
```

## Build from length

Every simulator that subclasses `BaseSimulation` exposes a `.from_length` constructor for the case where you do not care about the specific alphabet, just its size:

```python
import random
from gpmap.simulate import NKSimulation

sim = NKSimulation.from_length(
    length=5,
    alphabet_size=3,
    kind="DNA",
    K=2,
    rng=random.Random(0),
)
```

| Parameter | Type | Meaning |
|---|---|---|
| `length` | `int` | Number of sites |
| `alphabet_size` | `int` | Per-site alphabet size (>= 2) |
| `kind` | `"AA" \| "DNA" \| "RNA" \| "BINARY"` | Source alphabet to sample from |

The wildtype is the first letter of each site's randomly-chosen alphabet.

## Masking: subsample a map

`gpmap.simulate.mask` returns a `(true_fraction, GenotypePhenotypeMap)` named tuple holding a fraction of the input genotypes, sampled uniformly without replacement:

```python
from gpmap.simulate import mask

masked = mask(sim, fraction=0.5, rng=np.random.default_rng(0))
masked.fraction      # actual fraction kept (rounded)
masked.gpm           # subsampled GenotypePhenotypeMap
```

!!! tip "Named-tuple return"

    v1 returned `(float, GPM)`. v2 returns a `MaskedGPM` named tuple so you can index by name (`.fraction`, `.gpm`) instead of remembering the order.

## Adding noise to a simulator

Simulators expose a `.set_stdeviations(sigma)` method so you can attach a uniform measurement uncertainty after construction:

```python
sim = NKSimulation(wildtype="AAAA", mutations={i: ["A", "T"] for i in range(4)})
sim.set_stdeviations(0.05)
sim.stdeviations  # array([0.05, 0.05, ...])
```

For per-genotype noise, build the simulator, copy `gpm.phenotypes`, perturb, and feed the result back through `GenotypePhenotypeMap(...)`.
