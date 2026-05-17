---
title: "Error transforms"
description: "upper_transform, lower_transform, StandardDeviationMap, StandardErrorMap. Helpers for matplotlib error bars."
---

# `gpmap.errors`

Error-bar transforms and `StandardDeviationMap` / `StandardErrorMap` views. Pinned by `SCHEMA.md` section 5.

## `upper_transform`

```python
def upper_transform(
    bar_y: np.ndarray,
    upper: np.ndarray,
    *,
    log_transform: bool = False,
    logbase: float = 10.0,
) -> np.ndarray
```

Matplotlib-style positive offset from `bar_y` up to the upper bound.

- **Linear** (`log_transform=False`): returns `abs(upper - bar_y)`.
- **Log** (`log_transform=True`): returns `log_b(upper) - log_b(bar_y)`, clipped at 0.

`matplotlib.pyplot.errorbar` expects positive offsets, so this returns the magnitude rather than the signed difference.

## `lower_transform`

```python
def lower_transform(
    bar_y: np.ndarray,
    lower: np.ndarray,
    *,
    log_transform: bool = False,
    logbase: float = 10.0,
) -> np.ndarray
```

Matplotlib-style positive offset from `bar_y` down to the lower bound.

- **Linear**: returns `abs(bar_y - lower)`.
- **Log**: returns `log_b(bar_y) - log_b(lower)`, clipped at 0.

!!! note "v1 bug fixed"

    v1's `errors.py` had identical implementations of `upper_transform` and `lower_transform` (copy-paste bug). v2's `lower_transform` is genuinely the lower-bound distance.

## `StandardDeviationMap`

```python
class StandardDeviationMap:
    def __init__(self, gpm: GenotypePhenotypeMap) -> None: ...

    @property
    def upper(self) -> np.ndarray: ...   # raw stdeviations
    @property
    def lower(self) -> np.ndarray: ...   # raw stdeviations
```

A thin view over a `GenotypePhenotypeMap` that exposes the raw `stdeviations` as both upper and lower error bars. Both bounds are the same value: this is just a name for symmetric error bars at `+/- sigma`.

Available as `gpm.stdeviation_map`.

## `StandardErrorMap`

```python
class StandardErrorMap:
    def __init__(self, gpm: GenotypePhenotypeMap) -> None: ...

    @property
    def upper(self) -> np.ndarray: ...   # stdeviations / sqrt(n_replicates)
    @property
    def lower(self) -> np.ndarray: ...   # stdeviations / sqrt(n_replicates)
```

A view that returns standard errors (`std / sqrt(n_replicates)`) as symmetric error bars. Available as `gpm.standard_error_map`.

Use this when you want error bars on a mean rather than on individual measurements.

## Pattern: plotting with matplotlib

```python
import matplotlib.pyplot as plt
from gpmap import upper_transform, lower_transform

bar_y = gpm.phenotypes
upper = gpm.phenotypes + gpm.stdeviations
lower = gpm.phenotypes - gpm.stdeviations

yerr = [lower_transform(bar_y, lower), upper_transform(bar_y, upper)]
plt.errorbar(range(len(bar_y)), bar_y, yerr=yerr, fmt="o")
```

Or, equivalently, lean on the `stdeviation_map` view:

```python
plt.errorbar(
    range(gpm.n),
    gpm.phenotypes,
    yerr=[gpm.stdeviation_map.lower, gpm.stdeviation_map.upper],
    fmt="o",
)
```
