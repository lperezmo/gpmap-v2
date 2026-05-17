---
title: "Statistics"
description: "Unbiased standard deviation, variance, and standard error with c4 correction."
---

# `gpmap.stats`

Statistical utilities for unbiased estimators. v1 had two bugs (`unbiased_var` ignored its `axis` kwarg, plus dead `coverage()` stub). v2 honors `axis` and uses the Bessel `c4` correction for unbiased standard deviation.

## `c4_correction`

```python
def c4_correction(n_samples: int) -> float
```

Bessel `c4` factor for unbiased standard deviation estimation. Numerically stable via `lgamma` for every `n`.

- `n <= 1`: returns `1.0` (no correction possible).
- `n > 170`: uses `lgamma` to avoid factorial overflow.

## `unbiased_std`

```python
def unbiased_std(x: np.ndarray, axis: int | None = None) -> np.ndarray | float
```

Unbiased standard deviation: `np.std(x, ddof=1) / c4(n)` where `n = x.shape[axis]` or `x.size` if `axis is None`.

## `unbiased_var`

```python
def unbiased_var(x: np.ndarray, axis: int | None = None) -> np.ndarray | float
```

Unbiased variance: `unbiased_std(x, axis) ** 2`.

!!! note "v1 bug fixed"

    v1 had `unbiased_var(x, axis=None)` that always hardcoded `axis=1` internally. v2 honors the `axis` kwarg correctly.

## `unbiased_sterror`

```python
def unbiased_sterror(x: np.ndarray, axis: int | None = None) -> np.ndarray | float
```

Unbiased standard error: `unbiased_std(x, axis) / sqrt(n)`.

## `corrected_std` and `corrected_sterror`

```python
def corrected_std(biased_std: float | np.ndarray, n_samples: int) -> float | np.ndarray

def corrected_sterror(biased_sterror: float | np.ndarray, n_samples: int) -> float | np.ndarray
```

Given a biased (ddof=0 style) estimator and the sample count, recover the c4-corrected unbiased estimator. `corrected_sterror` is a thin alias for `corrected_std`; both apply the same algebraic correction.

## Use case

The simulator subclasses store per-genotype `n_replicates`. When you want to recover an unbiased per-genotype standard error from a biased measurement, pass it through `corrected_std`:

```python
import numpy as np
from gpmap.stats import corrected_std

biased = np.std(samples, ddof=0)
unbiased = corrected_std(biased, n_samples=samples.size)
```

For an array of biased stds, the function broadcasts.
