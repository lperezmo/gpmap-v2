use numpy::{IntoPyArray, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub fn enumerate_genotypes<'py>(
    py: Python<'py>,
    alphabet_sizes: Vec<usize>,
    max_genotypes: usize,
    allow_huge: bool,
) -> PyResult<Bound<'py, PyArray2<u8>>> {
    let l = alphabet_sizes.len();
    if l == 0 {
        return Err(PyValueError::new_err("alphabet_sizes is empty"));
    }
    for (i, &s) in alphabet_sizes.iter().enumerate() {
        if s == 0 {
            return Err(PyValueError::new_err(format!(
                "alphabet_sizes[{}] is 0",
                i
            )));
        }
        if s > 255 {
            return Err(PyValueError::new_err(format!(
                "alphabet_sizes[{}] = {} exceeds u8 max (255)",
                i, s
            )));
        }
    }

    // Count genotypes with overflow detection.
    let mut total: u128 = 1;
    for &s in &alphabet_sizes {
        total = total.checked_mul(s as u128).ok_or_else(|| {
            PyValueError::new_err("alphabet size product overflowed u128")
        })?;
    }
    if !allow_huge && total > max_genotypes as u128 {
        return Err(PyValueError::new_err(format!(
            "enumerate_genotypes would produce {} rows, exceeds max_genotypes = {}; \
             pass allow_huge=True to override",
            total, max_genotypes
        )));
    }
    if total > usize::MAX as u128 {
        return Err(PyValueError::new_err(
            "alphabet size product exceeds usize::MAX even with allow_huge",
        ));
    }
    let n = total as usize;

    // Precompute strides, last site changes fastest.
    let mut strides: Vec<usize> = vec![1; l];
    for i in (0..l - 1).rev() {
        strides[i] = strides[i + 1] * alphabet_sizes[i + 1];
    }

    let out_flat: Vec<u8> = py.allow_threads(|| {
        let mut out = vec![0u8; n * l];
        for row in 0..n {
            let base = row * l;
            let mut rem = row;
            for site in 0..l {
                let letter = rem / strides[site];
                rem %= strides[site];
                out[base + site] = letter as u8;
            }
        }
        out
    });

    let arr = ndarray::Array2::from_shape_vec((n, l), out_flat)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(arr.into_pyarray_bound(py))
}
