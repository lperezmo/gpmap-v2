//! gpmap-core: Rust hot-path primitives exposed to Python via PyO3.
//!
//! This crate is loaded as `gpmap._rust` inside the Python package.
//! See `SCHEMA.md` at the repo root for the contract.

use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::exceptions::{PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyList};
use rayon::prelude::*;

mod encoding;
mod enumerate;
mod hamming;

/// Convert iterable of genotype strings to packed uint8 2D array.
///
/// Parameters
/// ----------
/// genotypes : list[str]
/// site_offsets : np.ndarray[int64] shape (L + 1,)
///     site i covers letters at positions site_offsets[i]..site_offsets[i+1].
/// letter_to_bits_flat : np.ndarray[uint8] shape (sum n_alphabet_i, n_bits)
/// letter_to_bits_offsets : np.ndarray[int64] shape (L + 1,)
///     site i's letter rows live at rows [letter_to_bits_offsets[i], letter_to_bits_offsets[i+1]).
/// letter_index : list[dict[str,int]]
///     per-site mapping from single-character letter to local letter-row index.
///
/// Returns
/// -------
/// np.ndarray[uint8] shape (n_genotypes, n_bits).
#[pyfunction]
#[pyo3(signature = (genotypes, n_sites, n_bits, letter_to_bits_flat, letter_to_bits_offsets, letter_index))]
fn genotypes_to_binary_packed<'py>(
    py: Python<'py>,
    genotypes: &Bound<'py, PyList>,
    n_sites: usize,
    n_bits: usize,
    letter_to_bits_flat: PyReadonlyArray2<'py, u8>,
    letter_to_bits_offsets: numpy::PyReadonlyArray1<'py, i64>,
    letter_index: &Bound<'py, PyList>,
) -> PyResult<Bound<'py, PyArray2<u8>>> {
    let n = genotypes.len();

    // Copy Python strings to Rust Vec<Vec<u8>> so the parallel loop does not
    // touch the GIL.
    let mut byte_rows: Vec<Vec<u8>> = Vec::with_capacity(n);
    for item in genotypes.iter() {
        let s: &str = item.extract()?;
        if s.len() != n_sites {
            return Err(PyValueError::new_err(format!(
                "genotype {:?} has length {}, expected {}",
                s,
                s.len(),
                n_sites
            )));
        }
        byte_rows.push(s.as_bytes().to_vec());
    }

    // Lift encoding lookups to Rust-native structures so the parallel loop
    // is GIL-free.
    let offsets_view = letter_to_bits_offsets.as_array();
    let offsets: Vec<i64> = offsets_view.to_vec();
    if offsets.len() != n_sites + 1 {
        return Err(PyValueError::new_err(format!(
            "letter_to_bits_offsets length {} != n_sites + 1 = {}",
            offsets.len(),
            n_sites + 1
        )));
    }
    let bits_view = letter_to_bits_flat.as_array();
    let bits_rows = bits_view.shape()[0];
    let bits_cols = bits_view.shape()[1];
    if bits_cols != n_bits {
        return Err(PyValueError::new_err(format!(
            "letter_to_bits_flat has {} cols, expected n_bits = {}",
            bits_cols, n_bits
        )));
    }
    let bits_flat: Vec<u8> = bits_view.iter().copied().collect();
    if bits_flat.len() != bits_rows * bits_cols {
        return Err(PyValueError::new_err(
            "letter_to_bits_flat is not contiguous or has unexpected shape",
        ));
    }

    // Per-site: Vec<(letter_byte -> local_row)>. We flatten each dict to a
    // 256-entry lookup array (i32, -1 = invalid) so the hot loop is O(1).
    let mut per_site_lookup: Vec<[i32; 256]> = vec![[-1i32; 256]; n_sites];
    if letter_index.len() != n_sites {
        return Err(PyValueError::new_err(format!(
            "letter_index length {} != n_sites = {}",
            letter_index.len(),
            n_sites
        )));
    }
    for (site_i, entry) in letter_index.iter().enumerate() {
        let d: std::collections::HashMap<String, i32> = entry.extract()?;
        for (letter, idx) in d {
            let b = letter.as_bytes();
            if b.len() != 1 {
                return Err(PyValueError::new_err(format!(
                    "letter {:?} at site {} is not a single ASCII byte",
                    letter, site_i
                )));
            }
            per_site_lookup[site_i][b[0] as usize] = idx;
        }
    }

    // Release the GIL and run the hot loop in parallel.
    let output: Vec<u8> = py.allow_threads(|| -> Result<Vec<u8>, String> {
        if n_bits == 0 {
            // Nothing to pack (all sites frozen). Still validate letters.
            for (i, row) in byte_rows.iter().enumerate() {
                for site_i in 0..n_sites {
                    let byte = row[site_i];
                    if per_site_lookup[site_i][byte as usize] < 0 {
                        return Err(format!(
                            "unknown letter {:?} at site {} (row {})",
                            byte as char, site_i, i
                        ));
                    }
                }
            }
            return Ok(Vec::new());
        }
        let mut out = vec![0u8; n * n_bits];
        out.par_chunks_mut(n_bits)
            .zip(byte_rows.par_iter())
            .try_for_each(|(dst, row)| -> Result<(), String> {
                for site_i in 0..n_sites {
                    let byte = row[site_i];
                    let local = per_site_lookup[site_i][byte as usize];
                    if local < 0 {
                        return Err(format!(
                            "unknown letter {:?} at site {}",
                            byte as char, site_i
                        ));
                    }
                    let letter_row = (offsets[site_i] as usize) + (local as usize);
                    let src = &bits_flat[letter_row * n_bits..(letter_row + 1) * n_bits];
                    // XOR into dst: WT letters are all-zero, mutant letters have a single 1 bit
                    // at their slot. Each site's bits occupy disjoint column ranges, so XOR == OR.
                    for k in 0..n_bits {
                        dst[k] |= src[k];
                    }
                }
                Ok(())
            })?;
        Ok(out)
    }).map_err(PyValueError::new_err)?;

    // Hand ownership back to numpy.
    let arr = ndarray::Array2::from_shape_vec((n, n_bits), output)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(arr.into_pyarray_bound(py))
}

/// Hamming distance from each row of `genotypes_int` to `reference`.
#[pyfunction]
fn hamming_to_reference<'py>(
    py: Python<'py>,
    genotypes_int: PyReadonlyArray2<'py, u8>,
    reference: numpy::PyReadonlyArray1<'py, u8>,
) -> PyResult<Bound<'py, numpy::PyArray1<u32>>> {
    hamming::hamming_to_reference(py, genotypes_int, reference)
}

/// Enumerate the Cartesian product of per-site alphabet indices as a 2D u8 array.
///
/// Refuses to run if the product of `alphabet_sizes` exceeds `max_genotypes`
/// unless `allow_huge` is True.
#[pyfunction]
#[pyo3(signature = (alphabet_sizes, max_genotypes=1usize << 28, allow_huge=false))]
fn enumerate_genotypes<'py>(
    py: Python<'py>,
    alphabet_sizes: Vec<usize>,
    max_genotypes: usize,
    allow_huge: bool,
) -> PyResult<Bound<'py, PyArray2<u8>>> {
    enumerate::enumerate_genotypes(py, alphabet_sizes, max_genotypes, allow_huge)
}

#[pymodule]
fn _rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(genotypes_to_binary_packed, m)?)?;
    m.add_function(wrap_pyfunction!(hamming_to_reference, m)?)?;
    m.add_function(wrap_pyfunction!(enumerate_genotypes, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}

// Silence dead_code on encoding mod until we wire more entry points.
#[allow(dead_code)]
fn _touch_encoding_mod() {
    encoding::_noop();
}
