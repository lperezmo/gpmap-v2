use numpy::{IntoPyArray, PyArray1, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rayon::prelude::*;

pub fn hamming_to_reference<'py>(
    py: Python<'py>,
    genotypes_int: PyReadonlyArray2<'py, u8>,
    reference: PyReadonlyArray1<'py, u8>,
) -> PyResult<Bound<'py, PyArray1<u32>>> {
    let g = genotypes_int.as_array();
    let r = reference.as_array();
    let n_sites = g.shape()[1];
    if r.len() != n_sites {
        return Err(PyValueError::new_err(format!(
            "reference length {} does not match genotypes columns {}",
            r.len(),
            n_sites
        )));
    }
    let n = g.shape()[0];
    let g_flat: Vec<u8> = g.iter().copied().collect();
    let r_vec: Vec<u8> = r.iter().copied().collect();

    let out: Vec<u32> = py.allow_threads(|| {
        let mut out = vec![0u32; n];
        out.par_iter_mut().enumerate().for_each(|(i, slot)| {
            let row = &g_flat[i * n_sites..(i + 1) * n_sites];
            let mut d = 0u32;
            for j in 0..n_sites {
                if row[j] != r_vec[j] {
                    d += 1;
                }
            }
            *slot = d;
        });
        out
    });
    Ok(ndarray::Array1::from_vec(out).into_pyarray_bound(py))
}
