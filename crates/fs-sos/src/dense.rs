//! Small dense kernels for the SDP solver and the Lyapunov solve. Row-major
//! `n × n` storage in a flat `Vec<f64>`. Sizes here are Gram/Schur matrices of
//! low-dimensional polynomial programs (tens to a few hundred rows), so plain
//! cache-friendly loops are the honest choice; nothing here claims a roofline.

// Dense index kernels: explicit `i, j, k` loops are the readable form.
#![allow(clippy::needless_range_loop)]

/// Cholesky factor `L` (lower, row-major) of a symmetric positive definite
/// matrix, or `None` if a pivot is not strictly positive and finite.
pub(crate) fn cholesky(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0; n * n];
    for j in 0..n {
        let mut s = a[j * n + j];
        for k in 0..j {
            s -= l[j * n + k] * l[j * n + k];
        }
        if !(s > 0.0) || !s.is_finite() {
            return None;
        }
        let d = s.sqrt();
        l[j * n + j] = d;
        for i in (j + 1)..n {
            let mut t = a[i * n + j];
            for k in 0..j {
                t -= l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = t / d;
        }
    }
    Some(l)
}

/// Solve `L Lᵀ x = b` in place.
pub(crate) fn chol_solve(l: &[f64], n: usize, b: &mut [f64]) {
    for i in 0..n {
        let mut s = b[i];
        for k in 0..i {
            s -= l[i * n + k] * b[k];
        }
        b[i] = s / l[i * n + i];
    }
    for i in (0..n).rev() {
        let mut s = b[i];
        for k in (i + 1)..n {
            s -= l[k * n + i] * b[k];
        }
        b[i] = s / l[i * n + i];
    }
}

/// Inverse of `L Lᵀ` (symmetric), column by column.
pub(crate) fn chol_inverse(l: &[f64], n: usize) -> Vec<f64> {
    let mut inv = vec![0.0; n * n];
    let mut col = vec![0.0; n];
    for j in 0..n {
        col.fill(0.0);
        col[j] = 1.0;
        chol_solve(l, n, &mut col);
        for i in 0..n {
            inv[i * n + j] = col[i];
        }
    }
    // Exact symmetry: average the two triangles.
    for i in 0..n {
        for j in (i + 1)..n {
            let m = 0.5 * (inv[i * n + j] + inv[j * n + i]);
            inv[i * n + j] = m;
            inv[j * n + i] = m;
        }
    }
    inv
}

/// `A · B` for square row-major matrices.
pub(crate) fn matmul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i * n + j] += aik * b[k * n + j];
            }
        }
    }
    c
}

/// `(A + Aᵀ)/2` in place.
pub(crate) fn symmetrize(a: &mut [f64], n: usize) {
    for i in 0..n {
        for j in (i + 1)..n {
            let m = 0.5 * (a[i * n + j] + a[j * n + i]);
            a[i * n + j] = m;
            a[j * n + i] = m;
        }
    }
}

/// `L⁻¹ S L⁻ᵀ` for symmetric `S` and lower-triangular `L`.
pub(crate) fn congruence_inv(l: &[f64], s: &[f64], n: usize) -> Vec<f64> {
    // W = L⁻¹ S (forward substitution per column).
    let mut w = s.to_vec();
    for j in 0..n {
        for i in 0..n {
            let mut t = w[i * n + j];
            for k in 0..i {
                t -= l[i * n + k] * w[k * n + j];
            }
            w[i * n + j] = t / l[i * n + i];
        }
    }
    // T = L⁻¹ Wᵀ, then T is symmetric (= L⁻¹ S L⁻ᵀ).
    let mut t = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let mut v = w[j * n + i];
            for k in 0..i {
                v -= l[i * n + k] * t[k * n + j];
            }
            t[i * n + j] = v / l[i * n + i];
        }
    }
    symmetrize(&mut t, n);
    t
}

/// All eigenvalues of a symmetric matrix by cyclic Jacobi (ascending).
pub(crate) fn sym_eigenvalues(a: &[f64], n: usize) -> Vec<f64> {
    sym_eigen(a, n, false).0
}

/// Eigen-decomposition of a symmetric matrix by cyclic Jacobi: ascending
/// eigenvalues and (when requested) the matching orthonormal eigenvectors as
/// COLUMNS of a row-major matrix.
pub(crate) fn sym_eigen(a: &[f64], n: usize, vectors: bool) -> (Vec<f64>, Vec<f64>) {
    let mut m = a.to_vec();
    symmetrize(&mut m, n);
    let mut v = if vectors {
        let mut id = vec![0.0; n * n];
        for i in 0..n {
            id[i * n + i] = 1.0;
        }
        id
    } else {
        Vec::new()
    };
    let scale = m
        .iter()
        .fold(0.0f64, |s, x| s.max(x.abs()))
        .max(f64::MIN_POSITIVE);
    for _sweep in 0..60 {
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += m[i * n + j] * m[i * n + j];
            }
        }
        if off.sqrt() <= 1e-15 * scale {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = m[p * n + q];
                if apq.abs() <= 1e-300 {
                    continue;
                }
                let theta = (m[q * n + q] - m[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (mkp, mkq) = (m[k * n + p], m[k * n + q]);
                    m[k * n + p] = c * mkp - s * mkq;
                    m[k * n + q] = s * mkp + c * mkq;
                }
                for k in 0..n {
                    let (mpk, mqk) = (m[p * n + k], m[q * n + k]);
                    m[p * n + k] = c * mpk - s * mqk;
                    m[q * n + k] = s * mpk + c * mqk;
                }
                if vectors {
                    for k in 0..n {
                        let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                        v[k * n + p] = c * vkp - s * vkq;
                        v[k * n + q] = s * vkp + c * vkq;
                    }
                }
            }
        }
    }
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| m[i * n + i].total_cmp(&m[j * n + j]));
    let vals: Vec<f64> = idx.iter().map(|&i| m[i * n + i]).collect();
    let vecs = if vectors {
        let mut out = vec![0.0; n * n];
        for (newc, &oldc) in idx.iter().enumerate() {
            for r in 0..n {
                out[r * n + newc] = v[r * n + oldc];
            }
        }
        out
    } else {
        Vec::new()
    };
    (vals, vecs)
}

/// Solve the general square system `A x = b` by Gaussian elimination with
/// partial pivoting; `None` if singular to working precision.
pub(crate) fn lu_solve(a: &[f64], n: usize, b: &[f64]) -> Option<Vec<f64>> {
    let mut m = a.to_vec();
    let mut x = b.to_vec();
    let scale = m.iter().fold(0.0f64, |s, v| s.max(v.abs()));
    if !(scale > 0.0) {
        return None;
    }
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if m[r * n + col].abs() > m[piv * n + col].abs() {
                piv = r;
            }
        }
        if m[piv * n + col].abs() <= 1e-14 * scale {
            return None;
        }
        if piv != col {
            for k in 0..n {
                m.swap(col * n + k, piv * n + k);
            }
            x.swap(col, piv);
        }
        let d = m[col * n + col];
        for r in (col + 1)..n {
            let f = m[r * n + col] / d;
            if f == 0.0 {
                continue;
            }
            for k in col..n {
                m[r * n + k] -= f * m[col * n + k];
            }
            x[r] -= f * x[col];
        }
    }
    for r in (0..n).rev() {
        let mut s = x[r];
        for k in (r + 1)..n {
            s -= m[r * n + k] * x[k];
        }
        x[r] = s / m[r * n + r];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cholesky_inverse_round_trip() {
        let a = vec![4.0, 2.0, 0.4, 2.0, 5.0, 1.0, 0.4, 1.0, 3.0];
        let l = cholesky(&a, 3).expect("spd");
        let inv = chol_inverse(&l, 3);
        let id = matmul(&a, &inv, 3);
        for i in 0..3 {
            for j in 0..3 {
                let e = if i == j { 1.0 } else { 0.0 };
                assert!((id[i * 3 + j] - e).abs() < 1e-13);
            }
        }
        assert!(cholesky(&[1.0, 2.0, 2.0, 1.0], 2).is_none());
    }

    #[test]
    fn jacobi_eigenvalues_and_vectors() {
        let a = vec![2.0, 1.0, 1.0, 2.0];
        let (vals, vecs) = sym_eigen(&a, 2, true);
        assert!((vals[0] - 1.0).abs() < 1e-14 && (vals[1] - 3.0).abs() < 1e-14);
        // A v = λ v for the top eigenvector.
        let (v0, v1) = (vecs[1], vecs[3]);
        assert!((2.0 * v0 + v1 - 3.0 * v0).abs() < 1e-13);
    }

    #[test]
    fn congruence_matches_definition() {
        let x = vec![4.0, 1.0, 1.0, 3.0];
        let l = cholesky(&x, 2).unwrap();
        let s = vec![1.0, -2.0, -2.0, 0.5];
        let t = congruence_inv(&l, &s, 2);
        // L T Lᵀ == S
        let lt = matmul(&l, &t, 2);
        let ltr = vec![l[0], l[2], l[1], l[3]];
        let back = matmul(&lt, &ltr, 2);
        for k in 0..4 {
            assert!((back[k] - s[k]).abs() < 1e-13);
        }
    }

    #[test]
    fn lu_solves_and_detects_singularity() {
        let a = vec![0.0, 2.0, 1.0, 1.0];
        let x = lu_solve(&a, 2, &[2.0, 3.0]).unwrap();
        assert!((x[0] - 2.0).abs() < 1e-14 && (x[1] - 1.0).abs() < 1e-14);
        assert!(lu_solve(&[1.0, 2.0, 2.0, 4.0], 2, &[1.0, 2.0]).is_none());
    }
}
