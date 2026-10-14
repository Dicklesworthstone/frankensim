//! Small dense symmetric kernels for the COLD voice-set change.
//!
//! Why this exists: the offline reference uses `fs_modal::eigh_gen_dense` for
//! every small pencil, and so does this crate at construction. Inside the
//! browser, however, that routine (and everything built on `f64::mul_add`)
//! runs on a software fused multiply-add: measured 1.79 ms per 25x25 pencil
//! in wasm against 0.085 ms natively, and a voice-set change needs four of
//! them. The routines below use plain multiply/add only. They compute the
//! SAME mathematical objects (generalized symmetric eigen-decomposition, SPD
//! inverse, matrix exponential of a small symmetric matrix); the native
//! validation compares the resulting motion against the reference engine.
//!
//! Row-major `n x n` slices throughout.

/// Lower Cholesky factor `L` of a symmetric positive definite matrix.
pub fn cholesky(a: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let mut l = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i * n + j];
            for k in 0..j {
                sum -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if !(sum > 0.0) || !sum.is_finite() {
                    return Err("matrix is not positive definite".into());
                }
                l[i * n + i] = sum.sqrt();
            } else {
                l[i * n + j] = sum / l[j * n + j];
            }
        }
    }
    Ok(l)
}

/// Inverse of a symmetric positive definite matrix through its Cholesky factor.
pub fn spd_inverse(a: &[f64], n: usize, what: &str) -> Result<Vec<f64>, String> {
    let l = cholesky(a, n).map_err(|_| format!("{what} is not positive definite"))?;
    // W = L^-1 (lower), then A^-1 = W^T W.
    let mut w = vec![0.0; n * n];
    for j in 0..n {
        w[j * n + j] = 1.0 / l[j * n + j];
        for i in j + 1..n {
            let mut sum = 0.0;
            for k in j..i {
                sum -= l[i * n + k] * w[k * n + j];
            }
            w[i * n + j] = sum / l[i * n + i];
        }
    }
    let mut inverse = vec![0.0; n * n];
    for i in 0..n {
        for j in i..n {
            let mut sum = 0.0;
            for k in j..n {
                sum += w[k * n + i] * w[k * n + j];
            }
            inverse[i * n + j] = sum;
            inverse[j * n + i] = sum;
        }
    }
    if inverse.iter().any(|x| !x.is_finite()) {
        return Err(format!("{what} inverse overflow"));
    }
    Ok(inverse)
}

/// Cyclic Jacobi eigen-decomposition of a symmetric matrix (destroyed).
/// Returns eigenvalues (unsorted) and the eigenvector matrix (row-major: row `p` is eigenvector `p`).
fn jacobi(a: &mut [f64], n: usize) -> Result<(Vec<f64>, Vec<f64>), String> {
    let mut vt = vec![0.0; n * n];
    for i in 0..n {
        vt[i * n + i] = 1.0;
    }
    for _sweep in 0..60 {
        let mut rotated = false;
        for p in 0..n.saturating_sub(1) {
            for q in p + 1..n {
                let apq = a[p * n + q];
                let (app, aqq) = (a[p * n + p], a[q * n + q]);
                if apq == 0.0 || apq.abs() <= 1.0e-17 * (app.abs() + aqq.abs()) {
                    continue;
                }
                rotated = true;
                let theta = (aqq - app) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                a[p * n + q] = 0.0;
                a[q * n + p] = 0.0;
                let (head, tail) = vt.split_at_mut(q * n);
                let vp = &mut head[p * n..(p + 1) * n];
                let vq = &mut tail[..n];
                for (vkp, vkq) in vp.iter_mut().zip(vq.iter_mut()) {
                    let (old_p, old_q) = (*vkp, *vkq);
                    *vkp = c * old_p - s * old_q;
                    *vkq = s * old_p + c * old_q;
                }
            }
        }
        if !rotated {
            let values = (0..n).map(|i| a[i * n + i]).collect();
            return Ok((values, vt));
        }
    }
    Err("symmetric eigen-solve did not converge".into())
}

/// Generalized symmetric-definite eigenproblem `K x = lambda M x`.
/// Returns ascending `lambda` and `phi[i*n+j]` = component `i` of mode `j`,
/// mass-normalized: `Phi^T M Phi = I`.
pub fn eigh_gen(k: &[f64], m: &[f64], n: usize) -> Result<(Vec<f64>, Vec<f64>), String> {
    let l = cholesky(m, n).map_err(|_| "loaded board mass is not positive definite".to_string())?;
    // X = L^-1 K, then C = X L^-T = (L^-1 X^T)^T.
    let mut x = vec![0.0; n * n];
    for col in 0..n {
        for i in 0..n {
            let mut sum = k[i * n + col];
            for j in 0..i {
                sum -= l[i * n + j] * x[j * n + col];
            }
            x[i * n + col] = sum / l[i * n + i];
        }
    }
    let mut c = vec![0.0; n * n];
    for row in 0..n {
        // Solve L y = (row `row` of X)^T; y is row `row` of C.
        for i in 0..n {
            let mut sum = x[row * n + i];
            for j in 0..i {
                sum -= l[i * n + j] * c[row * n + j];
            }
            c[row * n + i] = sum / l[i * n + i];
        }
    }
    for i in 0..n {
        for j in i + 1..n {
            let mean = 0.5 * (c[i * n + j] + c[j * n + i]);
            c[i * n + j] = mean;
            c[j * n + i] = mean;
        }
    }
    let (values, yt) = jacobi(&mut c, n)?;
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|a, b| values[*a].total_cmp(&values[*b]));
    let mut phi = vec![0.0; n * n];
    let mut lambda = vec![0.0; n];
    for (j, &source) in order.iter().enumerate() {
        lambda[j] = values[source];
        let y_row = &yt[source * n..(source + 1) * n];
        // Back-substitute L^T phi = y.
        for i in (0..n).rev() {
            let mut sum = y_row[i];
            for kk in i + 1..n {
                sum -= l[kk * n + i] * phi[kk * n + j];
            }
            phi[i * n + j] = sum / l[i * n + i];
        }
    }
    if lambda.iter().chain(&phi).any(|v| !v.is_finite()) {
        return Err("loaded board eigen-solve overflow".into());
    }
    Ok((lambda, phi))
}

/// `exp(-X)` of a small symmetric matrix by its Taylor series. Refuses a
/// matrix whose norm is too large for the declared term budget.
pub fn exp_neg_small(x: &[f64], n: usize) -> Result<Vec<f64>, String> {
    let norm = (0..n).map(|i| x[i * n..(i + 1) * n].iter().map(|v| v.abs()).sum::<f64>()).fold(0.0_f64, f64::max);
    if !norm.is_finite() || norm > 0.25 {
        return Err("board damping per half step is too large for the series exponential".into());
    }
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        result[i * n + i] = 1.0;
    }
    let mut term = result.clone();
    let mut next = vec![0.0; n * n];
    for order in 1..=24 {
        let inv_order = -1.0 / order as f64;
        // Both `term` and `x` are symmetric: (term * (-X) / order)[i, j]
        // uses two contiguous rows `term[i, ..]` and `x[j, ..]` for j >= i.
        for i in 0..n {
            let row_i = &term[i * n..(i + 1) * n];
            for j in i..n {
                let row_j = &x[j * n..(j + 1) * n];
                let mut sum = 0.0;
                for (&tik, &xjk) in row_i.iter().zip(row_j) {
                    sum += tik * xjk;
                }
                let val = sum * inv_order;
                next[i * n + j] = val;
                next[j * n + i] = val;
            }
        }
        std::mem::swap(&mut term, &mut next);
        let size = term.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        for (r, t) in result.iter_mut().zip(&term) {
            *r += t;
        }
        if size < 1.0e-19 {
            break;
        }
    }
    for i in 0..n {
        for j in i + 1..n {
            let mean = 0.5 * (result[i * n + j] + result[j * n + i]);
            result[i * n + j] = mean;
            result[j * n + i] = mean;
        }
    }
    Ok(result)
}
