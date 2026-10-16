//! Small dense kernels shared by the Koopman/DMD and DEIM reduced-order
//! models. Row-major flat storage. Sizes are reduced-order (tens to a few
//! hundred), so plain loops are the honest choice; no roofline is claimed.

// Dense index kernels: explicit `i, j, k` loops are the readable form.
#![allow(clippy::needless_range_loop)]

/// Gaussian elimination with partial pivoting: solve `A x = b` for one RHS.
/// `None` when `A` is singular to working precision.
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

/// Inverse of a square matrix by repeated [`lu_solve`].
pub(crate) fn inverse(a: &[f64], n: usize) -> Option<Vec<f64>> {
    let mut inv = vec![0.0; n * n];
    let mut e = vec![0.0; n];
    for j in 0..n {
        e.fill(0.0);
        e[j] = 1.0;
        let col = lu_solve(a, n, &e)?;
        for i in 0..n {
            inv[i * n + j] = col[i];
        }
    }
    Some(inv)
}

/// Cholesky solve of `(G + ridge·I) x = b` for SPD-ish `G`, escalating the
/// ridge by ×10 until the factorization succeeds (at most 12 times).
pub(crate) fn ridge_solve(g: &[f64], n: usize, b: &[f64], ridge: f64) -> Option<Vec<f64>> {
    let dmax = (0..n)
        .fold(0.0f64, |s, i| s.max(g[i * n + i].abs()))
        .max(1e-300);
    let mut lam = ridge.max(0.0);
    for _ in 0..12 {
        let mut a = g.to_vec();
        for i in 0..n {
            a[i * n + i] += lam * dmax;
        }
        if let Some(l) = cholesky(&a, n) {
            let mut x = b.to_vec();
            for i in 0..n {
                let mut s = x[i];
                for k in 0..i {
                    s -= l[i * n + k] * x[k];
                }
                x[i] = s / l[i * n + i];
            }
            for i in (0..n).rev() {
                let mut s = x[i];
                for k in (i + 1)..n {
                    s -= l[k * n + i] * x[k];
                }
                x[i] = s / l[i * n + i];
            }
            return Some(x);
        }
        lam = if lam == 0.0 { 1e-14 } else { lam * 10.0 };
    }
    None
}

fn cholesky(a: &[f64], n: usize) -> Option<Vec<f64>> {
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

/// Eigenvalues `(re, im)` of a real square matrix: Householder reduction to
/// upper Hessenberg form, then the Francis double-shift QR iteration (the
/// classical `hqr` scheme of Wilkinson–Reinsch / Numerical Recipes §11.6,
/// with exceptional shifts at iterations 10 and 20). Complex eigenvalues
/// come in conjugate pairs. `None` if an eigenvalue fails to converge in 60
/// iterations.
///
/// The negligibility tests `|h| + s == s` are deliberate exact float
/// comparisons — the standard deflation criterion.
#[allow(
    clippy::float_cmp,
    clippy::too_many_lines,
    clippy::many_single_char_names
)]
pub(crate) fn eigenvalues_real(a_in: &[f64], n: usize) -> Option<Vec<(f64, f64)>> {
    if n == 0 {
        return Some(Vec::new());
    }
    let mut a = a_in.to_vec();
    // Householder → upper Hessenberg.
    for k in 0..n.saturating_sub(2) {
        let len = n - k - 1;
        let mut v: Vec<f64> = (0..len).map(|i| a[(k + 1 + i) * n + k]).collect();
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if v[0] > 0.0 { -norm } else { norm };
        v[0] -= alpha;
        let vn = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if vn == 0.0 {
            continue;
        }
        for x in &mut v {
            *x /= vn;
        }
        for j in 0..n {
            let s: f64 = (0..len).map(|i| v[i] * a[(k + 1 + i) * n + j]).sum();
            for i in 0..len {
                a[(k + 1 + i) * n + j] -= 2.0 * v[i] * s;
            }
        }
        for i in 0..n {
            let s: f64 = (0..len).map(|j| a[i * n + k + 1 + j] * v[j]).sum();
            for j in 0..len {
                a[i * n + k + 1 + j] -= 2.0 * s * v[j];
            }
        }
    }
    // 1-based accessor over the Hessenberg matrix.
    let ix = |i: isize, j: isize| -> usize { ((i - 1) as usize) * n + (j - 1) as usize };
    let ni = n as isize;
    let mut wr = vec![0.0; n + 1];
    let mut wi = vec![0.0; n + 1];
    let mut anorm = 0.0;
    for i in 1..=ni {
        for j in (i - 1).max(1)..=ni {
            anorm += a[ix(i, j)].abs();
        }
    }
    let mut nn = ni;
    let mut t = 0.0;
    while nn >= 1 {
        let mut its = 0;
        loop {
            let mut l = nn;
            while l >= 2 {
                let mut s = a[ix(l - 1, l - 1)].abs() + a[ix(l, l)].abs();
                if s == 0.0 {
                    s = anorm;
                }
                if a[ix(l, l - 1)].abs() + s == s {
                    let k = ix(l, l - 1);
                    a[k] = 0.0;
                    break;
                }
                l -= 1;
            }
            let mut x = a[ix(nn, nn)];
            if l == nn {
                wr[nn as usize] = x + t;
                wi[nn as usize] = 0.0;
                nn -= 1;
            } else {
                let mut y = a[ix(nn - 1, nn - 1)];
                let mut w = a[ix(nn, nn - 1)] * a[ix(nn - 1, nn)];
                if l == nn - 1 {
                    let p = 0.5 * (y - x);
                    let q = p * p + w;
                    let mut z = q.abs().sqrt();
                    x += t;
                    if q >= 0.0 {
                        z = p + if p >= 0.0 { z } else { -z };
                        wr[(nn - 1) as usize] = x + z;
                        wr[nn as usize] = x + z;
                        if z != 0.0 {
                            wr[nn as usize] = x - w / z;
                        }
                        wi[(nn - 1) as usize] = 0.0;
                        wi[nn as usize] = 0.0;
                    } else {
                        wr[(nn - 1) as usize] = x + p;
                        wr[nn as usize] = x + p;
                        wi[(nn - 1) as usize] = -z;
                        wi[nn as usize] = z;
                    }
                    nn -= 2;
                } else {
                    if its == 60 {
                        return None;
                    }
                    if its == 10 || its == 20 {
                        t += x;
                        for i in 1..=nn {
                            let k = ix(i, i);
                            a[k] -= x;
                        }
                        let s = a[ix(nn, nn - 1)].abs() + a[ix(nn - 1, nn - 2)].abs();
                        x = 0.75 * s;
                        y = x;
                        w = -0.4375 * s * s;
                    }
                    its += 1;
                    let mut m = nn - 2;
                    let (mut p, mut q, mut r);
                    loop {
                        let z = a[ix(m, m)];
                        let rr = x - z;
                        let ss = y - z;
                        p = (rr * ss - w) / a[ix(m + 1, m)] + a[ix(m, m + 1)];
                        q = a[ix(m + 1, m + 1)] - z - rr - ss;
                        r = a[ix(m + 2, m + 1)];
                        let s = p.abs() + q.abs() + r.abs();
                        p /= s;
                        q /= s;
                        r /= s;
                        if m == l {
                            break;
                        }
                        let u = a[ix(m, m - 1)].abs() * (q.abs() + r.abs());
                        let v = p.abs()
                            * (a[ix(m - 1, m - 1)].abs() + z.abs() + a[ix(m + 1, m + 1)].abs());
                        if u + v == v {
                            break;
                        }
                        m -= 1;
                    }
                    for i in (m + 2)..=nn {
                        let k = ix(i, i - 2);
                        a[k] = 0.0;
                        if i != m + 2 {
                            let k = ix(i, i - 3);
                            a[k] = 0.0;
                        }
                    }
                    let mut k = m;
                    while k < nn {
                        if k != m {
                            p = a[ix(k, k - 1)];
                            q = a[ix(k + 1, k - 1)];
                            r = 0.0;
                            if k != nn - 1 {
                                r = a[ix(k + 2, k - 1)];
                            }
                            x = p.abs() + q.abs() + r.abs();
                            if x != 0.0 {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }
                        let s0 = (p * p + q * q + r * r).sqrt();
                        let s = if p >= 0.0 { s0 } else { -s0 };
                        if s != 0.0 {
                            if k == m {
                                if l != m {
                                    let kk = ix(k, k - 1);
                                    a[kk] = -a[kk];
                                }
                            } else {
                                let kk = ix(k, k - 1);
                                a[kk] = -s * x;
                            }
                            p += s;
                            x = p / s;
                            y = q / s;
                            let z = r / s;
                            q /= p;
                            r /= p;
                            for j in k..=nn {
                                let mut pp = a[ix(k, j)] + q * a[ix(k + 1, j)];
                                if k != nn - 1 {
                                    pp += r * a[ix(k + 2, j)];
                                    let kk = ix(k + 2, j);
                                    a[kk] -= pp * z;
                                }
                                let kk = ix(k + 1, j);
                                a[kk] -= pp * y;
                                let kk = ix(k, j);
                                a[kk] -= pp * x;
                            }
                            let mmin = if nn < k + 3 { nn } else { k + 3 };
                            for i in l..=mmin {
                                let mut pp = x * a[ix(i, k)] + y * a[ix(i, k + 1)];
                                if k != nn - 1 {
                                    pp += z * a[ix(i, k + 2)];
                                    let kk = ix(i, k + 2);
                                    a[kk] -= pp * r;
                                }
                                let kk = ix(i, k + 1);
                                a[kk] -= pp * q;
                                let kk = ix(i, k);
                                a[kk] -= pp;
                            }
                        }
                        k += 1;
                    }
                }
            }
            if l >= nn - 1 {
                break;
            }
        }
    }
    let mut out: Vec<(f64, f64)> = (1..=n).map(|i| (wr[i], wi[i])).collect();
    // Deterministic order: descending modulus, then by real, then imaginary.
    out.sort_by(|a, b| {
        let ma = a.0.hypot(a.1);
        let mb = b.0.hypot(b.1);
        mb.total_cmp(&ma)
            .then_with(|| b.0.total_cmp(&a.0))
            .then_with(|| b.1.total_cmp(&a.1))
    });
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hqr_recovers_real_and_complex_spectra() {
        // Block diag(rotation-scaling, 0.5) conjugated by a dense matrix.
        let (r, th) = (0.9f64, 0.3f64);
        let d = [
            r * th.cos(),
            -r * th.sin(),
            0.0,
            r * th.sin(),
            r * th.cos(),
            0.0,
            0.0,
            0.0,
            0.5,
        ];
        let s = [2.0, 1.0, 0.0, 0.5, 1.0, 1.0, 0.0, 0.3, 1.0];
        let sinv = inverse(&s, 3).unwrap();
        let mut sd = [0.0; 9];
        let mut a = vec![0.0; 9];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    sd[i * 3 + j] += s[i * 3 + k] * d[k * 3 + j];
                }
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    a[i * 3 + j] += sd[i * 3 + k] * sinv[k * 3 + j];
                }
            }
        }
        let ev = eigenvalues_real(&a, 3).unwrap();
        assert!((ev[0].0 - r * th.cos()).abs() < 1e-12);
        assert!((ev[0].1.abs() - r * th.sin()).abs() < 1e-12);
        assert!((ev[0].1 + ev[1].1).abs() < 1e-12, "conjugate pair");
        assert!((ev[2].0 - 0.5).abs() < 1e-12 && ev[2].1 == 0.0);
    }

    #[test]
    fn solvers_round_trip() {
        let a = vec![4.0, 1.0, 1.0, 3.0];
        let x = lu_solve(&a, 2, &[1.0, 2.0]).unwrap();
        assert!((4.0 * x[0] + x[1] - 1.0).abs() < 1e-14);
        let y = ridge_solve(&a, 2, &[1.0, 2.0], 0.0).unwrap();
        assert!((x[0] - y[0]).abs() < 1e-14 && (x[1] - y[1]).abs() < 1e-14);
        assert!(lu_solve(&[1.0, 2.0, 2.0, 4.0], 2, &[1.0, 1.0]).is_none());
    }
}
