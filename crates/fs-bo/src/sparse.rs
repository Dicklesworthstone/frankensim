//! Inducing-point sparse GPs (DTC/SoR predictive) with the TITSIAS
//! variational ELBO as the honesty instrument: the ELBO provably
//! lower-bounds the exact log marginal likelihood, is TIGHT when the
//! inducing set equals the data (Q_XX = K_XX), and its trace slack
//! (1/2σ²)·tr(K_XX − Q_XX) reports exactly how much the approximation
//! discards — a number, not a vibe. All linear algebra through fs-la
//! Cholesky; O(n·m²) fit against the exact GP's O(n³).

use crate::gp::Kernel;
use fs_la::factor::{Cholesky, cholesky};

/// A fitted sparse GP.
pub struct SparseGp {
    /// Kernel (shared with the exact machinery).
    pub kernel: Kernel,
    /// Observation-noise variance.
    pub noise: f64,
    /// Inducing locations.
    pub z: Vec<Vec<f64>>,
    /// Cholesky of K_ZZ (+ jitter).
    kzz_chol: Cholesky,
    /// Cholesky of A = K_ZZ + σ⁻²·K_ZX·K_XZ.
    a_chol: Cholesky,
    /// σ⁻²·A⁻¹·K_ZX·y (the predictive-mean weights).
    beta: Vec<f64>,
    /// The Titsias evidence lower bound.
    pub elbo: f64,
}

fn lsolve(chol: &Cholesky, v: &mut [f64]) {
    let n = v.len();
    for i in 0..n {
        let mut acc = v[i];
        for (j, vj) in v.iter().enumerate().take(i) {
            acc = (-chol.l(i, j)).mul_add(*vj, acc);
        }
        v[i] = acc / chol.l(i, i);
    }
}

impl SparseGp {
    /// Fit at fixed inducing locations `z`.
    ///
    /// # Panics
    /// If K_ZZ (jittered) or A lose positive-definiteness — degenerate
    /// inducing sets are a caller error at fixture scale.
    #[must_use]
    pub fn fit(
        x: &[Vec<f64>],
        y: &[f64],
        kernel: Kernel,
        noise: f64,
        z: Vec<Vec<f64>>,
    ) -> SparseGp {
        let n = x.len();
        let m = z.len();
        assert_eq!(n, y.len());
        assert!(m >= 1 && noise > 0.0);
        let s2inv = 1.0 / noise;
        // K_ZZ (+ jitter) and K_ZX.
        let mut kzz = vec![0.0f64; m * m];
        for i in 0..m {
            for j in 0..=i {
                let v = kernel.eval(&z[i], &z[j]);
                kzz[i * m + j] = v;
                kzz[j * m + i] = v;
            }
            kzz[i * m + i] += 1e-10;
        }
        let mut kzx = vec![0.0f64; m * n];
        for i in 0..m {
            for (j, xj) in x.iter().enumerate() {
                kzx[i * n + j] = kernel.eval(&z[i], xj);
            }
        }
        let kzz_chol = cholesky(&kzz, m).expect("K_ZZ must be PD (jittered)");
        // A = K_ZZ + σ⁻²·K_ZX·K_XZ.
        let mut a = kzz;
        for i in 0..m {
            for j in 0..=i {
                let mut acc = 0.0f64;
                for t in 0..n {
                    acc = kzx[i * n + t].mul_add(kzx[j * n + t], acc);
                }
                a[i * m + j] += s2inv * acc;
                a[j * m + i] = a[i * m + j];
            }
        }
        let a_chol = cholesky(&a, m).expect("A must be PD");
        // beta = σ⁻²·A⁻¹·K_ZX·y.
        let mut kzx_y = vec![0.0f64; m];
        for i in 0..m {
            let mut acc = 0.0f64;
            for (t, &yt) in y.iter().enumerate() {
                acc = kzx[i * n + t].mul_add(yt, acc);
            }
            kzx_y[i] = acc;
        }
        let mut beta = kzx_y.clone();
        a_chol.solve(&mut beta);
        for b in &mut beta {
            *b *= s2inv;
        }
        // ELBO = −n/2·log 2π − ½·log|Q+σ²I| − ½·yᵀ(Q+σ²I)⁻¹y
        //        − (1/2σ²)·tr(K_XX − Q_XX)
        // log|Q+σ²I| = log|A| − log|K_ZZ| + n·log σ².
        let logdet = |c: &Cholesky, k: usize| -> f64 {
            (0..k).map(|i| 2.0 * fs_math::det::ln(c.l(i, i))).sum()
        };
        let log_q = logdet(&a_chol, m) - logdet(&kzz_chol, m)
            + n as f64 * fs_math::det::ln(noise);
        // yᵀ(Q+σ²I)⁻¹y = σ⁻²·(yᵀy − σ⁻²·(K_ZX y)ᵀ A⁻¹ (K_ZX y)).
        let yty: f64 = y.iter().map(|v| v * v).sum();
        let mut a_inv_kzxy = kzx_y.clone();
        a_chol.solve(&mut a_inv_kzxy);
        let quad_inner: f64 = kzx_y.iter().zip(&a_inv_kzxy).map(|(p, q)| p * q).sum();
        let quad = s2inv * s2inv.mul_add(-quad_inner, yty);
        // tr(K_XX − Q_XX) = Σ k(xᵢ,xᵢ) − Σ ‖L_zz⁻¹·k_Z(xᵢ)‖².
        let mut trace_slack = 0.0f64;
        for (j, xj) in x.iter().enumerate() {
            trace_slack += kernel.eval(xj, xj);
            let mut col: Vec<f64> = (0..m).map(|i| kzx[i * n + j]).collect();
            lsolve(&kzz_chol, &mut col);
            trace_slack -= col.iter().map(|v| v * v).sum::<f64>();
        }
        trace_slack = trace_slack.max(0.0);
        let elbo = (-0.5f64).mul_add(
            quad,
            (-0.5f64).mul_add(
                log_q,
                (-0.5 * n as f64).mul_add(
                    fs_math::det::ln(2.0 * core::f64::consts::PI),
                    -(0.5 * s2inv * trace_slack),
                ),
            ),
        );
        SparseGp {
            kernel,
            noise,
            z,
            kzz_chol,
            a_chol,
            beta,
            elbo,
        }
    }

    /// DTC predictive mean and variance at a point.
    #[must_use]
    pub fn predict(&self, xs: &[f64]) -> (f64, f64) {
        let m = self.z.len();
        let kstar: Vec<f64> = self.z.iter().map(|zi| self.kernel.eval(zi, xs)).collect();
        let mean: f64 = kstar.iter().zip(&self.beta).map(|(a, b)| a * b).sum();
        // Q_** = ‖L_zz⁻¹·k_*‖²; explained = ‖L_A⁻¹·k_*‖².
        let mut v1 = kstar.clone();
        lsolve(&self.kzz_chol, &mut v1);
        let q_ss: f64 = v1.iter().map(|t| t * t).sum();
        let mut v2 = kstar;
        lsolve(&self.a_chol, &mut v2);
        let expl: f64 = v2.iter().map(|t| t * t).sum();
        let kss = self.kernel.eval(xs, xs);
        let var = (kss - q_ss + expl).max(0.0);
        let _ = m;
        (mean, var)
    }
}

/// Deterministic farthest-point inducing selection: start at index 0,
/// greedily add the point maximizing the minimum distance to the
/// chosen set (index tie-break).
#[must_use]
pub fn farthest_point_inducing(x: &[Vec<f64>], m: usize) -> Vec<Vec<f64>> {
    let n = x.len();
    assert!(m >= 1 && m <= n);
    let mut chosen = vec![0usize];
    let d2 = |a: &[f64], b: &[f64]| -> f64 {
        a.iter().zip(b).map(|(p, q)| (p - q) * (p - q)).sum()
    };
    let mut min_d: Vec<f64> = x.iter().map(|xi| d2(xi, &x[0])).collect();
    while chosen.len() < m {
        let mut best = (0usize, f64::NEG_INFINITY);
        for (i, &d) in min_d.iter().enumerate() {
            if d > best.1 {
                best = (i, d);
            }
        }
        chosen.push(best.0);
        for (i, md) in min_d.iter_mut().enumerate() {
            *md = md.min(d2(&x[i], &x[best.0]));
        }
    }
    chosen.into_iter().map(|i| x[i].clone()).collect()
}
