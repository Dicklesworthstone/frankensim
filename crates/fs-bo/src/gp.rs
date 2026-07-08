//! Gaussian-process regression: Matérn-family kernels (½, 3⁄2, 5⁄2)
//! with ARD lengthscales, EXACT inference through fs-la Cholesky, log
//! marginal likelihood, and hyperparameter fitting by fs-ascent
//! L-BFGS with QMC-SEEDED MULTISTART (scrambled Sobol over the
//! log-parameter box — the plan names this detail). All math through
//! the strict fs-math kernels: posteriors are bit-deterministic.

use fs_la::factor::{Cholesky, cholesky};

/// Matérn smoothness family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Matern {
    /// ν = ½ (exponential; rough).
    Half,
    /// ν = 3⁄2 (once-differentiable).
    ThreeHalves,
    /// ν = 5⁄2 (twice-differentiable; the BO default).
    FiveHalves,
}

/// Kernel: signal variance σ² and per-dimension (ARD) lengthscales.
#[derive(Debug, Clone)]
pub struct Kernel {
    /// Smoothness family.
    pub family: Matern,
    /// Signal variance σ².
    pub signal: f64,
    /// ARD lengthscales (one per input dimension).
    pub lengthscales: Vec<f64>,
}

impl Kernel {
    /// Scaled distance r = √Σ((xᵢ−yᵢ)/ℓᵢ)².
    fn scaled_dist(&self, x: &[f64], y: &[f64]) -> f64 {
        let mut acc = 0.0f64;
        for ((xi, yi), l) in x.iter().zip(y).zip(&self.lengthscales) {
            let d = (xi - yi) / l;
            acc = d.mul_add(d, acc);
        }
        fs_math::det::sqrt(acc)
    }

    /// k(x, y).
    #[must_use]
    pub fn eval(&self, x: &[f64], y: &[f64]) -> f64 {
        let r = self.scaled_dist(x, y);
        // Guard the r → ∞ limit (degenerate lengthscales during
        // hyperparameter search): the Matérn polynomial×exp forms hit
        // inf·0 = NaN there, but the true limit is 0 — and a NaN here
        // poisons the Cholesky as a fake NotSpd (the first draft
        // crashed exactly this way mid-L-BFGS).
        if !r.is_finite() || r > 1e8 {
            return 0.0;
        }
        let core = match self.family {
            Matern::Half => fs_math::det::exp(-r),
            Matern::ThreeHalves => {
                let a = fs_math::det::sqrt(3.0) * r;
                (1.0 + a) * fs_math::det::exp(-a)
            }
            Matern::FiveHalves => {
                let a = fs_math::det::sqrt(5.0) * r;
                (1.0 + a + a * a / 3.0) * fs_math::det::exp(-a)
            }
        };
        self.signal * core
    }
}

/// A fitted exact GP (zero prior mean).
pub struct Gp {
    /// Kernel used.
    pub kernel: Kernel,
    /// Observation-noise variance σ_n².
    pub noise: f64,
    x: Vec<Vec<f64>>,
    /// α = (K + σ_n²I)⁻¹·y.
    alpha: Vec<f64>,
    chol: Cholesky,
    /// Log marginal likelihood of the fit.
    pub lml: f64,
}

impl Gp {
    /// Exact fit: Cholesky of K + σ_n²I.
    ///
    /// # Panics
    /// If the kernel matrix is not SPD (pathological hyperparameters —
    /// a modeling error when calling directly; the hyperparameter
    /// SEARCH uses [`Gp::try_fit`], where near-duplicate exploitation
    /// points with tiny noise make NotSpd a rejectable candidate, not
    /// a crash).
    #[must_use]
    pub fn fit(x: &[Vec<f64>], y: &[f64], kernel: Kernel, noise: f64) -> Gp {
        Gp::try_fit(x, y, kernel, noise).expect("kernel matrix must be SPD")
    }

    /// Fallible exact fit: `None` when K + σ_n²I is not SPD.
    #[must_use]
    pub fn try_fit(x: &[Vec<f64>], y: &[f64], kernel: Kernel, noise: f64) -> Option<Gp> {
        let n = x.len();
        assert_eq!(n, y.len());
        let mut k = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = kernel.eval(&x[i], &x[j]);
                k[i * n + j] = v;
                k[j * n + i] = v;
            }
            k[i * n + i] += noise;
        }
        let chol = cholesky(&k, n).ok()?;
        let mut alpha = y.to_vec();
        chol.solve(&mut alpha);
        // LML = −½yᵀα − Σᵢ log Lᵢᵢ − (n/2)·log 2π.
        let mut logdet_half = 0.0f64;
        for i in 0..n {
            logdet_half += fs_math::det::ln(chol.l(i, i));
        }
        let yta: f64 = y.iter().zip(&alpha).map(|(a, b)| a * b).sum();
        let lml = (-0.5f64).mul_add(
            yta,
            -logdet_half - 0.5 * n as f64 * fs_math::det::ln(2.0 * core::f64::consts::PI),
        );
        Some(Gp {
            kernel,
            noise,
            x: x.to_vec(),
            alpha,
            chol,
            lml,
        })
    }

    /// Posterior mean and variance at a point (latent f, no noise).
    #[must_use]
    pub fn predict(&self, xs: &[f64]) -> (f64, f64) {
        let n = self.x.len();
        let kstar: Vec<f64> = self.x.iter().map(|xi| self.kernel.eval(xi, xs)).collect();
        let mean: f64 = kstar.iter().zip(&self.alpha).map(|(a, b)| a * b).sum();
        // v = L⁻¹k*: forward substitution.
        let mut v = kstar;
        forward_sub(&self.chol, &mut v, n);
        let kss = self.kernel.eval(xs, xs);
        let var = (kss - v.iter().map(|t| t * t).sum::<f64>()).max(0.0);
        (mean, var)
    }

    /// Posterior mean vector, variance vector, and the CROSS-covariance
    /// Cholesky factor for q points (the q-EI reparameterization
    /// substrate): returns (μ, L_post) with Σ_post = L_postL_postᵀ.
    #[must_use]
    pub fn predict_joint(&self, xs: &[Vec<f64>]) -> (Vec<f64>, Vec<f64>) {
        let q = xs.len();
        let n = self.x.len();
        let mut mu = vec![0.0f64; q];
        // V column b = L⁻¹ k*(x_b) — reused for the posterior covariance.
        let mut vcols = vec![vec![0.0f64; n]; q];
        for (b, xb) in xs.iter().enumerate() {
            let kstar: Vec<f64> = self.x.iter().map(|xi| self.kernel.eval(xi, xb)).collect();
            mu[b] = kstar.iter().zip(&self.alpha).map(|(a, c)| a * c).sum();
            let v = &mut vcols[b];
            v.copy_from_slice(&kstar);
            forward_sub(&self.chol, v, n);
        }
        let mut sigma = vec![0.0f64; q * q];
        for a in 0..q {
            for b in 0..=a {
                let prior = self.kernel.eval(&xs[a], &xs[b]);
                let dot: f64 = vcols[a].iter().zip(&vcols[b]).map(|(p, r)| p * r).sum();
                let val = prior - dot;
                sigma[a * q + b] = val;
                sigma[b * q + a] = val;
            }
            // Jitter for the reparameterization factor.
            sigma[a * q + a] = (sigma[a * q + a]).max(0.0) + 1e-10;
        }
        let lp = cholesky(&sigma, q).expect("posterior covariance must be SPD (jittered)");
        let mut lflat = vec![0.0f64; q * q];
        for i in 0..q {
            for j in 0..=i {
                lflat[i * q + j] = lp.l(i, j);
            }
        }
        (mu, lflat)
    }
}

/// Forward substitution v ← L⁻¹v (shared by the predictive paths).
fn forward_sub(chol: &Cholesky, v: &mut [f64], n: usize) {
    for i in 0..n {
        let mut acc = v[i];
        for (j, vj) in v.iter().enumerate().take(i) {
            acc = (-chol.l(i, j)).mul_add(*vj, acc);
        }
        v[i] = acc / chol.l(i, i);
    }
}

/// Hyperparameter fit: maximize LML over log-parameters
/// (lengthscales ×D, signal, noise) by fs-ascent L-BFGS with
/// FD gradients (≤ D+2 parameters — the FD cost is trivial at
/// fixture scale; analytic LML gradients are recorded follow-up),
/// QMC-multistarted from `starts` scrambled-Sobol points in
/// `log_box = [lo, hi]` per parameter. Deterministic per seed.
#[must_use]
pub fn fit_hyperparams(
    x: &[Vec<f64>],
    y: &[f64],
    family: Matern,
    log_box: (f64, f64),
    starts: usize,
    seed: u64,
) -> Gp {
    let d = x[0].len();
    let np = d + 2;
    let sobol = fs_rand::qmc::Sobol::scrambled(np, seed);
    let nll = |p: &[f64]| -> f64 {
        let kernel = Kernel {
            family,
            signal: fs_math::det::exp(p[d].clamp(-25.0, 25.0)),
            lengthscales: p[..d]
                .iter()
                .map(|v| fs_math::det::exp(v.clamp(-12.0, 12.0)))
                .collect(),
        };
        let noise = fs_math::det::exp(p[d + 1].clamp(-16.0, 7.0)).max(1e-8);
        Gp::try_fit(x, y, kernel, noise).map_or(f64::INFINITY, |g| -g.lml)
    };
    let (lo, hi) = log_box;
    let mut best: Option<(f64, Vec<f64>)> = None;
    let mut pt = vec![0.0f64; np];
    for s in 0..starts {
        sobol.point(u32::try_from(s).expect("few starts"), &mut pt);
        let p0: Vec<f64> = pt.iter().map(|u| (hi - lo).mul_add(*u, lo)).collect();
        let mut fg = |p: &[f64]| -> (f64, Vec<f64>) {
            let f0 = nll(p);
            let mut g = vec![0.0f64; np];
            let eps = 1e-5;
            for i in 0..np {
                let mut pp = p.to_vec();
                pp[i] += eps;
                let mut pm = p.to_vec();
                pm[i] -= eps;
                g[i] = (nll(&pp) - nll(&pm)) / (2.0 * eps);
            }
            (f0, g)
        };
        let mut st = fs_ascent::LbfgsState::new(&p0, 8, &mut fg);
        let rep = st.run(&mut fg, &fs_ascent::StopRule::GradNorm(1e-6), 60);
        let better = best.as_ref().is_none_or(|(bf, _)| rep.f < *bf);
        if better {
            best = Some((rep.f, st.x.clone()));
        }
    }
    let (_, p) = best.expect("at least one start");
    let kernel = Kernel {
        family,
        signal: fs_math::det::exp(p[d].clamp(-25.0, 25.0)),
        lengthscales: p[..d]
            .iter()
            .map(|v| fs_math::det::exp(v.clamp(-12.0, 12.0)))
            .collect(),
    };
    let noise = fs_math::det::exp(p[d + 1].clamp(-16.0, 7.0)).max(1e-8);
    // The winning candidate was SPD when scored; re-fit is safe.
    Gp::fit(x, y, kernel, noise)
}
