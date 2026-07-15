//! Multi-fidelity Bayesian optimization: a two-fidelity joint GP via
//! the intrinsic-coregionalization (ICM) kernel
//! K((x,m),(x',m')) = B[m][m']·k_x(x,x') with B = LLᵀ (2×2,
//! Cholesky-parameterized — PSD by construction; the between-fidelity
//! correlation is a LEARNED hyperparameter), plus the MFEI-class
//! cost-aware acquisition: evaluate cheap when the posterior
//! between-fidelity correlation times the cost ratio beats one.

use crate::acq::{phi_cdf, phi_pdf};
use fs_la::factor::{Cholesky, cholesky};

/// ICM hyperparameters.
#[derive(Debug, Clone)]
pub struct MfKernel {
    /// ARD lengthscales for the spatial Matérn-5⁄2 part.
    pub lengthscales: Vec<f64>,
    /// Fidelity-covariance Cholesky factor entries (l11, l21, l22):
    /// B = LLᵀ with L = [[l11, 0], [l21, l22]].
    pub l_fid: [f64; 3],
    /// Observation-noise variance.
    pub noise: f64,
}

impl MfKernel {
    fn assert_dims(&self, x: &[f64], y: &[f64]) {
        assert_eq!(
            x.len(),
            y.len(),
            "MF kernel inputs must have the same dimension"
        );
        assert_eq!(
            x.len(),
            self.lengthscales.len(),
            "MF kernel input dimension must match ARD lengthscales"
        );
    }

    /// The 2×2 fidelity covariance B = LLᵀ.
    #[must_use]
    pub fn b(&self) -> [[f64; 2]; 2] {
        let [l11, l21, l22] = self.l_fid;
        [
            [l11 * l11, l11 * l21],
            [l11 * l21, l21.mul_add(l21, l22 * l22)],
        ]
    }

    /// Learned between-fidelity correlation ρ = B01/√(B00·B11).
    #[must_use]
    pub fn correlation(&self) -> f64 {
        let b = self.b();
        b[0][1] / fs_math::det::sqrt(b[0][0] * b[1][1]).max(1e-30)
    }

    fn k_x(&self, x: &[f64], y: &[f64]) -> f64 {
        self.assert_dims(x, y);
        let mut acc = 0.0f64;
        for ((xi, yi), l) in x.iter().zip(y).zip(&self.lengthscales) {
            let d = (xi - yi) / l;
            acc = d.mul_add(d, acc);
        }
        let r = fs_math::det::sqrt(acc);
        if !r.is_finite() || r > 1e8 {
            return 0.0;
        }
        let a = fs_math::det::sqrt(5.0) * r;
        (1.0 + a + a * a / 3.0) * fs_math::det::exp(-a)
    }

    /// Full kernel value between (x, fidelity m) and (y, fidelity n).
    #[must_use]
    pub fn eval(&self, x: &[f64], m: usize, y: &[f64], n: usize) -> f64 {
        assert!(m < 2 && n < 2, "fidelity index must be 0 or 1");
        self.b()[m][n] * self.k_x(x, y)
    }
}

/// A fitted two-fidelity GP (zero prior mean; fidelity 1 = HIGH).
pub struct MfGp {
    /// Kernel used.
    pub kernel: MfKernel,
    x: Vec<Vec<f64>>,
    fid: Vec<usize>,
    alpha: Vec<f64>,
    chol: Cholesky,
    /// Log marginal likelihood.
    pub lml: f64,
}

impl MfGp {
    /// Fallible exact fit.
    #[must_use]
    pub fn try_fit(x: &[Vec<f64>], fid: &[usize], y: &[f64], kernel: MfKernel) -> Option<MfGp> {
        let n = x.len();
        assert_eq!(n, fid.len());
        assert_eq!(n, y.len());
        let mut k = vec![0.0f64; n * n];
        for i in 0..n {
            for j in 0..=i {
                let v = kernel.eval(&x[i], fid[i], &x[j], fid[j]);
                k[i * n + j] = v;
                k[j * n + i] = v;
            }
            k[i * n + i] += kernel.noise;
        }
        let chol = cholesky(&k, n).ok()?;
        let mut alpha = y.to_vec();
        chol.solve(&mut alpha);
        let mut logdet_half = 0.0f64;
        for i in 0..n {
            logdet_half += fs_math::det::ln(chol.l(i, i));
        }
        let yta: f64 = y.iter().zip(&alpha).map(|(a, b)| a * b).sum();
        let lml = (-0.5f64).mul_add(
            yta,
            -logdet_half - 0.5 * n as f64 * fs_math::det::ln(2.0 * core::f64::consts::PI),
        );
        Some(MfGp {
            kernel,
            x: x.to_vec(),
            fid: fid.to_vec(),
            alpha,
            chol,
            lml,
        })
    }

    fn kstar(&self, xs: &[f64], m: usize) -> Vec<f64> {
        self.x
            .iter()
            .zip(&self.fid)
            .map(|(xi, &fi)| self.kernel.eval(xi, fi, xs, m))
            .collect()
    }

    fn lsolve(&self, v: &mut [f64]) {
        let n = v.len();
        for i in 0..n {
            let mut acc = v[i];
            for (j, vj) in v.iter().enumerate().take(i) {
                acc = (-self.chol.l(i, j)).mul_add(*vj, acc);
            }
            v[i] = acc / self.chol.l(i, i);
        }
    }

    /// Posterior mean and variance at (x, fidelity m).
    #[must_use]
    pub fn predict(&self, xs: &[f64], m: usize) -> (f64, f64) {
        let kstar = self.kstar(xs, m);
        let mean: f64 = kstar.iter().zip(&self.alpha).map(|(a, b)| a * b).sum();
        let mut v = kstar;
        self.lsolve(&mut v);
        let kss = self.kernel.eval(xs, m, xs, m);
        let var = (kss - v.iter().map(|t| t * t).sum::<f64>()).max(0.0);
        (mean, var)
    }

    /// Posterior CORRELATION between f_low(x) and f_high(x) — the
    /// fidelity-selection signal.
    #[must_use]
    pub fn posterior_fid_correlation(&self, xs: &[f64]) -> f64 {
        let klo = self.kstar(xs, 0);
        let khi = self.kstar(xs, 1);
        let mut vlo = klo;
        let mut vhi = khi;
        self.lsolve(&mut vlo);
        self.lsolve(&mut vhi);
        let cross_prior = self.kernel.eval(xs, 0, xs, 1);
        let cross: f64 = cross_prior - vlo.iter().zip(&vhi).map(|(a, b)| a * b).sum::<f64>();
        let (_, var_lo) = self.predict(xs, 0);
        let (_, var_hi) = self.predict(xs, 1);
        cross / fs_math::det::sqrt(var_lo * var_hi).max(1e-30)
    }
}

/// Fit ICM hyperparameters by L-BFGS on −LML with hybrid QMC
/// multistart (the house pattern). Log-parameterization: D
/// lengthscales, 3 fidelity-Cholesky entries (l21 signed via tanh of
/// the raw parameter times the l11 scale), noise.
#[must_use]
pub fn fit_mf(
    x: &[Vec<f64>],
    fid: &[usize],
    y: &[f64],
    log_box: (f64, f64),
    starts: usize,
    seed: u64,
) -> MfGp {
    let d = x[0].len();
    let np = d + 4;
    let kq = np.min(fs_rand::qmc::MAX_SOBOL_DIM);
    let sobol = fs_rand::qmc::Sobol::scrambled(kq, seed);
    let build = |p: &[f64]| -> MfKernel {
        let lengthscales: Vec<f64> = p[..d]
            .iter()
            .map(|v| fs_math::det::exp(v.clamp(-8.0, 8.0)))
            .collect();
        let l11 = fs_math::det::exp(p[d].clamp(-8.0, 8.0));
        // l21 signed: raw in R, squashed to (−2, 2)·l11.
        let t = p[d + 1].clamp(-6.0, 6.0);
        let e2 = fs_math::det::exp(2.0 * t);
        let l21 = 2.0 * l11 * ((e2 - 1.0) / (e2 + 1.0));
        let l22 = fs_math::det::exp(p[d + 2].clamp(-8.0, 8.0));
        let noise = fs_math::det::exp(p[d + 3].clamp(-16.0, 4.0)).max(1e-8);
        MfKernel {
            lengthscales,
            l_fid: [l11, l21, l22],
            noise,
        }
    };
    let nll =
        |p: &[f64]| -> f64 { MfGp::try_fit(x, fid, y, build(p)).map_or(f64::INFINITY, |g| -g.lml) };
    let (lo, hi) = log_box;
    let mut best: Option<(f64, Vec<f64>)> = None;
    let mut pt = vec![0.0f64; kq];
    for s in 0..starts {
        sobol.point(u32::try_from(s + 1).expect("few starts"), &mut pt);
        let mut tail = fs_rand::StreamKey {
            seed,
            kernel: 0x30F1,
            tile: u32::try_from(s).expect("few starts"),
        }
        .stream();
        let p0: Vec<f64> = (0..np)
            .map(|i| {
                let u = if i < kq { pt[i] } else { tail.next_f64() };
                (hi - lo).mul_add(u, lo)
            })
            .collect();
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
        let rep = st.run(&mut fg, &fs_ascent::StopRule::GradNorm(1e-6), 50);
        if best.as_ref().is_none_or(|(bf, _)| rep.f < *bf) {
            best = Some((rep.f, st.x.clone()));
        }
    }
    let (_, p) = best.expect("at least one start");
    MfGp::try_fit(x, fid, y, build(&p)).expect("winning candidate was SPD when scored")
}

/// Multi-fidelity BO configuration.
#[derive(Debug, Clone)]
pub struct MfConfig {
    /// Search box.
    pub bounds: (f64, f64),
    /// Hyperparameter log-box.
    pub log_box: (f64, f64),
    /// Hyperparameter multistarts.
    pub hyper_starts: usize,
    /// Initial low-fidelity points.
    pub n_init_low: usize,
    /// Initial high-fidelity points.
    pub n_init_high: usize,
    /// Cost of one low-fidelity evaluation.
    pub cost_low: f64,
    /// Cost of one high-fidelity evaluation.
    pub cost_high: f64,
    /// Total cost budget.
    pub budget: f64,
    /// Acquisition CMA-ES restarts.
    pub acq_starts: usize,
    /// Acquisition evaluations per restart.
    pub acq_evals: usize,
    /// Root seed for all derived study streams.
    pub seed: u64,
}

/// Outcome of an MF-BO run.
#[derive(Debug, Clone)]
pub struct MfReport {
    /// Best HIGH-fidelity value observed.
    pub f_best_high: f64,
    /// Low-fidelity evaluations.
    pub evals_low: usize,
    /// High-fidelity evaluations.
    pub evals_high: usize,
    /// Total cost spent.
    pub cost: f64,
    /// (cost, best-high-so-far) curve — the ledgered evidence.
    pub trace: Vec<(f64, f64)>,
    /// Final learned between-fidelity correlation.
    pub learned_correlation: f64,
}

/// Run cost-aware two-fidelity BO for MINIMIZATION of the HIGH
/// fidelity. `f(x, m)` evaluates fidelity m ∈ {0 = low, 1 = high}.
#[allow(clippy::too_many_lines)] // one coherent cost-aware BO loop
pub fn mf_minimize(
    f: &mut dyn FnMut(&[f64], usize) -> f64,
    dim: usize,
    config: &MfConfig,
) -> MfReport {
    let (lo, hi) = config.bounds;
    let span = hi - lo;
    let kq = dim.min(fs_rand::qmc::MAX_SOBOL_DIM);
    let sobol = fs_rand::qmc::Sobol::scrambled(kq, config.seed);
    let mut tail = fs_rand::StreamKey {
        seed: config.seed,
        kernel: 0x30F2,
        tile: 0,
    }
    .stream();
    let mut pt = vec![0.0f64; kq];
    let mut xs: Vec<Vec<f64>> = Vec::new();
    let mut fid: Vec<usize> = Vec::new();
    let mut ys: Vec<f64> = Vec::new();
    let mut cost = 0.0f64;
    let mut trace: Vec<(f64, f64)> = Vec::new();
    let mut best_high = f64::INFINITY;
    let draw = |s: usize, pt: &mut [f64], tail: &mut fs_rand::Stream| -> Vec<f64> {
        sobol.point(u32::try_from(s + 1).expect("few inits"), pt);
        (0..dim)
            .map(|i| {
                let u = if i < kq { pt[i] } else { tail.next_f64() };
                span.mul_add(u, lo)
            })
            .collect()
    };
    for s in 0..config.n_init_low {
        let x = draw(s, &mut pt, &mut tail);
        let y = f(&x, 0);
        cost += config.cost_low;
        xs.push(x);
        fid.push(0);
        ys.push(y);
        trace.push((cost, best_high));
    }
    for s in 0..config.n_init_high {
        let x = draw(config.n_init_low + s, &mut pt, &mut tail);
        let y = f(&x, 1);
        cost += config.cost_high;
        best_high = best_high.min(y);
        xs.push(x);
        fid.push(1);
        ys.push(y);
        trace.push((cost, best_high));
    }
    let mut learned_correlation = 0.0f64;
    let mut it = 0u64;
    let mut evals_low = config.n_init_low;
    let mut evals_high = config.n_init_high;
    while cost + config.cost_low <= config.budget {
        it += 1;
        // Standardize on the joint set.
        let n = ys.len() as f64;
        let mean = ys.iter().sum::<f64>() / n;
        let var = ys.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n;
        let scale = fs_math::det::sqrt(var.max(1e-30));
        let ys_std: Vec<f64> = ys.iter().map(|v| (v - mean) / scale).collect();
        let gp = fit_mf(
            &xs,
            &fid,
            &ys_std,
            config.log_box,
            config.hyper_starts,
            config.seed ^ it.wrapping_mul(0x9e37_79b9),
        );
        learned_correlation = gp.kernel.correlation();
        let f_best_std = (best_high - mean) / scale;
        // EI on the HIGH-fidelity posterior, maximized by CMA-ES.
        let ei = |x: &[f64]| -> f64 {
            let (mu, var) = gp.predict(x, 1);
            let sigma = fs_math::det::sqrt(var.max(1e-18));
            let delta = f_best_std - mu;
            let z = delta / sigma;
            delta.mul_add(phi_cdf(z), sigma * phi_pdf(z)).max(0.0)
        };
        let cand_sobol = fs_rand::qmc::Sobol::scrambled(kq, config.seed ^ 0x5EED ^ (it << 8));
        let mut best_x: Option<Vec<f64>> = None;
        let mut best_v = f64::NEG_INFINITY;
        for s in 0..config.acq_starts {
            cand_sobol.point(u32::try_from(s + 1).expect("few starts"), &mut pt);
            let x0: Vec<f64> = (0..dim)
                .map(|i| {
                    let u = if i < kq { pt[i] } else { tail.next_f64() };
                    span.mul_add(u, lo)
                })
                .collect();
            let mut obj = |x: &[f64]| -> f64 {
                let xc: Vec<f64> = x.iter().map(|v| v.clamp(lo, hi)).collect();
                -ei(&xc)
            };
            let params =
                fs_dfo::CmaParams::standard(dim, 0.2 * span, config.acq_evals, f64::NEG_INFINITY);
            let rep = fs_dfo::cmaes(&mut obj, &x0, &params, config.seed ^ (s as u64) << 4);
            if -rep.f_best > best_v {
                best_v = -rep.f_best;
                best_x = Some(rep.x_best.iter().map(|v| v.clamp(lo, hi)).collect());
            }
        }
        let x_new = best_x.expect("at least one acquisition start");
        // MFEI-class fidelity choice: cheap wins when
        // corr²·(cost_high/cost_low) > 1 (variance-reduction per cost).
        let corr = gp.posterior_fid_correlation(&x_new).clamp(-1.0, 1.0);
        let gain_low = corr * corr * (config.cost_high / config.cost_low);
        let take_low = gain_low > 1.0 && cost + config.cost_low <= config.budget;
        if take_low {
            let y = f(&x_new, 0);
            cost += config.cost_low;
            evals_low += 1;
            xs.push(x_new);
            fid.push(0);
            ys.push(y);
        } else {
            if cost + config.cost_high > config.budget {
                break;
            }
            let y = f(&x_new, 1);
            cost += config.cost_high;
            evals_high += 1;
            best_high = best_high.min(y);
            xs.push(x_new);
            fid.push(1);
            ys.push(y);
        }
        trace.push((cost, best_high));
    }
    MfReport {
        f_best_high: best_high,
        evals_low,
        evals_high,
        cost,
        trace,
        learned_correlation,
    }
}
