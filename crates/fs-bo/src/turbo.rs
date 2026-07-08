//! TuRBO-class trust-region Bayesian optimization — the honest answer
//! to BO's dimensionality ceiling: a LOCAL GP inside an adaptive
//! hyperrectangle (sides weighted by ARD lengthscales), Thompson
//! sampling over Sobol candidates through the joint-posterior Cholesky
//! with a FIXED normal bank (common random numbers — every run replays
//! bitwise), success/failure counters driving expand/shrink, and
//! restarts on collapse that keep the global best.

use crate::gp::{Gp, Matern, fit_hyperparams};

/// TuRBO configuration.
#[derive(Debug, Clone)]
pub struct TurboConfig {
    /// Box bounds (same per dimension).
    pub bounds: (f64, f64),
    /// Kernel family.
    pub family: Matern,
    /// Hyperparameter log-box.
    pub log_box: (f64, f64),
    /// Hyperparameter multistarts.
    pub hyper_starts: usize,
    /// Initial points per (re)start.
    pub n_init: usize,
    /// Candidates per iteration.
    pub candidates: usize,
    /// Total evaluation budget.
    pub max_evals: usize,
    /// Initial TR side (fraction of the box).
    pub l_init: f64,
    /// Collapse threshold.
    pub l_min: f64,
    /// Expansion cap.
    pub l_max: f64,
    /// Successes before doubling.
    pub succ_tol: usize,
    /// Failures before halving.
    pub fail_tol: usize,
    /// Master seed.
    pub seed: u64,
    /// Local-data cap: the GP trains on at most this many points
    /// nearest the region center (O(n³) inference — unbounded local
    /// sets measurably stalled the 30-d battery).
    pub max_local: usize,
    /// Refit hyperparameters every this many iterations (reusing the
    /// last kernel in between — the standard TuRBO economy).
    pub refit_every: u64,
}

/// Outcome of a TuRBO run.
#[derive(Debug, Clone)]
pub struct TurboReport {
    /// Best point found.
    pub x_best: Vec<f64>,
    /// Best value found.
    pub f_best: f64,
    /// Evaluations spent.
    pub evals: usize,
    /// Restarts performed.
    pub restarts: usize,
    /// Best-so-far trace (per evaluation — the ledgered curve).
    pub trace: Vec<f64>,
}

/// Run TuRBO for MINIMIZATION.
#[allow(clippy::too_many_lines)] // one coherent restart/TR loop
pub fn turbo_minimize(
    f: &mut dyn FnMut(&[f64]) -> f64,
    dim: usize,
    config: &TurboConfig,
) -> TurboReport {
    let (lo, hi) = config.bounds;
    assert!(dim > 0, "TuRBO requires dim > 0");
    assert!(
        lo.is_finite() && hi.is_finite() && hi > lo,
        "TuRBO bounds must be finite and ordered"
    );
    assert!(config.hyper_starts > 0, "TuRBO requires hyper_starts > 0");
    assert!(config.n_init > 0, "TuRBO requires n_init > 0");
    assert!(config.candidates > 0, "TuRBO requires candidates > 0");
    assert!(config.max_local > 0, "TuRBO requires max_local > 0");
    assert!(config.refit_every > 0, "TuRBO requires refit_every > 0");
    assert!(config.succ_tol > 0, "TuRBO requires succ_tol > 0");
    assert!(config.fail_tol > 0, "TuRBO requires fail_tol > 0");
    assert!(
        config.l_min.is_finite()
            && config.l_init.is_finite()
            && config.l_max.is_finite()
            && config.l_min > 0.0
            && config.l_init > 0.0
            && config.l_max >= config.l_init,
        "TuRBO trust-region lengths must be positive finite values with l_max >= l_init"
    );
    let span = hi - lo;
    let mut evals = 0usize;
    let mut restarts = 0usize;
    let mut global_best = f64::INFINITY;
    let mut global_x = vec![0.0f64; dim];
    let mut trace = Vec::with_capacity(config.max_evals);
    let mut round = 0u64;
    'restart: while evals < config.max_evals {
        // Fresh local history per (re)start.
        let init_sobol = fs_rand::qmc::Sobol::scrambled(
            dim.min(fs_rand::qmc::MAX_SOBOL_DIM),
            config.seed ^ round,
        );
        let mut tail = fs_rand::StreamKey {
            seed: config.seed ^ 0x7B0,
            kernel: 0x7B00,
            tile: u32::try_from(round & 0xFFFF).expect("fits"),
        }
        .stream();
        let kq = dim.min(fs_rand::qmc::MAX_SOBOL_DIM);
        let mut pt = vec![0.0f64; kq];
        let mut xs: Vec<Vec<f64>> = Vec::new();
        let mut ys: Vec<f64> = Vec::new();
        for s in 0..config.n_init {
            if evals >= config.max_evals {
                break;
            }
            init_sobol.point(u32::try_from(s + 1).expect("few inits"), &mut pt);
            let mut x: Vec<f64> = pt.iter().map(|u| span.mul_add(*u, lo)).collect();
            while x.len() < dim {
                x.push(span.mul_add(tail.next_f64(), lo));
            }
            let y = f(&x);
            evals += 1;
            if y < global_best {
                global_best = y;
                global_x.clone_from(&x);
            }
            trace.push(global_best);
            xs.push(x);
            ys.push(y);
        }
        let mut length = config.l_init;
        let mut succ = 0usize;
        let mut fail = 0usize;
        let mut iter_in_round = 0u64;
        let mut kernel_cache: Option<crate::gp::Kernel> = None;
        let mut noise_cache = 1e-6f64;
        while evals < config.max_evals {
            iter_in_round += 1;
            // Local best anchors the region.
            let (bi, _) = ys
                .iter()
                .enumerate()
                .min_by(|a, b| a.1.total_cmp(b.1).then(a.0.cmp(&b.0)))
                .expect("nonempty");
            let center = xs[bi].clone();
            let f_local_best = ys[bi];
            // Fit the local GP on in-region points (fall back to all
            // points when the region is data-poor).
            let half = 0.5 * length * span;
            let in_region: Vec<usize> = (0..xs.len())
                .filter(|&i| {
                    xs[i]
                        .iter()
                        .zip(&center)
                        .all(|(a, c)| (a - c).abs() <= half + 1e-12)
                })
                .collect();
            let mut use_idx: Vec<usize> = if in_region.len() >= config.n_init.min(2 * dim) {
                in_region
            } else {
                // Data-poor region: use everything.
                (0..xs.len()).collect()
            };
            // Cap by distance to the center (deterministic tie-break).
            if use_idx.len() > config.max_local {
                let d2 = |i: usize| -> f64 {
                    xs[i]
                        .iter()
                        .zip(&center)
                        .map(|(a, c)| (a - c) * (a - c))
                        .sum()
                };
                use_idx.sort_by(|&a, &b| d2(a).total_cmp(&d2(b)).then(a.cmp(&b)));
                use_idx.truncate(config.max_local);
            }
            let xl: Vec<Vec<f64>> = use_idx.iter().map(|&i| xs[i].clone()).collect();
            let yl_raw: Vec<f64> = use_idx.iter().map(|&i| ys[i]).collect();
            let n = yl_raw.len() as f64;
            let mean = yl_raw.iter().sum::<f64>() / n;
            let var = yl_raw.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n;
            let scale = fs_math::det::sqrt(var.max(1e-30));
            let yl: Vec<f64> = yl_raw.iter().map(|v| (v - mean) / scale).collect();
            let refit =
                kernel_cache.is_none() || (iter_in_round - 1).is_multiple_of(config.refit_every);
            let gp = if refit {
                let g = fit_hyperparams(
                    &xl,
                    &yl,
                    config.family,
                    config.log_box,
                    config.hyper_starts,
                    config.seed ^ round.wrapping_mul(0x9e37_79b9) ^ iter_in_round,
                );
                kernel_cache = Some(g.kernel.clone());
                noise_cache = g.noise;
                g
            } else {
                let kernel = kernel_cache.clone().expect("cache primed");
                Gp::try_fit(&xl, &yl, kernel, noise_cache).unwrap_or_else(|| {
                    let g = fit_hyperparams(
                        &xl,
                        &yl,
                        config.family,
                        config.log_box,
                        config.hyper_starts,
                        config.seed ^ round.wrapping_mul(0x9e37_79b9) ^ iter_in_round,
                    );
                    kernel_cache = Some(g.kernel.clone());
                    noise_cache = g.noise;
                    g
                })
            };
            // ARD-weighted TR sides: side_d ∝ ℓ_d, normalized to the
            // geometric mean (long lengthscales get wide sides).
            let ell = &gp.kernel.lengthscales;
            let log_gm: f64 = ell.iter().map(|l| fs_math::det::ln(*l)).sum::<f64>() / dim as f64;
            let gm = fs_math::det::exp(log_gm);
            let sides: Vec<f64> = ell
                .iter()
                .map(|l| (l / gm).clamp(0.3, 3.0) * length * span)
                .collect();
            // Sobol candidates in the weighted box (clamped to bounds).
            let cand_sobol = fs_rand::qmc::Sobol::scrambled(
                kq,
                config.seed ^ 0xCAFE ^ (round << 20) ^ iter_in_round,
            );
            let mut cands: Vec<Vec<f64>> = Vec::with_capacity(config.candidates);
            for s in 0..config.candidates {
                cand_sobol.point(u32::try_from(s + 1).expect("few candidates"), &mut pt);
                let mut x: Vec<f64> = Vec::with_capacity(dim);
                for d in 0..dim {
                    let u = if d < kq { pt[d] } else { tail.next_f64() };
                    let v = (u - 0.5).mul_add(sides[d], center[d]);
                    x.push(v.clamp(lo, hi));
                }
                cands.push(x);
            }
            // Thompson sample over the candidate set. The noise bank
            // comes from a Philox stream (candidate counts exceed the
            // Sobol table's 10-dim cap, and Thompson noise gains
            // nothing from QMC) — fixed per (seed, round, iteration):
            // common random numbers, bitwise-replayable.
            let (mu, lchol) = gp.predict_joint(&cands);
            let q = cands.len();
            let mut zs = fs_rand::StreamKey {
                seed: config.seed ^ 0x7541,
                kernel: 0x754B,
                tile: u32::try_from((round << 16 | iter_in_round) & 0xFFFF_FFFF).expect("fits"),
            }
            .stream();
            let bank: Vec<f64> = (0..q).map(|_| zs.next_normal()).collect();
            let mut best_c = 0usize;
            let mut best_v = f64::INFINITY;
            for i in 0..q {
                let mut v = mu[i];
                for j in 0..=i {
                    v = lchol[i * q + j].mul_add(bank[j], v);
                }
                if v < best_v {
                    best_v = v;
                    best_c = i;
                }
            }
            let x_new = cands[best_c].clone();
            let y_new = f(&x_new);
            evals += 1;
            if y_new < global_best {
                global_best = y_new;
                global_x.clone_from(&x_new);
            }
            trace.push(global_best);
            // Counters: success = improved the LOCAL best.
            if y_new < f_local_best - 1e-12 * f_local_best.abs() {
                succ += 1;
                fail = 0;
            } else {
                fail += 1;
                succ = 0;
            }
            if succ >= config.succ_tol {
                length = (2.0 * length).min(config.l_max);
                succ = 0;
            }
            if fail >= config.fail_tol {
                length *= 0.5;
                fail = 0;
            }
            xs.push(x_new);
            ys.push(y_new);
            if length < config.l_min {
                restarts += 1;
                round += 1;
                continue 'restart;
            }
        }
        break;
    }
    TurboReport {
        x_best: global_x,
        f_best: global_best,
        evals,
        restarts,
        trace,
    }
}
