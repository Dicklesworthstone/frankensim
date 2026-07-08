//! The BO loop: fit (QMC-multistart hyperparameters) → optimize the
//! acquisition (fs-dfo CMA-ES multistarted from Sobol points — a
//! deterministic derivative-free inner optimizer; the FrankenTorch
//! reparameterized-gradient tape is the recorded follow-up lane) →
//! evaluate → repeat. Boxes are [lo, hi]^d; every run is
//! deterministic per seed and replays bitwise.

use crate::acq::{expected_improvement, normal_bank, q_expected_improvement};
use crate::gp::{Matern, fit_hyperparams};

/// BO configuration.
#[derive(Debug, Clone)]
pub struct BoConfig {
    /// Search box (lo, hi) per dimension (same for all dims).
    pub bounds: (f64, f64),
    /// Kernel family.
    pub family: Matern,
    /// Hyperparameter log-box.
    pub log_box: (f64, f64),
    /// Hyperparameter multistarts.
    pub hyper_starts: usize,
    /// Acquisition CMA-ES restarts.
    pub acq_starts: usize,
    /// Acquisition evaluation budget per restart.
    pub acq_evals: usize,
    /// Batch size q (1 = sequential EI).
    pub q: usize,
    /// MC samples for q-EI.
    pub mc_samples: usize,
    /// Master seed.
    pub seed: u64,
}

/// Outcome of a BO run.
#[derive(Debug, Clone)]
pub struct BoReport {
    /// All evaluated points.
    pub x: Vec<Vec<f64>>,
    /// All observed values.
    pub y: Vec<f64>,
    /// Best value per iteration (the anytime curve — the ledgered
    /// evidence for baseline comparisons).
    pub best_trace: Vec<f64>,
}

fn clamp_box(x: &mut [f64], lo: f64, hi: f64) {
    for v in x {
        *v = v.clamp(lo, hi);
    }
}

/// Run BO for MINIMIZATION from `n_init` Sobol seeds, `iters` batches
/// of `q` acquisitions each.
pub fn minimize(
    f: &mut dyn FnMut(&[f64]) -> f64,
    dim: usize,
    n_init: usize,
    iters: usize,
    config: &BoConfig,
) -> BoReport {
    let (lo, hi) = config.bounds;
    let sobol = fs_rand::qmc::Sobol::scrambled(dim, config.seed);
    let mut xs: Vec<Vec<f64>> = Vec::new();
    let mut ys: Vec<f64> = Vec::new();
    let mut pt = vec![0.0f64; dim];
    for s in 0..n_init {
        sobol.point(u32::try_from(s + 1).expect("few inits"), &mut pt);
        let x: Vec<f64> = pt.iter().map(|u| (hi - lo).mul_add(*u, lo)).collect();
        ys.push(f(&x));
        xs.push(x);
    }
    let mut best_trace: Vec<f64> = vec![ys.iter().copied().fold(f64::INFINITY, f64::min)];
    for it in 0..iters {
        // STANDARDIZE observations (zero mean, unit variance) before
        // fitting: EI is invariant under consistent affine maps of y,
        // and without this the signal-variance search box cannot span
        // objectives at arbitrary scale (measured: raw-y BO LOST to
        // random on Branin, whose values span O(300)).
        let n = ys.len() as f64;
        let mean = ys.iter().sum::<f64>() / n;
        let var = ys.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n;
        let scale = fs_math::det::sqrt(var.max(1e-30));
        let ys_std: Vec<f64> = ys.iter().map(|v| (v - mean) / scale).collect();
        let gp = fit_hyperparams(
            &xs,
            &ys_std,
            config.family,
            config.log_box,
            config.hyper_starts,
            config.seed ^ (it as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15),
        );
        let f_best = ys_std.iter().copied().fold(f64::INFINITY, f64::min);
        let batch: Vec<Vec<f64>> = if config.q <= 1 {
            vec![argmax_acq(
                &|x: &[f64]| expected_improvement(&gp, x, f_best, 0.0),
                dim,
                config,
                it,
            )]
        } else {
            // Greedy q-EI construction: grow the batch one point at a
            // time, each maximizing the JOINT q-EI with the already-
            // chosen points fixed (the standard sequential-greedy
            // batch heuristic; joint q-EI over the fixed z-bank keeps
            // it deterministic).
            let bank = normal_bank(config.mc_samples, config.q, config.seed ^ 0xACC5);
            let mut chosen: Vec<Vec<f64>> = Vec::new();
            for _ in 0..config.q {
                let cand = argmax_acq(
                    &|x: &[f64]| {
                        let mut trial: Vec<Vec<f64>> = chosen.clone();
                        trial.push(x.to_vec());
                        let sub_bank_cols = trial.len();
                        // Use the leading columns of the bank for the
                        // current batch size (fixed common random
                        // numbers across the greedy growth).
                        let samples = config.mc_samples;
                        let mut sub = vec![0.0f64; samples * sub_bank_cols];
                        for s in 0..samples {
                            sub[s * sub_bank_cols..(s + 1) * sub_bank_cols]
                                .copy_from_slice(&bank[s * config.q..s * config.q + sub_bank_cols]);
                        }
                        q_expected_improvement(&gp, &trial, f_best, &sub)
                    },
                    dim,
                    config,
                    it,
                );
                chosen.push(cand);
            }
            chosen
        };
        for x in batch {
            ys.push(f(&x));
            xs.push(x);
        }
        best_trace.push(ys.iter().copied().fold(f64::INFINITY, f64::min));
    }
    BoReport {
        x: xs,
        y: ys,
        best_trace,
    }
}

/// Maximize an acquisition over the box by CMA-ES restarts from Sobol
/// starts (deterministic per seed/iteration).
fn argmax_acq(
    acq: &dyn Fn(&[f64]) -> f64,
    dim: usize,
    config: &BoConfig,
    iteration: usize,
) -> Vec<f64> {
    let (lo, hi) = config.bounds;
    let sobol = fs_rand::qmc::Sobol::scrambled(dim, config.seed ^ 0x5EED ^ (iteration as u64) << 8);
    let mut best_x: Option<Vec<f64>> = None;
    let mut best_v = f64::NEG_INFINITY;
    let mut pt = vec![0.0f64; dim];
    for s in 0..config.acq_starts {
        sobol.point(u32::try_from(s + 1).expect("few starts"), &mut pt);
        let x0: Vec<f64> = pt.iter().map(|u| (hi - lo).mul_add(*u, lo)).collect();
        let mut obj = |x: &[f64]| -> f64 {
            let mut xc = x.to_vec();
            clamp_box(&mut xc, lo, hi);
            -acq(&xc)
        };
        let params =
            fs_dfo::CmaParams::standard(dim, 0.2 * (hi - lo), config.acq_evals, f64::NEG_INFINITY);
        let rep = fs_dfo::cmaes(&mut obj, &x0, &params, config.seed ^ (s as u64) << 4);
        let v = -rep.f_best;
        if v > best_v {
            let mut xb = rep.x_best;
            clamp_box(&mut xb, lo, hi);
            best_v = v;
            best_x = Some(xb);
        }
    }
    best_x.expect("at least one acquisition start")
}
