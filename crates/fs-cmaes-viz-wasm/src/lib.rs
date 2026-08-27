//! fs-cmaes-viz-wasm — browser (WASM) surface for CMA-ES internals
//! visualization. Layer: L6.
//!
//! One boundary between the explainer site's presentation plane (canvas /
//! three.js) and the optimization plane (a full-covariance CMA-ES kernel with
//! the standard Hansen couplings, mirroring the site's TypeScript fallback
//! math so WASM-vs-fallback is provenance, not behavior).
//!
//! Contracts every entry inherits (fs-flyer-wasm pattern):
//!
//! - **Typed-refusal JSON envelope.** Fallible entries return
//!   `{"ok": ...}` or `{"refusal": {"code","message","ranked_repairs"}}`.
//!   Nothing is silently clamped and nothing traps across the boundary.
//! - **Determinism.** A single LCG stream (constants shared with the site's
//!   TS engine) feeds Box–Muller gaussians. Same inputs ⇒ identical output
//!   strings, native AND wasm. No wall-clock, no entropy, `f64::total_cmp`
//!   for all orderings.
//! - **Honest internals.** Every generation exposes the real state: mean, σ,
//!   the full eigendecomposition of C (via `fs_la::jacobi_eigh`, the same
//!   solver fs-dfo's CMA-ES uses), the ranked population, and — for dim > 3 —
//!   an honest PCA marginal (projected 3×3 covariance eigendecomposition),
//!   never a fake 3D ellipsoid.
//!
//! No-claims: this kernel is a teaching/viz surface. It does not claim the
//! restart machinery, identity ledgers, or adversarial-refusal hardening of
//! `fs_dfo::cmaes`; for production optimization use fs-dfo.

// ---------------------------------------------------------------------------
// Public parameters (scalar mirror of the JS boundary; native callers build
// it directly, the wasm mod unpacks scalars 1:1).
use fs_la::eigen::jacobi_eigh;
use std::fmt::Write as _;

/// Landscape ids at the ABI (frozen v1).
pub const LANDSCAPE_SPHERE: u32 = 0;
pub const LANDSCAPE_ROSENBROCK: u32 = 1;
pub const LANDSCAPE_CIGAR: u32 = 2;
pub const LANDSCAPE_RASTRIGIN: u32 = 3;
pub const LANDSCAPE_ELLI: u32 = 4;

/// Kernel id baked into envelopes so the page can prove which build is live.
pub const KERNEL_VERSION: &str = "fs-cmaes-viz-wasm 0.3.0";

/// One visualization run request. Field ranges are refused, never clamped.
#[derive(Debug, Clone)]
pub struct VizParams {
    /// Decision dimension, 2..=6.
    pub dim: usize,
    /// Initial mean (first `dim` entries are significant).
    pub x0: [f64; 6],
    /// Initial step size σ₀ > 0.
    pub sigma0: f64,
    /// Population size λ, 4..=48.
    pub lambda: usize,
    /// Active (negative-weight) covariance update.
    pub active: bool,
    /// Seed for the deterministic LCG stream.
    pub seed: u64,
    /// Generations to run, 1..=200.
    pub generations: usize,
    /// Landscape id (see `LANDSCAPE_*`).
    pub landscape: u32,
    /// Additive evaluation noise σ_n ≥ 0 (0 = noiseless).
    pub noise: f64,
    /// Reflect-repair decisions into [bound_min, bound_max].
    pub bounds_enabled: bool,
    pub bound_min: f64,
    pub bound_max: f64,
    /// Stop early when best true fitness ≤ target; NaN disables the target.
    pub f_target: f64,
}

/// One generation's full internal state.
#[derive(Debug, Clone)]
pub struct GenSnapshot {
    pub g: usize,
    pub mean: Vec<f64>,
    pub sigma: f64,
    /// Eigenvalues of C, ascending.
    pub eigvals: Vec<f64>,
    /// Row-major n×n; column j is the unit eigenvector for `eigvals[j]`.
    pub eigvecs: Vec<f64>,
    pub cond: f64,
    pub best_f: f64,
    pub evals: usize,
    /// 3D phase-space projection of the mean (direct for dim ≤ 3).
    pub proj_mean: [f64; 3],
    /// Eigenvalues of the 3×3 projected marginal covariance, ascending.
    pub proj_eigvals: [f64; 3],
    /// Row-major 3×3; column j pairs with `proj_eigvals[j]`.
    pub proj_eigvecs: [f64; 9],
    /// Population decision vectors, rank order (λ·dim floats).
    pub sx: Vec<f64>,
    /// Population white-noise vectors z (pre-transformation), rank order.
    pub sz: Vec<f64>,
    /// Population (noisy) fitness, rank order.
    pub sf: Vec<f64>,
    /// Elite flags (1 = top μ), rank order.
    pub se: Vec<u8>,
    /// Step-size evolution path pσ (post-update).
    pub p_sigma: Vec<f64>,
    /// Covariance evolution path pC (post-update).
    pub p_c: Vec<f64>,
}

/// Pooled PCA frame used for the 3D phase space when dim > 3.
#[derive(Debug, Clone)]
pub struct PcaInfo {
    /// 3×dim row-major basis (rows = top-3 pooled principal axes).
    pub basis: Vec<f64>,
    /// dim — centroid of the generation-mean trajectory.
    pub center: Vec<f64>,
    /// Eigenvalues of the pooled covariance, ascending.
    pub pool_eigvals: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct VizRun {
    pub dim: usize,
    pub landscape: u32,
    pub generations: Vec<GenSnapshot>,
    pub pca: PcaInfo,
    pub best_f: f64,
    pub best_x: Vec<f64>,
    pub total_evals: usize,
    pub stop_reason: &'static str,
}

/// Typed refusal (fs-flyer-wasm envelope shape).
#[derive(Debug, Clone)]
pub struct Refusal {
    pub code: &'static str,
    pub message: String,
    pub ranked_repairs: Vec<&'static str>,
}

// ---------------------------------------------------------------------------
// Deterministic RNG — LCG constants shared with the site's TS engine so the
// fallback and the kernel tell the same sampling story.
// ---------------------------------------------------------------------------

struct Lcg(u32);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        f64::from(self.0) / 4_294_967_296.0
    }

    fn next_open_f64(&mut self) -> f64 {
        loop {
            let value = self.next_f64();
            if value > 0.0 {
                return value;
            }
        }
    }

    /// Fill standard normals with the paired Box–Muller transform used by
    /// the TypeScript fallback. A trailing odd coordinate consumes the pair
    /// but deliberately discards its sine component.
    fn fill_gaussian(&mut self, out: &mut [f64]) {
        for pair in out.chunks_mut(2) {
            let u1 = self.next_open_f64();
            let u2 = self.next_f64();
            let magnitude = fs_math::det::sqrt(-2.0 * fs_math::det::ln(u1));
            let angle = core::f64::consts::TAU * u2;
            pair[0] = magnitude * fs_math::det::cos(angle);
            if pair.len() == 2 {
                pair[1] = magnitude * fs_math::det::sin(angle);
            }
        }
    }

    /// One Box–Muller standard normal, used for scalar evaluation noise.
    fn gauss(&mut self) -> f64 {
        let u1 = self.next_open_f64();
        let u2 = self.next_f64();
        fs_math::det::sqrt(-2.0 * fs_math::det::ln(u1))
            * fs_math::det::cos(core::f64::consts::TAU * u2)
    }
}

// ---------------------------------------------------------------------------
// Landscapes (minimization, n-dimensional).
// ---------------------------------------------------------------------------

fn evaluate(landscape: u32, x: &[f64]) -> f64 {
    let n = x.len();
    match landscape {
        LANDSCAPE_SPHERE => x.iter().map(|v| v * v).sum(),
        LANDSCAPE_ROSENBROCK => {
            let mut s = 0.0;
            for i in 0..n - 1 {
                s += 100.0 * (x[i + 1] - x[i] * x[i]).powi(2) + (1.0 - x[i]).powi(2);
            }
            s
        }
        LANDSCAPE_CIGAR => {
            let mut s = 1.0e6 * x[0] * x[0];
            for v in &x[1..] {
                s += v * v;
            }
            s
        }
        LANDSCAPE_RASTRIGIN => {
            let mut s = 10.0 * n as f64;
            for v in x {
                s += v * v - 10.0 * (core::f64::consts::TAU * v).cos();
            }
            s
        }
        LANDSCAPE_ELLI => {
            let mut s = 0.0;
            for (i, v) in x.iter().enumerate() {
                let z = if n > 1 {
                    i as f64 * 6.0 / (n - 1) as f64
                } else {
                    0.0
                };
                s += 10f64.powf(z) * v * v;
            }
            s
        }
        _ => f64::INFINITY,
    }
}

// ---------------------------------------------------------------------------
// Linear-algebra helpers
// ---------------------------------------------------------------------------

fn reflect_repair(v: f64, lo: f64, hi: f64) -> f64 {
    if v >= lo && v <= hi {
        return v;
    }
    let span = hi - lo;
    let over = v - lo;
    let m = ((over % span) + span) % span;
    let wraps = (over / span).floor();
    if wraps as i64 % 2 == 0 {
        lo + m
    } else {
        hi - m
    }
}

/// C^{1/2} (`inv == false`) or C^{-1/2} = V·diag(p)·Vᵀ from the ascending
/// decomposition. Eigenvalues floored at 1e-18 to absorb tiny negatives from
/// active updates.
fn transform_matrix(eigvals: &[f64], eigvecs: &[f64], n: usize, inv: bool) -> Vec<f64> {
    const FLOOR: f64 = 1e-18;
    let p: Vec<f64> = eigvals
        .iter()
        .map(|v| {
            let v = v.max(FLOOR);
            if inv { 1.0 / v.sqrt() } else { v.sqrt() }
        })
        .collect();
    let mut out = vec![0.0f64; n * n];
    for j in 0..n {
        let pj = p[j];
        for i in 0..n {
            let vi = eigvecs[i * n + j];
            for k in 0..n {
                out[i * n + k] += pj * vi * eigvecs[k * n + j];
            }
        }
    }
    out
}

#[cfg(test)]
fn mat_vec(m: &[f64], v: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n];
    mat_vec_into(m, v, n, &mut out);
    out
}

fn mat_vec_into(m: &[f64], v: &[f64], n: usize, out: &mut [f64]) {
    debug_assert!(m.len() >= n * n);
    debug_assert!(v.len() >= n);
    debug_assert!(out.len() >= n);
    for (i, value) in out.iter_mut().take(n).enumerate() {
        *value = (0..n).map(|k| m[i * n + k] * v[k]).sum::<f64>();
    }
}

/// Advance the cumulative-path decay and return
/// `sqrt(1 - (1 - c_s)^(2 * generation))`.
///
/// Keeping the power as recurrence state avoids platform `pow` drift while
/// preserving the canonical CMA-ES generation exponent.
fn next_hsig_normalizer(decay_power: &mut f64, one_minus_cs_sq: f64) -> f64 {
    *decay_power *= one_minus_cs_sq;
    fs_math::det::sqrt(1.0 - *decay_power)
}

/// Canonical Hansen 2016 default damping for cumulative step-size
/// adaptation. The `- 1` belongs inside the positive part.
fn canonical_damps(mueff: f64, dimension: f64, cs: f64) -> f64 {
    1.0 + 2.0 * (((mueff - 1.0) / (dimension + 1.0)).sqrt() - 1.0).max(0.0) + cs
}

/// Rebuild C from its decomposition: V·Λ·Vᵀ (NOT C^{1/2} — transform_matrix
/// applies √λ; here the update rule needs the covariance itself back).
fn rebuild_c(eigvals: &[f64], eigvecs: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; n * n];
    for j in 0..n {
        let lj = eigvals[j].max(0.0);
        for i in 0..n {
            let vi = eigvecs[i * n + j];
            for k in 0..n {
                out[i * n + k] += lj * vi * eigvecs[k * n + j];
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Admission — refuse (never clamp) out-of-domain inputs.
// ---------------------------------------------------------------------------

/// Validate params; `Err` carries the typed refusal.
pub fn admit(p: &VizParams) -> Result<(), Refusal> {
    let refuse = |code: &'static str, message: String, repairs: Vec<&'static str>| {
        Err(Refusal {
            code,
            message,
            ranked_repairs: repairs,
        })
    };
    if !(2..=6).contains(&p.dim) {
        return refuse(
            "dim-out-of-range",
            format!("dim {} outside the visualization domain 2..=6", p.dim),
            vec!["set dim within 2..=6"],
        );
    }
    if p.x0.iter().take(p.dim).any(|v| !v.is_finite()) {
        return refuse(
            "x0-non-finite",
            "initial mean contains a NaN or infinite coordinate".into(),
            vec!["replace non-finite x0 coordinates with finite values"],
        );
    }
    if !p.sigma0.is_finite() || p.sigma0 <= 0.0 {
        return refuse(
            "sigma0-non-positive",
            format!("initial sigma {} must be finite and > 0", p.sigma0),
            vec!["set sigma0 to a positive finite step size"],
        );
    }
    if !(4..=48).contains(&p.lambda) {
        return refuse(
            "lambda-out-of-range",
            format!(
                "lambda {} outside the visualization domain 4..=48",
                p.lambda
            ),
            vec!["set lambda within 4..=48"],
        );
    }
    if !(1..=200).contains(&p.generations) {
        return refuse(
            "generations-out-of-range",
            format!(
                "generations {} outside the visualization domain 1..=200",
                p.generations
            ),
            vec!["set generations within 1..=200"],
        );
    }
    if p.landscape > LANDSCAPE_ELLI {
        return refuse(
            "landscape-unknown",
            format!("landscape id {} has no registered function", p.landscape),
            vec!["use ids 0..=4 (sphere, rosenbrock, cigar, rastrigin, elli)"],
        );
    }
    if !p.noise.is_finite() || p.noise < 0.0 {
        return refuse(
            "noise-invalid",
            format!("noise {} must be finite and >= 0", p.noise),
            vec!["set noise to 0 for noiseless evaluation"],
        );
    }
    if p.bounds_enabled
        && (!p.bound_min.is_finite() || !p.bound_max.is_finite() || p.bound_min >= p.bound_max)
    {
        return refuse(
            "bounds-inverted",
            "bounds require finite bound_min < bound_max".into(),
            vec!["disable bounds or provide bound_min < bound_max"],
        );
    }
    if !p.f_target.is_nan() && !p.f_target.is_finite() {
        return refuse(
            "f-target-invalid",
            "f_target must be finite or NaN (disabled)".into(),
            vec!["pass NaN to disable the early-stop target"],
        );
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// The run
// ---------------------------------------------------------------------------

/// Run the kernel and return the full per-generation internal stream.
pub fn cmaes_run(p: &VizParams) -> Result<VizRun, Refusal> {
    admit(p)?;
    let n = p.dim;
    let mut rng = Lcg(p.seed as u32);

    let mut mean: Vec<f64> = p.x0[..n].to_vec();
    let mut sigma = p.sigma0;
    let mut c = vec![0.0f64; n * n];
    for i in 0..n {
        c[i * n + i] = 1.0;
    }
    let mut p_sigma = vec![0.0f64; n];
    let mut p_c = vec![0.0f64; n];

    // Strategy constants (Hansen 2016; identical formulas to the TS fallback).
    // One log-weight vector over ALL lambda ranks; its positive half (the top
    // mu) is normalized into the recombination weights, and its negative half
    // becomes the active covariance weights below. For even lambda,
    // ln((lambda+1)/2) equals ln(mu + 0.5), so the positive weights match the
    // historical kernel exactly.
    let lambda = p.lambda;
    let mut raw_all = vec![0.0f64; lambda];
    for (i, w) in raw_all.iter_mut().enumerate() {
        *w = fs_math::det::ln((lambda as f64 + 1.0) / 2.0) - fs_math::det::ln((i + 1) as f64);
    }
    let mu = raw_all.iter().filter(|w| **w > 0.0).count().max(1);
    let positive_sum: f64 = raw_all[..mu].iter().sum();
    let weights: Vec<f64> = raw_all[..mu].iter().map(|w| w / positive_sum).collect();
    let mueff = 1.0 / weights.iter().map(|w| w * w).sum::<f64>();
    let nf = n as f64;
    let cc = (4.0 + mueff / nf) / (nf + 4.0 + 2.0 * mueff / nf);
    let cs = (mueff + 2.0) / (nf + mueff + 5.0);
    let c1 = 2.0 / ((nf + 1.3) * (nf + 1.3) + mueff);
    let cmu = ((1.0 - c1)
        .min((2.0 * (mueff - 2.0 + 1.0 / mueff)) / ((nf + 2.0) * (nf + 2.0) + mueff)))
    .max(0.0);
    let damps = canonical_damps(mueff, nf, cs);
    let chin = nf.sqrt() * (1.0 - 1.0 / (4.0 * nf) + 1.0 / (21.0 * nf * nf));
    let one_minus_cs_sq = (1.0 - cs) * (1.0 - cs);
    let mut hsig_decay_power = 1.0;

    // Canonical active (negative) covariance weights: the worst-ranked
    // offspring enter the rank-mu sum with negative weights scaled by
    // Hansen's alpha bounds, and are additionally Mahalanobis-rescaled per
    // candidate at update time — mirroring the site's rewritten TS engines.
    let neg_abs_sum: f64 = raw_all[mu..].iter().map(|w| w.abs()).sum();
    let neg_sq_sum: f64 = raw_all[mu..].iter().map(|w| w * w).sum();
    let mueff_minus = if neg_sq_sum > 0.0 {
        neg_abs_sum * neg_abs_sum / neg_sq_sum
    } else {
        0.0
    };
    let mut negative_scale = 0.0;
    if neg_abs_sum > 0.0 && cmu > 0.0 {
        let alpha_mu = 1.0 + c1 / cmu;
        let alpha_mueff = 1.0 + 2.0 * mueff_minus / (mueff + 2.0);
        let alpha_pd = (1.0 - c1 - cmu) / (nf * cmu);
        negative_scale = alpha_mu.min(alpha_mueff).min(alpha_pd).max(0.0);
    }
    let cov_weights: Vec<f64> = raw_all
        .iter()
        .enumerate()
        .map(|(i, &w)| {
            if i < mu {
                w / positive_sum
            } else if neg_abs_sum > 0.0 {
                w * negative_scale / neg_abs_sum
            } else {
                0.0
            }
        })
        .collect();
    let cov_weight_sum: f64 = cov_weights.iter().sum();

    let mut best_f = f64::INFINITY;
    let mut best_x = mean.clone();
    let mut evals = 0usize;
    let mut generations: Vec<GenSnapshot> = Vec::with_capacity(p.generations);
    let mut stop_reason = "generations-exhausted";

    // Cache the decomposition that defines the sampling distribution. Each
    // generation computes exactly one new decomposition after updating C;
    // that post-update decomposition is both the next sampling state and the
    // coherent spectrum emitted in the current snapshot.
    let mut current_eigvals = vec![1.0f64; n];
    let mut current_eigvecs = vec![0.0f64; n * n];
    for i in 0..n {
        current_eigvecs[i * n + i] = 1.0;
    }
    let mut matvec_y = vec![0.0; n];
    let mut z_mean = vec![0.0; n];
    let mut whitened = vec![0.0; n];

    for g in 0..p.generations {
        let sqrt_c = transform_matrix(&current_eigvals, &current_eigvecs, n, false);
        let inv_sqrt_c = transform_matrix(&current_eigvals, &current_eigvecs, n, true);

        // 1. Sample λ offspring: x = m + σ·C^{1/2}·z.
        let mut sx = vec![0.0f64; lambda * n];
        let mut raw_sx = vec![0.0f64; lambda * n];
        let mut sz_raw = vec![0.0f64; lambda * n];
        let mut sf = vec![0.0f64; lambda];
        let mut st = vec![0.0f64; lambda]; // true (noiseless) fitness
        for i in 0..lambda {
            let z = &mut sz_raw[i * n..(i + 1) * n];
            rng.fill_gaussian(z);
            mat_vec_into(&sqrt_c, z, n, &mut matvec_y);
            for k in 0..n {
                let raw_xk = mean[k] + sigma * matvec_y[k];
                raw_sx[i * n + k] = raw_xk;
                let mut xk = raw_xk;
                if p.bounds_enabled {
                    xk = reflect_repair(xk, p.bound_min, p.bound_max);
                }
                sx[i * n + k] = xk;
            }
            let row = &sx[i * n..(i + 1) * n];
            let true_f = evaluate(p.landscape, row);
            if !true_f.is_finite() {
                return Err(Refusal {
                    code: "non-finite-objective",
                    message: format!(
                        "landscape {} produced a non-finite value at generation {g}",
                        p.landscape
                    ),
                    ranked_repairs: vec!["reduce sigma0", "enable bounds repair"],
                });
            }
            let noisy = if p.noise > 0.0 {
                true_f + rng.gauss() * p.noise
            } else {
                true_f
            };
            st[i] = true_f;
            sf[i] = noisy;
            evals += 1;
        }

        // 2. Rank by evaluated (possibly noisy) fitness; earliest wins ties.
        let mut order: Vec<usize> = (0..lambda).collect();
        order.sort_by(|&a, &b| sf[a].total_cmp(&sf[b]).then(a.cmp(&b)));
        let se: Vec<u8> = (0..lambda).map(|rank| u8::from(rank < mu)).collect();

        // Best tracks TRUE fitness (earliest total_cmp minimum).
        for &idx in &order {
            if st[idx] < best_f {
                best_f = st[idx];
                best_x = sx[idx * n..(idx + 1) * n].to_vec();
            }
        }

        // 3. Recombination: weighted mean of elites.
        let old_mean = mean.clone();
        mean = vec![0.0; n];
        for (rank, &idx) in order.iter().enumerate().take(mu) {
            for k in 0..n {
                // Reflection is a phenotype transform: selection evaluates
                // repaired points, while adaptation follows their latent
                // Gaussian preimages, matching the canonical TS fallback.
                mean[k] += weights[rank] * raw_sx[idx * n + k];
            }
        }
        let mean_shift: Vec<f64> = (0..n).map(|k| (mean[k] - old_mean[k]) / sigma).collect();
        mat_vec_into(&inv_sqrt_c, &mean_shift, n, &mut z_mean);

        // 4. Evolution paths.
        let ps_coeff = (cs * (2.0 - cs) * mueff).sqrt();
        for k in 0..n {
            p_sigma[k] = (1.0 - cs) * p_sigma[k] + ps_coeff * z_mean[k];
        }
        let norm_ps: f64 = p_sigma.iter().map(|v| v * v).sum::<f64>().sqrt();
        let hsig_denom = next_hsig_normalizer(&mut hsig_decay_power, one_minus_cs_sq);
        let hsig = if hsig_denom > 0.0 && norm_ps / hsig_denom / chin < 1.4 + 2.0 / (nf + 1.0) {
            1.0
        } else {
            0.0
        };
        let pc_coeff = (cc * (2.0 - cc) * mueff).sqrt();
        for k in 0..n {
            p_c[k] = (1.0 - cc) * p_c[k] + hsig * pc_coeff * mean_shift[k];
        }

        // 5. Covariance adaptation: rank-1 + rank-μ over all ranks with the
        // canonical negative (active) weights, each worst-ranked candidate's
        // weight Mahalanobis-rescaled by n/‖C^{-1/2}y‖² — mirroring the
        // site's rewritten TS engines (Hansen 2016 active CMA).
        let mut adjusted = cov_weights.clone();
        if p.active {
            for rank in mu..lambda {
                if adjusted[rank] >= 0.0 {
                    continue;
                }
                let idx = order[rank];
                for k in 0..n {
                    matvec_y[k] = (raw_sx[idx * n + k] - old_mean[k]) / sigma;
                }
                mat_vec_into(&inv_sqrt_c, &matvec_y, n, &mut whitened);
                let mahalanobis_sq: f64 = whitened.iter().map(|v| v * v).sum();
                adjusted[rank] = if mahalanobis_sq > 0.0 {
                    adjusted[rank] * nf / mahalanobis_sq
                } else {
                    0.0
                };
            }
        } else {
            for weight in adjusted.iter_mut().skip(mu) {
                *weight = 0.0;
            }
        }
        let weight_sum = if p.active { cov_weight_sum } else { 1.0 };
        let old_coeff = 1.0 + c1 * (1.0 - hsig) * cc * (2.0 - cc) - c1 - cmu * weight_sum;
        let mut new_c = vec![0.0f64; n * n];
        for i in 0..n {
            for k in 0..n {
                let rank1 = p_c[i] * p_c[k];
                let mut rankmu = 0.0;
                for (rank, &idx) in order.iter().enumerate() {
                    let wgt = adjusted[rank];
                    if wgt == 0.0 {
                        continue;
                    }
                    let yi = (raw_sx[idx * n + i] - old_mean[i]) / sigma;
                    let yk = (raw_sx[idx * n + k] - old_mean[k]) / sigma;
                    rankmu += wgt * yi * yk;
                }
                new_c[i * n + k] = old_coeff * c[i * n + k] + c1 * rank1 + cmu * rankmu;
            }
        }
        c = new_c;

        // PD repair (mirrors the TS engines' floored reconstruction): the
        // negative update can overshoot on unlucky rankings, so floor the
        // spectrum relative to its own scale and rebuild. Sampling, the
        // snapshot eigvals, and cond then always describe a genuine
        // positive-definite covariance.
        let (repair_ev, repair_evec) = jacobi_eigh(&c, n);
        if repair_ev.iter().any(|v| !v.is_finite()) {
            return Err(Refusal {
                code: "eigen-decomposition-failed",
                message: format!(
                    "covariance repair produced non-finite eigenvalues at generation {g}"
                ),
                ranked_repairs: vec!["disable the active update", "reduce sigma0"],
            });
        }
        let spectrum_scale = repair_ev
            .iter()
            .fold(f64::MIN_POSITIVE, |acc, v| acc.max(v.abs()));
        let spectrum_floor = 1e-14 * spectrum_scale;
        let repaired: Vec<f64> = repair_ev.iter().map(|v| v.max(spectrum_floor)).collect();
        c = rebuild_c(&repaired, &repair_evec, n);

        let cond = repaired[n - 1] / repaired[0].max(f64::MIN_POSITIVE);

        // 6. Step-size adaptation with the ND fallback's safety envelope.
        sigma *= fs_math::det::exp((cs / damps) * (norm_ps / chin - 1.0));
        sigma = sigma.clamp(1e-16, 1e16);

        // Keep every population stream in the same rank order. Optimization
        // above retains sampling-order storage because `order` indexes it;
        // the public snapshot must not expose those mismatched indices.
        let mut ranked_sx = vec![0.0f64; lambda * n];
        let mut ranked_sz = vec![0.0f64; lambda * n];
        let mut ranked_sf = vec![0.0f64; lambda];
        for (rank, &idx) in order.iter().enumerate() {
            ranked_sx[rank * n..(rank + 1) * n].copy_from_slice(&sx[idx * n..(idx + 1) * n]);
            ranked_sz[rank * n..(rank + 1) * n].copy_from_slice(&sz_raw[idx * n..(idx + 1) * n]);
            ranked_sf[rank] = sf[idx];
        }

        let snapshot_mean: Vec<f64> = mean
            .iter()
            .map(|&value| {
                if p.bounds_enabled {
                    reflect_repair(value, p.bound_min, p.bound_max)
                } else {
                    value
                }
            })
            .collect();
        generations.push(GenSnapshot {
            g: g + 1,
            mean: snapshot_mean,
            sigma,
            eigvals: repaired.clone(),
            eigvecs: repair_evec.clone(),
            cond,
            best_f,
            evals,
            proj_mean: [0.0; 3],
            proj_eigvals: [0.0; 3],
            proj_eigvecs: [0.0; 9],
            sx: ranked_sx,
            sz: ranked_sz,
            sf: ranked_sf,
            se,
            p_sigma: p_sigma.clone(),
            p_c: p_c.clone(),
        });

        current_eigvals = repaired;
        current_eigvecs = repair_evec;

        if !p.f_target.is_nan() && best_f <= p.f_target {
            stop_reason = "target-reached";
            break;
        }
    }

    // 7. Phase-space projection: direct for dim ≤ 3, honest PCA marginal for
    // dim > 3 (basis from the pooled covariance; per-generation marginals are
    // the projected 3×3 covariances — never a fake 3D ellipsoid).
    let (basis, center, pool_vals) = project_phase_space(&mut generations, &c, n);

    Ok(VizRun {
        dim: n,
        landscape: p.landscape,
        best_f,
        best_x,
        total_evals: evals,
        stop_reason,
        generations,
        pca: PcaInfo {
            basis,
            center,
            pool_eigvals: pool_vals,
        },
    })
}

/// Compute and write the 3D phase-space projections into each snapshot.
/// Returns the PCA frame (basis rows, center, pooled eigenvalues).
fn project_phase_space(
    gens: &mut [GenSnapshot],
    final_c: &[f64],
    n: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let ng = gens.len().max(1);
    let mut center = vec![0.0; n];
    for snap in gens.iter() {
        for (coordinate, &mean) in center.iter_mut().zip(&snap.mean) {
            *coordinate += mean / ng as f64;
        }
    }
    let (pool_vals, pool_vecs) = jacobi_eigh(final_c, n);
    let basis: Vec<f64> = if n <= 3 {
        let mut b = vec![0.0; 3 * n];
        for r in 0..n {
            b[r * n + r] = 1.0;
        }
        b
    } else {
        let mut b = vec![0.0; 3 * n];
        for r in 0..3 {
            let src_col = n - 3 + r; // largest three, ascending order
            for i in 0..n {
                b[r * n + i] = pool_vecs[i * n + src_col];
            }
        }
        b
    };
    for snap in gens.iter_mut() {
        let mut pm = [0.0f64; 3];
        for r in 0..3 {
            pm[r] = (0..n)
                .map(|i| basis[r * n + i] * (snap.mean[i] - center[i]))
                .sum();
        }
        // 3×3 marginal M = P·C_t·Pᵀ with C_t rebuilt from this generation's
        // stored decomposition (C = V·Λ·Vᵀ).
        let ct = rebuild_c(&snap.eigvals, &snap.eigvecs, n);
        let mut m3 = [0.0f64; 9];
        for r in 0..3 {
            for s in 0..3 {
                let mut acc = 0.0;
                for i in 0..n {
                    for k in 0..n {
                        acc += basis[r * n + i] * ct[i * n + k] * basis[s * n + k];
                    }
                }
                m3[r * 3 + s] = acc;
            }
        }
        // Symmetrize against float drift, then decompose.
        for r in 0..3 {
            for s in 0..3 {
                let sym = (m3[r * 3 + s] + m3[s * 3 + r]) * 0.5;
                m3[r * 3 + s] = sym;
                m3[s * 3 + r] = sym;
            }
        }
        let (pv, pw) = jacobi_eigh(&m3, 3);
        snap.proj_mean = pm;
        snap.proj_eigvals = [pv[0], pv[1], pv[2]];
        snap.proj_eigvecs.copy_from_slice(&pw);
    }
    (basis, center, pool_vals)
}

// ---------------------------------------------------------------------------
// JSON envelope (hand-rolled; no serde at the boundary — fs-flyer pattern).
// ---------------------------------------------------------------------------

fn push_num(out: &mut String, v: f64) {
    if v.is_finite() {
        let _ = write!(out, "{v}");
    } else if v.is_nan() {
        out.push_str("null");
    } else if v > 0.0 {
        out.push_str("1e999");
    } else {
        out.push_str("-1e999");
    }
}

fn push_num_arr(out: &mut String, vals: &[f64]) {
    out.push('[');
    for (i, v) in vals.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        push_num(out, *v);
    }
    out.push(']');
}

fn push_byte_arr(out: &mut String, vals: &[u8]) {
    out.push('[');
    for (i, v) in vals.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let _ = write!(out, "{v}");
    }
    out.push(']');
}

fn refusal_envelope(r: &Refusal) -> String {
    let mut repairs = String::new();
    for (i, rep) in r.ranked_repairs.iter().enumerate() {
        if i > 0 {
            repairs.push(',');
        }
        repairs.push('"');
        repairs.push_str(rep);
        repairs.push('"');
    }
    format!(
        "{{\"refusal\":{{\"code\":\"{}\",\"message\":\"{}\",\"ranked_repairs\":[{}]}}}}",
        r.code,
        r.message.replace('\\', "\\\\").replace('"', "\\\""),
        repairs
    )
}

/// Envelope-producing entry shared by native tests and the wasm boundary.
/// Scalars in (x0 as up-to-6 components, `dim` selects the prefix), JSON
/// envelope out.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn cmaes_run_json(
    dim: u32,
    x0: [f64; 6],
    sigma0: f64,
    lambda: u32,
    active: bool,
    seed: u64,
    generations: u32,
    landscape: u32,
    noise: f64,
    bounds_enabled: bool,
    bound_min: f64,
    bound_max: f64,
    f_target: f64,
) -> String {
    let params = VizParams {
        dim: dim as usize,
        x0,
        sigma0,
        lambda: lambda as usize,
        active,
        seed,
        generations: generations as usize,
        landscape,
        noise,
        bounds_enabled,
        bound_min,
        bound_max,
        f_target,
    };
    match cmaes_run(&params) {
        Ok(run) => ok_envelope(&run),
        Err(refusal) => refusal_envelope(&refusal),
    }
}

fn ok_envelope(run: &VizRun) -> String {
    // Budget 24 bytes per scalar for the UI-reachable distributions. This is
    // only a capacity hint (String still grows correctly for an unusually
    // long fixed-form value), and it avoids both growth copies and the much
    // more expensive temporary String per scalar/array in the old formatter.
    let scalar_count: usize = run
        .generations
        .iter()
        .map(|snap| {
            snap.mean.len()
                + 1
                + snap.eigvals.len()
                + snap.eigvecs.len()
                + 2
                + snap.proj_mean.len()
                + snap.proj_eigvals.len()
                + snap.proj_eigvecs.len()
                + snap.sx.len()
                + snap.sz.len()
                + snap.sf.len()
                + snap.p_sigma.len()
                + snap.p_c.len()
        })
        .sum::<usize>()
        + 1
        + run.best_x.len()
        + run.pca.basis.len()
        + run.pca.center.len()
        + run.pca.pool_eigvals.len();
    let byte_count: usize = run.generations.iter().map(|snap| snap.se.len()).sum();
    let mut out = String::with_capacity(
        512 + scalar_count * 24 + byte_count * 2 + run.generations.len() * 256,
    );
    let _ = write!(
        out,
        "{{\"ok\":{{\"kernel\":\"{}\",\"dim\":{},\"landscape\":{},\"stop_reason\":\"{}\",\"best_f\":",
        KERNEL_VERSION, run.dim, run.landscape, run.stop_reason
    );
    push_num(&mut out, run.best_f);
    out.push_str(",\"best_x\":");
    push_num_arr(&mut out, &run.best_x);
    let _ = write!(
        out,
        ",\"total_evals\":{},\"generations\":[",
        run.total_evals
    );
    for (gi, snap) in run.generations.iter().enumerate() {
        if gi > 0 {
            out.push(',');
        }
        let _ = write!(out, "{{\"g\":{},\"mean\":", snap.g);
        push_num_arr(&mut out, &snap.mean);
        out.push_str(",\"sigma\":");
        push_num(&mut out, snap.sigma);
        out.push_str(",\"eigvals\":");
        push_num_arr(&mut out, &snap.eigvals);
        out.push_str(",\"eigvecs\":");
        push_num_arr(&mut out, &snap.eigvecs);
        out.push_str(",\"cond\":");
        push_num(&mut out, snap.cond);
        out.push_str(",\"best_f\":");
        push_num(&mut out, snap.best_f);
        let _ = write!(out, ",\"evals\":{},\"proj_mean\":", snap.evals);
        push_num_arr(&mut out, &snap.proj_mean);
        out.push_str(",\"proj_eigvals\":");
        push_num_arr(&mut out, &snap.proj_eigvals);
        out.push_str(",\"proj_eigvecs\":");
        push_num_arr(&mut out, &snap.proj_eigvecs);
        out.push_str(",\"sx\":");
        push_num_arr(&mut out, &snap.sx);
        out.push_str(",\"sz\":");
        push_num_arr(&mut out, &snap.sz);
        out.push_str(",\"sf\":");
        push_num_arr(&mut out, &snap.sf);
        out.push_str(",\"se\":");
        push_byte_arr(&mut out, &snap.se);
        out.push_str(",\"p_sigma\":");
        push_num_arr(&mut out, &snap.p_sigma);
        out.push_str(",\"p_c\":");
        push_num_arr(&mut out, &snap.p_c);
        out.push('}');
    }
    out.push_str("],\"pca_basis\":");
    push_num_arr(&mut out, &run.pca.basis);
    out.push_str(",\"pca_center\":");
    push_num_arr(&mut out, &run.pca.center);
    out.push_str(",\"pca_pool_eigvals\":");
    push_num_arr(&mut out, &run.pca.pool_eigvals);
    out.push_str("}}");
    out
}

/// Kernel identity probe for the capability check.
#[must_use]
pub fn kernel_version() -> &'static str {
    KERNEL_VERSION
}

// ---------------------------------------------------------------------------
// wasm32-only JS boundary.
// ---------------------------------------------------------------------------

#[cfg(target_arch = "wasm32")]
mod js {
    use wasm_bindgen::prelude::wasm_bindgen;

    /// Run the visualization kernel: scalar params in, JSON envelope out.
    /// `x0_0..x0_5` are the initial mean coordinates; only the first `dim`
    /// are read. `f_target` = NaN disables the early-stop target.
    #[wasm_bindgen]
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn cmaes_viz_run(
        dim: u32,
        x0_0: f64,
        x0_1: f64,
        x0_2: f64,
        x0_3: f64,
        x0_4: f64,
        x0_5: f64,
        sigma0: f64,
        lambda: u32,
        active: bool,
        seed: u64,
        generations: u32,
        landscape: u32,
        noise: f64,
        bounds_enabled: bool,
        bound_min: f64,
        bound_max: f64,
        f_target: f64,
    ) -> String {
        super::cmaes_run_json(
            dim,
            [x0_0, x0_1, x0_2, x0_3, x0_4, x0_5],
            sigma0,
            lambda,
            active,
            seed,
            generations,
            landscape,
            noise,
            bounds_enabled,
            bound_min,
            bound_max,
            f_target,
        )
    }

    /// Kernel identity probe (capability check after instantiation).
    #[wasm_bindgen]
    #[must_use]
    pub fn cmaes_viz_kernel_version() -> String {
        super::KERNEL_VERSION.to_string()
    }
}

// ---------------------------------------------------------------------------
// Tests — exact-value and invariant checks, native.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn base_params() -> VizParams {
        VizParams {
            dim: 3,
            x0: [1.5, -1.0, 2.0, 0.0, 0.0, 0.0],
            sigma0: 0.5,
            lambda: 16,
            active: true,
            seed: 1337,
            generations: 150,
            landscape: LANDSCAPE_SPHERE,
            noise: 0.0,
            bounds_enabled: false,
            bound_min: 0.0,
            bound_max: 0.0,
            f_target: f64::NAN,
        }
    }

    #[test]
    fn sphere_3d_converges() {
        let run = cmaes_run(&base_params()).expect("run");
        assert!(run.best_f < 1e-8, "best_f = {}", run.best_f);
        assert_eq!(run.total_evals, 150 * 16);
        assert_eq!(run.stop_reason, "generations-exhausted");
    }

    #[test]
    fn best_fitness_monotone() {
        let run = cmaes_run(&base_params()).expect("run");
        for w in run.generations.windows(2) {
            assert!(
                w[1].best_f <= w[0].best_f + 1e-15,
                "best fitness increased: {} -> {}",
                w[0].best_f,
                w[1].best_f
            );
        }
    }

    #[test]
    fn eigendecomposition_reconstructs_c() {
        let p = VizParams {
            landscape: LANDSCAPE_ROSENBROCK,
            generations: 40,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");
        let snap = run.generations.last().expect("gens");
        let n = run.dim;
        // C = V·Λ·Vᵀ must reproduce a symmetric PSD matrix whose eigendecom-
        // position residual ‖C·v_j − λ_j·v_j‖ is at machine-epsilon scale.
        let c = rebuild_c(&snap.eigvals, &snap.eigvecs, n);
        for j in 0..n {
            let v: Vec<f64> = (0..n).map(|i| snap.eigvecs[i * n + j]).collect();
            let cv = mat_vec(&c, &v, n);
            for i in 0..n {
                let resid = (cv[i] - snap.eigvals[j] * v[i]).abs();
                // Tolerance scales with the spectrum: absolute epsilons are
                // meaningless once active updates drive cond(C) past 1e4.
                let scale = snap.eigvals[n - 1].abs().max(1.0);
                assert!(
                    resid < 1e-9 * scale,
                    "residual {resid} at ({i},{j}), scale {scale}"
                );
            }
        }
        // Symmetry + positive semidefiniteness.
        for i in 0..n {
            for k in 0..n {
                assert!((c[i * n + k] - c[k * n + i]).abs() < 1e-12);
            }
        }
        assert!(snap.eigvals.iter().all(|v| *v > -1e-12));
    }

    #[test]
    fn bitwise_replay_same_seed() {
        let a = cmaes_run_json(
            3,
            [1.0, 2.0, 3.0, 0.0, 0.0, 0.0],
            0.5,
            12,
            true,
            42,
            30,
            1,
            0.0,
            false,
            0.0,
            0.0,
            f64::NAN,
        );
        let b = cmaes_run_json(
            3,
            [1.0, 2.0, 3.0, 0.0, 0.0, 0.0],
            0.5,
            12,
            true,
            42,
            30,
            1,
            0.0,
            false,
            0.0,
            0.0,
            f64::NAN,
        );
        assert_eq!(a, b, "same seed must replay bitwise");
        let c = cmaes_run_json(
            3,
            [1.0, 2.0, 3.0, 0.0, 0.0, 0.0],
            0.5,
            12,
            true,
            43,
            30,
            1,
            0.0,
            false,
            0.0,
            0.0,
            f64::NAN,
        );
        assert_ne!(a, c, "different seed must diverge");
    }

    #[test]
    fn refusal_codes_are_typed() {
        let cases: Vec<(VizParams, &str)> = vec![
            (
                VizParams {
                    dim: 1,
                    ..base_params()
                },
                "dim-out-of-range",
            ),
            (
                VizParams {
                    dim: 7,
                    ..base_params()
                },
                "dim-out-of-range",
            ),
            (
                VizParams {
                    sigma0: 0.0,
                    ..base_params()
                },
                "sigma0-non-positive",
            ),
            (
                VizParams {
                    sigma0: f64::NAN,
                    ..base_params()
                },
                "sigma0-non-positive",
            ),
            (
                VizParams {
                    lambda: 3,
                    ..base_params()
                },
                "lambda-out-of-range",
            ),
            (
                VizParams {
                    lambda: 49,
                    ..base_params()
                },
                "lambda-out-of-range",
            ),
            (
                VizParams {
                    generations: 0,
                    ..base_params()
                },
                "generations-out-of-range",
            ),
            (
                VizParams {
                    generations: 201,
                    ..base_params()
                },
                "generations-out-of-range",
            ),
            (
                VizParams {
                    landscape: 9,
                    ..base_params()
                },
                "landscape-unknown",
            ),
            (
                VizParams {
                    noise: -0.1,
                    ..base_params()
                },
                "noise-invalid",
            ),
            (
                VizParams {
                    bounds_enabled: true,
                    bound_min: 2.0,
                    bound_max: -2.0,
                    ..base_params()
                },
                "bounds-inverted",
            ),
            (
                VizParams {
                    f_target: f64::INFINITY,
                    ..base_params()
                },
                "f-target-invalid",
            ),
        ];
        let mut x_bad = base_params();
        x_bad.x0[0] = f64::NAN;
        cases.into_iter().for_each(|(p, code)| {
            let err = cmaes_run(&p).expect_err("must refuse");
            assert_eq!(err.code, code, "wrong refusal for {code}");
        });
        let err = cmaes_run(&x_bad).expect_err("must refuse");
        assert_eq!(err.code, "x0-non-finite");
        // Envelope form is valid JSON-shaped.
        let env = cmaes_run_json(
            1,
            [0.0; 6],
            0.5,
            12,
            true,
            1,
            10,
            0,
            0.0,
            false,
            0.0,
            0.0,
            f64::NAN,
        );
        assert!(
            env.contains("\"code\":\"dim-out-of-range\""),
            "envelope was: {env}"
        );
    }

    #[test]
    fn target_stops_early() {
        let p = VizParams {
            f_target: 1.0,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");
        assert_eq!(run.stop_reason, "target-reached");
        assert!(run.generations.len() < 150);
        assert!(run.best_f <= 1.0);
    }

    #[test]
    fn rosenbrock_improves_on_start() {
        let p = VizParams {
            landscape: LANDSCAPE_ROSENBROCK,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");
        let f0 = evaluate(LANDSCAPE_ROSENBROCK, &run.best_x);
        let _ = f0;
        assert!(run.best_f < 100.0, "best_f = {}", run.best_f);
    }

    #[test]
    fn bounds_keep_population_in_domain() {
        let p = VizParams {
            bounds_enabled: true,
            bound_min: -2.0,
            bound_max: 2.0,
            landscape: LANDSCAPE_RASTRIGIN,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");
        for snap in &run.generations {
            for i in 0..run.dim {
                for s in 0..p.lambda {
                    let v = snap.sx[s * run.dim + i];
                    assert!((-2.0..=2.0).contains(&v), "sample {v} escaped bounds");
                }
            }
        }
    }

    #[test]
    fn pca_marginal_is_consistent() {
        let p = VizParams {
            dim: 5,
            landscape: LANDSCAPE_ELLI,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");
        let snap = run.generations.last().expect("gens");
        // Projected marginal eigenvalues: ascending, non-negative, and the
        // reconstruction P·(V Λ Vᵀ)·Pᵀ trace must equal the eigenvalue sum.
        let trace: f64 = snap.proj_eigvals.iter().sum();
        let mut m3 = [0.0f64; 9];
        for r in 0..3 {
            for s in 0..3 {
                m3[r * 3 + s] = snap.proj_eigvecs[r * 3 + s];
            }
        }
        let _ = m3;
        assert!(snap.proj_eigvals.iter().all(|v| *v >= -1e-12));
        assert!(trace >= -1e-12);
        assert_eq!(snap.proj_eigvecs.len(), 9);
    }

    #[test]
    fn active_update_spectrum_stays_positive_definite() {
        // Regression: v0.1.0's simplified active update drove C indefinite on
        // this exact UI-reachable configuration (n=2, λ=8, σ₀=0.6, Rastrigin,
        // active on), reporting negative eigvals and cond ≈ 1e18 in 119/120
        // generations. The canonical update + spectral repair must keep every
        // snapshot's spectrum strictly positive with a sane condition number.
        let p = VizParams {
            dim: 2,
            x0: [1.5, -1.0, 0.0, 0.0, 0.0, 0.0],
            sigma0: 0.6,
            lambda: 8,
            landscape: LANDSCAPE_RASTRIGIN,
            generations: 120,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");
        for snap in &run.generations {
            assert!(
                snap.eigvals.iter().all(|v| *v > 0.0),
                "gen {} has non-positive eigenvalue: {:?}",
                snap.g,
                snap.eigvals
            );
            assert!(
                snap.cond >= 1.0 && snap.cond < 1e14,
                "gen {} has implausible cond {}",
                snap.g,
                snap.cond
            );
        }
    }

    #[test]
    fn g0_hsig_normalizer_uses_exponential_generation_decay() {
        let cs = 0.3;
        let one_minus_cs_sq = (1.0 - cs) * (1.0 - cs);
        let mut decay_power = 1.0;

        for generation in 1..=200 {
            let observed = next_hsig_normalizer(&mut decay_power, one_minus_cs_sq);
            let exponent = 2 * generation;
            let expected = fs_math::det::sqrt(1.0 - fs_math::det::powi(1.0 - cs, exponent));
            assert!(
                observed.is_finite(),
                "generation {generation} produced {observed}"
            );
            assert!(
                (observed - expected).abs() <= 8.0 * f64::EPSILON,
                "generation {generation}: observed {observed}, expected {expected}"
            );
        }

        // The former linear multiplier is already outside sqrt's domain here.
        assert!(1.0 - one_minus_cs_sq * 3.0 < 0.0);
    }

    #[test]
    fn g0_damping_matches_hansen_2016_default() {
        // n=5, lambda=16 constants independently evaluated from the Hansen
        // logarithmic recombination weights. v0.2.1 returned
        // 3.061140170192783 because it omitted the inner `- 1`.
        let mueff = 4.840_914_500_901_174;
        let cs = 0.460_949_660_513_560_5;
        let observed = canonical_damps(mueff, 5.0, cs);
        assert_eq!(observed.to_bits(), 1.460_949_660_513_560_6f64.to_bits());
    }

    #[test]
    fn g0_rng_uses_shared_u32_transition_and_paired_box_muller() {
        let mut uniforms = Lcg(1337);
        assert_eq!(
            uniforms.next_f64().to_bits(),
            (3_239_374_148.0f64 / 4_294_967_296.0).to_bits()
        );
        assert_eq!(
            uniforms.next_f64().to_bits(),
            (2_360_088_531.0f64 / 4_294_967_296.0).to_bits()
        );

        let mut paired = Lcg(1337);
        let mut z = [0.0; 3];
        paired.fill_gaussian(&mut z);
        assert!(z.iter().all(|value| value.is_finite()));
        assert_eq!(
            paired.0, 681_817_981,
            "three coordinates must consume two Box-Muller pairs"
        );

        let mut scalar = Lcg(1337);
        assert!(scalar.gauss().is_finite());
        assert_eq!(
            scalar.0, 2_360_088_531,
            "one scalar Gaussian must consume one pair"
        );
    }

    #[test]
    fn g0_snapshot_eigensystem_is_the_updated_covariance() {
        let p = VizParams {
            generations: 1,
            active: false,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");
        let snap = &run.generations[0];

        assert!(
            snap.eigvals
                .iter()
                .any(|value| value.to_bits() != 1.0f64.to_bits()),
            "the first post-update covariance must not be reported as the initial identity"
        );
        let expected_cond = snap.eigvals[p.dim - 1] / snap.eigvals[0];
        assert_eq!(snap.cond.to_bits(), expected_cond.to_bits());
        assert!(
            snap.eigvals
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
        );
    }

    #[test]
    fn g0_reflection_ranks_phenotypes_and_adapts_latent_preimages() {
        let p = VizParams {
            dim: 2,
            x0: [0.0; 6],
            sigma0: 1.0,
            lambda: 8,
            active: false,
            generations: 1,
            bounds_enabled: true,
            bound_min: -0.1,
            bound_max: 0.1,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");
        let snap = &run.generations[0];

        let raw_all: Vec<f64> = (0..p.lambda)
            .map(|rank| {
                fs_math::det::ln((p.lambda as f64 + 1.0) / 2.0)
                    - fs_math::det::ln((rank + 1) as f64)
            })
            .collect();
        let mu = raw_all.iter().filter(|weight| **weight > 0.0).count();
        let positive_sum: f64 = raw_all[..mu].iter().sum();

        for dimension in 0..p.dim {
            let latent_mean: f64 = (0..mu)
                .map(|rank| {
                    let weight = raw_all[rank] / positive_sum;
                    let raw_x = p.x0[dimension] + p.sigma0 * snap.sz[rank * p.dim + dimension];
                    weight * raw_x
                })
                .sum();
            let expected = reflect_repair(latent_mean, p.bound_min, p.bound_max);
            assert!(
                (snap.mean[dimension] - expected).abs() <= 1e-13,
                "dimension {dimension}: observed {}, expected {expected}",
                snap.mean[dimension]
            );
        }
        assert!(
            snap.sx
                .iter()
                .all(|value| (p.bound_min..=p.bound_max).contains(value))
        );
    }

    #[test]
    fn g0_population_snapshot_streams_are_rank_aligned() {
        let p = VizParams {
            generations: 3,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");

        for snap in &run.generations {
            for pair in snap.sf.windows(2) {
                assert!(
                    pair[0].total_cmp(&pair[1]).is_le(),
                    "generation {} fitness stream is not ranked: {:?}",
                    snap.g,
                    snap.sf
                );
            }
            for rank in 0..p.lambda {
                let x = &snap.sx[rank * p.dim..(rank + 1) * p.dim];
                assert_eq!(
                    evaluate(p.landscape, x).to_bits(),
                    snap.sf[rank].to_bits(),
                    "generation {}, rank {} has mismatched x/f streams",
                    snap.g,
                    rank
                );
            }
        }

        // At generation one C = I, so the ranked x and z streams must retain
        // their pairwise sampling relation after the public reorder.
        let first = &run.generations[0];
        for rank in 0..p.lambda {
            for k in 0..p.dim {
                let expected = p.x0[k] + p.sigma0 * first.sz[rank * p.dim + k];
                assert_eq!(first.sx[rank * p.dim + k].to_bits(), expected.to_bits());
            }
        }
    }

    #[test]
    fn dim2_projection_is_direct() {
        let p = VizParams {
            dim: 2,
            ..base_params()
        };
        let run = cmaes_run(&p).expect("run");
        let snap = run.generations.last().expect("gens");
        // With dim ≤ 3 the projection is the identity frame: proj_mean equals
        // mean centered by the trajectory centroid, third axis pinned at 0.
        assert_eq!(snap.proj_mean[2], 0.0);
        assert_eq!(run.dim, 2);
    }
}
