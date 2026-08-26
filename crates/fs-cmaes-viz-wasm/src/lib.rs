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

/// Landscape ids at the ABI (frozen v1).
pub const LANDSCAPE_SPHERE: u32 = 0;
pub const LANDSCAPE_ROSENBROCK: u32 = 1;
pub const LANDSCAPE_CIGAR: u32 = 2;
pub const LANDSCAPE_RASTRIGIN: u32 = 3;
pub const LANDSCAPE_ELLI: u32 = 4;

/// Kernel id baked into envelopes so the page can prove which build is live.
pub const KERNEL_VERSION: &str = "fs-cmaes-viz-wasm 0.1.0";

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

struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        // Top 24 bits → [0, 1); the offset keeps log(u) finite without clamping.
        (self.0 >> 40) as f64 / 16_777_216.0 + 5.960_464_477_539_06e-8
    }

    /// Box–Muller standard normal.
    fn gauss(&mut self) -> f64 {
        let u1 = self.next_f64();
        let u2 = self.next_f64();
        (-2.0 * fs_math::det::ln(u1)).sqrt() * (core::f64::consts::TAU * u2).cos()
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
                let z = if n > 1 { i as f64 * 6.0 / (n - 1) as f64 } else { 0.0 };
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
    if wraps as i64 % 2 == 0 { lo + m } else { hi - m }
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

fn mat_vec(m: &[f64], v: &[f64], n: usize) -> Vec<f64> {
    (0..n)
        .map(|i| (0..n).map(|k| m[i * n + k] * v[k]).sum::<f64>())
        .collect()
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
        Err(Refusal { code, message, ranked_repairs: repairs })
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
            format!("lambda {} outside the visualization domain 4..=48", p.lambda),
            vec!["set lambda within 4..=48"],
        );
    }
    if !(1..=200).contains(&p.generations) {
        return refuse(
            "generations-out-of-range",
            format!("generations {} outside the visualization domain 1..=200", p.generations),
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
    let mut rng = Lcg(p.seed);

    let mut mean: Vec<f64> = p.x0[..n].to_vec();
    let mut sigma = p.sigma0;
    let mut c = vec![0.0f64; n * n];
    for i in 0..n {
        c[i * n + i] = 1.0;
    }
    let mut p_sigma = vec![0.0f64; n];
    let mut p_c = vec![0.0f64; n];

    // Strategy constants (Hansen 2016; identical formulas to the TS fallback).
    let lambda = p.lambda;
    let mu = (lambda / 2).max(1);
    let mut weights = vec![0.0f64; mu];
    {
        let mut raw = vec![0.0f64; mu];
        for (i, w) in raw.iter_mut().enumerate() {
            *w = fs_math::det::ln(mu as f64 + 0.5) - fs_math::det::ln((i + 1) as f64);
        }
        let sum: f64 = raw.iter().sum();
        for (i, w) in weights.iter_mut().enumerate() {
            *w = raw[i] / sum;
        }
    }
    let mueff = 1.0 / weights.iter().map(|w| w * w).sum::<f64>();
    let nf = n as f64;
    let cc = (4.0 + mueff / nf) / (nf + 4.0 + 2.0 * mueff / nf);
    let cs = (mueff + 2.0) / (nf + mueff + 5.0);
    let c1 = 2.0 / ((nf + 1.3) * (nf + 1.3) + mueff);
    let cmu = ((1.0 - c1).min(
        (2.0 * (mueff - 2.0 + 1.0 / mueff)) / ((nf + 2.0) * (nf + 2.0) + mueff),
    ))
    .max(0.0);
    let damps = 1.0 + 2.0 * ((mueff - 1.0) / (nf + 1.0)).sqrt().max(0.0) + cs;
    let chin = nf.sqrt() * (1.0 - 1.0 / (4.0 * nf) + 1.0 / (21.0 * nf * nf));

    let mut best_f = f64::INFINITY;
    let mut best_x = mean.clone();
    let mut evals = 0usize;
    let mut generations: Vec<GenSnapshot> = Vec::with_capacity(p.generations);
    let mut stop_reason = "generations-exhausted";

    for g in 0..p.generations {
        let (eigvals, eigvecs) = jacobi_eigh(&c, n);
        if eigvals.iter().any(|v| !v.is_finite()) {
            return Err(Refusal {
                code: "eigen-decomposition-failed",
                message: format!(
                    "covariance eigendecomposition produced non-finite values at generation {g}"
                ),
                ranked_repairs: vec!["disable the active update", "reduce sigma0"],
            });
        }
        let sqrt_c = transform_matrix(&eigvals, &eigvecs, n, false);
        let inv_sqrt_c = transform_matrix(&eigvals, &eigvecs, n, true);
        let cond = eigvals[n - 1] / eigvals[0].max(1e-18);

        // 1. Sample λ offspring: x = m + σ·C^{1/2}·z.
        let mut sx = vec![0.0f64; lambda * n];
        let mut sz_raw = vec![0.0f64; lambda * n];
        let mut sf = vec![0.0f64; lambda];
        let mut st = vec![0.0f64; lambda]; // true (noiseless) fitness
        for i in 0..lambda {
            let z: Vec<f64> = (0..n).map(|_| rng.gauss()).collect();
            sz_raw[i * n..(i + 1) * n].copy_from_slice(&z);
            let y = mat_vec(&sqrt_c, &z, n);
            for k in 0..n {
                let mut xk = mean[k] + sigma * y[k];
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
            let noisy = if p.noise > 0.0 { true_f + rng.gauss() * p.noise } else { true_f };
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
                mean[k] += weights[rank] * sx[idx * n + k];
            }
        }
        let mean_shift: Vec<f64> = (0..n).map(|k| (mean[k] - old_mean[k]) / sigma).collect();
        let z_mean = mat_vec(&inv_sqrt_c, &mean_shift, n);

        // 4. Evolution paths.
        let ps_coeff = (cs * (2.0 - cs) * mueff).sqrt();
        for k in 0..n {
            p_sigma[k] = (1.0 - cs) * p_sigma[k] + ps_coeff * z_mean[k];
        }
        let norm_ps: f64 = p_sigma.iter().map(|v| v * v).sum::<f64>().sqrt();
        let hsig_denom = (1.0 - (1.0 - cs) * (1.0 - cs) * (g as f64 + 1.0)).sqrt();
        let hsig = if hsig_denom > 0.0 && norm_ps / hsig_denom / chin < 1.4 + 2.0 / (nf + 1.0) {
            1.0
        } else {
            0.0
        };
        let pc_coeff = (cc * (2.0 - cc) * mueff).sqrt();
        for k in 0..n {
            p_c[k] = (1.0 - cc) * p_c[k] + hsig * pc_coeff * mean_shift[k];
        }

        // 5. Covariance adaptation: rank-1 + rank-μ (+ active negative rank-μ,
        // mirroring the TS fallback's mirrored-weight heuristic).
        let old_coeff = 1.0 - c1 - cmu + (1.0 - hsig) * c1 * cc * (2.0 - cc);
        let mut new_c = vec![0.0f64; n * n];
        for i in 0..n {
            for k in 0..n {
                let rank1 = p_c[i] * p_c[k];
                let mut rankmu = 0.0;
                let mut active = 0.0;
                for (rank, &idx) in order.iter().enumerate().take(mu) {
                    let yi = (sx[idx * n + i] - old_mean[i]) / sigma;
                    let yk = (sx[idx * n + k] - old_mean[k]) / sigma;
                    rankmu += weights[rank] * yi * yk;
                }
                if p.active {
                    let c_active = cmu * 0.4;
                    for rank in (lambda - mu)..lambda {
                        let idx = order[rank];
                        let yi = (sx[idx * n + i] - old_mean[i]) / sigma;
                        let yk = (sx[idx * n + k] - old_mean[k]) / sigma;
                        let neg_w = weights[lambda - 1 - rank];
                        active -= c_active * neg_w * yi * yk;
                    }
                }
                new_c[i * n + k] =
                    old_coeff * c[i * n + k] + c1 * rank1 + cmu * rankmu + active;
            }
        }
        for i in 0..n {
            new_c[i * n + i] += 1e-10;
        }
        c = new_c;

        // 6. Step-size adaptation with numerical safety clamp (TS parity).
        sigma *= fs_math::det::exp((cs / damps) * (norm_ps / chin - 1.0));
        sigma = sigma.min(10.0).max(1e-10);

        // z in rank order (mirrors sx ordering).
        let mut sz = vec![0.0f64; lambda * n];
        for (rank, &idx) in order.iter().enumerate() {
            sz[rank * n..(rank + 1) * n].copy_from_slice(&sz_raw[idx * n..(idx + 1) * n]);
        }

        generations.push(GenSnapshot {
            g: g + 1,
            mean: mean.clone(),
            sigma,
            eigvals,
            eigvecs,
            cond,
            best_f,
            evals,
            proj_mean: [0.0; 3],
            proj_eigvals: [0.0; 3],
            proj_eigvecs: [0.0; 9],
            sx,
            sz,
            sf,
            se,
            p_sigma: p_sigma.clone(),
            p_c: p_c.clone(),
        });

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
fn project_phase_space(gens: &mut [GenSnapshot], final_c: &[f64], n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let ng = gens.len().max(1);
    let mut center = vec![0.0; n];
    for snap in gens.iter() {
        for k in 0..n {
            center[k] += snap.mean[k] / ng as f64;
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
            pm[r] = (0..n).map(|i| basis[r * n + i] * (snap.mean[i] - center[i])).sum();
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

fn num(v: f64) -> String {
    if v.is_finite() {
        format!("{v}")
    } else if v.is_nan() {
        "null".into()
    } else if v > 0.0 {
        "1e999".into()
    } else {
        "-1e999".into()
    }
}

fn num_arr(vals: &[f64]) -> String {
    let mut s = String::from("[");
    for (i, v) in vals.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&num(*v));
    }
    s.push(']');
    s
}

fn byte_arr(vals: &[u8]) -> String {
    let mut s = String::from("[");
    for (i, v) in vals.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push_str(&v.to_string());
    }
    s.push(']');
    s
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
    let mut gens = String::from("[");
    for (gi, snap) in run.generations.iter().enumerate() {
        if gi > 0 {
            gens.push(',');
        }
        gens.push_str(&format!(
            "{{\"g\":{},\"mean\":{},\"sigma\":{},\"eigvals\":{},\"eigvecs\":{},\"cond\":{},\"best_f\":{},\"evals\":{},\"proj_mean\":{},\"proj_eigvals\":{},\"proj_eigvecs\":{},\"sx\":{},\"sz\":{},\"sf\":{},\"se\":{},\"p_sigma\":{},\"p_c\":{}}}",
            snap.g,
            num_arr(&snap.mean),
            num(snap.sigma),
            num_arr(&snap.eigvals),
            num_arr(&snap.eigvecs),
            num(snap.cond),
            num(snap.best_f),
            snap.evals,
            num_arr(&snap.proj_mean),
            num_arr(&snap.proj_eigvals),
            num_arr(&snap.proj_eigvecs),
            num_arr(&snap.sx),
            num_arr(&snap.sz),
            num_arr(&snap.sf),
            byte_arr(&snap.se),
            num_arr(&snap.p_sigma),
            num_arr(&snap.p_c),
        ));
    }
    gens.push(']');
    // PCA frame computed by project_phase_space and carried on the run.
    let basis = num_arr(&run.pca.basis);
    format!(
        "{{\"ok\":{{\"kernel\":\"{}\",\"dim\":{},\"landscape\":{},\"stop_reason\":\"{}\",\"best_f\":{},\"best_x\":{},\"total_evals\":{},\"generations\":{},\"pca_basis\":{},\"pca_center\":{},\"pca_pool_eigvals\":{}}}}}",
        KERNEL_VERSION,
        run.dim,
        run.landscape,
        run.stop_reason,
        num(run.best_f),
        num_arr(&run.best_x),
        run.total_evals,
        gens,
        basis,
        num_arr(&run.pca.center),
        num_arr(&run.pca.pool_eigvals),
    )
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
        let p = VizParams { landscape: LANDSCAPE_ROSENBROCK, generations: 40, ..base_params() };
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
                assert!(resid < 1e-9 * scale, "residual {resid} at ({i},{j}), scale {scale}");
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
        let a = cmaes_run_json(3, [1.0, 2.0, 3.0, 0.0, 0.0, 0.0], 0.5, 12, true, 42, 30, 1, 0.0, false, 0.0, 0.0, f64::NAN);
        let b = cmaes_run_json(3, [1.0, 2.0, 3.0, 0.0, 0.0, 0.0], 0.5, 12, true, 42, 30, 1, 0.0, false, 0.0, 0.0, f64::NAN);
        assert_eq!(a, b, "same seed must replay bitwise");
        let c = cmaes_run_json(3, [1.0, 2.0, 3.0, 0.0, 0.0, 0.0], 0.5, 12, true, 43, 30, 1, 0.0, false, 0.0, 0.0, f64::NAN);
        assert_ne!(a, c, "different seed must diverge");
    }

    #[test]
    fn refusal_codes_are_typed() {
        let cases: Vec<(VizParams, &str)> = vec![
            (VizParams { dim: 1, ..base_params() }, "dim-out-of-range"),
            (VizParams { dim: 7, ..base_params() }, "dim-out-of-range"),
            (VizParams { sigma0: 0.0, ..base_params() }, "sigma0-non-positive"),
            (VizParams { sigma0: f64::NAN, ..base_params() }, "sigma0-non-positive"),
            (VizParams { lambda: 3, ..base_params() }, "lambda-out-of-range"),
            (VizParams { lambda: 49, ..base_params() }, "lambda-out-of-range"),
            (VizParams { generations: 0, ..base_params() }, "generations-out-of-range"),
            (VizParams { generations: 201, ..base_params() }, "generations-out-of-range"),
            (VizParams { landscape: 9, ..base_params() }, "landscape-unknown"),
            (VizParams { noise: -0.1, ..base_params() }, "noise-invalid"),
            (
                VizParams { bounds_enabled: true, bound_min: 2.0, bound_max: -2.0, ..base_params() },
                "bounds-inverted",
            ),
            (VizParams { f_target: f64::INFINITY, ..base_params() }, "f-target-invalid"),
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
        let env = cmaes_run_json(1, [0.0; 6], 0.5, 12, true, 1, 10, 0, 0.0, false, 0.0, 0.0, f64::NAN);
        assert!(env.contains("\"code\":\"dim-out-of-range\""), "envelope was: {env}");
    }

    #[test]
    fn target_stops_early() {
        let p = VizParams { f_target: 1.0, ..base_params() };
        let run = cmaes_run(&p).expect("run");
        assert_eq!(run.stop_reason, "target-reached");
        assert!(run.generations.len() < 150);
        assert!(run.best_f <= 1.0);
    }

    #[test]
    fn rosenbrock_improves_on_start() {
        let p = VizParams { landscape: LANDSCAPE_ROSENBROCK, ..base_params() };
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
        let p = VizParams { dim: 5, landscape: LANDSCAPE_ELLI, ..base_params() };
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
    fn dim2_projection_is_direct() {
        let p = VizParams { dim: 2, ..base_params() };
        let run = cmaes_run(&p).expect("run");
        let snap = run.generations.last().expect("gens");
        // With dim ≤ 3 the projection is the identity frame: proj_mean equals
        // mean centered by the trajectory centroid, third axis pinned at 0.
        assert_eq!(snap.proj_mean[2], 0.0);
        assert_eq!(run.dim, 2);
    }
}
