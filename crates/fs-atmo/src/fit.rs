//! Complete Hermitian fit of the mode basis to the Mann target (bead
//! wf-root-guzez.4.6.2.2, E3.3b-ii-b). Round-3 Q2: diagonal-only may
//! INITIALIZE, never ACCEPT — and in this basis that law is a theorem:
//! with independent amplitude components the wall-compatible curl gives
//! E[u·w] ≡ 0 identically (phase averages kill every diagonal term),
//! so the surface shear stress uw < 0 is reachable ONLY through
//! amplitude cross-covariances.
//!
//! The FINITE basis changes the realizable tensor (the plan's
//! per-truncation clause): at the reference height this basis cannot in
//! general reach the full Mann uw at exact diagonals — measured, not
//! assumed. The fit is therefore CONSTRAINED: per-mode amplitude
//! covariances stay PSD (correlation budget ρ ≤ ρ_max), the uw channel
//! is capacity-limited, and the achieved errors are RECORDED as the
//! per-truncation artifact. Deterministic coordinate descent (fixed
//! sweep count, fixed step schedule) — bitwise repeatable.
//!
//! Per-mode ensemble closed forms (θ uniform, amplitudes zero-mean):
//!   E[u²]  = ch²·(k_y²·s_z² + k_h²·s_y²)/2
//!   E[v²]  = ch²·(k_h²·s_x² + k_x²·s_z²)/2
//!   E[w²]  = sh²·(k_y²·s_x² + k_x²·s_y²)/2
//! and per-mode uw capacity ρ·√s_z²·√(α²·s_x² + β²·s_y²) with
//! α = k_y²|ch·sh|/2, β = |k_x·k_y·ch·sh|/2 (sign-aligned crosses).

use crate::mann::{MannParams, energy_spectrum, stress_integrals};
use crate::{Mode, Refusal, StreamKey, TurbulenceField};
use fs_math::det;

/// Mode-count cap for the fit (matches the field cap).
pub const MAX_FIT_MODES: usize = crate::MAX_MODES;

/// Correlation budget cap (strict PSD margin).
pub const RHO_MAX: f64 = 0.95;

/// Coordinate-descent sweeps (fixed — deterministic).
pub const FIT_SWEEPS: usize = 120;

/// The fitted, PSD-constrained amplitude model.
#[derive(Clone, Debug, PartialEq)]
pub struct FittedAmplitudes {
    /// E[a_x²] scale (per unit spectral weight).
    pub sx2: f64,
    /// E[a_y²] scale.
    pub sy2: f64,
    /// E[a_z²] scale.
    pub sz2: f64,
    /// Cross-correlation budget used, ρ ∈ (0, RHO_MAX].
    pub rho: f64,
    /// Fraction of the target |uw| the truncation realizes ∈ (0, 1].
    pub uw_fraction: f64,
    /// Scale-normalized errors on (u2, v2, w2, uw) — the V-04b1 fit
    /// receipt AND the per-truncation artifact.
    pub errors: [f64; 4],
    /// Mode count this fit is FOR.
    pub n_modes: usize,
    /// Reference height [m].
    pub h_ref_m: f64,
}

/// Deterministic mode-geometry list for a seed (the same draws the
/// field builder makes — geometry only).
fn mode_geometry(seed: u64, n_modes: usize, length_m: f64) -> Vec<(f64, f64, f64)> {
    let base_k = core::f64::consts::TAU / (8.0 * length_m);
    let mut out = Vec::with_capacity(n_modes);
    for i in 0..n_modes {
        let mut s = StreamKey {
            seed,
            kernel: crate::ATMO_KERNEL,
            tile: i as u32,
        }
        .stream();
        let shell = base_k * det::exp(det::ln(64.0) * s.next_f64());
        let angle = core::f64::consts::TAU * s.next_f64();
        let kx = shell * det::cos(angle);
        let ky = shell * det::sin(angle);
        let kh = shell * (0.5 + s.next_f64());
        out.push((kx, ky, kh));
    }
    out
}

/// Per-mode spectral weight (von Kármán level at the mode's shell).
fn spectral_weight(p: &MannParams, kx: f64, ky: f64, kh: f64) -> f64 {
    let k = det::sqrt(kx * kx + ky * ky + kh * kh);
    energy_spectrum(p, k) / (k * k)
}

struct Precomp {
    /// Diagonal map rows: (u2|v2|w2) coefficients on (sx2, sy2, sz2).
    diag: [[f64; 3]; 3],
    /// Per-mode (w·α, w·β) for the uw capacity.
    uw_terms: Vec<(f64, f64)>,
}

fn precompute(geo: &[(f64, f64, f64)], weights: &[f64], h_ref: f64) -> Precomp {
    let mut diag = [[0.0f64; 3]; 3];
    let mut uw_terms = Vec::with_capacity(geo.len());
    for ((kx, ky, kh), w) in geo.iter().zip(weights) {
        let ch = det::cos(kh * h_ref);
        let sh = det::sin(kh * h_ref);
        let (c2, s2) = (ch * ch, sh * sh);
        diag[0][1] += w * c2 * kh * kh / 2.0;
        diag[0][2] += w * c2 * ky * ky / 2.0;
        diag[1][0] += w * c2 * kh * kh / 2.0;
        diag[1][2] += w * c2 * kx * kx / 2.0;
        diag[2][0] += w * s2 * ky * ky / 2.0;
        diag[2][1] += w * s2 * kx * kx / 2.0;
        let alpha = ky * ky * (ch * sh).abs() / 2.0;
        let beta = (kx * ky * ch * sh).abs() / 2.0;
        uw_terms.push((w * alpha, w * beta));
    }
    Precomp { diag, uw_terms }
}

/// Realized (u2, v2, w2, |uw|_capacity) for a parameter state.
fn realize(pc: &Precomp, s: &[f64; 3], rho: f64) -> [f64; 4] {
    let mut out = [0.0f64; 4];
    for r in 0..3 {
        out[r] = (0..3).map(|c| pc.diag[r][c] * s[c]).sum();
    }
    let cap: f64 = pc
        .uw_terms
        .iter()
        .map(|&(wa, wb)| det::sqrt(wa * wa * s[0] + wb * wb * s[1]))
        .sum::<f64>()
        * rho
        * det::sqrt(s[2]);
    out[3] = cap;
    out
}

/// Fit the constrained amplitude model for (seed, n_modes).
///
/// # Errors
/// `fit-params-invalid` (mode cap at cap AND cap+1, h_ref ≤ 0);
/// `fit-degenerate` (no positive-variance descent point).
pub fn fit_amplitudes(
    seed: u64,
    n_modes: usize,
    p: &MannParams,
    h_ref_m: f64,
) -> Result<FittedAmplitudes, Refusal> {
    p.admit()?;
    if n_modes == 0 || n_modes > MAX_FIT_MODES || !(h_ref_m > 0.0) {
        return Err(Refusal {
            code: "fit-params-invalid",
            message: format!("n_modes {n_modes}, h_ref {h_ref_m:?}"),
            ranked_repairs: vec![format!("1..={MAX_FIT_MODES} modes; h_ref > 0")],
        });
    }
    let t = stress_integrals(p)?;
    let target = [t[0][0], t[1][1], t[2][2], -t[0][2]]; // |uw| in slot 3
    let geo = mode_geometry(seed, n_modes, p.length_m);
    let weights: Vec<f64> = geo
        .iter()
        .map(|&(kx, ky, kh)| spectral_weight(p, kx, ky, kh))
        .collect();
    let pc = precompute(&geo, &weights, h_ref_m);
    // Objective: scale-normalized squared error; realized uw =
    // min(capacity, target) so excess capacity is never penalized.
    let cost = |s: &[f64; 3], rho: f64| -> f64 {
        let r = realize(&pc, s, rho);
        let mut j = 0.0;
        for c in 0..3 {
            let e = (r[c] - target[c]) / target[c];
            j += e * e;
        }
        let uw_real = r[3].min(target[3]);
        let e = (uw_real - target[3]) / target[3];
        j + e * e
    };
    // Deterministic start: equal variances scaled to the trace.
    let trace_coeff: f64 = (0..3)
        .map(|r| (0..3).map(|c| pc.diag[r][c]).sum::<f64>())
        .sum();
    let s0 = (target[0] + target[1] + target[2]) / trace_coeff.max(1e-300);
    if !(s0 > 0.0) {
        return Err(Refusal {
            code: "fit-degenerate",
            message: "zero diagonal capacity".into(),
            ranked_repairs: vec!["more modes".into()],
        });
    }
    let mut s = [s0, s0, s0];
    let mut rho = RHO_MAX;
    let mut j = cost(&s, rho);
    // Fixed multiplicative coordinate descent.
    let mut step = 0.30f64;
    for _sweep in 0..FIT_SWEEPS {
        for c in 0..3 {
            for dir in [1.0 + step, 1.0 / (1.0 + step)] {
                let mut cand = s;
                cand[c] *= dir;
                let jc = cost(&cand, rho);
                if jc < j {
                    s = cand;
                    j = jc;
                }
            }
        }
        for dir in [1.0 + step, 1.0 / (1.0 + step)] {
            let cand = (rho * dir).min(RHO_MAX);
            let jc = cost(&s, cand);
            if jc < j {
                rho = cand;
                j = jc;
            }
        }
        step *= 0.96;
    }
    let r = realize(&pc, &s, rho);
    let uw_real = r[3].min(target[3]);
    let errors = [
        (r[0] - target[0]) / target[0],
        (r[1] - target[1]) / target[1],
        (r[2] - target[2]) / target[2],
        (uw_real - target[3]) / target[3],
    ];
    Ok(FittedAmplitudes {
        sx2: s[0],
        sy2: s[1],
        sz2: s[2],
        rho,
        uw_fraction: uw_real / target[3],
        errors,
        n_modes,
        h_ref_m,
    })
}

/// Build a turbulence field whose amplitudes are drawn from the FITTED
/// covariance: per-mode sign-aligned crosses at the fitted budget,
/// uniformly scaled so the aggregate uw equals the REALIZED value.
///
/// # Errors
/// Parameter refusals pass through.
pub fn build_fitted_field(
    seed: u64,
    fit: &FittedAmplitudes,
    p: &MannParams,
    u_adv_mps: f64,
) -> Result<TurbulenceField, Refusal> {
    let geo = mode_geometry(seed, fit.n_modes, p.length_m);
    let weights: Vec<f64> = geo
        .iter()
        .map(|&(kx, ky, kh)| spectral_weight(p, kx, ky, kh))
        .collect();
    // Recompute the capacity scale so the builder carries no hidden
    // state: aggregate uw must land on uw_fraction × target.
    let pc = precompute(&geo, &weights, fit.h_ref_m);
    let s = [fit.sx2, fit.sy2, fit.sz2];
    let cap = realize(&pc, &s, fit.rho)[3];
    let t = stress_integrals(p)?;
    let uw_wanted = (-t[0][2]) * fit.uw_fraction;
    let scale = if cap > 0.0 {
        (uw_wanted / cap).min(1.0)
    } else {
        0.0
    };
    let mut modes = Vec::with_capacity(fit.n_modes);
    for (i, &(kx, ky, kh)) in geo.iter().enumerate() {
        let mut st = StreamKey {
            seed,
            kernel: crate::ATMO_KERNEL,
            tile: i as u32,
        }
        .stream();
        // Re-consume the geometry draws so the amplitude draws land on
        // the same counter offsets regardless of caller.
        let _ = st.next_f64();
        let _ = st.next_f64();
        let _ = st.next_f64();
        let w = weights[i];
        let ch = det::cos(kh * fit.h_ref_m);
        let sh = det::sin(kh * fit.h_ref_m);
        let alpha = ky * ky * (ch * sh).abs() / 2.0;
        let beta = (kx * ky * ch * sh).abs() / 2.0;
        // Per-mode optimal direction (γx ∝ sx2·α, γy ∝ sy2·β), budget
        // ρ·scale of the PSD ellipse.
        let norm = det::sqrt(alpha * alpha * fit.sx2 + beta * beta * fit.sy2);
        let (gx, gy) = if norm > 0.0 {
            let g = fit.rho * scale * det::sqrt(fit.sz2) / norm;
            (g * fit.sx2 * alpha, g * fit.sy2 * beta)
        } else {
            (0.0, 0.0)
        };
        let s1 = if ch * sh >= 0.0 { 1.0 } else { -1.0 };
        let s2 = if kx * ky * ch * sh >= 0.0 { 1.0 } else { -1.0 };
        let cxz_m = gx * w * s1;
        let cyz_m = -gy * w * s2;
        let l11 = det::sqrt(fit.sx2 * w);
        let l22 = det::sqrt(fit.sy2 * w);
        let l31 = if l11 > 0.0 { cxz_m / l11 } else { 0.0 };
        let l32 = if l22 > 0.0 { cyz_m / l22 } else { 0.0 };
        let rem = fit.sz2 * w - l31 * l31 - l32 * l32;
        let l33 = det::sqrt(rem.max(0.0));
        let (n1, n2, n3) = (st.next_normal(), st.next_normal(), st.next_normal());
        let ax = l11 * n1;
        let ay = l22 * n2;
        let az = l31 * n1 + l32 * n2 + l33 * n3;
        let phi = core::f64::consts::TAU * st.next_f64();
        modes.push(Mode {
            kx,
            ky,
            kh,
            a: [ax, ay, az],
            phi,
        });
    }
    Ok(TurbulenceField { modes, u_adv_mps })
}

/// EXACT per-realization stress: the θ-averaged one-point tensor of a
/// built field's DRAWN amplitudes at height h (closed form per mode —
/// no spatial window, no grid error). The V-04b2 unit statistic.
#[must_use]
pub fn realization_stress(field: &TurbulenceField, h_m: f64) -> [f64; 4] {
    let mut out = [0.0f64; 4];
    for m in &field.modes {
        let ch = det::cos(m.kh * h_m);
        let sh = det::sin(m.kh * h_m);
        let [ax, ay, az] = m.a;
        // E_theta with drawn amplitudes (st^2, ct^2 -> 1/2; st*ct -> 0).
        out[0] += ch * ch * (az * az * m.ky * m.ky + ay * ay * m.kh * m.kh) / 2.0;
        out[1] += ch * ch * (ax * ax * m.kh * m.kh + az * az * m.kx * m.kx) / 2.0;
        let cw = ax * m.ky - ay * m.kx;
        out[2] += sh * sh * cw * cw / 2.0;
        // E[ux*uz] = -ky*ch*sh*E[az*(ax*ky - ay*kx)]/2 with drawn amps.
        out[3] += -m.ky * ch * sh * az * (ax * m.ky - ay * m.kx) / 2.0;
    }
    out
}
