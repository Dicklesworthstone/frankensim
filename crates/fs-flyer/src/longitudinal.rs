//! Longitudinal small-perturbation linearization + the V-02a GATE (bead
//! wf-root-guzez.5.13.3, E4.6a-iii). Central finite differences of the
//! FULL force build-up about the fixed-throttle trim give the classical
//! 4-state body-axes system x = (u, w, q, θ):
//!
//!   u̇ = X_u u + X_w w + (X_q − w₀) q − g cosθ₀ θ
//!   ẇ = Z_u u + Z_w w + (Z_q + u₀) q − g sinθ₀ θ
//!   q̇ = (M_u u + M_w w + M_q q) / I_yy
//!   θ̇ = q
//!
//! Eigenvalues from the characteristic quartic (Leverrier–Faddeev
//! coefficients, fixed-iteration Durand–Kerner roots — deterministic).
//!
//! The GATE holds the model to the A4 (Culick/Jex) claims at the level
//! the dossier PERMITS (culick-a4-anchors-v1.json): pole STRUCTURE
//! (an unstable longitudinal mode exists), pitch time-to-double inside
//! the declared order band (the value is REPORTED, never a point claim),
//! and derivative SIGNS (M_α > 0 statically unstable, M_q < 0 damping).
//! Quantitative pole levels stay out of reach until the derivative
//! tables are transcribed (machine_readable_absence).

use crate::Refusal;
use crate::aircraft::{OpenLoopDesign, TrimResult};
use fs_math::det;

/// Pitch inertia, Culick/Jex mass-reconstruction class
/// (culick-a4-anchors-v1: ~289 slug·ft², Estimated ±25%).
pub const IYY_KG_M2: f64 = 392.0;

/// The A4 time-to-double order band [s] (anchors file).
pub const T2_BAND_S: (f64, f64) = (0.15, 3.0);

/// One complex eigenvalue.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pole {
    /// Real part [1/s].
    pub re: f64,
    /// Imaginary part [rad/s].
    pub im: f64,
}

/// Linearization + gate receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct LongitudinalReport {
    /// The 4×4 A matrix, rows (u̇, ẇ, q̇, θ̇).
    pub a: [[f64; 4]; 4],
    /// Dimensional derivatives (X_u, X_w, Z_u, Z_w, M_u, M_w, M_q).
    pub derivatives: [f64; 7],
    /// dM/dα at trim [N·m/rad] (the static-stability sign carrier).
    pub m_alpha_nm_per_rad: f64,
    /// dM/dq at trim [N·m·s/rad] (the damping sign carrier).
    pub m_q_nm_s_per_rad: f64,
    /// The four poles.
    pub poles: [Pole; 4],
    /// Largest real part [1/s].
    pub max_re: f64,
    /// Time to double of the dominant unstable mode [s] (ln2 / max_re).
    pub time_to_double_s: f64,
    /// The trim this was linearized about.
    pub trim: TrimResult,
}

/// V-02a gate verdict (per-clause, never a bare boolean).
#[derive(Clone, Debug, PartialEq)]
pub struct V02aVerdict {
    /// An unstable longitudinal mode exists.
    pub unstable_mode_present: bool,
    /// Time-to-double inside the A4 order band.
    pub t2_in_band: bool,
    /// M_α > 0 (statically unstable, matching A4).
    pub m_alpha_sign_ok: bool,
    /// M_q < 0 (pitch damping present).
    pub m_q_sign_ok: bool,
    /// Scalar gate score (lower = closer to the A4 claims; the
    /// anti-vacuity comparison quantity).
    pub score: f64,
}

impl V02aVerdict {
    /// All clauses green.
    #[must_use]
    pub fn passes(&self) -> bool {
        self.unstable_mode_present && self.t2_in_band && self.m_alpha_sign_ok && self.m_q_sign_ok
    }
}

/// Linearize the design about a trim (central differences; fixed,
/// deterministic steps).
///
/// # Errors
/// Build-up refusals pass through.
pub fn linearize(
    design: &OpenLoopDesign,
    trim: &TrimResult,
    rho_kg_m3: f64,
) -> Result<LongitudinalReport, Refusal> {
    let m = design.gross_mass_kg;
    let (v0, a0, dc, om) = (
        trim.v_mps,
        trim.alpha_rad,
        trim.delta_canard_rad,
        trim.omega_prop_rad_s,
    );
    let (u0, w0) = (v0 * det::cos(a0), v0 * det::sin(a0));
    let g = 9.80665;
    // Perturb (u, w) via the (V, alpha) map, q directly.
    let du = 0.05f64;
    let dw = 0.05f64;
    let dq = 0.01f64;
    let eval = |u: f64, w: f64, q: f64| -> Result<[f64; 3], Refusal> {
        let uu = u * u;
        let ww = w * w;
        let v = det::sqrt(uu + ww);
        let alpha = det::atan2(w, u);
        let b = design.force_buildup(v, alpha, dc, om, q, rho_kg_m3)?;
        Ok([b.force_n[0], b.force_n[2], b.moment_y_nm])
    };
    let fd = |p: [f64; 3], mmm: [f64; 3], h: f64| -> [f64; 3] {
        [
            (p[0] - mmm[0]) / (2.0 * h),
            (p[1] - mmm[1]) / (2.0 * h),
            (p[2] - mmm[2]) / (2.0 * h),
        ]
    };
    let d_u = fd(eval(u0 + du, w0, 0.0)?, eval(u0 - du, w0, 0.0)?, du);
    let d_w = fd(eval(u0, w0 + dw, 0.0)?, eval(u0, w0 - dw, 0.0)?, dw);
    let d_q = fd(eval(u0, w0, dq)?, eval(u0, w0, -dq)?, dq);
    let (x_u, z_u, m_u) = (d_u[0] / m, d_u[1] / m, d_u[2]);
    let (x_w, z_w, m_w) = (d_w[0] / m, d_w[1] / m, d_w[2]);
    let (x_q, z_q, m_q) = (d_q[0] / m, d_q[1] / m, d_q[2]);
    let a = [
        [x_u, x_w, x_q - w0, -g * det::cos(a0)],
        [z_u, z_w, z_q + u0, -g * det::sin(a0)],
        [m_u / IYY_KG_M2, m_w / IYY_KG_M2, m_q / IYY_KG_M2, 0.0],
        [0.0, 0.0, 1.0, 0.0],
    ];
    let poles = eig4(&a);
    let max_re = poles.iter().map(|p| p.re).fold(f64::NEG_INFINITY, f64::max);
    let time_to_double_s = if max_re > 0.0 {
        core::f64::consts::LN_2 / max_re
    } else {
        f64::INFINITY
    };
    // dM/dalpha at constant V (chain rule: dα at fixed V is dw = V·dα
    // with du = −w0·dα to first order; use the direct FD instead).
    let da = 0.002;
    let mp = design
        .force_buildup(v0, a0 + da, dc, om, 0.0, rho_kg_m3)?
        .moment_y_nm;
    let mm = design
        .force_buildup(v0, a0 - da, dc, om, 0.0, rho_kg_m3)?
        .moment_y_nm;
    let m_alpha = (mp - mm) / (2.0 * da);
    Ok(LongitudinalReport {
        a,
        derivatives: [x_u, x_w, z_u, z_w, m_u, m_w, m_q],
        m_alpha_nm_per_rad: m_alpha,
        m_q_nm_s_per_rad: m_q,
        poles,
        max_re,
        time_to_double_s,
        trim: trim.clone(),
    })
}

/// Apply the V-02a gate clauses (anchors file levels).
#[must_use]
pub fn v02a_gate(rep: &LongitudinalReport) -> V02aVerdict {
    let unstable = rep.max_re > 0.0;
    let t2 = rep.time_to_double_s;
    let t2_in_band = t2 >= T2_BAND_S.0 && t2 <= T2_BAND_S.1;
    let ma_ok = rep.m_alpha_nm_per_rad > 0.0;
    let mq_ok = rep.m_q_nm_s_per_rad < 0.0;
    // Score: log-distance of t2 from the band's geometric center plus
    // heavy penalties for structural misses (used by the anti-vacuity
    // comparison; lower is better).
    let center = det::sqrt(T2_BAND_S.0 * T2_BAND_S.1);
    let mut score = if t2.is_finite() {
        det::ln(t2 / center).abs()
    } else {
        10.0
    };
    if !unstable {
        score += 10.0;
    }
    if !ma_ok {
        score += 5.0;
    }
    if !mq_ok {
        score += 5.0;
    }
    V02aVerdict {
        unstable_mode_present: unstable,
        t2_in_band,
        m_alpha_sign_ok: ma_ok,
        m_q_sign_ok: mq_ok,
        score,
    }
}

/// Characteristic-polynomial coefficients via Leverrier–Faddeev, then
/// fixed-iteration Durand–Kerner (deterministic init and count).
/// Public so batteries can verify it against analytic fixtures.
#[must_use]
pub fn eig4(a: &[[f64; 4]; 4]) -> [Pole; 4] {
    // p(λ) = λ⁴ + c3 λ³ + c2 λ² + c1 λ + c0.
    let mut b = *a;
    let mut c = [0.0f64; 4]; // c[3] = c3 ... c[0] = c0
    let tr = |m: &[[f64; 4]; 4]| m[0][0] + m[1][1] + m[2][2] + m[3][3];
    let matmul = |x: &[[f64; 4]; 4], y: &[[f64; 4]; 4]| {
        let mut out = [[0.0f64; 4]; 4];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate() {
                *v = (0..4).map(|k| x[i][k] * y[k][j]).sum();
            }
        }
        out
    };
    let mut coeff = -tr(&b);
    c[3] = coeff;
    for step in 2..=4 {
        for (i, row) in b.iter_mut().enumerate() {
            row[i] += coeff;
        }
        b = matmul(a, &b);
        coeff = -tr(&b) / step as f64;
        c[4 - step] = coeff;
    }
    // Durand–Kerner on p with deterministic starts (0.4+0.9i)^k.
    let (mut zr, mut zi) = ([0.0f64; 4], [0.0f64; 4]);
    let (mut pr, mut pi) = (1.0f64, 0.0f64);
    for k in 0..4 {
        // Statement-split (guzez.7.2.1): complex mul without FMA.
        let nr_a = pr * 0.4;
        let nr_b = pi * 0.9;
        let nr = nr_a - nr_b;
        let ni_a = pr * 0.9;
        let ni_b = pi * 0.4;
        let ni = ni_a + ni_b;
        pr = nr;
        pi = ni;
        zr[k] = pr;
        zi[k] = pi;
    }
    let poly = |x: f64, y: f64| -> (f64, f64) {
        // Horner for p(z), z = x+iy, coeffs [1, c3, c2, c1, c0].
        let (mut ar, mut ai) = (1.0f64, 0.0f64);
        for &cc in &[c[3], c[2], c[1], c[0]] {
            let nr_a = ar * x;
            let nr_b = ai * y;
            let nr = (nr_a - nr_b) + cc;
            let ni_a = ar * y;
            let ni_b = ai * x;
            let ni = ni_a + ni_b;
            ar = nr;
            ai = ni;
        }
        (ar, ai)
    };
    for _ in 0..200 {
        for k in 0..4 {
            let (pvr, pvi) = poly(zr[k], zi[k]);
            // denominator: prod_{j!=k} (z_k - z_j)
            let (mut dr, mut di) = (1.0f64, 0.0f64);
            for j in 0..4 {
                if j != k {
                    let (ur, ui) = (zr[k] - zr[j], zi[k] - zi[j]);
                    let nr_a = dr * ur;
                    let nr_b = di * ui;
                    let nr = nr_a - nr_b;
                    let ni_a = dr * ui;
                    let ni_b = di * ur;
                    let ni = ni_a + ni_b;
                    dr = nr;
                    di = ni;
                }
            }
            let den_a = dr * dr;
            let den_b = di * di;
            let den = den_a + den_b;
            if den > 1e-300 {
                let qr_a = pvr * dr;
                let qr_b = pvi * di;
                let qr = (qr_a + qr_b) / den;
                let qi_a = pvi * dr;
                let qi_b = pvr * di;
                let qi = (qi_a - qi_b) / den;
                zr[k] -= qr;
                zi[k] -= qi;
            }
        }
    }
    // Canonical order: by (re, im) descending re then ascending im — a
    // FIXED tie rule so the pole list is deterministic.
    let mut idx = [0usize, 1, 2, 3];
    idx.sort_by(|&i, &j| {
        zr[j]
            .partial_cmp(&zr[i])
            .unwrap_or(core::cmp::Ordering::Equal)
            .then(
                zi[i]
                    .partial_cmp(&zi[j])
                    .unwrap_or(core::cmp::Ordering::Equal),
            )
    });
    let mut out = [Pole { re: 0.0, im: 0.0 }; 4];
    for (n, &i) in idx.iter().enumerate() {
        out[n] = Pole {
            re: zr[i],
            im: zi[i],
        };
    }
    out
}
