//! Stage 4: TRIM AND CERTIFY. A 2-state nonlinear pitch model per
//! candidate,
//!
//! ```text
//! θ̈ = −k(c)·θ·(1 − θ²/θ_s²) − d(c)·θ̇
//! ```
//!
//! with the restoring stiffness from the BEM lift slope, damping from
//! thickness and flapping (the smoke reduced model), and the cubic STALL
//! SOFTENING that makes the region of attraction genuinely finite: past
//! the stall angle `θ_s` the restoring moment reverses, the trims at
//! `θ = ±θ_s` are saddles, and their stable manifolds bound the basin.
//! (Koopman/DMD reduced models with per-trim conformal e-bands are the
//! recorded successor.)
//!
//! Stability is not a vibe: fs-sos PROVES a region of attraction — an
//! SOS Lyapunov/S-procedure certificate, verified with interval
//! arithmetic and interval Cholesky, that `V̇ ≤ −ε‖x‖²` on the ellipse
//! `{xᵀPx ≤ c}` — and the reported ROA volume is that ellipse's area. A
//! candidate whose certificate does not verify gets volume 0, never a
//! pretended basin. The linearization is cross-checked by the classic
//! `AᵀP + PA ≺ 0` Lyapunov test. The screening surrogate carries a
//! DISTRIBUTION-FREE conformal band (fs-surrogate) — certify-or-escalate,
//! gated on coverage.

use crate::param::OrnithCandidate;
use crate::screen::lift_to_drag;
use fs_sos::{
    MPoly, RoaOptions, SdpSettings, certify_roa_with, lyapunov_certifies_stability, solve_lyapunov,
};
use fs_surrogate::{ConformalBand, conformal_band};

/// Stall angle of the cubic pitch-softening model (rad).
pub const STALL_ANGLE: f64 = 0.35;

/// The stability certificate row.
#[derive(Debug, Clone)]
pub struct CertifyReport {
    /// Pitch dynamics linearization at trim (companion form).
    pub a: [[f64; 2]; 2],
    /// The Lyapunov matrix of the certified quadratic `V = xᵀPx`.
    pub p: [[f64; 2]; 2],
    /// SOS region-of-attraction certificate verified (and the
    /// linearization's Lyapunov inequality holds).
    pub certified: bool,
    /// The proved sublevel `c` of `{xᵀPx ≤ c}` (0.0 when uncertified).
    pub level: f64,
    /// Area of the certified ellipse in the (θ, θ̇) plane (0.0 when
    /// uncertified — never pretended).
    pub roa_volume: f64,
    /// Number of SOS feasibility solves the certificate took.
    pub sos_solves: usize,
    /// Maneuver proxy: control authority over pitch stiffness.
    pub maneuver: f64,
}

/// The candidate's pitch stiffness and damping: stiffness from the lift
/// slope at trim, damping from section thickness and flapping (thicker =
/// more damped, the smoke law).
#[must_use]
pub fn pitch_coefficients(c: &OrnithCandidate) -> (f64, f64) {
    let foil = c.section(crate::screen::PANELS);
    let dcl = fs_bem::panel2d::dcl_dalpha_adjoint(&foil, c.alpha)
        .expect("ornith candidate must satisfy the bounded adjoint contract");
    let k = 0.4 * dcl; // restoring stiffness ∝ lift slope
    let d = 8.0 * c.thickness + 0.4 * c.flap_amp; // damping
    (k, d)
}

/// The candidate's linearized pitch matrix at trim.
#[must_use]
pub fn pitch_model(c: &OrnithCandidate) -> [[f64; 2]; 2] {
    let (k, d) = pitch_coefficients(c);
    [[0.0, 1.0], [-k, -d]]
}

/// The nonlinear pitch vector field `(θ̇, −kθ + (k/θ_s²)θ³ − dθ̇)` in the
/// state `x = (θ, θ̇)`; its float coefficients DEFINE the model the
/// certificate is about.
#[must_use]
pub fn pitch_dynamics(c: &OrnithCandidate) -> Vec<MPoly> {
    let (k, d) = pitch_coefficients(c);
    pitch_field(k, d)
}

/// The pitch vector field for given stiffness and damping.
fn pitch_field(k: f64, d: f64) -> Vec<MPoly> {
    let th = MPoly::var(2, 0);
    let om = MPoly::var(2, 1);
    let k3 = k / (STALL_ANGLE * STALL_ANGLE);
    vec![
        om.clone(),
        th.scale(-k).add(&th.pow(3).scale(k3)).sub(&om.scale(d)),
    ]
}

/// Certify one candidate's trim state.
#[must_use]
pub fn certify(c: &OrnithCandidate) -> CertifyReport {
    // One BEM adjoint solve feeds both the linearization and the field.
    let (k, d) = pitch_coefficients(c);
    let a = [[0.0, 1.0], [-k, -d]];
    let maneuver = c.flap_amp * c.flap_freq / (k + 0.2);
    let roa = if k > 0.0 && d > 0.0 {
        // V = xᵀPx from the Lyapunov equation of the linearization (Q = I).
        // The post-stall saddles (±θ_s, 0) are equilibria, so no certified
        // ellipse can reach them: max θ on {xᵀPx ≤ c} is √(c·(P⁻¹)₁₁), giving
        // the strict ceiling c < θ_s²/(P⁻¹)₁₁. Starting the level search there
        // avoids paying for SDPs that must fail.
        let a_flat = [a[0][0], a[0][1], a[1][0], a[1][1]];
        solve_lyapunov(&a_flat, &[1.0, 0.0, 0.0, 1.0], 2).and_then(|p| {
            let det = p[0].mul_add(p[3], -(p[1] * p[2]));
            let pinv11 = p[3] / det;
            let ceiling = STALL_ANGLE * STALL_ANGLE / pinv11;
            let opts = RoaOptions {
                level_cap: ceiling,
                relative_tolerance: 1e-2,
                sdp: SdpSettings {
                    max_iter: 60,
                    ..SdpSettings::default()
                },
                ..RoaOptions::default()
            };
            (det > 0.0 && ceiling.is_finite() && ceiling > 0.0)
                .then(|| certify_roa_with(&pitch_field(k, d), &p, 1e-3, &opts).ok())
                .flatten()
        })
    } else {
        None
    };
    match roa {
        Some(r) => {
            let p = [
                [r.lyapunov[0], r.lyapunov[1]],
                [r.lyapunov[2], r.lyapunov[3]],
            ];
            let certified = lyapunov_certifies_stability(a, p);
            CertifyReport {
                a,
                p,
                certified,
                level: if certified { r.level } else { 0.0 },
                roa_volume: if certified { r.volume } else { 0.0 },
                sos_solves: r.sos_solves,
                maneuver,
            }
        }
        None => CertifyReport {
            a,
            p: [[0.0; 2]; 2],
            certified: false,
            level: 0.0,
            roa_volume: 0.0,
            sos_solves: 0,
            maneuver,
        },
    }
}

/// The screening surrogate with its conformal e-band: predict L/D from
/// the two dominant genes (thickness, alpha) with a fitted quadratic;
/// the band is split-conformal over held-out residuals — coverage is
/// GATED in the battery, and consumers must escalate outside the band.
pub struct LdSurrogate {
    coef: [f64; 6],
    /// The conformal band around predictions.
    pub band: ConformalBand,
}

impl LdSurrogate {
    /// Fit on a training set, calibrate the band on a held-out split.
    ///
    /// # Panics
    /// If fewer than 12 samples are supplied (6 coefficients + a
    /// calibration half need data).
    #[must_use]
    pub fn fit(samples: &[(OrnithCandidate, f64)], alpha: f64) -> LdSurrogate {
        assert!(samples.len() >= 12, "surrogate needs >= 12 samples");
        let half = samples.len() / 2;
        let (train, cal) = samples.split_at(half);
        // Least squares on [1, t, a, t², a², t·a] via normal equations.
        let feats = |c: &OrnithCandidate| -> [f64; 6] {
            let (t, a) = (c.thickness, c.alpha);
            [1.0, t, a, t * t, a * a, t * a]
        };
        let mut ata = [[0.0f64; 6]; 6];
        let mut atb = [0.0f64; 6];
        for (c, y) in train {
            let f = feats(c);
            for i in 0..6 {
                for j in 0..6 {
                    ata[i][j] += f[i] * f[j];
                }
                atb[i] += f[i] * y;
            }
        }
        // Ridge for conditioning (documented).
        for (i, row) in ata.iter_mut().enumerate() {
            row[i] += 1e-9;
        }
        let coef = solve6(&ata, &atb);
        let predict = |c: &OrnithCandidate| -> f64 {
            let f = feats(c);
            (0..6).map(|i| coef[i] * f[i]).sum()
        };
        let residuals: Vec<f64> = cal.iter().map(|(c, y)| y - predict(c)).collect();
        let band = conformal_band(&residuals, alpha);
        LdSurrogate { coef, band }
    }

    /// Predict L/D.
    #[must_use]
    pub fn predict(&self, c: &OrnithCandidate) -> f64 {
        let (t, a) = (c.thickness, c.alpha);
        let f = [1.0, t, a, t * t, a * a, t * a];
        (0..6).map(|i| self.coef[i] * f[i]).sum()
    }

    /// Empirical coverage of the band on fresh candidates.
    #[must_use]
    pub fn coverage(&self, fresh: &[OrnithCandidate]) -> f64 {
        let hits = fresh
            .iter()
            .filter(|c| self.band.covers(self.predict(c), lift_to_drag(c)))
            .count();
        hits as f64 / fresh.len().max(1) as f64
    }
}

/// Tiny dense 6×6 Gaussian elimination (fixture-scale).
fn solve6(a: &[[f64; 6]; 6], b: &[f64; 6]) -> [f64; 6] {
    let mut m = *a;
    let mut r = *b;
    for col in 0..6 {
        let mut piv = col;
        for row in col + 1..6 {
            if m[row][col].abs() > m[piv][col].abs() {
                piv = row;
            }
        }
        m.swap(col, piv);
        r.swap(col, piv);
        let d = m[col][col];
        assert!(d.abs() > 1e-30, "surrogate normal equations singular");
        let pivot_row = m[col];
        for row in col + 1..6 {
            let f = m[row][col] / d;
            for (cell, pivot) in m[row][col..].iter_mut().zip(pivot_row[col..].iter()) {
                *cell -= f * *pivot;
            }
            r[row] -= f * r[col];
        }
    }
    let mut x = [0.0f64; 6];
    for row in (0..6).rev() {
        let mut s = r[row];
        for (m_row_k, x_k) in m[row][row + 1..].iter().zip(x[row + 1..].iter()) {
            s -= *m_row_k * *x_k;
        }
        x[row] = s / m[row][row];
    }
    x
}
