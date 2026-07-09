//! Stage 4: TRIM AND CERTIFY. A 2-state pitch model per candidate
//! (θ̈ = −k(c)·θ − d(c)·θ̇, stiffness from the BEM lift slope,
//! damping from thickness — the smoke reduced model; Koopman/DMD
//! reduced models with per-trim conformal e-bands are the recorded
//! successor). Stability is not a vibe: fs-sos certifies the Lyapunov
//! matrix (A ᵀP + PA ≺ 0), and the CERTIFIED region-of-attraction
//! proxy is the volume of the P-ellipsoid level set under the trim
//! saturation bound. The screening surrogate carries a DISTRIBUTION-
//! FREE conformal band (fs-surrogate) — certify-or-escalate, gated on
//! coverage.

use crate::param::OrnithCandidate;
use crate::screen::lift_to_drag;
use fs_sos::lyapunov_certifies_stability;
use fs_surrogate::{ConformalBand, conformal_band};

/// The stability certificate row.
#[derive(Debug, Clone)]
pub struct CertifyReport {
    /// Pitch dynamics matrix (companion form).
    pub a: [[f64; 2]; 2],
    /// The certified Lyapunov matrix (when certified).
    pub p: [[f64; 2]; 2],
    /// SOS/Lyapunov certificate verified.
    pub certified: bool,
    /// Certified ROA proxy volume (P-ellipsoid area at the saturation
    /// level; 0.0 when uncertified — never pretended).
    pub roa_volume: f64,
    /// Maneuver proxy: control authority over pitch stiffness.
    pub maneuver: f64,
}

/// The candidate's pitch model: stiffness from the lift slope at trim,
/// damping from section thickness (thicker = more damped, the smoke
/// law).
#[must_use]
pub fn pitch_model(c: &OrnithCandidate) -> [[f64; 2]; 2] {
    let foil = c.section(crate::screen::PANELS);
    let dcl = fs_bem::panel2d::dcl_dalpha_adjoint(&foil, c.alpha);
    let k = 0.4 * dcl; // restoring stiffness ∝ lift slope
    let d = 8.0 * c.thickness + 0.4 * c.flap_amp; // damping
    [[0.0, 1.0], [-k, -d]]
}

/// Certify one candidate's trim state.
#[must_use]
pub fn certify(c: &OrnithCandidate) -> CertifyReport {
    let a = pitch_model(c);
    // Candidate Lyapunov matrix from the (2×2) Lyapunov equation with
    // Q = I, solved in closed form for the companion structure.
    let (k, d) = (-a[1][0], -a[1][1]);
    // For A = [[0,1],[−k,−d]], P = [[p11,p12],[p12,p22]] solving
    // AᵀP + PA = −I in closed form:
    //   −2k·p12 = −1            → p12 = 1/(2k)
    //   2(p12 − d·p22) = −1     → p22 = (p12 + 1/2)/d
    //   p11 = d·p12 + k·p22     (off-diagonal stationarity)
    let p12 = 1.0 / (2.0 * k);
    let p22 = (p12 + 0.5) / d;
    let p11 = d.mul_add(p12, k * p22);
    let p = [[p11, p12], [p12, p22]];
    let certified = k > 0.0 && d > 0.0 && lyapunov_certifies_stability(a, p);
    // ROA proxy: the ellipsoid {xᵀPx ≤ c*} inside the pitch saturation
    // |θ| ≤ 0.35 rad: c* = 0.35²/ (P⁻¹)₁₁-normalized bound; area =
    // π·c*/√det P.
    let roa_volume = if certified {
        let det = p11.mul_add(p22, -(p12 * p12));
        let cstar = 0.35 * 0.35 * det / p22; // sup |θ| on the ellipse
        std::f64::consts::PI * cstar / fs_math::det::sqrt(det)
    } else {
        0.0
    };
    let maneuver = c.flap_amp * c.flap_freq / (k + 0.2);
    CertifyReport {
        a,
        p,
        certified,
        roa_volume,
        maneuver,
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
