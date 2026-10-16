//! Mann-class spectral tensor (bead wf-root-guzez.4.6.2.1,
//! E3.3b-ii-a). Plan §5.4 Round-2: the atmosphere's statistical target
//! is a MANN-CLASS sheared von Kármán tensor — the shear distortion is
//! what produces the u–w cross-spectrum (uw < 0 surface stress) that a
//! diagonal-only model cannot represent (Round-3 Q2: diagonal-only may
//! initialize, never accept).
//!
//! Closed forms (Mann 1994, JFM 273, eqs. 3.20–3.24):
//!   E(k)  = α·ε^{2/3}·L^{5/3}·(kL)⁴ / (1 + (kL)²)^{17/6}
//!   β(k)  = Γ·(kL)^{−2/3} / √(₂F₁(1/3, 17/6; 4/3; −(kL)^{−2}))
//!   k₀    = (k₁, k₂, k₃ + β k₁),  ζ₁/ζ₂ from the C₁/C₂ algebra.
//! The hypergeometric is evaluated through the Pfaff transformation
//! (argument mapped to w = 1/(1+(kL)²) ∈ (0,1); convergent series,
//! deterministic fixed tolerance + term cap).

use crate::Refusal;
use fs_math::det;

/// Γ admission cap (published fits live near 3–4).
pub const MAX_GAMMA: f64 = 10.0;

/// kL admission floor (β and the series degrade below it).
pub const MIN_KL: f64 = 1.0e-4;

/// Hypergeometric term cap.
pub const HYPERGEOM_MAX_TERMS: usize = 4000;

/// The pinned Mann-class parameter tuple (enters the target artifact).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MannParams {
    /// α·ε^{2/3} [m^{4/3}/s²] (spectral level).
    pub alpha_eps23: f64,
    /// Length scale L [m].
    pub length_m: f64,
    /// Shear-distortion parameter Γ (0 = isotropic von Kármán).
    pub gamma: f64,
}

/// The registered neutral-surface-layer target at the Flyer reference
/// height class (mann-target-v1.json carries provenance; Estimated).
pub const MANN_TARGET_V1: MannParams = MannParams {
    alpha_eps23: 0.32,
    length_m: 1.77, // 0.59 * z_ref at z_ref = 3 m (neutral similarity)
    gamma: 3.9,
};

impl MannParams {
    /// Admit the tuple.
    ///
    /// # Errors
    /// `mann-params-invalid` (non-finite; non-positive level/length;
    /// Γ outside [0, MAX_GAMMA] — cap AND one ulp past).
    pub fn admit(&self) -> Result<(), Refusal> {
        let ok = self.alpha_eps23.is_finite()
            && self.alpha_eps23 > 0.0
            && self.length_m.is_finite()
            && self.length_m > 0.0
            && self.gamma.is_finite()
            && (0.0..=MAX_GAMMA).contains(&self.gamma);
        if !ok {
            return Err(Refusal {
                code: "mann-params-invalid",
                message: format!("{self:?}"),
                ranked_repairs: vec![format!("positive level/length; Gamma in [0, {MAX_GAMMA}]")],
            });
        }
        Ok(())
    }
}

/// ₂F₁(1/3, 17/6; 4/3; z) for z ≤ 0 via Pfaff:
/// (1−z)^{−1/3}·₂F₁(1/3, −3/2; 4/3; z/(z−1)); the mapped argument lies
/// in [0, 1) and the series converges (c − a − b' = 5/2 > 0 at w → 1).
///
/// # Errors
/// `hypergeom-did-not-converge` (term cap).
pub fn hypergeom_mann(z: f64) -> Result<f64, Refusal> {
    debug_assert!(z <= 0.0);
    let gauss = |a: f64, b: f64, c: f64, x: f64| -> Result<f64, Refusal> {
        let mut term = 1.0f64;
        let mut sum = 1.0f64;
        for n in 0..HYPERGEOM_MAX_TERMS {
            let nf = n as f64;
            term *= (a + nf) * (b + nf) / (c + nf) * x / (nf + 1.0);
            sum += term;
            if term.abs() < 1e-15 * sum.abs() {
                return Ok(sum);
            }
        }
        Err(Refusal {
            code: "hypergeom-did-not-converge",
            message: format!("series argument {x}, {HYPERGEOM_MAX_TERMS} terms"),
            ranked_repairs: vec!["argument mapping bug — both branches keep |x| <= 0.5".into()],
        })
    };
    let (a, b, c) = (1.0 / 3.0, 17.0 / 6.0, 4.0 / 3.0);
    if z >= -1.0 {
        // Pfaff: w = z/(z−1) ∈ [0, 0.5] here — fast.
        let w = z / (z - 1.0);
        Ok(det::pow(1.0 - z, -1.0 / 3.0) * gauss(a, c - b, c, w)?)
    } else {
        // z → 1/z connection formula (|1/z| < 1). Constants:
        // C1 = Γ(4/3)Γ(5/2)/(Γ(17/6)Γ(1)), C2 = Γ(4/3)Γ(−5/2)/(Γ(1/3)Γ(−3/2)) = −2/15
        // (cross-verified against direct series at z = −2, −50, −8870).
        const C1: f64 = 0.688_343_942_614_314_09;
        const C2: f64 = -2.0 / 15.0;
        let inv = 1.0 / z;
        let t1 = C1 * det::pow(-z, -a) * gauss(a, a - c + 1.0, a - b + 1.0, inv)?;
        let t2 = C2 * det::pow(-z, -b) * gauss(b, b - c + 1.0, b - a + 1.0, inv)?;
        Ok(t1 + t2)
    }
}

/// The von Kármán energy spectrum E(k).
#[must_use]
pub fn energy_spectrum(p: &MannParams, k: f64) -> f64 {
    let kl = k * p.length_m;
    p.alpha_eps23 * det::pow(p.length_m, 5.0 / 3.0) * det::pow(kl, 4.0)
        / det::pow(1.0 + kl * kl, 17.0 / 6.0)
}

/// The ISOTROPIC von Kármán tensor (independent oracle for the Γ → 0
/// limit — deliberately its own arithmetic, not the Γ = 0 path).
#[must_use]
pub fn isotropic_tensor(p: &MannParams, k: [f64; 3]) -> [[f64; 3]; 3] {
    let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
    let kn = det::sqrt(k2);
    let e = energy_spectrum(p, kn);
    let c = e / (4.0 * core::f64::consts::PI * k2 * k2);
    let mut phi = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { k2 } else { 0.0 };
            phi[i][j] = c * (delta - k[i] * k[j]);
        }
    }
    phi
}

/// The Mann sheared spectral tensor Φ_ij(k).
///
/// # Errors
/// Admission refusals; `mann-wavevector-invalid` (non-finite k or
/// kL below [`MIN_KL`] — floor AND one ulp under);
/// hypergeometric refusals pass through.
pub fn mann_tensor(p: &MannParams, k: [f64; 3]) -> Result<[[f64; 3]; 3], Refusal> {
    p.admit()?;
    let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
    let kn = det::sqrt(k2);
    let kl = kn * p.length_m;
    if !kn.is_finite() || kl < MIN_KL {
        return Err(Refusal {
            code: "mann-wavevector-invalid",
            message: format!("|k|L = {kl:?} below the floor {MIN_KL}"),
            ranked_repairs: vec!["evaluate above the admitted wavenumber floor".into()],
        });
    }
    // Eddy-lifetime shear parameter.
    let beta = if p.gamma == 0.0 {
        0.0
    } else {
        let f = hypergeom_mann(-1.0 / (kl * kl))?;
        p.gamma * det::pow(kl, -2.0 / 3.0) / det::sqrt(f)
    };
    let (k1, k2c, k3) = (k[0], k[1], k[2]);
    let k30 = k3 + beta * k1;
    let k0sq = k1 * k1 + k2c * k2c + k30 * k30;
    let k0n = det::sqrt(k0sq);
    let e0 = energy_spectrum(p, k0n);
    let kperp2 = k1 * k1 + k2c * k2c;
    let (zeta1, zeta2) = if beta == 0.0 {
        (0.0, 0.0)
    } else if k1.abs() < 1e-12 * kn {
        // Mann's k1 → 0 limit: ζ1 → −β, ζ2 → 0.
        (-beta, 0.0)
    } else {
        let c1 = beta * k1 * k1 * (k0sq - 2.0 * k30 * k30 + beta * k1 * k30) / (k2 * kperp2);
        let c2 = k2c * k0sq / det::pow(kperp2, 1.5)
            * det::atan2(beta * k1 * det::sqrt(kperp2), k0sq - k30 * k1 * beta);
        (c1 - k2c / k1 * c2, k2c / k1 * c1 + c2)
    };
    let common = e0 / (4.0 * core::f64::consts::PI * k0sq * k0sq);
    let phi11 = common * (k0sq - k1 * k1 - 2.0 * k1 * k30 * zeta1 + kperp2 * zeta1 * zeta1);
    let phi22 = common * (k0sq - k2c * k2c - 2.0 * k2c * k30 * zeta2 + kperp2 * zeta2 * zeta2);
    let phi33 = e0 / (4.0 * core::f64::consts::PI * k2 * k2) * kperp2;
    let phi12 =
        common * (-k1 * k2c - k1 * k30 * zeta2 - k2c * k30 * zeta1 + kperp2 * zeta1 * zeta2);
    let c13 = e0 / (4.0 * core::f64::consts::PI * k0sq * k2);
    let phi13 = c13 * (-k1 * k30 + kperp2 * zeta1);
    let phi23 = c13 * (-k2c * k30 + kperp2 * zeta2);
    Ok([
        [phi11, phi12, phi13],
        [phi12, phi22, phi23],
        [phi13, phi23, phi33],
    ])
}

/// Reynolds-stress integrals ∫Φ_ij d³k by deterministic log-spherical
/// quadrature (fixed grid — no randomness).
///
/// # Errors
/// Tensor refusals pass through.
pub fn stress_integrals(p: &MannParams) -> Result<[[f64; 3]; 3], Refusal> {
    p.admit()?;
    let (nr, nth, nph) = (96usize, 24usize, 24usize);
    let (r_lo, r_hi) = (1e-2 / p.length_m, 1e3 / p.length_m);
    let lr_lo = det::ln(r_lo);
    let lr_hi = det::ln(r_hi);
    let dlr = (lr_hi - lr_lo) / nr as f64;
    let dth = core::f64::consts::PI / nth as f64;
    let dph = core::f64::consts::TAU / nph as f64;
    let mut out = [[0.0f64; 3]; 3];
    for ir in 0..nr {
        let r = det::exp(lr_lo + (ir as f64 + 0.5) * dlr);
        for it in 0..nth {
            let th = (it as f64 + 0.5) * dth;
            let (sth, cth) = (det::sin(th), det::cos(th));
            for ip in 0..nph {
                let ph = (ip as f64 + 0.5) * dph;
                let k = [r * sth * det::cos(ph), r * sth * det::sin(ph), r * cth];
                let phi = mann_tensor(p, k)?;
                let w = r * r * r * sth * dlr * dth * dph; // r² sinθ · r dlr
                for i in 0..3 {
                    for j in 0..3 {
                        out[i][j] += phi[i][j] * w;
                    }
                }
            }
        }
    }
    Ok(out)
}
