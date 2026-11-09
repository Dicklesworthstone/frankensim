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
//! The hypergeometric uses Pfaff and inverse-argument transformations,
//! keeping both series away from their convergence boundaries, with a
//! deterministic fixed tolerance and term cap.

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

/// ₂F₁(1/3, 17/6; 4/3; z) for z ≤ 0. Pfaff maps z ≥ −2 to
/// w = z/(z−1) ∈ [0, 2/3]; below −2 the connection formula uses
/// 1/z ∈ (−1/2, 0). Neither series approaches its convergence boundary.
/// Transformations: DLMF 15.8.1 and 15.8.2, <https://dlmf.nist.gov/15.8>.
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
            ranked_repairs: vec![
                "check argument mapping: Pfaff x in [0, 2/3], inverse x in (-1/2, 0)".into(),
            ],
        })
    };
    let (a, b, c) = (1.0 / 3.0, 17.0 / 6.0, 4.0 / 3.0);
    if z >= -2.0 {
        // Keep Pfaff through z = -2: switching at -1 sent the inverse
        // series arbitrarily close to -1, refusing valid modes near kL = 1.
        let w = z / (z - 1.0);
        Ok(det::pow(1.0 - z, -1.0 / 3.0) * gauss(a, c - b, c, w)?)
    } else {
        // z → 1/z connection formula (|1/z| < 1/2). Constants:
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

/// Reynolds stresses ∫Φ_ij d³k over 0.01 ≤ |k|L ≤ 1000.
///
/// Fixed Gauss–Legendre quadrature in log(|k|L), polar angle and a
/// symmetry-reduced azimuth resolves the sheared energy-containing scales.
/// The integration is dimensionless: level and length enter only through
/// α·ε^(2/3)·L^(2/3). The finite radial band is not extrapolated or rescaled.
///
/// # Errors
/// Parameter and hypergeometric refusals pass through.
pub fn stress_integrals(p: &MannParams) -> Result<[[f64; 3]; 3], Refusal> {
    p.admit()?;
    let radial = gauss_legendre::<64>();
    let lr_lo = det::ln(1e-2);
    let lr_hi = det::ln(1e3);
    let lr_mid = 0.5 * (lr_hi + lr_lo);
    let lr_half = 0.5 * (lr_hi - lr_lo);
    let scale = p.alpha_eps23 * det::pow(p.length_m, 2.0 / 3.0);
    let mut stress = [0.0; 4];
    let mut compensation = [0.0; 4];

    // For Γ = 0, the angular integral is exactly (2/3)∫E(k) dk
    // on each diagonal. Do not manufacture small anisotropy or uw by
    // numerically integrating an analytically isotropic angular factor.
    if p.gamma == 0.0 {
        for (node, weight) in radial {
            let r = det::exp(lr_mid + lr_half * node);
            let r2 = r * r;
            let value =
                weight * lr_half * (2.0 / 3.0) * r2 * r2 * r / det::pow(1.0 + r2, 17.0 / 6.0);
            compensated_add(&mut stress[0], &mut compensation[0], value);
        }
        let variance = scale * stress[0];
        return Ok([
            [variance, 0.0, 0.0],
            [0.0, variance, 0.0],
            [0.0, 0.0, variance],
        ]);
    }

    let pi = core::f64::consts::PI;
    let polar = gauss_legendre::<48>().map(|(node, weight)| {
        let (sin, cos) = det::sin_cos(0.5 * pi * (node + 1.0));
        (sin, cos, weight * 0.5 * pi * sin)
    });
    let azimuth = gauss_legendre::<24>().map(|(node, weight)| {
        let (sin, cos) = det::sin_cos(0.25 * pi * (node + 1.0));
        // Φ11, Φ22, Φ33 and Φ13 are invariant under ky reflection and
        // full k inversion. These map all azimuths onto this quadrant.
        // Φ12 and Φ23 are odd under ky reflection and integrate to zero.
        (sin, cos, weight * pi)
    });
    for (node, weight) in radial {
        let r = det::exp(lr_mid + lr_half * node);
        let r2 = r * r;
        // Lifetime depends only on |k|L: evaluate its series once per
        // radius, not once per angular point. All dimensional powers
        // cancel before the spectral tensor is integrated.
        let beta = p.gamma * det::pow(r, -2.0 / 3.0) / det::sqrt(hypergeom_mann(-1.0 / r2)?);
        let radial_weight = weight * lr_half * r2 * r2 * r / (4.0 * pi);
        for (sth, n3, polar_weight) in polar {
            let np2 = sth * sth;
            for (sph, cph, azimuth_weight) in azimuth {
                let n1 = sth * cph;
                let n2 = sth * sph;
                let n30 = n3 + beta * n1;
                let n0sq = np2 + n30 * n30;
                // Same distortion algebra as mann_tensor, with |n| = 1.
                // Interior quadrature nodes avoid both polar and n1=0
                // singularities; sqrt(np2) is the positive polar sine.
                let c1 = beta * n1 * n1 * (n0sq - 2.0 * n30 * n30 + beta * n1 * n30) / np2;
                let c2 =
                    n2 * n0sq / (np2 * sth) * det::atan2(beta * n1 * sth, n0sq - n30 * n1 * beta);
                let z1 = c1 - n2 / n1 * c2;
                let z2 = n2 / n1 * c1 + c2;
                let tensor = [
                    n0sq - n1 * n1 - 2.0 * n1 * n30 * z1 + np2 * z1 * z1,
                    n0sq - n2 * n2 - 2.0 * n2 * n30 * z2 + np2 * z2 * z2,
                    n0sq * n0sq * np2,
                    n0sq * (-n1 * n30 + np2 * z1),
                ];
                // E(k0)/|k0|^4 cancels its (|k0|L)^4 numerator.
                // The remaining r^5 includes the r^3 log-radius Jacobian.
                let w = radial_weight * polar_weight * azimuth_weight
                    / det::pow(1.0 + r2 * n0sq, 17.0 / 6.0);
                for i in 0..4 {
                    compensated_add(&mut stress[i], &mut compensation[i], w * tensor[i]);
                }
            }
        }
    }
    let [u2, v2, w2, uw] = stress.map(|value| scale * value);
    Ok([[u2, 0.0, uw], [0.0, v2, 0.0], [uw, 0.0, w2]])
}

fn compensated_add(sum: &mut f64, compensation: &mut f64, value: f64) {
    let increment = value - *compensation;
    let next = *sum + increment;
    *compensation = (next - *sum) - increment;
    *sum = next;
}

/// Fixed-order Gauss–Legendre nodes and weights on [-1, 1]. Newton's
/// iteration is bounded, and each mirrored pair shares exactly one weight.
fn gauss_legendre<const N: usize>() -> [(f64, f64); N] {
    let polynomial = |x: f64| {
        let (mut previous, mut value) = (1.0, x);
        for degree in 2..=N {
            let d = degree as f64;
            let next = ((2.0 * d - 1.0) * x * value - (d - 1.0) * previous) / d;
            (previous, value) = (value, next);
        }
        let derivative = N as f64 * (x * value - previous) / (x * x - 1.0);
        (value, derivative)
    };
    let mut nodes = [(0.0, 0.0); N];
    for i in 0..N.div_ceil(2) {
        let mut root = det::cos(core::f64::consts::PI * (i as f64 + 0.75) / (N as f64 + 0.5));
        for _ in 0..16 {
            let (value, derivative) = polynomial(root);
            let correction = value / derivative;
            root -= correction;
            if correction.abs() <= 2.0 * f64::EPSILON {
                break;
            }
        }
        let (_, derivative) = polynomial(root);
        let weight = 2.0 / ((1.0 - root * root) * derivative * derivative);
        nodes[i] = (-root, weight);
        nodes[N - 1 - i] = (root, weight);
    }
    nodes
}
