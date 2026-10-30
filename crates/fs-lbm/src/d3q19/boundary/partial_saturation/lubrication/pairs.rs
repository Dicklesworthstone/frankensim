//! Objective unequal-sphere leading-logarithmic translation/rotation resistance.
//!
//! Townsend (2023), https://arxiv.org/abs/1802.08226, equations 19-25, 33-34a,
//! with the author's corrected coefficient functions G2ya/G2yb/G2yc/G4yc:
//! https://github.com/Pecnut/stokesian-dynamics/blob/HEAD/find_resistance_scalars/functions_nearfield.py
//! For lambda=b/a the dimensional leading coefficients are
//! YA11 = 6*pi*mu*a * 4*lambda*(2+lambda+2*lambda^2)/(15*(1+lambda)^3),
//! YB11 = 4*pi*mu*a^2 * -lambda*(4+lambda)/(5*(1+lambda)^2),
//! YC11 = 8*pi*mu*a^3 * 2*lambda/(5*(1+lambda)),
//! YC12 = 8*pi*mu*a^3 * lambda^2/(10*(1+lambda)).
//! These all multiply ln(1/xi). At contact their reciprocal block factors into
//! k*[[4,1],[1,4]] on radius-scaled spins relative to line-of-centers rotation.
//! The two eigenvalues are 5k and 3k: slip alone would omit a singular mode.
//!
//! Using the actual center distance in the frame rate extends that leading
//! block objectively to finite positive gaps. Its exact null rigid motions and
//! power-conjugate wrench preserve total angular momentum; departures from the
//! literal contact-limit coefficients are only omitted-order h*ln(h) terms.

use super::*;

pub(super) struct TangentialPairMode {
    pub(super) linear: [f64; 3],
    pub(super) angular_a: [f64; 3],
    pub(super) angular_b: [f64; 3],
    pub(super) resistance: f64,
}

pub(super) fn tangential_modes(
    radii: [f64; 2],
    distance: f64,
    normal: [f64; 3],
    gap: f64,
    scales: Option<(f64, f64)>,
    viscosity: f64,
) -> Result<Option<[TangentialPairMode; 4]>> {
    let Some((cutoff, roughness)) = scales else {
        return Ok(None);
    };
    let effective = gap.max(roughness);
    if effective >= cutoff {
        return Ok(None);
    }
    let small = radii[0].min(radii[1]);
    let large = radii[0].max(radii[1]);
    let radius = small / (1.0 + small / large);
    let k = (4.0 * core::f64::consts::PI / 5.0)
        * viscosity
        * radius
        * fs_math::det::ln(cutoff / effective);
    if !(k.is_finite() && k > 0.0) {
        return Err(invalid(
            "tangential sphere-pair resistance is unrepresentable",
        ));
    }
    let basis = tangent_basis(normal);
    let modes = core::array::from_fn(|i| {
        let tangent = basis[i / 2];
        let spin = cross(normal, tangent);
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        // The divided radii stay bounded for arbitrary admitted size ratios.
        // This is the derivative of the SAME actual line-of-centers frame rate
        // used in z_a,z_b, rather than a corrective torque added after solving.
        let linear =
            (radii[0] / distance + sign * (radii[1] / distance)) * core::f64::consts::FRAC_1_SQRT_2;
        TangentialPairMode {
            linear: tangent.map(|x| linear * x),
            angular_a: spin.map(|x| radii[0] * x * core::f64::consts::FRAC_1_SQRT_2),
            angular_b: spin.map(|x| sign * radii[1] * x * core::f64::consts::FRAC_1_SQRT_2),
            resistance: if sign > 0.0 { 5.0 * k } else { 3.0 * k },
        }
    });
    if modes.iter().any(|mode| {
        !mode.resistance.is_finite()
            || !mode
                .linear
                .into_iter()
                .chain(mode.angular_a)
                .chain(mode.angular_b)
                .all(f64::is_finite)
    }) {
        return Err(invalid("tangential sphere-pair modes overflowed"));
    }
    Ok(Some(modes))
}
