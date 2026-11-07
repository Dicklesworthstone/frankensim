//! Physical Maxwell coefficients in energy-normalized elastic-strain coordinates.
//! No material admission, numerical integrator or eigenvalue solver lives here.

use fs_math::c64::C64;

#[derive(Clone, Copy, Debug)]
pub(super) struct MemoryCoefficient {
    pub stiffness: f64,
    pub root: f64,
    pub rate: f64,
}

impl MemoryCoefficient {
    fn new(stiffness: f64, rate: f64) -> Option<Self> {
        let root = stiffness.sqrt();
        if [stiffness, root, rate]
            .iter()
            .all(|x| x.is_finite() && *x > 0.0)
        {
            Some(Self {
                stiffness,
                root,
                rate,
            })
        } else {
            None
        }
    }

    // Peak phasor observation: z = sqrt(k) * (extension - viscous_extension).
    // The energy is |z|^2 / 4 and heat rate is |z|^2 / (2*tau).
    pub fn observe(self, extension: C64, z: C64) -> Option<MemoryObservation> {
        let elastic = z.scale(1.0 / self.root);
        let viscous = extension - elastic;
        let force = z.scale(self.root);
        let energy = z.scale(0.5).norm_sq();
        let heat = z
            .scale(self.rate.sqrt() * core::f64::consts::FRAC_1_SQRT_2)
            .norm_sq();
        if [extension, z, elastic, viscous, force]
            .iter()
            .any(|x| !x.re.is_finite() || !x.im.is_finite())
            || !energy.is_finite()
            || !heat.is_finite()
        {
            return None;
        }
        Some(MemoryObservation {
            elastic,
            viscous,
            force,
            energy,
            heat,
        })
    }
}

pub(super) struct MemoryObservation {
    pub elastic: C64,
    pub viscous: C64,
    pub force: C64,
    pub energy: f64,
    pub heat: f64,
}

/// The original longitudinal delta_E contributes to three independent exact
/// reference strain traces: delta_E*A/L, 3*delta_E*I/L^3, delta_E*I/L^3.
/// A zero source arm is handled by the source owner, without inventing a state.
pub(super) fn longitudinal_coefficients(
    delta: f64,
    area: f64,
    second_moment: f64,
    length: f64,
    tau: f64,
) -> Option<[MemoryCoefficient; 3]> {
    longitudinal_coefficients_on_reference(delta, area, second_moment, length, tau, 1.0)
}

/// Same three geometric coordinates on the original Timoshenko interpolation.
/// Only symmetric bending curvature changes, by the original bending fraction.
pub(super) fn longitudinal_coefficients_on_reference(
    delta: f64,
    area: f64,
    second_moment: f64,
    length: f64,
    tau: f64,
    bending_fraction: f64,
) -> Option<[MemoryCoefficient; 3]> {
    if [delta, area, second_moment, length, tau]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || !bending_fraction.is_finite()
        || bending_fraction <= 0.0
        || bending_fraction > 1.0
    {
        return None;
    }
    let rate = 1.0 / tau;
    let axial = delta * area / length;
    let bending = delta * second_moment / length / length / length;
    Some([
        MemoryCoefficient::new(axial, rate)?,
        MemoryCoefficient::new(3.0 * bending * bending_fraction * bending_fraction, rate)?,
        MemoryCoefficient::new(bending, rate)?,
    ])
}

/// Constant engineering shear is -s*x1/(2*L) on the original reference
/// interpolation. Its integrated energy is kappa*delta_G*A*s^2*x1^2/(8*L).
/// The original shear source already owns the corrected-area rigidity [N].
pub(super) fn transverse_coefficient(
    rigidity: f64,
    length: f64,
    tau: f64,
    shear_fraction: f64,
) -> Option<MemoryCoefficient> {
    if [rigidity, length, tau, shear_fraction]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || shear_fraction > 1.0
    {
        return None;
    }
    MemoryCoefficient::new(
        (rigidity / length) * (0.25 * shear_fraction * shear_fraction),
        1.0 / tau,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn original_geometry_channels_and_source_rate_survive() {
        let coefficients = longitudinal_coefficients(5.0, 2.0, 3.0, 2.0, 0.25).unwrap();
        for (coefficient, stiffness) in coefficients.into_iter().zip([5.0, 5.625, 1.875]) {
            assert_eq!(coefficient.stiffness, stiffness);
            assert_eq!(coefficient.rate, 4.0);
            assert!((coefficient.root.powi(2) - stiffness).abs() < 2e-15);
        }
        for (delta, tau) in [(0.0, 1.0), (-1.0, 1.0), (1.0, 0.0), (1.0, f64::INFINITY)] {
            assert!(longitudinal_coefficients(delta, 2.0, 3.0, 2.0, tau).is_none());
        }
    }

    #[test]
    fn timoshenko_channels_use_one_fixed_strain_basis_and_independent_times() {
        let (young, shear, area, moment, length, correction) = (8.0, 2.0, 0.2, 0.01, 0.7, 0.8);
        let phi = 12.0 * young * moment / (correction * shear * area * length * length);
        let beta = 1.0 / (1.0 + phi);
        let s = phi / (1.0 + phi);
        let e =
            longitudinal_coefficients_on_reference(young, area, moment, length, 0.4, beta).unwrap();
        let g = transverse_coefficient(correction * shear * area, length, 0.23, s).unwrap();
        // Reference bending+shear energy reconstructs the native stiffness.
        let expected = 3.0 * young * moment * beta / length.powi(3);
        assert!((e[1].stiffness + g.stiffness - expected).abs() < 2e-16);
        assert_eq!(e[0].stiffness, young * area / length);
        assert!((e[2].stiffness - young * moment / length.powi(3)).abs() < 1e-16);
        assert_eq!(e[1].rate, 1.0 / 0.4);
        assert_eq!(g.rate, 1.0 / 0.23);
        for omega in [1e-6, 1.0, 1.0 / 0.23, 1e6] {
            let laplace = C64::new(0.0, -omega);
            let x = C64::new(0.03, -0.01);
            let mut heat = 0.0;
            let mut force = C64::ZERO;
            for coefficient in [e[1], g] {
                let z = laplace * x.scale(coefficient.root)
                    / (laplace + C64::from_re(coefficient.rate));
                let seen = coefficient.observe(x, z).unwrap();
                heat += seen.heat;
                force = force + seen.force;
            }
            assert!((heat - 0.5 * (force.conj() * (laplace * x)).re).abs() < 3e-14);
        }
        assert!(transverse_coefficient(1.0, length, 0.23, 0.0).is_none());
        assert!(transverse_coefficient(1.0, length, 0.23, 1.1).is_none());
        assert!(
            longitudinal_coefficients_on_reference(young, area, moment, length, 0.4, 0.0).is_none()
        );
    }

    #[test]
    fn relaxation_changes_storage_and_loss_without_a_frozen_damping_fit() {
        let arm = longitudinal_coefficients(5.0, 2.0, 3.0, 2.0, 0.25).unwrap()[0];
        for omega in [1e-8, 0.4, 4.0, 40.0, 1e8] {
            let s = C64::new(0.0, -omega);
            let extension = C64::new(0.7, -0.2);
            let z = s * extension.scale(arm.root) / (s + C64::from_re(arm.rate));
            let seen = arm.observe(extension, z).unwrap();
            let ratio = omega / arm.rate;
            let dynamic_stiffness = C64::new(
                arm.stiffness * ratio * ratio / (1.0 + ratio * ratio),
                -arm.stiffness * ratio / (1.0 + ratio * ratio),
            );
            assert!((seen.force - dynamic_stiffness * extension).abs() < 2e-14);
            assert!((seen.elastic + seen.viscous - extension).abs() < 1e-15);
            let supplied = 0.5 * (seen.force.conj() * (s * extension)).re;
            assert!((supplied - seen.heat).abs() < 1e-14_f64.max(2e-8 * seen.heat));
            assert!((seen.heat - 2.0 * arm.rate * seen.energy).abs() < 3e-15);
        }
    }

    #[test]
    fn original_maxwell_cubic_has_all_three_coupled_rates() {
        // M=1, K_infinity=4, C=4.5, delta_K=.25 and tau=2/3.
        // (s+1.5)*(s^2+4.5*s+4)+.25*s=(s+1)*(s+2)*(s+3).
        let arm = longitudinal_coefficients(0.25, 1.0, 1.0, 1.0, 2.0 / 3.0).unwrap()[0];
        for s in [-1.0, -2.0, -3.0] {
            let z = arm.root * s / (s + arm.rate);
            let physical_row = s * s + 4.5 * s + 4.0 + arm.root * z;
            assert!(physical_row.abs() < 1e-14);
            assert!((s * z + arm.rate * z - arm.root * s).abs() < 1e-14);
        }
    }
}
