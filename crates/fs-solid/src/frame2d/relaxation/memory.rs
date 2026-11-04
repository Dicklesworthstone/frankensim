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
    if [delta, area, second_moment, length, tau]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return None;
    }
    let rate = 1.0 / tau;
    let axial = delta * area / length;
    let bending = delta * second_moment / length / length / length;
    Some([
        MemoryCoefficient::new(axial, rate)?,
        MemoryCoefficient::new(3.0 * bending, rate)?,
        MemoryCoefficient::new(bending, rate)?,
    ])
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
