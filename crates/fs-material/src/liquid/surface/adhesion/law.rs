//! Ideal partial-wetting thermodynamics from explicit reversible adhesion work.
//!
//! Young--Dupre: W = sigma (1 + cos(theta)). The wetted-solid energy relative
//! to the dry substrate is (sigma-W) A_sl; liquid/gas energy is sigma A_lv.
//! This is an explicitly selected smooth, homogeneous, rigid-substrate model.
//! Fracture toughness, roughness, plastic peel work and advancing/receding
//! hysteresis measurements must not be substituted for equilibrium W.
//!
//! See Nagy, Scientific Reports 14, 29991 (2024), equations (1)-(3),
//! https://doi.org/10.1038/s41598-024-81710-6. No drop pull-off force is inferred.

use core::f64::consts::PI;
use fs_math::det;

/// Refusal of unsupported complete wetting/dewetting or nonrepresentable data.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum YoungDupreError {
    /// Surface tension and reversible work must be finite and positive.
    InvalidCoefficient,
    /// A finite-angle partial-wetting boundary requires 0 < W/sigma < 2.
    NoPartialWettingEquilibrium,
    /// The requested geometry or angle is outside the finite scalar domain.
    InvalidGeometry,
    /// A derived angle, force or energy cannot be represented as finite f64.
    Unrepresentable,
}
impl core::fmt::Display for YoungDupreError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Young-Dupre equilibrium refused: {self:?}")
    }
}
impl std::error::Error for YoungDupreError {}

/// A fixed-condition thermodynamic model, separate from its source receipts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YoungDupreLaw {
    work_j_m2: f64,
    tension_n_m: f64,
    angle_radians: f64,
}
impl YoungDupreLaw {
    /// Select the ideal model explicitly. Complete wetting/film formation and
    /// complete dewetting require different boundary owners; no clipping occurs.
    pub fn new(work_j_m2: f64, tension_n_m: f64) -> Result<Self, YoungDupreError> {
        if !work_j_m2.is_finite()
            || work_j_m2 <= 0.0
            || !tension_n_m.is_finite()
            || tension_n_m <= 0.0
        {
            return Err(YoungDupreError::InvalidCoefficient);
        }
        let ratio = work_j_m2 / tension_n_m;
        if !(ratio > 0.0 && ratio < 2.0) {
            return Err(YoungDupreError::NoPartialWettingEquilibrium);
        }
        // Half-angle form retains small angles near W=2 sigma and avoids
        // rounding W/sigma-1 to an endpoint before inverse trigonometry.
        let angle_radians = 2.0 * det::atan2((2.0 - ratio).sqrt(), ratio.sqrt());
        if !angle_radians.is_finite() || !(angle_radians > 0.0 && angle_radians < PI) {
            return Err(YoungDupreError::Unrepresentable);
        }
        Ok(Self {
            work_j_m2,
            tension_n_m,
            angle_radians,
        })
    }
    /// Original reversible work per unit area, J/m2.
    #[must_use]
    pub const fn work_of_adhesion_j_m2(self) -> f64 {
        self.work_j_m2
    }
    /// Original liquid/gas surface tension, N/m.
    #[must_use]
    pub const fn surface_tension_n_m(self) -> f64 {
        self.tension_n_m
    }
    /// Model-predicted equilibrium angle measured through liquid, radians.
    #[must_use]
    pub const fn contact_angle_radians(self) -> f64 {
        self.angle_radians
    }
    /// gamma_sl - gamma_sv, J/m2, relative to the initially dry substrate.
    #[must_use]
    pub fn wetted_surface_energy_difference_j_m2(self) -> f64 {
        self.tension_n_m - self.work_j_m2
    }
    /// Surface free energy [J] relative to the dry-substrate reference.
    /// Geometry owns the actual two areas; no VOF area or topology is inferred.
    pub fn surface_energy_j(
        self,
        liquid_gas_area_m2: f64,
        wetted_area_m2: f64,
    ) -> Result<f64, YoungDupreError> {
        if [liquid_gas_area_m2, wetted_area_m2]
            .into_iter()
            .any(|area| !area.is_finite() || area < 0.0)
        {
            return Err(YoungDupreError::InvalidGeometry);
        }
        let free = self.tension_n_m * liquid_gas_area_m2;
        let wet = self.wetted_surface_energy_difference_j_m2() * wetted_area_m2;
        let total = free + wet;
        if !free.is_finite() || !wet.is_finite() || !total.is_finite() {
            return Err(YoungDupreError::Unrepresentable);
        }
        Ok(total)
    }
    /// Unbalanced Young force per line length [N/m], positive for advancing
    /// liquid. This is a local thermodynamic driving force, not drop pull-off.
    pub fn advancing_force_per_length_n_m(
        self,
        geometric_angle_radians: f64,
    ) -> Result<f64, YoungDupreError> {
        if !geometric_angle_radians.is_finite() || !(0.0..=PI).contains(&geometric_angle_radians) {
            return Err(YoungDupreError::InvalidGeometry);
        }
        let force = self.work_j_m2 - self.tension_n_m * (1.0 + det::cos(geometric_angle_radians));
        if force.is_finite() {
            Ok(force)
        } else {
            Err(YoungDupreError::Unrepresentable)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol, "{a:e} != {b:e}");
    }
    #[test]
    fn adhesion_work_changes_angle_and_recovers_young_force_balance() {
        let mut previous = PI;
        for ratio in [0.02, 0.2, 0.7, 1.0, 1.6, 1.99] {
            let law = YoungDupreLaw::new(0.05 * ratio, 0.05).unwrap();
            let theta = law.contact_angle_radians();
            assert!(theta < previous);
            previous = theta;
            close(0.05 * (1.0 + theta.cos()), 0.05 * ratio, 2e-16);
            close(
                law.advancing_force_per_length_n_m(theta).unwrap(),
                0.0,
                2e-16,
            );
            assert!(law.advancing_force_per_length_n_m(0.5 * theta).unwrap() < 0.0);
            assert!(
                law.advancing_force_per_length_n_m(0.5 * (theta + PI))
                    .unwrap()
                    > 0.0
            );
        }
    }
    #[test]
    fn predicted_angle_minimizes_fixed_volume_spherical_cap_energy() {
        let law = YoungDupreLaw::new(0.08, 0.05).unwrap();
        let energy = |theta: f64| {
            let c = theta.cos();
            let radius = (3.0e-12 / (PI * (1.0 - c).powi(2) * (2.0 + c))).cbrt();
            law.surface_energy_j(
                2.0 * PI * radius * radius * (1.0 - c),
                PI * radius * radius * (1.0 - c * c),
            )
            .unwrap()
        };
        let theta = law.contact_angle_radians();
        let minimum = energy(theta);
        for offset in [-0.3, -0.1, 0.1, 0.3] {
            assert!(energy(theta + offset) > minimum);
        }
        let step = 1e-5;
        close(
            (energy(theta + step) - energy(theta - step)) / (2.0 * step),
            0.0,
            1e-18,
        );
    }
    #[test]
    fn surface_energy_scales_with_area_and_has_correct_wetting_derivative() {
        let law = YoungDupreLaw::new(0.08, 0.05).unwrap();
        let a = 2e-6;
        let b = 1e-6;
        close(law.surface_energy_j(a, b).unwrap(), 7e-8, 1e-22);
        close(
            law.surface_energy_j(7.0 * a, 7.0 * b).unwrap(),
            7.0 * law.surface_energy_j(a, b).unwrap(),
            1e-21,
        );
        let h = 1e-9;
        close(
            (law.surface_energy_j(a, b + h).unwrap() - law.surface_energy_j(a, b - h).unwrap())
                / (2.0 * h),
            -0.03,
            1e-14,
        );
    }
    #[test]
    fn unsupported_wetting_limits_and_nonfinite_geometry_refuse() {
        for (w, sigma) in [
            (0.0, 0.05),
            (-0.1, 0.05),
            (0.1, 0.05),
            (0.11, 0.05),
            (0.01, 0.0),
            (f64::NAN, 0.05),
            (0.01, f64::INFINITY),
        ] {
            assert!(YoungDupreLaw::new(w, sigma).is_err());
        }
        let law = YoungDupreLaw::new(0.08, 0.05).unwrap();
        assert!(law.surface_energy_j(-1.0, 0.0).is_err());
        assert!(law.surface_energy_j(0.0, f64::NAN).is_err());
        assert!(law.advancing_force_per_length_n_m(PI + 0.1).is_err());
        assert!(law.advancing_force_per_length_n_m(f64::NAN).is_err());
        let near = YoungDupreLaw::new(f64::from_bits(2.0_f64.to_bits() - 1), 1.0).unwrap();
        assert!(near.contact_angle_radians() > 0.0 && near.contact_angle_radians() < 1e-6);
    }
}
