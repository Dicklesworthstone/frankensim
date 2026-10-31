//! Volume-coordinate response of the ORIGINAL fixed BubbleGas inventory.
//!
//! The spherical reference radius only specifies its reference volume;
//! evaluating this law introduces no spherical liquid inertia or capillarity.
//! The work potential is minus the integral of gas pressure from V0 to V.
//! It is a barotropic mechanical potential, not generic gas internal energy,
//! heat exchange, latent heat or a finite-rate vaporization closure.

use super::{BubbleGas, CavitationError, finite, positive};

/// Gas pressure, positive compression stiffness and its conjugate work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasVolumeResponse {
    pub pressure_pa: f64,
    /// -dp_g/dV = kappa p_g/V [Pa/m³].
    pub stiffness_pa_per_m3: f64,
    /// -integral(V0..V) p_g(v) dv [J], with zero at the inventory's V0.
    pub mechanical_potential_j: f64,
}

impl BubbleGas {
    /// Author the SAME inventory with a physical reference volume [m³].
    /// No shape, extra surface pressure or gas thermal state is inferred.
    pub fn from_reference_volume_m3(
        volume_m3: f64,
        pressure_pa: f64,
        polytropic_exponent: f64,
    ) -> Result<Self, CavitationError> {
        if !positive(volume_m3) {
            return Err(CavitationError::InvalidInput(
                "positive gas reference volume",
            ));
        }
        let radius = (volume_m3 / (4.0 * core::f64::consts::PI / 3.0)).cbrt();
        let gas = Self::try_new(radius, pressure_pa, polytropic_exponent)?;
        gas.reference_volume_m3()?;
        Ok(gas)
    }

    /// Volume equivalent of the original inventory's reference radius.
    pub fn reference_volume_m3(self) -> Result<f64, CavitationError> {
        let radius = self.reference_radius_m;
        let volume = (radius * radius) * radius * (4.0 * core::f64::consts::PI / 3.0);
        if !positive(volume) {
            return Err(CavitationError::InvalidInput(
                "gas reference volume is not representable",
            ));
        }
        Ok(volume)
    }

    /// Exact polytropic pressure, analytic derivative, and analytic work.
    ///
    /// The kappa=1 limit is logarithmic. exp_m1 avoids subtracting nearly
    /// equal powers as kappa approaches one. Zero original gas inventory has
    /// exactly zero response; a positive inventory that underflows refuses.
    pub fn volume_response(self, volume_m3: f64) -> Result<GasVolumeResponse, CavitationError> {
        if !positive(volume_m3) {
            return Err(CavitationError::InvalidInput("positive gas pocket volume"));
        }
        if self.reference_pressure_pa == 0.0 {
            return Ok(GasVolumeResponse {
                pressure_pa: 0.0,
                stiffness_pa_per_m3: 0.0,
                mechanical_potential_j: 0.0,
            });
        }
        let reference = self.reference_volume_m3()?;
        let relative = (volume_m3 - reference) / reference;
        let logarithm = if relative.is_finite() && relative.abs() <= 0.5 {
            -relative.ln_1p()
        } else {
            reference.ln() - volume_m3.ln()
        };
        let exponent = self.polytropic_exponent * logarithm;
        let pressure = if exponent.abs() <= 100.0 {
            self.reference_pressure_pa * exponent.exp()
        } else {
            (self.reference_pressure_pa.ln() + exponent).exp()
        };
        let stiffness = (self.polytropic_exponent.ln() + pressure.ln() - volume_m3.ln()).exp();
        if !positive(pressure) || !positive(stiffness) {
            return Err(CavitationError::InvalidInput(
                "positive gas pressure or volume stiffness is not representable",
            ));
        }
        let reference_work = finite(
            self.reference_pressure_pa * reference,
            "reference gas pressure-volume work",
        )?;
        if reference_work == 0.0 {
            return Err(CavitationError::InvalidInput(
                "positive reference gas pressure-volume work underflows",
            ));
        }
        let delta = self.polytropic_exponent - 1.0;
        let integral_factor = if delta == 0.0 {
            logarithm
        } else {
            (delta * logarithm).exp_m1() / delta
        };
        let potential = finite(
            reference_work * integral_factor,
            "gas mechanical volume potential",
        )?;
        if logarithm != 0.0 && potential == 0.0 {
            return Err(CavitationError::InvalidInput(
                "nonzero gas mechanical volume work underflows",
            ));
        }
        Ok(GasVolumeResponse {
            pressure_pa: pressure,
            stiffness_pa_per_m3: stiffness,
            mechanical_potential_j: potential,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn near(a: f64, b: f64, tolerance: f64) {
        assert!((a - b).abs() <= tolerance, "{a:e} != {b:e} ({tolerance:e})");
    }

    #[test]
    fn radius_and_volume_are_the_same_gas_inventory() {
        let gas = BubbleGas::try_new(0.002, 80_000.0, 1.4).unwrap();
        for ratio in [0.4, 0.9, 1.0, 1.5, 4.0] {
            let radius = gas.reference_radius_m() * ratio;
            let volume = 4.0 * core::f64::consts::PI / 3.0 * radius.powi(3);
            near(
                gas.volume_response(volume).unwrap().pressure_pa,
                gas.pressure(radius).unwrap(),
                2e-8,
            );
        }
    }

    #[test]
    fn volume_stiffness_and_work_have_independent_finite_difference_oracles() {
        for kappa in [0.8, 1.0, 1.4, 2.0] {
            let gas = BubbleGas::from_reference_volume_m3(2e-6, 70_000.0, kappa).unwrap();
            let v = 1.3e-6;
            let h = 1e-11;
            let at = gas.volume_response(v).unwrap();
            let left = gas.volume_response(v - h).unwrap();
            let right = gas.volume_response(v + h).unwrap();
            near(
                (right.pressure_pa - left.pressure_pa) / (2.0 * h),
                -at.stiffness_pa_per_m3,
                3e-9 * at.stiffness_pa_per_m3,
            );
            near(
                (right.mechanical_potential_j - left.mechanical_potential_j) / (2.0 * h),
                -at.pressure_pa,
                3e-9 * at.pressure_pa,
            );
        }
    }

    #[test]
    fn near_isothermal_work_converges_without_cancellation() {
        let exact = BubbleGas::from_reference_volume_m3(1e-6, 1e5, 1.0).unwrap();
        let v0 = exact.reference_volume_m3().unwrap();
        for ratio in [0.2, 0.999999999, 1.000000001, 5.0] {
            let v = v0 * ratio;
            let expected = -1e5 * v0 * (v / v0).ln();
            // The independent ratio/log oracle and log differences round at the
            // reference pressure-volume energy scale, including near zero work.
            let roundoff = 8.0 * f64::EPSILON * (1e5 * v0).max(expected.abs());
            near(
                exact.volume_response(v).unwrap().mechanical_potential_j,
                expected,
                roundoff,
            );
            for kappa in [1.0 - 1e-10, 1.0 + 1e-10] {
                let gas = BubbleGas::from_reference_volume_m3(v0, 1e5, kappa).unwrap();
                near(
                    gas.volume_response(v).unwrap().mechanical_potential_j,
                    expected,
                    3e-11,
                );
            }
        }
    }

    #[test]
    fn zero_inventory_is_exact_and_positive_inventory_underflow_refuses() {
        let empty = BubbleGas::try_new(1.0, 0.0, 2.0).unwrap();
        for volume in [1e-300, 1.0, 1e300] {
            assert_eq!(
                empty.volume_response(volume).unwrap(),
                GasVolumeResponse {
                    pressure_pa: 0.0,
                    stiffness_pa_per_m3: 0.0,
                    mechanical_potential_j: 0.0
                }
            );
        }
        let tiny = BubbleGas::from_reference_volume_m3(1e-6, 1e-300, 1.4).unwrap();
        assert!(tiny.volume_response(1e300).is_err());
        let negligible_work = BubbleGas::from_reference_volume_m3(1e-30, 1e-300, 1.4).unwrap();
        assert!(negligible_work.volume_response(1e-30).is_err());
        assert!(empty.volume_response(0.0).is_err());
    }
}
