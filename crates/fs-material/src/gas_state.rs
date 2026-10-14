//! Shared admission of the numerical gas snapshot, independent of a gas name.
//!
//! Constitutive/source validity still belongs to the model that produced the
//! snapshot. In particular, water mole fraction is not relative humidity and
//! does not establish saturation, phase stability, or ISO 9613 applicability.

use crate::MaterialError;
use crate::fluid_transport::FluidTransportState;
use crate::gas::GasState;

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

fn close(actual: f64, expected: f64) -> bool {
    positive(expected)
        && (actual - expected).abs() <= 2.0e-12 * actual.abs().max(expected.abs())
}

impl GasState {
    /// Check one complete ideal-gas snapshot before passing it between solvers.
    ///
    /// All public fields must describe the SAME state: EOS density, frozen
    /// sound speed, local cp/gamma, transport groups and acoustic impedance
    /// cannot be independently edited. Re-evaluate the owning material model
    /// when temperature, pressure, composition or source data changes.
    ///
    /// This retains the existing GasState temperature/pressure envelope but
    /// admits any explicit water mole fraction in [0, 1]. Humid-air correlation
    /// limits belong to that correlation, not to generic gas mechanics. This
    /// numerical check grants no source authority or phase-stability claim.
    ///
    /// # Errors
    /// Refuses nonphysical, out-of-envelope, inconsistent or unrepresentable
    /// states. No coefficients are silently repaired or replaced by air data.
    pub fn validate_coherence(&self) -> Result<(), MaterialError> {
        if ![
            self.temperature,
            self.pressure,
            self.density,
            self.sound_speed,
            self.dynamic_viscosity,
            self.thermal_conductivity,
            self.gamma,
            self.specific_gas_constant,
            self.specific_heat_cp,
            self.prandtl,
            self.characteristic_impedance,
        ]
        .into_iter()
        .all(positive)
            || !(50.0..=2000.0).contains(&self.temperature)
            || self.pressure > 1.0e7
            || self.gamma <= 1.0
            || self.gamma > 5.0 / 3.0 + 1.0e-9
            || !self.water_mole_fraction.is_finite()
            || !(0.0..=1.0).contains(&self.water_mole_fraction)
        {
            return Err(MaterialError::Parameters {
                what: "gas snapshot is nonphysical or outside the GasState envelope".into(),
            });
        }
        let transport = FluidTransportState::try_new(
            self.density,
            self.specific_heat_cp,
            self.thermal_conductivity,
            self.dynamic_viscosity,
        )?;
        let rt = self.specific_gas_constant * self.temperature;
        if !close(self.density, self.pressure / rt)
            || !close(self.sound_speed, fs_math::det::sqrt(self.gamma * rt))
            || !close(
                self.specific_heat_cp,
                self.gamma * self.specific_gas_constant / (self.gamma - 1.0),
            )
            || !close(self.prandtl, transport.prandtl())
            || !close(self.characteristic_impedance, self.density * self.sound_speed)
        {
            return Err(MaterialError::Parameters {
                what: "gas snapshot contains inconsistent EOS, caloric, transport or acoustic coefficients".into(),
            });
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gas::GasSpec;
    use fs_qty::semantic::{Composition, CompositionBasis};

    fn air() -> GasState {
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100000.0).unwrap()
    }

    #[test]
    fn evaluated_dry_humid_and_steam_states_share_admission() {
        air().validate_coherence().unwrap();
        GasState::try_new_moist_air(300.0, 100000.0, 0.8)
            .unwrap().validate_coherence().unwrap();
        let composition = Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap();
        let steam = GasState::try_new_mixture(
            &[GasSpec::water_vapor_nist()], &composition, 450.0, 100000.0, Some(0),
        ).unwrap();
        assert_eq!(steam.water_mole_fraction, 1.0);
        steam.validate_coherence().unwrap();
    }

    #[test]
    fn changing_a_public_coefficient_cannot_create_two_different_media() {
        let corruptions: [fn(&mut GasState); 11] = [
            |s| s.temperature *= 2.0,
            |s| s.pressure *= 2.0,
            |s| s.density *= 2.0,
            |s| s.sound_speed *= 2.0,
            |s| s.dynamic_viscosity *= 2.0,
            |s| s.thermal_conductivity *= 2.0,
            |s| s.gamma = 1.2,
            |s| s.specific_gas_constant *= 2.0,
            |s| s.specific_heat_cp *= 2.0,
            |s| s.prandtl *= 2.0,
            |s| s.characteristic_impedance *= 2.0,
        ];
        for corrupt in corruptions {
            let mut gas = air();
            corrupt(&mut gas);
            assert!(gas.validate_coherence().is_err(), "accepted incoherent {gas:?}");
        }
    }

    #[test]
    fn invalid_composition_and_nonfinite_values_refuse() {
        for fraction in [-1.0, 1.1, f64::NAN, f64::INFINITY] {
            let mut gas = air();
            gas.water_mole_fraction = fraction;
            assert!(gas.validate_coherence().is_err());
        }
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            let mut gas = air();
            gas.density = value;
            assert!(gas.validate_coherence().is_err());
        }
    }

    #[test]
    fn coherent_pressure_change_retains_the_ideal_gas_sound_speed() {
        let first = air();
        let second = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 200000.0).unwrap();
        first.validate_coherence().unwrap();
        second.validate_coherence().unwrap();
        assert_eq!(second.density, 2.0 * first.density);
        assert_eq!(second.sound_speed, first.sound_speed);
        assert_eq!(second.characteristic_impedance, 2.0 * first.characteristic_impedance);
    }
}
