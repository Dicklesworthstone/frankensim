//! Lumped steady radiative balance for a declared incandescent filament.
//!
//! This is the cheapest admissible rung for a small filament in a highly
//! evacuated receiver: declared electrical power is balanced against
//! gray-body surface radiation. It intentionally does not infer resistance,
//! emissivity, geometry, lead conduction, gas conduction, or useful visible
//! output. Those remain explicit inputs or higher-fidelity work.

use crate::STEFAN_BOLTZMANN_W_M2_K4;

/// Inputs to the declared gray-body filament balance, all in SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IncandescentRadiativeInput {
    /// Applied terminal potential, V.
    pub voltage_v: f64,
    /// Declared operating-point filament resistance, ohm.
    pub hot_resistance_ohm: f64,
    /// Declared radiating surface area of the filament, m².
    pub radiating_area_m2: f64,
    /// Declared total hemispherical emissivity, dimensionless in `(0, 1]`.
    pub emissivity: f64,
    /// Declared radiative surroundings temperature, K.
    pub ambient_temperature_k: f64,
}

/// Admitted steady radiative balance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IncandescentRadiativeState {
    /// Applied terminal potential, V.
    pub voltage_v: f64,
    /// Current from Ohm's law at the declared resistance, A.
    pub current_a: f64,
    /// Electrical input `V²/R`, W.
    pub joule_power_w: f64,
    /// Gray-body equilibrium temperature, K.
    pub filament_temperature_k: f64,
    /// Re-evaluated outward gray-body surface radiation, W.
    pub radiative_power_w: f64,
    /// Relative closure `|P_electric - P_radiative| / max(P_electric, 1)`.
    pub relative_energy_closure: f64,
}

/// Typed refusal from the lumped filament boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IncandescentRadiativeError {
    /// At least one input was NaN or infinite.
    NonFinite,
    /// Voltage was negative.
    NegativeVoltage,
    /// Resistance was not strictly positive.
    NonPositiveResistance,
    /// Radiating area was not strictly positive.
    NonPositiveArea,
    /// Emissivity was outside `(0, 1]`.
    InvalidEmissivity,
    /// Ambient temperature was not strictly positive.
    NonPositiveAmbientTemperature,
}

/// Solve the steady balance
/// `V²/R = ε σ A (T⁴ - T_ambient⁴)`.
///
/// The returned temperature is a model result only for the declared inputs;
/// there is no silent material-property lookup or hidden convective loss.
pub fn solve_incandescent_radiative_balance(
    input: IncandescentRadiativeInput,
) -> Result<IncandescentRadiativeState, IncandescentRadiativeError> {
    let values = [
        input.voltage_v,
        input.hot_resistance_ohm,
        input.radiating_area_m2,
        input.emissivity,
        input.ambient_temperature_k,
    ];
    if values.iter().any(|value| !value.is_finite()) {
        return Err(IncandescentRadiativeError::NonFinite);
    }
    if input.voltage_v < 0.0 {
        return Err(IncandescentRadiativeError::NegativeVoltage);
    }
    if input.hot_resistance_ohm <= 0.0 {
        return Err(IncandescentRadiativeError::NonPositiveResistance);
    }
    if input.radiating_area_m2 <= 0.0 {
        return Err(IncandescentRadiativeError::NonPositiveArea);
    }
    if !(input.emissivity > 0.0 && input.emissivity <= 1.0) {
        return Err(IncandescentRadiativeError::InvalidEmissivity);
    }
    if input.ambient_temperature_k <= 0.0 {
        return Err(IncandescentRadiativeError::NonPositiveAmbientTemperature);
    }

    let current_a = input.voltage_v / input.hot_resistance_ohm;
    let joule_power_w = input.voltage_v * current_a;
    let denominator = input.emissivity * STEFAN_BOLTZMANN_W_M2_K4 * input.radiating_area_m2;
    let ambient_fourth = input.ambient_temperature_k.powi(4);
    let filament_temperature_k = (ambient_fourth + joule_power_w / denominator).sqrt().sqrt();
    let radiative_power_w =
        denominator * (filament_temperature_k.powi(4) - input.ambient_temperature_k.powi(4));
    let relative_energy_closure =
        (joule_power_w - radiative_power_w).abs() / joule_power_w.max(1.0);

    Ok(IncandescentRadiativeState {
        voltage_v: input.voltage_v,
        current_a,
        joule_power_w,
        filament_temperature_k,
        radiative_power_w,
        relative_energy_closure,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source_bounded_fixture() -> IncandescentRadiativeInput {
        IncandescentRadiativeInput {
            voltage_v: 110.0,
            hot_resistance_ohm: 150.0,
            // 22 cm long, 0.007 inch diameter cylindrical display filament.
            radiating_area_m2: core::f64::consts::PI * 0.000_177_8 * 0.22,
            emissivity: 0.8,
            ambient_temperature_k: 293.15,
        }
    }

    #[test]
    fn closes_declared_joule_to_gray_body_radiation() {
        let state = solve_incandescent_radiative_balance(source_bounded_fixture()).unwrap();
        assert!((state.current_a - 110.0 / 150.0).abs() < 1e-12);
        assert!(state.filament_temperature_k > 1_800.0);
        assert!(state.filament_temperature_k < 2_200.0);
        assert!(state.relative_energy_closure < 1e-12);
    }

    #[test]
    fn higher_voltage_raises_equilibrium_temperature() {
        let low = solve_incandescent_radiative_balance(IncandescentRadiativeInput {
            voltage_v: 70.0,
            ..source_bounded_fixture()
        })
        .unwrap();
        let high = solve_incandescent_radiative_balance(source_bounded_fixture()).unwrap();
        assert!(high.filament_temperature_k > low.filament_temperature_k);
        assert!(high.joule_power_w > low.joule_power_w);
    }

    #[test]
    fn refuses_unphysical_or_undeclared_domains() {
        assert_eq!(
            solve_incandescent_radiative_balance(IncandescentRadiativeInput {
                hot_resistance_ohm: 0.0,
                ..source_bounded_fixture()
            }),
            Err(IncandescentRadiativeError::NonPositiveResistance)
        );
        assert_eq!(
            solve_incandescent_radiative_balance(IncandescentRadiativeInput {
                emissivity: 1.2,
                ..source_bounded_fixture()
            }),
            Err(IncandescentRadiativeError::InvalidEmissivity)
        );
    }
}
