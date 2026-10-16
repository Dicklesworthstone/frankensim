//! Browser boundary for a declared US 223,898 filament radiative balance.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt::Write as _;
use fs_conduction::incandescent::{
    IncandescentRadiativeError, IncandescentRadiativeInput,
    solve_incandescent_radiative_balance,
};

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;

fn refusal_json(code: &str, message: &str, repair: &str) -> String {
    format!(
        "{{\"refusal\":{{\"code\":\"{code}\",\"message\":\"{message}\",\"repairs\":[\"{repair}\"]}}}}"
    )
}

/// Evaluate one steady declared Joule-to-gray-body filament balance.
#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[must_use]
pub fn edison_radiative_step(
    voltage_v: f64,
    hot_resistance_ohm: f64,
    radiating_area_m2: f64,
    emissivity: f64,
    ambient_temperature_k: f64,
) -> String {
    let state = match solve_incandescent_radiative_balance(IncandescentRadiativeInput {
        voltage_v,
        hot_resistance_ohm,
        radiating_area_m2,
        emissivity,
        ambient_temperature_k,
    }) {
        Ok(state) => state,
        Err(error) => {
            let (code, message, repair) = match error {
                IncandescentRadiativeError::NonFinite => (
                    "non-finite-input",
                    "every radiative-balance input must be finite",
                    "Supply finite SI values",
                ),
                IncandescentRadiativeError::NegativeVoltage => (
                    "negative-voltage",
                    "voltage magnitude cannot be negative",
                    "Supply a non-negative terminal-potential magnitude",
                ),
                IncandescentRadiativeError::NonPositiveResistance => (
                    "non-positive-resistance",
                    "hot filament resistance must be positive",
                    "Supply a declared positive operating-point resistance",
                ),
                IncandescentRadiativeError::NonPositiveArea => (
                    "non-positive-area",
                    "filament radiating area must be positive",
                    "Supply a declared positive surface area in square metres",
                ),
                IncandescentRadiativeError::InvalidEmissivity => (
                    "invalid-emissivity",
                    "gray-body emissivity must lie in (0, 1]",
                    "Supply a declared emissivity no greater than one",
                ),
                IncandescentRadiativeError::NonPositiveAmbientTemperature => (
                    "non-positive-ambient-temperature",
                    "radiative surroundings temperature must be positive",
                    "Supply ambient absolute temperature in kelvin",
                ),
            };
            return refusal_json(code, message, repair);
        }
    };

    let mut output = String::with_capacity(380);
    let _ = write!(
        output,
        "{{\"ok\":{{\
         \"voltage_v\":{},\
         \"current_a\":{},\
         \"joule_power_w\":{},\
         \"filament_temperature_k\":{},\
         \"radiative_power_w\":{},\
         \"relative_energy_closure\":{}\
         }}}}",
        state.voltage_v,
        state.current_a,
        state.joule_power_w,
        state.filament_temperature_k,
        state.radiative_power_w,
        state.relative_energy_closure,
    );
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_closed_radiative_balance() {
        let output = edison_radiative_step(
            110.0,
            150.0,
            core::f64::consts::PI * 0.000_177_8 * 0.22,
            0.8,
            293.15,
        );
        assert!(output.contains("\"ok\""));
        assert!(output.contains("\"filament_temperature_k\""));
        assert!(output.contains("\"relative_energy_closure\""));
    }

    #[test]
    fn voltage_changes_operating_point() {
        let low = edison_radiative_step(70.0, 150.0, 0.000_123, 0.8, 293.15);
        let high = edison_radiative_step(110.0, 150.0, 0.000_123, 0.8, 293.15);
        assert_ne!(low, high);
    }

    #[test]
    fn invalid_emissivity_refuses_without_ok_payload() {
        let output = edison_radiative_step(110.0, 150.0, 0.000_123, 1.2, 293.15);
        assert!(output.contains("\"code\":\"invalid-emissivity\""));
        assert!(!output.contains("\"ok\""));
    }
}
