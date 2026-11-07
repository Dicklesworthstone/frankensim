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
    /// A derived current, power, or temperature cannot be represented in `f64`.
    Unrepresentable,
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
    if !current_a.is_finite()
        || !joule_power_w.is_finite()
        || (input.voltage_v > 0.0 && joule_power_w == 0.0)
    {
        return Err(IncandescentRadiativeError::Unrepresentable);
    }

    let ambient = input.ambient_temperature_k;
    let (filament_temperature_k, radiative_power_w) = if joule_power_w == 0.0 {
        // Equilibrium needs no powers of temperature, even when T^4 itself
        // would overflow or underflow. Preserve the supplied bath exactly.
        (ambient, 0.0)
    } else {
        let root = coefficient_fourth_root(input.emissivity, input.radiating_area_m2);
        let heating_temperature = joule_power_w.sqrt().sqrt() / root;
        let temperature = if heating_temperature <= ambient {
            let ratio4 = fs_math::det::powi(heating_temperature / ambient, 4);
            let q = (1.0 + ratio4).sqrt().sqrt();
            // q - 1 = ratio4 / ((q + 1)(q^2 + 1)). Form the small rise
            // before adding it to the bath instead of subtracting near-ones.
            ambient + ambient * (ratio4 / ((q + 1.0) * (q * q + 1.0)))
        } else {
            let ratio4 = fs_math::det::powi(ambient / heating_temperature, 4);
            heating_temperature * (1.0 + ratio4).sqrt().sqrt()
        };
        (
            temperature,
            gray_radiative_power_w(root, temperature, ambient),
        )
    };
    let relative_energy_closure =
        (joule_power_w - radiative_power_w).abs() / joule_power_w.max(1.0);
    if !filament_temperature_k.is_finite()
        || !radiative_power_w.is_finite()
        || !relative_energy_closure.is_finite()
    {
        return Err(IncandescentRadiativeError::Unrepresentable);
    }

    Ok(IncandescentRadiativeState {
        voltage_v: input.voltage_v,
        current_a,
        joule_power_w,
        filament_temperature_k,
        radiative_power_w,
        relative_energy_closure,
    })
}

// Taking roots before multiplying preserves coefficients whose product is
// subnormal or smaller than f64 while their equilibrium temperature is finite.
fn coefficient_fourth_root(emissivity: f64, area_m2: f64) -> f64 {
    emissivity.sqrt().sqrt() * STEFAN_BOLTZMANN_W_M2_K4.sqrt().sqrt() * area_m2.sqrt().sqrt()
}

// k(T^4-Ta^4), with k distributed over four factors. The temperature
// difference stays in kelvin until after subtraction, preserving close baths.
fn gray_radiative_power_w(root: f64, temperature_k: f64, ambient_k: f64) -> f64 {
    if temperature_k == ambient_k {
        return 0.0;
    }
    let t = root * temperature_k;
    let a = root * ambient_k;
    (root * (temperature_k - ambient_k)) * (t + a) * (t * t + a * a)
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
    fn g0_zero_drive_preserves_bath_even_when_its_fourth_power_is_unrepresentable() {
        for ambient_temperature_k in [1.0e-200, 293.15, 1.0e100, f64::MAX] {
            let state = solve_incandescent_radiative_balance(IncandescentRadiativeInput {
                voltage_v: 0.0,
                ambient_temperature_k,
                ..source_bounded_fixture()
            })
            .unwrap();
            assert_eq!(state.filament_temperature_k, ambient_temperature_k);
            assert_eq!(state.current_a, 0.0);
            assert_eq!(state.joule_power_w, 0.0);
            assert_eq!(state.radiative_power_w, 0.0);
            assert_eq!(state.relative_energy_closure, 0.0);
        }
    }

    #[test]
    fn g1_near_equilibrium_retains_small_net_radiation() {
        let ambient = 300.0;
        let target = ambient + 1.0e-8;
        let power = STEFAN_BOLTZMANN_W_M2_K4
            * (target - ambient)
            * (target + ambient)
            * (target * target + ambient * ambient);
        let state = solve_incandescent_radiative_balance(IncandescentRadiativeInput {
            voltage_v: power.sqrt(),
            hot_resistance_ohm: 1.0,
            radiating_area_m2: 1.0,
            emissivity: 1.0,
            ambient_temperature_k: ambient,
        })
        .unwrap();
        assert_eq!(state.filament_temperature_k, target);
        assert!((state.radiative_power_w / power - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn g1_representable_balance_survives_intermediate_power_range_limits() {
        for (emissivity, radiating_area_m2) in [(1.0, 1.0e-305), (1.0e-200, 1.0e-200)] {
            let state = solve_incandescent_radiative_balance(IncandescentRadiativeInput {
                voltage_v: 1.0,
                hot_resistance_ohm: 1.0,
                radiating_area_m2,
                emissivity,
                ambient_temperature_k: 300.0,
            })
            .unwrap();
            // Independent logarithmic fourth-root reference. The bath's
            // contribution is far below f64 resolution in these fixtures.
            let expected = fs_math::det::exp(
                -0.25
                    * (fs_math::det::ln(emissivity)
                        + fs_math::det::ln(STEFAN_BOLTZMANN_W_M2_K4)
                        + fs_math::det::ln(radiating_area_m2)),
            );
            assert!((state.filament_temperature_k / expected - 1.0).abs() < 1.0e-12);
            assert!((state.radiative_power_w - 1.0).abs() < 1.0e-12);
        }
    }

    #[test]
    fn g3_equal_emitting_area_emissivity_products_preserve_temperature() {
        let input = source_bounded_fixture();
        let base = solve_incandescent_radiative_balance(input).unwrap();
        let scaled = solve_incandescent_radiative_balance(IncandescentRadiativeInput {
            emissivity: input.emissivity * 0.5,
            radiating_area_m2: input.radiating_area_m2 * 2.0,
            ..input
        })
        .unwrap();
        assert!(
            (base.filament_temperature_k / scaled.filament_temperature_k - 1.0).abs() < 1.0e-14
        );
        assert_eq!(base.joule_power_w, scaled.joule_power_w);
    }

    #[test]
    fn g0_unrepresentable_electrical_outputs_refuse() {
        for (voltage_v, hot_resistance_ohm) in [(f64::MAX, 0.5), (1.0e200, 1.0), (1.0e-200, 1.0)] {
            assert_eq!(
                solve_incandescent_radiative_balance(IncandescentRadiativeInput {
                    voltage_v,
                    hot_resistance_ohm,
                    ..source_bounded_fixture()
                }),
                Err(IncandescentRadiativeError::Unrepresentable)
            );
        }
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
