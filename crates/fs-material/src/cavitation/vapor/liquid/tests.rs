use super::*;
use crate::cavitation::vapor::SaturatedWaterVapor;

#[test]
fn liquid_storage_and_phase_transfer_retain_the_same_source_enthalpy() {
    for t in [273.16, 293.15, 373.15, 450.0, 623.15] {
        let pair = If97SaturationPair::at_temperature(Temperature::new(t)).unwrap();
        let liquid = SaturatedWaterLiquidThermal::try_from_saturation(pair).unwrap();
        let vapor = SaturatedWaterVapor::try_from_saturation(pair).unwrap();
        assert_eq!(*liquid.saturation(), pair);
        assert_eq!(liquid.enthalpy_j_kg(), vapor.liquid_enthalpy_j_kg());
        assert_eq!(liquid.density_kg_m3(), pair.liquid().density().value());
        assert_eq!(liquid.cp_j_kg_k(), pair.liquid().cp().value());
        assert_eq!(liquid.saturation_pressure_pa(), pair.pressure().value());
        assert_eq!(liquid.temperature_k(), t);
    }
}

#[test]
fn analytic_liquid_caloric_tangent_matches_independent_saturation_differences() {
    let dt = 1e-3;
    for t in [273.2, 293.15, 373.15, 450.0, 600.0] {
        let point = SaturatedWaterLiquidThermal::at_temperature_k(t).unwrap();
        let lo = If97SaturationPair::at_temperature(Temperature::new(t - dt)).unwrap();
        let hi = If97SaturationPair::at_temperature(Temperature::new(t + dt)).unwrap();
        let difference = (hi.liquid().specific_enthalpy().value()
            - lo.liquid().specific_enthalpy().value())
            / (2.0 * dt);
        assert!(
            (point.enthalpy_temperature_derivative_j_kg_k() - difference).abs()
                < 2e-7 * difference.abs()
        );
    }
}

#[test]
fn actual_industrial_transport_and_pressure_path_defect_remain_visible() {
    for t in [293.15, 373.15, 600.0] {
        let state = SaturatedWaterLiquidThermal::at_temperature_k(t).unwrap();
        let original =
            If97TransportState::from_state(If97FluidState::Gibbs(*state.saturation().liquid()))
                .unwrap();
        assert_eq!(*state.transport(), original);
        assert_eq!(
            state.conductivity_w_m_k(),
            original.thermal_conductivity().value()
        );
        assert_eq!(
            state.dynamic_viscosity_pa_s(),
            original.dynamic_viscosity().value()
        );
        assert_eq!(
            state.sound_speed_m_s(),
            state.saturation().liquid().sound_speed().value()
        );
        let relative = (state.enthalpy_temperature_derivative_j_kg_k() - state.cp_j_kg_k()).abs()
            / state.cp_j_kg_k();
        assert!((relative - state.pressure_path_capacity_fraction()).abs() < 5e-16);
        assert!(state.isothermal_compressibility_per_pa() > 0.0);
    }
    let cool = SaturatedWaterLiquidThermal::at_temperature_k(293.15).unwrap();
    let hot = SaturatedWaterLiquidThermal::at_temperature_k(600.0).unwrap();
    assert!(cool.pressure_path_capacity_fraction() < 1e-3);
    assert!(hot.pressure_path_capacity_fraction() > cool.pressure_path_capacity_fraction());
}

#[test]
fn original_temperature_domain_refuses_extrapolated_liquid_calorics() {
    for t in [f64::NAN, f64::INFINITY, 273.14, 273.15, 623.16, 647.096] {
        assert!(SaturatedWaterLiquidThermal::at_temperature_k(t).is_err());
    }
}
