use super::*;

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn retains_original_phase_endpoints_and_caloric_reference() {
    let pair = If97SaturationPair::at_temperature(Temperature::new(373.15)).unwrap();
    let phase = SaturatedWaterVapor::try_from_saturation(pair).unwrap();
    assert_eq!(*phase.saturation(), pair);
    assert_eq!(phase.pressure_pa(), pair.pressure().value());
    assert_eq!(phase.vapor_density_kg_m3(), pair.vapor().density().value());
    assert_eq!(
        phase.liquid_density_kg_m3(),
        pair.liquid().density().value()
    );
    assert_eq!(
        phase.vapor_internal_energy_j_kg(),
        pair.vapor().specific_internal_energy().value()
    );
    assert_eq!(
        phase.liquid_enthalpy_j_kg(),
        pair.liquid().specific_enthalpy().value()
    );
    assert_eq!(phase.latent_enthalpy_j_kg(), pair.latent_enthalpy().value());
    assert_eq!(
        phase.vapor_sound_speed_m_s(),
        pair.vapor().sound_speed().value()
    );
}

#[test]
fn analytic_saturation_tangents_match_independent_source_differences() {
    let dt = 1e-3;
    for t in [273.2, 293.15, 373.15, 523.15, 600.0] {
        let point = SaturatedWaterVapor::at_temperature_k(t).unwrap();
        let lower = If97SaturationPair::at_temperature(Temperature::new(t - dt)).unwrap();
        let upper = If97SaturationPair::at_temperature(Temperature::new(t + dt)).unwrap();
        relative(
            point.pressure_temperature_derivative_pa_k(),
            (upper.pressure().value() - lower.pressure().value()) / (2.0 * dt),
            2e-7,
        );
        relative(
            point.vapor_density_temperature_derivative_kg_m3_k(),
            (upper.vapor().density().value() - lower.vapor().density().value()) / (2.0 * dt),
            2e-7,
        );
        relative(
            point.vapor_internal_energy_temperature_derivative_j_kg_k(),
            (upper.vapor().specific_internal_energy().value()
                - lower.vapor().specific_internal_energy().value())
                / (2.0 * dt),
            2e-7,
        );
    }
}

#[test]
fn phase_capacity_matches_open_system_energy_change() {
    let dt = 1e-3;
    for t in [280.0, 373.15, 500.0, 600.0] {
        let point = SaturatedWaterVapor::at_temperature_k(t).unwrap();
        let lower = If97SaturationPair::at_temperature(Temperature::new(t - dt)).unwrap();
        let upper = If97SaturationPair::at_temperature(Temperature::new(t + dt)).unwrap();
        let (lo, hi) = (lower.vapor(), upper.vapor());
        let energy = hi.density().value() * hi.specific_internal_energy().value()
            - lo.density().value() * lo.specific_internal_energy().value();
        let carried = point.liquid_enthalpy_j_kg() * (hi.density().value() - lo.density().value());
        relative(
            point.phase_temperature_capacity_j_m3_k(),
            (energy - carried) / (2.0 * dt),
            2e-7,
        );
    }
}

#[test]
fn fixed_temperature_expansion_requires_original_latent_heat() {
    // One unit of added vapor volume: first law Q = u_v dm - h_l dm + p dV.
    for t in [293.15, 373.15, 573.15] {
        let point = SaturatedWaterVapor::at_temperature_k(t).unwrap();
        let heat_per_volume = (point.vapor_internal_energy_j_kg() - point.liquid_enthalpy_j_kg())
            * point.vapor_density_kg_m3()
            + point.pressure_pa();
        relative(
            heat_per_volume,
            point.vapor_density_kg_m3() * point.latent_enthalpy_j_kg(),
            5e-15,
        );
    }
}

#[test]
fn preserves_original_saturation_domain_without_critical_extrapolation() {
    for t in [273.15, 623.15] {
        let point = SaturatedWaterVapor::at_temperature_k(t).unwrap();
        assert!(point.phase_temperature_capacity_j_m3_k() > 0.0);
        assert_eq!(point.temperature_k(), t);
    }
    for t in [f64::NAN, f64::INFINITY, 273.14, 623.16, 647.096] {
        assert!(SaturatedWaterVapor::at_temperature_k(t).is_err());
    }
}
