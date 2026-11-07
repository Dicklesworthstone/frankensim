use super::*;
use crate::gas_species::input::read_gas_record;

fn law(text: &str) -> IsobaricGasLaw {
    IsobaricGasLaw::try_from_record(
        read_gas_record(&mut text.split_ascii_whitespace(), None, None).unwrap(),
    )
    .unwrap()
}

fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

fn values(gas: &GasRecord) -> [f64; 4] {
    let state = gas.state();
    [
        state.density,
        state.specific_heat_cp,
        state.thermal_conductivity,
        state.dynamic_viscosity,
    ]
}

fn derivatives(response: &IsobaricGasTransportTangent) -> [f64; 4] {
    [
        response.density_derivative_kg_m3_k,
        response.specific_heat_derivative_j_kg_k2,
        response.conductivity_derivative_w_m_k2,
        response.dynamic_viscosity_derivative_pa_s_k,
    ]
}

#[test]
fn g0_dry_air_preserves_values_and_original_sutherland_and_ideal_gas_slopes() {
    let law = law("300 100000 0");
    let spec = GasSpec::dry_air_ussa1976();
    for t in [50.0, 100.0, 300.0, 1000.0, 2000.0] {
        let response = law
            .transport_temperature_tangent(Temperature::new(t))
            .unwrap();
        assert_eq!(
            response.gas,
            law.at_temperature(Temperature::new(t)).unwrap()
        );
        let expected_mu = spec.sutherland_beta * t.sqrt() * (0.5 * t + 1.5 * spec.sutherland_s)
            / (t + spec.sutherland_s).powi(2);
        near(
            response.dynamic_viscosity_derivative_pa_s_k,
            expected_mu,
            4e-15,
        );
        near(
            response.density_derivative_kg_m3_k,
            -100000.0 * spec.molar_mass / (R_USSA_1976 * t * t),
            4e-15,
        );
        assert_eq!(response.specific_heat_derivative_j_kg_k2, 0.0);
        assert!(response.conductivity_derivative_w_m_k2 > 0.0);
    }
}

#[test]
fn g1_analytic_transport_response_matches_original_value_evaluations() {
    // Five-point differences are a test oracle only. In particular, the humid
    // samples retain water fraction, and Wilke/WMS denominators move with T.
    for (text, temperatures) in [
        ("300 100000 0", [100.0, 300.0, 1000.0]),
        ("300 100000 0.5", [300.0, 310.0, 315.0]),
        ("mixture 300 100000 mole 1 He 1", [280.0, 300.0, 315.0]),
        (
            "mixture 300 100000 mole 3 He 0.55 CO2 0.4 N2 0.05",
            [280.0, 300.0, 315.0],
        ),
        (
            "mixture 300 100000 mole 6 N2 0.60 O2 0.20 Ar 0.05 CO2 0.10 H2O 0.005 He 0.045",
            [285.0, 300.0, 315.0],
        ),
    ] {
        let law = law(text);
        for t in temperatures {
            let h = 0.02;
            let samples = [-2.0, -1.0, 1.0, 2.0]
                .map(|scale| values(&law.at_temperature(Temperature::new(t + scale * h)).unwrap()));
            let response = law
                .transport_temperature_tangent(Temperature::new(t))
                .unwrap();
            assert_eq!(
                response.gas,
                law.at_temperature(Temperature::new(t)).unwrap()
            );
            let exact = derivatives(&response);
            for i in [0, 2, 3] {
                let difference = (samples[0][i] - 8.0 * samples[1][i] + 8.0 * samples[2][i]
                    - samples[3][i])
                    / (12.0 * h);
                near(exact[i], difference, 3e-8);
            }
            assert_eq!(exact[1], 0.0);
            assert_eq!(
                response.gas.state().water_mole_fraction,
                law.at_temperature(Temperature::new(t))
                    .unwrap()
                    .state()
                    .water_mole_fraction
            );
        }
    }
}

#[test]
fn g3_pressure_and_species_order_keep_the_original_transport_meaning() {
    let low = law("mixture 300 50000 mole 3 He 0.55 CO2 0.4 N2 0.05");
    let high = law("mixture 300 100000 mole 3 He 0.55 CO2 0.4 N2 0.05");
    let permuted = law("mixture 300 100000 mole 3 N2 0.05 He 0.55 CO2 0.4");
    let a = low
        .transport_temperature_tangent(Temperature::new(310.0))
        .unwrap();
    let b = high
        .transport_temperature_tangent(Temperature::new(310.0))
        .unwrap();
    let c = permuted
        .transport_temperature_tangent(Temperature::new(310.0))
        .unwrap();
    assert_eq!(
        2.0 * a.density_derivative_kg_m3_k,
        b.density_derivative_kg_m3_k
    );
    assert_eq!(
        a.dynamic_viscosity_derivative_pa_s_k,
        b.dynamic_viscosity_derivative_pa_s_k
    );
    assert_eq!(
        a.conductivity_derivative_w_m_k2,
        b.conductivity_derivative_w_m_k2
    );
    for i in [0, 2, 3] {
        near(derivatives(&b)[i], derivatives(&c)[i], 3e-14);
    }
    let helium = law("mixture 300 100000 mole 1 He 1")
        .transport_temperature_tangent(Temperature::new(300.0))
        .unwrap();
    near(
        helium.conductivity_derivative_w_m_k2,
        3.75 * helium.gas.state().specific_gas_constant
            * helium.dynamic_viscosity_derivative_pa_s_k,
        4e-15,
    );
}

#[test]
fn g0_temperature_phase_and_original_source_limits_are_retained() {
    let humid = law("300 100000 0.9");
    let original = humid.clone();
    for t in [280.0, 0.0, f64::NAN, f64::INFINITY] {
        assert!(
            humid
                .transport_temperature_tangent(Temperature::new(t))
                .is_err()
        );
    }
    assert_eq!(humid, original);
    let named = law("mixture 300 100000 mole 1 N2 1");
    assert!(
        named
            .transport_temperature_tangent(Temperature::new(330.0))
            .is_err()
    );
    // No numerical perturbation leaves a source boundary: its inward branch
    // tangent remains available exactly at the model's temperature endpoints.
    for t in [273.15, 323.15] {
        named
            .transport_temperature_tangent(Temperature::new(t))
            .unwrap();
    }
    let saturated = law("300 50001 1");
    saturated
        .transport_temperature_tangent(Temperature::new(300.0))
        .unwrap();
    assert!(
        saturated
            .transport_temperature_tangent(Temperature::new(299.0))
            .is_err()
    );
}

#[test]
fn g0_thermochemical_values_do_not_imply_unsupported_transport_derivatives() {
    use crate::thermochemical_gas::tests::{model_in_domains, point};
    let model = model_in_domains(
        vec![0.8, 0.2],
        0.001,
        [250.0, 400.0],
        [290.0, 330.0],
        [50000.0, 150000.0],
    );
    let law = IsobaricGasLaw::try_from_model(model, point(300.0, 100000.0)).unwrap();
    let original = law.clone();
    law.at_temperature(Temperature::new(315.0)).unwrap();
    assert!(
        matches!(law.transport_temperature_tangent(Temperature::new(315.0)),
        Err(MaterialError::State { what }) if what.contains("original source temperature derivatives"))
    );
    assert_eq!(law, original);
}
