use super::*;
use crate::{
    ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1, NASA9_STATE_SCHEMA_VERSION_V1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure};
use std::collections::BTreeMap;

fn region(bounds: [f64; 2], coefficients: [f64; 9]) -> Nasa9RegionV1 {
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    Nasa9RegionV1::from_card(ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", bounds[0], bounds[1]),
        sources: vec![],
        provenance: Provenance {
            source: "manufactured Cp polynomial".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    })
    .unwrap()
}

fn model(regions: Vec<Nasa9RegionV1>) -> Nasa9StandardStateModelV1 {
    Nasa9StandardStateModelV1::new(
        SpeciesId::new("A").unwrap(),
        MolarMass::new(0.028),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("manufactured-reference").unwrap(),
        )
        .unwrap(),
        regions,
    )
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 5e-9 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn original_nasa_terms_have_analytic_temperature_tangents_without_formation_reference_work() {
    let coefficients = [700.0, -20.0, 3.5, 0.003, -2e-6, 4e-9, 1e-12, 8e9, -3e7];
    let source = model(vec![region([200.0, 2000.0], coefficients)]);
    for t in [200.0, 300.0, 700.0, 2000.0] {
        let derivative = source
            .evaluate_heat_capacity_temperature_derivative(Temperature::new(t))
            .unwrap();
        let expected = R
            * (-1400.0 / t.powi(3) + 20.0 / t.powi(2) + 0.003 - 4e-6 * t
                + 12e-9 * t.powi(2)
                + 4e-12 * t.powi(3));
        near(derivative.value(), expected);
    }
    for t in [300.0, 700.0, 1500.0] {
        let h = 0.01;
        let values = [-2.0, -1.0, 1.0, 2.0].map(|scale| {
            source
                .evaluate_heat_capacity(Temperature::new(t + scale * h))
                .unwrap()
                .value()
        });
        let difference = (values[0] - 8.0 * values[1] + 8.0 * values[2] - values[3]) / (12.0 * h);
        near(
            source
                .evaluate_heat_capacity_temperature_derivative(Temperature::new(t))
                .unwrap()
                .value(),
            difference,
        );
    }
    let mut shifted = coefficients;
    shifted[7] *= -100.0;
    shifted[8] *= 200.0;
    let shifted = model(vec![region([200.0, 2000.0], shifted)]);
    assert_eq!(
        source.evaluate_heat_capacity_temperature_derivative(Temperature::new(600.0)),
        shifted.evaluate_heat_capacity_temperature_derivative(Temperature::new(600.0))
    );
}

#[test]
fn source_endpoints_and_smooth_seams_admit_but_cp_jumps_corners_and_gaps_refuse() {
    let a = [0.0, 0.0, 3.5, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0];
    let smooth = model(vec![region([200.0, 400.0], a), region([400.0, 800.0], a)]);
    for t in [200.0, 400.0, 800.0] {
        near(
            smooth
                .evaluate_heat_capacity_temperature_derivative(Temperature::new(t))
                .unwrap()
                .value(),
            R * 0.001,
        );
    }
    let mut jump = a;
    jump[2] += 0.25;
    let mut corner = a;
    corner[2] = 3.1;
    corner[3] = 0.002; // Same Cp at 400 K, unequal derivatives.
    for b in [jump, corner] {
        let source = model(vec![region([200.0, 400.0], a), region([400.0, 800.0], b)]);
        source
            .evaluate_heat_capacity(Temperature::new(400.0))
            .unwrap();
        assert!(matches!(
            source.evaluate_heat_capacity_temperature_derivative(Temperature::new(400.0)),
            Err(Nasa9TemperatureDerivativeErrorV1::NonDifferentiable {
                lower_region: 0,
                upper_region: 1,
                ..
            })
        ));
        source
            .evaluate_heat_capacity_temperature_derivative(Temperature::new(399.0))
            .unwrap();
        source
            .evaluate_heat_capacity_temperature_derivative(Temperature::new(401.0))
            .unwrap();
    }
    let gap = model(vec![region([200.0, 400.0], a), region([500.0, 800.0], a)]);
    for t in [0.0, f64::NAN, 199.0, 450.0, 801.0] {
        assert!(matches!(
            gap.evaluate_heat_capacity_temperature_derivative(Temperature::new(t)),
            Err(Nasa9TemperatureDerivativeErrorV1::Source(_))
        ));
    }
}
