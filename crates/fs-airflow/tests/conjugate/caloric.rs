//! Original source Cp(T), analytic channel heat, and native coupled admission.

#[path = "caloric/linearization.rs"]
mod linearization;

use std::sync::Arc;

use fs_airflow::conjugate::caloric::CaloricAirPath;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};
use fs_qty::{Pressure, QuantitySpec};

use super::*;

fn source(knots: &[(f64, f64)], datum: f64) -> Arc<IsobaricHeatCapacityLaw> {
    let key = PropertyKey::new("Cp", SPECIFIC_HEAT_CAPACITY_DIMS);
    let low = knots[0].0;
    let high = knots.last().unwrap().0;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: knots.to_vec(),
                dims: SPECIFIC_HEAT_CAPACITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), low, high)
                .with_quantity(
                    "p",
                    QuantitySpec::dimensional(Pressure::DIMS),
                    100_000.0,
                    100_000.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "analytic Cp channel fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured coolant".into(),
            phase: "single-phase".into(),
            process: "fixed pressure".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
            .unwrap()
            .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), 100_000.0)
            .unwrap()
    };
    Arc::new(
        IsobaricHeatCapacityLaw::try_new(
            &card,
            IsobaricHeatCapacityQuery {
                specific_heat_cp: &key,
                temperature_axis: "T",
                pressure_axis: "p",
                lower: &point(low),
                upper: &point(high),
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            },
            datum,
        )
        .unwrap(),
    )
}

fn path(law: Arc<IsobaricHeatCapacityLaw>, inlet: f64, count: usize) -> CaloricAirPath {
    let segments = (0..count)
        .map(|i| AirSegment::new(&format!("region-{i}"), 0.2 / count as f64, 200.0).unwrap())
        .collect();
    CaloricAirPath::new(inlet, 0.02, law, segments).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:e}"
    );
}

#[test]
fn original_constant_cp_recovers_native_channel_in_heating_and_cooling() {
    let law = source(&[(250.0, 1005.0), (650.0, 1005.0)], 300.0);
    for (inlet, wall) in [(300.0, 450.0), (450.0, 300.0), (300.0, 300.0)] {
        let path = path(law.clone(), inlet, 1);
        let native = AirPath::new(inlet, 0.02, 1005.0, path.segments().to_vec()).unwrap();
        let result = path.march(&[wall]).unwrap();
        let expected = native.march(&[wall]).unwrap();
        near(
            result.outlet_temperature_k,
            expected.outlet_temperature_k,
            2e-11,
        );
        near(result.total_heat_rate_w, expected.total_heat_rate_w, 2e-9);
        near(
            result.segments[0].reference_temperature_k,
            expected.segments[0].reference_temperature_k,
            2e-11,
        );
        assert_eq!(path.caloric_law().identity(), law.identity());
    }
}

#[test]
fn weak_exchange_retains_heat_and_robin_mean_before_rounding_absolute_temperature() {
    for (intercept, slope) in [(1000.0, 0.0), (800.0, 2.0)] {
        let law = source(
            &[
                (250.0, intercept + slope * 250.0),
                (650.0, intercept + slope * 650.0),
            ],
            300.0,
        );
        for (inlet, wall) in [(300.0, 400.0), (550.0, 280.0)] {
            let cp = intercept + slope * inlet;
            let driving = wall - inlet;
            let curvature = slope * driving / cp;
            for ntu in [1e-9, 1e-12] {
                let conductance = cp * ntu;
                let path = CaloricAirPath::new(
                    inlet,
                    1.0,
                    law.clone(),
                    vec![AirSegment::new("weak wall", 1.0, conductance).unwrap()],
                )
                .unwrap();
                let result = path.march(&[wall]).unwrap();
                let state = &result.segments[0];

                // Independent local solution of dz/dn=(1-z)/(1+b*z), where
                // z=(T-Tin)/(Tw-Tin), n=UA/(mass_flow*Cp(Tin)), and
                // b=Cp'(Tin)*(Tw-Tin)/Cp(Tin). The omitted relative heat
                // term is O(n^3), far below f64 rounding at these exposures.
                let expected_heat =
                    conductance * driving * (1.0 - 0.5 * ntu + (1.0 + curvature) * ntu * ntu / 6.0);
                let expected_mean = inlet
                    + driving
                        * ntu
                        * (0.5 - (1.0 + curvature) * ntu / 6.0
                            + (1.0 + curvature) * (1.0 + 3.0 * curvature) * ntu * ntu / 24.0);
                let expected_outlet = inlet
                    + driving
                        * ntu
                        * (1.0 - 0.5 * (1.0 + curvature) * ntu
                            + (1.0 + curvature) * (1.0 + 3.0 * curvature) * ntu * ntu / 6.0);
                let temperature_rounding = 8.0 * f64::EPSILON * inlet.max(wall);
                let heat_rounding = 2e-12 * expected_heat.abs();
                near(state.heat_rate_w, expected_heat, heat_rounding);
                near(result.total_heat_rate_w, expected_heat, heat_rounding);
                near(
                    state.reference_temperature_k,
                    expected_mean,
                    temperature_rounding,
                );
                near(
                    state.outlet_temperature_k,
                    expected_outlet,
                    temperature_rounding,
                );
                near(
                    conductance * (wall - state.reference_temperature_k),
                    state.heat_rate_w,
                    heat_rounding + conductance * temperature_rounding,
                );
                assert!(state.reference_temperature_k >= inlet.min(state.outlet_temperature_k));
                assert!(state.reference_temperature_k <= inlet.max(state.outlet_temperature_k));
            }
        }
    }
}

#[test]
fn affine_cp_matches_independent_channel_integral_and_conserves_heat() {
    // Cp(T) = 800 + 2T. The exact exchange primitive is
    // (800 + 2Tw) ln(|Tw-Tin|/|Tw-Tout|) - 2(Tout-Tin).
    let law = source(&[(250.0, 1300.0), (400.0, 1600.0), (650.0, 2100.0)], 300.0);
    for (inlet, wall) in [(300.0, 500.0), (550.0, 280.0)] {
        let original = path(law.clone(), inlet, 1);
        let result = original.march(&[wall]).unwrap();
        let t = result.outlet_temperature_k;
        let exact_exposure = (800.0 + 2.0 * wall) * ((wall - inlet).abs() / (wall - t).abs()).ln()
            - 2.0 * (t - inlet);
        near(exact_exposure, 40.0 / 0.02, 2e-9);
        let heat = 0.02 * ((800.0 + 2.0 * inlet) * (t - inlet) + (t - inlet).powi(2));
        near(result.total_heat_rate_w, heat, 2e-10);
        near(
            40.0 * (wall - result.segments[0].reference_temperature_k),
            heat,
            2e-10,
        );
        assert!(t >= inlet.min(wall) && t <= inlet.max(wall));
        let frozen = AirPath::new(
            inlet,
            0.02,
            800.0 + 2.0 * inlet,
            original.segments().to_vec(),
        )
        .unwrap()
        .march(&[wall])
        .unwrap();
        assert!(
            (t - frozen.outlet_temperature_k).abs() > 1.0,
            "the nonlinear source must affect the physical outlet"
        );
    }
}

#[test]
fn original_knots_refinement_and_caloric_datum_preserve_the_same_channel() {
    let knots = [
        (250.0, 1200.0),
        (330.0, 1500.0),
        (410.0, 1100.0),
        (650.0, 2200.0),
    ];
    for (inlet, wall) in [(290.0, 620.0), (610.0, 260.0)] {
        let (one, one_response) = path(source(&knots, 300.0), inlet, 1)
            .march_with_inlet_response(&[wall])
            .unwrap();
        for count in [2, 9, 32] {
            let (split, split_response) = path(source(&knots, 600.0), inlet, count)
                .march_with_inlet_response(&vec![wall; count])
                .unwrap();
            near(split.outlet_temperature_k, one.outlet_temperature_k, 2e-9);
            near(split.total_heat_rate_w, one.total_heat_rate_w, 2e-7);
            near(
                split_response.enthalpy_gain,
                one_response.enthalpy_gain,
                2e-12,
            );
            near(
                split_response.temperature_gain,
                one_response.temperature_gain,
                2e-12,
            );
        }
    }
}

#[test]
fn analytic_inlet_response_recovers_constant_cp_in_heating_cooling_and_equilibrium() {
    let cp = 1005.0;
    let law = source(&[(250.0, cp), (650.0, cp)], 300.0);
    let expected = (-40.0_f64 / (0.02 * cp)).exp();
    for (inlet, wall) in [(300.0, 450.0), (450.0, 300.0), (300.0, 300.0)] {
        for count in [1, 4] {
            let original = path(law.clone(), inlet, count);
            let walls = vec![wall; count];
            let (march, response) = original.march_with_inlet_response(&walls).unwrap();
            assert_eq!(march, original.march(&walls).unwrap());
            near(response.enthalpy_gain, expected, 2e-14);
            near(response.temperature_gain, expected, 2e-14);
        }
    }
}

#[test]
fn affine_cp_inlet_response_matches_analytic_transfer_and_independent_differences() {
    // An independently prescribed outlet gives the complete transfer exposure
    // and its analytic endpoint derivatives, without reading the solver result.
    let knots = [(250.0, 1300.0), (400.0, 1600.0), (650.0, 2100.0)];
    for (inlet, wall, outlet) in [(300.0_f64, 500.0_f64, 425.0_f64), (550.0, 280.0, 370.0)] {
        let exposure =
            (800.0 + 2.0 * wall) * ((wall - inlet) / (wall - outlet)).ln() - 2.0 * (outlet - inlet);
        let segments = vec![AirSegment::new("original source", 1.0, 0.02 * exposure).unwrap()];
        let law = source(&knots, 300.0);
        let original = CaloricAirPath::new(inlet, 0.02, law.clone(), segments.clone()).unwrap();
        let (march, response) = original.march_with_inlet_response(&[wall]).unwrap();
        let expected_enthalpy_gain = (wall - outlet) / (wall - inlet);
        let expected_temperature_gain =
            expected_enthalpy_gain * (800.0 + 2.0 * inlet) / (800.0 + 2.0 * outlet);
        near(march.outlet_temperature_k, outlet, 2e-11);
        near(response.enthalpy_gain, expected_enthalpy_gain, 2e-12);
        near(response.temperature_gain, expected_temperature_gain, 2e-12);

        let shift = 1e-3;
        let outputs = [-shift, shift].map(|offset| {
            CaloricAirPath::new(inlet + offset, 0.02, law.clone(), segments.clone())
                .unwrap()
                .march(&[wall])
                .unwrap()
                .outlet_temperature_k
        });
        near(
            response.temperature_gain,
            (outputs[1] - outputs[0]) / (2.0 * shift),
            2e-8,
        );
        let output_heat = law
            .specific_enthalpy_change_j_kg(outputs[0], outputs[1])
            .unwrap();
        let inlet_heat = law
            .specific_enthalpy_change_j_kg(inlet - shift, inlet + shift)
            .unwrap();
        near(response.enthalpy_gain, output_heat / inlet_heat, 2e-8);

        let changed_datum = CaloricAirPath::new(inlet, 0.02, source(&knots, 600.0), segments)
            .unwrap()
            .march_with_inlet_response(&[wall])
            .unwrap();
        assert_eq!((march, response), changed_datum);
    }
}

#[test]
fn inlet_response_retains_tiny_gain_when_absolute_outlet_rounds_to_wall() {
    let law = source(&[(250.0, 1300.0), (650.0, 2100.0)], 300.0);
    for (inlet, wall) in [(300.0, 500.0), (550.0, 280.0)] {
        // Cp(T)=800+2T; selecting log attenuation 50 leaves an unrepresentable
        // absolute temperature gap, but both derivatives remain representable.
        let cp_wall = 800.0 + 2.0 * wall;
        let exposure = cp_wall * 50.0 - 2.0 * (wall - inlet);
        let original = CaloricAirPath::new(
            inlet,
            1.0,
            law.clone(),
            vec![AirSegment::new("strong exchange", 1.0, exposure).unwrap()],
        )
        .unwrap();
        let (march, response) = original.march_with_inlet_response(&[wall]).unwrap();
        assert_eq!(march.outlet_temperature_k, wall);
        let expected = (-50.0_f64).exp();
        near(response.enthalpy_gain, expected, 2e-12 * expected);
        near(
            response.temperature_gain,
            expected * (800.0 + 2.0 * inlet) / cp_wall,
            2e-12 * expected,
        );
        assert!(response.enthalpy_gain > 0.0);
    }
}

#[test]
fn inlet_response_preserves_equilibrium_limits_and_native_resolution_refusal() {
    let inlet: f64 = 400.0;
    let law = source(&[(250.0, 1300.0), (650.0, 2100.0)], 300.0);
    let cp = law.specific_heat_cp_j_kg_k(inlet).unwrap();
    for wall in [inlet.next_down(), inlet, inlet.next_up()] {
        let original = CaloricAirPath::new(
            inlet,
            1.0,
            law.clone(),
            vec![AirSegment::new("near equilibrium", 1.0, cp).unwrap()],
        )
        .unwrap();
        let (_, response) = original.march_with_inlet_response(&[wall]).unwrap();
        near(response.enthalpy_gain, (-1.0_f64).exp(), 2e-14);
        near(response.temperature_gain, (-1.0_f64).exp(), 2e-14);
    }
    let unresolved = CaloricAirPath::new(
        inlet,
        1.0,
        law,
        vec![AirSegment::new("unresolved", 1.0, cp * 1e-20).unwrap()],
    )
    .unwrap();
    assert!(matches!(
        unresolved.march(&[500.0]),
        Err(AirflowError::CaloricSource(_))
    ));
    assert!(matches!(
        unresolved.march_with_inlet_response(&[500.0]),
        Err(AirflowError::CaloricSource(_))
    ));
}

#[test]
fn decreasing_cp_can_amplify_temperature_response_without_wall_extrapolation() {
    // Cp(T)=6400-18T on the actual fluid interval. Its extrapolation to the
    // hotter wall would be negative and must never be used as caloric data.
    let (inlet, wall, outlet) = (300.0_f64, 1000.0_f64, 340.0_f64);
    let law = source(&[(290.0, 1180.0), (350.0, 100.0)], inlet);
    assert!(law.specific_heat_cp_j_kg_k(wall).is_err());
    let exposure =
        (6400.0 - 18.0 * wall) * ((wall - inlet) / (wall - outlet)).ln() + 18.0 * (outlet - inlet);
    let original = CaloricAirPath::new(
        inlet,
        1.0,
        law,
        vec![AirSegment::new("external wall", 1.0, exposure).unwrap()],
    )
    .unwrap();
    let (march, response) = original.march_with_inlet_response(&[wall]).unwrap();
    let gain = (wall - outlet) / (wall - inlet);
    near(march.outlet_temperature_k, outlet, 3e-11);
    near(response.enthalpy_gain, gain, 2e-12);
    near(response.temperature_gain, gain * 1000.0 / 280.0, 2e-11);
    assert!(response.enthalpy_gain < 1.0);
    assert!(response.temperature_gain > 1.0);
}

#[test]
fn actual_fluid_support_is_enforced_without_requiring_wall_calorimetry() {
    let law = source(&[(290.0, 1000.0), (350.0, 1000.0)], 300.0);
    let supported = CaloricAirPath::new(
        300.0,
        1.0,
        law.clone(),
        vec![AirSegment::new("wall", 0.001, 1.0).unwrap()],
    )
    .unwrap();
    assert!(supported.march(&[1000.0]).unwrap().outlet_temperature_k < 350.0);
    let crossing = path(law, 300.0, 1);
    assert!(matches!(
        crossing.march(&[1000.0]),
        Err(AirflowError::CaloricSource(_))
    ));
    for wall in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(supported.march(&[wall]).is_err());
    }
    assert!(supported.march(&[]).is_err());
}

#[test]
fn native_coupled_solver_delivers_power_into_the_original_variable_cp_stream() {
    let law = source(&[(250.0, 1300.0), (650.0, 2100.0)], 300.0);
    let path = path(law, 300.0, 1);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1e-10,
        balance_tolerance_w: 1e-7,
        max_iterations: 200,
        ..ConjugateConfig::default()
    };
    let power = 800.0;
    with_cx(|cx| {
        let result = path
            .solve(cx, &config, |_, references| {
                Ok(vec![SolidRegionState {
                    region: "region-0".into(),
                    area_m2: 0.2,
                    mean_wall_temperature_k: references[0] + power / 40.0,
                    heat_rate_w: power,
                    mean_reference_temperature_k: Some(references[0]),
                }])
            })
            .unwrap();
        // Global calorimetry: 1400*dT + dT^2 = power / mass_flow.
        let rise =
            2.0 * (power / 0.02) / (1400.0 + (1400.0_f64.powi(2) + 4.0 * power / 0.02).sqrt());
        near(result.march.outlet_temperature_k, 300.0 + rise, 1e-8);
        near(result.march.total_heat_rate_w, power, 1e-7);
        assert!(result.balance.max_region_imbalance_w < 1e-7);
        let repeated = path
            .march(&[result.solid[0].mean_wall_temperature_k])
            .unwrap();
        assert_eq!(repeated, result.march);
        let accelerated = path
            .solve_iqn(
                cx,
                &config,
                fs_airflow::conjugate::IqnIlsConfig::default(),
                |_, references| {
                    Ok(vec![SolidRegionState {
                        region: "region-0".into(),
                        area_m2: 0.2,
                        mean_wall_temperature_k: references[0] + power / 40.0,
                        heat_rate_w: power,
                        mean_reference_temperature_k: Some(references[0]),
                    }])
                },
            )
            .unwrap();
        near(accelerated.march.outlet_temperature_k, 300.0 + rise, 1e-8);
        near(accelerated.march.total_heat_rate_w, power, 1e-7);
    });
}

#[test]
fn caloric_coupling_keeps_native_cancellation_before_any_solid_work() {
    let path = path(source(&[(250.0, 1300.0), (650.0, 2100.0)], 300.0), 300.0, 1);
    let gate = CancelGate::new();
    gate.request();
    with_gated_cx(&gate, |cx| {
        let result = path.solve(cx, &ConjugateConfig::default(), |_, _| {
            panic!("cancelled solve must not call the solid")
        });
        assert!(matches!(result, Err(AirflowError::Cancelled { .. })));
    });
}
