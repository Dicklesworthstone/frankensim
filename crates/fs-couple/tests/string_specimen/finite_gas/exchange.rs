//! Synthetic constant-capacity G1 and source-bound reciprocal G3/G4 checks.
use super::*;
use fs_couple::cavity_phs::piston::thermal::{
    ConductanceSample, HeatExchangeSolve, layers::MaterialConductionLayers,
};

fn exchange_solve() -> HeatExchangeSolve {
    HeatExchangeSolve {
        heat_bracket_j: [-8e-5, 8e-5],
        source_tolerance_j: 1e-13,
        energy_tolerance_j: 2e-10,
        maximum_iterations: 64,
    }
}

fn warm_gas(temperature_k: f64) -> IsochoricGasParcel {
    let mut gas = parcel(0.028, 100000.0);
    let cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
    gas.add_heat(
        gas.mass_kg() * cv * (temperature_k - 300.0),
        solve().gas,
        || false,
    )
    .unwrap();
    gas
}

// Reuse the existing synthetic mechanical/enthalpy fixture and append a
// source-owned conductivity curve. Tiny cp is intentional test data.
pub(super) fn quiet_body(
    cx: &fs_exec::Cx<'_>,
    gas: &fs_material::gas::GasState,
    temperature_k: f64,
) -> (
    ThermalMaterialStringRuntime,
    MaterialCard,
    QueryPoint,
    PropertyKey,
) {
    let (old, original, curve) =
        heated_string_parts(cx, 48000, false, InterpolationPolicy::LinearInside);
    let key = PropertyKey::new(
        "thermal_conductivity",
        fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS,
    );
    let mut claims = original.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: key.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: fs_qty::Temperature::DIMS,
                knots: vec![(300.0, 20.0), (340.0, 40.0)],
                dims: fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 340.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic reciprocal string conductor".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic reciprocal string".into(),
            phase: "solid".into(),
            process: "synthetic".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with("T", temperature_k)
        .unwrap()
        .with("omega", 1.0)
        .unwrap();
    let mut requirements: Vec<_> = old
        .specimen()
        .material()
        .properties()
        .iter()
        .map(|p| p.requirement().clone())
        .collect();
    requirements.push(
        ScalarPropertyRequirement::try_with_key(&key, ScalarAdmissibility::NonNegative).unwrap(),
    );
    let state = resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let bound = old.specimen();
    let specimen = with_uniform_circular_material_and_constraints(
        bound.string(),
        bound.geometry_constraint(),
        &state,
        bound.prestress(),
    )
    .unwrap()
    .with_kelvin_voigt_bending_loss()
    .unwrap();
    let mechanical =
        LinearMaterialStringRuntime::try_new_with_gas(cx, specimen, None, gas, 1.0, 48000).unwrap();
    let curve = fs_material::phase::EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        curve.knots().to_vec(),
    )
    .unwrap();
    let live = ThermalMaterialStringRuntime::try_new(
        cx,
        mechanical,
        card.clone(),
        curve,
        0.001 * (temperature_k - 300.0),
        100.0,
    )
    .unwrap();
    (live, card, point, key)
}

#[test]
fn g1_finite_gas_exchange_matches_two_capacity_implicit_solution_in_both_directions() {
    with_string_cx(|cx, _| {
        for (solid_k, gas_k, solid_source_k, gas_source_k) in [
            (300.0, 330.0, 2.0, 1.0),
            // Source heating reverses the initial temperature ordering.
            (320.0, 330.0, 10.0, -15.0),
        ] {
            let mut gas = warm_gas(gas_k);
            let (mut live, _, _, _) = quiet_body(cx, gas.state().state(), solid_k);
            let initial = live.clone();
            let gas_before = gas.clone();
            let cs = live.mechanical().specimen().mass_kg() * 0.001;
            let cv = 2.5 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028;
            let cg = gas.mass_kg() * cv;
            let g = 0.1;
            let dt_g = g / 48000.0;
            let heat = FiniteGasStringHeat {
                string_heat_j: cs * solid_source_k,
                gas_heat_j: cg * gas_source_k,
            };
            let source_s = live.thermal().temperature_k() + heat.string_heat_j / cs;
            let source_g = gas.state().state().temperature + heat.gas_heat_j / cg;
            let expected_q = dt_g * (source_s - source_g) / (1.0 + dt_g * (1.0 / cs + 1.0 / cg));
            let frame = live
                .step_with_isochoric_gas_exchange(
                    cx,
                    0,
                    &[0.0],
                    &mut gas,
                    heat,
                    solve(),
                    exchange_solve(),
                    |_, _| {
                        Ok::<_, &'static str>(ConductanceSample {
                            conductance_w_k: g,
                            source: "synthetic declared constant G",
                        })
                    },
                )
                .unwrap();
            let gas_tolerance = solve().gas.energy_tolerance_j_kg / cv;
            let q_tolerance = exchange_solve().source_tolerance_j + dt_g * gas_tolerance;
            assert!((frame.solid_to_gas_heat_j - expected_q).abs() <= 2.0 * q_tolerance);
            assert!(
                (live.thermal().temperature_k() - (source_s - expected_q / cs)).abs()
                    <= 2.0 * q_tolerance / cs + 1e-11
            );
            assert!(
                (gas.state().state().temperature - (source_g + expected_q / cg)).abs()
                    <= 2.0 * gas_tolerance + 2.0 * q_tolerance / cg + 1e-11
            );
            assert!(
                frame.solid_to_gas_heat_j
                    * (live.thermal().temperature_k() - gas.state().state().temperature)
                    >= 0.0
            );
            assert!(live.thermal().temperature_k() >= source_s.min(source_g));
            assert!(live.thermal().temperature_k() <= source_s.max(source_g));
            assert!(gas.state().state().temperature >= source_s.min(source_g) - gas_tolerance);
            assert!(gas.state().state().temperature <= source_s.max(source_g) + gas_tolerance);
            assert_eq!(frame.external_string_heat_j, heat.string_heat_j);
            assert_eq!(frame.coupled.external_gas_heat_j, heat.gas_heat_j);
            assert_eq!(frame.coupled.drag_heat_j, 0.0);
            assert_eq!(
                frame.coupled.string.external_heat_j,
                heat.string_heat_j - frame.solid_to_gas_heat_j
            );
            assert_eq!(
                frame.coupled.gas.supplied_heat_j,
                heat.gas_heat_j + frame.solid_to_gas_heat_j
            );
            let solid_storage = live.mechanical().specimen().mass_kg()
                * (live.thermal().specific_enthalpy_j_kg()
                    - initial.thermal().specific_enthalpy_j_kg());
            let gas_storage = gas.mass_kg()
                * (gas.state().specific_internal_energy_j_kg()
                    - gas_before.state().specific_internal_energy_j_kg());
            assert!(
                (solid_storage + gas_storage - heat.string_heat_j - heat.gas_heat_j).abs()
                    <= exchange_solve().energy_tolerance_j
            );
            assert!(frame.source_residual_j.abs() <= exchange_solve().source_tolerance_j);
            assert!(frame.thermal_energy_residual_j.abs() <= exchange_solve().energy_tolerance_j);
            assert!(
                frame.coupled.energy_balance_residual_j.abs() <= frame.coupled.energy_tolerance_j
            );
            assert_eq!(live.mechanical().ambient(), gas.state().state());
            assert_eq!(live.mechanical().epoch(), 1);
            assert_eq!(live.mechanical().accepted_samples(), 1);
        }
    });
}

#[test]
fn g3_finite_gas_conduction_uses_endpoint_sources_and_actual_geometry() {
    with_string_cx(|cx, _| {
        let gas = warm_gas(330.0);
        let (original, card, point, key) = quiet_body(cx, gas.state().state(), 300.0);
        let specimen = original.mechanical().specimen();
        let radius = specimen.radius_m();
        let area = 2.0 * core::f64::consts::PI * radius * (specimen.string().length_m + radius);
        let mut layers = MaterialConductionLayers {
            material: &card,
            context: &point,
            temperature_axis: "T",
            liquid_fraction_axis: None,
            conductivity: &key,
            selection: MaterialPropertySelection::SingleClaimOnly,
            area_m2: area,
            gas_distance_m: 0.001,
            body_distance_m: radius / 2.0,
        };
        let mut transfers = Vec::new();
        for distance in [0.001, 0.002] {
            layers.gas_distance_m = distance;
            let mut live = original.clone();
            let mut next_gas = gas.clone();
            let frame = live
                .step_with_isochoric_gas_conduction(
                    cx,
                    0,
                    &[0.0],
                    &mut next_gas,
                    FiniteGasStringHeat::default(),
                    solve(),
                    exchange_solve(),
                    &layers,
                )
                .unwrap();
            let kb = 20.0 + 0.5 * (live.thermal().temperature_k() - 300.0);
            let kg = next_gas.state().state().thermal_conductivity;
            close(
                frame.conductance_w_k,
                area / (distance / kg + layers.body_distance_m / kb),
            );
            close(frame.source.property_by_key(&key).unwrap().value_si(), kb);
            assert_eq!(frame.source.card_identity(), card.content_hash());
            assert_eq!(
                frame.source.query_point().iter().find(|(axis, _)| axis == "T").map(|(_, value)| value),
                Some(&live.thermal().temperature_k())
            );
            assert_eq!(
                frame.source.query_point().iter().find(|(axis, _)| axis == "omega").map(|(_, value)| value),
                Some(&1.0)
            );
            assert_eq!(live.mechanical().ambient(), next_gas.state().state());
            assert!(
                live.mechanical().modes()[0].damping_ratio
                    > original.mechanical().modes()[0].damping_ratio
            );
            assert_eq!(live.mechanical().states(), original.mechanical().states());
            transfers.push(-frame.solid_to_gas_heat_j);
        }
        assert!(transfers[0] > transfers[1] && transfers[1] > 0.0);
        for geometry in 0..3 {
            let mut live = original.clone();
            let mut next_gas = gas.clone();
            let wrong_point = point.clone().with("omega", 2.0).unwrap();
            let rejected = MaterialConductionLayers {
                material: &card,
                context: if geometry == 2 { &wrong_point } else { &point },
                temperature_axis: "T",
                liquid_fraction_axis: None,
                conductivity: &key,
                selection: MaterialPropertySelection::SingleClaimOnly,
                area_m2: if geometry == 0 { area * 1.01 } else { area },
                gas_distance_m: if geometry == 1 { 1e-9 } else { 0.001 },
                body_distance_m: if geometry == 1 { 1e-9 } else { radius / 2.0 },
            };
            assert!(
                live.step_with_isochoric_gas_conduction(
                    cx,
                    0,
                    &[0.0],
                    &mut next_gas,
                    FiniteGasStringHeat::default(),
                    solve(),
                    exchange_solve(),
                    &rejected,
                )
                .is_err()
            );
            assert_eq!(next_gas, gas);
            assert_eq!(live.thermal(), original.thermal());
            assert_eq!(live.mechanical().epoch(), 0);
        }
    });
}

#[test]
fn g3_zero_exchange_preserves_existing_finite_gas_motion_and_observer_history() {
    with_string_cx(|cx, _| {
        let mut gas = parcel(0.028, 100000.0);
        let mut live = runtime(cx, gas.state().state());
        let mut baseline = live.clone();
        let mut baseline_gas = gas.clone();
        for epoch in 0..3 {
            let frame = live
                .step_with_isochoric_gas_exchange(
                    cx,
                    epoch,
                    &[0.0; 2],
                    &mut gas,
                    FiniteGasStringHeat::default(),
                    solve(),
                    exchange_solve(),
                    |_, _| {
                        Ok::<_, &'static str>(ConductanceSample {
                            conductance_w_k: 0.0,
                            source: "explicit insulation",
                        })
                    },
                )
                .unwrap();
            let old = baseline
                .step_with_isochoric_gas(
                    cx,
                    epoch,
                    &[0.0; 2],
                    &mut baseline_gas,
                    FiniteGasStringHeat::default(),
                    solve(),
                )
                .unwrap();
            assert_eq!(frame.solid_to_gas_heat_j, 0.0);
            assert_eq!(frame.trials, 1);
            assert_eq!(frame.coupled.string.vibration, old.string.vibration);
            assert_eq!(gas, baseline_gas);
            assert_eq!(live.thermal(), baseline.thermal());
            assert_eq!(live.mechanical().states(), baseline.mechanical().states());
            assert_eq!(
                live.mechanical().specimen(),
                baseline.mechanical().specimen()
            );
        }
        assert_eq!(
            live.step(cx, 3, &[0.0; 2], 0.0).unwrap(),
            baseline.step(cx, 3, &[0.0; 2], 0.0).unwrap()
        );
    });
}

#[test]
fn g4_finite_gas_exchange_refusals_and_late_cancellation_are_atomic() {
    with_string_cx(|cx, gate| {
        let gas = warm_gas(330.0);
        let (pristine, _, _, _) = quiet_body(cx, gas.state().state(), 300.0);
        let cs = pristine.mechanical().specimen().mass_kg() * 0.001;
        for failure in 0..5 {
            let mut live = pristine.clone();
            let mut next_gas = gas.clone();
            let mut limits = exchange_solve();
            if failure == 0 {
                limits.heat_bracket_j = [-1e-8, 1e-8];
            }
            let heat = if failure == 1 {
                FiniteGasStringHeat {
                    string_heat_j: cs * 50.0,
                    gas_heat_j: 0.0,
                }
            } else if failure == 2 {
                FiniteGasStringHeat {
                    string_heat_j: 0.0,
                    gas_heat_j: 1000.0,
                }
            } else {
                FiniteGasStringHeat::default()
            };
            let result = live.step_with_isochoric_gas_exchange(
                cx,
                0,
                &[0.0],
                &mut next_gas,
                heat,
                solve(),
                limits,
                |_, _| {
                    if failure == 3 {
                        return Err("source unavailable");
                    }
                    Ok(ConductanceSample {
                        conductance_w_k: if failure == 1 {
                            0.0
                        } else if failure == 4 {
                            -1.0
                        } else {
                            0.1
                        },
                        source: "synthetic refusal probe",
                    })
                },
            );
            assert!(result.is_err(), "failure case {failure} was accepted");
            assert_eq!(next_gas, gas);
            assert_eq!(live.thermal(), pristine.thermal());
            assert_eq!(live.mechanical().states(), pristine.mechanical().states());
            assert_eq!(live.mechanical().ambient(), pristine.mechanical().ambient());
            assert_eq!(live.mechanical().epoch(), 0);
            let mut control = pristine.clone();
            assert_eq!(
                live.step(cx, 0, &[0.0], 0.0).unwrap(),
                control.step(cx, 0, &[0.0], 0.0).unwrap()
            );
        }
        let mut live = pristine.clone();
        let mut next_gas = gas.clone();
        let mut reached_source = false;
        let result = live.step_with_isochoric_gas_exchange(
            cx,
            0,
            &[0.0],
            &mut next_gas,
            FiniteGasStringHeat::default(),
            solve(),
            exchange_solve(),
            |_, _| {
                reached_source = true;
                gate.request();
                Ok::<_, &'static str>(ConductanceSample {
                    conductance_w_k: 0.1,
                    source: "cancel after thermal proposal",
                })
            },
        );
        assert!(reached_source && result.is_err());
        assert_eq!(next_gas, gas);
        assert_eq!(live.thermal(), pristine.thermal());
        assert_eq!(live.mechanical().states(), pristine.mechanical().states());
        assert_eq!(live.mechanical().epoch(), 0);
    });
}
