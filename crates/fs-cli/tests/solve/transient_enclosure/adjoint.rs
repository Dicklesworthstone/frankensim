//! Complete native enclosure history derivatives against fresh physical runs.
use super::*;
use fs_project::TransientInitialTemperature;

const FINAL: &str = "temperature-final-adjoint";
const PEAK: &str = "temperature-sampled-peak-adjoint";
const SECOND_EPSILON: f64 = 0.6;

fn base() -> (ProjectSpec, CardPackSet) {
    let (cards, first) = radiation_product::emissivity_cards(&fixture_cards(), EPSILON);
    let (cards, second) = radiation_product::emissivity_cards(&cards, SECOND_EPSILON);
    let mut spec = project(&first);
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .enclosure_radiation.as_mut().unwrap().surfaces[1].card = second;
    time(&mut spec).initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
    (spec, cards)
}

fn evaluate(spec: &ProjectSpec, cards: &CardPackSet, output: Option<&str>) -> (J, Vec<f64>) {
    let mut spec = spec.clone();
    if let Some(name) = output {
        spec.outputs.as_mut().unwrap().push(OutputRequest {
            name: name.into(), kind: "report".into(), region: None,
        });
    }
    let solved = run(&spec, cards);
    let receipt = J::parse(&solved.conduction).unwrap();
    let field = J::parse(&solved.field).unwrap();
    (receipt, values(&field, "temperature"))
}

fn maximum(values: &[f64]) -> f64 {
    values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
}

fn gradient(receipt: &J) -> &J {
    let row = receipt.get("transient").unwrap().get("adjoint").unwrap();
    assert_eq!(row.str_field("schema"), Some("fs-cli-native-enthalpy-adjoint-v1"));
    row
}

fn first<'a>(row: &'a J, key: &str) -> &'a J {
    &row.get(key).unwrap().as_array().unwrap()[0]
}

fn derivative_close(actual: f64, expected: f64) {
    close(actual, expected, 5e-6 + 3e-4 * expected.abs());
}

fn radiation_surface(row: &J, surface: usize) -> &J {
    &row.get("radiation").unwrap().get("surfaces").unwrap().as_array().unwrap()[surface]
}

fn emissivity_difference(
    spec: &ProjectSpec, cards: &CardPackSet, surface: usize, value: f64,
) -> f64 {
    let step = 0.0005;
    let trial = |delta: f64| {
        let (changed_cards, card) = radiation_product::emissivity_cards(cards, value + delta);
        let mut changed = spec.clone();
        changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
            .enclosure_radiation.as_mut().unwrap().surfaces[surface].card = card;
        maximum(&evaluate(&changed, &changed_cards, None).1)
    };
    (trial(step) - trial(-step)) / (2.0 * step)
}

#[test]
fn g1_native_enclosure_history_controls_match_independent_runs_without_changing_forward_state() {
    let (spec, cards) = base();
    let (forward, expected_field) = evaluate(&spec, &cards, None);
    let (receipt, actual_field) = evaluate(&spec, &cards, Some(FINAL));
    assert_eq!(actual_field, expected_field, "reverse must retain the accepted forward field");
    let g = gradient(&receipt);
    close(n(g, "value_k"), maximum(&actual_field), 0.0);
    close(n(g, "time_s"), 1.0, 0.0);
    close(n(g, "revalidated_solid_endpoints"), 6.0, 0.0);
    assert!(n(g, "maximum_dual_relative_residual") < 1e-10);
    assert!(n(g, "adjoint_krylov_iterations") <= n(g, "max_adjoint_krylov_iterations"));
    assert_eq!(g.get("radiation").unwrap().str_field("method"),
        Some("complete-closed-gray-diffuse-enthalpy-tangent"));
    for index in 0..2 {
        let row = radiation_surface(g, index);
        assert_eq!(row.str_field("target"), Some(if index == 0 { "wall-x" } else { "wall-z" }));
        assert_eq!(row.str_field("card"), radiation_surface(&forward, index).str_field("card"));
        close(n(row, "emissivity"), if index == 0 { EPSILON } else { SECOND_EPSILON }, 0.0);
        let derivative = n(row, "dtemperature_demissivity_k");
        assert!(derivative.abs() > 1e-8, "unequal patch temperatures must exercise reflection");
        derivative_close(derivative, emissivity_difference(
            &spec, &cards, index, if index == 0 { EPSILON } else { SECOND_EPSILON },
        ));
    }

    let initial = n(g, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j");
    let density = n(first(g, "regional_storage"), "dtemperature_dreference_density_k_m3_per_kg");
    let first_pulse: f64 = g.get("intervals").unwrap().as_array().unwrap().iter()
        .filter(|row| n(row, "end_s") <= 0.3)
        .map(|row| n(first(row, "regional_power"), "dtemperature_ddelivered_watts_k_per_w"))
        .sum();
    let trial = |delta: f64| {
        let mut changed = spec.clone();
        let policy = time(&mut changed);
        let material = &mut policy.enthalpy.as_mut().unwrap().materials[0];
        material.initial_specific_enthalpy.value += 2.0 * delta;
        material.reference_density.value += 0.4 * delta;
        policy.power_schedules[0].steps[0].watts.value += 3.0 * delta;
        maximum(&evaluate(&changed, &cards, None).1)
    };
    derivative_close(2.0 * initial + 0.4 * density + 3.0 * first_pulse,
        (trial(0.02) - trial(-0.02)) / 0.04);

    let robin = g.get("prescribed_convection").unwrap();
    close(n(robin, "control_steps"), 6.0, 0.0);
    let row = first(robin, "boundaries");
    assert_eq!(row.str_field("target"), Some("wall-x"));
    for coefficient in [false, true] {
        let trial = |delta: f64| {
            let mut changed = spec.clone();
            let ThermalBoundaryCondition::Convection {
                coefficient: h, reference_temperature,
            } = &mut changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                .boundaries[0].condition else { panic!("prescribed fixture"); };
            if coefficient { h.value += delta; } else { reference_temperature.value += delta; }
            maximum(&evaluate(&changed, &cards, None).1)
        };
        let key = if coefficient { "dtemperature_dcoefficient_k_per_w_m2_k" }
            else { "dtemperature_dreference_temperature_k_per_k" };
        derivative_close(n(row, key), (trial(0.01) - trial(-0.01)) / 0.02);
    }
}

#[test]
fn g1_native_enclosure_peak_reverse_keeps_initial_and_selected_prefix_semantics() {
    let (mut spec, cards) = base();
    time(&mut spec).enthalpy.as_mut().unwrap().materials[0]
        .initial_specific_enthalpy.value = 600.0;
    let (receipt, field) = evaluate(&spec, &cards, Some(PEAK));
    let g = gradient(&receipt);
    close(n(g, "state_index"), 0.0, 0.0);
    close(n(g, "value_k"), INITIAL_K + 600.0 / CP, 1e-12);
    assert!(maximum(&field) < n(g, "value_k"));
    close(n(g, "adjoint_krylov_iterations"), 0.0, 0.0);
    close(n(g, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j"), 1.0 / CP, 1e-14);
    close(n(g, "dtemperature_dreference_density_multiplier_k"), 0.0, 0.0);
    for surface in 0..2 { close(n(radiation_surface(g, surface), "dtemperature_demissivity_k"), 0.0, 0.0); }
    for row in g.get("prescribed_convection").unwrap().get("boundaries").unwrap().as_array().unwrap() {
        close(n(row, "dtemperature_dreference_temperature_k_per_k"), 0.0, 0.0);
        close(n(row, "dtemperature_dcoefficient_k_per_w_m2_k"), 0.0, 0.0);
    }

    time(&mut spec).enthalpy.as_mut().unwrap().materials[0].initial_specific_enthalpy.value = 0.0;
    for (index, step) in time(&mut spec).power_schedules[0].steps.iter_mut().enumerate() {
        step.watts.value = if index == 0 { 200.0 } else { 0.0 };
    }
    for boundary in &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
        if let ThermalBoundaryCondition::Convection { reference_temperature, .. } = &mut boundary.condition {
            reference_temperature.value = 290.0;
        }
    }
    let (receipt, _) = evaluate(&spec, &cards, Some(PEAK));
    let g = gradient(&receipt);
    close(n(g, "time_s"), 0.3, 0.0);
    close(n(g, "state_index"), 2.0, 0.0);
    close(n(g, "revalidated_solid_endpoints"), 2.0, 0.0);
    for interval in g.get("intervals").unwrap().as_array().unwrap() {
        if n(interval, "start_s") >= 0.3 {
            close(n(first(interval, "regional_power"), "dtemperature_ddelivered_watts_k_per_w"), 0.0, 0.0);
        }
    }
    let peak = receipt.get("transient").unwrap().get("sampled_peak").unwrap();
    close(n(g, "value_k"), n(peak, "value_k"), 0.0);
    close(n(g, "time_s"), n(peak, "time_s"), 0.0);
}

#[test]
fn g3_native_mixed_air_natural_enclosure_reverse_keeps_complete_boundary_feedback() {
    let (mut spec, cards) = base();
    spec.envelope.as_mut().unwrap().ambient_lo.value = 299.0;
    spec.envelope.as_mut().unwrap().ambient_hi.value = 351.0;
    time(&mut spec).enthalpy.as_mut().unwrap().materials[0].initial_specific_enthalpy.value = CP * (330.0 - INITIAL_K);
    let prototype = conjugate_fixture_project(7, &tetra_stl(), "convection.gnielinski");
    let mut air = prototype.cooling.as_ref().unwrap().conduction.as_ref().unwrap().boundaries[0].condition.clone();
    if let ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. } = &mut air {
        inlet_temperature.value = 300.0;
    }
    let setup = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.boundaries[0].condition = air.clone();
    setup.boundaries[1].condition = ThermalBoundaryCondition::NaturalConvection {
        characteristic_length: QtyAny::new(0.1, dims::LENGTH),
        ambient_temperature: QtyAny::new(350.0, dims::TEMPERATURE),
        correlation: "convection.churchill-chu-vertical-plate".into(),
    };
    if let ThermalBoundaryCondition::AirflowConvection { order, .. } = &mut air { *order = 1; }
    setup.boundaries.push(ThermalBoundary { target: "air".into(), condition: air });
    // Bound the nested fixture while retaining the same tight physical endpoint gates.
    setup.enclosure_radiation.as_mut().unwrap().max_iterations = 64;
    let (receipt, _) = evaluate(&spec, &cards, Some(FINAL));
    let g = gradient(&receipt);
    let convection = g.get("convection").unwrap();
    close(n(convection, "natural_ports"), 1.0, 0.0);
    close(n(convection, "air_ports"), 2.0, 0.0);
    assert!(n(g, "maximum_dual_relative_residual") < 1e-10);
    let first_air = first(convection, "air_inlets");
    let natural = first(convection, "natural_reservoirs");
    for inlet in [true, false] {
        let trial = |delta: f64| {
            let mut changed = spec.clone();
            for row in &mut changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
                match &mut row.condition {
                    ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. } if inlet => {
                        inlet_temperature.value += delta;
                    }
                    ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. } if !inlet => {
                        ambient_temperature.value += delta;
                    }
                    _ => {},
                }
            }
            maximum(&evaluate(&changed, &cards, None).1)
        };
        let expected = if inlet {
            n(first_air, "dtemperature_dinlet_temperature_k_per_k")
        } else {
            n(natural, "dtemperature_dambient_temperature_k_per_k")
        };
        derivative_close(expected, (trial(0.02) - trial(-0.02)) / 0.04);
    }
    derivative_close(n(radiation_surface(g, 0), "dtemperature_demissivity_k"),
        emissivity_difference(&spec, &cards, 0, EPSILON));
}
