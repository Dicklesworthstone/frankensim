//! The native reverse follows spatial emission, full reflection and original heat.
use super::*;
use super::adjoint::{base, evaluate, maximum, gradient, first, derivative_close,
    radiation_surface, emissivity_difference};

const FINAL: &str = "temperature-final-adjoint";
const PEAK: &str = "temperature-sampled-peak-adjoint";

fn pointwise() -> (ProjectSpec, CardPackSet) {
    let (mut spec, cards) = base();
    let declaration = spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .enclosure_radiation.as_mut().unwrap();
    declaration.spatial_model = Some(fs_project::RadiationSpatialModel::Pointwise);
    declaration.max_iterations = 1;
    (spec, cards)
}

fn check_pointwise(receipt: &J) {
    let g = gradient(receipt);
    assert!(g.get("radiation").unwrap().str_field("method").unwrap().contains("pointwise"));
    close(n(g.get("radiation").unwrap(), "control_steps"), n(g, "state_index"), 0.0);
    assert!(n(g, "maximum_dual_relative_residual") < 1e-9);
    assert!(n(g, "adjoint_krylov_iterations") <= n(g, "max_adjoint_krylov_iterations"));
    assert_eq!(receipt.get("radiation").unwrap().str_field("spatial_model"), Some("pointwise"));
}

#[test]
fn g1_pointwise_vacuum_enclosure_reverse_preserves_physical_history_and_reflected_controls() {
    let (mut spec, cards) = pointwise();
    for (index, boundary) in spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .boundaries.iter_mut().enumerate()
    {
        boundary.condition = ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(if index == 0 { -12.0 } else { 2.0 }, dims::HEAT_FLUX),
        };
    }
    let (forward, expected_field) = evaluate(&spec, &cards, None);
    let (receipt, field) = evaluate(&spec, &cards, Some(FINAL));
    assert_eq!(field, expected_field);
    for key in ["energy", "radiation"] { assert_eq!(receipt.get(key), forward.get(key)); }
    check_pointwise(&receipt);
    let g = gradient(&receipt);
    close(n(g, "active_ties"), 1.0, 0.0);
    assert!(g.get("prescribed_convection").is_none());
    close(n(receipt.get("energy").unwrap(), "neumann_out_w"), -5.0, 1e-12);
    close(n(receipt.get("transient").unwrap().get("energy").unwrap(), "stored_change_j"),
        14.0, 1e-5);
    let patches = receipt.get("radiation").unwrap().get("surfaces").unwrap().as_array().unwrap();
    assert!(n(&patches[0], "outward_heat_w").abs() > 1e-5);
    assert!(n(&patches[0], "maximum_temperature_k") - n(&patches[0], "minimum_temperature_k") > 1e-4);
    for (index, epsilon) in [(0, EPSILON), (1, 0.6)] {
        let derivative = n(radiation_surface(g, index), "dtemperature_demissivity_k");
        assert!(derivative.abs() > 1e-8);
        derivative_close(derivative, emissivity_difference(&spec, &cards, index, epsilon));
    }
    let initial = n(g, "dtemperature_duniform_initial_specific_enthalpy_k_kg_per_j");
    let density = n(first(g, "regional_storage"), "dtemperature_dreference_density_k_m3_per_kg");
    let first_pulse: f64 = g.get("intervals").unwrap().as_array().unwrap().iter()
        .filter(|row| n(row, "end_s") <= 0.3)
        .map(|row| n(first(row, "regional_power"), "dtemperature_ddelivered_watts_k_per_w")).sum();
    let vertex = n(g, "active_vertex") as usize;
    let probe = |delta: f64| {
        let mut changed = spec.clone();
        let policy = time(&mut changed);
        let material = &mut policy.enthalpy.as_mut().unwrap().materials[0];
        material.initial_specific_enthalpy.value += 2.0 * delta;
        material.reference_density.value += 0.4 * delta;
        policy.power_schedules[0].steps[0].watts.value += 3.0 * delta;
        evaluate(&changed, &cards, None).1[vertex]
    };
    derivative_close(2.0 * initial + 0.4 * density + 3.0 * first_pulse,
        (probe(0.02) - probe(-0.02)) / 0.04);
}

#[test]
fn g1_pointwise_enclosure_peak_differentiates_only_the_selected_physical_prefix() {
    let (mut spec, cards) = pointwise();
    for (index, step) in time(&mut spec).power_schedules[0].steps.iter_mut().enumerate() {
        step.watts.value = if index == 0 { 200.0 } else { 0.0 };
    }
    for boundary in &mut spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
        if let ThermalBoundaryCondition::Convection { reference_temperature, .. } = &mut boundary.condition {
            reference_temperature.value = 290.0;
        }
    }
    let (receipt, final_field) = evaluate(&spec, &cards, Some(PEAK));
    check_pointwise(&receipt);
    let g = gradient(&receipt);
    close(n(g, "time_s"), 0.3, 0.0);
    close(n(g, "state_index"), 2.0, 0.0);
    close(n(g, "revalidated_solid_endpoints"), 2.0, 0.0);
    assert!(n(g, "value_k") > maximum(&final_field));
    let peak = receipt.get("transient").unwrap().get("sampled_peak").unwrap();
    close(n(g, "value_k"), n(peak, "value_k"), 0.0);
    for interval in g.get("intervals").unwrap().as_array().unwrap() {
        if n(interval, "start_s") >= 0.3 {
            close(n(first(interval, "regional_power"), "dtemperature_ddelivered_watts_k_per_w"),
                0.0, 0.0);
        }
    }
    let first_pulse: f64 = g.get("intervals").unwrap().as_array().unwrap().iter()
        .filter(|row| n(row, "end_s") <= 0.3)
        .map(|row| n(first(row, "regional_power"), "dtemperature_ddelivered_watts_k_per_w")).sum();
    let probe = |delta: f64| {
        let mut changed = spec.clone();
        time(&mut changed).power_schedules[0].steps[0].watts.value += delta;
        let receipt = evaluate(&changed, &cards, None).0;
        n(receipt.get("transient").unwrap().get("sampled_peak").unwrap(), "value_k")
    };
    derivative_close(first_pulse, (probe(0.02) - probe(-0.02)) / 0.04);
}

#[test]
fn g3_pointwise_enclosure_reverse_keeps_streamwise_air_and_natural_reservoir_feedback() {
    let (mut spec, cards) = pointwise();
    spec.envelope.as_mut().unwrap().ambient_lo.value = 299.0;
    spec.envelope.as_mut().unwrap().ambient_hi.value = 351.0;
    time(&mut spec).enthalpy.as_mut().unwrap().materials[0]
        .initial_specific_enthalpy.value = CP * (330.0 - INITIAL_K);
    let prototype = conjugate_fixture_project(7, &tetra_stl(), "convection.gnielinski");
    let mut air = prototype.cooling.as_ref().unwrap().conduction.as_ref().unwrap()
        .boundaries[0].condition.clone();
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
    let (forward, expected) = evaluate(&spec, &cards, None);
    let (receipt, actual) = evaluate(&spec, &cards, Some(FINAL));
    assert_eq!(actual, expected);
    assert_eq!(receipt.get("energy"), forward.get("energy"));
    check_pointwise(&receipt);
    let g = gradient(&receipt);
    let vertex = n(g, "active_vertex") as usize;
    let controls = g.get("convection").unwrap();
    close(n(controls, "natural_ports"), 1.0, 0.0);
    close(n(controls, "air_ports"), 2.0, 0.0);
    for inlet in [true, false] {
        let probe = |delta: f64| {
            let mut changed = spec.clone();
            for row in &mut changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries {
                match &mut row.condition {
                    ThermalBoundaryCondition::AirflowConvection { inlet_temperature, .. } if inlet =>
                        inlet_temperature.value += delta,
                    ThermalBoundaryCondition::NaturalConvection { ambient_temperature, .. } if !inlet =>
                        ambient_temperature.value += delta,
                    _ => {},
                }
            }
            evaluate(&changed, &cards, None).1[vertex]
        };
        let expected = if inlet {
            n(first(controls, "air_inlets"), "dtemperature_dinlet_temperature_k_per_k")
        } else {
            n(first(controls, "natural_reservoirs"), "dtemperature_dambient_temperature_k_per_k")
        };
        derivative_close(expected, (probe(0.02) - probe(-0.02)) / 0.04);
    }
}
