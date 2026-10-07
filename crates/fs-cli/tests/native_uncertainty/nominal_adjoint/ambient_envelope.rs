//! Real native ambient-envelope re-solves must not replace radiation sources.
use super::*;

fn fixture_project(natural: bool, lo: f64, hi: f64) -> fs_project::ProjectSpec {
    let mut spec = project(natural);
    let cooling = spec.cooling.as_mut().unwrap();
    // No hydraulic uncertainty in this test of the fluid-temperature envelope.
    cooling.fans.clear();
    cooling.vents.clear();
    cooling.fan_system = None;
    cooling.airflow_leakage = None;
    let setup = cooling.conduction.as_mut().unwrap();
    setup.radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value =
        if natural { 350.0 } else { 280.0 };
    spec.envelope.as_mut().unwrap().ambient_lo.value = lo;
    spec.envelope.as_mut().unwrap().ambient_hi.value = hi;
    spec
}

fn maximum(fixture: &Fixture, receipt: &JsonValue) -> f64 {
    field(fixture, receipt).as_array().unwrap().iter()
        .map(|v| v.as_f64().unwrap()).reduce(f64::max).unwrap()
}

#[test]
fn exact_fluid_ambient_does_not_invent_radiative_boundary_uncertainty() {
    for natural in [false, true] {
        let fixture = Fixture::new();
        let spec = fixture_project(natural, 300.0, 300.0);
        let result = run(&fixture, &spec, 700);
        let propagation = result.get("propagation").unwrap();
        let nominal = propagation.f64_field("nominal_base_k").unwrap();
        let boundary = propagation.get("boundary_conditions").unwrap();
        assert_eq!(boundary.str_field("state"), Some("measured"));
        assert_eq!(boundary.f64_field("half_width_k"), Some(0.0),
            "a zero-width fluid envelope cannot replace an independently fixed reservoir");
        let vertices = boundary.get("vertices").unwrap().as_array().unwrap();
        assert_eq!(vertices.len(), 2);
        for vertex in vertices {
            assert_eq!(vertex.f64_field("t_max_k").unwrap().to_bits(), nominal.to_bits());
        }
        assert!(boundary.str_field("detail").unwrap().contains("radiative reservoirs remain fixed"));
    }
}

#[test]
fn ambient_budget_vertices_equal_isolated_physical_resolves_not_reservoir_substitutions() {
    for natural in [false, true] {
        let fixture = Fixture::new();
        let spec = fixture_project(natural, 295.0, 305.0);
        let nominal = run(&fixture, &spec, 710);
        let boundary = nominal.get("propagation").unwrap().get("boundary_conditions").unwrap();
        assert_eq!(boundary.str_field("state"), Some("measured"));
        let vertices = boundary.get("vertices").unwrap().as_array().unwrap();
        assert_eq!(vertices.len(), 2);
        let ambient = if natural { "natural-convection-ambient" } else { "convection-temperature" };
        for (i, temperature) in [295.0, 305.0].into_iter().enumerate() {
            let mut changed = spec.clone();
            perturb(&mut changed, ambient, temperature - 300.0);
            let actual = run(&fixture, &changed, 711 + 2 * i);
            let expected = maximum(&fixture, &actual);
            let retained = vertices[i].f64_field("t_max_k").unwrap();
            assert!((retained - expected).abs() < 1e-8 * expected.abs().max(1.0),
                "natural={natural}, ambient={temperature}: budget {retained} != physical {expected}");
            // Independent negative physical control: the previous code also
            // changed the radiation reservoir, which is a different experiment.
            changed.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                .radiation.as_mut().unwrap().surfaces[0].reservoir_temperature.value = temperature;
            let wrong = run(&fixture, &changed, 712 + 2 * i);
            assert!((maximum(&fixture, &wrong) - expected).abs() > 1e-3,
                "fixture must distinguish a reservoir substitution from an ambient perturbation");
        }
    }
}

#[path = "joint_envelope.rs"]
mod joint_envelope;
