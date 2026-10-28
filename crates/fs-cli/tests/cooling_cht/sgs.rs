//! The real CLI selects the shared transient SGS owner and publishes its
//! modeled momentum/thermal fields. These are workflow, not turbulence-DNS tests.

use super::{DUCT, J, f, run, scratch};

const EXAMPLE: &str = include_str!("../../../../examples/cooling-cht/sgs-heated-obstacle.json");
const VREMAN: &str = include_str!("../../../../examples/cooling-cht/vreman-graded-obstacle.json");

fn short_scene(name: &str, model: &str, coefficient: f64) -> String {
    EXAMPLE
        .replace("\"steps\": 40", "\"steps\": 2")
        .replace("\"model\": \"wale\"", &format!("\"model\": \"{model}\""))
        .replace(
            "\"coefficient\": 0.5",
            &format!("\"coefficient\": {coefficient}"),
        )
        .replace("sgs-fields.vtr", &format!("{name}.vtr"))
}

#[test]
fn declared_models_run_coupled_cooling_and_export_the_final_closure() {
    for (model, coefficient) in [("wale", 0.5), ("smagorinsky", 0.17), ("vreman", 0.07)] {
        let name = format!("sgs-{model}");
        let scene = short_scene(&name, model, coefficient);
        let (code, result, error) = run(&scratch(&format!("{name}.json"), &scene));
        assert_eq!(code, 0, "{error}");
        assert_eq!(result.str_field("solver"), Some("fv-simplec-unsteady-sgs"));
        assert_eq!(
            result.path(&["flow", "turbulence"]).and_then(J::as_str),
            Some(model)
        );
        let sgs = result.path(&["flow", "sgs"]).unwrap();
        assert_eq!(sgs.str_field("model"), Some(model));
        assert_eq!(sgs.str_field("usage"), Some("exploration"));
        assert_eq!(f(sgs, &["coefficient"]), coefficient);
        assert_eq!(f(sgs, &["maximum_eddy_to_molecular"]), 50.0);
        assert!((f(sgs, &["min_filter_width_m"]) - 0.002).abs() < 1e-15);
        assert!((f(sgs, &["max_filter_width_m"]) - 0.002).abs() < 1e-15);
        let maximum = f(&result, &["flow", "max_eddy_viscosity_ratio"]);
        let mean = f(sgs, &["final_volume_mean_eddy_viscosity_ratio"]);
        assert!(maximum > 0.0 && maximum <= 50.0, "{maximum}");
        assert!(mean > 0.0 && mean <= maximum, "mean={mean}, max={maximum}");
        let fluid = result.path(&["material_inputs", "fluid"]).unwrap();
        let expected_k = f(fluid, &["density_kg_m3"])
            * f(fluid, &["specific_heat_j_kg_k"])
            * f(sgs, &["final_max_eddy_viscosity_m2_s"])
            / f(sgs, &["thermal_prandtl"]);
        assert!(
            (f(sgs, &["final_max_turbulent_conductivity_w_m_k"]) / expected_k - 1.0).abs() < 1e-12
        );
        assert!(f(&result, &["transient", "worst_step_closure_j"]) < 1e-8);
        assert!(f(&result, &["max_solid_temperature_k"]) > 300.0);
        assert_eq!(result.get("probes").and_then(J::as_array).unwrap().len(), 2);
        let fields = std::fs::read_to_string(result.str_field("vtk").unwrap()).unwrap();
        for name in [
            "temperature_k",
            "velocity_m_s",
            "eddy_viscosity_m2_s",
            "turbulent_conductivity_w_m_k",
        ] {
            assert!(
                fields.contains(&format!("Name=\"{name}\"")),
                "missing {name}"
            );
        }
    }
}

#[test]
fn vreman_graded_scene_reports_actual_directional_filter_widths() {
    let scene = VREMAN.replace("\"steps\": 40", "\"steps\": 2");
    let (code, result, error) = run(&scratch("vreman-graded.json", &scene));
    assert_eq!(code, 0, "{error}");
    let sgs = result.path(&["flow", "sgs"]).unwrap();
    assert_eq!(sgs.str_field("model"), Some("vreman"));
    assert_eq!(sgs.str_field("filter_width_definition"), Some("directional-cell-widths"));
    assert_eq!(f(sgs, &["coefficient"]), 0.07);
    for (key, expected) in [
        ("min_filter_width_by_axis_m", [0.002, 0.001, 0.002]),
        ("max_filter_width_by_axis_m", [0.004, 0.003, 0.002]),
    ] {
        let axes = sgs.get(key).and_then(J::as_array).unwrap();
        assert_eq!(axes.len(), 3);
        for (value, width) in axes.iter().zip(expected) {
            let actual = value.as_f64().unwrap();
            assert!((actual - width).abs() < 1e-14, "{key}: {actual} != {width}");
        }
    }
    assert!((f(sgs, &["min_filter_width_m"]) - 0.001).abs() < 1e-14);
    assert!((f(sgs, &["max_filter_width_m"]) - 0.004).abs() < 1e-14);
    assert!(f(sgs, &["final_max_eddy_viscosity_m2_s"]) > 0.0);
    assert!(f(&result, &["transient", "worst_step_closure_j"]) < 1e-8);
}

#[test]
fn sgs_controls_refuse_incomplete_or_inconsistent_solver_ownership() {
    let base = short_scene("sgs-refusal", "wale", 0.5);
    for (i, (from, to, reason)) in [
        (
            "\"transient\":",
            "\"unused_transient\":",
            "requires transient.flow",
        ),
        (
            "\"flow\": \"unsteady\"",
            "\"flow\": \"frozen\"",
            "requires transient.flow",
        ),
        (
            "\"energy\": \"march\"",
            "\"energy\": \"steady-on-mean-flow\"",
            "velocity-temperature correlation",
        ),
        (
            "\"solver\": {",
            "\"solver\": {\"turbulence\": \"lvel\",",
            "cannot be combined",
        ),
        (
            "\"usage\": \"exploration\"",
            "\"usage\": \"acceptance-render\"",
            "AcceptanceRenderForbidden",
        ),
        (
            "\"model\": \"wale\"",
            "\"model\": \"unimplemented\"",
            "solver.sgs.model",
        ),
        (
            "\"coefficient\": 0.5",
            "\"coefficient\": -0.5",
            "InvalidCoefficient",
        ),
        (
            "\"filter_width_factor\": 1.0",
            "\"filter_width_factor\": 0.0",
            "InvalidFilterWidth",
        ),
        (
            "\"maximum_eddy_to_molecular\": 50.0",
            "\"maximum_eddy_to_molecular\": -1.0",
            "InvalidViscosityRatio",
        ),
        (
            "\"usage\": \"exploration\"",
            "\"use\": \"exploration\"",
            "unsupported solver.sgs field",
        ),
        (
            "\"sgs\": {",
            "\"sgs\": null, \"unused_sgs\": {",
            "solver.sgs must be an object",
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let input = base.replace(from, to);
        assert_ne!(input, base, "fixture mutation {i}");
        let (code, result, error) = run(&scratch(&format!("sgs-refusal-{i}.json"), &input));
        assert_ne!(code, 0, "unexpected acceptance of {input}");
        assert_eq!(result.str_field("code"), Some("cooling-cht-input"));
        assert!(error.contains(reason), "{error}");
    }
    // A valid cap of zero is accepted as a control, then the actual nonzero
    // closure refuses at the cell where it exceeds that cap; it is not clipped.
    let capped = base.replace(
        "\"maximum_eddy_to_molecular\": 50.0",
        "\"maximum_eddy_to_molecular\": 0.0",
    );
    let (code, _, error) = run(&scratch("sgs-cap.json", &capped));
    assert_ne!(code, 0);
    assert!(error.contains("eddy/molecular ratio"), "{error}");
}

#[test]
fn explicit_zero_coefficient_keeps_the_existing_cli_trajectory() {
    let base = DUCT.replace(
        "\"conductivity_w_m_k\": 167.0}",
        "\"conductivity_w_m_k\": 167.0, \"volumetric_heat_capacity_j_m3_k\": 2.4e6}",
    ).replace(
        "\"solver\": {\"tolerance\": 1e-8}",
        "\"solver\": {\"tolerance\": 1e-8}, \"transient\": {\"flow\": \"unsteady\", \"time_step_s\": 0.01, \"steps\": 2}",
    );
    let disabled = base.replace(
        "\"solver\": {\"tolerance\": 1e-8}",
        "\"solver\": {\"tolerance\": 1e-8, \"sgs\": {\"model\": \"smagorinsky\", \"coefficient\": 0.0, \"filter_width_factor\": 1.0, \"maximum_eddy_to_molecular\": 0.0, \"usage\": \"exploration\"}}",
    );
    assert_ne!(base, disabled);
    let (code, laminar, error) = run(&scratch("sgs-zero-base.json", &base));
    assert_eq!(code, 0, "{error}");
    for model in ["smagorinsky", "vreman"] {
        let scene = disabled.replace("\"model\": \"smagorinsky\"", &format!("\"model\": \"{model}\""));
        let (code, zero, error) = run(&scratch(&format!("sgs-zero-{model}.json"), &scene));
        assert_eq!(code, 0, "{error}");
        assert_eq!(laminar.get("transient"), zero.get("transient"));
        assert_eq!(laminar.get("sources"), zero.get("sources"));
        for key in [
            "final_kinetic_energy_j",
            "final_mass_residual",
            "final_momentum_residual",
            "inflow_m3_s",
        ] {
            assert_eq!(laminar.path(&["flow", key]), zero.path(&["flow", key]));
        }
        assert_eq!(f(&zero, &["flow", "max_eddy_viscosity_ratio"]), 0.0);
    }
}
