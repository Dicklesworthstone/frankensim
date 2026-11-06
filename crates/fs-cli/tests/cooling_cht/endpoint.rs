//! Exercise endpoint coupling through the same scene parser and real binary.
use super::{J, f, run, scratch, unsteady_material_atmosphere};

const CONTROLS: &str = r#""endpoint_coupling":{
    "maximum_passes_per_step":80,"maximum_momentum_solves":500,
    "maximum_energy_solves":10000,"maximum_source_refreshes":1001,
    "temperature_tolerance_k":1e-8,"transport_relative_tolerance":1e-8,
    "force_tolerance_n_m3":1e-9,"temperature_relaxation":0.8},"#;

fn implicit(scene: &str) -> String {
    scene.replacen(
        r#""transient":{"#,
        &format!(r#""transient":{{{CONTROLS}"#),
        1,
    )
}

fn assert_closed(result: &J, steps: usize) {
    let report = result.get("endpoint_coupling").unwrap();
    assert_eq!(report.str_field("mode"), Some("implicit-endpoint"));
    assert_eq!(f(report, &["accepted_steps"]), steps as f64);
    assert_eq!(f(result, &["transient", "steps"]), steps as f64);
    assert!(f(report, &["momentum_solves"]) >= steps as f64);
    assert!(f(report, &["energy_solves"]) >= f(report, &["momentum_solves"]));
    for (defect, tolerance) in [
        ("maximum_temperature_defect_k", "temperature_tolerance_k"),
        (
            "maximum_transport_relative_defect",
            "transport_relative_tolerance",
        ),
        ("maximum_force_defect_n_m3", "force_tolerance_n_m3"),
    ] {
        assert!(f(report, &[defect]) <= f(report, &["controls", tolerance]));
    }
    assert!(f(result, &["transient", "worst_step_closure_j"]) < 1e-7);
}

#[test]
fn native_buoyancy_closes_the_first_heated_interval_from_the_scene() {
    let native = unsteady_material_atmosphere("300 100000 0", -9.81)
        .replace(
            r#""fluid":{"mode":"temperature-dependent","gas":"300 100000 0"}"#,
            r#""fluid":"dry-air-300k","reference_temperature_k":300"#,
        )
        .replace(r#""steps":4"#, r#""steps":1"#);
    let (code, lagged, error) = run(&scratch("endpoint-native-lagged.json", &native));
    assert_eq!(code, 0, "{error}");
    assert_eq!(f(&lagged, &["flow", "final_kinetic_energy_j"]), 0.0);
    let (code, coupled, error) = run(&scratch("endpoint-native-coupled.json", &implicit(&native)));
    assert_eq!(code, 0, "{error}");
    assert_eq!(
        coupled.str_field("solver"),
        Some("fv-simplec-unsteady-boussinesq-implicit")
    );
    assert!(f(&coupled, &["flow", "final_kinetic_energy_j"]) > 0.0);
    assert_closed(&coupled, 1);
    let records = coupled
        .path(&["endpoint_coupling", "records"])
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(records.len(), 1);
    assert!(f(&records[0], &["passes"]) >= 2.0);
}

#[test]
fn source_material_endpoint_coupling_composes_with_calorics_and_sgs() {
    for (index, caloric, sgs) in [(0, false, false), (1, true, false), (2, true, true)] {
        let mut scene = implicit(&unsteady_material_atmosphere("300 100000 0", -9.81));
        if caloric {
            scene = scene.replacen(
                '{',
                r#"{"fluid_enthalpy":{
                    "minimum_temperature_k":300,"maximum_temperature_k":310,
                    "maximum_flux_defect_w":1e-9,
                    "transient_controls":{"maximum_passes_per_step":80,"maximum_energy_solves":10000,
                        "temperature_tolerance_k":1e-9,"maximum_storage_defect_w":1e-9,
                        "maximum_energy_closure_j":1e-8,"temperature_relaxation":1}},"#,
                1,
            ).replace(r#""solver":{"#, r#""solver":{"energy_convection":"enthalpy-minmod","#);
        }
        if sgs {
            scene = scene.replace(
                r#""solver":{"#,
                r#""solver":{"sgs":{
                "model":"wale","coefficient":0.17,"filter_width_factor":1,
                "maximum_eddy_to_molecular":50,"usage":"exploration"},"turbulent_prandtl":0.7,"#,
            );
        }
        let (code, result, error) =
            run(&scratch(&format!("endpoint-material-{index}.json"), &scene));
        assert_eq!(code, 0, "{error}");
        assert_closed(&result, 4);
        let material = result.get("material_feedback").unwrap();
        assert_eq!(material.str_field("coupling"), Some("implicit-endpoint"));
        assert!((f(material, &["applied_time_s"]) - 0.1).abs() < 1e-14);
        assert!((f(material, &["endpoint_time_s"]) - 0.1).abs() < 1e-14);
        assert!(f(&result, &["flow", "final_kinetic_energy_j"]) > 0.0);
        if caloric {
            assert_eq!(
                material
                    .path(&["fluid_enthalpy", "advection"])
                    .and_then(J::as_str),
                Some("source-enthalpy-minmod")
            );
        }
        if sgs {
            assert_eq!(
                result.str_field("solver"),
                Some("fv-simplec-unsteady-material-sgs")
            );
            assert!(result.path(&["flow", "sgs"]).is_some());
        }
    }
}

#[test]
fn endpoint_controls_cannot_be_ignored_and_budget_exhaustion_publishes_no_march() {
    let scene = implicit(&unsteady_material_atmosphere("300 100000 0", -9.81));
    for (index, invalid) in [
        scene.replace(r#""flow":"unsteady""#, r#""flow":"frozen""#),
        scene.replace(
            r#""flow":"unsteady""#,
            r#""flow":"unsteady","energy":"steady-on-mean-flow""#,
        ),
        scene.replace(
            r#""maximum_passes_per_step":80"#,
            r#""maximum_passes_per_step":1"#,
        ),
        scene.replace(
            r#""maximum_momentum_solves":500"#,
            r#""maximum_momentum_solves":1"#,
        ),
        scene.replace(
            r#""temperature_relaxation":0.8"#,
            r#""temperature_relaxation":0"#,
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let (code, result, error) = run(&scratch(
            &format!("endpoint-refusal-{index}.json"),
            &invalid,
        ));
        assert_ne!(code, 0, "{error}");
        assert!(result.get("flow").is_none());
        assert!(result.get("endpoint_coupling").is_none());
        assert!(result.get("material_feedback").is_none());
    }
}
