//! G0/G1/G3/G4/G5: manufactured joint inputs consume the same source chain
//! through both accepted capacity and total-enthalpy trajectories.
use super::*;
use fs_project::{
    ConductionTransient, TransientEnthalpy, TransientEnthalpyKnot,
    TransientInitialTemperature, TransientRegionCapacity, TransientRegionEnthalpy, spec::dims,
};
use fs_qty::QtyAny;

const FINAL_MATERIAL: &str = "temperature-final-material-adjoint";

fn fixture(enthalpy: bool) -> Fixture {
    let fixture = Fixture::new();
    let path = fixture.inputs.join("contact.fsim");
    let mut project = fs_project::parse_sexpr(&std::fs::read_to_string(&path).unwrap())
        .unwrap().spec;
    let materials = project.materials.as_ref().unwrap();
    let setup = project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let mut time = ConductionTransient {
        initial_temperature: TransientInitialTemperature::Uniform(QtyAny::new(
            293.15, dims::TEMPERATURE)),
        horizon: QtyAny::new(2.0, dims::TIME),
        max_step: QtyAny::new(0.5, dims::TIME),
        max_steps: 12,
        energy_tolerance: QtyAny::new(1e-6, dims::ENERGY),
        capacities: setup.regions.iter().map(|region| TransientRegionCapacity {
            region: region.region.clone(),
            volumetric_heat_capacity: QtyAny::new(5.0, dims::VOLUMETRIC_HEAT_CAPACITY),
            source: "synthetic contact-control storage".into(),
        }).collect(),
        enthalpy: None,
        power_schedules: Vec::new(),
    };
    if enthalpy {
        time.capacities.clear();
        time.enthalpy = Some(TransientEnthalpy {
            max_newton_iterations: 32,
            max_backtracks: 24,
            materials: materials.iter().map(|binding| TransientRegionEnthalpy {
                region: binding.region.clone(),
                material_card_identity: ContentHash::from_hex(&binding.card).unwrap(),
                source: "synthetic sensible enthalpy for contact uncertainty".into(),
                reference_density: QtyAny::new(1.0, dims::DENSITY),
                initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
                knots: [(-250.0, 243.15), (500.0, 393.15)].into_iter()
                    .map(|(h, temperature)| TransientEnthalpyKnot {
                        specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                        temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                        liquid_mass_fraction: 0.0,
                    }).collect(),
                phase_conductivity: None,
            }).collect(),
        });
    }
    setup.transient = Some(time);
    project.solver.as_mut().unwrap().tolerance_rel = 1e-11;
    std::fs::write(path, fs_project::print_sexpr(&project).unwrap()).unwrap();
    fixture
}

fn controlled() -> String {
    SOURCE.replace(":version 1", ":version 3 :mean-control (nominal-adjoint :max-solves 1)")
}

#[test]
fn g1_g3_contact_mean_controls_follow_capacity_and_enthalpy_trajectories() {
    for enthalpy in [false, true] {
        let f = fixture(enthalpy);
        let (raw, raw_path) = f.run(SOURCE, "raw-transient", crate::exit::SUCCESS);
        let (calibrated, path) = f.run(&controlled(), "controlled-transient", crate::exit::SUCCESS);
        let raw = Fixture::report(&raw, &raw_path);
        let calibrated = Fixture::report(&calibrated, &path);
        assert_eq!(calibrated.get("observations"), raw.get("observations"));
        assert_eq!(calibrated.get("statistics"), raw.get("statistics"));
        let control = calibrated.get("mean_control").unwrap();
        assert_eq!(control.f64_field("probe_solves_planned"), Some(1.0));
        assert_eq!(control.f64_field("probe_solves_completed"), Some(1.0));
        assert_eq!(control.f64_field("probe_solves_attempted"), Some(1.0));
        let gradient: Vec<_> = control.get("gradient").unwrap().as_array().unwrap()
            .iter().map(|value| value.as_f64().unwrap()).collect();
        assert!(gradient[0] < 0.0 && gradient[1] > 0.0);
        let model = Model::load(&f.inputs.join("controlled-transient.fsim")).unwrap();
        let view = model.nominal_view().unwrap();
        let outputs = view.bound.base().outputs.as_ref().unwrap();
        assert_eq!(outputs.iter().filter(|output| output.name.ends_with("-adjoint")).count(), 1);
        assert!(outputs.iter().any(|output| output.name == FINAL_MATERIAL));
        assert!(model.bound.base().outputs.as_ref().unwrap().iter()
            .all(|output| !output.name.ends_with("-adjoint")));
        let ledger = Ledger::open(path.to_str().unwrap()).unwrap();
        let gate = CancelGate::new_clock_free();
        // The physical mean uses the 0.20 -> 0.10 source segment. The original
        // project deliberately uses a different segment at 1.5 MPa.
        for (direction, epsilon) in [([1.0, 0.0], 100.0), ([2000.0, 0.2], 0.02)] {
            let mean = [900000.0, 5.0];
            let lo: Vec<_> = mean.iter().zip(direction).map(|(m, d)| m - epsilon * d).collect();
            let hi: Vec<_> = mean.iter().zip(direction).map(|(m, d)| m + epsilon * d).collect();
            let low = model.sample(&ledger, &gate, &lo, 120.0).unwrap().unwrap();
            let high = model.sample(&ledger, &gate, &hi, 120.0).unwrap().unwrap();
            let expected = (high.value_k - low.value_k) / (2.0 * epsilon);
            let actual: f64 = gradient.iter().zip(direction).map(|(g, d)| g * d).sum();
            assert!((actual - expected).abs() <= 5e-4 * expected.abs().max(1e-9),
                "enthalpy={enthalpy}, direction={direction:?}: adjoint={actual:e}, FD={expected:e}");
        }
    }
}

#[test]
fn g4_g5_contact_history_calibration_resumes_and_preserves_source_slope_refusals() {
    let f = fixture(true);
    let source = controlled();
    let (full, full_path) = f.run(&source, "full-transient", crate::exit::SUCCESS);
    let expected = Fixture::report(&full, &full_path);
    let path = f.inputs.join("resume-transient.fsim");
    std::fs::write(&path, &source).unwrap();
    let ledger = f.root.join("resume-transient.db");
    let prefix = crate::study::study_path(&path, &ledger, Some("1"), crate::OutputMode::Json);
    assert_eq!(prefix.exit_code, crate::exit::BUDGET, "{}", prefix.stderr);
    let prefix = J::parse(&prefix.stdout).unwrap();
    let report = Fixture::report(&prefix, &ledger);
    assert!(report.get("observations").unwrap().as_array().unwrap().is_empty());
    assert_eq!(report.get("mean_control").unwrap().f64_field("probe_solves_completed"), Some(1.0));
    let kink = source.replace("700000Pa", "1100000Pa").replace("1100000Pa)", "1300000Pa)");
    let (failed, failed_path) = f.run(&kink, "kink-transient", crate::exit::REFUSED);
    let failed = Fixture::report(&failed, &failed_path);
    assert!(failed.get("observations").unwrap().as_array().unwrap().is_empty());
    assert!(failed.str_field("failure").unwrap().contains("unequal-slope"));
    std::fs::rename(&f.inputs, f.root.join("retained-transient-inputs")).unwrap();
    let finished = crate::study::resume_path(prefix.str_field("run").unwrap(),
        &ledger, None, crate::OutputMode::Json);
    assert_eq!(finished.exit_code, crate::exit::SUCCESS, "{}", finished.stderr);
    let resumed = Fixture::report(&J::parse(&finished.stdout).unwrap(), &ledger);
    assert_eq!(resumed.get("mean_control"), expected.get("mean_control"));
    assert_eq!(resumed.get("observations"), expected.get("observations"));
}

#[test]
fn g0_contact_history_rows_bind_the_complete_prefix_scale_units_and_owner() {
    for schema in ["fs-cli-native-transient-adjoint-v1", "fs-cli-native-enthalpy-adjoint-v1"] {
        let row = r#"{"interface":"joint","interface_card":"card","mapped":false,"face_pairs":2,"reference_value":1,"log_reference_value":0,"dtemperature_dresistance_multiplier_k":0.4}"#;
        let source = format!(r#"{{"schema":"{schema}","output":"{FINAL_MATERIAL}","state_index":8,"material_controls":{{"control_steps":8,"units":"K","parameter_scope":"constant-through-selected-trajectory-independent-regional-effective-conductivity-and-named-matching-resistance-multipliers","contact_interfaces":[{row}]}}}}"#);
        assert_eq!(transient_coefficient(&J::parse(&source).unwrap(), "joint").unwrap(), 0.4);
        for bad in [
            source.replace(FINAL_MATERIAL, "temperature-final-adjoint"),
            source.replace("\"control_steps\":8", "\"control_steps\":7"),
            source.replace("\"units\":\"K\"", "\"units\":\"1\""),
            source.replace("constant-through-selected-trajectory", "one-endpoint"),
            source.replace("\"joint\"", "\"another-joint\""),
            source.replace("\"reference_value\":1", "\"reference_value\":2"),
            source.replace("\"log_reference_value\":0", "\"log_reference_value\":1"),
            source.replace("\"mapped\":false", "\"mapped\":true"),
            source.replace("\"face_pairs\":2", "\"face_pairs\":0"),
            source.replace(row, &format!("{row},{row}")),
        ] {
            assert!(transient_coefficient(&J::parse(&bad).unwrap(), "joint").is_err(), "{bad}");
        }
    }
}
