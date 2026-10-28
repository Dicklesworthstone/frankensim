//! Actual native capacity/enthalpy trajectories and one-probe amplitude control.
use super::*;

fn fixture(enthalpy: bool) -> Fixture {
    let mut fixture = pulse_fixture();
    if enthalpy {
        use fs_project::{TransientEnthalpy, TransientEnthalpyKnot, TransientRegionEnthalpy};
        let card = fs_blake3::ContentHash::from_hex(&fixture.project.materials.as_ref().unwrap()[0].card).unwrap();
        let time = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap();
        time.capacities.clear();
        time.enthalpy = Some(TransientEnthalpy {
            max_newton_iterations: 32, max_backtracks: 24,
            materials: vec![TransientRegionEnthalpy {
                region: "air".into(), material_card_identity: card,
                source: "synthetic sensible enthalpy for amplitude-chain integration".into(),
                reference_density: QtyAny::new(10.0, dims::DENSITY),
                initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
                knots: [(-500.0, 243.15), (1000.0, 393.15)].into_iter().map(|(h, t)| TransientEnthalpyKnot {
                    specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                    temperature: QtyAny::new(t, dims::TEMPERATURE), liquid_mass_fraction: 0.0,
                }).collect(), phase_conductivity: None,
            }],
        });
        write_project(&fixture);
    }
    fixture
}
fn source(peak: bool, controlled: bool) -> String {
    let mut source = study_source();
    if !peak { source = source.replace("temperature-peak", "temperature-max"); }
    if controlled {
        source = source.replace(":version 1", ":version 3 :mean-control (nominal-adjoint :max-solves 1)");
    }
    source
}
fn solve(fixture: &Fixture, project: &fs_project::ProjectSpec, ordinal: usize) -> JsonValue {
    let path = fixture.dir.join(format!("schedule-probe-{ordinal}.fsim"));
    std::fs::write(&path, fs_project::print_sexpr(project).unwrap()).unwrap();
    let ledger = fixture.dir.join("schedule-adjoint-reference.db");
    command(&["--json", "import", path.to_str().unwrap(),
        fixture.sources.join("plate.stl").to_str().unwrap(), ledger.to_str().unwrap(),
        "--unit", "m", "--max-hole-edges", "0"], fs_cli::exit::SUCCESS);
    let result = command(&["--json", "solve", path.to_str().unwrap(), ledger.to_str().unwrap(),
        "--materials", fixture.sources.join("aa6061.fsmcdpk").to_str().unwrap()], fs_cli::exit::SUCCESS);
    let receipt = json_artifact(&ledger, result.str_field("run_receipt").unwrap());
    let conduction = receipt.get("stages").unwrap().as_array().unwrap().iter()
        .find(|stage| stage.str_field("stage") == Some("conduction")).unwrap();
    json_artifact(&ledger, conduction.str_field("receipt").unwrap())
}
fn value(receipt: &JsonValue, peak: bool) -> f64 {
    if peak {
        receipt.path(&["transient", "sampled_peak", "value_k"]).unwrap().as_f64().unwrap()
    } else {
        receipt.path(&["transient", "temporal_error", "fine_final_k"]).unwrap().as_f64().unwrap()
    }
}

#[test]
fn g1_g3_native_amplitude_adjoint_matches_physical_schedule_probes_and_leaves_samples_unchanged() {
    for enthalpy in [false, true] {
        for peak in [false, true] {
            let fixture = fixture(enthalpy);
            std::fs::write(fixture.sources.join("study.fsim"), source(peak, false)).unwrap();
            let raw_result = fixture.study(None, fs_cli::exit::SUCCESS);
            let (_, raw) = fixture.report(&raw_result);
            std::fs::write(fixture.sources.join("study.fsim"), source(peak, true)).unwrap();
            let result = fixture.study(None, fs_cli::exit::SUCCESS);
            let (_, report) = fixture.report(&result);
            assert_eq!(rows(&report), rows(&raw));
            assert_eq!(report.get("statistics"), raw.get("statistics"));
            let control = report.get("mean_control").unwrap();
            assert_eq!(control.f64_field("probe_solves_completed"), Some(1.0));
            let gradient = control.get("gradient").unwrap().as_array().unwrap();
            let derivative = gradient[0].as_f64().unwrap();
            assert!(derivative > 0.0);
            assert_eq!(gradient[1].as_f64(), Some(0.0));
            let h = 0.01;
            let low = solve(&fixture, &scaled(&fixture, 1.0 - h), 10);
            let high = solve(&fixture, &scaled(&fixture, 1.0 + h), 11);
            let difference = (value(&high, peak) - value(&low, peak)) / (2.0 * h);
            assert!((derivative - difference).abs() <= 1e-5 * difference.abs().max(1.0),
                "enthalpy={enthalpy}, peak={peak}: adjoint={derivative} FD={difference}");
            assert!(control.get("estimate").unwrap().f64_field("variance_ratio").unwrap() < 1e-4,
                "this linear fixture's pulse response is affine in amplitude");
        }
    }
}

#[test]
fn g4_g5_amplitude_calibration_recovers_the_original_workload_without_reopening_sources() {
    let full = fixture(true);
    std::fs::write(full.sources.join("study.fsim"), source(true, true)).unwrap();
    let completed = full.study(None, fs_cli::exit::SUCCESS);
    let expected = full.report(&completed).0;
    let split = fixture(true);
    std::fs::write(split.sources.join("study.fsim"), source(true, true)).unwrap();
    let prefix = split.study(Some("1"), fs_cli::exit::BUDGET);
    let (_, initial) = split.report(&prefix);
    assert!(rows(&initial).is_empty(), "the one physical adjoint probe is not a probability sample");
    let calibration = initial.path(&["mean_control", "calibration"]).unwrap().clone();
    std::fs::rename(&split.sources, split.dir.join("relocated-amplitude-inputs")).unwrap();
    let partial = split.resume(prefix.str_field("run").unwrap(), Some("1"), fs_cli::exit::BUDGET);
    assert_eq!(split.report(&partial).1.path(&["mean_control", "calibration"]), Some(&calibration));
    let completed = split.resume(partial.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    assert_eq!(split.report(&completed).0, expected);
}
