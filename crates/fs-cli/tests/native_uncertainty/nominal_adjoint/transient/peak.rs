//! Native sampled-peak calibration, physical coefficient checks and retained resume.
use super::*;

fn peak_fixture(enthalpy: bool, initially_hot: bool, pulse: bool) -> Fixture {
    let mut fixture = if enthalpy { enthalpy_fixture() } else { fixture() };
    let time = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap();
    if initially_hot {
        if let Some(storage) = &mut time.enthalpy {
            // The declared sensible chart is 10 J/(kg K), with h=0 at 293.15 K.
            storage.materials[0].initial_specific_enthalpy.value = 368.5;
            time.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
        } else {
            time.initial_temperature = TransientInitialTemperature::Uniform(
                QtyAny::new(330.0, dims::TEMPERATURE));
        }
    }
    if pulse {
        time.power_schedules.push(TransientRegionPower {
            region: "air".into(),
            source: "synthetic finite heating pulse for independent peak controls".into(),
            steps: vec![
                TransientPowerStep {
                    until: QtyAny::new(0.5, dims::TIME),
                    watts: QtyAny::new(1000.0, dims::POWER),
                },
                TransientPowerStep {
                    until: QtyAny::new(2.0, dims::TIME),
                    watts: QtyAny::new(0.0, dims::POWER),
                },
            ],
        });
    }
    write_project(&fixture);
    fixture
}

fn peak_study(fixture: &Fixture, controlled: bool, qmc: bool) {
    let source = study(fixture, controlled, qmc)
        .replace(":qoi \"temperature-max\"", ":qoi \"temperature-peak\"");
    assert!(source.contains(":qoi \"temperature-peak\""));
    std::fs::write(fixture.sources.join("study.fsim"), source).unwrap();
}

fn branch(report: &JsonValue) -> &JsonValue {
    report.get("mean_control").unwrap().get("sampled_peak_branch").unwrap()
}

#[test]
fn g1_g3_initial_tied_peak_retains_zero_future_controls_and_identical_raw_samples() {
    for enthalpy in [false, true] {
        let raw = peak_fixture(enthalpy, true, false);
        peak_study(&raw, false, false);
        let raw_result = raw.study(None, fs_cli::exit::SUCCESS);
        let (_, raw_report) = raw.report(&raw_result);
        let controlled = peak_fixture(enthalpy, true, false);
        peak_study(&controlled, true, false);
        let result = controlled.study(None, fs_cli::exit::SUCCESS);
        let (_, report) = controlled.report(&result);
        assert_eq!(rows(&report), rows(&raw_report));
        assert_eq!(report.get("statistics"), raw_report.get("statistics"));
        assert_eq!(gradient(&report), [0.0, 0.0],
            "the genuine initial branch has no future power or reservoir influence");
        let selected = branch(&report);
        assert_eq!(selected.f64_field("state_index"), Some(0.0));
        assert_eq!(selected.f64_field("control_prefix_steps"), Some(0.0));
        assert_eq!(selected.f64_field("time_s"), Some(0.0));
        assert!(selected.f64_field("active_ties").unwrap() > 1.0);
        assert_eq!(selected.get("unique_sampled_maximum"), Some(&JsonValue::Bool(false)));
        assert_eq!(selected.get("unique_maximum_derivative_claim"), Some(&JsonValue::Bool(false)));
        assert_eq!(selected.str_field("coefficient_semantics"), Some("frozen-selected-branch"));
        let control = report.get("mean_control").unwrap();
        assert!(control.str_field("scope").unwrap().contains("fixed control coefficient"));
        assert_eq!(control.get("calibration").unwrap().get("sampled_peak_branch"), Some(selected));
        assert!(rows(&report).iter().all(|row|
            (row.f64_field("value_k").unwrap() - 330.0).abs() < 1e-10));
    }
}

#[test]
fn g1_g3_interior_peak_controls_match_fresh_native_peak_resolves_for_both_storage_laws() {
    for enthalpy in [false, true] {
        let raw = peak_fixture(enthalpy, false, true);
        peak_study(&raw, false, true);
        let raw_result = raw.study(None, fs_cli::exit::SUCCESS);
        let (_, raw_report) = raw.report(&raw_result);
        let controlled = peak_fixture(enthalpy, false, true);
        peak_study(&controlled, true, true);
        let result = controlled.study(None, fs_cli::exit::SUCCESS);
        let (_, report) = controlled.report(&result);
        assert_eq!(rows(&report), rows(&raw_report));
        assert_eq!(report.get("statistics"), raw_report.get("statistics"));
        assert_eq!(report.get("qmc"), raw_report.get("qmc"));
        let selected = branch(&report);
        assert_eq!(selected.f64_field("time_s"), Some(0.5));
        assert_eq!(selected.f64_field("state_index"), Some(2.0));
        assert_eq!(selected.f64_field("control_prefix_steps"), Some(2.0));
        assert_eq!(selected.f64_field("fine_steps"), Some(8.0));
        let coefficients = gradient(&report);
        assert_eq!(coefficients[0], 0.0, "scheduled absolute watts override static power");
        assert!(coefficients[1] > 0.0);
        let epsilon = 0.01;
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut project = controlled.project.clone();
            let fs_project::ThermalBoundaryCondition::Convection { reference_temperature, .. } =
                &mut project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
                    .boundaries[0].condition
            else { panic!("prescribed reference fixture"); };
            reference_temperature.value = 297.0 + sign * epsilon;
            let receipt = solve(&controlled, &project, 300 + side);
            values.push(receipt.path(&["transient", "sampled_peak"]).unwrap()
                .f64_field("value_k").unwrap());
        }
        let difference = (values[1] - values[0]) / (2.0 * epsilon);
        assert!((coefficients[1] - difference).abs() <= 2e-5 * difference.abs().max(1.0),
            "enthalpy={enthalpy}: peak branch {} versus independent peak FD {difference}",
            coefficients[1]);
    }
}

#[test]
fn g4_g5_peak_calibration_resumes_with_the_same_branch_without_original_sources() {
    let full = peak_fixture(true, false, true);
    peak_study(&full, true, true);
    let result = full.study(None, fs_cli::exit::SUCCESS);
    let (expected, report) = full.report(&result);
    let expected_branch = branch(&report).clone();

    let split = peak_fixture(true, false, true);
    peak_study(&split, true, true);
    let prefix = split.study(Some("1"), fs_cli::exit::BUDGET);
    let (_, prefix_report) = split.report(&prefix);
    assert!(rows(&prefix_report).is_empty(), "calibration is not a probability observation");
    assert_eq!(branch(&prefix_report), &expected_branch);
    let frozen = prefix_report.get("mean_control").unwrap().get("calibration").unwrap().clone();
    std::fs::rename(&split.sources, split.dir.join("retained-peak-inputs")).unwrap();
    let completed = split.resume(prefix.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    let (actual, resumed) = split.report(&completed);
    assert_eq!(branch(&resumed), &expected_branch);
    assert_eq!(resumed.get("mean_control").unwrap().get("calibration"), Some(&frozen));
    assert_eq!(actual, expected, "split calibration and sampling reproduce uninterrupted bytes");
}
