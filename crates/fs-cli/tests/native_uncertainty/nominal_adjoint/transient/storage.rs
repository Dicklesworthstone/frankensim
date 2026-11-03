//! Original heat-storage declarations feed actual samples and native trajectory duals.
use super::*;

fn storage_fixture(enthalpy: bool, initially_hot: bool) -> Fixture {
    let mut result = if enthalpy {
        enthalpy_fixture()
    } else {
        fixture()
    };
    if initially_hot {
        let time = result
            .project
            .cooling
            .as_mut()
            .unwrap()
            .conduction
            .as_mut()
            .unwrap()
            .transient
            .as_mut()
            .unwrap();
        if let Some(storage) = &mut time.enthalpy {
            time.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
            storage.materials[0].initial_specific_enthalpy.value = 368.5;
        } else {
            time.initial_temperature =
                TransientInitialTemperature::Uniform(QtyAny::new(330.0, dims::TEMPERATURE));
        }
        write_project(&result);
    }
    result
}

fn storage_study(fixture: &Fixture, enthalpy: bool, controlled: bool, peak: bool) -> String {
    let coordinate = if enthalpy {
        "enthalpy-reference-density :entity \"air\" :low 10kg/m^3 :high 20kg/m^3"
    } else {
        "thermal-capacity-scale :entity \"thermal-storage\" :low 1 :high 2"
    };
    let mut source = study_text(false, false, controlled)
        .replace(":wall-time 120s", ":wall-time 240s")
        .replace(
            "convection-temperature :entity \"air\" :low 294K :high 300K",
            coordinate,
        );
    if peak {
        source = source.replace(":qoi \"temperature-max\"", ":qoi \"temperature-peak\"");
    }
    std::fs::write(fixture.sources.join("study.fsim"), &source).unwrap();
    source
}

#[test]
fn g1_storage_controls_match_fresh_native_trajectories_and_keep_raw_samples() {
    for enthalpy in [false, true] {
        let raw = storage_fixture(enthalpy, false);
        storage_study(&raw, enthalpy, false, false);
        let result = raw.study(None, fs_cli::exit::SUCCESS);
        let (_, raw_report) = raw.report(&result);
        let full = storage_fixture(enthalpy, false);
        let source = storage_study(&full, enthalpy, true, false);
        let result = full.study(None, fs_cli::exit::SUCCESS);
        let (_, report) = full.report(&result);
        assert_eq!(rows(&report), rows(&raw_report));
        assert_eq!(report.get("statistics"), raw_report.get("statistics"));
        assert_eq!(
            report
                .get("mean_control")
                .unwrap()
                .f64_field("probe_solves_completed"),
            Some(1.0)
        );
        let coefficients = gradient(&report);
        assert!(
            coefficients[1] < 0.0,
            "more storage slows this originally uniform body's heating"
        );

        let study = fs_project::uncertainty::UncertaintyStudy::parse(&source).unwrap();
        let nominal = study
            .mean_control()
            .unwrap()
            .probe(study.parameters(), 0)
            .unwrap();
        assert_eq!(
            nominal[1],
            if enthalpy { 15.0 } else { 1.5 },
            "a nonunit capacity scale exercises the normalized-to-authored chain rule"
        );
        let bound = study.bind(&full.project).unwrap();
        let epsilon = if enthalpy { 0.01 } else { 0.001 };
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut point = nominal.clone();
            point[1] += sign * epsilon;
            let project = bound.sample_project(&point).unwrap();
            values.push(full.solve_project(430 + side, &project).2);
        }
        let difference = (values[1] - values[0]) / (2.0 * epsilon);
        assert!(
            (coefficients[1] - difference).abs() <= 1e-5 * difference.abs().max(1.0),
            "enthalpy={enthalpy}: original storage dual {} versus FD {difference}",
            coefficients[1]
        );
    }
}

#[test]
fn g1_initial_peaks_have_zero_storage_sensitivity_when_the_original_h_t_chart_is_fixed() {
    for enthalpy in [false, true] {
        let full = storage_fixture(enthalpy, true);
        storage_study(&full, enthalpy, true, true);
        let result = full.study(None, fs_cli::exit::SUCCESS);
        let (_, report) = full.report(&result);
        assert_eq!(
            gradient(&report),
            [0.0, 0.0],
            "storage coefficients and future power cannot move the authored initial temperature"
        );
        let branch = report
            .get("mean_control")
            .unwrap()
            .get("sampled_peak_branch")
            .unwrap();
        assert_eq!(branch.f64_field("state_index"), Some(0.0));
        assert_eq!(branch.f64_field("control_prefix_steps"), Some(0.0));
        for row in rows(&report) {
            assert!((row.f64_field("value_k").unwrap() - 330.0).abs() < 1e-10);
        }
    }
}
