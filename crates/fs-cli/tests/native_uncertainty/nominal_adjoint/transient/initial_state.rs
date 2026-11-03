//! Initial-state uncertainty reaches native trajectories and their original duals.
use super::*;

fn initial_fixture(enthalpy: bool, initially_hot: bool) -> Fixture {
    let mut result = if enthalpy { enthalpy_fixture() } else { fixture() };
    let time = result.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap();
    if let Some(storage) = &mut time.enthalpy {
        time.initial_temperature = TransientInitialTemperature::FromRegionalEnthalpy;
        storage.materials[0].initial_specific_enthalpy.value = if initially_hot { 368.5 } else { 0.0 };
    } else if initially_hot {
        time.initial_temperature = TransientInitialTemperature::Uniform(QtyAny::new(330.0, dims::TEMPERATURE));
    }
    write_project(&result);
    result
}

fn initial_study(fixture: &Fixture, enthalpy: bool, controlled: bool, peak: bool) -> String {
    let coordinate = match (enthalpy, peak) {
        (false, false) => "initial-temperature :entity \"initial-state\" :low 288.15K :high 298.15K",
        (false, true) => "initial-temperature :entity \"initial-state\" :low 325K :high 335K",
        (true, false) => "initial-specific-enthalpy :entity \"air\" :low -50J/kg :high 50J/kg",
        (true, true) => "initial-specific-enthalpy :entity \"air\" :low 318.5J/kg :high 418.5J/kg",
    };
    let mut source = study_text(false, false, controlled)
        .replace(":samples 4", ":samples 2")
        .replace(":wall-time 120s", ":wall-time 240s")
        .replace("convection-temperature :entity \"air\" :low 294K :high 300K", coordinate);
    if peak { source = source.replace(":qoi \"temperature-max\"", ":qoi \"temperature-peak\""); }
    std::fs::write(fixture.sources.join("study.fsim"), &source).unwrap();
    source
}

#[test]
fn g1_initial_state_mean_controls_match_native_resolves_and_preserve_raw_samples() {
    for enthalpy in [false, true] {
        let raw = initial_fixture(enthalpy, false);
        initial_study(&raw, enthalpy, false, false);
        let result = raw.study(None, fs_cli::exit::SUCCESS);
        let (_, raw_report) = raw.report(&result);
        let full = initial_fixture(enthalpy, false);
        let source = initial_study(&full, enthalpy, true, false);
        let result = full.study(None, fs_cli::exit::SUCCESS);
        let (_, report) = full.report(&result);
        assert_eq!(rows(&report), rows(&raw_report));
        assert_eq!(report.get("statistics"), raw_report.get("statistics"));
        assert_eq!(report.get("mean_control").unwrap().f64_field("probe_solves_completed"), Some(1.0));
        let coefficients = gradient(&report);
        assert!(coefficients[1] > 0.0, "the final field retains initial stored energy");

        // Fresh original projects come through the real sampler, then the
        // ordinary native solve; no derivative report or analytic stand-in.
        let study = fs_project::uncertainty::UncertaintyStudy::parse(&source).unwrap();
        let nominal = study.mean_control().unwrap().probe(study.parameters(), 0).unwrap();
        let bound = study.bind(&full.project).unwrap();
        let epsilon = if enthalpy { 0.1 } else { 0.01 };
        let mut values = Vec::new();
        for (side, sign) in [-1.0, 1.0].into_iter().enumerate() {
            let mut point = nominal.clone();
            point[1] += sign * epsilon;
            let project = bound.sample_project(&point).unwrap();
            values.push(full.solve_project(410 + side, &project).2);
        }
        let difference = (values[1] - values[0]) / (2.0 * epsilon);
        assert!((coefficients[1] - difference).abs() <= 1e-5 * difference.abs().max(1.0),
            "enthalpy={enthalpy}: original initial-state dual {} versus FD {difference}", coefficients[1]);
    }
}

#[test]
fn g1_initial_sampled_peaks_keep_direct_state_sensitivity_and_zero_future_power() {
    for enthalpy in [false, true] {
        let full = initial_fixture(enthalpy, true);
        initial_study(&full, enthalpy, true, true);
        let result = full.study(None, fs_cli::exit::SUCCESS);
        let (_, report) = full.report(&result);
        let coefficients = gradient(&report);
        assert_eq!(coefficients[0], 0.0);
        let expected = if enthalpy { 0.1 } else { 1.0 };
        assert!((coefficients[1] - expected).abs() < 1e-12,
            "initial-state sensitivity survives the empty reverse prefix");
        let branch = report.get("mean_control").unwrap().get("sampled_peak_branch").unwrap();
        assert_eq!(branch.f64_field("state_index"), Some(0.0));
        assert_eq!(branch.f64_field("control_prefix_steps"), Some(0.0));
        for row in rows(&report) {
            let input = row.get("parameters").unwrap().as_array().unwrap()[1].as_f64().unwrap();
            let initial = if enthalpy { 293.15 + 0.1 * input } else { input };
            assert!((row.f64_field("value_k").unwrap() - initial).abs() < 1e-10);
        }
    }
}
