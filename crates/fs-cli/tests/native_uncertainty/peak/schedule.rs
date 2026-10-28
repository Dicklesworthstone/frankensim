//! Actual native schedule sampling, pulse sizing, project export and recovery.
use super::*;
use fs_project::{TransientPowerStep, TransientRegionPower, ThermalBoundaryCondition};

fn pulse_fixture() -> Fixture {
    let mut fixture = super::fixture();
    fixture.project.power.as_mut().unwrap()[0].duty = 0.25;
    let setup = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    let time = setup.transient.as_mut().unwrap();
    time.initial_temperature = TransientInitialTemperature::Uniform(QtyAny::new(293.15, dims::TEMPERATURE));
    time.max_steps = 24;
    time.power_schedules = vec![TransientRegionPower {
        region: "air".into(), source: "synthetic 2 J startup pulse".into(),
        steps: [(0.25, 8.0), (2.0, 0.0)].into_iter().map(|(until, watts)| TransientPowerStep {
            until: QtyAny::new(until, dims::TIME), watts: QtyAny::new(watts, dims::POWER),
        }).collect(),
    }];
    for boundary in &mut setup.boundaries {
        if let ThermalBoundaryCondition::Convection { reference_temperature, .. } = &mut boundary.condition {
            reference_temperature.value = 293.15;
        }
    }
    write_project(&fixture);
    fixture
}
fn write_project(fixture: &Fixture) {
    std::fs::write(fixture.sources.join("cooling-reference.fsim"),
        fs_project::print_sexpr(&fixture.project).unwrap()).unwrap();
}
fn study_source() -> String {
    STUDY.replace("temperature-max", "temperature-peak")
        .replace(":target power ", ":target power-schedule-scale ")
        .replace(":low 4W :high 6W", ":low 0 :high 2")
        .replace(":low 294K :high 300K", ":low 293.15K :high 293.15K")
}
fn scaled(fixture: &Fixture, amplitude: f64) -> fs_project::ProjectSpec {
    // Independent physical edit, not BoundStudy::sample_project.
    let mut project = fixture.project.clone();
    for step in &mut project.cooling.as_mut().unwrap().conduction.as_mut().unwrap()
        .transient.as_mut().unwrap().power_schedules[0].steps {
        step.watts.value *= amplitude;
    }
    project
}

#[test]
fn g1_schedule_samples_match_independently_scaled_native_trajectories() {
    let fixture = pulse_fixture();
    let source = study_source();
    std::fs::write(fixture.sources.join("study.fsim"), source.replace("temperature-peak", "temperature-max")).unwrap();
    let final_result = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, final_report) = fixture.report(&final_result);
    std::fs::write(fixture.sources.join("study.fsim"), source).unwrap();
    let peak_result = fixture.study(None, fs_cli::exit::SUCCESS);
    let (_, peak_report) = fixture.report(&peak_result);
    for (index, (peak, final_row)) in rows(&peak_report).iter().zip(rows(&final_report)).enumerate() {
        let amplitude = peak.get("parameters").unwrap().as_array().unwrap()[0].as_f64().unwrap();
        let (hash, run, value) = fixture.solve_project(900 + index, &scaled(&fixture, amplitude));
        assert_eq!(peak.str_field("project_hash"), Some(hash.as_str()));
        assert_eq!(peak.str_field("run"), Some(run.as_str()));
        assert_eq!(peak.str_field("run"), final_row.str_field("run"));
        assert_eq!(final_row.f64_field("value_k").unwrap().to_bits(), value.to_bits());
        assert!(peak.f64_field("value_k").unwrap() > value, "cooling must not erase the startup peak");
        assert_eq!(fixture.project.power.as_ref().unwrap()[0].duty, 0.25);
    }
}

const DESIGN: &str = r#"(fsim-native-design-study
    :version 1 :project "cooling-reference.fsim" :evaluations 32 :seed 29
    :wall-time 240s :qoi "temperature-max"
    :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
    :materials ("aa6061.fsmcdpk") :interfaces ()
    :control (power-schedule-scale :entity "air" :minimum 0 :maximum 2 :tolerance 0.001)
    :temperature-tolerance 0.001K :assumption nondecreasing-temperature-with-schedule-scale)"#;

fn design_fixture() -> Fixture {
    let mut fixture = pulse_fixture();
    // Observe the original pulse independently of the design search. The
    // singleton probability coordinate leaves its physical schedule unchanged.
    let reference = study_source().replace(":low 0 :high 2", ":low 1 :high 1");
    std::fs::write(fixture.sources.join("study.fsim"), reference).unwrap();
    let result = fixture.study(None, fs_cli::exit::SUCCESS);
    let value = rows(&fixture.report(&result).1)[0].f64_field("value_k").unwrap();
    fixture.project.requirements.as_mut().unwrap()[0].limit.value = value + 1.0;
    fixture.project.requirements.as_mut().unwrap()[0].margin.value = 1.0;
    write_project(&fixture);
    std::fs::write(fixture.sources.join("study.fsim"), DESIGN).unwrap();
    fixture
}
fn report(fixture: &Fixture, result: &JsonValue) -> JsonValue {
    let receipt = result.get("receipt").unwrap();
    assert_eq!(receipt.str_field("driver"), Some("native-cooling-schedule-design-v1"));
    json_artifact(&fixture.ledger, receipt.str_field("report_json").unwrap())
}

#[test]
fn g1_g4_g5_pulse_design_exports_the_evaluated_schedule_and_resumes_without_sources() {
    let full = design_fixture();
    let result = full.study(None, fs_cli::exit::SUCCESS);
    let expected = report(&full, &result);
    assert_eq!(expected.str_field("termination"), Some("target-bracketed"));
    assert_eq!(expected.str_field("objective"), Some("maximize-regional-power-schedule-scale"));
    let selected = expected.get("selected").unwrap();
    let amplitude = selected.f64_field("power_schedule_scale").unwrap();
    assert!((amplitude - 1.0).abs() <= 0.001);
    assert!(selected.f64_field("temperature_k").unwrap() <= expected.f64_field("temperature_limit_k").unwrap());
    assert!(selected.f64_field("final_temperature_k").unwrap() < selected.f64_field("temperature_k").unwrap());
    assert!(expected.f64_field("bracket_width_scale").unwrap() <= 0.001);
    let bytes = artifact(&full.ledger, result.path(&["receipt", "selected_project"]).unwrap().as_str().unwrap());
    let project = fs_project::parse_sexpr(std::str::from_utf8(&bytes).unwrap()).unwrap();
    assert_eq!(project.spec, scaled(&full, amplitude));
    let (hash, run, value) = full.solve_project(950, &project.spec);
    assert_eq!(selected.str_field("project_hash"), Some(hash.as_str()));
    assert_eq!(selected.str_field("run"), Some(run.as_str()));
    assert_eq!(selected.f64_field("final_temperature_k").unwrap().to_bits(), value.to_bits());
    let split = design_fixture();
    let prefix = split.study(Some("1"), fs_cli::exit::BUDGET);
    assert_eq!(report(&split, &prefix).path(&["selected", "power_schedule_scale"]).unwrap().as_f64(), Some(0.0));
    std::fs::rename(&split.sources, split.dir.join("relocated-pulse-sources")).unwrap();
    let finished = split.resume(prefix.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    assert_eq!(report(&split, &finished), expected);
}

#[test]
fn g3_zero_pulse_cannot_repair_initial_overheating() {
    let mut fixture = pulse_fixture();
    fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap().transient.as_mut().unwrap()
        .initial_temperature = TransientInitialTemperature::Uniform(QtyAny::new(330.0, dims::TEMPERATURE));
    write_project(&fixture);
    std::fs::write(fixture.sources.join("study.fsim"), DESIGN).unwrap();
    let result = fixture.study(None, fs_cli::exit::REFUSED);
    let report = report(&fixture, &result);
    assert_eq!(report.str_field("termination"), Some("minimum-schedule-amplitude-infeasible"));
    assert_eq!(report.get("selected"), Some(&JsonValue::Null));
    let initial = &report.get("history").unwrap().as_array().unwrap()[0];
    assert_eq!(initial.f64_field("power_schedule_scale"), Some(0.0));
    assert_eq!(initial.f64_field("time_s"), Some(0.0));
    assert_eq!(initial.f64_field("temperature_k"), Some(330.0));
}
