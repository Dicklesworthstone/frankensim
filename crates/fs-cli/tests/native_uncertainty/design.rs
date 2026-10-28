//! G1/G3/G4 native power sizing: actual projects, actual QoIs, retained recovery.
mod cases;
use super::*;
use fs_project::{ConductionTransient, TransientInitialTemperature, TransientRegionCapacity, spec::dims};
use fs_qty::QtyAny;

const DESIGN: &str = r#"(fsim-native-design-study
    :version 1 :project "cooling-reference.fsim" :evaluations 32 :seed 29
    :wall-time 240s :qoi "temperature-max"
    :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
    :materials ("aa6061.fsmcdpk") :interfaces ()
    :control (regional-power :entity "air" :minimum 0W :maximum 16W :tolerance 0.05W)
    :temperature-tolerance 0.01K :assumption nondecreasing-temperature-with-power)"#;

fn write(fixture: &Fixture) {
    std::fs::write(fixture.sources.join("cooling-reference.fsim"),
        fs_project::print_sexpr(&fixture.project).unwrap()).unwrap();
    std::fs::write(fixture.sources.join("study.fsim"), DESIGN).unwrap();
}

fn add_transient(fixture: &mut Fixture, initial: f64) {
    let setup = fixture.project.cooling.as_mut().unwrap().conduction.as_mut().unwrap();
    setup.transient = Some(ConductionTransient {
        initial_temperature: TransientInitialTemperature::Uniform(QtyAny::new(initial, dims::TEMPERATURE)),
        horizon: QtyAny::new(2.0, dims::TIME),
        max_step: QtyAny::new(0.5, dims::TIME),
        max_steps: 12,
        energy_tolerance: QtyAny::new(1e-6, dims::ENERGY),
        capacities: setup.regions.iter().map(|row| TransientRegionCapacity {
            region: row.region.clone(),
            volumetric_heat_capacity: QtyAny::new(100.0, dims::VOLUMETRIC_HEAT_CAPACITY),
            source: "synthetic storage for native power design integration".into(),
        }).collect(),
        enthalpy: None,
        power_schedules: Vec::new(),
    });
}

fn fixture(transient: bool) -> Fixture {
    let mut fixture = Fixture::new();
    fixture.project.power.as_mut().unwrap()[0].duty = 0.5;
    fixture.project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
    if transient { add_transient(&mut fixture, 293.15); }
    // Establish the limit with an independently authored ordinary project,
    // not the design search or its coordinate-mutation helper.
    let mut oracle = fixture.project.clone();
    oracle.power.as_mut().unwrap()[0].watts.value = 16.0; // 8 delivered W.
    let (_, _, target) = fixture.solve_project(700, &oracle);
    let requirement = &mut fixture.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = target + 1.0;
    requirement.margin.value = 1.0;
    write(&fixture);
    fixture
}

fn report(fixture: &Fixture, result: &JsonValue) -> JsonValue {
    let receipt = result.get("receipt").unwrap();
    assert_eq!(receipt.str_field("driver"), Some("native-cooling-power-design-v1"));
    assert_eq!(receipt.str_field("status"), result.str_field("status"));
    json_artifact(&fixture.ledger, receipt.str_field("report_json").unwrap())
}

#[test]
fn g1_g3_native_power_design_selects_a_feasible_replayable_project() {
    for transient in [false, true] {
        let fixture = fixture(transient);
        let result = fixture.study(None, fs_cli::exit::SUCCESS);
        let report = report(&fixture, &result);
        assert_eq!(report.str_field("status"), Some("completed"));
        assert_eq!(report.str_field("termination"), Some("target-bracketed"));
        assert_eq!(report.str_field("temperature_observation"),
            Some(if transient { "fine-grid-sampled-peak-including-initial" } else { "steady" }));
        let selected = report.get("selected").unwrap();
        let power = selected.f64_field("power_w").unwrap();
        let temperature = selected.f64_field("temperature_k").unwrap();
        let limit = report.f64_field("temperature_limit_k").unwrap();
        assert!(temperature <= limit);
        assert!((power - 8.0).abs() <= 0.05, "selected {power} W");
        assert!(report.f64_field("bracket_width_w").unwrap() <= 0.05);
        let bytes = artifact(&fixture.ledger,
            result.path(&["receipt", "selected_project"]).unwrap().as_str().unwrap());
        let project = fs_project::parse_sexpr(std::str::from_utf8(&bytes).unwrap()).unwrap();
        let authored = &project.spec.power.as_ref().unwrap()[0];
        assert_eq!(authored.duty, 0.5);
        assert_eq!(authored.watts.value * authored.duty, power);
        let (hash, run, value) = fixture.solve_project(701, &project.spec);
        assert_eq!(selected.str_field("project_hash"), Some(hash.as_str()));
        assert_eq!(selected.str_field("run"), Some(run.as_str()));
        assert_eq!(selected.f64_field("final_temperature_k").unwrap().to_bits(), value.to_bits());
        if transient {
            assert!(selected.f64_field("time_s").unwrap() <= 2.0);
            assert!(selected.f64_field("state").unwrap() <= 8.0);
            assert!(temperature >= value);
        } else {
            assert_eq!(temperature.to_bits(), value.to_bits());
        }
        let exported = command(&["--json", "report", result.str_field("run").unwrap(),
            fixture.ledger.to_str().unwrap()], fs_cli::exit::SUCCESS);
        let path = exported.str_field("selected_project").expect("selected native project exported");
        assert_eq!(std::fs::read(path).unwrap(), bytes);
    }
}

#[test]
fn g4_g5_power_design_retains_baseline_and_resumes_without_original_sources() {
    let full = fixture(false);
    let complete = full.study(None, fs_cli::exit::SUCCESS);
    let expected = report(&full, &complete);
    let split = fixture(false);
    let first = split.study(Some("2"), fs_cli::exit::BUDGET);
    let partial = report(&split, &first);
    assert_eq!(partial.str_field("status"), Some("budget-exhausted"));
    assert_eq!(partial.str_field("termination"), Some("invocation-evaluation-budget"));
    assert_eq!(partial.f64_field("evaluations_completed"), Some(2.0));
    assert_eq!(partial.path(&["selected", "power_w"]).unwrap().as_f64(), Some(0.0));
    assert!(first.path(&["receipt", "selected_project"]).unwrap().as_str().is_some());
    std::fs::rename(&split.sources, split.dir.join("relocated-inputs")).unwrap();
    let finished = split.resume(first.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    let actual = report(&split, &finished);
    assert_eq!(actual, expected, "resumption must reproduce actual candidate history and selected project");
}

#[test]
fn g1_transient_design_refuses_initial_overheating_even_when_the_final_field_passes() {
    let mut fixture = Fixture::new();
    fixture.project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
    add_transient(&mut fixture, 330.0);
    let requirement = &mut fixture.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = 321.0;
    requirement.margin.value = 1.0;
    write(&fixture);
    let result = fixture.study(None, fs_cli::exit::REFUSED);
    let report = report(&fixture, &result);
    assert_eq!(report.str_field("termination"), Some("minimum-power-infeasible"));
    assert_eq!(report.get("selected"), Some(&JsonValue::Null));
    assert_eq!(result.path(&["receipt", "selected_project"]), Some(&JsonValue::Null));
    let row = &report.get("history").unwrap().as_array().unwrap()[0];
    assert_eq!(row.f64_field("temperature_k"), Some(330.0));
    assert_eq!(row.f64_field("time_s"), Some(0.0));
    assert_eq!(row.f64_field("state"), Some(0.0));
    assert!(row.f64_field("final_temperature_k").unwrap() < 320.0);
}

fn fan_source() -> String {
    DESIGN.replace("regional-power :entity \"air\" :minimum 0W :maximum 16W :tolerance 0.05W",
        "fan-speed-ratio :entity \"fixture-bank\" :minimum 0.8 :maximum 1.2 :tolerance 0.001")
        .replace("nondecreasing-temperature-with-power", "nonincreasing-temperature-with-fan-speed")
}

fn fan_fixture(transient: bool) -> Fixture {
    let mut fixture = super::fan_speed::fixture();
    fixture.project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
    if transient { add_transient(&mut fixture, 294.0); }
    // Independently set the absolute bank ratio rather than using the new
    // design coordinate binder or its bracketing implementation.
    let mut oracle = fixture.project.clone();
    oracle.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks[0].speed_ratio = 1.0;
    let (_, _, target) = fixture.solve_project(800, &oracle);
    let requirement = &mut fixture.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = target + 1.0;
    requirement.margin.value = 1.0;
    write(&fixture);
    std::fs::write(fixture.sources.join("study.fsim"), fan_source()).unwrap();
    fixture
}

fn fan_report(fixture: &Fixture, result: &JsonValue) -> JsonValue {
    let receipt = result.get("receipt").unwrap();
    assert_eq!(receipt.str_field("driver"), Some("native-cooling-fan-design-v1"));
    let report = json_artifact(&fixture.ledger, receipt.str_field("report_json").unwrap());
    assert_eq!(report.str_field("objective"), Some("minimize-fan-speed-ratio"));
    assert_eq!(report.str_field("fan_bank"), Some("fixture-bank"));
    report
}

#[test]
fn g1_g3_native_fan_design_selects_the_minimum_feasible_speed_and_exports_its_actual_project() {
    for transient in [false, true] {
        let fixture = fan_fixture(transient);
        let result = fixture.study(None, fs_cli::exit::SUCCESS);
        let report = fan_report(&fixture, &result);
        assert_eq!(report.str_field("termination"), Some("target-bracketed"));
        assert_eq!(report.str_field("temperature_observation"),
            Some(if transient { "fine-grid-sampled-peak-including-initial" } else { "steady" }));
        let selected = report.get("selected").unwrap();
        let speed = selected.f64_field("fan_speed_ratio").unwrap();
        assert!((1.0 - 1e-12..=1.001).contains(&speed), "selected absolute ratio {speed}");
        assert!(report.f64_field("bracket_width_ratio").unwrap() <= 0.001);
        let limit = report.f64_field("temperature_limit_k").unwrap();
        assert!(selected.f64_field("temperature_k").unwrap() <= limit);
        let history = report.get("history").unwrap().as_array().unwrap();
        assert_eq!(history[0].f64_field("fan_speed_ratio"), Some(1.2), "start with a measured feasible baseline");
        assert_eq!(history[1].f64_field("fan_speed_ratio"), Some(0.8));
        for row in history {
            if row.f64_field("fan_speed_ratio").unwrap() < speed {
                assert!(row.f64_field("temperature_k").unwrap() > limit,
                    "the selected project must be the lowest-speed observed feasible candidate");
            }
        }
        let bytes = artifact(&fixture.ledger,
            result.path(&["receipt", "selected_project"]).unwrap().as_str().unwrap());
        let exported = fs_project::parse_sexpr(std::str::from_utf8(&bytes).unwrap()).unwrap();
        let mut expected = fixture.project.clone();
        expected.cooling.as_mut().unwrap().fan_system.as_mut().unwrap().banks[0].speed_ratio = speed;
        assert_eq!(exported.spec, expected, "only the chosen bank ratio may change");
        let (hash, run, final_temperature) = fixture.solve_project(801, &expected);
        assert_eq!(selected.str_field("project_hash"), Some(hash.as_str()));
        assert_eq!(selected.str_field("run"), Some(run.as_str()));
        assert_eq!(selected.f64_field("final_temperature_k").unwrap().to_bits(), final_temperature.to_bits());
        assert!(selected.f64_field("temperature_k").unwrap() >= final_temperature);
        let paths = command(&["--json", "report", result.str_field("run").unwrap(),
            fixture.ledger.to_str().unwrap()], fs_cli::exit::SUCCESS);
        assert_eq!(std::fs::read(paths.str_field("selected_project").unwrap()).unwrap(), bytes);
    }
}

#[test]
fn g4_g5_native_fan_design_retains_safe_speed_and_resumes_from_original_sources() {
    let whole = fan_fixture(false);
    let complete = whole.study(None, fs_cli::exit::SUCCESS);
    let expected = fan_report(&whole, &complete);
    let split = fan_fixture(false);
    let partial = split.study(Some("2"), fs_cli::exit::BUDGET);
    let before = fan_report(&split, &partial);
    assert_eq!(before.str_field("termination"), Some("invocation-evaluation-budget"));
    assert_eq!(before.path(&["selected", "fan_speed_ratio"]).unwrap().as_f64(), Some(1.2));
    std::fs::rename(&split.sources, split.dir.join("relocated-fan-design-inputs")).unwrap();
    let result = split.resume(partial.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    assert_eq!(fan_report(&split, &result), expected);
    // Terminal recovery must route the fan producer without returning a
    // probability-study authority label or spending more evaluations.
    let again = split.resume(result.str_field("run").unwrap(), None, fs_cli::exit::SUCCESS);
    assert_eq!(again, result);
}

#[test]
fn g1_native_fan_design_refuses_initial_overheating_at_maximum_speed() {
    let mut fixture = super::fan_speed::fixture();
    fixture.project.solver.as_mut().unwrap().tolerance_rel = 1e-10;
    add_transient(&mut fixture, 330.0);
    let requirement = &mut fixture.project.requirements.as_mut().unwrap()[0];
    requirement.limit.value = 321.0;
    requirement.margin.value = 1.0;
    write(&fixture);
    std::fs::write(fixture.sources.join("study.fsim"), fan_source()).unwrap();
    let result = fixture.study(None, fs_cli::exit::REFUSED);
    let report = fan_report(&fixture, &result);
    assert_eq!(report.str_field("termination"), Some("maximum-fan-speed-infeasible"));
    assert_eq!(report.get("selected"), Some(&JsonValue::Null));
    let row = &report.get("history").unwrap().as_array().unwrap()[0];
    assert_eq!(row.f64_field("fan_speed_ratio"), Some(1.2));
    assert_eq!(row.f64_field("temperature_k"), Some(330.0));
    assert_eq!(row.f64_field("time_s"), Some(0.0));
    assert!(row.f64_field("final_temperature_k").unwrap() < 320.0);
}
