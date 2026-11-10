//! `frankensim system`: scenes from examples/system run through the real
//! binary, and their reports agree with references computed here.

use std::path::PathBuf;
use std::process::Command;

use fs_system::elements::Waveform;

fn example(name: &str) -> String {
    format!(
        "{}/../../examples/system/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

fn run_binary(args: &[&str]) -> (i32, String, String) {
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(args)
        .output()
        .expect("start binary");
    (
        output.status.code().unwrap_or(-1),
        String::from_utf8(output.stdout).unwrap(),
        String::from_utf8(output.stderr).unwrap(),
    )
}

/// The number after `"key":` that follows `anchor` in a JSON report.
fn value(json: &str, anchor: &str, key: &str) -> f64 {
    let from = json
        .find(anchor)
        .unwrap_or_else(|| panic!("{anchor} not in report"));
    let tail = &json[from..];
    let field = format!("\"{key}\":");
    let start = tail
        .find(&field)
        .unwrap_or_else(|| panic!("{key} not after {anchor}"))
        + field.len();
    let digits: String = tail[start..]
        .chars()
        .take_while(|c| c.is_ascii_digit() || matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .collect();
    digits
        .parse()
        .unwrap_or_else(|_| panic!("{key} after {anchor} is not a number: {digits:?}"))
}

fn element(name: &str) -> String {
    format!("{{\"name\":\"{name}\"")
}

#[test]
fn the_driven_oscillator_matches_the_monolithic_midpoint_rule() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("driven_oscillator.json")]);
    assert_eq!(code, 0, "{stderr}");
    assert!(stderr.is_empty());
    assert!(json.contains("\"schema\":\"frankensim.system.result.v0\""));

    // m x'' = -k x - c x' + F(t), F the exact step mean of the sine.
    let (m, k, c) = (0.25_f64, 400.0_f64, 1.5_f64);
    let force = Waveform::sine(0.8, 28.64788975654116, 0.0, 0.0).unwrap();
    let (dt, steps) = (2.0e-4_f64, 10_000_u32);
    let (mut x, mut v, mut loss, mut peak) = (0.002_f64, -0.3_f64, 0.0_f64, 0.002_f64);
    let half = 0.5 * dt / m;
    for i in 0..steps {
        let f = force.mean(f64::from(i) * dt, dt);
        let mean = (v + half * (f - k * x)) / (1.0 + half * c + 0.5 * dt * half * k);
        x += dt * mean;
        v = 2.0 * mean - v;
        loss += dt * c * mean * mean;
        peak = peak.max(x.abs());
    }
    let counted = value(&json, "\"steps\"", "steps");
    assert!((counted - f64::from(steps)).abs() < 0.5, "{counted} steps");
    let block = element("block");
    let position = value(&json, &block, "position_m");
    let velocity = value(&json, &block, "velocity_m_s");
    assert!((position - x).abs() <= 1.0e-9 * peak, "{position} vs {x}");
    assert!((velocity - v).abs() <= 1.0e-9 * peak * (k / m).sqrt(), "{velocity} vs {v}");
    // The spring's extension is minus the block's position, to the
    // interface tolerance.
    let extension = value(&json, &element("spring"), "extension_m");
    assert!((extension + position).abs() <= 1.0e-12);

    // Every joule the damper absorbed went to the housing, and from there
    // partly to the ambient boundary.
    let damper = element("damper");
    let dissipated = value(&json, &damper, "dissipated_j");
    assert!((dissipated - loss).abs() <= 1.0e-9 * loss, "{dissipated} vs {loss}");
    let housing = value(&json, &element("housing"), "stored_change_j");
    let ambient = value(&json, "\"port\":\"skin.b\"", "energy_in_j");
    assert!(ambient < 0.0, "heat leaves through the skin");
    assert!((housing - ambient - dissipated).abs() <= 1.0e-10 * dissipated);
    let residual = value(&json, "\"totals\"", "energy_residual_j");
    assert!(residual.abs() <= 1.0e-12, "{residual}");
    assert!(value(&json, "\"totals\"", "entropy_generated_j_per_k") > 0.0);
    // Samples every 0.01 s from 0 to 2 s.
    assert_eq!(json.matches("{\"time_s\":").count(), 201);
}

#[test]
fn the_heated_body_approaches_its_steady_temperature_by_the_discrete_factor() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("heated_body.json")]);
    assert_eq!(code, 0, "{stderr}");
    let (c, g, power, ambient, dt, steps) = (500.0_f64, 2.0_f64, 30.0_f64, 293.15_f64, 5.0_f64, 400);
    let steady = ambient + power / g;
    let factor = (1.0 - 0.5 * g * dt / c) / (1.0 + 0.5 * g * dt / c);
    let expected = steady + (ambient - steady) * factor.powi(steps);
    let reached = value(&json, &element("body"), "temperature_k");
    assert!((reached - expected).abs() <= 1.0e-9 * steady, "{reached} vs {expected}");
    // The heater's account is the energy it gave.
    let given = -value(&json, &element("heater"), "energy_in_j");
    assert!((given - power * dt * f64::from(steps)).abs() <= 1.0e-9 * given);
}

#[test]
fn the_quarter_car_returns_the_road_work_as_motion_and_heat_split_seventy_thirty() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("quarter_car.json")]);
    assert_eq!(code, 0, "{stderr}");
    let road = value(&json, "\"port\":\"tyre.a\"", "energy_in_j");
    let shock = element("shock");
    let dissipated = value(&json, &shock, "dissipated_j");
    assert!(dissipated > 0.0);
    let warmed = value(&json, &element("shock_body"), "stored_change_j");
    let air = value(&json, "\"name\":\"air\"", "heat_j");
    assert!((warmed - 0.7 * dissipated).abs() <= 1.0e-12 * dissipated, "{warmed} vs {dissipated}");
    assert!((air - 0.3 * dissipated).abs() <= 1.0e-12 * dissipated, "{air} vs {dissipated}");
    // The road's work is what is left as motion and spring energy plus all
    // the heat.
    let mechanical: f64 = ["body", "wheel", "suspension", "tyre"]
        .iter()
        .map(|name| value(&json, &element(name), "stored_change_j"))
        .sum();
    assert!(
        (road - mechanical - dissipated).abs() <= 1.0e-9 * road.abs().max(dissipated),
        "road {road}, mechanical {mechanical}, heat {dissipated}"
    );
    // The bump is over by 0.7 s and the shock has nearly stopped the body
    // by 3 s.
    assert!(mechanical.abs() < 0.05 * dissipated, "{mechanical} vs {dissipated}");
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-9 * dissipated);
}

#[test]
fn library_and_binary_give_the_same_bytes_and_text_mode_reports_the_ledger() {
    let scene = example("quarter_car.json");
    let library = fs_cli::run(["--json", "system", &scene].map(str::to_owned));
    assert_eq!(library.exit_code, fs_cli::exit::SUCCESS, "{}", library.stderr);
    let (code, json, _) = run_binary(&["--json", "system", &scene]);
    assert_eq!(code, 0);
    assert_eq!(library.stdout, json, "the run is deterministic");
    let (code, text, _) = run_binary(&["system", &scene]);
    assert_eq!(code, 0);
    assert!(text.contains("shock_body"));
    assert!(text.contains("sink air: heat"));
    assert!(text.contains("No claim:"));
    let (code, help, _) = run_binary(&["system", "--help"]);
    assert_eq!(code, 0);
    assert!(help.contains("frankensim.system.v0"));
    let (code, _, _) = run_binary(&["system"]);
    assert_eq!(code, 2);
}

fn scratch(name: &str, text: &str) -> PathBuf {
    let directory = std::env::temp_dir().join(format!("fs-cli-system-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join(name);
    std::fs::write(&path, text).unwrap();
    path
}

fn refused(name: &str, text: &str) -> (i32, String) {
    let path = scratch(name, text);
    let (code, stdout, stderr) = run_binary(&["system", path.to_str().unwrap()]);
    assert!(stdout.is_empty(), "{stdout}");
    (code, stderr)
}

const OSCILLATOR: &str = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.01, "step_s": 0.001,
  "elements": {
    "m": {"kind": "mass", "mass_kg": 1.0, "velocity_m_s": 1.0},
    "s": {"kind": "spring", "stiffness_n_m": 100.0}
  },
  "connections": [
    {"join": ["m.p", "s.a"]},
    {"prescribe": "s.b", "velocity_m_s": 0}
  ]
}"#;

#[test]
fn scenes_that_cannot_run_are_refused_with_the_reason() {
    // The template itself runs.
    let path = scratch("ok.json", OSCILLATOR);
    let (code, _, stderr) = run_binary(&["system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");

    let (code, stderr) = refused("kind.json", &OSCILLATOR.replace("\"spring\"", "\"sprung\""));
    assert_eq!(code, 3);
    assert!(stderr.contains("elements.s: unknown kind \"sprung\""), "{stderr}");

    let (code, stderr) = refused("typo.json", &OSCILLATOR.replace("\"mass_kg\"", "\"mass_g\""));
    assert_eq!(code, 3);
    assert!(stderr.contains("unsupported field \"mass_g\""), "{stderr}");

    let (code, stderr) = refused("schema.json", &OSCILLATOR.replace("system.v0", "system.v9"));
    assert_eq!(code, 3);
    assert!(stderr.contains("schema"), "{stderr}");

    // An end left unconnected is refused by the kernel, by name.
    let loose = OSCILLATOR.replace(",\n    {\"prescribe\": \"s.b\", \"velocity_m_s\": 0}", "");
    let (code, stderr) = refused("loose.json", &loose);
    assert_eq!(code, 4);
    assert!(stderr.contains("s.b"), "{stderr}");

    // A port that does not exist names the ones that do.
    let (code, stderr) = refused("port.json", &OSCILLATOR.replace("\"s.b\"", "\"s.c\""));
    assert_eq!(code, 3);
    assert!(stderr.contains("a spring has ports a, b"), "{stderr}");

    // Two masses cannot share a node: nothing sets it.
    let two = OSCILLATOR.replace(
        "\"s\": {",
        "\"n\": {\"kind\": \"mass\", \"mass_kg\": 1.0},\n    \"s\": {",
    )
    .replace("[\"m.p\", \"s.a\"]", "[\"m.p\", \"n.p\", \"s.a\"]");
    let (code, stderr) = refused("two.json", &two);
    assert_eq!(code, 3);
    assert!(stderr.contains("this one has 2"), "{stderr}");

    // A velocity cannot be prescribed on a mass; its port takes a force.
    let wrong = OSCILLATOR.replace(
        "{\"prescribe\": \"s.b\", \"velocity_m_s\": 0}",
        "{\"prescribe\": \"s.b\", \"force_n\": 0}",
    );
    let (code, stderr) = refused("quantity.json", &wrong);
    assert_eq!(code, 3);
    assert!(stderr.contains("exactly one quantity, velocity_m_s"), "{stderr}");

    // Too many steps for the cap.
    let (code, stderr) = refused("long.json", &OSCILLATOR.replace("\"duration_s\": 0.01", "\"duration_s\": 1.0e6"));
    assert_eq!(code, 3);
    assert!(stderr.contains("steps"), "{stderr}");
}
