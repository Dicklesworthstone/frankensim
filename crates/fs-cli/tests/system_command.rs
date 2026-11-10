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
fn the_dc_motor_scene_matches_the_monolithic_midpoint_rule_across_three_domains() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("dc_motor.json")]);
    assert_eq!(code, 0, "{stderr}");

    // L di/dt = V - R i - k w,  J dw/dt = k i - b w, by the midpoint rule:
    // two linear equations for the step means of i and w.
    let (volts, l, r, k, j, b) = (12.0_f64, 2.0e-3_f64, 1.2_f64, 0.08_f64, 3.0e-4_f64, 2.0e-4_f64);
    let (dt, steps) = (2.0e-5_f64, 20_000_u32);
    let (a11, a12, a21, a22) = (2.0 * l / dt + r, k, -k, 2.0 * j / dt + b);
    let determinant = a11 * a22 - a12 * a21;
    let (mut current, mut speed) = (0.0_f64, 0.0_f64);
    let (mut supplied, mut heat) = (0.0_f64, 0.0_f64);
    for _ in 0..steps {
        let (b1, b2) = (volts + 2.0 * l * current / dt, 2.0 * j * speed / dt);
        let mean_i = (b1 * a22 - a12 * b2) / determinant;
        let mean_w = (a11 * b2 - a21 * b1) / determinant;
        current = 2.0 * mean_i - current;
        speed = 2.0 * mean_w - speed;
        supplied += dt * volts * mean_i;
        heat += dt * (r * mean_i * mean_i + b * mean_w * mean_w);
    }
    let reached_i = value(&json, &element("winding_l"), "current_a");
    let reached_w = value(&json, &element("rotor"), "angular_velocity_rad_s");
    assert!((reached_i - current).abs() <= 1.0e-9 * volts / r, "{reached_i} vs {current}");
    assert!((reached_w - speed).abs() <= 1.0e-9 * speed, "{reached_w} vs {speed}");
    // Near k V / (k^2 + R b) after many time constants.
    let steady = k * volts / (k * k + r * b);
    assert!((reached_w - steady).abs() <= 0.02 * steady, "{reached_w} vs {steady}");

    // The supply's account, both heats in the one case, and a gyrator that
    // keeps nothing.
    let delivered = -value(&json, &element("supply"), "energy_in_j");
    assert!((delivered - supplied).abs() <= 1.0e-9 * supplied, "{delivered} vs {supplied}");
    let warmed = value(&json, &element("case"), "stored_change_j");
    assert!((warmed - heat).abs() <= 1.0e-9 * supplied, "{warmed} vs {heat}");
    let converted = value(&json, &element("machine"), "port_energy_in_j");
    assert!(converted.abs() <= 1.0e-12 * supplied, "{converted}");
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-11 * supplied);

    // A rotor port joined into the current loop is refused by the kernel's
    // port kinds, by name.
    let scene = std::fs::read_to_string(example("dc_motor.json")).unwrap();
    let crossed = scene
        .replace("\"machine.a\"]", "\"machine.a\", \"bearing.b\"]")
        .replace("    {\"prescribe\": \"bearing.b\", \"angular_velocity_rad_s\": 0},\n", "");
    assert_ne!(crossed, scene);
    let (code, stderr) = refused("motor_crossed.json", &crossed);
    assert_eq!(code, 4);
    assert!(stderr.contains("bearing.b"), "{stderr}");
    // A current cannot be prescribed on a shaft, and a gyrator names its sides.
    let (code, stderr) = refused(
        "motor_quantity.json",
        &scene.replace("\"angular_velocity_rad_s\": 0", "\"current_a\": 0"),
    );
    assert_eq!(code, 3);
    assert!(stderr.contains("exactly one quantity, angular_velocity_rad_s"), "{stderr}");
    let (code, stderr) = refused("motor_side.json", &scene.replace("\"b\": \"rotation\"", "\"b\": \"shaft\""));
    assert_eq!(code, 3);
    assert!(stderr.contains("b is one of translation, rotation, electrical"), "{stderr}");
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

/// The midpoint rule for `I f' = k e - r f + effort`, `e' = -f`: the flow of
/// a node and the extension of the storage element anchored beside it.
fn loop_reference(
    (inertance, stiffness, resistance, effort): (f64, f64, f64, f64),
    extension: f64,
    dt: f64,
    steps: u32,
) -> (f64, f64) {
    let (mut flow, mut extension) = (0.0_f64, extension);
    let pivot = 2.0 * inertance / dt + 0.5 * stiffness * dt + resistance;
    for _ in 0..steps {
        let mean = (2.0 * inertance * flow / dt + stiffness * extension + effort) / pivot;
        flow = 2.0 * mean - flow;
        extension -= dt * mean;
    }
    (flow, extension)
}

#[test]
fn an_rlc_loop_and_a_driven_torsion_pendulum_follow_the_same_law_in_their_own_units() {
    // A charged capacitor rings down through an inductor and a resistor.
    let (henry, farad, ohm, coulomb) = (10.0e-3_f64, 100.0e-6_f64, 2.0_f64, 1.0e-3_f64);
    let scene = format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 0.02, "step_s": 1.0e-5,
  "elements": {{
    "coil": {{"kind": "inductor", "inductance_h": {henry}}},
    "cap": {{"kind": "capacitor", "capacitance_f": {farad}, "charge_c": {coulomb}}},
    "load": {{"kind": "resistor", "resistance_ohm": {ohm}}},
    "air": {{"kind": "heat_capacity", "capacity_j_k": 1.0, "temperature_k": 300.0}}
  }},
  "connections": [
    {{"join": ["coil.p", "cap.a", "load.a"]}},
    {{"prescribe": "cap.b", "current_a": 0}},
    {{"prescribe": "load.b", "current_a": 0}},
    {{"route": "load.heat", "to": [["air.heat", 1.0]]}}
  ]
}}"#
    );
    let path = scratch("rlc.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (current, charge) = loop_reference((henry, 1.0 / farad, ohm, 0.0), coulomb, 1.0e-5, 2_000);
    let peak = coulomb / (henry * farad).sqrt();
    let reached = value(&json, &element("coil"), "current_a");
    assert!((reached - current).abs() <= 1.0e-9 * peak, "{reached} vs {current}");
    let left = value(&json, &element("cap"), "charge_c");
    assert!((left - charge).abs() <= 1.0e-9 * coulomb, "{left} vs {charge}");
    // The energy the loop lost is the heat in the air.
    let start = 0.5 * coulomb * coulomb / farad;
    let now = 0.5 * left * left / farad + 0.5 * henry * reached * reached;
    let warmed = value(&json, &element("air"), "stored_change_j");
    assert!(warmed > 0.5 * start, "most of the charge's energy is heat by now");
    assert!((start - now - warmed).abs() <= 1.0e-10 * start, "{start} - {now} vs {warmed}");

    // A shaft on a torsion spring, turned by a constant torque.
    let (inertia, kappa, torque) = (2.0e-3_f64, 5.0_f64, 0.4_f64);
    let scene = format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 0.5, "step_s": 1.0e-4,
  "elements": {{
    "wheel": {{"kind": "inertia", "inertia_kg_m2": {inertia}}},
    "wire": {{"kind": "torsion_spring", "stiffness_n_m_rad": {kappa}}},
    "motor": {{"kind": "torque", "torque_n_m": {torque}}}
  }},
  "connections": [
    {{"join": ["wheel.p", "wire.a", "motor.p"]}},
    {{"prescribe": "wire.b", "angular_velocity_rad_s": 0}}
  ]
}}"#
    );
    let path = scratch("pendulum.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (speed, twist) = loop_reference((inertia, kappa, 0.0, torque), 0.0, 1.0e-4, 5_000);
    let swing = torque / kappa;
    let reached = value(&json, &element("wheel"), "angular_velocity_rad_s");
    assert!((reached - speed).abs() <= 1.0e-9 * swing * (kappa / inertia).sqrt(), "{reached} vs {speed}");
    let wound = value(&json, &element("wire"), "twist_rad");
    assert!((wound - twist).abs() <= 1.0e-9 * swing, "{wound} vs {twist}");
    // The wheel's angle is minus the wire's twist, and it swings about
    // torque / stiffness without loss.
    let angle = value(&json, &element("wheel"), "angle_rad");
    assert!((angle + wound).abs() <= 1.0e-12, "{angle} vs {wound}");
    assert!(angle > 0.0 && angle < 2.0 * swing * (1.0 + 1.0e-9), "{angle} vs {swing}");
    let work = -value(&json, &element("motor"), "energy_in_j");
    let stored = 0.5 * kappa * wound * wound + 0.5 * inertia * reached * reached;
    assert!((work - stored).abs() <= 1.0e-10 * torque * swing, "{work} vs {stored}");
}

#[test]
fn a_prescribed_end_moves_the_way_it_is_written() {
    // The far end of a spring is pulled at +1 m/s: the spring stretches and
    // the mass follows in the same direction.
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.1, "step_s": 1.0e-3,
  "elements": {
    "m": {"kind": "mass", "mass_kg": 1.0},
    "s": {"kind": "spring", "stiffness_n_m": 100.0}
  },
  "connections": [
    {"join": ["m.p", "s.a"]},
    {"prescribe": "s.b", "velocity_m_s": 1.0}
  ]
}"#;
    let path = scratch("pulled.json", scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let position = value(&json, &element("m"), "position_m");
    let velocity = value(&json, &element("m"), "velocity_m_s");
    let extension = value(&json, &element("s"), "extension_m");
    assert!(position > 0.0 && velocity > 0.0, "{position}, {velocity}");
    // The end has moved 0.1 m; the spring holds what the mass has not yet covered.
    assert!((extension - (0.1 - position)).abs() <= 1.0e-12, "{extension} vs {position}");
    // The pull did positive work, all of it now in the spring and the mass.
    let work = value(&json, "\"port\":\"s.b\"", "energy_in_j");
    let stored = 50.0 * extension * extension + 0.5 * velocity * velocity;
    assert!(work > 0.0 && (work - stored).abs() <= 1.0e-12, "{work} vs {stored}");
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
