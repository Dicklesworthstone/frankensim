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
fn adaptive_oscillator_tracks_the_analytic_phase_more_accurately_with_tighter_tolerances() {
    let scene = std::fs::read_to_string(example("adaptive_oscillator.json")).unwrap();
    let mut errors = Vec::new();
    let mut accepted_steps = Vec::new();
    for tight in [false, true] {
        let text = if tight {
            scene.replace("\"absolute\": 1e-5", "\"absolute\": 1e-7")
                .replace("\"absolute\": 1e-6", "\"absolute\": 1e-8")
        } else { scene.clone() };
        let path = scratch(if tight { "adaptive_tight.json" } else { "adaptive_loose.json" }, &text);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        let time = value(&json, "\"schema\"", "time_s");
        assert!((time - 4.0).abs() < 1.0e-11);
        // x'' + 4x = 0.4, x(0)=0.2, v(0)=0: the equilibrium is 0.1 m
        // and the natural frequency is 2 rad/s, independently of stepping.
        let expected_x = 0.1 + 0.1 * (2.0 * time).cos();
        let expected_v = -0.2 * (2.0 * time).sin();
        let x = value(&json, &element("block"), "position_m");
        let v = value(&json, &element("block"), "velocity_m_s");
        errors.push((x - expected_x).hypot(0.5 * (v - expected_v)));
        let steps = value(&json, "\"schema\"", "steps");
        accepted_steps.push(steps);
        assert!(value(&json, "\"adaptive\"", "attempts") > 0.5 * steps);
        assert!(value(&json, "\"adaptive\"", "trial_steps") > steps);
        assert!(value(&json, "\"adaptive\"", "max_error_norm") <= 1.0);
        assert!(value(&json, "\"adaptive\"", "last_error_norm") <= 1.0);

        // Trial work must not survive rejection. Only the actual accepted
        // displacement enters the source's cumulative work account.
        let source_energy = value(&json, &element("drive"), "energy_in_j");
        assert!((source_energy + 0.4 * (x - 0.2)).abs() < 1.0e-11);
        assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-10);
    }
    assert!(errors[1] < 0.15 * errors[0], "loose {}, tight {}", errors[0], errors[1]);
    assert!(errors[1] < 5.0e-5, "tight error {}", errors[1]);
    assert!(accepted_steps[1] > accepted_steps[0]);
}

#[test]
fn an_exhausted_adaptive_trial_reports_the_original_state_time_after_rollback() {
    let scene = std::fs::read_to_string(example("adaptive_oscillator.json")).unwrap()
        .replace("\"max_attempts\": 32", "\"max_attempts\": 1")
        .replace("\"absolute\": 1e-5", "\"absolute\": 1e-20")
        .replace("\"absolute\": 1e-6", "\"absolute\": 1e-20");
    let (code, stderr) = refused("adaptive_attempt_exhaustion.json", &scene);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("system-adaptive"), "{stderr}");
    assert!(stderr.contains("step 1 at t = 0 s"), "{stderr}");
}

#[test]
fn adaptive_tolerances_are_complete_and_a_short_final_horizon_can_cross_the_step_floor() {
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.03125, "step_s": 0.125,
  "adaptive": {"min_step_s": 0.0625, "tolerances": {
    "cap.charge_c": {"absolute": 1e-12, "relative": 0}
  }},
  "elements": {"cap": {"kind": "node_capacitor", "capacitance_f": 2.0, "charge_c": 0.25}},
  "connections": [{"prescribe": "cap.p", "current_a": 0.5}]
}"#;
    let path = scratch("adaptive_final_interval.json", scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    assert!((value(&json, "\"schema\"", "time_s") - 0.03125).abs() < 1.0e-15);
    assert!((value(&json, &element("cap"), "charge_c") - 0.265625).abs() < 1.0e-15);
    assert!((value(&json, "\"schema\"", "steps") - 2.0).abs() < 0.5);
    for (name, changed, expected) in [
        ("adaptive_missing_account.json", scene.replace("\"cap.charge_c\": {\"absolute\": 1e-12, \"relative\": 0}", ""), "cap.charge_c"),
        ("adaptive_derived_account.json", scene.replace("\"cap.charge_c\"", "\"cap.voltage_v\""), "cap.voltage_v"),
        ("adaptive_zero_tolerance.json", scene.replace("\"absolute\": 1e-12", "\"absolute\": 0"), "absolute"),
    ] {
        let (code, stderr) = refused(name, &changed);
        assert_eq!(code, 3, "{stderr}");
        assert!(stderr.contains(expected), "{stderr}");
    }
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

    // The assertions above follow from the kernel's audit and the route
    // fractions for any dynamics. This pins the dynamics: the midpoint rule
    // for the body (1) and the wheel (2) over the same road,
    //   m1 v1' = -k1 (x1 - x2) - c (v1 - v2),
    //   m2 v2' =  k1 (x1 - x2) + c (v1 - v2) - k2 (x2 - r).
    let (m1, m2, k1, c, k2) = (300.0_f64, 40.0_f64, 2.0e4_f64, 1500.0_f64, 1.8e5_f64);
    let road = Waveform::table(vec![(0.5, 0.0), (0.55, 0.6), (0.6, 0.0), (0.65, -0.6), (0.7, 0.0)])
        .unwrap();
    let (dt, steps) = (5.0e-4_f64, 6_000_u32);
    let (alpha, beta) = (0.5 * k1 * dt + c, 0.5 * k2 * dt);
    let (a11, a12, a22) = (2.0 * m1 / dt + alpha, -alpha, 2.0 * m2 / dt + alpha + beta);
    let determinant = a11 * a22 - a12 * a12;
    let (mut x1, mut v1, mut x2, mut v2, mut r) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let (mut heat, mut crest) = (0.0_f64, 0.0_f64);
    for i in 0..steps {
        let u = road.mean(f64::from(i) * dt, dt);
        let b1 = 2.0 * m1 * v1 / dt - k1 * (x1 - x2);
        let b2 = 2.0 * m2 * v2 / dt + k1 * (x1 - x2) - k2 * (x2 - r) + beta * u;
        let mean1 = (b1 * a22 - a12 * b2) / determinant;
        let mean2 = (a11 * b2 - a12 * b1) / determinant;
        x1 += dt * mean1;
        x2 += dt * mean2;
        v1 = 2.0 * mean1 - v1;
        v2 = 2.0 * mean2 - v2;
        r += dt * u;
        heat += dt * c * (mean1 - mean2) * (mean1 - mean2);
        crest = crest.max(r);
    }
    // The road rose 30 mm and came back.
    assert!((crest - 0.03).abs() <= 1.0e-12 && r.abs() <= 1.0e-12, "{crest}, {r}");
    let close = |name: &str, key: &str, expected: f64, scale: f64| {
        let actual = value(&json, &element(name), key);
        assert!((actual - expected).abs() <= 1.0e-9 * scale, "{name}.{key}: {actual} vs {expected}");
    };
    close("body", "position_m", x1, crest);
    close("wheel", "position_m", x2, crest);
    close("body", "velocity_m_s", v1, 0.6);
    close("wheel", "velocity_m_s", v2, 0.6);
    close("suspension", "extension_m", x1 - x2, crest);
    close("tyre", "extension_m", x2 - r, crest);
    assert!((dissipated - heat).abs() <= 1.0e-9 * heat, "{dissipated} vs {heat}");
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
fn the_parallel_rlc_scene_matches_node_equations_source_work_and_joule_heat() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("parallel_rlc.json")]);
    assert_eq!(code, 0, "{stderr}");
    assert!(stderr.is_empty());

    // C V' = J(t) - G V - I, L I' = V. The source waveform is -J:
    // positive flow enters each owner, including the prescribed source.
    let (c, l, g, v0, i0) = (0.01_f64, 0.2_f64, 0.05_f64, 2.0_f64, -0.1_f64);
    let (dt, steps) = (0.001_f64, 2_000_u32);
    let (mut voltage, mut current) = (v0, i0);
    let (mut supplied_j, mut heat_j) = (0.0_f64, 0.0_f64);
    for step in 0..steps {
        let omega = 6.0 * std::f64::consts::PI;
        let phase = omega * f64::from(step) * dt + 0.3;
        // Integrate the sine independently of the library's waveform helper.
        let supplied = 0.4 + 0.2 * (phase.cos() - (phase + omega * dt).cos()) / (omega * dt);
        let half = 0.5 * dt;
        let mean_v = (voltage + half / c * (supplied - current))
            / (1.0 + half * g / c + half * half / (c * l));
        voltage = 2.0 * mean_v - voltage;
        current += dt / l * mean_v;
        supplied_j += dt * supplied * mean_v;
        heat_j += dt * g * mean_v * mean_v;
    }
    let capacitor = element("capacitor");
    let reached_v = value(&json, &capacitor, "voltage_v");
    let reached_q = value(&json, &capacitor, "charge_c");
    let reached_i = value(&json, &element("coil"), "current_a");
    assert!((reached_v - voltage).abs() <= 1.0e-9, "{reached_v} vs {voltage}");
    assert!((reached_q - c * voltage).abs() <= 1.0e-11, "{reached_q} vs {}", c * voltage);
    assert!((reached_i - current).abs() <= 1.0e-10, "{reached_i} vs {current}");
    // The reported voltage is the endpoint q/C, not the last step's mean.
    assert!((reached_v - reached_q / c).abs() <= 1.0e-13);

    let delivered = -value(&json, &element("source"), "energy_in_j");
    let dissipated = value(&json, &element("load"), "dissipated_j");
    let warmed = value(&json, &element("case"), "stored_change_j");
    let temperature = value(&json, &element("case"), "temperature_k");
    assert!((delivered - supplied_j).abs() <= 1.0e-10, "{delivered} vs {supplied_j}");
    assert!((dissipated - heat_j).abs() <= 1.0e-10, "{dissipated} vs {heat_j}");
    assert!((warmed - heat_j).abs() <= 1.0e-10, "{warmed} vs {heat_j}");
    assert!((10.0 * (temperature - 293.15) - heat_j).abs() <= 2.0e-9);
    let stored_change = 0.5 * (c * (reached_v * reached_v - v0 * v0)
        + l * (reached_i * reached_i - i0 * i0));
    assert!((delivered - stored_change - warmed).abs() <= 1.0e-10);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-10);
    assert!(value(&json, "\"totals\"", "entropy_generated_j_per_k") > 0.0);
    assert_eq!(json.matches("{\"time_s\":").count(), 101);

    let scene = std::fs::read_to_string(example("parallel_rlc.json")).unwrap();
    for (file, original, replacement) in [
        ("parallel_bad_capacitance.json", "\"capacitance_f\": 0.01", "\"capacitance_f\": 0"),
        ("parallel_bad_resistance.json", "\"resistance_ohm\": 20.0", "\"resistance_ohm\": 0"),
    ] {
        let changed = scene.replace(original, replacement);
        assert_ne!(changed, scene);
        let (code, stderr) = refused(file, &changed);
        assert_eq!(code, 3, "{stderr}");
        assert!(stderr.contains("positive"), "{stderr}");
    }
    let wrong_ground = scene.replace("\"voltage_v\": 0", "\"current_a\": 0");
    let (code, stderr) = refused("parallel_ground_quantity.json", &wrong_ground);
    assert_eq!(code, 3, "{stderr}");
    assert!(stderr.contains("exactly one quantity, voltage_v"), "{stderr}");
}

#[test]
fn a_prescribed_current_enters_the_node_capacitor() {
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.1, "step_s": 0.001,
  "elements": {
    "cap": {"kind": "node_capacitor", "capacitance_f": 0.2, "charge_c": 0.03}
  },
  "connections": [{"prescribe": "cap.p", "current_a": 0.4}]
}"#;
    let path = scratch("charged_node.json", scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let charge = value(&json, &element("cap"), "charge_c");
    let voltage = value(&json, &element("cap"), "voltage_v");
    assert!((charge - 0.07).abs() <= 1.0e-12, "{charge}");
    assert!((voltage - 0.35).abs() <= 1.0e-12, "{voltage}");
    let entered = value(&json, "\"port\":\"cap.p\"", "energy_in_j");
    let expected = (0.07_f64.powi(2) - 0.03_f64.powi(2)) / (2.0 * 0.2);
    assert!((entered - expected).abs() <= 1.0e-12, "{entered} vs {expected}");
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-12);
}

// Independent piecewise laws for nonlinear_rlc.json: the outer slope is
// twice the inner slope, and energy integrates the physical conjugate value.
fn nonlinear_effort(coordinate: f64, knee: f64, at_knee: f64) -> f64 {
    let x = coordinate / knee;
    at_knee * if x.abs() <= 1.0 { x } else { 2.0 * x - x.signum() }
}

fn nonlinear_energy(coordinate: f64, knee: f64, at_knee: f64) -> f64 {
    let x = (coordinate / knee).abs();
    knee * at_knee * if x <= 1.0 { 0.5 * x * x } else { x * x - x + 0.5 }
}

fn nonlinear_rhs([charge, flux, _, _]: [f64; 4]) -> [f64; 4] {
    let voltage = nonlinear_effort(charge, 0.02, 2.0);
    let current = nonlinear_effort(flux, 0.05, 1.0);
    // q' = supplied current - resistor current - inductor current;
    // flux' = node voltage. The other states integrate heat and source work.
    [0.05 - 0.1 * voltage - current, voltage, 0.1 * voltage * voltage, 0.05 * voltage]
}

#[test]
fn nonlinear_storage_scene_follows_the_differential_laws_and_closes_all_energy_accounts() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("nonlinear_rlc.json")]);
    assert_eq!(code, 0, "{stderr}");
    let (mut reference, dt) = ([0.032, -0.04, 0.0, 0.0], 5.0e-6);
    // Independent continuous-time RK4, 100x finer than the scene. This
    // reference does not call the native curves or their discrete gradients.
    for _ in 0..100_000 {
        let a = nonlinear_rhs(reference);
        let b = nonlinear_rhs(std::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]));
        let c = nonlinear_rhs(std::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]));
        let d = nonlinear_rhs(std::array::from_fn(|i| reference[i] + dt * c[i]));
        reference = std::array::from_fn(|i| reference[i] + dt * (a[i] + 2.0*b[i] + 2.0*c[i] + d[i]) / 6.0);
    }
    let charge = value(&json, &element("capacitor"), "charge_c");
    let flux = value(&json, &element("coil"), "flux_wb");
    assert!((charge - reference[0]).abs() <= 3.0e-6, "{charge} vs {}", reference[0]);
    assert!((flux - reference[1]).abs() <= 8.0e-6, "{flux} vs {}", reference[1]);
    // Probes report endpoint efforts, not the last step's integral means.
    let voltage = value(&json, &element("capacitor"), "voltage_v");
    let current = value(&json, &element("coil"), "current_a");
    assert!((voltage - nonlinear_effort(charge, 0.02, 2.0)).abs() < 1.0e-13);
    assert!((current - nonlinear_effort(flux, 0.05, 1.0)).abs() < 1.0e-13);

    let heat = value(&json, &element("load"), "dissipated_j");
    let supplied = -value(&json, &element("source"), "energy_in_j");
    assert!((heat - reference[2]).abs() < 1.0e-6);
    assert!((supplied - reference[3]).abs() < 1.0e-6);
    let electrical_change = nonlinear_energy(charge, 0.02, 2.0)
        + nonlinear_energy(flux, 0.05, 1.0)
        - nonlinear_energy(0.032, 0.02, 2.0)
        - nonlinear_energy(-0.04, 0.05, 1.0);
    assert!((electrical_change + heat - supplied).abs() < 1.0e-10);
    let warmed = value(&json, &element("case"), "stored_change_j");
    assert!((warmed - heat).abs() < 1.0e-10);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-10);
    assert!(value(&json, "\"totals\"", "entropy_generated_j_per_k") > 0.0);
    assert_eq!(json.matches("{\"time_s\":").count(), 101);
}

#[test]
fn prescribed_current_and_voltage_advance_the_original_charge_and_flux_curves() {
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.3, "step_s": 0.03,
  "elements": {
    "cap": {"kind": "nonlinear_capacitor", "charge_voltage_curve": [[-0.04, -6.0], [-0.02, -2.0], [0.0, 0.0], [0.02, 2.0], [0.04, 6.0]], "charge_c": 0.012},
    "coil": {"kind": "nonlinear_inductor", "flux_current_curve": [[-0.1, -3.0], [-0.05, -1.0], [0.0, 0.0], [0.05, 1.0], [0.1, 3.0]], "flux_wb": 0.03}
  },
  "connections": [{"prescribe": "cap.p", "current_a": 0.08}, {"prescribe": "coil.p", "voltage_v": 0.2}]
}"#;
    let path = scratch("prescribed_nonlinear_storage.json", scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    assert!((value(&json, &element("cap"), "charge_c") - 0.036).abs() < 1.0e-14);
    assert!((value(&json, &element("cap"), "voltage_v") - 5.2).abs() < 1.0e-12);
    assert!((value(&json, &element("coil"), "flux_wb") - 0.09).abs() < 1.0e-14);
    assert!((value(&json, &element("coil"), "current_a") - 2.6).abs() < 1.0e-12);
    let capacitor_work = value(&json, "\"port\":\"cap.p\"", "energy_in_j");
    let inductor_work = value(&json, "\"port\":\"coil.p\"", "energy_in_j");
    assert!((capacitor_work - (nonlinear_energy(0.036, 0.02, 2.0) - nonlinear_energy(0.012, 0.02, 2.0))).abs() < 1.0e-12);
    assert!((inductor_work - (nonlinear_energy(0.09, 0.05, 1.0) - nonlinear_energy(0.03, 0.05, 1.0))).abs() < 1.0e-12);

    for (file, original, replacement, named_field) in [
        ("bad_charge_curve.json", "[0.02, 2.0]", "[0.02, -2.0]", "charge_voltage_curve"),
        ("bad_flux_curve.json", "[0.05, 1.0]", "[0.05, 0.0]", "flux_current_curve"),
        ("unsupported_initial_charge.json", "\"charge_c\": 0.012", "\"charge_c\": 0.5", "charge_c"),
    ] {
        let changed = scene.replace(original, replacement);
        assert_ne!(changed, scene);
        let (code, stderr) = refused(file, &changed);
        assert_eq!(code, 3, "{stderr}");
        assert!(stderr.contains(named_field), "{stderr}");
    }
    let changed = scene.replace("\"current_a\": 0.08", "\"current_a\": 10.0");
    let (code, stderr) = refused("nonlinear_storage_support_exit.json", &changed);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("cap") && stderr.contains("outside the supplied curve support"), "{stderr}");
}

#[test]
fn the_transformer_scene_has_reflected_stiffness_with_either_signed_ratio() {
    let scene = std::fs::read_to_string(example("geared_oscillator.json")).unwrap();
    for ratio in [2.5_f64, -2.5] {
        let changed = if ratio < 0.0 {
            scene.replace("\"ratio\": 2.5", "\"ratio\": -2.5")
                .replace("\"extension_m\": -0.25", "\"extension_m\": 0.25")
        } else {
            scene.clone()
        };
        let path = scratch(if ratio > 0.0 { "gear_forward.json" } else { "gear_reverse.json" }, &changed);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        let (m, k, x0, v0, dt, steps) = (2.0_f64, 80.0_f64, 0.1_f64, -0.4_f64, 0.002_f64, 1_500_u32);
        let omega = (k * ratio * ratio / m).sqrt();
        let angle = f64::from(steps) * 2.0 * (0.5 * omega * dt).atan();
        let expected_x = x0 * angle.cos() + v0 / omega * angle.sin();
        let expected_v = v0 * angle.cos() - omega * x0 * angle.sin();
        let reached_x = value(&json, &element("mass"), "position_m");
        let reached_v = value(&json, &element("mass"), "velocity_m_s");
        let extension = value(&json, &element("spring"), "extension_m");
        assert!((reached_x - expected_x).abs() <= 1.0e-10, "{reached_x} vs {expected_x}");
        assert!((reached_v - expected_v).abs() <= 1.0e-9, "{reached_v} vs {expected_v}");
        assert!((extension + ratio * reached_x).abs() <= 1.0e-12);
        let initial = 0.5 * (m * v0 * v0 + k * (ratio * x0).powi(2));
        let final_energy = 0.5 * (m * reached_v * reached_v + k * extension * extension);
        assert!((final_energy - initial).abs() <= 1.0e-10 * initial);
        assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-11);
    }
    let zero_ratio = scene.replace("\"ratio\": 2.5", "\"ratio\": 0");
    let (code, stderr) = refused("gear_zero_ratio.json", &zero_ratio);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("ratio"), "{stderr}");
}

#[test]
fn a_rectangular_power_map_suspends_a_rigid_body_without_constraining_transverse_motion() {
    let scene = std::fs::read_to_string(example("projected_suspension.json")).unwrap();
    let (mass, stiffness, dt, steps) = (2.0_f64, 18.0_f64, 0.01_f64, 100_u32);
    let direction = [0.6_f64, -0.8, 0.0];
    let initial_x = [0.2_f64, -0.1, 0.4];
    let initial_v = [0.3_f64, 0.2, -0.15];
    let (q0, u0) = (0.2_f64, 0.02_f64);
    for (name, ratio) in [("forward", 1.0_f64), ("reversed", -1.0), ("scaled", 2.0)] {
        let changed = scene
            .replace(
                "[[0.6, -0.8, 0.0]]",
                &format!(
                    "[[{}, {}, 0.0]]",
                    ratio * direction[0],
                    ratio * direction[1]
                ),
            )
            .replace(
                "\"extension_m\": -0.2",
                &format!("\"extension_m\": {}", -ratio * q0),
            );
        let path = scratch(&format!("projection_{name}.json"), &changed);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        assert!(stderr.is_empty());
        let time = f64::from(steps) * dt;
        assert!((value(&json, "\"schema\"", "time_s") - time).abs() <= 1.0e-12);

        // With n a unit vector, x = n q + x_perp and q'' = -k r^2 q / m.
        // The native body under a held force and the spring's discrete
        // gradient give the midpoint phase; transverse motion stays free.
        let omega = (stiffness * ratio * ratio / mass).sqrt();
        let phase = f64::from(steps) * 2.0 * (0.5 * omega * dt).atan();
        let q = q0 * phase.cos() + u0 / omega * phase.sin();
        let u = u0 * phase.cos() - omega * q0 * phase.sin();
        let mut kinetic = 0.0;
        for (axis, suffix) in ["x", "y", "z"].into_iter().enumerate() {
            let transverse_v = initial_v[axis] - direction[axis] * u0;
            let expected_x = initial_x[axis] + time * transverse_v + direction[axis] * (q - q0);
            let expected_v = transverse_v + direction[axis] * u;
            let position = value(
                &json,
                &element("carriage"),
                &format!("body0.position_{suffix}_m"),
            );
            let velocity = value(
                &json,
                &element("carriage"),
                &format!("body0.velocity_{suffix}_m_s"),
            );
            assert!(
                (position - expected_x).abs() <= 1.0e-9,
                "{name} {suffix}: {position} vs {expected_x}"
            );
            assert!(
                (velocity - expected_v).abs() <= 1.0e-9,
                "{name} {suffix}: {velocity} vs {expected_v}"
            );
            kinetic += 0.5 * mass * velocity * velocity;
        }
        let extension = value(&json, &element("spring"), "extension_m");
        assert!(
            (extension + ratio * q).abs() <= 1.0e-9,
            "{name}: {extension} vs {}",
            -ratio * q
        );
        let initial_energy = 0.5 * mass * initial_v.iter().map(|v| v * v).sum::<f64>()
            + 0.5 * stiffness * (ratio * q0).powi(2);
        assert!(
            (kinetic + 0.5 * stiffness * extension * extension - initial_energy).abs() <= 1.0e-9
        );
        assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-9);
    }
}

#[test]
fn power_map_scenes_require_rectangular_finite_matrices_and_valid_ports() {
    let scene = std::fs::read_to_string(example("projected_suspension.json")).unwrap();
    let matrix = "[[0.6, -0.8, 0.0]]";
    for (name, changed, reason) in [
        ("scalar", "0.6", "matrix row count must be 1"),
        ("empty", "[]", "matrix row count must be 1"),
        (
            "transposed",
            "[[0.6], [-0.8], [0.0]]",
            "matrix row count must be 1",
        ),
        ("short", "[[0.6, -0.8]]", "matrix[0] column count must be 3"),
        ("row_scalar", "[0.6]", "matrix[0] column count must be 3"),
        (
            "text",
            "[[0.6, \"west\", 0.0]]",
            "matrix[0][1] must be a finite number",
        ),
        ("nonfinite", "[[0.6, 1e999, 0.0]]", "not finite"),
    ] {
        let (code, stderr) = refused(&format!("map_{name}.json"), &scene.replace(matrix, changed));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
    for (name, from, to, reason) in [
        (
            "missing_matrix",
            ", \"matrix\": [[0.6, -0.8, 0.0]]",
            "",
            "matrix row count must be 1",
        ),
        (
            "unknown_field",
            "\"matrix\":",
            "\"weights\":",
            "unsupported field",
        ),
        (
            "missing_port",
            "[\"carriage.force0\", \"spring.a\"]",
            "[\"carriage.force0\"]",
            "power_map is",
        ),
        (
            "port_type",
            "\"carriage.force0\"",
            "7",
            "a port is a string",
        ),
        (
            "port_spelling",
            "\"carriage.force0\"",
            "\"carriage\"",
            "write it as element.port",
        ),
        (
            "unknown_element",
            "\"carriage.force0\"",
            "\"missing.force0\"",
            "no element named",
        ),
        (
            "unknown_port",
            "\"carriage.force0\"",
            "\"carriage.force1\"",
            "a rigid_network has ports",
        ),
    ] {
        assert!(scene.contains(from), "{name}: {from:?}");
        let (code, stderr) = refused(&format!("map_{name}.json"), &scene.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
    // Valid matrix syntax still goes through the kernel's role, carrier
    // and connection-ownership checks.
    let backwards = scene
        .replace(
            "[\"carriage.force0\", \"spring.a\"]",
            "[\"spring.a\", \"carriage.force0\"]",
        )
        .replace(matrix, "[[0.6], [-0.8], [0.0]]");
    let (code, stderr) = refused("map_roles.json", &backwards);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("spring.a"), "{stderr}");
    let heat = scene
        .replace(
            "\"spring\": {",
            "\"case\": {\"kind\": \"heat_capacity\", \"capacity_j_k\": 5, \"temperature_k\": 300},\n    \"spring\": {",
        )
        .replace("[\"carriage.force0\", \"spring.a\"]", "[\"carriage.force0\", \"case.heat\"]");
    let (code, stderr) = refused("map_heat.json", &heat);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("case.heat"), "{stderr}");
    let duplicate = scene.replace(
        "{\"prescribe\": \"spring.b\"",
        "{\"prescribe\": \"carriage.force0\", \"force_n\": [0, 0, 0]},\n    {\"prescribe\": \"spring.b\"",
    );
    let (code, stderr) = refused("map_connected_twice.json", &duplicate);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("carriage.force0"), "{stderr}");
}

#[test]
fn a_resistance_matrix_couples_rigid_body_motion_and_routes_all_lost_energy_to_heat() {
    let (code, json, stderr) =
        run_binary(&["--json", "system", &example("anisotropic_damping.json")]);
    assert_eq!(code, 0, "{stderr}");
    assert!(stderr.is_empty());
    let (dt, steps) = (0.01_f64, 100_i32);
    // R/m has the independent modes (x+y), (x-y), z with rates 2.5, 1.5, 1.
    // The held-force body and the resistance close at midpoint velocity.
    let decay = |rate: f64| ((1.0 - 0.5 * rate * dt) / (1.0 + 0.5 * rate * dt)).powi(steps);
    let (fast, slow, vertical) = (decay(2.5), decay(1.5), decay(1.0));
    let expected_v = [0.5 * (fast + slow), 0.5 * (fast - slow), 0.5 * vertical];
    let expected_x = [
        0.5 * ((1.0 - fast) / 2.5 + (1.0 - slow) / 1.5),
        0.5 * ((1.0 - fast) / 2.5 - (1.0 - slow) / 1.5),
        0.5 * (1.0 - vertical),
    ];
    assert!(
        expected_v[1] < 0.0,
        "off-diagonal damping excites transverse motion"
    );
    let mut kinetic = 0.0;
    for (axis, suffix) in ["x", "y", "z"].into_iter().enumerate() {
        let position = value(
            &json,
            &element("carriage"),
            &format!("body0.position_{suffix}_m"),
        );
        let velocity = value(
            &json,
            &element("carriage"),
            &format!("body0.velocity_{suffix}_m_s"),
        );
        assert!(
            (position - expected_x[axis]).abs() <= 1.0e-9,
            "{suffix}: {position} vs {}",
            expected_x[axis]
        );
        assert!(
            (velocity - expected_v[axis]).abs() <= 1.0e-9,
            "{suffix}: {velocity} vs {}",
            expected_v[axis]
        );
        kinetic += velocity * velocity; // m = 2 kg
    }
    let heat = 1.25 - expected_v.iter().map(|v| v * v).sum::<f64>();
    let dissipated = value(&json, &element("damping"), "dissipated_j");
    let stored_heat = value(&json, &element("case"), "stored_change_j");
    let temperature = value(&json, &element("case"), "temperature_k");
    assert!(
        (dissipated - heat).abs() <= 1.0e-9,
        "{dissipated} vs {heat}"
    );
    assert!(
        (stored_heat - heat).abs() <= 1.0e-9,
        "{stored_heat} vs {heat}"
    );
    assert!((temperature - 300.0 - heat / 10.0).abs() <= 1.0e-9);
    assert!((kinetic + stored_heat - 1.25).abs() <= 1.0e-9);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-9);
    assert!(value(&json, "\"totals\"", "entropy_generated_j_per_k") > 0.0);

    // Probes are means of the last step, not the final endpoint force.
    let fast_mean = fast / (1.0 - 0.5 * 2.5 * dt);
    let slow_mean = slow / (1.0 - 0.5 * 1.5 * dt);
    let mean = [
        0.5 * (fast_mean + slow_mean),
        0.5 * (fast_mean - slow_mean),
        0.5 * vertical / (1.0 - 0.5 * dt),
    ];
    let force = [
        -4.0 * mean[0] - mean[1],
        -mean[0] - 4.0 * mean[1],
        -2.0 * mean[2],
    ];
    for (row, expected) in force.into_iter().enumerate() {
        assert!(
            (value(&json, &element("damping"), &format!("u{row}.force_n")) - expected).abs()
                <= 1.0e-9
        );
    }
    let heat_w: f64 = force
        .iter()
        .zip(mean)
        .map(|(force, flow)| -force * flow)
        .sum();
    assert!((value(&json, &element("damping"), "heat_w") - heat_w).abs() <= 1.0e-9);
}

#[test]
fn resistance_matrices_use_the_declared_domain_and_need_no_adaptive_state_account() {
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.25, "step_s": 0.25,
  "adaptive": {"order": 2, "tolerances": {
    "case.temperature_k": {"absolute": 1e-8, "relative": 0}
  }},
  "elements": {
    "loss": {"kind": "resistance_matrix", "domain": "translation", "matrix": [[3.0]]},
    "case": {"kind": "heat_capacity", "capacity_j_k": 10, "temperature_k": 300}
  },
  "connections": [
    {"prescribe": "loss.u", "velocity_m_s": 2},
    {"route": "loss.heat", "to": [["case.heat", 1]]}
  ]
}"#;
    for (domain, flow, effort) in [
        ("translation", "velocity_m_s", "force_n"),
        ("rotation", "angular_velocity_rad_s", "torque_n_m"),
        ("electrical", "current_a", "voltage_v"),
    ] {
        let changed = scene
            .replace("translation", domain)
            .replace("velocity_m_s", flow);
        let path = scratch(&format!("resistance_{domain}.json"), &changed);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{domain}: {stderr}");
        // A prescribed attachment flow enters the element with the opposite sign.
        assert!((value(&json, &element("loss"), &format!("u0.{effort}")) + 6.0).abs() <= 1.0e-12);
        assert!((value(&json, &element("loss"), "heat_w") - 12.0).abs() <= 1.0e-12);
        assert!((value(&json, &element("loss"), "dissipated_j") - 3.0).abs() <= 1.0e-12);
        assert!((value(&json, &element("case"), "temperature_k") - 300.3).abs() <= 1.0e-11);
        assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-10);
    }
}

#[test]
fn resistance_matrix_scenes_refuse_invalid_laws_and_unconnected_heat() {
    let scene = std::fs::read_to_string(example("anisotropic_damping.json")).unwrap();
    let matrix = "[[4.0, 1.0, 0.0], [1.0, 4.0, 0.0], [0.0, 0.0, 2.0]]";
    for (name, changed, reason) in [
        ("empty", "[]", "matrix must have 1 to 64 square rows"),
        ("scalar", "2", "matrix must have 1 to 64 square rows"),
        (
            "ragged",
            "[[4, 1], [1, 4, 0], [0, 0, 2]]",
            "matrix[0] column count must be 3",
        ),
        (
            "row_scalar",
            "[4, 4, 2]",
            "matrix[0] column count must be 3",
        ),
        (
            "text",
            "[[4, 1, 0], [1, \"four\", 0], [0, 0, 2]]",
            "matrix[1][1] must be a finite number",
        ),
        (
            "nonfinite",
            "[[4, 1, 0], [1, 1e999, 0], [0, 0, 2]]",
            "not finite",
        ),
        (
            "asymmetric",
            "[[4, 1, 0], [2, 4, 0], [0, 0, 2]]",
            "differs from its mirror",
        ),
        (
            "semidefinite",
            "[[1, 1, 0], [1, 1, 0], [0, 0, 2]]",
            "not positive definite",
        ),
        (
            "indefinite",
            "[[1, 2, 0], [2, 1, 0], [0, 0, 2]]",
            "not positive definite",
        ),
    ] {
        let (code, stderr) = refused(
            &format!("resistance_{name}.json"),
            &scene.replace(matrix, changed),
        );
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
    let oversized = format!("[{}]", vec!["[1]"; 65].join(","));
    let (code, stderr) = refused(
        "resistance_oversized.json",
        &scene.replace(matrix, &oversized),
    );
    assert_eq!(code, 3, "{stderr}");
    assert!(
        stderr.contains("matrix must have 1 to 64 square rows"),
        "{stderr}"
    );
    for (name, from, to, reason) in [
        (
            "domain",
            "\"domain\": \"translation\"",
            "\"domain\": \"heat\"",
            "domain is one of",
        ),
        (
            "missing_domain",
            "\"domain\": \"translation\",",
            "",
            "domain is one of",
        ),
        (
            "unknown_field",
            "\"domain\":",
            "\"units\":",
            "unsupported field",
        ),
    ] {
        assert!(scene.contains(from));
        let (code, stderr) = refused(&format!("resistance_{name}.json"), &scene.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
    let no_heat = scene.replace(
        "{\"route\": \"damping.heat\", \"to\": [[\"case.heat\", 1.0]]}",
        "{\"prescribe\": \"case.heat\", \"heat_rate_w\": 0}",
    );
    let (code, stderr) = refused("resistance_unconnected_heat.json", &no_heat);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("damping.heat"), "{stderr}");
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

/// A scratch file of the calling test. Tests run side by side and several
/// use the same file names, so each test has a directory of its own, named
/// after its thread, which the test harness names after the test.
fn scratch(name: &str, text: &str) -> PathBuf {
    let thread = std::thread::current();
    let test: String = thread
        .name()
        .unwrap_or("unnamed")
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect();
    let directory = std::env::temp_dir()
        .join(format!("fs-cli-system-{}", std::process::id()))
        .join(test);
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

#[test]
fn a_motor_drives_a_load_through_a_compliant_shaft_and_a_run_reaches_rest() {
    // Second independent review, 2026-11-10: this scene was refused at step
    // 6912, because the shaft's two ends turn together at 140 rad/s while its
    // twist is tiny, and the spring's energy tolerance did not cover that.
    let scene = std::fs::read_to_string(example("dc_motor.json"))
        .unwrap()
        .replace(
            "    \"case\": {",
            "    \"shaft\": {\"kind\": \"torsion_spring\", \"stiffness_n_m_rad\": 500.0},\n    \
             \"load\": {\"kind\": \"inertia\", \"inertia_kg_m2\": 5.0e-4},\n    \
             \"load_bearing\": {\"kind\": \"rotary_damper\", \"damping_n_m_s_rad\": 1.0e-4},\n    \
             \"case\": {",
        )
        .replace(
            "\"bearing.a\"]},",
            "\"bearing.a\", \"shaft.a\"]},\n    \
             {\"join\": [\"load.p\", \"shaft.b\", \"load_bearing.a\"]},\n    \
             {\"prescribe\": \"load_bearing.b\", \"angular_velocity_rad_s\": 0},\n    \
             {\"route\": \"load_bearing.heat\", \"to\": [[\"case.heat\", 1.0]]},",
        );
    assert!(scene.contains("load_bearing.heat") && scene.contains("\"shaft\""));
    let path = scratch("motor_shaft.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    // The two inertias turn together, the shaft carrying the load's share.
    let rotor = value(&json, &element("rotor"), "angular_velocity_rad_s");
    let load = value(&json, &element("load"), "angular_velocity_rad_s");
    assert!(rotor > 50.0 && (rotor - load).abs() < 0.05 * rotor, "{rotor} vs {load}");
    let supplied = -value(&json, &element("supply"), "energy_in_j");
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-11 * supplied);

    // Two bodies through a wall, run until nothing more can change: the
    // kernel used to refuse at step 192 of this scene.
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 1000.0, "step_s": 0.5,
  "elements": {
    "x1": {"kind": "heat_capacity", "capacity_j_k": 10.0, "temperature_k": 400.0},
    "x2": {"kind": "heat_capacity", "capacity_j_k": 12.0, "temperature_k": 350.0},
    "wall": {"kind": "conductance", "conductance_w_k": 0.7}
  },
  "connections": [
    {"join": ["x1.heat", "wall.a"]},
    {"join": ["x2.heat", "wall.b"]}
  ]
}"#;
    let path = scratch("equilibrate.json", scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let mixed = (10.0 * 400.0 + 12.0 * 350.0) / 22.0;
    for name in ["x1", "x2"] {
        let reached = value(&json, &element(name), "temperature_k");
        assert!((reached - mixed).abs() <= 1.0e-12 * mixed, "{name}: {reached} vs {mixed}");
    }
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

#[test]
fn the_damped_pair_scene_runs_the_rigid_body_owner_and_matches_the_closed_form() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("damped_pair.json")]);
    assert_eq!(code, 0, "{stderr}");
    // Two 2 kg bodies closing at 2 m/s on an axial connector of 0.7 N s/m:
    // the reduced mass is 1 kg, so the relative speed decays as exp(-0.7 t)
    // and each body keeps its share of it; the heat is the relative kinetic
    // energy that has gone, 2 (1 - exp(-1.4 t)) J.
    let t = 0.5_f64;
    let decay = (-0.7 * t).exp();
    let heat = 2.0 * -(-1.4 * t).exp_m1();
    let pair = element("pair");
    assert!((value(&json, &pair, "body0.velocity_x_m_s") - decay).abs() <= 1.0e-13);
    assert!((value(&json, &pair, "body1.velocity_x_m_s") + decay).abs() <= 1.0e-13);
    // The centre of mass of the pair does not move. The owner moves a body
    // at the velocity it has after half a step's damping, so body 0 has
    // travelled the sum of those four drifts; that is second order in the
    // step and 1.3e-4 m short of the integral of the velocity.
    let dt = 0.125_f64;
    let drifted: f64 = (0..4).map(|n| dt * (-0.7 * (f64::from(n) + 0.5) * dt).exp()).sum();
    let position = value(&json, &pair, "body0.position_x_m");
    assert!((position - drifted).abs() <= 1.0e-13, "{position} vs {drifted}");
    let travelled = -(-0.7 * t).exp_m1() / 0.7;
    assert!((1.0e-4..2.0e-4).contains(&(travelled - position)), "{}", travelled - position);
    let apart = value(&json, &pair, "body1.position_x_m") + value(&json, &pair, "body0.position_x_m");
    assert!((apart - 1.0).abs() <= 1.0e-13, "{apart}");
    // Nothing leaves the line or turns: exactly.
    for quiet in ["body0.velocity_y_m_s", "body0.angular_velocity_z_rad_s", "body1.position_z_m"] {
        assert!(value(&json, &pair, quiet).abs() <= 0.0, "{quiet}");
    }
    assert!((value(&json, &pair, "body0.orientation_w") - 1.0).abs() <= 0.0);
    assert!((value(&json, &pair, "dissipated_j") - heat).abs() <= 1.0e-13, "{heat}");
    // A quarter of the heat to the shell (5 J/K), half to the oil (8 J/K),
    // the rest to the air.
    let shell = value(&json, &element("shell"), "temperature_k");
    let oil = value(&json, &element("oil"), "temperature_k");
    assert!((shell - 300.0 - 0.25 * heat / 5.0).abs() <= 1.0e-12, "{shell}");
    assert!((oil - 300.0 - 0.5 * heat / 8.0).abs() <= 1.0e-12, "{oil}");
    assert!((value(&json, "\"name\":\"air\"", "heat_j") - 0.25 * heat).abs() <= 1.0e-13);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-12);
    assert!(value(&json, "\"totals\"", "boundary_energy_in_j").abs() <= 0.0);
    // The last step's mean heat rate of the connector is reported.
    let last = 2.0 * ((-1.4 * (t - 0.125)).exp() - (-1.4 * t).exp()) / 0.125;
    assert!((value(&json, &pair, "connector0.heat_w") - last).abs() <= 1.0e-12);
    assert!(json.contains("elements stepped by the fs-mbd network owners"));
    let (code, text, _) = run_binary(&["system", &example("damped_pair.json")]);
    assert_eq!(code, 0);
    assert!(text.contains("rigid_network"), "{text}");
}

const ROTOR: &str = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.5, "step_s": 0.0125,
  "elements": {
    "rotor": {
      "kind": "rigid_network",
      "gravity_m_s2": [0, 0, -9.81],
      "bodies": [{"mass_kg": 2.0, "principal_inertia_kg_m2": [2.0, 3.0, 4.0],
                  "position_m": [0, 0, 10.0]}],
      "connectors": [],
      "torque_ports": [0]
    }
  },
  "connections": [
    {"prescribe": "rotor.force0", "force_n": [0, 0, 0]},
    {"prescribe": "rotor.torque0", "torque_n_m": [0, 0, 0.5]}
  ]
}"#;

#[test]
fn a_falling_body_spun_by_a_prescribed_torque_follows_its_two_equations() {
    let path = scratch("rotor.json", ROTOR);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    // Free fall from rest: v = -g t and z = 10 - g t^2 / 2. A torque of
    // 0.5 N m about the principal axis of 4 kg m^2: omega = tau t / I, and
    // the torque's work is tau^2 t^2 / (2 I), all of it rotation.
    let (t, g) = (0.5_f64, 9.81_f64);
    let rotor = element("rotor");
    assert!((value(&json, &rotor, "body0.velocity_z_m_s") + g * t).abs() <= 1.0e-12);
    assert!((value(&json, &rotor, "body0.position_z_m") - (10.0 - 0.5 * g * t * t)).abs() <= 1.0e-12);
    assert!((value(&json, &rotor, "body0.angular_velocity_z_rad_s") - 0.5 * t / 4.0).abs() <= 1.0e-14);
    assert!(value(&json, &rotor, "body0.angular_velocity_x_rad_s").abs() <= 0.0);
    let work = 0.5 * 0.5 * t * t / (2.0 * 4.0);
    assert!((value(&json, "rotor.torque0", "energy_in_j") - work).abs() <= 1.0e-14);
    assert!(value(&json, "rotor.force0", "energy_in_j").abs() <= 0.0);
    // Kinetic energy of the fall is paid by the potential: the stored
    // change is the torque's work alone.
    assert!((value(&json, &rotor, "stored_change_j") - work).abs() <= 1.0e-11);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-11);
}

#[test]
fn rigid_network_scenes_that_cannot_run_are_refused_with_the_reason() {
    let pair = std::fs::read_to_string(example("damped_pair.json")).unwrap();
    // A force port carries three components.
    let scalar = pair.replace("\"pair.force0\", \"force_n\": [0, 0, 0]", "\"pair.force0\", \"force_n\": 0");
    assert_ne!(scalar, pair);
    let (code, stderr) = refused("scalar.json", &scalar);
    assert_eq!(code, 3);
    assert!(stderr.contains("force_n for this port is an array of 3 waveforms"), "{stderr}");
    // A force port left unconnected is refused by the kernel, by name.
    let loose = pair.replace("{\"prescribe\": \"pair.force1\", \"force_n\": [0, 0, 0]},", "");
    assert_ne!(loose, pair);
    let (code, stderr) = refused("free.json", &loose);
    assert_eq!(code, 4);
    assert!(stderr.contains("pair.force1"), "{stderr}");
    // A connector joins two different bodies of the network.
    for ends in ["[0, 0]", "[0, 2]", "[0.5, 1]"] {
        let bad = pair.replace("\"bodies\": [0, 1]", &format!("\"bodies\": {ends}"));
        assert_ne!(bad, pair);
        let (code, stderr) = refused("ends.json", &bad);
        assert_eq!(code, 3);
        assert!(stderr.contains("connectors[0]: bodies is two different indices"), "{stderr}");
    }
    // Damping that is not passive, a field that does not exist, a mass
    // that is not one.
    let (code, stderr) = refused("negative.json", &pair.replace("[0.7, 0.0, 0.0]", "[-0.7, 0.0, 0.0]"));
    assert_eq!(code, 3);
    assert!(stderr.contains("connectors[0]: damping is not a passive form"), "{stderr}");
    let (code, stderr) = refused("field.json", &pair.replace("\"damping_n_s_m\"", "\"damping\""));
    assert_eq!(code, 3);
    assert!(stderr.contains("unsupported field \"damping\""), "{stderr}");
    let (code, stderr) = refused("mass.json", &pair.replacen("\"mass_kg\": 2.0", "\"mass_kg\": 0", 1));
    assert_eq!(code, 3);
    assert!(stderr.contains("bodies[0]: mass_kg must be positive"), "{stderr}");
    // A torque port on a body the network does not have, and one asked twice.
    for (ports, said) in [("[1]", "torque_ports is an array of body indices"), ("[0, 0]", "torque_ports")] {
        let bad = ROTOR.replace("\"torque_ports\": [0]", &format!("\"torque_ports\": {ports}"));
        let (code, stderr) = refused("torque.json", &bad);
        assert_eq!(code, 3);
        assert!(stderr.contains(said), "{stderr}");
    }
    // Error-controlled stepping has no account for the bodies.
    let adaptive = ROTOR.replace(
        "\"elements\"",
        "\"adaptive\": {\"tolerances\": {}},\n  \"elements\"",
    );
    let (code, stderr) = refused("adaptive.json", &adaptive);
    assert_eq!(code, 3);
    assert!(stderr.contains("rotor is a rigid_network"), "{stderr}");
}

/// Two 2 kg bodies on one axial connector with viscous damping `C`, one
/// Maxwell branch of stiffness `K` relaxing in `TAU` and a held force `f`
/// (RK4, 16384 steps). With u the relative speed (reduced mass 1 kg) and x
/// the branch's elastic extension: `u' = f - C u - K x`, `x' = u - x / TAU`,
/// heat `C u^2 + K x^2 / TAU`, actuator work `f u`. Returns
/// `[u, x, heat, work]` at `t`.
fn mount_reference(f: f64, u: f64, x: f64, t: f64) -> [f64; 4] {
    const C: f64 = 0.4;
    const K: f64 = 1.69;
    const TAU: f64 = 0.2;
    let rhs = |s: [f64; 4]| {
        [
            f - C * s[0] - K * s[1],
            s[0] - s[1] / TAU,
            C * s[0] * s[0] + K * s[1] * s[1] / TAU,
            f * s[0],
        ]
    };
    let add = |a: [f64; 4], b: [f64; 4], h: f64| core::array::from_fn(|i| a[i] + h * b[i]);
    let mut state = [u, x, 0.0, 0.0];
    let h = t / 16384.0;
    for _ in 0..16384 {
        let a = rhs(state);
        let b = rhs(add(state, a, 0.5 * h));
        let c = rhs(add(state, b, 0.5 * h));
        let d = rhs(add(state, c, h));
        for i in 0..4 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn the_viscoelastic_mount_scene_releases_its_stored_energy_as_motion_and_heat() {
    let (code, json, stderr) =
        run_binary(&["--json", "system", &example("viscoelastic_mount.json")]);
    assert_eq!(code, 0, "{stderr}");
    // At rest with the branch stretched by 0.5 m: 1.69 * 0.25 / 2 J stored.
    let [u, x, heat, _] = mount_reference(0.0, 0.0, 0.5, 0.25);
    let mount = element("mount");
    let speed = value(&json, &mount, "body1.velocity_x_m_s") - value(&json, &mount, "body0.velocity_x_m_s");
    assert!(u.abs() > 0.1, "the branch sets the bodies moving: {u} m/s");
    assert!((speed - u).abs() <= 1.0e-9, "{speed} vs {u}");
    let stored = value(&json, &mount, "connector0.branch0.energy_j");
    assert!((stored - 0.5 * 1.69 * x * x).abs() <= 1.0e-9, "{stored}");
    assert!((value(&json, &mount, "dissipated_j") - heat).abs() <= 1.0e-9, "{heat}");
    let casing = value(&json, &element("casing"), "temperature_k");
    assert!((casing - 300.0 - heat / 5.0).abs() <= 1.0e-9, "{casing}");
    // What left the branch is motion or heat.
    let kinetic = 0.5 * u * u;
    assert!((0.5 * 1.69 * 0.25 - stored - kinetic - heat).abs() <= 1.0e-9);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-12);
    // The first sample holds the starting state, before any step.
    let first = value(&json, "\"samples\"", "mount.connector0.branch0.energy_j");
    assert!((first - 0.5 * 1.69 * 0.25).abs() <= 1.0e-15, "{first}");
}

#[test]
fn an_actuator_on_a_viscoelastic_connector_does_the_work_of_its_own_equations() {
    // The mount relaxed and at rest, pushed apart by 3 N held on the
    // connector, in the first body's axes and in world axes (nothing turns,
    // so the two are the same run).
    let mount = std::fs::read_to_string(example("viscoelastic_mount.json")).unwrap();
    let mut reports = Vec::new();
    for frame in ["body", "world"] {
        let driven = mount
            .replace("\"extension_m\": [0.5, 0.0, 0.0]}", "\"extension_m\": [0.0, 0.0, 0.0]}")
            .replace(
                "\"damping_n_s_m\": [0.4, 0.0, 0.0],",
                &format!("\"damping_n_s_m\": [0.4, 0.0, 0.0], \"drive\": \"{frame}\","),
            )
            .replace(
                "{\"route\": \"mount.heat0\"",
                "{\"prescribe\": \"mount.drive0\", \"force_n\": [3.0, 0, 0]},\n    \
                 {\"prescribe\": \"mount.couple0\", \"torque_n_m\": [0, 0, 0]},\n    \
                 {\"route\": \"mount.heat0\"",
            );
        assert!(driven.contains("\"drive\"") && driven.contains("mount.couple0"));
        let path = scratch("driven_mount.json", &driven);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        let [u, x, heat, work] = mount_reference(3.0, 0.0, 0.0, 0.25);
        let name = element("mount");
        let speed = value(&json, &name, "body1.velocity_x_m_s") - value(&json, &name, "body0.velocity_x_m_s");
        assert!((speed - u).abs() <= 1.0e-9, "{frame}: {speed} vs {u}");
        let stored = value(&json, &name, "connector0.branch0.energy_j");
        assert!((stored - 0.5 * 1.69 * x * x).abs() <= 1.0e-9, "{frame}: {stored}");
        assert!((value(&json, &name, "dissipated_j") - heat).abs() <= 1.0e-9, "{frame}");
        // The actuator's work is its boundary's account, and nothing else
        // brought energy in.
        assert!(work > 0.05 && heat > 0.01, "work {work} J and heat {heat} J to measure");
        assert!((value(&json, "mount.drive0", "energy_in_j") - work).abs() <= 1.0e-9, "{frame}");
        assert!(value(&json, "mount.couple0", "energy_in_j").abs() <= 0.0);
        assert!((value(&json, "\"totals\"", "boundary_energy_in_j") - work).abs() <= 1.0e-9);
        assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-12);
        reports.push(json.replace("driven_mount.json", ""));
    }
    assert_eq!(reports[0], reports[1], "body and world axes agree while nothing turns");
}

#[test]
fn viscoelastic_scenes_that_cannot_run_are_refused_with_the_reason() {
    let mount = std::fs::read_to_string(example("viscoelastic_mount.json")).unwrap();
    let check = |name: &str, from: &str, to: &str, said: &str| {
        let bad = mount.replace(from, to);
        assert_ne!(bad, mount, "{name}");
        let (code, stderr) = refused(name, &bad);
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(said), "{name}: {stderr}");
    };
    check(
        "negative.json",
        "\"stiffness_n_m\": [1.69, 0.0, 0.0]",
        "\"stiffness_n_m\": [-1.69, 0.0, 0.0]",
        "branches[0]: a stiffness is never negative",
    );
    check(
        "axis.json",
        "\"extension_m\": [0.5, 0.0, 0.0]",
        "\"extension_m\": [0.5, 0.1, 0.0]",
        "an extension or twist needs a stiffness on the same axis",
    );
    check(
        "time.json",
        "\"relaxation_time_s\": 0.2",
        "\"relaxation_time_s\": 0",
        "relaxation_time_s must be positive",
    );
    check("field.json", "\"relaxation_time_s\"", "\"relaxation_s\"", "unsupported field");
    check(
        "drive.json",
        "\"damping_n_s_m\": [0.4, 0.0, 0.0],",
        "\"damping_n_s_m\": [0.4, 0.0, 0.0], \"drive\": \"sideways\",",
        "drive is \"body\" or \"world\"",
    );
    // A rigid_network has no branches: the field is refused there.
    check(
        "kind.json",
        "\"kind\": \"viscoelastic_network\"",
        "\"kind\": \"rigid_network\"",
        "unsupported field \"branches\"",
    );
    // An actuator's ports must be connected like any other.
    let loose = mount.replace(
        "\"damping_n_s_m\": [0.4, 0.0, 0.0],",
        "\"damping_n_s_m\": [0.4, 0.0, 0.0], \"drive\": \"body\",",
    );
    let (code, stderr) = refused("loose.json", &loose);
    assert_eq!(code, 4);
    assert!(stderr.contains("mount.drive0") || stderr.contains("mount.couple0"), "{stderr}");
}

/// The self-heating pair at step `dt`: relative speed `u` (reduced mass
/// 1 kg), damping `0.7 a(T)` with `a = exp(20000 (1 / T - 1 / 300))`, all
/// heat into 1 J/K. The scheme holds the damping over a step at the
/// temperature the step began with, and the owner's connector source is
/// exact for a held damping; so a step is `u -> u exp(-0.7 a(T) dt)` and the
/// kinetic energy lost is the heat. Returns `[u, T]` after `steps` steps.
fn self_heating_recursion(dt: f64, steps: u32) -> [f64; 2] {
    let (mut u, mut kelvin) = (-2.0_f64, 300.0_f64);
    for _ in 0..steps {
        let factor = (20_000.0 * (1.0 / kelvin - 1.0 / 300.0)).exp();
        let next = u * (-0.7 * factor * dt).exp();
        kelvin += 0.5 * (u * u - next * next);
        u = next;
    }
    [u, kelvin]
}

/// The coupled equations themselves, `u' = -0.7 a(T) u` and
/// `T' = 0.7 a(T) u^2` (RK4, 16384 steps). Returns `[u, T]` at `t`.
fn self_heating_reference(t: f64) -> [f64; 2] {
    let rhs = |s: [f64; 2]| {
        let c = 0.7 * (20_000.0 * (1.0 / s[1] - 1.0 / 300.0)).exp();
        [-c * s[0], c * s[0] * s[0]]
    };
    let add = |a: [f64; 2], b: [f64; 2], h: f64| [a[0] + h * b[0], a[1] + h * b[1]];
    let mut state = [-2.0, 300.0];
    let h = t / 16384.0;
    for _ in 0..16384 {
        let a = rhs(state);
        let b = rhs(add(state, a, 0.5 * h));
        let c = rhs(add(state, b, 0.5 * h));
        let d = rhs(add(state, c, h));
        for i in 0..2 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    state
}

#[test]
fn the_self_heating_damper_scene_thins_as_its_fluid_warms() {
    let scene = std::fs::read_to_string(example("self_heating_damper.json")).unwrap();
    let run = |text: &str| {
        let path = scratch("self_heating.json", text);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        let pair = element("pair");
        let speed = value(&json, &pair, "body1.velocity_x_m_s") - value(&json, &pair, "body0.velocity_x_m_s");
        (speed, value(&json, &element("fluid"), "temperature_k"), json)
    };
    // The scene as shipped: sixteen steps of 1/32 s, against the scheme's
    // own recursion.
    let (speed, kelvin, json) = run(&scene);
    let [u, expected_k] = self_heating_recursion(0.03125, 16);
    assert!((speed - u).abs() <= 1.0e-11, "{speed} vs {u}");
    assert!((kelvin - expected_k).abs() <= 1.0e-11, "{kelvin} vs {expected_k}");
    // The feedback is real: the fluid warms by 0.93 K and the damping ends
    // at 0.81 of its reference value, so the pair has slowed less (to
    // 1.461 m/s) than it would without the shift (1.409 m/s).
    let thinned = (20_000.0 * (1.0 / kelvin - 1.0 / 300.0)).exp();
    assert!(kelvin > 300.9 && thinned < 0.82, "{kelvin} K, factor {thinned}");
    assert!(speed.abs() > 2.0 * (-0.7_f64 * 0.5).exp() + 0.05, "{speed}");
    // The heat is what the motion lost, whatever the feedback did.
    let heat = value(&json, &element("pair"), "dissipated_j");
    assert!((heat - 0.5 * (4.0 - speed * speed)).abs() <= 1.0e-12, "{heat}");
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-12);

    // First order in the step: against the coupled equations, half the step
    // leaves half the distance.
    let [exact_u, exact_k] = self_heating_reference(0.5);
    let halved = scene.replace("\"step_s\": 0.03125", "\"step_s\": 0.015625");
    assert_ne!(halved, scene);
    let (fine_speed, fine_k, _) = run(&halved);
    let [fine_u, _] = self_heating_recursion(0.015625, 32);
    assert!((fine_speed - fine_u).abs() <= 1.0e-11);
    let (coarse, fine) = ((kelvin - exact_k).abs(), (fine_k - exact_k).abs());
    assert!(coarse > 1.0e-3, "a distance to measure: {coarse} K");
    assert!((1.7..=2.3).contains(&(coarse / fine)), "{coarse} K then {fine} K");
    assert!((speed - exact_u).abs() > (fine_speed - exact_u).abs());

    // The same scene with the temperature held at 300 K from outside is the
    // damped pair without a shift.
    let held = scene.replace(
        "{\"signal\": \"fluid.temperature\", \"to\": [\"pair.temperature0\"]}",
        "{\"prescribe\": \"pair.temperature0\", \"temperature_k\": 300.0}",
    );
    assert_ne!(held, scene);
    let (held_speed, _, _) = run(&held);
    assert!((held_speed + 2.0 * (-0.7_f64 * 0.5).exp()).abs() <= 1.0e-13, "{held_speed}");
}

#[test]
fn a_viscoelastic_mount_held_warm_relaxes_faster_by_its_shift() {
    // The mount of viscoelastic_mount.json with a shift, held at 330 K: its
    // viscosity and its relaxation time are both multiplied by
    // exp(20000 (1 / 330 - 1 / 300)) = 0.0023, its stiffness is not.
    let mount = std::fs::read_to_string(example("viscoelastic_mount.json")).unwrap();
    let warm = mount
        .replace(
            "\"connectors\": [",
            "\"shift\": {\"arrhenius\": {\"activation_energy_j_mol\": 166289.2523630648, \
             \"reference_k\": 300.0, \"valid_k\": [280.0, 360.0]}},\n      \"connectors\": [",
        )
        .replace(
            "{\"route\": \"mount.heat0\"",
            "{\"prescribe\": \"mount.temperature0\", \"temperature_k\": 330.0},\n    \
             {\"route\": \"mount.heat0\"",
        );
    assert!(warm.contains("\"shift\"") && warm.contains("mount.temperature0"));
    let path = scratch("warm_mount.json", &warm);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let factor = (20_000.0_f64 * (1.0 / 330.0 - 1.0 / 300.0)).exp();
    // u' = -c a u - K x, x' = u - x / (tau a), from rest with x = 0.5 m.
    let (c, k, tau) = (0.4 * factor, 1.69_f64, 0.2 * factor);
    let rhs = |s: [f64; 3]| {
        [
            -c * s[0] - k * s[1],
            s[0] - s[1] / tau,
            c * s[0] * s[0] + k * s[1] * s[1] / tau,
        ]
    };
    let add = |a: [f64; 3], b: [f64; 3], h: f64| core::array::from_fn(|i| a[i] + h * b[i]);
    let mut state = [0.0, 0.5, 0.0];
    let h = 0.25 / 1_048_576.0;
    for _ in 0..1_048_576 {
        let a = rhs(state);
        let b = rhs(add(state, a, 0.5 * h));
        let c = rhs(add(state, b, 0.5 * h));
        let d = rhs(add(state, c, h));
        for i in 0..3 {
            state[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let [u, x, heat] = state;
    let name = element("mount");
    let speed = value(&json, &name, "body1.velocity_x_m_s") - value(&json, &name, "body0.velocity_x_m_s");
    assert!((speed - u).abs() <= 1.0e-9, "{speed} vs {u}");
    let stored = value(&json, &name, "connector0.branch0.energy_j");
    assert!((stored - 0.5 * k * x * x).abs() <= 1.0e-9, "{stored}");
    assert!((value(&json, &name, "dissipated_j") - heat).abs() <= 1.0e-9, "{heat}");
    // At 330 K the branch has let go of nearly all it held (0.21 J), which
    // the mount at its reference temperature has not.
    assert!(heat > 0.2 && stored < 1.0e-6, "heat {heat} J, stored {stored} J");
}

#[test]
fn temperature_feedback_that_cannot_run_is_refused_with_the_reason() {
    let scene = std::fs::read_to_string(example("self_heating_damper.json")).unwrap();
    let check = |name: &str, from: &str, to: &str, code_expected: i32, said: &str| {
        let bad = scene.replace(from, to);
        assert_ne!(bad, scene, "{name}");
        let (code, stderr) = refused(name, &bad);
        assert_eq!(code, code_expected, "{name}: {stderr}");
        assert!(stderr.contains(said), "{name}: {stderr}");
    };
    // A shift is an Arrhenius or a WLF calibration, with its band.
    check("model.json", "\"arrhenius\"", "\"exponential\"", 3, "expected {\"arrhenius\": ...} or {\"wlf\": ...}");
    check("band.json", "\"valid_k\": [280.0, 360.0]", "\"valid_k\": [310.0, 360.0]", 3, "pair.shift");
    check("field.json", "\"reference_k\"", "\"reference\"", 3, "unsupported field");
    // The signal input of a shifted network must be fed.
    check(
        "unfed.json",
        ",\n    {\"signal\": \"fluid.temperature\", \"to\": [\"pair.temperature0\"]}",
        "",
        4,
        "pair.temperature0",
    );
    // A signal goes from an output to inputs.
    check(
        "backwards.json",
        "{\"signal\": \"fluid.temperature\", \"to\": [\"pair.temperature0\"]}",
        "{\"signal\": \"pair.temperature0\", \"to\": [\"fluid.temperature\"]}",
        4,
        "signal",
    );
    check("nowhere.json", "\"to\": [\"pair.temperature0\"]", "\"to\": []", 3, "to is an array of the ports that read the signal");
    // A network without a shift has no temperature input.
    let (code, stderr) = refused(
        "noshift.json",
        &std::fs::read_to_string(example("damped_pair.json")).unwrap().replace(
            "{\"prescribe\": \"pair.force0\"",
            "{\"signal\": \"shell.temperature\", \"to\": [\"pair.temperature0\"]},\n    {\"prescribe\": \"pair.force0\"",
        ),
    );
    assert_eq!(code, 3);
    assert!(stderr.contains("a rigid_network has ports"), "{stderr}");
    // A temperature that leaves the shift's band refuses the step, by name.
    let (code, stderr) = refused(
        "hot.json",
        &scene.replace("\"valid_k\": [280.0, 360.0]", "\"valid_k\": [280.0, 300.5]"),
    );
    assert_eq!(code, 4);
    assert!(stderr.contains("component `pair` refused"), "{stderr}");
    assert!(stderr.contains("temperature leaves the declared shift calibration band"), "{stderr}");
}

#[test]
fn the_rack_and_pinion_scene_carries_a_motor_and_a_rack_at_the_pinion_radius_on_one_inertia() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("rack_and_pinion.json")]);
    assert_eq!(code, 0, "{stderr}");

    // J w' = T + r F - (b + r² c) w - r² k a and a' = w, by the midpoint
    // rule: the rack's spring, friction and load reach the pinion through
    // its radius, the motor and the bearing directly.
    let (inertia, torque, bearing) = (2.0e-3_f64, 0.6_f64, 1.0e-3_f64);
    let (radius, stiffness, friction, load) = (0.02_f64, 4000.0_f64, 40.0_f64, -12.0_f64);
    let (dt, steps) = (5.0e-4_f64, 2_000_u32);
    let drive = torque + radius * load;
    let damping = bearing + radius * radius * friction;
    let reflected = radius * radius * stiffness;
    let (mut angle, mut speed) = (0.0_f64, 0.0_f64);
    let (mut motor_work, mut heat) = (0.0_f64, 0.0_f64);
    for _ in 0..steps {
        let mean = (drive - reflected * angle + 2.0 * inertia * speed / dt)
            / (2.0 * inertia / dt + damping + 0.5 * reflected * dt);
        angle += dt * mean;
        speed = 2.0 * mean - speed;
        motor_work += dt * torque * mean;
        heat += dt * damping * mean * mean;
    }
    let reached_angle = value(&json, &element("pinion"), "angle_rad");
    let reached_speed = value(&json, &element("pinion"), "angular_velocity_rad_s");
    assert!((reached_angle - angle).abs() <= 1.0e-10 * angle, "{reached_angle} vs {angle}");
    assert!((reached_speed - speed).abs() <= 1.0e-10, "{reached_speed} vs {speed}");
    // Close to rest at (T + r F) / (r² k) after four decay times.
    let rest = drive / reflected;
    assert!((reached_angle - rest).abs() <= 0.03 * rest, "{reached_angle} vs {rest}");

    // The spring is stretched by minus the rack's travel. The motor
    // delivers its work; the load, pushed back along the rack, receives
    // -F r a; both dampers warm the housing.
    let extension = value(&json, &element("rack_spring"), "extension_m");
    assert!((extension + radius * reached_angle).abs() <= 1.0e-13, "{extension}");
    let delivered = -value(&json, &element("motor"), "energy_in_j");
    assert!((delivered - motor_work).abs() <= 1.0e-10 * motor_work, "{delivered} vs {motor_work}");
    let lifted = value(&json, &element("load"), "energy_in_j");
    let expected = -load * radius * reached_angle;
    assert!((lifted - expected).abs() <= 1.0e-10 * expected, "{lifted} vs {expected}");
    let warmed = value(&json, &element("housing"), "stored_change_j");
    assert!((warmed - heat).abs() <= 1.0e-10 * motor_work, "{warmed} vs {heat}");
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-11 * motor_work);

    let scene = std::fs::read_to_string(example("rack_and_pinion.json")).unwrap();
    let geared_load = "{\"port\": \"load.p\", \"ratio\": 0.02}";
    let check = |name: &str, from: &str, to: &str, code: i32, needle: &str| {
        let changed = scene.replace(from, to);
        assert_ne!(changed, scene, "{name}");
        let (got, stderr) = refused(name, &changed);
        assert_eq!(got, code, "{name}: {stderr}");
        assert!(stderr.contains(needle), "{name}: {stderr}");
    };
    // The node sets the shared velocity; only a member sits behind a ratio.
    check(
        "rack_node.json",
        "\"pinion.p\",",
        "{\"port\": \"pinion.p\", \"ratio\": 2.0},",
        3,
        "a ratio belongs to a member",
    );
    check("rack_zero.json", geared_load, "{\"port\": \"load.p\", \"ratio\": 0}", 4, "ratio");
    check(
        "rack_field.json",
        geared_load,
        "{\"port\": \"load.p\", \"ratio\": 0.02, \"gain\": 1}",
        3,
        "unsupported field \"gain\"",
    );
    check("rack_missing.json", geared_load, "{\"port\": \"load.p\"}", 3, "missing ratio");
    // With no ratio anywhere the join is a plain one, and the rack's ports
    // are not the pinion's kind.
    let plain = scene
        .replace("{\"port\": \"rack_spring.a\", \"ratio\": 0.02}", "\"rack_spring.a\"")
        .replace("{\"port\": \"rack_friction.a\", \"ratio\": 0.02}", "\"rack_friction.a\"")
        .replace(geared_load, "\"load.p\"");
    let (code, stderr) = refused("rack_plain.json", &plain);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("rack_spring.a"), "{stderr}");
}

/// The pneumatic cylinder as `[p_tank, T_tank, p_chamber, T_chamber, x, v,
/// T_barrel]` at `t_end`, by classical RK4 on its own equations: two volumes
/// of perfect gas, a valve whose flow carries the enthalpy of the side it
/// leaves, a piston on the chamber, a spring, a seal whose friction warms the
/// barrel, and a film between the chamber's gas and the barrel.
fn cylinder_reference(t_end: f64, steps: u32) -> [f64; 7] {
    let (gas_constant, ratio, start_k) = (287.0_f64, 1.4_f64, 293.15_f64);
    let capacity = gas_constant / (ratio - 1.0);
    let (tank_m3, area, x0) = (2.0e-3_f64, 2.0e-3_f64, 0.02_f64);
    let (conductance, mass, outside) = (1.0e-8_f64, 0.5_f64, 200.0_f64);
    let (stiffness, damping, film, barrel) = (8000.0_f64, 60.0_f64, 0.5_f64, 150.0_f64);
    let held = |pressure: f64, volume: f64| {
        [pressure * volume / (gas_constant * start_k), pressure * volume / (ratio - 1.0)]
    };
    let ([tank_kg, tank_j], [chamber_kg, chamber_j]) = (held(6.0e5, tank_m3), held(1.0e5, area * x0));
    // m and U of the tank, m and U of the chamber, x, v, barrel temperature.
    let mut y = [tank_kg, tank_j, chamber_kg, chamber_j, x0, 0.0, start_k];
    let pressures = |y: &[f64; 7]| {
        [(ratio - 1.0) * y[1] / tank_m3, (ratio - 1.0) * y[3] / (area * y[4])]
    };
    let rhs = |y: &[f64; 7]| {
        let [p_tank, p_chamber] = pressures(y);
        let rate = conductance * (p_tank - p_chamber);
        let enthalpy = if rate >= 0.0 { ratio * y[1] / y[0] } else { ratio * y[3] / y[2] };
        let heat = film * (y[6] - y[3] / (y[2] * capacity));
        [
            -rate,
            -rate * enthalpy,
            rate,
            rate * enthalpy - p_chamber * area * y[5] + heat,
            y[5],
            (p_chamber * area - outside - stiffness * (y[4] - x0) - damping * y[5]) / mass,
            (damping * y[5] * y[5] - heat) / barrel,
        ]
    };
    let h = t_end / f64::from(steps);
    let moved = |y: &[f64; 7], d: &[f64; 7], by: f64| -> [f64; 7] {
        core::array::from_fn(|index| y[index] + by * d[index])
    };
    for _ in 0..steps {
        let k1 = rhs(&y);
        let k2 = rhs(&moved(&y, &k1, 0.5 * h));
        let k3 = rhs(&moved(&y, &k2, 0.5 * h));
        let k4 = rhs(&moved(&y, &k3, h));
        for index in 0..7 {
            y[index] += h * (k1[index] + 2.0 * k2[index] + 2.0 * k3[index] + k4[index]) / 6.0;
        }
    }
    let [p_tank, p_chamber] = pressures(&y);
    [
        p_tank,
        y[1] / (y[0] * capacity),
        p_chamber,
        y[3] / (y[2] * capacity),
        y[4],
        y[5],
        y[6],
    ]
}

/// The same seven numbers from a scene report.
fn cylinder_report(json: &str) -> [f64; 7] {
    [
        value(json, &element("tank"), "pressure_pa"),
        value(json, &element("tank"), "temperature_k"),
        value(json, &element("chamber"), "pressure_pa"),
        value(json, &element("chamber"), "temperature_k"),
        value(json, &element("piston"), "position_m"),
        value(json, &element("piston"), "velocity_m_s"),
        value(json, &element("barrel"), "temperature_k"),
    ]
}

/// Largest difference, each relative to the reference's value except the
/// velocity, which is taken against 1 m/s.
fn cylinder_error(got: &[f64; 7], want: &[f64; 7]) -> f64 {
    (0..7).fold(0.0_f64, |worst, index| {
        let scale = if index == 5 { 1.0 } else { want[index].abs() };
        worst.max((got[index] - want[index]).abs() / scale)
    })
}

#[test]
fn the_pneumatic_cylinder_scene_fills_from_its_tank_and_matches_its_own_equations() {
    let scene = std::fs::read_to_string(example("pneumatic_cylinder.json")).unwrap();
    // In mid-stroke, at 0.3 s, at two steps: second order against the
    // independent integration.
    let reference = cylinder_reference(0.3, 30_000);
    let mut errors = Vec::new();
    for step in ["2.0e-4", "1.0e-4"] {
        let changed = scene
            .replace("\"duration_s\": 3.0", "\"duration_s\": 0.3")
            .replace("\"step_s\": 2.0e-4", &format!("\"step_s\": {step}"));
        let path = scratch(&format!("cylinder_{step}.json"), &changed);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        errors.push(cylinder_error(&cylinder_report(&json), &reference));
    }
    // There the piston is 75 mm out and moving at 89 mm/s, and the chamber
    // is at 4.0 bar and 325 K against the tank's 5.5 bar. Measured errors:
    // 3.3e-8 and 8.1e-9.
    assert!(reference[5] > 0.05 && reference[3] > 320.0, "{reference:?}");
    let ratio = errors[0] / errors[1];
    assert!(errors[1] < 2.0e-8 && (3.8..4.2).contains(&ratio), "{errors:?}");

    // The scene as written, to 3 s.
    let (code, json, stderr) = run_binary(&["--json", "system", &example("pneumatic_cylinder.json")]);
    assert_eq!(code, 0, "{stderr}");
    let now = cylinder_report(&json);
    let reference = cylinder_reference(3.0, 300_000);
    // Measured 4.9e-9.
    let error = cylinder_error(&now, &reference);
    assert!(error < 2.0e-8, "{error}");
    // The piston has gone out more than 100 mm, the tank has cooled as it
    // emptied and the chamber warmed as it filled.
    assert!(now[4] > 0.12 && now[1] < 285.0 && now[3] > 300.0, "{now:?}");
    // No step was shortened, and the gas the tank lost is in the chamber.
    assert!(value(&json, "\"scene\"", "shortenings").abs() <= 0.0);
    let total = value(&json, &element("tank"), "mass_kg") + value(&json, &element("chamber"), "mass_kg");
    let initial = (6.0e5 * 2.0e-3 + 1.0e5 * 4.0e-5) / (287.0 * 293.15);
    assert!((total - initial).abs() <= 1.0e-12 * initial, "{total} vs {initial}");
    // The chamber's volume is the piston's face times its position, the
    // spring is compressed by the travel, and the outside has received the
    // work of 1 bar on the volume swept.
    let travel = now[4] - 0.02;
    let volume = value(&json, &element("chamber"), "volume_m3");
    assert!((volume - 2.0e-3 * now[4]).abs() <= 1.0e-12 * volume, "{volume}");
    assert!((value(&json, &element("load"), "extension_m") + travel).abs() <= 1.0e-13);
    let pushed = value(&json, &element("outside"), "energy_in_j");
    assert!((pushed - 200.0 * travel).abs() <= 1.0e-10 * pushed, "{pushed}");
    // Nearly at rest, the gas holds the spring and the outside pressure.
    let force = (now[2] - 1.0e5) * 2.0e-3;
    assert!((force - 8000.0 * travel).abs() <= 1.0e-3 * force, "{force}");
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-9);
}

/// A 20 l vessel of air at 1 bar and 300 K carrying one quantity, filled at
/// 5 g/s for 4 s from a line at 300 K that is all of that quantity. The
/// enthalpy rate is the mass rate times c_p T = 1004.5 J/(kg K) * 300 K.
const FILLED_VESSEL: &str = r#"{
  "schema": "frankensim.system.v0", "duration_s": 4.0, "step_s": 0.01,
  "elements": {
    "vessel": {"kind": "gas_volume", "gas_constant_j_kg_k": 287.0, "heat_capacity_ratio": 1.4,
               "volume_m3": 0.02, "pressure_pa": 1.0e5, "temperature_k": 300.0, "carried": [0.0]}
  },
  "connections": [
    {"prescribe": "vessel.stream", "stream": [5.0e-3, 1506.75, 5.0e-3]},
    {"prescribe": "vessel.heat", "heat_rate_w": 0}
  ]
}"#;

#[test]
fn gas_held_from_outside_fills_a_vessel_and_a_supply_sets_its_pressure() {
    // Filled by a prescribed stream: the gas pushed in does flow work, so
    // the vessel ends at 300 K * (m0 + ratio m_in) / (m0 + m_in), hotter
    // than the line, and the carried fraction is m_in / (m0 + m_in).
    let path = scratch("vessel_filled.json", FILLED_VESSEL);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (held, added) = (1.0e5 * 0.02 / (287.0 * 300.0), 0.02_f64);
    let kelvin = 300.0 * (held + 1.4 * added) / (held + added);
    let reached = value(&json, &element("vessel"), "temperature_k");
    assert!(kelvin > 350.0 && (reached - kelvin).abs() <= 1.0e-10 * kelvin, "{reached} vs {kelvin}");
    let fraction = value(&json, &element("vessel"), "carried0");
    assert!((fraction - added / (held + added)).abs() <= 1.0e-12, "{fraction}");
    let mass = value(&json, &element("vessel"), "mass_kg");
    assert!((mass - held - added).abs() <= 1.0e-12 * mass, "{mass}");
    // The boundary delivered m_in c_p T.
    let delivered = value(&json, "\"boundaries\"", "energy_in_j");
    assert!((delivered - added * 1004.5 * 300.0).abs() <= 1.0e-10 * delivered, "{delivered}");

    // Fed through a restriction from a supply held at 2 bar and 300 K: the
    // pressure ends at the supply's, and the energy that came in fixes the
    // temperature: U = p V / (ratio - 1) rose by m_in c_p T, which gives
    // 350 K here.
    let supplied = FILLED_VESSEL
        .replace(", \"carried\": [0.0]}", "},\n    \"inlet\": {\"kind\": \"restriction\", \"conductance_kg_s_pa\": 2.0e-7}")
        .replace(
            "{\"prescribe\": \"vessel.stream\", \"stream\": [5.0e-3, 1506.75, 5.0e-3]}",
            "{\"join\": [\"vessel.stream\", \"inlet.b\"]},\n    {\"prescribe\": \"inlet.a\", \"contents\": [2.0e5, 301350.0]}",
        )
        .replace("\"duration_s\": 4.0", "\"duration_s\": 40.0");
    let path = scratch("vessel_supplied.json", &supplied);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let pressure = value(&json, &element("vessel"), "pressure_pa");
    assert!((pressure - 2.0e5).abs() <= 1.0e-9 * 2.0e5, "{pressure}");
    let entered = 1.0e5 * 0.02 / (0.4 * 301_350.0);
    let kelvin = (2.0e5 * 0.02 / 0.4) / ((held + entered) * 717.5);
    let reached = value(&json, &element("vessel"), "temperature_k");
    assert!((reached - kelvin).abs() <= 1.0e-8 * kelvin, "{reached} vs {kelvin}");
    assert!((kelvin - 350.0).abs() < 0.01, "{kelvin}");

    // Refusals. A stream port takes its own quantity, with one waveform for
    // each thing it carries.
    let check = |name: &str, scene: &str, from: &str, to: &str, code: i32, needle: &str| {
        let changed = scene.replace(from, to);
        assert_ne!(changed, scene, "{name}");
        let (got, stderr) = refused(name, &changed);
        assert_eq!(got, code, "{name}: {stderr}");
        assert!(stderr.contains(needle), "{name}: {stderr}");
    };
    check("vessel_short.json", FILLED_VESSEL, "[5.0e-3, 1506.75, 5.0e-3]", "[5.0e-3, 1506.75]", 3, "stream for this port is an array of 3 waveforms");
    check("vessel_quantity.json", FILLED_VESSEL, "\"stream\": [5.0e-3, 1506.75, 5.0e-3]", "\"heat_rate_w\": 0", 3, "exactly one quantity, stream");
    check("vessel_contents.json", &supplied, "\"contents\": [2.0e5, 301350.0]", "\"stream\": [0, 0]", 3, "exactly one quantity, contents");
    check("vessel_ratio.json", FILLED_VESSEL, "\"heat_capacity_ratio\": 1.4", "\"heat_capacity_ratio\": 1.0", 3, "heat_capacity_ratio");
    check("vessel_fraction.json", FILLED_VESSEL, "\"carried\": [0.0]", "\"carried\": [1.5]", 3, "carried fraction");
    check("vessel_list.json", FILLED_VESSEL, "\"carried\": [0.0]", "\"carried\": 1", 3, "carried is an array");
    check("vessel_field.json", FILLED_VESSEL, "\"carried\": [0.0]", "\"carried\": [0.0], \"lid\": 1", 3, "unsupported field \"lid\"");
    check("vessel_wall.json", FILLED_VESSEL, "\"carried\": [0.0]", "\"carried\": [0.0], \"wall\": 1", 3, "wall is true or false");
    // A wall left unconnected is refused by the kernel, by name.
    check("vessel_loose.json", FILLED_VESSEL, "\"carried\": [0.0]", "\"carried\": [0.0], \"wall\": true", 4, "vessel.p");
    check("inlet_carried.json", &supplied, "2.0e-7}", "2.0e-7, \"carried\": 1.5}", 3, "carried is the number of carried quantities");
    // A restriction that carries something the vessel does not is refused
    // by the kernel's port match.
    check("inlet_mismatch.json", &supplied, "2.0e-7}", "2.0e-7, \"carried\": 1}", 4, "vessel.stream");
    // Two volumes have nothing between them, and one flow is not handed to
    // two of them.
    let cylinder = std::fs::read_to_string(example("pneumatic_cylinder.json")).unwrap();
    check("cylinder_direct.json", &cylinder, "[\"tank.stream\", \"valve.a\"]", "[\"tank.stream\", \"chamber.stream\"]", 3, "this one has 2");
    check("cylinder_open.json", &cylinder, "    {\"join\": [\"chamber.stream\", \"valve.b\"]},\n", "", 4, "valve.b");
    // A gas volume has no error account for adaptive stepping.
    check(
        "vessel_adaptive.json",
        FILLED_VESSEL,
        "\"step_s\": 0.01,",
        "\"step_s\": 0.01, \"adaptive\": {\"order\": 2, \"min_step_s\": 1e-9, \"max_attempts\": 8, \"tolerances\": {}},",
        3,
        "whose state has no error account here",
    );
}

/// One step of the midpoint rule for two bodies on a plane of no heat
/// capacity that carries a heat source: `C_k (T_k' - T_k) = dt G_k (T_n -
/// m_k)` with `m_k = (T_k + T_k') / 2`, and the plane's balance
/// `P + sum G_k (m_k - T_n) = 0`. Returns the bodies' end temperatures and
/// the plane's temperature.
fn plane_step(t: [f64; 2], c: [f64; 2], g: [f64; 2], power: f64, dt: f64) -> ([f64; 2], f64) {
    let total = g[0] + g[1];
    let a = [dt * g[0] / c[0], dt * g[1] / c[1]];
    // (2 + a_k) m_k - a_k T_n = 2 T_k, with T_n = (P + G_1 m_1 + G_2 m_2) / G.
    let rows = [
        [2.0 + a[0] - a[0] * g[0] / total, -a[0] * g[1] / total],
        [-a[1] * g[0] / total, 2.0 + a[1] - a[1] * g[1] / total],
    ];
    let rhs = [2.0 * t[0] + a[0] * power / total, 2.0 * t[1] + a[1] * power / total];
    let det = rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0];
    let mean = [
        (rhs[0] * rows[1][1] - rows[0][1] * rhs[1]) / det,
        (rows[0][0] * rhs[1] - rows[1][0] * rhs[0]) / det,
    ];
    let plane = (power + g[0] * mean[0] + g[1] * mean[1]) / total;
    ([2.0 * mean[0] - t[0], 2.0 * mean[1] - t[1]], plane)
}

#[test]
fn the_heated_plane_scene_divides_its_heat_by_where_the_bodies_put_it() {
    let scene = std::fs::read_to_string(example("heated_plane.json")).unwrap();
    let (c, g, power, dt, steps) = ([50.0_f64, 120.0], [2.0_f64, 6.0], 30.0_f64, 0.5_f64, 400_u32);
    let (mut t, mut plane, mut disc_share) = ([400.0_f64, 300.0], f64::NAN, f64::NAN);
    let mut pad_after_one_step = f64::NAN;
    for step in 0..steps {
        let before = t[1];
        (t, plane) = plane_step(t, c, g, power, dt);
        disc_share = c[1] * (t[1] - before) / (power * dt);
        if step == 0 {
            pad_after_one_step = t[0];
        }
    }
    let given_j = power * dt * f64::from(steps);

    // The heat source as a heater on the plane, and as a heat rate
    // prescribed on the plane's source port with no heater and no starting
    // temperature: one device, written two ways.
    let prescribed = scene
        .replace(", \"temperature_k\": 325.0}", "}")
        .replace(",\n    \"rubbing\": {\"kind\": \"heater\", \"power_w\": 30.0}", "")
        .replace(
            "{\"join\": [\"plane.source\", \"rubbing.heat\"]}",
            "{\"prescribe\": \"plane.source\", \"heat_rate_w\": 30.0}",
        );
    assert!(!prescribed.contains("rubbing") && !prescribed.contains("325.0"), "{prescribed}");
    let path = scratch("prescribed.json", &prescribed);
    for (scene, heater) in [(example("heated_plane.json"), true), (path.to_str().unwrap().to_owned(), false)] {
        let (code, json, stderr) = run_binary(&["--json", "system", &scene]);
        assert_eq!(code, 0, "{stderr}");
        let reached = [
            value(&json, &element("pad"), "temperature_k"),
            value(&json, &element("disc"), "temperature_k"),
        ];
        for body in 0..2 {
            assert!(
                (reached[body] - t[body]).abs() <= 1.0e-9 * t[body],
                "{} vs {}",
                reached[body],
                t[body]
            );
        }
        // The plane reports the temperature its balance gave it in the
        // last step: above both bodies, which it is heating.
        let node = value(&json, &element("plane"), "temperature_k");
        assert!((node - plane).abs() <= 1.0e-9 * plane, "{node} vs {plane}");
        assert!(node > reached[0] && node > reached[1], "{node} {reached:?}");
        // Every joule of the source is in the two bodies.
        let stored_j = c[0] * (reached[0] - 400.0) + c[1] * (reached[1] - 300.0);
        assert!((stored_j - given_j).abs() <= 1.0e-9 * given_j, "{stored_j} vs {given_j}");
        if heater {
            let account_j = -value(&json, &element("rubbing"), "energy_in_j");
            assert!((account_j - given_j).abs() <= 1.0e-9 * given_j, "{account_j}");
        }
    }
    // No fraction is prescribed. In the first step the hot pad cools: it
    // gives heat up through the plane, which is at 328.75 K before any of
    // it warms. By the end the two warm at nearly one rate and the disc
    // takes nearly its share of the heat capacity, 120 / 170. Both are
    // read from the binary, on the scene cut to one step and to its last
    // two, and agree with the recursion above.
    let stopped_after = |duration_s: &str| {
        let cut = scene
            .replace("\"duration_s\": 200.0", &format!("\"duration_s\": {duration_s}"))
            .replace("\"sample_every_s\": 20.0", "\"sample_every_s\": 0.5");
        assert!(cut.contains(&format!("\"duration_s\": {duration_s},")), "{cut}");
        let path = scratch(&format!("stopped_{duration_s}.json"), &cut);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        [
            value(&json, &element("pad"), "temperature_k"),
            value(&json, &element("disc"), "temperature_k"),
        ]
    };
    let [pad_then, _] = stopped_after("0.5");
    assert!(pad_then < 400.0, "{pad_then}");
    assert!((pad_then - pad_after_one_step).abs() <= 1.0e-9 * pad_then, "{pad_then}");
    let ([_, disc_before], [_, disc_last]) = (stopped_after("199.5"), stopped_after("200.0"));
    let share = c[1] * (disc_last - disc_before) / (power * dt);
    assert!((share - 120.0 / 170.0).abs() < 2.0e-3, "{share}");
    assert!((share - disc_share).abs() <= 1.0e-9, "{share} vs {disc_share}");
}

#[test]
fn a_heat_node_with_a_held_side_brings_its_body_to_where_it_exchanges_nothing() {
    // The heated plane with its disc taken away and that side of the plane
    // held at 300 K: `temperature_k` prescribed on a body port. The pad
    // then relaxes to the temperature at which it exchanges nothing,
    // 300 K + P / G1 = 305 K, at the rate G0 G1 / (G0 + G1) / C = 0.03 /s,
    // by the midpoint rule's factor each step.
    let scene = std::fs::read_to_string(example("heated_plane.json")).unwrap();
    let held = scene
        .replace(
            "    \"disc\": {\"kind\": \"heat_capacity\", \"capacity_j_k\": 120.0, \"temperature_k\": 300.0},\n",
            "",
        )
        .replace(
            "{\"join\": [\"disc.heat\", \"plane.body1\"]}",
            "{\"prescribe\": \"plane.body1\", \"temperature_k\": 300.0}",
        );
    assert!(!held.contains("disc") && held.contains("\"prescribe\": \"plane.body1\""), "{held}");
    let (rest_k, half_rate_step) = (300.0 + 30.0 / 6.0, 0.5 * 0.5 * (2.0 * 6.0 / 8.0 / 50.0));
    let factor: f64 = (1.0 - half_rate_step) / (1.0 + half_rate_step);

    let path = scratch("held.json", &held);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let pad = value(&json, &element("pad"), "temperature_k");
    let expected = rest_k + (400.0 - rest_k) * factor.powi(400);
    assert!((pad - expected).abs() <= 1.0e-12 * expected, "{pad} vs {expected}");
    // Still cooling at 200 s: the plane sits between the pad and where
    // both are going.
    let plane = value(&json, &element("plane"), "temperature_k");
    assert!(plane > rest_k && plane < pad, "{plane} {pad}");

    // Ten times as long, the pad has arrived and goes on stepping there.
    // (At rest the heat between pad and plane is a rounding-sized
    // difference of large terms across a junction that carries nothing;
    // the kernel refused that at step 1661 until its audit was told the
    // rounding of a step taken with no update.)
    let long = held
        .replace("\"duration_s\": 200.0", "\"duration_s\": 2000.0")
        .replace("\"sample_every_s\": 20.0", "\"sample_every_s\": 200.0");
    assert!(long.contains("\"duration_s\": 2000.0,"), "{long}");
    let path = scratch("held_long.json", &long);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    assert!(json.contains("\"steps\":4000,"), "{json}");
    for name in ["pad", "plane"] {
        let reached = value(&json, &element(name), "temperature_k");
        assert!((reached - rest_k).abs() <= 1.0e-9, "{name}: {reached}");
    }
}

#[test]
fn where_a_heat_node_s_solve_starts_does_not_change_the_heated_plane_s_report() {
    // `temperature_k` on a heat node is where the solve of its loop
    // starts, not a state. Every part of this scene's loop is linear, so
    // one update lands on the same numbers from any start: the reports
    // from five starts between 250 K and 5000 K are the same bytes. (A
    // loop with a part that is not linear would agree to the solve's
    // tolerance only.) Each is written to the same path, which the report
    // names.
    let scene = std::fs::read_to_string(example("heated_plane.json")).unwrap();
    let report_from = |start_k: &str| {
        let moved = scene.replace("\"temperature_k\": 325.0", &format!("\"temperature_k\": {start_k}"));
        assert!(moved.contains(&format!("\"temperature_k\": {start_k}}}")), "{moved}");
        let path = scratch("start.json", &moved);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        json
    };
    let usual = report_from("325.0");
    assert!(usual.contains("\"steps\":400,"), "{usual}");
    for start_k in ["250.0", "301.0", "600.0", "5000.0"] {
        assert_eq!(report_from(start_k), usual, "from {start_k} K");
    }
}

#[test]
fn heat_node_scenes_that_cannot_run_are_refused_with_the_reason() {
    let scene = std::fs::read_to_string(example("heated_plane.json")).unwrap();
    let check = |name: &str, from: &str, to: &str, code: i32, reason: &str| {
        assert!(scene.contains(from), "{name}: {from:?} is not in the scene");
        let (exit, stderr) = refused(name, &scene.replace(from, to));
        assert_eq!(exit, code, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    };
    let list = "\"conductances_w_k\": [2.0, 6.0]";
    check("node_scalar.json", list, "\"conductances_w_k\": 2.0", 3, "conductances_w_k is an array of 1 to 64");
    check("node_empty.json", list, "\"conductances_w_k\": []", 3, "conductances_w_k is an array of 1 to 64");
    check("node_text.json", list, "\"conductances_w_k\": [2.0, \"six\"]", 3, "conductances_w_k must be numbers");
    check("node_negative.json", list, "\"conductances_w_k\": [2.0, -6.0]", 3, "plane");
    check("node_field.json", list, "\"conductances_w_k\": [2.0, 6.0], \"capacity_j_k\": 1.0", 3, "unsupported field \"capacity_j_k\"");
    check("node_cold.json", "\"temperature_k\": 325.0", "\"temperature_k\": 0.0", 3, "temperature_k must be positive and finite");
    // A third conductance is a third body port, and the kernel refuses a
    // port left unconnected, by name.
    check("node_loose.json", list, "\"conductances_w_k\": [2.0, 6.0, 1.0]", 4, "plane.body2");
    // Two bodies cannot meet at one body port: nothing there sets a
    // temperature for both.
    check(
        "node_two.json",
        "{\"join\": [\"pad.heat\", \"plane.body0\"]},\n    {\"join\": [\"disc.heat\", \"plane.body1\"]}",
        "{\"join\": [\"pad.heat\", \"disc.heat\", \"plane.body0\"]},\n    {\"prescribe\": \"plane.body1\", \"temperature_k\": 300.0}",
        3,
        "this one has 2",
    );
    // With a heater on its source the temperature there is solved in a
    // loop, and the solve needs somewhere to start.
    // The scene is admitted and its first step fails: exit 4, the step's
    // code, where a scene refused on reading exits 3.
    let (exit, stderr) = refused("node_start.json", &scene.replace(", \"temperature_k\": 325.0}", "}"));
    assert_eq!(exit, 4, "{stderr}");
    assert!(
        stderr.contains("system-step")
            && stderr.contains("rubbing.heat")
            && stderr.contains("needs a positive starting value"),
        "{stderr}"
    );
}

#[test]
fn the_friction_belt_holds_breaks_away_and_catches_with_the_heat_in_its_lining() {
    // G1: piecewise constant acceleration, independent of the discrete law.
    // A 2 kg block needs 8 N to follow the first ramp (below the 10 N
    // static limit), then 16 N on the second. It slides at 4 N, hence
    // 2 m/s2, and catches the 3 m/s belt at 1.25 s. Its final displacement
    // is 2.875 m; the belt travelled 3.625 m. The 0.75 m of slip makes
    // 3 J of heat, while 9 J becomes kinetic energy: 12 J from the belt.
    let scene = std::fs::read_to_string(example("friction_belt.json")).unwrap();
    for (name, duration, position, velocity, slip, stuck, heat) in [
        ("held.json", "0.25", 0.125, 1.0, 0.0, 1.0, 0.0),
        ("sliding.json", "0.5", 0.4375, 1.5, -1.5, 0.0, 0.75),
        ("caught.json", "1.5", 2.875, 3.0, 0.0, 1.0, 3.0),
    ] {
        let path = scratch(
            name,
            &scene.replace(
                "\"duration_s\": 1.5",
                &format!("\"duration_s\": {duration}"),
            ),
        );
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        assert!(stderr.is_empty());
        for (owner, field, expected, tolerance) in [
            ("block", "position_m", position, 2.0e-8),
            ("block", "velocity_m_s", velocity, 2.0e-9),
            ("friction", "slip_m_s", slip, 2.0e-9),
            ("friction", "stuck", stuck, 0.0),
            ("friction", "dissipated_j", heat, 2.0e-8),
            ("lining", "temperature_k", 300.0 + heat / 50.0, 1.0e-9),
            ("belt", "energy_in_j", -(velocity * velocity + heat), 4.0e-8),
        ] {
            let found = value(&json, &element(owner), field);
            assert!(
                (found - expected).abs() <= tolerance,
                "{name} {owner}.{field}: {found} vs {expected}"
            );
        }
        assert!(
            value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
            "{json}"
        );
        assert!(json.contains("dry friction uses constant static/kinetic coefficients"));
    }
}

#[test]
fn a_slider_arrests_a_thrown_block_in_either_direction_and_returns_its_energy_as_heat() {
    // G1/G3: v^2/(2 mu_k N/m) is the stopping distance and m v^2/2
    // is the heat, unchanged when the motion is reversed. Initial slip
    // uses the scene's mover sign, not the opposite kernel port sign.
    for sign in [-1.0, 1.0] {
        let scene = format!(
            r#"{{
  "schema": "frankensim.system.v0", "duration_s": 2.0, "step_s": 0.005,
  "elements": {{
    "block": {{"kind": "mass", "mass_kg": 2.0, "velocity_m_s": {velocity}}},
    "friction": {{"kind": "slider", "normal_n": 20.0, "static_coefficient": 0.5,
                 "kinetic_coefficient": 0.2, "slip_m_s": {velocity}}},
    "lining": {{"kind": "heat_capacity", "capacity_j_k": 50.0, "temperature_k": 300.0}}
  }},
  "connections": [
    {{"join": ["block.p", "friction.a"]}},
    {{"prescribe": "friction.b", "velocity_m_s": 0.0}},
    {{"join": ["lining.heat", "friction.heat"]}}
  ]
}}"#,
            velocity = sign * 3.0
        );
        let path = scratch(
            if sign > 0.0 {
                "forward.json"
            } else {
                "backward.json"
            },
            &scene,
        );
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        assert!((value(&json, &element("block"), "position_m") - sign * 2.25).abs() < 1.0e-8);
        assert!(value(&json, &element("block"), "velocity_m_s").abs() < 1.0e-10);
        assert!(value(&json, &element("friction"), "slip_m_s").abs() < 1.0e-10);
        assert_eq!(value(&json, &element("friction"), "stuck"), 1.0);
        assert!((value(&json, &element("friction"), "dissipated_j") - 9.0).abs() < 1.0e-9);
        assert!((value(&json, &element("lining"), "temperature_k") - 300.18).abs() < 1.0e-9);
        assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9);
    }
}

#[test]
fn a_prescribed_slider_span_uses_the_waveform_endpoints_and_can_reverse() {
    // A sine's step mean is not either endpoint. Both ways of supplying
    // the same motion must carry the block through its reversals (max
    // force 2 * 0.4 * 2 pi < 10 N). The curved velocity retains the
    // kernel's finite-step position, endpoint-speed and heat errors.
    let scene = std::fs::read_to_string(example("friction_belt.json")).unwrap();
    let wave = r#"{"sine": {"amplitude": 0.4, "frequency_hz": 1.0}}"#;
    let table = r#"{"table": [[0.0, 0.0], [0.25, 1.0], [0.5, 3.0], [1.5, 3.0]]}"#;
    let source = scene.replace(table, wave);
    let boundary = source
        .replace(
            &format!(
                "    \"belt\": {{\"kind\": \"velocity_source\", \"velocity_m_s\": {wave}}},\n"
            ),
            "",
        )
        .replace(
            r#"{"join": ["belt.p", "friction.b"]}"#,
            &format!(r#"{{"prescribe": "friction.b", "velocity_m_s": {wave}}}"#),
        )
        .replace(
            r#"{"signal": "belt.span", "to": ["friction.span_b"]}"#,
            &format!(r#"{{"prescribe": "friction.span_b", "velocity_m_s": {wave}}}"#),
        );
    let mut projections = Vec::new();
    for (name, scene) in [("source.json", source), ("boundary.json", boundary)] {
        let path = scratch(name, &scene);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let position = value(&json, &element("block"), "position_m");
        assert!(
            (position - 0.4 / std::f64::consts::PI).abs() < 2.0e-5,
            "{name}: {position}"
        );
        let velocity = value(&json, &element("block"), "velocity_m_s");
        let heat = value(&json, &element("friction"), "dissipated_j");
        assert!(velocity.abs() < 1.0e-4, "{name}: {velocity}");
        assert_eq!(value(&json, &element("friction"), "stuck"), 1.0);
        for sample in json.split("{\"time_s\":").skip(1) {
            assert_eq!(value(sample, "\"values\"", "friction.stuck"), 1.0, "{name}");
        }
        assert!((0.0..1.0e-3).contains(&heat), "{name}: {heat}");
        projections.push([position, velocity, heat]);
    }
    for (source, boundary) in projections[0].iter().zip(&projections[1]) {
        assert!(
            (source - boundary).abs() < 1.0e-10,
            "{source} vs {boundary}"
        );
    }
}

#[test]
fn slider_scenes_refuse_invalid_laws_missing_spans_and_unsupported_error_control() {
    let scene = std::fs::read_to_string(example("friction_belt.json")).unwrap();
    for (name, from, to, wanted, reason) in [
        (
            "normal.json",
            "\"normal_n\": 20.0",
            "\"normal_n\": 0.0",
            3,
            "normal_n",
        ),
        (
            "coefficients.json",
            "\"kinetic_coefficient\": 0.2",
            "\"kinetic_coefficient\": 0.6",
            3,
            "kinetic_coefficient",
        ),
        (
            "slip.json",
            "\"prescribed_ends\": [\"b\"]",
            "\"prescribed_ends\": [\"b\"], \"slip_m_s\": \"moving\"",
            3,
            "slip_m_s must be a number",
        ),
        (
            "repeat.json",
            "\"prescribed_ends\": [\"b\"]",
            "\"prescribed_ends\": [\"b\", \"b\"]",
            3,
            "an end appears twice",
        ),
        (
            "end.json",
            "\"prescribed_ends\": [\"b\"]",
            "\"prescribed_ends\": [\"c\"]",
            3,
            "each end must be a or b",
        ),
        (
            "missing_span.json",
            "    {\"signal\": \"belt.span\", \"to\": [\"friction.span_b\"]},\n",
            "",
            4,
            "friction.span_b",
        ),
        (
            "adaptive.json",
            "\"duration_s\": 1.5",
            "\"adaptive\": {\"tolerances\": {}}, \"duration_s\": 1.5",
            3,
            "friction is a slider",
        ),
    ] {
        assert!(scene.contains(from), "{name}: test fixture changed");
        let (code, stderr) = refused(name, &scene.replace(from, to));
        assert_eq!(code, wanted, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
}
