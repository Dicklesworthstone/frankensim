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
