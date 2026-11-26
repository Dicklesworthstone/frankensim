//! `frankensim system`: scenes from examples/system run through the real
//! binary, and their reports agree with references computed here.

use std::path::PathBuf;
use std::process::Command;

use fs_system::elements::Waveform;

#[path = "system_command/checkpoints.rs"]
mod checkpoints;
#[path = "system_command/adaptive_networks.rs"]
mod adaptive_networks;
#[path = "system_command/timeline.rs"]
mod timeline;
#[path = "system_command/modal_network.rs"]
mod modal_network;
#[path = "system_command/contact.rs"]
mod contact;
#[path = "system_command/conduction.rs"]
mod conduction;
#[path = "system_command/solids.rs"]
mod solids;
#[path = "system_command/ohmic.rs"]
mod ohmic;
#[path = "system_command/laminar_pipe.rs"]
mod laminar_pipe;
#[path = "system_command/feedback.rs"]
mod feedback;
#[path = "system_command/gas.rs"]
mod gas;
#[path = "system_command/electrical_switching.rs"]
mod electrical_switching;
#[path = "system_command/dirac.rs"]
mod dirac;
#[path = "system_command/radiation.rs"]
mod radiation;
#[path = "system_command/dissipation.rs"]
mod dissipation;
#[path = "system_command/gap_field.rs"]
mod gap_field;
#[path = "system_command/battery.rs"]
mod battery;
#[path = "system_command/circuit.rs"]
mod circuit;
#[path = "system_command/motor_drive.rs"]
mod motor_drive;

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

fn scene_context(body: impl FnOnce(&fs_exec::Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let key = fs_exec::StreamKey {
            seed: 0,
            kernel_id: 0x5953,
            tile: 0,
            iteration: 0,
        };
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            key,
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        body(&cx, &gate);
    });
}

#[test]
fn controlled_sources_drive_mechanical_and_electrical_states_with_signed_work() {
    let path = example("controlled_sources.json");
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    // Independent midpoint recurrence on the authored piecewise affine commands.
    // The nominal step crosses a knot, so the actual schedule must resolve it.
    let mut time = 0.0_f64;
    let (mut position, mut velocity, mut charge) = (0.0, 1.0, 0.25);
    let force = |t: f64| 4.0 - 16.0 * (t - 0.25).clamp(0.0, 0.5);
    let current = |t: f64| -0.5 + 1.5 * (t - 0.25).clamp(0.0, 0.5);
    let mut steps = 0_u32;
    while time < 1.25 {
        let event = [0.25, 0.75, 1.25]
            .into_iter()
            .find(|knot| *knot > time)
            .unwrap();
        let dt = 0.375_f64.min(event - time);
        let next = time + dt;
        let next_velocity = velocity + dt * 0.5 * (force(time) + force(next)) / 2.0;
        position += dt * 0.5 * (velocity + next_velocity);
        charge -= dt * 0.5 * (current(time) + current(next));
        (time, velocity) = (next, next_velocity);
        steps += 1;
    }
    let drive_work = 1.0 - velocity * velocity;
    let supply_work = 0.25_f64.powi(2) - charge * charge;
    assert!(
        drive_work > 0.0,
        "the drive recovered the block's kinetic energy"
    );
    assert!(
        supply_work < 0.0,
        "the current source supplied stored electrical energy"
    );
    for (name, field, expected) in [
        ("block", "position_m", position),
        ("block", "velocity_m_s", velocity),
        ("capacitor", "charge_c", charge),
        ("drive", "energy_in_j", drive_work),
        ("supply", "energy_in_j", supply_work),
        ("drive", "force_n", -4.0),
        ("supply", "current_a", 0.25),
    ] {
        let actual = value(&json, &element(name), field);
        assert!(
            (actual - expected).abs() < 1e-11,
            "{name}.{field}: {actual} vs {expected}"
        );
    }
    assert_eq!(value(&json, "\"schema\"", "steps"), f64::from(steps));
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-11);
    assert!(json.contains("scalar commands carry no energy"));
    let text = std::fs::read_to_string(&path).unwrap();
    let mut native = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    assert_eq!(native.probe("drive.force_n").unwrap(), None);
    assert_eq!(native.probe("supply.energy_in_j").unwrap(), Some(0.0));
    scene_context(|cx, _| native.run(cx).unwrap());
    assert_eq!(native.render_json(), json);
}

#[test]
fn controlled_sources_have_adaptive_work_accounts_and_restore_the_accepted_prefix() {
    let text = std::fs::read_to_string(example("controlled_sources.json"))
        .unwrap()
        .replace(
            "\"step_s\": 0.375,",
            r#""step_s": 0.375,
  "adaptive": {"order": 2, "tolerances": {
    "block.position_m": {"absolute": 1e-7, "relative": 0},
    "block.velocity_m_s": {"absolute": 1e-7, "relative": 0},
    "capacitor.charge_c": {"absolute": 1e-7, "relative": 0},
    "drive.energy_in_j": {"absolute": 1e-7, "relative": 0},
    "supply.energy_in_j": {"absolute": 1e-7, "relative": 0}
  }},"#,
        );
    let path = scratch("adaptive_controlled_sources.json", &text);
    let mut uninterrupted = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    scene_context(|cx, _| {
        assert!(uninterrupted.advance(f64::MAX, cx).unwrap());
        assert!(uninterrupted.time_s() > 0.0 && !uninterrupted.is_finished());
        let (bytes, root) = uninterrupted.checkpoint().unwrap();
        resumed.resume_checkpoint(&bytes, root).unwrap();
        assert_eq!(resumed.render_json(), uninterrupted.render_json());
        uninterrupted.run(cx).unwrap();
        resumed.run(cx).unwrap();
    });
    assert_eq!(resumed.render_json(), uninterrupted.render_json());
    // Exact integrals of the authored commands, independent of the adaptive schedule.
    assert!((resumed.probe("block.velocity_m_s").unwrap().unwrap() - 0.5).abs() < 1e-10);
    assert!((resumed.probe("block.position_m").unwrap().unwrap() - 79.0 / 48.0).abs() < 1e-5);
    assert!((resumed.probe("capacitor.charge_c").unwrap().unwrap() - 0.3125).abs() < 1e-10);
    assert!((resumed.probe("drive.energy_in_j").unwrap().unwrap() - 0.75).abs() < 1e-10);
    assert!((resumed.probe("supply.energy_in_j").unwrap().unwrap() + 0.03515625).abs() < 1e-10);
    assert!(resumed.ledger().energy_residual_j().abs() < 1e-10);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(json, uninterrupted.render_json());

    let missing_work = text.replace(
        "\"drive.energy_in_j\": {\"absolute\": 1e-7, \"relative\": 0},",
        "",
    );
    let error = fs_cli::SystemScene::from_text(&missing_work, "missing-work.json")
        .err()
        .unwrap();
    assert!(error.message.contains("drive.energy_in_j"), "{error}");
}

#[test]
fn controlled_sources_refuse_unknown_domains_wrong_commands_and_missing_wiring() {
    let text = std::fs::read_to_string(example("controlled_sources.json")).unwrap();
    for (name, from, to, reason) in [
        (
            "domain.json",
            "\"domain\": \"translation\"",
            "\"domain\": \"temperature\"",
            "elements.drive: domain",
        ),
        (
            "field.json",
            "\"domain\": \"translation\"",
            "\"domain\": \"translation\", \"force_n\": 4",
            "unsupported field",
        ),
        (
            "effort-unit.json",
            "\"force_n\": {\"table\"",
            "\"voltage_v\": {\"table\"",
            "this port takes exactly one quantity, force_n",
        ),
        (
            "flow-unit.json",
            "\"current_a\": {\"table\"",
            "\"temperature_k\": {\"table\"",
            "this port takes exactly one quantity, current_a",
        ),
    ] {
        let changed = text.replace(from, to);
        assert_ne!(changed, text);
        let (code, message) = refused(name, &changed);
        assert_eq!(code, i32::from(fs_cli::exit::INPUT), "{message}");
        assert!(message.contains(reason), "{message}");
    }
    let disconnected = text.replace(
        "    {\"prescribe\": \"drive.command\", \"force_n\": {\"table\": [[0, 4], [0.25, 4], [0.75, -4], [1.25, -4]]}},\n",
        "",
    );
    assert_ne!(disconnected, text);
    let (code, message) = refused("missing-command.json", &disconnected);
    assert_eq!(code, i32::from(fs_cli::exit::REFUSED), "{message}");
    assert!(message.contains("drive.command"), "{message}");
}

#[test]
fn speed_controller_closes_the_loop_with_substep_invariant_sampling_and_signed_work() {
    let example = std::fs::read_to_string(example("closed_loop_speed.json")).unwrap();
    // The authored step is four controller periods; the second run takes four
    // plant substeps per period. Both must follow the same digital recurrence.
    for (name, step, substeps) in [
        ("speed_coarse.json", "0.5", 1_u32),
        ("speed_fine.json", "0.03125", 4_u32),
    ] {
        let text = example.replace("\"step_s\": 0.5", &format!("\"step_s\": {step}"));
        let path = scratch(name, &text);
        let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        let (mut angle, mut speed, mut integral, mut work) = (0.0_f64, 0.0, 0.0, 0.0);
        let mut supplied = false;
        let mut recovered = false;
        scene_context(|cx, _| {
            for sample in 0..120_u32 {
                let time = f64::from(sample) * 0.125;
                let target = if time < 2.0 {
                    10.0
                } else if time < 4.0 {
                    10.0 - 4.0 * (time - 2.0)
                } else if time < 8.0 {
                    2.0
                } else if time < 10.0 {
                    2.0 - 1.5 * (time - 8.0)
                } else {
                    -1.0
                };
                // Independent fixed-period PI and exact constant-torque inertia
                // recurrence; do not call the controller or plant owners here.
                let error = target - speed;
                let predicted = integral + 0.25 * error;
                let unlimited = error + predicted;
                let torque = unlimited.clamp(-2.0, 3.0);
                integral = predicted + (torque - unlimited) / 3.0;
                assert!((torque * speed).abs() < 100.0, "power limits are inactive");
                let next_speed = speed + 0.125 * (torque - 0.25);
                let displacement = 0.125 * 0.5 * (speed + next_speed);
                angle += displacement;
                speed = next_speed;
                let source_work = -torque * displacement;
                work += source_work;
                supplied |= source_work < -1e-6;
                recovered |= source_work > 1e-6;
                for _ in 0..substeps {
                    assert!(scene.advance(f64::MAX, cx).unwrap());
                }
                for (coordinate, expected) in [
                    ("rotor.angle_rad", angle),
                    ("rotor.angular_velocity_rad_s", speed),
                    ("control.integral_torque_nm", integral),
                    ("control.held_torque_nm", torque),
                    ("control.sample_count", f64::from(sample + 1)),
                    ("control.sample_phase_s", 0.0),
                    ("drive.energy_in_j", work),
                    ("load.energy_in_j", 0.25 * angle),
                ] {
                    let actual = scene.probe(coordinate).unwrap().unwrap();
                    assert!(
                        (actual - expected).abs() < 2e-10,
                        "{name} sample {sample}: {coordinate}: {actual} vs {expected}"
                    );
                }
                assert!(scene.ledger().energy_residual_j().abs() < 2e-10);
            }
        });
        assert!(
            supplied && recovered,
            "the physical drive supplies and recovers work"
        );
        assert!(scene.is_finished());
        assert_eq!(scene.steps(), u64::from(120 * substeps));
        assert!(
            (speed + 1.0).abs() < 0.08,
            "reversal settles toward the command"
        );
        assert!((0.5 * speed * speed + work + 0.25 * angle).abs() < 2e-10);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        assert_eq!(json, scene.render_json());
    }
}

#[test]
fn speed_controller_refuses_invalid_parameters_units_and_adaptive_policy() {
    let text = std::fs::read_to_string(example("closed_loop_speed.json")).unwrap();
    for (from, to, reason) in [
        (
            "\"publish_state\": true",
            "\"publish_state\": 1",
            "publish_state: expected a boolean",
        ),
        (
            "\"sample_period_s\": 0.125",
            "\"sample_period_s\": 0",
            "sample period must be finite and positive",
        ),
        (
            "\"relative_speed_rad_s\":",
            "\"angular_velocity_rad_s\":",
            "this port takes exactly one quantity, relative_speed_rad_s",
        ),
        (
            "\"duration_s\": 15.0,",
            "\"duration_s\": 15.0, \"adaptive\": {\"tolerances\": {}},",
            "discrete policy, not a continuous error coordinate",
        ),
    ] {
        let changed = text.replace(from, to);
        assert_ne!(changed, text);
        let error = fs_cli::SystemScene::from_text(&changed, "invalid-controller.json")
            .err()
            .unwrap();
        assert!(error.message.contains(reason), "{error}");
    }
}

#[test]
fn native_scene_runs_match_cli_and_incremental_mechanical_and_phase_runs() {
    for name in ["variable_transmission.json", "phase_change_store.json"] {
        let path = example(name);
        let mut bytes = std::fs::read(&path).unwrap();
        let mut whole = fs_cli::SystemScene::from_bytes(&bytes, &path).unwrap();
        let mut incremental =
            fs_cli::SystemScene::from_text(std::str::from_utf8(&bytes).unwrap(), &path).unwrap();
        assert_eq!(whole.source_path(), std::path::Path::new(&path));
        assert_eq!(whole.time_s(), 0.0);
        assert_eq!(whole.steps(), 0);
        assert!(!whole.is_finished());
        if name == "variable_transmission.json" {
            assert_eq!(whole.probe("transmission.effort_a").unwrap(), None);
            assert_eq!(
                whole.probe("rotor.angular_velocity_rad_s").unwrap(),
                Some(1.0)
            );
        }
        // All runtime input is owned after admission; caller bytes are reusable.
        bytes.fill(b'?');
        scene_context(|cx, _| {
            whole.run(cx).unwrap();
            while incremental.advance(f64::MAX, cx).unwrap() {}
            let finished = incremental.render_json();
            incremental.run(cx).unwrap();
            assert!(!incremental.advance(1.0, cx).unwrap());
            assert_eq!(incremental.render_json(), finished);
        });
        let (code, json, stderr) = run_binary(&["--json", "system", &path]);
        assert_eq!(code, 0, "{name}: {stderr}");
        assert_eq!(whole.render_json(), json, "{name}");
        assert_eq!(incremental.render_json(), json, "{name}");
        assert_eq!(whole.render_text(), incremental.render_text());
        assert_eq!(whole.time_s(), whole.duration_s());
        assert!(whole.is_finished());
        if name == "phase_change_store.json" {
            assert!(
                (whole
                    .probe("store.specific_enthalpy_j_kg")
                    .unwrap()
                    .unwrap()
                    - 240_000.0)
                    .abs()
                    < 1e-7
            );
            assert!(
                (whole.probe("store.liquid_mass_fraction").unwrap().unwrap() - 1.0).abs() < 1e-12
            );
            assert!(whole.ledger().energy_residual_j().abs() < 1e-7);
        }
    }
}

#[test]
fn native_adaptive_scene_resumes_the_same_step_and_sample_policy_as_cli() {
    let path = example("adaptive_oscillator.json");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    scene_context(|cx, _| {
        for _ in 0..3 {
            assert!(scene.advance(f64::MAX, cx).unwrap());
        }
        assert!(scene.time_s() > 0.0 && scene.time_s() < scene.duration_s());
        scene.run(cx).unwrap();
        let finished = scene.render_json();
        scene.run(cx).unwrap();
        assert_eq!(scene.render_json(), finished);
    });
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(scene.render_json(), json);
}

#[test]
fn native_scene_cancellation_preserves_the_initial_and_accepted_reports() {
    for adaptive in [false, true] {
        let text = phase_boundary_scene(0.0, 10_000.0);
        let text = if adaptive {
            text.replace(
                "\"step_s\": 0.25",
                r#""step_s": 0.25,
              "adaptive": {"tolerances": {
                "store.specific_enthalpy_j_kg": {"absolute": 1e-6, "relative": 0}
              }}"#,
            )
        } else {
            text
        };
        let path = std::path::Path::new("retained/native-scene.json");
        let mut scene = fs_cli::SystemScene::from_text(&text, path).unwrap();
        let initial = scene.render_json();
        scene_context(|cx, gate| {
            gate.request();
            let error = scene.advance(0.25, cx).unwrap_err();
            assert_eq!(error.exit_code, fs_cli::exit::CANCELLED);
            assert_eq!(scene.render_json(), initial);
        });
        scene_context(|cx, gate| {
            assert!(scene.advance(0.25, cx).unwrap());
            let prefix = scene.render_json();
            gate.request();
            let error = scene.run(cx).unwrap_err();
            assert_eq!(error.exit_code, fs_cli::exit::CANCELLED);
            assert_eq!(scene.render_json(), prefix);
        });
        let mut reference = fs_cli::SystemScene::from_text(&text, path).unwrap();
        scene_context(|cx, _| {
            scene.run(cx).unwrap();
            reference.run(cx).unwrap();
        });
        assert_eq!(scene.render_json(), reference.render_json());
        assert_eq!(scene.render_json().matches("{\"time_s\":").count(), 13);
    }
    let text = phase_boundary_scene(239_999.0, 10_000.0);
    let mut refused = fs_cli::SystemScene::from_text(&text, "retained/outside-chart.json").unwrap();
    let original = refused.render_json();
    scene_context(|cx, _| {
        let error = refused.advance(0.25, cx).unwrap_err();
        assert_eq!(error.exit_code, fs_cli::exit::REFUSED);
        assert!(error.message.contains("OutsideEnthalpyDomain"), "{error}");
        assert_eq!(refused.render_json(), original);
    });
}

#[test]
fn native_scene_admission_probes_and_caller_horizons_are_bounded() {
    let oversized = vec![b' '; 1024 * 1024 + 1];
    for bytes in [&b"{"[..], &b"\xff"[..], &oversized[..]] {
        let error = fs_cli::SystemScene::from_bytes(bytes, "never-read.json")
            .err()
            .unwrap();
        assert_eq!(error.code, "system-input");
        assert_eq!(error.exit_code, fs_cli::exit::INPUT);
    }
    let text =
        phase_boundary_scene(0.0, 10_000.0).replace("\"duration_s\": 24.0", "\"duration_s\": 0.3");
    for adaptive in [false, true] {
        let text = if adaptive {
            text.replace(
                "\"step_s\": 0.25",
                r#""step_s": 0.25,
              "adaptive": {"min_step_s": 0.125, "tolerances": {
                "store.specific_enthalpy_j_kg": {"absolute": 1e-6, "relative": 0}
              }}"#,
            )
        } else {
            text.clone()
        };
        let mut scene = fs_cli::SystemScene::from_text(&text, "never-read.json").unwrap();
        for path in ["store", "missing.temperature_k", "store.unsupported"] {
            let error = scene.probe(path).unwrap_err();
            assert_eq!(error.code, "system-probe");
            assert!(error.message.contains(path));
        }
        let initial = scene.render_json();
        scene_context(|cx, _| {
            for interval in [0.0, -1.0, f64::NAN, f64::INFINITY] {
                assert_eq!(
                    scene.advance(interval, cx).unwrap_err().code,
                    "system-interval"
                );
                assert_eq!(scene.render_json(), initial);
            }
            // A caller ceiling smaller than the adaptive floor is a requested
            // physical horizon, not a failed accuracy retry.
            assert!(scene.advance(0.03125, cx).unwrap());
            assert_eq!(scene.time_s(), 0.03125);
            assert_eq!(
                scene.probe("store.specific_enthalpy_j_kg").unwrap(),
                Some(312.5)
            );
            while scene.advance(f64::MAX, cx).unwrap() {
                assert!(scene.time_s() <= 0.3);
            }
            assert_eq!(scene.time_s(), 0.3);
            assert!(
                (scene
                    .probe("store.specific_enthalpy_j_kg")
                    .unwrap()
                    .unwrap()
                    - 3000.0)
                    .abs()
                    < 1e-9
            );
        });
        assert!(scene.is_finished());
        assert_eq!(scene.render_json().matches("{\"time_s\":").count(), 2);
        assert!(scene.ledger().energy_residual_j().abs() < 1e-9);
    }
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
fn an_enthalpy_body_heats_through_latent_storage_without_an_apparent_capacity() {
    let (code, json, stderr) =
        run_binary(&["--json", "system", &example("phase_change_store.json")]);
    assert_eq!(code, 0, "{stderr}");
    // G1: h = Q t / m. The prescribed 200 kJ/kg latent interval consumes
    // heat for twenty seconds while the temperature stays at 300 K.
    for sample in json.split("{\"time_s\":").skip(1) {
        let time: f64 = sample.split_once(',').unwrap().0.parse().unwrap();
        let h = 10_000.0 * time;
        let temperature = if h <= 20_000.0 {
            290.0 + h / 2000.0
        } else if h <= 220_000.0 {
            300.0
        } else {
            300.0 + (h - 220_000.0) / 2000.0
        };
        let fraction = ((h - 20_000.0) / 200_000.0).clamp(0.0, 1.0);
        assert!((value(sample, "", "store.specific_enthalpy_j_kg") - h).abs() < 1e-7);
        assert!((value(sample, "", "store.temperature_k") - temperature).abs() < 1e-10);
        assert!((value(sample, "", "store.liquid_mass_fraction") - fraction).abs() < 1e-12);
    }
    assert_eq!(json.matches("{\"time_s\":").count(), 13);
    assert!((value(&json, &element("store"), "stored_change_j") - 240_000.0).abs() < 1e-7);
    assert!((value(&json, &element("heater"), "energy_in_j") + 240_000.0).abs() < 1e-7);
    assert!((value(&json, &element("store"), "bulk_density_kg_m3") - 890.0).abs() < 1e-10);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-7);
    assert!(json.contains("authored bounded sensible and latent storage"));
    assert!(!json.contains("phase change or step-size error control"));
}

fn phase_boundary_scene(initial_h: f64, heat_w: f64) -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 24.0, "step_s": 0.25,
  "sample_every_s": 2.0,
  "elements": {{"store": {{
    "kind": "enthalpy_body", "mass_kg": 1.0,
    "specific_enthalpy_j_kg": {initial_h}, "phase": "solid-liquid",
    "enthalpy_curve": [[0.0, 290.0, 0.0, 1000.0], [20000.0, 300.0, 0.0, 990.0],
                       [220000.0, 300.0, 1.0, 900.0], [240000.0, 310.0, 1.0, 890.0]]
  }}}},
  "connections": [{{"prescribe": "store.heat", "heat_rate_w": {heat_w}}}]
}}"#
    )
}

#[test]
fn an_enthalpy_body_freezes_on_the_same_chart_and_exports_its_entropy_change() {
    let path = scratch(
        "freezing_store.json",
        &phase_boundary_scene(240_000.0, -10_000.0),
    );
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    assert!(value(&json, &element("store"), "specific_enthalpy_j_kg").abs() < 1e-7);
    assert!((value(&json, &element("store"), "temperature_k") - 290.0).abs() < 1e-10);
    assert!(value(&json, &element("store"), "liquid_mass_fraction").abs() < 1e-12);
    // G1/G3: cooling removes sensible entropy on both sides and L/T at
    // coexistence; no second latent source or fitted capacity appears.
    let entropy =
        2000.0 * (300.0_f64 / 290.0).ln() + 200_000.0 / 300.0 + 2000.0 * (310.0_f64 / 300.0).ln();
    let boundary = "\"port\":\"store.heat\"";
    assert!((value(&json, boundary, "energy_in_j") + 240_000.0).abs() < 1e-7);
    assert!((value(&json, boundary, "entropy_in_j_per_k") + entropy).abs() < 1e-8);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-7);
    let middle = json
        .split("{\"time_s\":")
        .skip(1)
        .find(|sample| sample.split_once(',').unwrap().0.parse::<f64>().unwrap() == 12.0)
        .unwrap();
    let expected_density = (0.5_f64 / 990.0 + 0.5 / 900.0).recip();
    assert!((value(middle, "", "store.bulk_density_kg_m3") - expected_density).abs() < 1e-9);
}

#[test]
fn enthalpy_bodies_in_a_conductive_loop_reach_the_energy_determined_equilibrium() {
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 200.0, "step_s": 0.5,
  "elements": {
    "hot": {"kind": "enthalpy_body", "mass_kg": 1.0, "specific_enthalpy_j_kg": 40000.0,
      "phase": "solid", "enthalpy_curve": [[0.0, 290.0, 0.0, 1000.0], [80000.0, 370.0, 0.0, 1000.0]]},
    "cold": {"kind": "enthalpy_body", "mass_kg": 2.0, "specific_enthalpy_j_kg": 0.0,
      "phase": "liquid", "enthalpy_curve": [[0.0, 290.0, 1.0, 900.0], [80000.0, 370.0, 1.0, 900.0]]},
    "contact": {"kind": "conductance", "conductance_w_k": 100.0}
  },
  "connections": [{"join": ["hot.heat", "contact.a"]}, {"join": ["cold.heat", "contact.b"]}]
}"#;
    let path = scratch("enthalpy_contact.json", scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    // Both authored charts have Cp=1000 J/(kg K). After thirty time
    // constants, the endpoint is within 1e-6 K of the independent limit.
    let equilibrium = (330.0 + 2.0 * 290.0) / 3.0;
    for name in ["hot", "cold"] {
        assert!((value(&json, &element(name), "temperature_k") - equilibrium).abs() < 1e-6);
    }
    let hot = value(&json, &element("hot"), "specific_enthalpy_j_kg");
    let cold = value(&json, &element("cold"), "specific_enthalpy_j_kg");
    assert!((hot + 2.0 * cold - 40_000.0).abs() < 1e-7);
    let entropy = 1000.0 * (equilibrium / 330.0_f64).ln() + 2000.0 * (equilibrium / 290.0_f64).ln();
    assert!((value(&json, "\"totals\"", "entropy_generated_j_per_k") - entropy).abs() < 1e-8);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-7);
}

#[test]
fn adaptive_enthalpy_storage_keeps_latent_energy_visible_and_refuses_derived_accounts() {
    let scene = phase_boundary_scene(0.0, 10_000.0)
        .replace("\"duration_s\": 24.0", "\"duration_s\": 12.0")
        .replace(
            "\"step_s\": 0.25",
            r#""step_s": 2.0,
  "adaptive": {"tolerances": {"store.specific_enthalpy_j_kg": {"absolute": 1e-6, "relative": 0}}}"#,
        );
    let path = scratch("adaptive_latent_store.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    assert!((value(&json, &element("store"), "specific_enthalpy_j_kg") - 120_000.0).abs() < 1e-7);
    assert!((value(&json, &element("store"), "temperature_k") - 300.0).abs() < 1e-10);
    assert!((value(&json, &element("store"), "liquid_mass_fraction") - 0.5).abs() < 1e-12);
    assert!(value(&json, "\"adaptive\"", "max_error_norm") <= 1.0);
    for field in [
        "temperature_k",
        "liquid_mass_fraction",
        "bulk_density_kg_m3",
    ] {
        let changed = scene.replace("store.specific_enthalpy_j_kg", &format!("store.{field}"));
        let (code, stderr) = refused(&format!("adaptive_{field}.json"), &changed);
        assert_eq!(code, 3, "{stderr}");
        assert!(stderr.contains(&format!("store.{field}")), "{stderr}");
    }
}

#[test]
fn enthalpy_scenes_refuse_bad_charts_and_source_extrapolation_by_name() {
    let scene = phase_boundary_scene(0.0, 10_000.0);
    for (name, from, to, expected) in [
        ("mass", "\"mass_kg\": 1.0", "\"mass_kg\": 0.0", "mass_kg"),
        (
            "phase",
            "\"phase\": \"solid-liquid\"",
            "\"phase\": \"gas\"",
            "phase",
        ),
        (
            "shape",
            "[0.0, 290.0, 0.0, 1000.0]",
            "[0.0, 290.0, 0.0]",
            "enthalpy_curve[0]",
        ),
        (
            "order",
            "[20000.0, 300.0",
            "[0.0, 300.0",
            "strictly increasing",
        ),
        (
            "temperature",
            "[20000.0, 300.0",
            "[20000.0, 280.0",
            "temperature must be nondecreasing",
        ),
        (
            "density",
            "[0.0, 290.0, 0.0, 1000.0]",
            "[0.0, 290.0, 0.0, 0.0]",
            "positive temperature/density",
        ),
        (
            "extra",
            "\"mass_kg\": 1.0",
            "\"mass_kg\": 1.0, \"temperature_k\": 290.0",
            "unsupported field",
        ),
    ] {
        assert!(scene.contains(from));
        let (code, stderr) = refused(&format!("enthalpy_{name}.json"), &scene.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(expected), "{name}: {stderr}");
    }
    let (code, stderr) = refused(
        "enthalpy_initial_support.json",
        &phase_boundary_scene(-1.0, 0.0),
    );
    assert_eq!(code, 3, "{stderr}");
    assert!(stderr.contains("specific_enthalpy_j_kg") && stderr.contains("OutsideEnthalpyDomain"));
    let (code, stderr) = refused(
        "enthalpy_step_support.json",
        &phase_boundary_scene(0.0, 10_001.0),
    );
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("OutsideEnthalpyDomain"), "{stderr}");
}

fn prescribed_layer_scene(geometry: &str, a: &str, b: &str, resistance: f64) -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 2.0, "step_s": 0.5,
  "elements": {{"wall": {{
    "kind": "layer_conductor", "conductivity_a": {a}, "conductivity_b": {b},
    "geometry": {geometry}, "contact_resistance_m2_k_w": {resistance}
  }}}},
  "connections": [{{"prescribe": "wall.a", "temperature_k": 400.0}},
                  {{"prescribe": "wall.b", "temperature_k": 300.0}}]
}}"#
    )
}

#[test]
fn layer_conductors_retain_planar_radial_and_contact_series_resistances() {
    let pi = std::f64::consts::PI;
    let contact = 0.03;
    // G1: independent Fourier integrals dr/(k A(r)), plus R''/A_contact.
    let cases = [
        (
            "planar",
            r#"{"kind":"planar","distances_m":[0.1,0.2],"area_m2":0.5}"#,
            0.1 / (2.0 * 0.5) + 0.2 / (4.0 * 0.5) + contact / 0.5,
        ),
        (
            "cylindrical",
            r#"{"kind":"cylindrical","radii_m":[0.01,0.04],"contact_radius_m":0.02,"length_m":0.3}"#,
            2.0_f64.ln() / (2.0 * pi * 0.3 * 2.0)
                + 2.0_f64.ln() / (2.0 * pi * 0.3 * 4.0)
                + contact / (2.0 * pi * 0.02 * 0.3),
        ),
        (
            "spherical",
            r#"{"kind":"spherical","radii_m":[0.01,0.04],"contact_radius_m":0.02}"#,
            (100.0 - 50.0) / (4.0 * pi * 2.0)
                + (50.0 - 25.0) / (4.0 * pi * 4.0)
                + contact / (4.0 * pi * 0.02 * 0.02),
        ),
    ];
    for (name, geometry, resistance) in cases {
        let scene =
            prescribed_layer_scene(geometry, "[[250,2],[450,2]]", "[[250,4],[450,4]]", contact);
        let path = scratch(&format!("layer_{name}.json"), &scene);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let expected = 100.0 / resistance;
        let actual = value(&json, &element("wall"), "heat_rate_w");
        assert!(
            (actual - expected).abs() < 1e-8 * expected,
            "{name}: {actual} vs {expected}"
        );
        assert!((value(&json, "\"port\":\"wall.a\"", "energy_in_j") - 2.0 * expected).abs() < 1e-7);
        assert!((value(&json, "\"port\":\"wall.b\"", "energy_in_j") + 2.0 * expected).abs() < 1e-7);
        let entropy = 2.0 * expected * (1.0 / 300.0 - 1.0 / 400.0);
        assert!((value(&json, "\"totals\"", "entropy_generated_j_per_k") - entropy).abs() < 1e-9);
        assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-9);
    }
}

#[test]
fn layer_conductors_integrate_affine_conductivity_instead_of_freezing_it() {
    let geometry = r#"{"kind":"planar","distances_m":[0.1,0.2],"area_m2":0.5}"#;
    let curve = "[[250,1.5],[450,3.5]]"; // k(T) = 0.01 T - 1.
    let nonlinear = prescribed_layer_scene(geometry, curve, curve, 0.0).replace(
        "\"contact_resistance_m2_k_w\": 0",
        r#""contact_resistance_m2_k_w": 0,
        "local_solver": {"max_iterations": 100, "absolute_heat_tolerance_w": 1e-10,
                         "relative_heat_tolerance": 1e-12}"#,
    );
    let constant = prescribed_layer_scene(geometry, "[[250,2],[450,2]]", "[[250,2],[450,2]]", 0.0);
    let mut powers = Vec::new();
    for (name, scene) in [("affine", nonlinear), ("constant", constant)] {
        let path = scratch(&format!("layer_{name}.json"), &scene);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        powers.push(value(&json, &element("wall"), "heat_rate_w"));
    }
    // G1: A/(La+Lb) integral_300^400 (0.01 T - 1) dT.
    let expected = 0.5 / 0.3 * (0.005 * (400.0 * 400.0 - 300.0 * 300.0) - 100.0);
    assert!(
        (powers[0] - expected).abs() < 1e-7,
        "{} vs {expected}",
        powers[0]
    );
    assert!((powers[1] - 0.5 / 0.3 * 2.0 * 100.0).abs() < 1e-7);
    assert!(powers[0] > 1.2 * powers[1]);
}

#[test]
fn layer_scenes_refuse_invalid_laws_geometry_and_local_budgets_by_name() {
    let geometry = r#"{"kind":"planar","distances_m":[0.1,0.2],"area_m2":0.5}"#;
    let scene = prescribed_layer_scene(geometry, "[[250,2],[450,2]]", "[[250,4],[450,4]]", 0.03);
    for (name, from, to, expected) in [
        ("short", "[[250,2],[450,2]]", "[[250,2]]", "conductivity_a"),
        (
            "shape",
            "[[250,2],[450,2]]",
            "[[250,2,3],[450,2]]",
            "conductivity_a[0]",
        ),
        (
            "order",
            "[[250,2],[450,2]]",
            "[[450,2],[250,2]]",
            "conductivity_a",
        ),
        (
            "zero_k",
            "[[250,4],[450,4]]",
            "[[250,0],[450,4]]",
            "conductivity_b",
        ),
        (
            "negative_t",
            "[[250,4],[450,4]]",
            "[[-250,4],[450,4]]",
            "conductivity_b",
        ),
        ("distance", "[0.1,0.2]", "[0.0,0.2]", "geometry"),
        ("pair", "[0.1,0.2]", "[0.1]", "geometry.distances_m"),
        ("area", "\"area_m2\":0.5", "\"area_m2\":0", "area_m2"),
        (
            "kind",
            "\"kind\":\"planar\"",
            "\"kind\":\"cone\"",
            "geometry.kind",
        ),
        (
            "extra",
            "\"area_m2\":0.5",
            "\"area_m2\":0.5,\"radius_m\":0.1",
            "unsupported field",
        ),
        (
            "contact",
            "\"contact_resistance_m2_k_w\": 0.03",
            "\"contact_resistance_m2_k_w\": -1",
            "contact_resistance_m2_k_w",
        ),
        (
            "iterations",
            "\"contact_resistance_m2_k_w\": 0.03",
            "\"local_solver\":{\"max_iterations\":1001}",
            "local_solver.max_iterations",
        ),
        (
            "tolerance",
            "\"contact_resistance_m2_k_w\": 0.03",
            "\"local_solver\":{\"relative_heat_tolerance\":-1}",
            "local_solver.relative_heat_tolerance",
        ),
    ] {
        assert!(scene.contains(from));
        let (code, stderr) = refused(
            &format!("invalid_layer_{name}.json"),
            &scene.replace(from, to),
        );
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(
            stderr.contains("wall") && stderr.contains(expected),
            "{name}: {stderr}"
        );
    }
    let radial = r#"{"kind":"spherical","radii_m":[0.03,0.04],"contact_radius_m":0.02}"#;
    let (code, stderr) = refused("layer_radial_order.json", &scene.replace(geometry, radial));
    assert_eq!(code, 3, "{stderr}");
    assert!(
        stderr.contains("wall.geometry") && stderr.contains("straddle"),
        "{stderr}"
    );
    let (code, stderr) = refused(
        "layer_support.json",
        &scene.replace("\"temperature_k\": 400.0", "\"temperature_k\": 500.0"),
    );
    assert_eq!(code, 4, "{stderr}");
    assert!(
        stderr.contains("wall") && stderr.contains("MaterialDomain"),
        "{stderr}"
    );
}

#[test]
fn nonlinear_layers_heat_a_finite_phase_change_store_and_close_every_sample_energy() {
    let (code, json, stderr) =
        run_binary(&["--json", "system", &example("phase_change_conduction.json")]);
    assert_eq!(code, 0, "{stderr}");
    let mut previous_h = -1.0;
    let mut plateau_samples = 0;
    for sample in json.split("{\"time_s\":").skip(1) {
        let hot = value(sample, "", "hot.temperature_k");
        let h = value(sample, "", "store.specific_enthalpy_j_kg");
        let fraction = value(sample, "", "store.liquid_mass_fraction");
        // G3: neither the quasisteady layers nor their contact store heat.
        assert!((2000.0 * (hot - 400.0) + h).abs() < 1e-5);
        assert!(h >= previous_h && fraction <= 0.9 + 1e-10);
        if h > 20_000.0 {
            plateau_samples += 1;
            assert!((value(sample, "", "store.temperature_k") - 300.0).abs() < 1e-10);
            assert!((fraction - (h - 20_000.0) / 200_000.0).abs() < 1e-12);
        }
        previous_h = h;
    }
    assert!(plateau_samples >= 10);
    let fraction = value(&json, &element("store"), "liquid_mass_fraction");
    assert!(fraction > 0.75 && fraction < 0.9, "{fraction}");
    let hot_change = value(&json, &element("hot"), "stored_change_j");
    let pcm_change = value(&json, &element("store"), "stored_change_j");
    assert!((hot_change + pcm_change).abs() < 1e-5);
    assert_eq!(value(&json, &element("wall"), "stored_change_j"), 0.0);
    assert!(value(&json, &element("wall"), "heat_rate_w") > 0.0);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-5);
    assert!(value(&json, "\"totals\"", "entropy_generated_j_per_k") > 0.0);
    assert!(json.contains("authored bounded k(T) charts"));
    assert!(json.contains("no verified material-pack provenance"));
}

#[test]
fn adaptive_phase_change_networks_do_not_require_an_account_for_stateless_layers() {
    let scene = std::fs::read_to_string(example("phase_change_conduction.json"))
        .unwrap()
        .replace("\"duration_s\": 3000.0", "\"duration_s\": 300.0")
        .replace(
            "\"step_s\": 5.0",
            r#""step_s": 5.0,
  "adaptive": {"tolerances": {
    "hot.temperature_k": {"absolute": 0.001, "relative": 1e-6},
    "store.specific_enthalpy_j_kg": {"absolute": 1.0, "relative": 1e-6}
  }}"#,
        );
    let path = scratch("adaptive_phase_change_layers.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    assert!(value(&json, "\"adaptive\"", "max_error_norm") <= 1.0);
    assert!(value(&json, &element("store"), "liquid_mass_fraction") > 0.0);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-5);
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

/// One parallel node with storage, a resistor to ground and a supply, in
/// electrical or hydraulic names: the same numbers in either.
fn parallel_rc(hydraulic: bool) -> String {
    let [
        storage,
        branch,
        source,
        capacity,
        content,
        resistance,
        flow,
        ground,
    ] = if hydraulic {
        [
            "accumulator",
            "hydraulic_resistor",
            "flow_source",
            "capacitance_m3_pa",
            "volume_m3",
            "resistance_pa_s_m3",
            "flow_m3_s",
            "pressure_pa",
        ]
    } else {
        [
            "node_capacitor",
            "resistor_branch",
            "current_source",
            "capacitance_f",
            "charge_c",
            "resistance_ohm",
            "current_a",
            "voltage_v",
        ]
    };
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 2.0, "step_s": 0.001,
  "elements": {{
    "store": {{"kind": "{storage}", "{capacity}": 1e-9, "{content}": 2e-4}},
    "drain": {{"kind": "{branch}", "{resistance}": 1e9}},
    "supply": {{"kind": "{source}", "{flow}": -1e-4}},
    "oil": {{"kind": "heat_capacity", "capacity_j_k": 50.0, "temperature_k": 300.0}}
  }},
  "connections": [
    {{"common_effort": ["store.p", "drain.a", "supply.p"]}},
    {{"prescribe": "drain.b", "{ground}": 0}},
    {{"route": "drain.heat", "to": [["oil.heat", 1.0]]}}
  ]
}}"#
    )
}

#[test]
fn an_accumulator_drains_through_a_hydraulic_resistor_as_the_electrical_node_does() {
    // The hydraulic scene is the example, run where it lies.
    let hydraulic = example("hydraulic_accumulator.json");
    let committed = std::fs::read_to_string(&hydraulic).unwrap();
    assert_eq!(committed.trim_end(), parallel_rc(true));
    let electrical = scratch("electrical_rc.json", &parallel_rc(false));
    let (code, json, stderr) = run_binary(&["--json", "system", &hydraulic]);
    assert_eq!(code, 0, "{stderr}");
    let (code, twin, stderr) = run_binary(&["--json", "system", electrical.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");

    // C p' = J - p / R with J = 1e-4 m^3/s supplied, by the midpoint rule.
    let (c, r, supplied, dt) = (1.0e-9_f64, 1.0e9_f64, 1.0e-4_f64, 0.001_f64);
    let (mut pressure, mut work_j, mut heat_j) = (2.0e-4 / c, 0.0_f64, 0.0_f64);
    for _ in 0..2_000 {
        let half = 0.5 * dt;
        let mean = (pressure + half / c * supplied) / (1.0 + half / (r * c));
        pressure = 2.0 * mean - pressure;
        work_j += dt * supplied * mean;
        heat_j += dt * mean * mean / r;
    }
    let store = element("store");
    let reached = value(&json, &store, "pressure_pa");
    let volume = value(&json, &store, "volume_m3");
    assert!(
        (reached - pressure).abs() <= 1.0e-9 * pressure,
        "{reached} vs {pressure}"
    );
    assert!(
        (volume - c * pressure).abs() <= 1.0e-9 * c * pressure,
        "{volume}"
    );
    let delivered = -value(&json, &element("supply"), "energy_in_j");
    let dissipated = value(&json, &element("drain"), "dissipated_j");
    let warmed = value(&json, &element("oil"), "stored_change_j");
    assert!(
        (delivered - work_j).abs() <= 1.0e-9 * work_j,
        "{delivered} vs {work_j}"
    );
    assert!(
        (dissipated - heat_j).abs() <= 1.0e-9 * heat_j,
        "{dissipated} vs {heat_j}"
    );
    assert!(
        (warmed - heat_j).abs() <= 1.0e-9 * heat_j,
        "{warmed} vs {heat_j}"
    );
    let stored_change = (volume * volume - 2.0e-4 * 2.0e-4) / (2.0 * c);
    assert!((delivered - stored_change - warmed).abs() <= 1.0e-9 * heat_j);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-9 * heat_j);

    // The hydraulic names carry the electrical node's numbers, bit for bit.
    for (anchor, hydraulic_key, electrical_key) in [
        ("store", "volume_m3", "charge_c"),
        ("store", "pressure_pa", "voltage_v"),
        ("drain", "flow_m3_s", "current_a"),
        ("drain", "heat_w", "heat_w"),
        ("drain", "dissipated_j", "dissipated_j"),
        ("supply", "energy_in_j", "energy_in_j"),
        ("oil", "temperature_k", "temperature_k"),
    ] {
        let a = value(&json, &element(anchor), hydraulic_key);
        let b = value(&twin, &element(anchor), electrical_key);
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "{anchor}.{hydraulic_key} {a} vs {b}"
        );
    }

    let scene = parallel_rc(true);
    for (file, original, replacement, expected) in [
        (
            "accumulator_zero.json",
            "\"capacitance_m3_pa\": 1e-9",
            "\"capacitance_m3_pa\": 0",
            "positive",
        ),
        (
            "resistor_zero.json",
            "\"resistance_pa_s_m3\": 1e9",
            "\"resistance_pa_s_m3\": 0",
            "positive",
        ),
        (
            "accumulator_farads.json",
            "\"capacitance_m3_pa\"",
            "\"capacitance_f\"",
            "capacitance_f",
        ),
        (
            "ground_volts.json",
            "\"pressure_pa\": 0",
            "\"voltage_v\": 0",
            "exactly one quantity, pressure_pa",
        ),
    ] {
        let changed = scene.replace(original, replacement);
        assert_ne!(changed, scene);
        let (code, stderr) = refused(file, &changed);
        assert_eq!(code, 3, "{stderr}");
        assert!(stderr.contains(expected), "{stderr}");
    }
}

const GAS_CONSTANT_J_MOL_K: f64 = 8.314_462_618_153_24;

/// exp(E/R (1/T - 1/T_ref)): the Arrhenius factor on the viscosity.
fn arrhenius(activation_j_mol: f64, reference_k: f64, temperature_k: f64) -> f64 {
    (activation_j_mol / GAS_CONSTANT_J_MOL_K * (1.0 / temperature_k - 1.0 / reference_k)).exp()
}

/// Laminar flow of a Carreau-Yasuda liquid in a pipe by the
/// Rabinowitsch-Mooney integral, Q = pi r^3 / s_w^3 int_0^s_w s^2 g(s) ds:
/// the shear rate g at each stress s found by bisection on s = mu(g) g,
/// Simpson's rule over the stress. Independent of the fs-material solver.
fn carreau_flow(
    law: [f64; 5],
    factor: f64,
    rate_too: bool,
    radius_m: f64,
    length_m: f64,
    drop_pa: f64,
) -> f64 {
    let [mu0, mu_inf, time, n, a] = law;
    let stress = |rate: f64| {
        let reduced = if rate_too { factor * rate } else { rate };
        factor
            * (mu_inf + (mu0 - mu_inf) * (1.0 + (time * reduced).powf(a)).powf((n - 1.0) / a))
            * rate
    };
    let rate_at = |target: f64| {
        let mut high = 1.0;
        while stress(high) < target {
            high *= 2.0;
        }
        let mut low = 0.0;
        for _ in 0..200 {
            let middle = 0.5 * (low + high);
            if stress(middle) < target {
                low = middle;
            } else {
                high = middle;
            }
        }
        0.5 * (low + high)
    };
    let wall = drop_pa * radius_m / (2.0 * length_m);
    let panels = 4_000_u32;
    let h = wall / f64::from(panels);
    let sum: f64 = (0..=panels)
        .map(|i| {
            let s = h * f64::from(i);
            let weight = if i == 0 || i == panels {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            weight * s * s * rate_at(s)
        })
        .sum();
    std::f64::consts::PI * radius_m.powi(3) / wall.powi(3) * sum * h / 3.0
}

/// A `carreau_pipe` of radius 5 mm and length 1 m from a pressure supply to
/// a drain held at 0 Pa; `oil` adds elements and `heat` connects its heat
/// port. Arrhenius shift of 30 kJ/mol about 300 K.
fn pipe_scene(
    law: [f64; 5],
    mode: &str,
    drop_pa: f64,
    ceiling: f64,
    oil: &str,
    heat: &str,
) -> String {
    let [mu0, mu_inf, time, n, a] = law;
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 0.1, "step_s": 0.05,
  "elements": {{
    "supply": {{"kind": "pressure", "pressure_pa": {drop_pa}}},
    "pipe": {{"kind": "carreau_pipe",
      "viscosity_zero_pa_s": {mu0}, "viscosity_infinite_pa_s": {mu_inf},
      "time_constant_s": {time}, "power_index": {n}, "transition": {a},
      "rate_support_per_s": [0, 100000],
      "shift": {{"arrhenius": {{"activation_energy_j_mol": 30000, "reference_k": 300, "valid_k": [280, 360]}}}},
      "shift_mode": "{mode}",
      "radius_m": 0.005, "length_m": 1.0, "maximum_wall_rate_per_s": {ceiling}}}{oil}
  }},
  "connections": [
    {{"common_effort": ["supply.p", "pipe.a"]}},
    {{"prescribe": "pipe.b", "pressure_pa": 0}},
    {heat}
  ]
}}"#
    )
}

const AT_320_K: &str = r#"{"prescribe": "pipe.heat", "temperature_k": 320}"#;

#[test]
fn a_newtonian_carreau_pipe_at_a_shifted_temperature_passes_the_poiseuille_flow() {
    let scene = pipe_scene(
        [0.5, 0.5, 0.5, 0.4, 2.0],
        "vertical",
        1.0e4,
        5.0e4,
        "",
        AT_320_K,
    );
    let path = scratch("poiseuille.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    // Equal plateaus: a constant viscosity 0.5 Pa s at 300 K, times a_T.
    let viscosity = 0.5 * arrhenius(30_000.0, 300.0, 320.0);
    let radius: f64 = 0.005;
    let flow = std::f64::consts::PI * radius.powi(4) * 1.0e4 / (8.0 * viscosity * 1.0);
    let pipe = element("pipe");
    let reached = value(&json, &pipe, "flow_m3_s");
    assert!(
        (reached - flow).abs() <= 1.0e-9 * flow,
        "{reached} vs Poiseuille {flow}"
    );
    assert!((value(&json, &pipe, "heat_w") - flow * 1.0e4).abs() <= 1.0e-9 * flow * 1.0e4);
    assert!((value(&json, &pipe, "temperature_k") - 320.0).abs() <= 0.0);
    // All the supply's work leaves as heat through the held heat port.
    let delivered = -value(&json, &element("supply"), "energy_in_j");
    let expected = 0.1 * flow * 1.0e4;
    assert!(
        (delivered - expected).abs() <= 1.0e-9 * expected,
        "{delivered} vs {expected}"
    );
    let left = -value(&json, "\"port\":\"pipe.heat\"", "energy_in_j");
    assert!(
        (left - expected).abs() <= 1.0e-9 * expected,
        "{left} vs {expected}"
    );
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-9 * expected);

    // A ceiling below the wall rate (4 Q / pi r^3, about 106 1/s) refuses the step.
    let ceiling = 4.0 * flow / (std::f64::consts::PI * radius.powi(3));
    assert!(ceiling > 100.0 && ceiling < 110.0, "{ceiling}");
    let low = pipe_scene(
        [0.5, 0.5, 0.5, 0.4, 2.0],
        "vertical",
        1.0e4,
        50.0,
        "",
        AT_320_K,
    );
    let (code, stderr) = refused("pipe_over_ceiling.json", &low);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("wall-rate ceiling"), "{stderr}");
}

#[test]
fn a_shear_thinning_carreau_pipe_follows_the_rabinowitsch_mooney_flow_in_both_shift_modes() {
    let law = [2.0, 0.01, 0.5, 0.4, 2.0];
    let factor = arrhenius(30_000.0, 300.0, 320.0);
    let mut flows = Vec::new();
    for (mode, rate_too) in [("vertical", false), ("vertical-and-rate", true)] {
        let scene = pipe_scene(law, mode, 1.0e4, 5.0e4, "", AT_320_K);
        let path = scratch(&format!("thinning_{mode}.json"), &scene);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        let expected = carreau_flow(law, factor, rate_too, 0.005, 1.0, 1.0e4);
        let reached = value(&json, &element("pipe"), "flow_m3_s");
        assert!(
            (reached - expected).abs() <= 1.0e-7 * expected,
            "{mode}: {reached} vs {expected}"
        );
        let heat = value(&json, &element("pipe"), "heat_w");
        assert!(
            (heat - reached * 1.0e4).abs() <= 1.0e-9 * heat,
            "{mode}: {heat}"
        );
        flows.push(reached);
    }
    // The shift's rate factor (a_T < 1 at 320 K) moves the curve's knee,
    // so the two modes differ by much more than the tolerance.
    assert!(flows[0] > 1.5 * flows[1], "{flows:?}");

    let scene = pipe_scene(law, "vertical", 1.0e4, 5.0e4, "", AT_320_K);
    for (file, original, replacement, expected) in [
        (
            "pipe_no_shift.json",
            "\"shift\": {\"arrhenius\": {\"activation_energy_j_mol\": 30000, \"reference_k\": 300, \"valid_k\": [280, 360]}},",
            "",
            "missing shift",
        ),
        (
            "pipe_mode.json",
            "\"shift_mode\": \"vertical\"",
            "\"shift_mode\": \"both\"",
            "shift_mode",
        ),
        (
            "pipe_index.json",
            "\"power_index\": 0.4",
            "\"power_index\": 1.5",
            "Carreau-Yasuda",
        ),
        (
            "pipe_support.json",
            "[0, 100000]",
            "[1, 100000]",
            "calibrated reference curve",
        ),
        (
            "pipe_radius.json",
            "\"radius_m\": 0.005",
            "\"radius_m\": 0",
            "radius_m",
        ),
        (
            "pipe_field.json",
            "\"length_m\": 1",
            "\"diameter_m\": 0.01, \"length_m\": 1",
            "diameter_m",
        ),
    ] {
        let changed = scene.replace(original, replacement);
        assert_ne!(changed, scene, "{file}");
        let (code, stderr) = refused(file, &changed);
        assert_eq!(code, 3, "{stderr}");
        assert!(stderr.contains(expected), "{file}: {stderr}");
    }
    // Outside the shift's band the step refuses.
    let hot = scene.replace("\"temperature_k\": 320", "\"temperature_k\": 370");
    let (code, stderr) = refused("pipe_too_hot.json", &hot);
    assert_eq!(code, 4, "{stderr}");
}

#[test]
fn the_carreau_pipe_example_heats_its_oil_by_the_supply_work_and_thins_as_it_warms() {
    let (code, json, stderr) = run_binary(&["--json", "system", &example("carreau_pipe.json")]);
    assert_eq!(code, 0, "{stderr}");
    let law = [2.0, 0.01, 0.5, 0.4, 2.0];
    let pipe = element("pipe");
    // The last step's flow is the pipe's law at the temperature it read.
    let temperature = value(&json, &pipe, "temperature_k");
    let factor = arrhenius(30_000.0, 300.0, temperature);
    let expected = carreau_flow(law, factor, false, 0.005, 1.0, 1.0e4);
    let reached = value(&json, &pipe, "flow_m3_s");
    assert!(
        (reached - expected).abs() <= 1.0e-7 * expected,
        "{reached} vs {expected}"
    );
    // Warmer oil is thinner: the flow has risen above the 300 K flow.
    let cold = carreau_flow(law, 1.0, false, 0.005, 1.0, 1.0e4);
    assert!(
        temperature > 300.5 && reached > 1.05 * cold,
        "{temperature} K, {reached} vs {cold}"
    );
    // Every joule of the supply's work is the pipe's heat, held by the oil.
    let delivered = -value(&json, &element("supply"), "energy_in_j");
    let dissipated = value(&json, &pipe, "dissipated_j");
    let warmed = value(&json, &element("oil"), "stored_change_j");
    let oil_k = value(&json, &element("oil"), "temperature_k");
    assert!(delivered > 0.1, "{delivered}");
    assert!(
        (dissipated - delivered).abs() <= 1.0e-9 * delivered,
        "{dissipated} vs {delivered}"
    );
    assert!(
        (warmed - delivered).abs() <= 1.0e-9 * delivered,
        "{warmed} vs {delivered}"
    );
    assert!(
        (0.2 * (oil_k - 300.0) - warmed).abs() <= 1.0e-9 * warmed,
        "{oil_k}"
    );
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-9 * delivered);
    assert!(value(&json, "\"totals\"", "entropy_generated_j_per_k") > 0.0);
    assert!(json.contains("authored Carreau-Yasuda pipe laws"), "{json}");
    assert!(
        json.contains("no fluid inertia, compressibility, entrance loss"),
        "{json}"
    );
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
fn a_signal_modulated_transmission_matches_the_changing_ratio_rotation() {
    let (code, json, stderr) =
        run_binary(&["--json", "system", &example("variable_transmission.json")]);
    assert_eq!(code, 0, "{stderr}");
    let (inertia, stiffness, dt) = (0.02_f64, 400.0_f64, 1.0 / 256.0);
    let (initial_w, initial_q) = (1.0, 0.01);
    let (mut phase, mut angle, mut omega, mut extension) = (0.0_f64, 0.0_f64, initial_w, initial_q);
    let (mut last_ratio, mut last_effort, mut last_flow) = (0.0, 0.0, 0.0);
    for i in 0..512 {
        // Independently integrate the authored affine ratio on each interval.
        // The zero plateau disconnects the power exchange; the last ramp
        // reverses it. In energy coordinates each interval is a rotation.
        let midpoint = (f64::from(i) + 0.5) * dt;
        let ratio = if midpoint < 0.5 {
            0.02 * (1.0 - 2.0 * midpoint)
        } else if midpoint < 1.0 {
            0.0
        } else {
            -0.02 * (midpoint - 1.0)
        };
        phase += 2.0 * (0.5 * ratio * (stiffness / inertia).sqrt() * dt).atan();
        let next_w =
            initial_w * phase.cos() + initial_q * (stiffness / inertia).sqrt() * phase.sin();
        let next_q =
            initial_q * phase.cos() - initial_w * (inertia / stiffness).sqrt() * phase.sin();
        angle += dt * 0.5 * (omega + next_w);
        last_ratio = ratio;
        last_effort = ratio * stiffness * 0.5 * (extension + next_q);
        last_flow = ratio * 0.5 * (omega + next_w);
        (omega, extension) = (next_w, next_q);
    }
    for (name, field, expected) in [
        ("rotor", "angle_rad", angle),
        ("rotor", "angular_velocity_rad_s", omega),
        ("spring", "extension_m", extension),
        ("transmission", "ratio", last_ratio),
        ("transmission", "effort_a", last_effort),
        ("transmission", "flow_b", last_flow),
    ] {
        let actual = value(&json, &element(name), field);
        assert!(
            (actual - expected).abs() <= 1.0e-10,
            "{name}.{field}: {actual} vs {expected}"
        );
    }
    let reached_w = value(&json, &element("rotor"), "angular_velocity_rad_s");
    let reached_q = value(&json, &element("spring"), "extension_m");
    let energy = 0.5 * (inertia * reached_w.powi(2) + stiffness * reached_q.powi(2));
    assert!((energy - 0.03).abs() <= 1.0e-11);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-11);
    assert_eq!(value(&json, &element("transmission"), "dissipated_j"), 0.0);
    assert!(json.contains("no crank/cam geometry"), "{json}");
}

#[test]
fn a_modulated_transformer_needs_no_adaptive_state_account() {
    let scene = std::fs::read_to_string(example("variable_transmission.json"))
        .unwrap()
        .replace(
            "{\"table\": [[0, 0.02], [0.5, 0], [1, 0], [2, -0.02]]}",
            "0",
        )
        .replace(
            "\"step_s\": 0.00390625",
            r#""step_s": 0.125,
  "adaptive": {"order": 2, "tolerances": {
    "rotor.angle_rad": {"absolute": 1e-8, "relative": 0},
    "rotor.angular_velocity_rad_s": {"absolute": 1e-8, "relative": 0},
    "spring.extension_m": {"absolute": 1e-8, "relative": 0}
  }}"#,
        );
    let path = scratch("adaptive_transmission.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    assert!((value(&json, &element("rotor"), "angle_rad") - 2.0).abs() <= 1.0e-12);
    assert!((value(&json, &element("spring"), "extension_m") - 0.01).abs() <= 1.0e-12);
    assert_eq!(value(&json, &element("transmission"), "ratio"), 0.0);
    assert!(value(&json, "\"adaptive\"", "max_error_norm") <= 1.0);
}

#[test]
fn modulated_transformer_scenes_refuse_named_domains_fields_and_missing_ratios() {
    let scene = std::fs::read_to_string(example("variable_transmission.json")).unwrap();
    for (file, from, to, reason) in [
        (
            "transmission_domain.json",
            "\"a\": \"rotation\"",
            "\"a\": \"angle\"",
            "elements.transmission: a",
        ),
        (
            "transmission_field.json",
            "\"b\": \"translation\"",
            "\"b\": \"translation\", \"ratio\": 2",
            "unsupported field \"ratio\"",
        ),
        (
            "transmission_quantity.json",
            "\"ratio\": {\"table\"",
            "\"temperature_k\": {\"table\"",
            "this port takes exactly one quantity, ratio",
        ),
        (
            "transmission_waveform.json",
            "{\"table\": [[0, 0.02], [0.5, 0], [1, 0], [2, -0.02]]}",
            "\"fast\"",
            "connections[3].ratio",
        ),
    ] {
        let changed = scene.replace(from, to);
        assert_ne!(changed, scene);
        let (code, stderr) = refused(file, &changed);
        assert_eq!(code, 3, "{stderr}");
        assert!(stderr.contains(reason), "{stderr}");
    }
    let disconnected = scene.replace(
        ",\n    {\"prescribe\": \"transmission.ratio\", \"ratio\": {\"table\": [[0, 0.02], [0.5, 0], [1, 0], [2, -0.02]]}}",
        "",
    );
    assert_ne!(disconnected, scene);
    let (code, stderr) = refused("transmission_unconnected_ratio.json", &disconnected);
    assert_eq!(code, 4, "{stderr}");
    assert!(stderr.contains("transmission.ratio"), "{stderr}");
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

/// The closed form of `m q'' + c q' + k q = F` from rest, returned as
/// `[momentum, displacement]`, for port_hamiltonian_oscillator.json:
/// m = 0.25 kg, k = 400 N/m, c = 1.5 N s/m and F = 2 N.
fn held_force_oscillator(time_s: f64) -> [f64; 2] {
    let (mass, stiffness, damping, force) = (0.25_f64, 400.0_f64, 1.5_f64, 2.0_f64);
    let decay = damping / (2.0 * mass);
    let frequency = (stiffness / mass - decay * decay).sqrt();
    let envelope = (-decay * time_s).exp();
    let (sine, cosine) = (frequency * time_s).sin_cos();
    [
        force / frequency * envelope * sine,
        force / stiffness * (1.0 - envelope * (cosine + decay / frequency * sine)),
    ]
}

/// The work the held 2 N force has done, less the energy still stored.
fn held_force_heat([momentum, displacement]: [f64; 2]) -> f64 {
    2.0 * displacement - (2.0 * momentum * momentum + 200.0 * displacement * displacement)
}

#[test]
fn a_port_hamiltonian_oscillator_steps_a_held_force_exactly_and_heats_its_case_by_the_lost_work() {
    let (code, json, stderr) = run_binary(&[
        "--json",
        "system",
        &example("port_hamiltonian_oscillator.json"),
    ]);
    assert_eq!(code, 0, "{stderr}");
    assert!(stderr.is_empty(), "{stderr}");
    let (end, dt) = (0.5, 1.0e-3);
    let (momentum_scale, displacement_scale) = (0.05, 0.005);
    let exact = held_force_oscillator(end);
    let [momentum, displacement] = exact;
    let reported = |name: &str, key: &str| value(&json, &element(name), key);
    // The owner's map is exact for a held input: 500 steps agree with the
    // closed form to rounding (about 2e-16 on aarch64). The midpoint rule
    // is about 3e-6 m off (below).
    let x0 = reported("oscillator", "x0");
    let x1 = reported("oscillator", "x1");
    assert!(
        (x0 - momentum).abs() <= 1.0e-12 * momentum_scale,
        "{x0} vs {momentum}"
    );
    assert!(
        (x1 - displacement).abs() <= 1.0e-12 * displacement_scale,
        "{x1} vs {displacement}"
    );

    // A source's account is the energy that entered it: minus its work.
    let supplied = 2.0 * displacement;
    let heat = held_force_heat(exact);
    assert!(heat > 0.0);
    for (name, key, expected) in [
        ("shaker", "energy_in_j", -supplied),
        ("oscillator", "stored_change_j", supplied - heat),
        ("oscillator", "dissipated_j", heat),
        ("case", "stored_change_j", heat),
    ] {
        let actual = reported(name, key);
        assert!(
            (actual - expected).abs() <= 1.0e-12,
            "{name}.{key}: {actual} vs {expected}"
        );
    }
    let temperature = reported("case", "temperature_k");
    assert!((temperature - 300.0 - heat / 10.0).abs() <= 1.0e-11);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() <= 1.0e-12);

    // Last-step means: the flow the port returned and the heat it sent.
    let before = held_force_oscillator(end - dt);
    let velocity = (displacement - before[1]) / dt;
    let heat_w = (heat - held_force_heat(before)) / dt;
    let mean = reported("oscillator", "p.velocity_m_s");
    assert!(
        (mean - velocity).abs() <= 1.0e-8 * velocity.abs(),
        "{mean} vs {velocity}"
    );
    let rate = reported("oscillator", "heat_w");
    assert!(
        (rate - heat_w).abs() <= 1.0e-6 * heat_w,
        "{rate} vs {heat_w}"
    );

    // The same law built from kernel elements takes the midpoint step: the
    // same motion, but not exact at this step.
    let elements = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.5, "step_s": 1.0e-3,
  "elements": {
    "block": {"kind": "mass", "mass_kg": 0.25},
    "spring": {"kind": "spring", "stiffness_n_m": 400.0},
    "damper": {"kind": "damper", "damping_n_s_m": 1.5},
    "shaker": {"kind": "force", "force_n": 2.0},
    "case": {"kind": "heat_capacity", "capacity_j_k": 10.0, "temperature_k": 300.0}
  },
  "connections": [
    {"join": ["block.p", "spring.a", "damper.a", "shaker.p"]},
    {"prescribe": "spring.b", "velocity_m_s": 0},
    {"prescribe": "damper.b", "velocity_m_s": 0},
    {"route": "damper.heat", "to": [["case.heat", 1.0]]}
  ]
}"#;
    let path = scratch("element_oscillator.json", elements);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let midpoint = value(&json, &element("block"), "position_m");
    let error = (midpoint - displacement).abs();
    assert!(
        error > 1.0e-6 * displacement_scale && error < 0.1 * displacement_scale,
        "{midpoint} vs {displacement}"
    );
}

#[test]
fn port_hamiltonian_states_are_adaptive_error_accounts_by_index() {
    let scene = std::fs::read_to_string(example("port_hamiltonian_oscillator.json")).unwrap();
    let adaptive = r#""step_s": 1.0e-3,
  "adaptive": {"order": 2, "tolerances": {
    "oscillator.x0": {"absolute": 1e-9, "relative": 0},
    "oscillator.x1": {"absolute": 1e-10, "relative": 0},
    "shaker.energy_in_j": {"absolute": 1e-10, "relative": 0},
    "case.temperature_k": {"absolute": 1e-9, "relative": 0}
  }},"#;
    let changed = scene.replace("\"step_s\": 1.0e-3,", adaptive);
    assert_ne!(changed, scene);
    let path = scratch("phs_adaptive.json", &changed);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let [momentum, displacement] = held_force_oscillator(0.5);
    let x0 = value(&json, &element("oscillator"), "x0");
    let x1 = value(&json, &element("oscillator"), "x1");
    assert!((x0 - momentum).abs() <= 1.0e-9 * 0.05, "{x0} vs {momentum}");
    assert!(
        (x1 - displacement).abs() <= 1.0e-9 * 0.005,
        "{x1} vs {displacement}"
    );

    // Every state coordinate needs its own tolerance, under its index.
    let missing = changed.replace(
        "    \"oscillator.x1\": {\"absolute\": 1e-10, \"relative\": 0},\n",
        "",
    );
    assert_ne!(missing, changed);
    let (code, stderr) = refused("phs_adaptive_missing.json", &missing);
    assert_eq!(code, 3, "{stderr}");
    assert!(
        stderr.contains("adaptive.tolerances.oscillator.x1: missing explicit tolerance"),
        "{stderr}"
    );
    let (code, stderr) = refused(
        "phs_adaptive_unknown.json",
        &changed.replace("oscillator.x1", "oscillator.x2"),
    );
    assert_eq!(code, 3, "{stderr}");
    assert!(
        stderr.contains("adaptive.tolerances.oscillator.x2: not a dynamic state coordinate"),
        "{stderr}"
    );
}

#[test]
fn a_flow_in_port_hamiltonian_spring_follows_the_spring_element_under_the_same_end_motion() {
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.25, "step_s": 0.05,
  "elements": {
    "matrix_spring": {"kind": "port_hamiltonian",
      "energy": [[400.0]], "structure": [[0.0]], "dissipation": [[0.0]],
      "input": [[1.0]], "state": [0.001],
      "ports": [{"name": "a", "domain": "translation", "role": "flow_in"}]},
    "spring": {"kind": "spring", "stiffness_n_m": 400.0, "extension_m": 0.001},
    "case": {"kind": "heat_capacity", "capacity_j_k": 10.0, "temperature_k": 300.0}
  },
  "connections": [
    {"prescribe": "matrix_spring.a", "velocity_m_s": 0.2},
    {"prescribe": "spring.a", "velocity_m_s": 0.2},
    {"prescribe": "spring.b", "velocity_m_s": 0},
    {"route": "matrix_spring.heat", "to": [["case.heat", 1.0]]}
  ]
}"#;
    let path = scratch("phs_spring.json", scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let extension = value(&json, &element("spring"), "extension_m");
    let coordinate = value(&json, &element("matrix_spring"), "x0");
    // 0.2 m/s held on one end for 0.25 s moves the extension by 0.05 m.
    assert!(
        ((extension - 0.001).abs() - 0.05).abs() <= 1.0e-12,
        "{extension}"
    );
    assert!(
        (coordinate - extension).abs() <= 1.0e-12,
        "{coordinate} vs {extension}"
    );
    let stored = 200.0 * (extension * extension - 0.001 * 0.001);
    for name in ["spring", "matrix_spring"] {
        let change = value(&json, &element(name), "stored_change_j");
        assert!(
            (change - stored).abs() <= 1.0e-12,
            "{name}: {change} vs {stored}"
        );
    }
    // The returned effort is the step mean of k x over the last step.
    let force = value(&json, &element("matrix_spring"), "a.force_n");
    let mean = 400.0 * (extension + 0.5 * (0.001 - extension) / 5.0);
    assert!((force - mean).abs() <= 1.0e-9, "{force} vs {mean}");
    assert!(value(&json, &element("matrix_spring"), "dissipated_j").abs() <= 1.0e-15);
}

/// Each case replaces `from` in port_hamiltonian_oscillator.json by `to`;
/// the scene must then refuse as an input, naming `reason`.
fn port_hamiltonian_refusals(cases: &[(&str, &str, &str, &str)]) {
    let scene = std::fs::read_to_string(example("port_hamiltonian_oscillator.json")).unwrap();
    for &(name, from, to, reason) in cases {
        assert!(scene.contains(from), "{name}");
        let (code, stderr) = refused(&format!("phs_{name}.json"), &scene.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
}

#[test]
fn port_hamiltonian_scenes_refuse_invalid_matrices_by_name() {
    port_hamiltonian_refusals(&[
        (
            "skew",
            "[[0.0, -1.0], [1.0, 0.0]]",
            "[[0.0, -1.0], [2.0, 0.0]]",
            "elements.oscillator: J violates its symmetry class",
        ),
        (
            "negative_dissipation",
            "[[1.5, 0.0], [0.0, 0.0]]",
            "[[-1.5, 0.0], [0.0, 0.0]]",
            "elements.oscillator: R is not positive semidefinite",
        ),
        (
            "asymmetric_dissipation",
            "[[1.5, 0.0], [0.0, 0.0]]",
            "[[1.5, 0.5], [0.0, 0.0]]",
            "elements.oscillator: R violates its symmetry class",
        ),
        (
            "indefinite_energy",
            "[[4.0, 0.0], [0.0, 400.0]]",
            "[[4.0, 0.0], [0.0, -400.0]]",
            "elements.oscillator.energy: Q is not positive semidefinite",
        ),
        (
            "asymmetric_energy",
            "[[4.0, 0.0], [0.0, 400.0]]",
            "[[4.0, 1.0], [0.0, 400.0]]",
            "elements.oscillator.energy: Q violates its symmetry class",
        ),
        (
            "structure_rows",
            "[[0.0, -1.0], [1.0, 0.0]]",
            "[[0.0, -1.0]]",
            "elements.oscillator.structure: matrix row count must be 2",
        ),
        (
            "input_columns",
            "[[1.0], [0.0]]",
            "[[1.0, 0.0], [0.0, 0.0]]",
            "elements.oscillator.input: matrix[0] column count must be 1",
        ),
        (
            "text",
            "[[1.5, 0.0], [0.0, 0.0]]",
            "[[1.5, 0.0], [0.0, \"none\"]]",
            "elements.oscillator.dissipation: matrix[1][1] must be a finite number",
        ),
        (
            "nonfinite",
            "[[1.5, 0.0], [0.0, 0.0]]",
            "[[1.5, 0.0], [0.0, 1e999]]",
            "not finite",
        ),
    ]);
    let scene = std::fs::read_to_string(example("port_hamiltonian_oscillator.json")).unwrap();
    let oversized = format!("[{}]", vec!["[1]"; 65].join(","));
    let (code, stderr) = refused(
        "phs_oversized.json",
        &scene.replace("[[4.0, 0.0], [0.0, 400.0]]", &oversized),
    );
    assert_eq!(code, 3, "{stderr}");
    assert!(
        stderr.contains("elements.oscillator.energy: Q must have 1 to 64 square rows"),
        "{stderr}"
    );
}

#[test]
fn port_hamiltonian_scenes_refuse_invalid_ports_state_and_fields_by_name() {
    let port = "[{\"name\": \"p\", \"domain\": \"translation\", \"role\": \"effort_in\"}]";
    port_hamiltonian_refusals(&[
        (
            "role",
            "\"effort_in\"",
            "\"effort\"",
            "elements.oscillator.ports[0]: role is",
        ),
        (
            "domain",
            "\"translation\"",
            "\"hydraulic\"",
            "elements.oscillator.ports[0]: domain is one of",
        ),
        (
            "reserved_heat",
            "{\"name\": \"p\"",
            "{\"name\": \"heat\"",
            "heat is reserved",
        ),
        (
            "state_length",
            "\"state\": [0.0, 0.0]",
            "\"state\": [0.0]",
            "elements.oscillator.state: 2 finite numbers",
        ),
        (
            "no_ports",
            port,
            "[]",
            "elements.oscillator.ports: name 1 to 64 power ports",
        ),
        (
            "unknown_field",
            "\"dissipation\":",
            "\"damping\":",
            "unsupported field",
        ),
    ]);
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
    // A charged capacitor rings down through an inductor and a resistor:
    // `examples/system/rlc_loop.json`, which holds these four values.
    let (henry, farad, ohm, coulomb) = (10.0e-3_f64, 100.0e-6_f64, 2.0_f64, 1.0e-3_f64);
    let (code, json, stderr) = run_binary(&["--json", "system", &example("rlc_loop.json")]);
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

    // A shaft on a torsion spring, turned by a constant torque:
    // `examples/system/torsion_pendulum.json`, which holds these three values.
    let (inertia, kappa, torque) = (2.0e-3_f64, 5.0_f64, 0.4_f64);
    let (code, json, stderr) = run_binary(&["--json", "system", &example("torsion_pendulum.json")]);
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
    // Every native body account needs an explicit accuracy tolerance.
    let adaptive = ROTOR.replace(
        "\"elements\"",
        "\"adaptive\": {\"tolerances\": {}},\n  \"elements\"",
    );
    let (code, stderr) = refused("adaptive.json", &adaptive);
    assert_eq!(code, 3);
    assert!(stderr.contains("rotor.body0.position_m: missing explicit tolerance"), "{stderr}");
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
    // What a gas volume holds are error accounts of adaptive stepping, and
    // each needs its tolerance.
    check(
        "vessel_adaptive.json",
        FILLED_VESSEL,
        "\"step_s\": 0.01,",
        "\"step_s\": 0.01, \"adaptive\": {\"order\": 2, \"min_step_s\": 1e-9, \"max_attempts\": 8, \"tolerances\": {}},",
        3,
        "adaptive.tolerances.vessel.mass_kg: missing explicit tolerance",
    );
}

/// The step entries of a scene under error control: a greatest step, and an
/// `adaptive` object of order 2 with an absolute tolerance for each account.
fn error_controlled(max_step_s: f64, tolerances: &[(&str, f64)]) -> String {
    let accounts: Vec<String> = tolerances
        .iter()
        .map(|(name, absolute)| format!(r#""{name}": {{"absolute": {absolute:e}, "relative": 0}}"#))
        .collect();
    format!(
        r#""step_s": {max_step_s:e}, "adaptive": {{"order": 2, "min_step_s": 1e-9, "max_attempts": 40, "tolerances": {{{}}}}},"#,
        accounts.join(", ")
    )
}

/// An error-controlled run that reaches its end: its report, its accepted
/// steps and the comparisons its driver rejected. The driver counts every
/// comparison of a full step with two half steps that it attempts, and an
/// accepted one leaves its two half steps as accepted steps.
fn controlled_run(name: &str, scene: &str) -> (String, f64, f64) {
    let path = scratch(name, scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{name}: {stderr}");
    assert!(
        value(&json, "\"adaptive\"", "max_error_norm") <= 1.0,
        "{name}"
    );
    let steps = value(&json, "\"schema\"", "steps");
    let attempts = value(&json, "\"adaptive\"", "attempts");
    (json, steps, attempts - 0.5 * steps)
}

/// The pneumatic cylinder over 0.5 s with one carried quantity, which the
/// tank holds alone and the valve takes to the chamber; `steps` are its step
/// entries.
fn carrying_cylinder(steps: &str) -> String {
    let mut scene = std::fs::read_to_string(example("pneumatic_cylinder.json")).unwrap();
    for (from, to) in [
        ("\"duration_s\": 3.0", "\"duration_s\": 0.5"),
        ("\"step_s\": 2.0e-4,", steps),
        (
            "\"temperature_k\": 293.15},",
            "\"temperature_k\": 293.15, \"carried\": [1.0]},",
        ),
        (
            "\"conductance_kg_s_pa\": 1.0e-8}",
            "\"conductance_kg_s_pa\": 1.0e-8, \"carried\": 1}",
        ),
        (
            "\"temperature_k\": 293.15, \"wall\": true}",
            "\"temperature_k\": 293.15, \"wall\": true, \"carried\": [0.0]}",
        ),
    ] {
        assert_eq!(scene.matches(from).count(), 1, "{from}");
        scene = scene.replace(from, to);
    }
    scene
}

/// Claims, declared before the first run: with `adaptive`, what a gas volume
/// holds are error accounts: its mass, its energy, its volume and each
/// carried mass. On the pneumatic cylinder over 0.5 s, with every other
/// tolerance too loose to bind, tightening one of those four alone on the
/// two volumes takes more accepted steps and rejects a comparison, and
/// tightening all four ends closer to a run at a fixed small step. A scene
/// without the chamber's carried mass is refused by that account's name.
#[test]
fn what_a_gas_volume_holds_are_adaptive_error_accounts() {
    // The size of what each account of a volume measures.
    let sizes = [
        [("tank.mass_kg", 1.4e-2), ("chamber.mass_kg", 4.8e-5)],
        [("tank.energy_j", 3.0e3), ("chamber.energy_j", 10.0)],
        [("tank.volume_m3", 2.0e-3), ("chamber.volume_m3", 4.0e-5)],
        [
            ("tank.carried0_kg", 1.4e-2),
            ("chamber.carried0_kg", 4.8e-5),
        ],
    ];
    // Tolerances of 1e-9 of the size for the accounts `tight` picks and of a
    // tenth of it for the other gas accounts; the rest never bind.
    let accounts = |tight: &dyn Fn(usize) -> bool| {
        let mut all = vec![
            ("piston.position_m", 1.0e3),
            ("piston.velocity_m_s", 1.0e3),
            ("outside.energy_in_j", 1.0e3),
            ("load.extension_m", 1.0e3),
            ("barrel.temperature_k", 1.0e3),
        ];
        for (index, pair) in sizes.iter().enumerate() {
            let part = if tight(index) { 1.0e-9 } else { 1.0e-1 };
            all.extend(pair.iter().map(|(name, size)| (*name, part * size)));
        }
        all
    };
    let run = |name: &str, tight: &dyn Fn(usize) -> bool| {
        controlled_run(
            name,
            &carrying_cylinder(&error_controlled(5.0e-3, &accounts(tight))),
        )
    };
    let (loose, loose_steps, _) = run("cylinder_loose.json", &|_| false);
    for (index, pair) in sizes.iter().enumerate() {
        let name = format!("cylinder_tight_{index}.json");
        let (_, steps, rejected) = run(&name, &|account| account == index);
        eprintln!(
            "cylinder: {} alone, {steps} steps and {rejected} rejected, loose {loose_steps}",
            pair[1].0
        );
        assert!(
            steps > loose_steps,
            "{}: {steps} vs {loose_steps}",
            pair[1].0
        );
        assert!(rejected >= 1.0, "{}: no comparison rejected", pair[1].0);
    }
    let (tight, _, _) = run("cylinder_tight.json", &|_| true);

    let ends = |json: &str| {
        [
            ("tank", "mass_kg"),
            ("tank", "pressure_pa"),
            ("chamber", "mass_kg"),
            ("chamber", "pressure_pa"),
        ]
        .map(|(name, field)| value(json, &element(name), field))
    };
    let path = scratch(
        "cylinder_fixed.json",
        &carrying_cylinder("\"step_s\": 2.0e-5,"),
    );
    let (code, fixed, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let reference = ends(&fixed);
    let distance = |json: &str| {
        let parts = ends(json).into_iter().zip(reference);
        parts
            .map(|(end, fixed)| (end / fixed - 1.0).abs())
            .fold(0.0, f64::max)
    };
    let (far, near) = (distance(&loose), distance(&tight));
    eprintln!("cylinder: from the fixed run, loose {far:e}, tight {near:e}");
    assert!(near < 0.1 * far, "tight {near}, loose {far}");

    let partial = accounts(&|_| true);
    let without: Vec<_> = partial
        .iter()
        .copied()
        .filter(|(name, _)| *name != "chamber.carried0_kg")
        .collect();
    let (code, stderr) = refused(
        "cylinder_partial.json",
        &carrying_cylinder(&error_controlled(5.0e-3, &without)),
    );
    assert_eq!(code, 3, "{stderr}");
    let missing = "adaptive.tolerances.chamber.carried0_kg: missing explicit tolerance";
    assert!(stderr.contains(missing), "{stderr}");

    // A gas volume steps by the trapezoidal rule: an order above two is refused.
    let third = error_controlled(5.0e-3, &partial).replace("\"order\": 2", "\"order\": 3");
    let (code, stderr) = refused("cylinder_order.json", &carrying_cylinder(&third));
    assert_eq!(code, 3, "{stderr}");
    let order = "adaptive.order: tank is a gas_volume, whose step is of order 2; \
                 the scene declares order 3";
    assert!(stderr.contains(order), "{stderr}");
}

/// Two modes of a bar, at 200 Hz and 510 Hz, driven at 150 Hz by a force
/// through one attachment over 0.02 s; `steps` are its step entries. The force
/// varies within a step, so the structure's step, which is exact for a held
/// force, is not exact here.
fn driven_bar(steps: &str) -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0",
  "duration_s": 0.02,
  {steps}
  "elements": {{
    "bar": {{"kind": "modal_structure", "sample_rate_hz": 48000,
      "modes": [{{"frequency_hz": 200.0, "damping_ratio": 0.01}},
                {{"frequency_hz": 510.0, "damping_ratio": 0.01}}],
      "attachments": {{"drive": {{"domain": "translation", "weights": [[1.0, 0.7]]}}}}}},
    "shaker": {{"kind": "force",
      "force_n": {{"sine": {{"amplitude": 1.0, "frequency_hz": 150.0}}}}}},
    "case": {{"kind": "heat_capacity", "capacity_j_k": 1.0, "temperature_k": 300.0}}
  }},
  "connections": [
    {{"join": ["shaker.p", "bar.drive"]}},
    {{"join": ["case.heat", "bar.heat"]}}
  ]
}}"#
    )
}

/// Claims, declared before the first run: with `adaptive`, a modal
/// structure's whole state is one error account, `energy_norm_sqrt_j`. On the
/// driven bar, with every other tolerance too loose to bind, tightening it
/// takes more accepted steps, rejects a comparison and ends closer to a run at
/// a fixed small step. A scene without it is refused by its name, and a scene
/// that declares an order above the structure's two is refused by the
/// structure's name.
#[test]
fn a_modal_structures_state_is_one_adaptive_error_account_in_its_energy_norm() {
    let accounts = |bar: f64| {
        [
            ("bar.energy_norm_sqrt_j", bar),
            ("shaker.energy_in_j", 1.0e3),
            ("case.temperature_k", 1.0e3),
        ]
    };
    let run = |name: &str, bar: f64| {
        controlled_run(
            name,
            &driven_bar(&error_controlled(4.0e-4, &accounts(bar))),
        )
    };
    let (loose, loose_steps, _) = run("bar_loose.json", 1.0e3);
    let (tight, steps, rejected) = run("bar_tight.json", 2.0e-9);
    eprintln!("bar: tight {steps} steps and {rejected} rejected, loose {loose_steps}");
    assert!(steps > loose_steps, "{steps} vs {loose_steps}");
    assert!(rejected >= 1.0, "no comparison rejected");

    let ends = |json: &str| {
        [
            "mode0.displacement",
            "mode0.velocity",
            "mode1.displacement",
            "mode1.velocity",
        ]
        .map(|field| value(json, &element("bar"), field))
    };
    let path = scratch("bar_fixed.json", &driven_bar("\"step_s\": 2.5e-6,"));
    let (code, fixed, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let reference = ends(&fixed);
    let distance = |json: &str| {
        let parts = ends(json).into_iter().zip(reference);
        parts
            .map(|(end, fixed)| (end / fixed - 1.0).abs())
            .fold(0.0, f64::max)
    };
    let (far, near) = (distance(&loose), distance(&tight));
    eprintln!("bar: from the fixed run, loose {far:e}, tight {near:e}");
    assert!(near < 0.1 * far, "tight {near}, loose {far}");

    let (code, stderr) = refused(
        "bar_partial.json",
        &driven_bar(&error_controlled(4.0e-4, &accounts(1.0)[1..])),
    );
    assert_eq!(code, 3, "{stderr}");
    let missing = "adaptive.tolerances.bar.energy_norm_sqrt_j: missing explicit tolerance";
    assert!(stderr.contains(missing), "{stderr}");

    // Exact for a held force, and a force that varies is held: order two.
    let third = error_controlled(4.0e-4, &accounts(1.0)).replace("\"order\": 2", "\"order\": 3");
    let (code, stderr) = refused("bar_order.json", &driven_bar(&third));
    assert_eq!(code, 3, "{stderr}");
    let order = "adaptive.order: bar is a modal_structure, whose step is of order 2; \
                 the scene declares order 3";
    assert!(stderr.contains(order), "{stderr}");
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

/// A 2 kg block thrown at `velocity` along a still surface, through a
/// friction interface of `law` under 20 N whose heat goes to a 50 J/K
/// lining at 300 K.
fn thrown_block(law: &str, velocity: f64, step_s: f64, duration_s: f64) -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": {duration_s}, "step_s": {step_s},
  "elements": {{
    "block": {{"kind": "mass", "mass_kg": 2.0, "velocity_m_s": {velocity}}},
    "friction": {{"kind": "friction_interface", "normal_n": 20.0, "law": {law},
                 "slip_m_s": {velocity}}},
    "lining": {{"kind": "heat_capacity", "capacity_j_k": 50.0, "temperature_k": 300.0}}
  }},
  "connections": [
    {{"join": ["block.p", "friction.a"]}},
    {{"prescribe": "friction.b", "velocity_m_s": 0.0}},
    {{"join": ["lining.heat", "friction.heat"]}}
  ]
}}"#
    )
}

#[test]
fn a_coulomb_friction_interface_arrests_a_thrown_block_in_either_direction() {
    // G1: fs-tribo's Coulomb rung, kinetic coefficient 0.2 under 20 N,
    // decelerates the 2 kg block at 2 m/s2: from 3 m/s it stops at 1.5 s
    // after 2.25 m, and its 9 J of kinetic energy are the lining's heat.
    for sign in [-1.0, 1.0] {
        let law = r#"{"coulomb": {"static_mu": 0.5, "kinetic_mu": 0.2}}"#;
        let path = scratch(
            if sign > 0.0 {
                "forward.json"
            } else {
                "backward.json"
            },
            &thrown_block(law, sign * 3.0, 0.005, 2.0),
        );
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        assert!(stderr.is_empty(), "{stderr}");
        for (owner, field, expected, tolerance) in [
            ("block", "position_m", sign * 2.25, 1.0e-8),
            ("block", "velocity_m_s", 0.0, 1.0e-10),
            ("friction", "slip_m_s", 0.0, 1.0e-10),
            ("friction", "stuck", 1.0, 0.0),
            ("friction", "dissipated_j", 9.0, 1.0e-9),
            ("lining", "temperature_k", 300.18, 1.0e-9),
        ] {
            let found = value(&json, &element(owner), field);
            assert!(
                (found - expected).abs() <= tolerance,
                "{sign} {owner}.{field}: {found} vs {expected}"
            );
        }
        assert!(
            value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
            "{json}"
        );
        assert!(
            json.contains("friction interfaces step an fs-tribo dry law"),
            "{json}"
        );
    }
}

#[test]
fn a_belt_at_constant_speed_drags_a_block_until_its_friction_interface_catches_it() {
    // G1: the block starts at rest on a belt already moving at 1 m/s, a
    // slip of -1 m/s in the movers' sign. Kinetic friction, 0.2 x 20 N,
    // accelerates the 2 kg block at 2 m/s2: it is caught at 0.5 s after
    // 0.25 m while the belt moved 0.5 m, and static friction (10 N
    // available, none needed) carries it on. At 1 s: 0.75 m and 1 m/s,
    // 4 N x 0.25 m = 1 J of heat, 1 J of kinetic energy, 2 J from the belt.
    // The belt is the interface's prescribed end, and then a velocity
    // source joined to the end, in 200 whole steps each. (Bead `.5.8.1`:
    // the kernel differenced the source's flow, and every step of the stuck
    // interface was cut to the stop of the kernel's own perturbation.)
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 1.0, "step_s": 0.005,
  "elements": {
    "block": {"kind": "mass", "mass_kg": 2.0},
    "friction": {"kind": "friction_interface", "normal_n": 20.0,
                 "law": {"coulomb": {"static_mu": 0.5, "kinetic_mu": 0.2}}, "slip_m_s": -1.0},
    "lining": {"kind": "heat_capacity", "capacity_j_k": 50.0, "temperature_k": 300.0}
  },
  "connections": [
    {"join": ["block.p", "friction.a"]},
    {"prescribe": "friction.b", "velocity_m_s": 1.0},
    {"join": ["lining.heat", "friction.heat"]}
  ]
}"#;
    let prescribed = r#"{"prescribe": "friction.b", "velocity_m_s": 1.0}"#;
    let mass = r#""block": {"kind": "mass", "mass_kg": 2.0},"#;
    assert!(scene.contains(prescribed) && scene.contains(mass));
    let driven = scene
        .replace(prescribed, r#"{"join": ["belt.p", "friction.b"]}"#)
        .replace(
            mass,
            &format!(r#"{mass} "belt": {{"kind": "velocity_source", "velocity_m_s": 1.0}},"#),
        );
    for (name, scene) in [
        ("belt.json", scene.to_owned()),
        ("driven_belt.json", driven),
    ] {
        let path = scratch(name, &scene);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        for (owner, field, expected, tolerance) in [
            ("block", "position_m", 0.75, 1.0e-8),
            ("block", "velocity_m_s", 1.0, 1.0e-9),
            ("friction", "slip_m_s", 0.0, 1.0e-9),
            ("friction", "stuck", 1.0, 0.0),
            ("friction", "dissipated_j", 1.0, 1.0e-8),
            ("lining", "temperature_k", 300.02, 1.0e-9),
        ] {
            let found = value(&json, &element(owner), field);
            assert!(
                (found - expected).abs() <= tolerance,
                "{name}: {owner}.{field}: {found} vs {expected}"
            );
        }
        let steps = value(&json, "\"schema\"", "steps");
        let shortenings = value(&json, "\"schema\"", "shortenings");
        assert!(
            (steps - 200.0).abs() <= 0.0 && shortenings.abs() <= 0.0,
            "{name}: {steps} steps, {shortenings} shortenings"
        );
        // What crossed the prescribed end, or what the source gave.
        let belt = if name == "belt.json" {
            value(&json, "\"port\":\"friction.b\"", "energy_in_j")
        } else {
            -value(&json, &element("belt"), "energy_in_j")
        };
        assert!((belt - 2.0).abs() <= 4.0e-8, "{name}: belt {belt}");
        assert!(
            value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
            "{name}: {json}"
        );
    }
}

/// `(m / N) * integral from 0 to v0 of v / mu(v) dv`: where a block of mass
/// `m` thrown at `v0` stops under `mu(v) N`, by composite Simpson's rule.
fn stopping_distance(mu: impl Fn(f64) -> f64, mass_kg: f64, normal_n: f64, v0: f64) -> f64 {
    let intervals = 20_000;
    let h = v0 / f64::from(intervals);
    let integrand = |v: f64| v / mu(v);
    let mut sum = integrand(0.0) + integrand(v0);
    for i in 1..intervals {
        let weight = if i % 2 == 1 { 4.0 } else { 2.0 };
        sum += weight * integrand(f64::from(i) * h);
    }
    mass_kg / normal_n * sum * h / 3.0
}

#[test]
fn velocity_weakening_friction_interfaces_stop_a_thrown_block_where_the_law_s_integral_says() {
    // G1/G3: a block of mass m thrown at v0 under mu(v) N stops after
    // (m / N) * integral of v / mu(v) dv from 0 to v0, all of its kinetic
    // energy heat. The interface asks its law at the step's mean slip, so
    // the kernel's stop converges at second order: halving the step
    // divides the error by about four. The Stribeck block is the example
    // scene; under 20 N both laws decelerate the block at 2 m/s2 or more,
    // so a duration of v0 seconds is well past the stop.
    let example_scene = std::fs::read_to_string(example("stribeck_block.json")).unwrap();
    let stribeck = |v: f64| 0.2 + 0.3 * (-(v / 0.5) * (v / 0.5)).exp() + 0.05 * v;
    stops_at_second_order("stribeck", &example_scene, stribeck, 3.0, ["0.01", "0.005"]);
    let rosin_scene = thrown_block(r#""rosin_woodhouse_galluzzo_2025_eq13""#, 0.3, 1.0e-3, 0.3);
    let rosin = |v: f64| 0.2 * (v + 0.011).powf(-0.4);
    stops_at_second_order("rosin", &rosin_scene, rosin, 0.3, ["0.001", "0.0005"]);
}

/// Runs a thrown-block `scene` written at the `coarse` step and at the
/// `fine` one, and asserts the block stops, all of its kinetic energy heat,
/// within 1e-3 of the law's distance, the error divided by 3 to 5.
fn stops_at_second_order(
    name: &str,
    scene: &str,
    mu: impl Fn(f64) -> f64,
    v0: f64,
    [coarse, fine]: [&str; 2],
) {
    let reference = stopping_distance(mu, 2.0, 20.0, v0);
    let coarse_step = format!("\"step_s\": {coarse}");
    assert!(scene.contains(&coarse_step), "{name}: test fixture changed");
    let mut errors = Vec::new();
    for (index, step_s) in [coarse, fine].into_iter().enumerate() {
        let path = scratch(
            &format!("{name}{index}.json"),
            &scene.replace(&coarse_step, &format!("\"step_s\": {step_s}")),
        );
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name} {step_s}: {stderr}");
        for (owner, field, expected, tolerance) in [
            ("friction", "stuck", 1.0, 0.0),
            ("block", "velocity_m_s", 0.0, 1.0e-9),
            ("friction", "dissipated_j", v0 * v0, 1.0e-9 * v0 * v0),
        ] {
            let found = value(&json, &element(owner), field);
            assert!(
                (found - expected).abs() <= tolerance,
                "{name} {step_s}: {owner}.{field} {found} vs {expected}"
            );
        }
        assert!(
            value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
            "{json}"
        );
        let position = value(&json, &element("block"), "position_m");
        eprintln!("{name}: step {step_s} s, stops at {position} m, reference {reference} m");
        errors.push((position - reference).abs());
    }
    let ratio = errors[0] / errors[1];
    assert!(
        (3.0..5.0).contains(&ratio) && errors[1] < 1.0e-3 * reference,
        "{name}: errors {errors:?}, ratio {ratio}"
    );
}

#[test]
fn friction_interface_scenes_refuse_invalid_laws_budgets_spans_and_error_control() {
    let scene = thrown_block(
        r#"{"coulomb": {"static_mu": 0.5, "kinetic_mu": 0.2}}"#,
        3.0,
        0.005,
        2.0,
    );
    let coulomb = r#""law": {"coulomb": {"static_mu": 0.5, "kinetic_mu": 0.2}}"#;
    let no_iterations = format!(r#"{coulomb}, "solve": {{"max_iterations": 0}}"#);
    let no_impulse = format!(r#"{coulomb}, "solve": {{"impulse_tolerance_n_s": 0.0}}"#);
    // Every one is an input refusal, exit 3.
    for (name, from, to, reason) in [
        (
            "rung.json",
            coulomb,
            r#""law": {"lubricated": {"static_mu": 0.5}}"#,
            "unknown law \"lubricated\"",
        ),
        (
            "named.json",
            coulomb,
            r#""law": "rosin""#,
            "the named law is",
        ),
        (
            "two.json",
            coulomb,
            r#""law": {"coulomb": {"static_mu": 0.5, "kinetic_mu": 0.2}, "hyperbolic": {}}"#,
            "an object with one key",
        ),
        (
            "kinetic.json",
            coulomb,
            r#""law": {"coulomb": {"static_mu": 0.2, "kinetic_mu": 0.5}}"#,
            "kinetic_mu > static_mu",
        ),
        (
            "coefficient.json",
            coulomb,
            r#""law": {"coulomb": {"static_mu": 0.5, "kinetic_mu": 0.2, "rolling_mu": 0.1}}"#,
            "unsupported field \"rolling_mu\"",
        ),
        (
            "speed.json",
            coulomb,
            r#""law": {"hyperbolic": {"static_mu": 0.5, "kinetic_mu": 0.2, "characteristic_speed_m_s": 0.0}}"#,
            "characteristic_speed",
        ),
        (
            "missing.json",
            coulomb,
            r#""damping": 1.0"#,
            "unsupported field \"damping\"",
        ),
        (
            "iterations.json",
            coulomb,
            no_iterations.as_str(),
            "max_iterations is a whole number from 1 to 128",
        ),
        (
            "budget.json",
            coulomb,
            no_impulse.as_str(),
            "impulse_tolerance_n_s",
        ),
        (
            "normal.json",
            "\"normal_n\": 20.0",
            "\"normal_n\": -20.0",
            "normal_force",
        ),
        (
            "span.json",
            "{\"prescribe\": \"friction.b\", \"velocity_m_s\": 0.0}",
            "{\"prescribe\": \"friction.b\", \"velocity_m_s\": 0.0},\n    \
             {\"prescribe\": \"friction.span_b\", \"velocity_m_s\": 0.0}",
            "a friction_interface has ports a, b, heat",
        ),
        (
            "adaptive.json",
            "\"duration_s\": 2",
            "\"adaptive\": {\"tolerances\": {}}, \"duration_s\": 2",
            "friction is a friction_interface",
        ),
    ] {
        assert!(scene.contains(from), "{name}: test fixture changed");
        let (code, stderr) = refused(name, &scene.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
    let missing = scene.replace(&format!("{coulomb},"), "");
    let (code, stderr) = refused("missing_law.json", &missing);
    assert_eq!(code, 3, "{stderr}");
    assert!(stderr.contains("missing law"), "{stderr}");
}

/// One 5 Hz mode with 2 % damping, seen at its tip through the weight
/// `1 / sqrt(2 kg)`, so that the tip moves as a 2 kg mass on a spring and a
/// damper. A constant 1 N pushes the tip from rest for 1 s; the modal heat
/// warms a 10 J/K case.
const MODAL_TIP: &str = r#"{
  "schema": "frankensim.system.v0", "duration_s": 1.0, "step_s": 0.001,
  "elements": {
    "rod": {"kind": "modal_structure", "sample_rate_hz": 1000,
            "modes": [{"frequency_hz": 5.0, "damping_ratio": 0.02}],
            "attachments": {"tip": {"domain": "translation",
                                    "weights": [[0.7071067811865476]]}}},
    "push": {"kind": "force", "force_n": 1.0},
    "case": {"kind": "heat_capacity", "capacity_j_k": 10.0, "temperature_k": 300.0}
  },
  "connections": [
    {"join": ["rod.tip", "push.p"]},
    {"join": ["case.heat", "rod.heat"]}
  ]
}"#;

#[test]
fn a_modal_structure_pushed_by_a_constant_force_follows_the_exact_damped_step_response() {
    // G1: the owner steps each mode exactly under the effort held over the
    // step, and a constant force is held exactly. The mode q obeys
    // q'' + 2 zeta w q' + w^2 q = s F with q(0) = q'(0) = 0, so
    // q = s F / w^2 (1 - e^(-zeta w t) (cos wd t + zeta w / wd sin wd t)) and
    // q' = s F / wd e^(-zeta w t) sin wd t. The force's work is F s q, of
    // which the modal energy (q'^2 + w^2 q^2) / 2 is stored and the rest is
    // the case's heat.
    let path = scratch("modal_tip.json", MODAL_TIP);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (s, force, t) = (0.5_f64.sqrt(), 1.0, 1.0);
    let (w, zeta) = (core::f64::consts::TAU * 5.0, 0.02_f64);
    let wd = w * (1.0 - zeta * zeta).sqrt();
    let decay = (-zeta * w * t).exp();
    let q = s * force / (w * w) * (1.0 - decay * ((wd * t).cos() + zeta * w / wd * (wd * t).sin()));
    let v = s * force / wd * decay * (wd * t).sin();
    let work = force * s * q;
    let heat = work - 0.5 * (v * v + w * w * q * q);
    for (owner, field, expected, tolerance) in [
        ("rod", "mode0.displacement", q, 1.0e-12),
        ("rod", "mode0.velocity", v, 1.0e-11),
        ("rod", "tip.velocity_m_s", s * v, 1.0e-11),
        ("push", "energy_in_j", -work, 1.0e-12),
        ("rod", "dissipated_j", heat, 1.0e-12),
        // 1000 additions to a 300 K accumulator, each rounded to its ulp.
        ("case", "temperature_k", 300.0 + heat / 10.0, 1.0e-10),
    ] {
        let found = value(&json, &element(owner), field);
        assert!(
            (found - expected).abs() <= tolerance,
            "{owner}.{field}: {found} vs {expected}"
        );
    }
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
        "{json}"
    );
    assert!(
        json.contains("modal structures step fs-nlmodal's damped modes"),
        "{json}"
    );
}

#[test]
fn three_struck_modes_ring_down_as_their_closed_forms_and_heat_the_case() {
    // G1: three modes of a bar set ringing by a strike (initial modal
    // velocities), with nothing acting on the bar after it: each rings down
    // as q = v0 / wd e^(-zeta w t) sin wd t, and the energy the modes lose
    // is the case's heat. The strike point's flow is the weighted sum of
    // the modes' velocities.
    // (frequency [Hz], initial modal velocity [m sqrt(kg)/s], damping ratio)
    let modes = [
        (331.125_927_8, 0.4, 1.0e-3_f64),
        (1_324.505_956, -0.2, 2.0e-3),
        (2_972.946_707, 0.1, 5.0e-3),
    ];
    let weights = [2.1, -1.4, 0.9];
    let rows: Vec<String> = modes
        .iter()
        .map(|(hz, v0, zeta)| {
            format!(r#"{{"frequency_hz": {hz}, "damping_ratio": {zeta}, "velocity": {v0}}}"#)
        })
        .collect();
    let scene = format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 0.05, "step_s": 2.5e-5,
  "elements": {{
    "bar": {{"kind": "modal_structure", "sample_rate_hz": 40000,
            "modes": [{}],
            "attachments": {{"strike": {{"domain": "translation",
                                        "weights": [[2.1, -1.4, 0.9]]}}}}}},
    "case": {{"kind": "heat_capacity", "capacity_j_k": 0.01, "temperature_k": 300.0}}
  }},
  "connections": [
    {{"prescribe": "bar.strike", "force_n": 0}},
    {{"join": ["case.heat", "bar.heat"]}}
  ]
}}"#,
        rows.join(", ")
    );
    let path = scratch("struck_bar.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let t = 0.05;
    let (mut start_j, mut end_j, mut strike) = (0.0, 0.0, 0.0);
    for (index, &(hz, v0, zeta)) in modes.iter().enumerate() {
        let w = core::f64::consts::TAU * hz;
        let wd = w * (1.0 - zeta * zeta).sqrt();
        let decay = (-zeta * w * t).exp();
        let q = v0 / wd * decay * (wd * t).sin();
        let v = v0 * decay * ((wd * t).cos() - zeta * w / wd * (wd * t).sin());
        for (field, expected) in [("displacement", q), ("velocity", v)] {
            let found = value(&json, &element("bar"), &format!("mode{index}.{field}"));
            assert!(
                (found - expected).abs() <= 1.0e-9 * v0.abs(),
                "mode{index}.{field}: {found} vs {expected}"
            );
        }
        start_j += 0.5 * v0 * v0;
        end_j += 0.5 * (v * v + w * w * q * q);
        strike += weights[index] * v;
    }
    let heat = start_j - end_j;
    let found = value(&json, &element("bar"), "strike.velocity_m_s");
    assert!(
        (found - strike).abs() <= 1.0e-9,
        "strike: {found} vs {strike}"
    );
    let found = value(&json, &element("case"), "temperature_k");
    assert!(
        (found - 300.0 - heat / 0.01).abs() <= 1.0e-9 * heat / 0.01,
        "case: {found} vs heat {heat}"
    );
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
        "{json}"
    );
}

#[test]
fn a_bowed_mode_slides_under_a_fast_bow_and_catches_a_slow_one_as_the_closed_forms_say() {
    // G1: one undamped 2 Hz mode of unit modal mass (weight 1), so that its
    // contact moves as 1 kg on a spring of k = w^2, rubbed by a bow at speed
    // V through Coulomb friction under 10 N (static 5 N, kinetic 3 N), the
    // interface reading the end of the string's flow. From rest, kinetic
    // friction pushes the mode: x = Fk / k (1 - cos w t), v = Fk / w sin w t.
    // A bow faster than Fk / w (0.239 m/s) is never caught: at 0.5 m/s the
    // mode slides for the whole 0.9 s and the heat is Fk (V t - x). A bow at
    // 0.1 m/s is caught when v = V, at sin w t1 = V w / Fk (the step the
    // catch falls in ends at it); the contact then rides with the bow while
    // the spring needs less than the static 5 N (until x = 5 N / k, at
    // 0.33 s), with the heat of the slide before the catch. Stuck, the
    // interface's effort is held over each step and the contact's flow is
    // the bow's at both ends of it, so over a step of length s the mode
    // moves 2 V tan(w s / 2) / w, not V s, under F = k x + V w tan(w s / 2):
    // a held force where the riding one grows along the step. The stuck
    // steps are whole from the catch, the last cut at the end; at 0.2 s
    // the mode is 2.2e-7 m past x(t1) + V (0.2 - t1), a gap that falls as
    // the step squared. The held force's work on that gap,
    // -sum F (moved - V s), is the interface's numerical loss. In both
    // runs the bow, the interface's prescribed end, gives the heat, the
    // mode's energy (v^2 + k x^2) / 2 and that loss.
    let scene = std::fs::read_to_string(example("bowed_mode.json")).unwrap();
    let w = core::f64::consts::TAU * 2.0;
    let (k, fk, dt) = (w * w, 3.0, 0.001);
    for (name, speed, t) in [("slow", 0.1, 0.2), ("fast", 0.5, 0.9)] {
        let mut text = scene.clone();
        for (from, to) in [
            (
                "\"velocity_m_s\": 0.1",
                format!("\"velocity_m_s\": {speed}"),
            ),
            ("\"slip_m_s\": -0.1", format!("\"slip_m_s\": {}", -speed)),
            ("\"duration_s\": 0.2", format!("\"duration_s\": {t}")),
        ] {
            assert!(text.contains(from), "{name}: test fixture changed");
            text = text.replace(from, &to);
        }
        let path = scratch(&format!("bowed_{name}.json"), &text);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let (x, v, heat, stuck, loss) = if speed < fk / w {
            let t1 = (speed * w / fk).asin() / w;
            let x1 = fk / k * (1.0 - (w * t1).cos());
            let (mut x, mut at, mut loss) = (x1, t1, 0.0);
            while at < t {
                let s = (t - at).min(dt);
                let tan = (w * s / 2.0).tan();
                let moved = 2.0 * speed * tan / w;
                loss -= (k * x + speed * w * tan) * (moved - speed * s);
                x += moved;
                at += s;
            }
            (x, speed, fk * (speed * t1 - x1), 1.0, loss)
        } else {
            let x = fk / k * (1.0 - (w * t).cos());
            (x, fk / w * (w * t).sin(), fk * (speed * t - x), 0.0, 0.0)
        };
        let given = heat + 0.5 * (v * v + k * x * x) + loss;
        for (owner, field, expected, tolerance) in [
            ("string", "mode0.displacement", x, 1.0e-9),
            ("string", "contact.velocity_m_s", v, 1.0e-9),
            ("string", "dissipated_j", 0.0, 0.0),
            ("rubbing", "stuck", stuck, 0.0),
            ("rubbing", "dissipated_j", heat, 1.0e-9),
            ("rubbing", "numerical_loss_j", loss, 1.0e-12),
        ] {
            let found = value(&json, &element(owner), field);
            assert!(
                (found - expected).abs() <= tolerance,
                "{name} {owner}.{field}: {found} vs {expected}"
            );
        }
        let bow = value(&json, "\"port\":\"rubbing.b\"", "energy_in_j");
        assert!(
            (bow - given).abs() <= 1.0e-9,
            "{name} bow: {bow} vs {given}"
        );
        assert!(
            value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-9,
            "{json}"
        );
    }
}

#[test]
fn modal_structure_scenes_refuse_bad_clocks_modes_weights_and_error_control() {
    for (name, from, to, reason) in [
        (
            "rate.json",
            "\"sample_rate_hz\": 1000",
            "\"sample_rate_hz\": 999.5",
            "sample_rate_hz is a whole number of hertz",
        ),
        (
            "empty.json",
            "[{\"frequency_hz\": 5.0, \"damping_ratio\": 0.02}]",
            "[]",
            "modes is a nonempty array",
        ),
        (
            "field.json",
            "\"damping_ratio\": 0.02}",
            "\"damping_ratio\": 0.02, \"mass_kg\": 2.0}",
            "unsupported field \"mass_kg\"",
        ),
        (
            "weights.json",
            "[[0.7071067811865476]]",
            "[[0.7071067811865476, 0.5]]",
            "one weight per mode, 1 here",
        ),
        (
            "domain.json",
            "\"domain\": \"translation\"",
            "\"domain\": \"thermal\"",
            "domain is one of translation, rotation, electrical",
        ),
        (
            "nyquist.json",
            "\"frequency_hz\": 5.0",
            "\"frequency_hz\": 460.0",
            "FS-COUPLE-MODAL-TIME-NYQUIST",
        ),
        ("heat.json", "{\"tip\": {", "{\"heat\": {", "port refused"),
        (
            "adaptive.json",
            "\"duration_s\": 1.0",
            "\"adaptive\": {\"tolerances\": {\
             \"push.energy_in_j\": {\"absolute\": 1, \"relative\": 0}, \
             \"case.temperature_k\": {\"absolute\": 1, \"relative\": 0}}}, \"duration_s\": 1.0",
            "adaptive.tolerances.rod.energy_norm_sqrt_j: missing explicit tolerance",
        ),
    ] {
        assert!(MODAL_TIP.contains(from), "{name}: test fixture changed");
        let (code, stderr) = refused(name, &MODAL_TIP.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
}

/// The seed N2 packs, compiled by `xtask matdb-pack` from data/matdb/seed-v1:
/// species, McBride, Zehe and Gordon (2002) NASA-9, Lemmon and Jacobsen
/// (2004) dilute transport.
fn n2_sources() -> String {
    let pack = |name: &str| example(&format!("n2/{name}"));
    format!(
        r#"[{{"species": "{}", "nasa9": "{}", "transport": "{}"}}]"#,
        pack("species.pack"),
        pack("nasa9.pack"),
        pack("transport.fsmcdpk")
    )
}

fn seed(path: &str) -> String {
    std::fs::read_to_string(format!(
        "{}/../../data/matdb/seed-v1/{path}",
        env!("CARGO_MANIFEST_DIR")
    ))
    .unwrap()
}

/// N2's 200 to 1000 K coefficients a0..a8, read from the seed source the
/// packs were compiled from, not from the packs.
fn n2_nasa9_low() -> [f64; 9] {
    let mut a = [f64::NAN; 9];
    for line in seed("nitrogen-nasa9-mcbride2002/nasa9.tsv").lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() == 6 && fields[0] == "coefficient" && fields[2] == "low" {
            a[fields[3][1..].parse::<usize>().unwrap()] = fields[4].parse().unwrap();
        }
    }
    assert!(a.iter().all(|c| c.is_finite()), "{a:?}");
    a
}

/// N2's specific gas constant [J/(kg K)], from the seed species source.
fn n2_gas_constant() -> f64 {
    let species = seed("nitrogen/species.tsv");
    let row: Vec<&str> = species
        .lines()
        .find(|line| line.starts_with("species\t"))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(row[3], "g/mol");
    GAS_CONSTANT_J_MOL_K / (row[2].parse::<f64>().unwrap() * 1.0e-3)
}

/// The NASA-9 molar internal energy over R, u/R = h/R - T [K], less the
/// formation constant a7 (it cancels from every difference).
fn n2_energy(a: &[f64; 9], t: f64) -> f64 {
    -a[0] / t
        + a[1] * t.ln()
        + a[2] * t
        + a[3] * t.powi(2) / 2.0
        + a[4] * t.powi(3) / 3.0
        + a[5] * t.powi(4) / 4.0
        + a[6] * t.powi(5) / 5.0
        - t
}

/// The NASA-9 standard molar entropy over R, less its constant a8.
fn n2_entropy(a: &[f64; 9], t: f64) -> f64 {
    -a[0] / (2.0 * t * t) - a[1] / t
        + a[2] * t.ln()
        + a[3] * t
        + a[4] * t.powi(2) / 2.0
        + a[5] * t.powi(3) / 3.0
        + a[6] * t.powi(4) / 4.0
}

/// The NASA-9 cp/R.
fn n2_heat_capacity(a: &[f64; 9], t: f64) -> f64 {
    a[0] / (t * t)
        + a[1] / t
        + a[2]
        + a[3] * t
        + a[4] * t.powi(2)
        + a[5] * t.powi(3)
        + a[6] * t.powi(4)
}

/// The root of an increasing `f` between `lo` and `hi`, to the last bit.
fn bisect(mut lo: f64, mut hi: f64, f: impl Fn(f64) -> f64) -> f64 {
    assert!(f(lo) < 0.0 && f(hi) > 0.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if mid <= lo || mid >= hi {
            break;
        }
        if f(mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// 1 L of N2 at 250 K and 80 kPa, held at its volume and heated at 40 W
/// for 1 s.
fn heated_n2() -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 1.0, "step_s": 0.01,
  "elements": {{
    "cavity": {{"kind": "gas_cavity", "gas": "thermochemical 250 80000 mole 1 N2 1",
               "sources": {}, "volume_m3": 1.0e-3,
               "temperature_min_k": 200, "temperature_max_k": 900}}
  }},
  "connections": [
    {{"prescribe": "cavity.p", "flow_m3_s": 0}},
    {{"prescribe": "cavity.heat", "heat_rate_w": 40}}
  ]
}}"#,
        n2_sources()
    )
}

#[test]
fn a_gas_cavity_heated_at_fixed_volume_ends_where_the_nasa9_energy_integral_puts_it() {
    // G1: at fixed volume the heat is all internal energy, so the cavity
    // ends at the temperature where the McBride (2002) enthalpy integral
    // gives u(T1) - u(T0) = Q / m, at the pressure m R T1 / V. The oracle
    // reads the seed coefficients and molar mass the packs were compiled
    // from and integrates nothing; the owner's caloric solve stops within
    // 1e-7 J/kg a step, about 1e-8 K over the 100 steps.
    let path = scratch("heated.json", &heated_n2());
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (a, r) = (n2_nasa9_low(), n2_gas_constant());
    let (t0, p0, volume, heat) = (250.0, 80_000.0, 1.0e-3, 40.0);
    let mass = p0 * volume / (r * t0);
    let target = n2_energy(&a, t0) + heat / (mass * r);
    let t1 = bisect(t0, 400.0, |t| n2_energy(&a, t) - target);
    for (field, expected, tolerance) in [
        ("temperature_k", t1, 1.0e-7),
        ("pressure_pa", mass * r * t1 / volume, 1.0e-4),
        ("mass_kg", mass, 1.0e-15),
        ("volume_m3", volume, 0.0),
    ] {
        let found = value(&json, &element("cavity"), field);
        assert!(
            (found - expected).abs() <= tolerance,
            "cavity.{field}: {found} vs {expected}"
        );
    }
    let boundary = "\"port\":\"cavity.heat\"";
    assert!(
        (value(&json, boundary, "energy_in_j") - heat).abs() < 1.0e-9,
        "{json}"
    );
    assert!(
        (value(&json, "\"totals\"", "stored_change_j") - heat).abs() < 1.0e-6,
        "{json}"
    );
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-6,
        "{json}"
    );
    assert!(
        json.contains("a gas_cavity holds a closed, uniform, frozen ideal-gas mixture"),
        "{json}"
    );
    // The narrowed base sentence: a scene with a gas cavity resolves packs.
    assert!(
        json.contains("no material-card resolution besides gas cavities' packs, "),
        "{json}"
    );
}

#[test]
fn an_insulated_gas_cavity_expanded_by_a_prescribed_flow_keeps_to_the_nasa9_isentrope() {
    // G1: the owner's discrete-gradient step changes the entropy only by
    // the heat over its mean temperature, so an insulated cavity keeps its
    // entropy at every step, however long: it ends on the McBride (2002)
    // isentrope s0(T1)/R - s0(T0)/R - ln(T1/T0) + ln(V1/V0) = 0, and the
    // work its port took out is its internal energy lost. A negative
    // compression flow enlarges the cavity; 100 kPa falls to 82 kPa, inside
    // the transport card's window.
    let scene = format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 1.0, "step_s": 0.01,
  "elements": {{
    "cavity": {{"kind": "gas_cavity", "gas": "thermochemical 300 100000 mole 1 N2 1",
               "sources": {}, "volume_m3": 1.0e-3, "insulated": true,
               "temperature_min_k": 200, "temperature_max_k": 900}}
  }},
  "connections": [{{"prescribe": "cavity.p", "flow_m3_s": -1.5e-4}}]
}}"#,
        n2_sources()
    );
    let path = scratch("expanded.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (a, r) = (n2_nasa9_low(), n2_gas_constant());
    let (t0, p0, v0, v1) = (300.0, 100_000.0, 1.0e-3, 1.15e-3);
    let mass = p0 * v0 / (r * t0);
    let start = n2_entropy(&a, t0) - (v1 / v0).ln();
    let t1 = bisect(250.0, t0, |t| n2_entropy(&a, t) - (t / t0).ln() - start);
    let work = mass * r * (n2_energy(&a, t1) - n2_energy(&a, t0));
    for (field, expected, tolerance) in [
        ("temperature_k", t1, 1.0e-7),
        ("pressure_pa", mass * r * t1 / v1, 1.0e-4),
        ("volume_m3", v1, 1.0e-15),
    ] {
        let found = value(&json, &element("cavity"), field);
        assert!(
            (found - expected).abs() <= tolerance,
            "cavity.{field}: {found} vs {expected}"
        );
    }
    let boundary = "\"port\":\"cavity.p\"";
    assert!(
        (value(&json, boundary, "energy_in_j") - work).abs() < 1.0e-6,
        "{json}"
    );
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-6,
        "{json}"
    );
}

#[test]
fn a_piston_on_an_n2_gas_spring_rings_at_the_adiabatic_frequency_and_keeps_its_energy() {
    // examples/system/n2_gas_spring.json: a 1 kg piston of 10 cm2 on 0.1 L
    // of insulated N2 at 90 kPa and 300 K, held there by 90 N pressing the
    // way that shrinks it, set moving at 5 cm/s. G3: no energy is lost, so the ledger closes and the
    // piston's energy plus the gas's equals the outside force's work.
    // G1: small, it rings at w^2 = gamma p A^2 / (m V), gamma = cp/cv of
    // the McBride (2002) N2 at 300 K. Its swing a = v A / (w V), 1.4 % of
    // the volume, lengthens the period by (gamma + 1)(2 gamma - 1) a^2 / 48
    // (Landau and Lifshitz, Mechanics, section 28, for p V^gamma held),
    // 1.8e-5, and the midpoint step by (w dt)^2 / 12, 1.0e-6. The terms left
    // out are near a^4, about 4e-8; the probe run agreed to 1.6e-8.
    let (code, json, stderr) = run_binary(&["--json", "system", &example("n2_gas_spring.json")]);
    assert_eq!(code, 0, "{stderr}");
    let a = n2_nasa9_low();
    let cp = n2_heat_capacity(&a, 300.0);
    let gamma = cp / (cp - 1.0);
    let (mass, area, p0, v0, speed, dt) = (1.0, 1.0e-3, 90_000.0, 1.0e-4, 0.05, 1.0e-4);
    let w = (gamma * p0 * area * area / (mass * v0)).sqrt();
    let swing = speed * area / (w * v0);
    let period = core::f64::consts::TAU / w
        * (1.0
            + (gamma + 1.0) * (2.0 * gamma - 1.0) * swing * swing / 48.0
            + (w * dt).powi(2) / 12.0);
    let mut crossings = Vec::new();
    let mut previous: Option<(f64, f64)> = None;
    for sample in json.split("{\"time_s\":").skip(1) {
        let t: f64 = sample.split_once(',').unwrap().0.parse().unwrap();
        let v = value(sample, "", "piston.velocity_m_s");
        if let Some((t_old, v_old)) = previous
            && v_old > 0.0
            && v <= 0.0
        {
            crossings.push(t_old + (t - t_old) * v_old / (v_old - v));
        }
        previous = Some((t, v));
    }
    assert!(crossings.len() >= 5, "{crossings:?}");
    let measured = (crossings[crossings.len() - 1] - crossings[0]) / (crossings.len() - 1) as f64;
    assert!(
        (measured / period - 1.0).abs() < 1.0e-6,
        "period {measured} vs {period}"
    );
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-6,
        "{json}"
    );
    let stored = value(&json, "\"totals\"", "stored_change_j");
    let given = value(&json, "\"totals\"", "boundary_energy_in_j");
    assert!((stored - given).abs() < 1.0e-6, "{stored} vs {given}");
}

/// Claims, declared before the first run: with `adaptive`, a gas cavity's
/// volume and the temperature of its gas are error accounts. On the N2 gas
/// spring over 0.15 s, with every other tolerance too loose to bind,
/// tightening either of the two alone takes more accepted steps and rejects
/// a comparison, and tightening both ends closer to a run at a fixed small
/// step.
/// A scene without the temperature's tolerance is refused by that account's
/// name.
#[test]
fn a_gas_cavitys_volume_and_temperature_are_adaptive_error_accounts() {
    // The device of examples/system/n2_gas_spring.json.
    let spring = |steps: &str| {
        format!(
            r#"{{
  "schema": "frankensim.system.v0", "duration_s": 0.15, {steps}
  "elements": {{
    "spring": {{"kind": "gas_cavity", "gas": "thermochemical 300 90000 mole 1 N2 1",
               "sources": {sources}, "volume_m3": 1.0e-4, "insulated": true,
               "temperature_min_k": 200, "temperature_max_k": 900}},
    "piston": {{"kind": "mass", "mass_kg": 1.0, "velocity_m_s": 0.05}},
    "outside": {{"kind": "force", "force_n": -90.0}}
  }},
  "connections": [
    {{"join": ["piston.p", {{"port": "spring.p", "ratio": 1.0e-3}}, "outside.p"]}}
  ]
}}"#,
            sources = n2_sources()
        )
    };
    // The cavity's two accounts with the size of what each measures. A tight
    // one has 1e-8 of that for its tolerance and a loose one a hundredth.
    let sizes = [
        ("spring.volume_m3", 1.0e-4),
        ("spring.temperature_k", 300.0),
    ];
    let accounts = |tight: &dyn Fn(usize) -> bool| {
        let mut all = vec![
            ("piston.position_m", 1.0e3),
            ("piston.velocity_m_s", 1.0e3),
            ("outside.energy_in_j", 1.0e3),
        ];
        for (index, (name, size)) in sizes.iter().enumerate() {
            let part = if tight(index) { 1.0e-8 } else { 1.0e-2 };
            all.push((*name, part * size));
        }
        all
    };
    let run = |name: &str, tight: &dyn Fn(usize) -> bool| {
        controlled_run(name, &spring(&error_controlled(5.0e-3, &accounts(tight))))
    };
    let (loose, loose_steps, _) = run("gas_spring_loose.json", &|_| false);
    for (index, (account, _)) in sizes.iter().enumerate() {
        let name = format!("gas_spring_tight_{index}.json");
        let (_, steps, rejected) = run(&name, &|which| which == index);
        eprintln!(
            "gas spring: {account} alone, {steps} steps and {rejected} rejected, loose {loose_steps}"
        );
        assert!(steps > loose_steps, "{account}: {steps} vs {loose_steps}");
        assert!(rejected >= 1.0, "{account}: no comparison rejected");
    }
    let (tight, _, _) = run("gas_spring_tight.json", &|_| true);

    let ends = |json: &str| {
        ["volume_m3", "temperature_k"].map(|field| value(json, &element("spring"), field))
    };
    let path = scratch("gas_spring_fixed.json", &spring("\"step_s\": 2.5e-5,"));
    let (code, fixed, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let reference = ends(&fixed);
    let distance = |json: &str| {
        let parts = ends(json).into_iter().zip(reference);
        parts
            .map(|(end, fixed)| (end / fixed - 1.0).abs())
            .fold(0.0, f64::max)
    };
    let (far, near) = (distance(&loose), distance(&tight));
    eprintln!("gas spring: from the fixed run, loose {far:e}, tight {near:e}");
    assert!(near < 0.1 * far, "tight {near}, loose {far}");

    let partial = accounts(&|_| true);
    let without: Vec<_> = partial
        .iter()
        .copied()
        .filter(|(name, _)| !name.ends_with("_k"))
        .collect();
    let (code, stderr) = refused(
        "gas_spring_partial.json",
        &spring(&error_controlled(5.0e-3, &without)),
    );
    assert_eq!(code, 3, "{stderr}");
    let missing = "adaptive.tolerances.spring.temperature_k: missing explicit tolerance";
    assert!(stderr.contains(missing), "{stderr}");

    // A gas cavity's step measures order two: an order above it is refused.
    let third = error_controlled(5.0e-3, &partial).replace("\"order\": 2", "\"order\": 3");
    let (code, stderr) = refused("gas_spring_order.json", &spring(&third));
    assert_eq!(code, 3, "{stderr}");
    let order = "adaptive.order: spring is a gas_cavity, whose step is of order 2; \
                 the scene declares order 3";
    assert!(stderr.contains(order), "{stderr}");
}

#[test]
fn gas_cavity_scenes_refuse_bad_records_packs_and_controls_and_gas_outside_the_cards() {
    let base = heated_n2();
    for (name, from, to, reason) in [
        (
            "record.json",
            "thermochemical 250 80000",
            "mixture 250 80000",
            "require a thermochemical record",
        ),
        (
            "pack.json",
            "species.pack",
            "absent.pack",
            "sources[0].species",
        ),
        (
            "insulated.json",
            "\"volume_m3\": 1.0e-3",
            "\"volume_m3\": 1.0e-3, \"insulated\": 1",
            "insulated is true or false",
        ),
        (
            "iterations.json",
            "\"volume_m3\": 1.0e-3",
            "\"volume_m3\": 1.0e-3, \"max_iterations\": 0",
            "max_iterations is a whole number from 1 to 128",
        ),
        (
            "field.json",
            "\"volume_m3\": 1.0e-3",
            "\"volume_m3\": 1.0e-3, \"wall\": true",
            "unsupported field \"wall\"",
        ),
        // 150 kPa is outside the Lemmon (2004) transport card's 75 to
        // 101.325 kPa, which the parcel resolves though it uses no transport.
        (
            "window.json",
            "thermochemical 250 80000",
            "thermochemical 250 150000",
            "NoClaimInDomain",
        ),
        (
            "adaptive.json",
            "\"duration_s\": 1.0",
            "\"adaptive\": {\"tolerances\": {}}, \"duration_s\": 1.0",
            "adaptive.tolerances.cavity.volume_m3: missing explicit tolerance",
        ),
    ] {
        assert!(base.contains(from), "{name}: test fixture changed");
        let (code, stderr) = refused(name, &base.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
    // In the run: compressing the gas past 101.325 kPa leaves the transport
    // card, and a caloric bracket reaching 1000 K spans both McBride charts;
    // the step that would do either is refused.
    for (name, scene, reason) in [
        (
            "squeezed.json",
            base.replace("thermochemical 250 80000", "thermochemical 250 100000")
                .replace("\"flow_m3_s\": 0", "\"flow_m3_s\": 1.0e-4"),
            "NoClaimInDomain",
        ),
        (
            "bracket.json",
            base.replace("\"temperature_max_k\": 900", "\"temperature_max_k\": 1000"),
            "within one NASA-9 region",
        ),
    ] {
        assert!(scene != base, "{name}: test fixture changed");
        let (code, stderr) = refused(name, &scene);
        assert_eq!(code, 4, "{name}: {stderr}");
        assert!(stderr.contains("refused in Propose"), "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
}

/// A sliding interface under 10 N at the slip rate `slip` [m/s] with `law`,
/// its heat into a 2 J/K body at 300 K, for 1 s in steps of 0.01 s.
fn rubbed_body(law: &str, slip: &str) -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": 1.0, "step_s": 0.01,
  "elements": {{
    "rub": {{"kind": "sliding_interface", "normal_n": 10.0, "law": {law},
            "slip_rate_m_s": {slip}}},
    "body": {{"kind": "heat_capacity", "capacity_j_k": 2.0, "temperature_k": 300.0}}
  }},
  "connections": [{{"join": ["body.heat", "rub.heat"]}}]
}}"#
    )
}

#[test]
fn sliding_interfaces_turn_the_work_of_a_prescribed_slip_into_the_heat_a_body_keeps() {
    // G1/G3: at a constant slip v under N, the law's traction is mu(|v|) N
    // along the slip, and each step dissipates mu N |v| dt, all of it into
    // the body: T1 = T0 + mu N |v| t / C, and the interface's work is the
    // body's heat. Coulomb slides at its kinetic coefficient; the cited
    // rosin fit, at 0.05 m/s along -y, at 0.2 (v + 0.011)^-0.4.
    for (name, law, slip, rate, mu) in [
        (
            "coulomb.json",
            r#"{"coulomb": {"static_mu": 0.4, "kinetic_mu": 0.3}}"#,
            "[0.2, 0.0]",
            [0.2, 0.0],
            0.3,
        ),
        (
            "rosin.json",
            r#""rosin_woodhouse_galluzzo_2025_eq13""#,
            "[0.0, -0.05]",
            [0.0, -0.05],
            0.2 * (0.05_f64 + 0.011).powf(-0.4),
        ),
    ] {
        let path = scratch(name, &rubbed_body(law, slip));
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let speed = f64::hypot(rate[0], rate[1]);
        let power = mu * 10.0 * speed;
        let rub = element("rub");
        for (field, expected) in [
            ("work_j", power),
            ("heat_w", power),
            ("traction_x_n", mu * 10.0 * rate[0] / speed),
            ("traction_y_n", mu * 10.0 * rate[1] / speed),
        ] {
            let found = value(&json, &rub, field);
            assert!(
                (found - expected).abs() <= 1.0e-12,
                "{name}: rub.{field}: {found} vs {expected}"
            );
        }
        let kelvin = value(&json, &element("body"), "temperature_k");
        let expected = 300.0 + power / 2.0;
        assert!(
            (kelvin - expected).abs() <= 1.0e-10,
            "{name}: body {kelvin} vs {expected}"
        );
        assert!(
            value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-12,
            "{name}: {json}"
        );
        assert!(
            value(&json, "\"totals\"", "stored_change_j").abs() < 1.0e-12,
            "{name}: {json}"
        );
    }
}

#[test]
fn a_dragged_lining_heats_its_pad_as_the_stribeck_law_and_the_midpoint_rule_say() {
    // examples/system/dragged_lining.json: a Stribeck lining, with the
    // coefficients of the friction_interface tests, dragged at [0.3, 0.4]
    // m/s (0.5 m/s) under 40 N. G1: every step dissipates P = mu(0.5) N 0.5,
    // and the traction is mu N along the slip. The heat warms a 20 J/K pad
    // that loses heat through 0.5 W/K to 293.15 K. The midpoint rule shrinks
    // T - T* by (1 - g) / (1 + g) a step, with g = G dt / (2 C) and
    // T* = Ta + P / G. G3: the pad's heat plus the wall's is the work.
    let (code, json, stderr) = run_binary(&["--json", "system", &example("dragged_lining.json")]);
    assert_eq!(code, 0, "{stderr}");
    let (normal, speed, capacity, conductance, ambient, dt): (f64, f64, f64, f64, f64, f64) =
        (40.0, 0.5, 20.0, 0.5, 293.15, 0.5);
    let mu = 0.2 + 0.3 * (-(speed / 0.5) * (speed / 0.5)).exp() + 0.05 * speed;
    let power = mu * normal * speed;
    let g = conductance * dt / (2.0 * capacity);
    let settled = ambient + power / conductance;
    let pad = settled + (ambient - settled) * ((1.0 - g) / (1.0 + g)).powi(120);
    let lining = element("lining");
    for (field, expected, tolerance) in [
        ("work_j", power * 60.0, 1.0e-10),
        ("heat_w", power, 1.0e-12),
        ("traction_x_n", mu * normal * 0.6, 1.0e-12),
        ("traction_y_n", mu * normal * 0.8, 1.0e-12),
    ] {
        let found = value(&json, &lining, field);
        assert!(
            (found - expected).abs() <= tolerance,
            "lining.{field}: {found} vs {expected}"
        );
    }
    let kelvin = value(&json, &element("pad"), "temperature_k");
    assert!((kelvin - pad).abs() <= 1.0e-10, "pad {kelvin} vs {pad}");
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-10,
        "{json}"
    );
    // The wall's heat left through its held end: the work the pad did not keep.
    let kept = capacity * (kelvin - ambient);
    let left = -value(&json, "\"totals\"", "boundary_energy_in_j");
    assert!(
        (kept + left - power * 60.0).abs() < 1.0e-9,
        "{kept} + {left} vs {}",
        power * 60.0
    );
}

#[test]
fn sliding_interface_scenes_refuse_bad_slips_fields_laws_loads_and_error_control() {
    let base = rubbed_body(
        r#"{"coulomb": {"static_mu": 0.4, "kinetic_mu": 0.3}}"#,
        "[0.2, 0.0]",
    );
    for (name, from, to, reason) in [
        (
            "rate.json",
            "\"slip_rate_m_s\": [0.2, 0.0]",
            "\"slip_rate_m_s\": 0.2",
            "slip_rate_m_s is two finite numbers",
        ),
        (
            "axes.json",
            "\"slip_rate_m_s\": [0.2, 0.0]",
            "\"slip_rate_m_s\": [0.2, 0.0, 0.0]",
            "slip_rate_m_s is two finite numbers",
        ),
        (
            "field.json",
            "\"normal_n\": 10.0",
            "\"normal_n\": 10.0, \"slip_m_s\": 0.0",
            "unsupported field \"slip_m_s\"",
        ),
        (
            "impulse.json",
            "\"normal_n\": 10.0",
            "\"normal_n\": 10.0, \"solve\": {\"impulse_tolerance_n_s\": 1e-10}",
            "unsupported field \"impulse_tolerance_n_s\"",
        ),
        (
            "law.json",
            r#"{"coulomb": {"static_mu": 0.4, "kinetic_mu": 0.3}}"#,
            "\"rosin\"",
            "the named law is",
        ),
        (
            "normal.json",
            "\"normal_n\": 10.0",
            "\"normal_n\": -1.0",
            "invalid or non-finite physical input: normal_force",
        ),
        (
            "adaptive.json",
            "\"duration_s\": 1.0",
            "\"adaptive\": {\"tolerances\": {}}, \"duration_s\": 1.0",
            "rub is a sliding_interface",
        ),
    ] {
        assert!(base.contains(from), "{name}: test fixture changed");
        let (code, stderr) = refused(name, &base.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
}

/// A 2 kg puck (principal inertias 0.01, 0.01, 0.02 kg m2) on a Hertz patch
/// of 1 cm under 10 N, static 0.4 and kinetic 0.3, at its centre of mass,
/// with the velocity `velocity` [m/s] and the spin `spin` [rad/s] and the
/// patch's starting motion `motion`; the patch's heat into a 2 J/K table at
/// 300 K, for `duration` s in steps of 0.01 s. The map gives the patch the
/// puck's velocity along x and y and the radius times its spin about z.
fn patched_puck(velocity: &str, spin: &str, motion: &str, duration: f64) -> String {
    format!(
        r#"{{
  "schema": "frankensim.system.v0", "duration_s": {duration}, "step_s": 0.01,
  "elements": {{
    "puck": {{"kind": "rigid_network", "torque_ports": [0], "connectors": [],
             "bodies": [{{"mass_kg": 2.0, "principal_inertia_kg_m2": [0.01, 0.01, 0.02],
                         "position_m": [0.0, 0.0, 0.0], "velocity_m_s": {velocity},
                         "angular_velocity_rad_s": {spin}}}]}},
    "patch": {{"kind": "hertz_patch", "normal_n": 10.0, "radius_m": 0.01,
              "static_mu": 0.4, "kinetic_mu": 0.3, "motion_m_s": {motion}}},
    "table": {{"kind": "heat_capacity", "capacity_j_k": 2.0, "temperature_k": 300.0}}
  }},
  "connections": [
    {{"power_map": [[{{"port": "puck.force0",
                     "matrix": [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 0.0]]}},
                    {{"port": "puck.torque0",
                     "matrix": [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -0.01]]}}],
                   "patch.wrench"]}},
    {{"join": ["table.heat", "patch.heat"]}}
  ]
}}"#
    )
}

#[test]
fn a_hertz_patch_slows_a_sliding_puck_by_coulomb_and_a_spinning_one_by_the_hertz_spin_moment() {
    // G1: a pure slide on a patch of uniform kinetic coefficient is resisted
    // by mu N against the slip, so the 2 kg puck slows at 1.5 m/s2 from
    // 1 m/s: after 0.5 s it moves at 0.25 m/s and has gone 0.3125 m. A pure
    // spin is resisted by the moment of mu p over the Hertz pressure p =
    // 3 N / (2 pi a2) sqrt(1 - r2 / a2), mu N a 3 pi / 16, so from 20 rad/s
    // it slows at M / I. The midpoint rule is exact for a constant wrench.
    // G3: every joule the puck loses warms the table: T = T0 + mu N x / C,
    // or T0 + M theta / C.
    let (mu, normal, radius, mass, inertia, capacity, time): (f64, f64, f64, f64, f64, f64, f64) =
        (0.3, 10.0, 0.01, 2.0, 0.02, 2.0, 0.5);
    let path = scratch(
        "slide.json",
        &patched_puck(
            "[1.0, 0.0, 0.0]",
            "[0.0, 0.0, 0.0]",
            "[1.0, 0.0, 0.0]",
            time,
        ),
    );
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "slide: {stderr}");
    let slowing = mu * normal / mass;
    let distance = time - 0.5 * slowing * time * time;
    let (patch, puck) = (element("patch"), element("puck"));
    for (anchor, field, expected) in [
        (&patch, "force_x_n", -mu * normal),
        (&patch, "force_y_n", 0.0),
        (&patch, "spin_moment_nm", 0.0),
        (&patch, "slip_x_m_s", 1.0 - slowing * time),
        (&patch, "owner_steps", 1.0),
        (&puck, "body0.velocity_x_m_s", 1.0 - slowing * time),
        (&puck, "body0.position_x_m", distance),
        (
            &element("table"),
            "temperature_k",
            300.0 + mu * normal * distance / capacity,
        ),
    ] {
        let found = value(&json, anchor, field);
        assert!(
            (found - expected).abs() <= 1.0e-9,
            "slide: {anchor} {field}: {found} vs {expected}"
        );
    }
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-10,
        "slide: {json}"
    );
    let path = scratch(
        "spin.json",
        &patched_puck(
            "[0.0, 0.0, 0.0]",
            "[0.0, 0.0, 20.0]",
            "[0.0, 0.0, 0.2]",
            time,
        ),
    );
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "spin: {stderr}");
    let moment = mu * normal * radius * 3.0 * std::f64::consts::PI / 16.0;
    let angle = 20.0 * time - 0.5 * moment / inertia * time * time;
    for (anchor, field, expected) in [
        (&patch, "spin_moment_nm", -moment),
        (&patch, "force_x_n", 0.0),
        (&patch, "force_y_n", 0.0),
        (
            &patch,
            "spin_speed_m_s",
            radius * (20.0 - moment / inertia * time),
        ),
        (
            &puck,
            "body0.angular_velocity_z_rad_s",
            20.0 - moment / inertia * time,
        ),
        (&puck, "body0.velocity_x_m_s", 0.0),
        (
            &element("table"),
            "temperature_k",
            300.0 + moment * angle / capacity,
        ),
    ] {
        let found = value(&json, anchor, field);
        assert!(
            (found - expected).abs() <= 1.0e-9,
            "spin: {anchor} {field}: {found} vs {expected}"
        );
    }
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-10,
        "spin: {json}"
    );
}

/// The wrench on side a of local Coulomb friction mu over the Hertz
/// pressure of a patch of radius a under N, at the slip `slip` [m/s] and the
/// spin `spin` [rad/s]: the slip at the point rho is slip + spin z x rho.
/// Midpoint quadrature on a polar grid.
fn hertz_patch_wrench(mu: f64, normal: f64, radius: f64, slip: [f64; 2], spin: f64) -> [f64; 3] {
    let (rings, spokes) = (400, 800);
    let mut wrench = [0.0; 3];
    for ring in 0..rings {
        let r = (f64::from(ring) + 0.5) * radius / f64::from(rings);
        let pressure = 3.0 * normal / (2.0 * std::f64::consts::PI * radius * radius)
            * (1.0 - (r / radius) * (r / radius)).sqrt();
        let area = r * (radius / f64::from(rings)) * (std::f64::consts::TAU / f64::from(spokes));
        for spoke in 0..spokes {
            let angle = (f64::from(spoke) + 0.5) * std::f64::consts::TAU / f64::from(spokes);
            let (x, y) = (r * angle.cos(), r * angle.sin());
            let (sx, sy) = (slip[0] - spin * y, slip[1] + spin * x);
            let speed = sx.hypot(sy);
            let (fx, fy) = (
                -mu * pressure * area * sx / speed,
                -mu * pressure * area * sy / speed,
            );
            wrench[0] += fx;
            wrench[1] += fy;
            wrench[2] += x * fy - y * fx;
        }
    }
    wrench
}

#[test]
fn a_spun_puck_feels_the_wrench_local_coulomb_friction_over_the_hertz_pressure_gives() {
    // examples/system/spun_puck.json: a 0.5 kg puck sliding at 2 m/s and
    // spinning at 100 rad/s on a 1 cm patch under 5 N, so the patch's two
    // motions are alike (a w = 1 m/s). G1: the owner's kinetic wrench is
    // local Coulomb friction mu p against the local slip, integrated over
    // the Hertz pressure; at the last step's mean motion, the mean of its
    // two end samples, a 400 x 800 polar quadrature of that integral gives
    // the force and moment to 1e-3. The joint motion lowers both below their
    // pure-slide and pure-spin values. G3: the puck's lost kinetic energy is
    // the ice's heat.
    let (code, json, stderr) = run_binary(&["--json", "system", &example("spun_puck.json")]);
    assert_eq!(code, 0, "{stderr}");
    let samples: Vec<&str> = json.split("{\"time_s\":").skip(1).collect();
    // A sample at the start and one after each of the 60 steps.
    assert_eq!(samples.len(), 61, "{json}");
    let motion = |sample: &str| {
        ["slip_x_m_s", "slip_y_m_s", "spin_speed_m_s"]
            .map(|field| value(sample, "\"values\"", &format!("patch.{field}")))
    };
    let (end, start) = (motion(samples[60]), motion(samples[59]));
    let mean: [f64; 3] = std::array::from_fn(|i| 0.5 * (end[i] + start[i]));
    let (mu, normal, radius): (f64, f64, f64) = (0.3, 5.0, 0.01);
    let expected = hertz_patch_wrench(mu, normal, radius, [mean[0], mean[1]], mean[2] / radius);
    let patch = element("patch");
    let found =
        ["force_x_n", "force_y_n", "spin_moment_nm"].map(|field| value(&json, &patch, field));
    let force_scale = mu * normal;
    let moment_scale = mu * normal * radius * 3.0 * std::f64::consts::PI / 16.0;
    for (index, scale) in [(0, force_scale), (1, force_scale), (2, moment_scale)] {
        assert!(
            (found[index] - expected[index]).abs() <= 1.0e-3 * scale,
            "wrench {index}: {found:?} vs {expected:?} at {mean:?}"
        );
    }
    assert!(found[0] < 0.0 && found[2] < 0.0, "{found:?}");
    assert!(found[0].hypot(found[1]) < 0.99 * force_scale, "{found:?}");
    assert!(-found[2] < 0.99 * moment_scale, "{found:?}");
    let puck = element("puck");
    let speed = |axis: &str| value(&json, &puck, &format!("body0.velocity_{axis}_m_s"));
    let spin = |axis: &str| {
        value(
            &json,
            &puck,
            &format!("body0.angular_velocity_{axis}_rad_s"),
        )
    };
    let before = 0.5 * 0.5 * 4.0 + 0.5 * 0.00025 * 100.0 * 100.0;
    let after = 0.5 * 0.5 * (speed("x").powi(2) + speed("y").powi(2) + speed("z").powi(2))
        + 0.5
            * (0.00015 * spin("x").powi(2)
                + 0.00015 * spin("y").powi(2)
                + 0.00025 * spin("z").powi(2));
    let heat = 1.0 * (value(&json, &element("ice"), "temperature_k") - 268.15);
    assert!(
        (before - after - heat).abs() < 1.0e-9,
        "lost {} vs heat {heat}",
        before - after
    );
    assert!(
        value(&json, "\"totals\"", "energy_residual_j").abs() < 1.0e-10,
        "{json}"
    );
}

#[test]
fn hertz_patch_scenes_refuse_bad_patches_maps_and_budgets_and_a_reversing_slip() {
    let base = patched_puck("[1.0, 0.0, 0.0]", "[0.0, 0.0, 0.0]", "[1.0, 0.0, 0.0]", 0.5);
    for (name, from, to, reason) in [
        (
            "order.json",
            "\"static_mu\": 0.4",
            "\"static_mu\": 0.2",
            "static_mu_below_kinetic_mu",
        ),
        (
            "radius.json",
            "\"radius_m\": 0.01",
            "\"radius_m\": 0.0",
            "invalid or non-finite physical input: contact_radius_m",
        ),
        (
            "normal.json",
            "\"normal_n\": 10.0",
            "\"normal_n\": -1.0",
            "invalid or non-finite physical input: normal_force_n",
        ),
        (
            "motion.json",
            "\"motion_m_s\": [1.0, 0.0, 0.0]",
            "\"motion_m_s\": [1.0, 0.0]",
            "motion_m_s is three finite numbers",
        ),
        (
            "budget.json",
            "\"motion_m_s\": [1.0, 0.0, 0.0]",
            "\"motion_m_s\": [1.0, 0.0, 0.0], \"solve\": {\"maximum_iterations\": 0}",
            "maximum_iterations is a whole number from 1 to 4096",
        ),
        (
            "field.json",
            "\"motion_m_s\": [1.0, 0.0, 0.0]",
            "\"motion_m_s\": [1.0, 0.0, 0.0], \"solve\": {\"tolerance\": 1.0}",
            "unsupported field \"tolerance\"",
        ),
        (
            "shape.json",
            "[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, -0.01]]",
            "[[0.0, 0.0, -0.01]]",
            "power_map[0][1]: matrix row count must be 3",
        ),
        (
            "matrix.json",
            "\"patch.wrench\"]}",
            "\"patch.wrench\"], \"matrix\": [[1.0]]}",
            "a power_map of members has a matrix per member",
        ),
        (
            "adaptive.json",
            "\"duration_s\": 0.5",
            "\"adaptive\": {\"tolerances\": {}}, \"duration_s\": 0.5",
            "patch is a hertz_patch",
        ),
    ] {
        assert!(base.contains(from), "{name}: test fixture changed");
        let (code, stderr) = refused(name, &base.replace(from, to));
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
    // The puck would stop at 2/3 s: the step in which the slip reverses is
    // refused, not shortened.
    let (code, stderr) = refused(
        "arrest.json",
        &patched_puck("[1.0, 0.0, 0.0]", "[0.0, 0.0, 0.0]", "[1.0, 0.0, 0.0]", 1.0),
    );
    assert_eq!(code, 4, "arrest: {stderr}");
    assert!(stderr.contains("refused in Propose"), "arrest: {stderr}");
    assert!(stderr.contains("hertz_patch_reversal"), "arrest: {stderr}");
}

fn array_values(json: &str, anchor: &str, key: &str) -> Vec<f64> {
    let tail = &json[json.find(anchor).unwrap()..];
    let field = format!("\"{key}\":[");
    let entries = tail
        .split_once(&field)
        .unwrap()
        .1
        .split_once(']')
        .unwrap()
        .0;
    if entries.is_empty() {
        Vec::new()
    } else {
        entries
            .split(',')
            .map(|entry| entry.parse().unwrap())
            .collect()
    }
}

#[test]
fn a_gas_boundary_reports_measured_mass_and_constituent_transport_with_full_coverage() {
    // G1: 5 g/s for four seconds delivers 20 g of gas and of its sole
    // constituent. This checks measured storage and transport, not a
    // residual synthesized from the receipt's energy account.
    let path = scratch("vessel_accounts.json", FILLED_VESSEL);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let steps = value(&json, "\"ledger_steps\"", "ledger_steps");
    assert!(steps > 0.0);
    assert_eq!(steps, value(&json, "\"steps\"", "steps"));
    let mass = &json[json.find("\"mass\":{").unwrap()..];
    let carried = &json[json.find("\"carried_kg\":[{").unwrap()..];
    for account in [mass, carried] {
        for field in ["port_in_kg", "stored_change_kg"] {
            assert!(
                (value(account, "", field) - 0.02).abs() < 1.0e-12,
                "{field}: {account}"
            );
        }
        assert!(value(account, "", "balance_defect_kg").abs() < 1.0e-12);
        assert_eq!(value(account, "", "audited_steps"), steps);
    }
    let boundary = "{\"port\":\"vessel.stream\"";
    assert!((value(&json, boundary, "mass_in_kg") - 0.02).abs() < 1.0e-12);
    let carried_in = array_values(&json, boundary, "carried_in_kg");
    assert_eq!(carried_in.len(), 1);
    assert!((carried_in[0] - 0.02).abs() < 1.0e-12);
    assert!(value(&json, "\"totals\"", "mass_residual_kg").abs() < 1.0e-12);
    let residuals = array_values(&json, "\"totals\"", "carried_residual_kg");
    assert_eq!(residuals.len(), 1);
    assert!(residuals[0].abs() < 1.0e-12);

    let (code, text, stderr) = run_binary(&["system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    for expected in [
        "stream ledger [kg]",
        "vessel mass:",
        "vessel carried[0]:",
        "boundary vessel.stream: mass in",
        "stream mass residual:",
        "carried[0] residual:",
    ] {
        assert!(text.contains(expected), "missing {expected}: {text}");
    }
    assert!(text.contains(&format!("{steps:.0}/{steps:.0} steps")));
}

#[test]
fn a_restriction_reports_each_carried_coordinate_and_the_connection_defect() {
    // G0/G3: two differently weighted inlet constituents retain their
    // coordinate order and separately close through the stateless valve.
    let scene = FILLED_VESSEL
        .replace("\"carried\": [0.0]}",
            "\"carried\": [0.0, 0.0]},\n    \"inlet\": {\"kind\": \"restriction\", \"conductance_kg_s_pa\": 2.0e-7, \"carried\": 2}")
        .replace("{\"prescribe\": \"vessel.stream\", \"stream\": [5.0e-3, 1506.75, 5.0e-3]}",
            "{\"join\": [\"vessel.stream\", \"inlet.b\"]},\n    {\"prescribe\": \"inlet.a\", \"contents\": [2.0e5, 301350.0, 0.25, 0.75]}");
    let path = scratch("constituent_accounts.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let initial_mass = 1.0e5 * 0.02 / (287.0 * 300.0);
    let added = value(&json, &element("vessel"), "mass_kg") - initial_mass;
    assert!(added > 0.0);
    let boundary = "{\"port\":\"inlet.a\"";
    assert!((value(&json, boundary, "mass_in_kg") - added).abs() < 1.0e-12);
    let crossing = array_values(&json, boundary, "carried_in_kg");
    assert_eq!(crossing.len(), 2);
    let vessel = &json[json.find(&element("vessel")).unwrap()..];
    let accounts = vessel
        .split_once("\"carried_kg\":[")
        .unwrap()
        .1
        .split_once(']')
        .unwrap()
        .0
        .split("},{")
        .collect::<Vec<_>>();
    assert_eq!(accounts.len(), 2);
    for (index, fraction) in [0.25, 0.75].into_iter().enumerate() {
        assert!((crossing[index] - fraction * added).abs() < 1.0e-12);
        let stored = value(&json, &element("vessel"), &format!("carried{index}"))
            * value(&json, &element("vessel"), "mass_kg");
        assert!((stored - crossing[index]).abs() < 1.0e-12);
        assert!((value(accounts[index], "", "stored_change_kg") - crossing[index]).abs() < 1.0e-12);
        assert_eq!(
            value(accounts[index], "", "audited_steps"),
            value(&json, "\"ledger_steps\"", "ledger_steps")
        );
    }
    let inlet = &json[json.find(&element("inlet")).unwrap()..];
    assert!(value(inlet, "\"mass\"", "stored_change_kg").abs() < 1.0e-12);
    assert_eq!(
        value(inlet, "\"mass\"", "audited_steps"),
        value(&json, "\"ledger_steps\"", "ledger_steps")
    );
    let connection = "{\"kernel_index\":0,";
    assert!(value(&json, connection, "mass_defect_kg").abs() < 1.0e-12);
    let defects = array_values(&json, connection, "carried_defect_kg");
    assert_eq!(defects.len(), 2);
    assert!(defects.iter().all(|value| value.abs() < 1.0e-12));
    let residuals = array_values(&json, "\"totals\"", "carried_residual_kg");
    assert_eq!(residuals.len(), 2);
    assert!(residuals.iter().all(|value| value.abs() < 1.0e-12));
    assert!(value(&json, "\"totals\"", "mass_residual_kg").abs() < 1.0e-12);
}

#[test]
fn a_scene_without_stream_measurements_never_reports_zero_mass_residual() {
    let (code, json, stderr) =
        run_binary(&["--json", "system", &example("driven_oscillator.json")]);
    assert_eq!(code, 0, "{stderr}");
    assert!(json.contains("\"mass_residual_kg\":null"));
    assert!(json.contains("\"has_streams\":false,\"mass\":null,\"carried_kg\":[]"));
    assert!(array_values(&json, "\"totals\"", "carried_residual_kg").is_empty());
}

// Original source laws used by these scenes are authored analytic fixtures,
// not observations of a physical material. All execution uses the real CLI.
fn system_material_pack(cp_slope: f64, semantic_k: bool, pressure_required: bool) -> Vec<u8> {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
        ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
        UncertaintyModel,
    };
    use fs_material::state_point::{SPECIFIC_HEAT_CAPACITY_DIMS, THERMAL_CONDUCTIVITY_DIMS};
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

    let mut claims = ClaimSet::new();
    let source = "authored affine thermal laws for generic System scenes; not measured";
    let provenance = Provenance {
        source: source.into(),
        license: "CC0: authored synthetic test data".into(),
        artifact: None,
    };
    for (name, dims, at_290, slope) in [
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            1000.0,
            cp_slope,
        ),
        ("density", fs_qty::Density::DIMS, 1000.0, -0.5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 2.0, 0.01),
    ] {
        let key = if semantic_k && name == "thermal_conductivity" {
            PropertyKey::with_quantity(
                name,
                fs_qty::QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity,
                    ValueForm::Static,
                )),
            )
        } else {
            PropertyKey::new(name, dims)
        };
        let mut validity = ValidityDomain::unconstrained().with("T", 290.0, 400.0);
        if pressure_required {
            validity = validity.with("pressure", 101_325.0, 101_325.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(290.0, at_290), (400.0, at_290 + 110.0 * slope)],
                    dims,
                },
                validity,
                uncertainty: if name == "specific_heat_capacity" {
                    UncertaintyModel::RelativeHalfWidth {
                        fraction: 0.01,
                        confidence: 0.95,
                    }
                } else {
                    UncertaintyModel::Unstated
                },
                interpolation: InterpolationPolicy::LinearInside,
                observations: Vec::new(),
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    // The pack retains an authored source artifact. Its claims deliberately
    // have no observation link, so they cannot become observation-backed.
    claims
        .register_observation(ObservationDataset {
            specimen: "synthetic affine material".into(),
            method: "analytic test fixture, not an experiment".into(),
            artifact: fs_blake3::hash_bytes(source.as_bytes()),
            caveats: "not measured or experimentally qualified".into(),
            provenance,
        })
        .unwrap();
    let pack = NormalizedPack::new(
        "generic-system-thermal-fixture",
        "test-fixture-v1",
        fs_blake3::hash_bytes(source.as_bytes()),
        "CC0: authored synthetic test data",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "synthetic affine material".into(),
            phase: "solid".into(),
            process: "authored for generic composition tests".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap()
    .to_bytes()
}

fn system_material_declaration(pack: &str, conductivity: bool) -> String {
    let conductivity = if conductivity {
        r#", "conductivity_property":"thermal_conductivity""#
    } else {
        ""
    };
    format!(
        r#"{{
      "pack":"{pack}", "phase":"solid", "temperature_axis":"T",
      "temperature_kind":"legacy", "state":[],
      "temperature_min_k":290.0, "temperature_max_k":400.0,
      "reference_specific_enthalpy_j_kg":0.0,
      "heat_capacity_property":"specific_heat_capacity", "density_property":"density"{conductivity}
    }}"#
    )
}

fn source_body_scene(initial_h: f64, power_w: f64, duration_s: f64) -> String {
    let material = system_material_declaration("material.fsmcdpk", true);
    format!(
        r#"{{
      "schema":"frankensim.system.v0", "duration_s":{duration_s}, "step_s":0.25,
      "materials":{{"metal":{material}}},
      "elements":{{"body":{{"kind":"enthalpy_body", "material":"metal",
        "mass_kg":2.0, "specific_enthalpy_j_kg":{initial_h}}}}},
      "connections":[{{"prescribe":"body.heat", "heat_rate_w":{power_w}}}]
    }}"#
    )
}

fn source_layer_scene() -> String {
    let material = system_material_declaration("material.fsmcdpk", true);
    format!(
        r#"{{
      "schema":"frankensim.system.v0", "duration_s":2.0, "step_s":0.5,
      "materials":{{"metal":{material}}},
      "elements":{{"wall":{{"kind":"layer_conductor", "materials":["metal","metal"],
        "geometry":{{"kind":"planar", "distances_m":[0.1,0.2], "area_m2":0.5}}}}}},
      "connections":[{{"prescribe":"wall.a", "temperature_k":400.0}},
                     {{"prescribe":"wall.b", "temperature_k":300.0}}]
    }}"#
    )
}

fn source_scene_file(name: &str, scene: &str, pack: &[u8]) -> PathBuf {
    let path = scratch(name, scene);
    std::fs::write(path.parent().unwrap().join("material.fsmcdpk"), pack).unwrap();
    path
}

fn reported_material_hash(json: &str, field: &str) -> String {
    let prefix = format!("\"{field}\":\"");
    let hash = json
        .split_once(&prefix)
        .unwrap()
        .1
        .split_once('"')
        .unwrap()
        .0;
    assert_eq!(hash.len(), 64, "{field}: {hash}");
    assert!(hash.bytes().all(|byte| byte.is_ascii_hexdigit()));
    hash.into()
}

#[test]
fn material_card_calorics_drive_live_enthalpy_temperature_and_density() {
    let pack = system_material_pack(2.0, false, false);
    let scene = source_body_scene(0.0, 2000.0, 10.0);
    let path = source_scene_file("card_calorics.json", &scene, &pack);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    // G1: integrate Cp=1000+2(T-290), then solve h=1000x+x².
    let h = 2000.0 * 10.0 / 2.0;
    let x = 2.0 * h / (1000.0 + (1_000_000.0_f64 + 4.0 * h).sqrt());
    let t = 290.0 + x;
    assert!((value(&json, &element("body"), "specific_enthalpy_j_kg") - h).abs() < 1e-8);
    assert!((value(&json, &element("body"), "temperature_k") - t).abs() < 1e-9);
    assert!(
        (value(&json, &element("body"), "bulk_density_kg_m3") - (1000.0 - 0.5 * x)).abs() < 1e-9
    );
    assert!((value(&json, &element("body"), "stored_change_j") - 20_000.0).abs() < 1e-8);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-8);
    // Cp(T)=420+2T, so this is the independent body entropy increase.
    let entropy = 2.0 * (420.0 * (t / 290.0).ln() + 2.0 * x);
    assert!((value(&json, "\"port\":\"body.heat\"", "entropy_in_j_per_k") - entropy).abs() < 1e-9);
    assert!(json.contains("\"uncertainty_propagated\":false"));
    assert!(json.contains("\"kind\":\"relative-half-width\""));
    assert!(json.contains("\"kind\":\"unstated\""));
    assert!(json.contains("\"observation_backed\":false"));
    assert!(json.contains("\"source_receipts\":["));
    let decoded = fs_matdb::NormalizedMaterialCardPack::from_bytes(&pack).unwrap();
    assert_eq!(
        reported_material_hash(&json, "pack_hash"),
        decoded.content_hash().to_hex()
    );
    assert_eq!(
        reported_material_hash(&json, "card_hash"),
        decoded.card().content_hash().to_hex()
    );
}

#[test]
fn material_card_layers_preserve_the_original_affine_conductivity_integral() {
    for semantic in [false, true] {
        let pack = system_material_pack(2.0, semantic, false);
        let mut scene = source_layer_scene();
        if semantic {
            scene = scene.replace("\"conductivity_property\":\"thermal_conductivity\"",
                "\"conductivity_property\":\"thermal_conductivity\",\"conductivity_kind\":\"thermal-conductivity\"");
        }
        let path = source_scene_file("card_conduction.json", &scene, &pack);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "semantic {semantic}: {stderr}");
        // G1: k(T)=0.01T-0.9; Q=A/L integral_300^400 k(T)dT.
        let q = 0.5 / 0.3 * (0.005 * (400.0 * 400.0 - 300.0 * 300.0) - 0.9 * 100.0);
        assert!((value(&json, &element("wall"), "heat_rate_w") - q).abs() < 1e-7);
        assert!((value(&json, "\"port\":\"wall.a\"", "energy_in_j") - 2.0 * q).abs() < 1e-7);
        assert!((value(&json, "\"port\":\"wall.b\"", "energy_in_j") + 2.0 * q).abs() < 1e-7);
        let entropy = 2.0 * q * (1.0 / 300.0 - 1.0 / 400.0);
        assert!((value(&json, "\"totals\"", "entropy_generated_j_per_k") - entropy).abs() < 1e-9);
        assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-9);
    }
}

#[test]
fn source_backed_closed_thermal_network_reaches_its_energy_determined_equilibrium() {
    let material = system_material_declaration("material.fsmcdpk", true);
    let scene = format!(
        r#"{{
      "schema":"frankensim.system.v0", "duration_s":500.0, "step_s":0.5, "sample_every_s":50.0,
      "materials":{{"metal":{material}}},
      "elements":{{
        "hot":{{"kind":"enthalpy_body", "material":"metal", "mass_kg":1.0,
          "specific_enthalpy_j_kg":63600.0}},
        "cold":{{"kind":"enthalpy_body", "material":"metal", "mass_kg":2.0,
          "specific_enthalpy_j_kg":10100.0}},
        "wall":{{"kind":"layer_conductor", "materials":["metal","metal"],
          "geometry":{{"kind":"planar", "distances_m":[0.01,0.02], "area_m2":0.5}}}}
      }},
      "connections":[{{"join":["hot.heat","wall.a"]}},{{"join":["cold.heat","wall.b"]}}]
    }}"#
    );
    let path = source_scene_file(
        "card_thermal_equilibrium.json",
        &scene,
        &system_material_pack(2.0, false, false),
    );
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let initial_energy = 63_600.0 + 2.0 * 10_100.0;
    let h_eq = initial_energy / 3.0;
    let equilibrium = 290.0 + 2.0 * h_eq / (1000.0 + (1_000_000.0_f64 + 4.0 * h_eq).sqrt());
    for name in ["hot", "cold"] {
        assert!((value(&json, &element(name), "temperature_k") - equilibrium).abs() < 1e-6);
    }
    let mut samples = 0;
    for sample in json.split("{\"time_s\":").skip(1) {
        samples += 1;
        let energy = value(sample, "", "hot.specific_enthalpy_j_kg")
            + 2.0 * value(sample, "", "cold.specific_enthalpy_j_kg");
        assert!(
            (energy - initial_energy).abs() < 1e-6,
            "{energy} vs {initial_energy}"
        );
    }
    assert_eq!(samples, 11);
    let entropy_change = |t0: f64| 420.0 * (equilibrium / t0).ln() + 2.0 * (equilibrium - t0);
    let entropy = entropy_change(350.0) + 2.0 * entropy_change(300.0);
    assert!((value(&json, "\"totals\"", "entropy_generated_j_per_k") - entropy).abs() < 1e-7);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-6);
    assert_eq!(value(&json, &element("wall"), "stored_change_j"), 0.0);
}

#[test]
fn material_identity_follows_pack_content_and_admitted_laws_not_filenames() {
    let scene = source_body_scene(0.0, 2000.0, 10.0);
    let bytes = system_material_pack(2.0, false, false);
    let path = source_scene_file("card_identity.json", &scene, &bytes);
    let (code, original, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    std::fs::write(path.parent().unwrap().join("renamed.fsmcdpk"), &bytes).unwrap();
    let renamed = scratch(
        "card_renamed.json",
        &scene.replace("material.fsmcdpk", "renamed.fsmcdpk"),
    );
    let (code, renamed, stderr) = run_binary(&["--json", "system", renamed.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    for field in ["pack_hash", "pack_bytes_hash", "card_hash", "caloric_hash"] {
        assert_eq!(
            reported_material_hash(&original, field),
            reported_material_hash(&renamed, field)
        );
    }
    assert_eq!(
        value(&original, &element("body"), "temperature_k"),
        value(&renamed, &element("body"), "temperature_k")
    );
    // The same path now names a different physical law, not a cached source.
    std::fs::write(
        path.parent().unwrap().join("material.fsmcdpk"),
        system_material_pack(4.0, false, false),
    )
    .unwrap();
    let (code, changed, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    for field in ["pack_hash", "card_hash", "caloric_hash"] {
        assert_ne!(
            reported_material_hash(&original, field),
            reported_material_hash(&changed, field)
        );
    }
    assert!(
        value(&changed, &element("body"), "temperature_k")
            < value(&original, &element("body"), "temperature_k")
    );
}

#[test]
fn material_scene_admission_refuses_unknown_mixed_and_unsupported_sources_by_name() {
    let body = source_body_scene(0.0, 2000.0, 10.0);
    let layer = source_layer_scene();
    let bytes = system_material_pack(2.0, false, false);
    source_scene_file("card_refusal_setup.json", &body, &bytes);
    let cases = [
        ("unknown_body", body.replace("\"material\":\"metal\"", "\"material\":\"missing\""), "missing"),
        ("unknown_layer", layer.replace("[\"metal\",\"metal\"]", "[\"metal\",\"missing\"]"), "missing"),
        ("mixed_body", body.replace("\"material\":\"metal\"", "\"material\":\"metal\",\"phase\":\"solid\""), "body"),
        ("mixed_layer", layer.replace("\"materials\":[\"metal\",\"metal\"]",
            "\"materials\":[\"metal\",\"metal\"],\"conductivity_a\":[[290,1],[400,1]]"), "wall"),
        ("missing_k", layer.replace(", \"conductivity_property\":\"thermal_conductivity\"", ""), "conductivity"),
        ("unknown_cp", body.replace("\"heat_capacity_property\":\"specific_heat_capacity\"",
            "\"heat_capacity_property\":\"missing_cp\""), "missing_cp"),
        ("source_range", body.replace("\"temperature_max_k\":400.0", "\"temperature_max_k\":401.0"), "calorics"),
        ("phase", body.replace("\"phase\":\"solid\"", "\"phase\":\"solid-liquid\""), "phase"),
        ("temperature_kind", body.replace("\"temperature_kind\":\"legacy\"", "\"temperature_kind\":\"temperature-difference\""), "temperature_kind"),
        ("conductivity_kind", layer.replace("\"conductivity_property\":\"thermal_conductivity\"",
            "\"conductivity_property\":\"thermal_conductivity\",\"conductivity_kind\":\"pressure\""), "conductivity_kind"),
        ("temperature_in_context", body.replace("\"state\":[]", "\"state\":[{\"name\":\"T\",\"value\":\"300 K\",\"kind\":\"legacy\"}]"), "state"),
        ("initial_domain", source_body_scene(122_101.0, 0.0, 1.0), "OutsideEnthalpyDomain"),
    ];
    for (name, scene, expected) in cases {
        let (code, stderr) = refused(&format!("card_{name}.json"), &scene);
        assert_eq!(code, 3, "{name}: {stderr}");
        assert!(stderr.contains(expected), "{name}: {stderr}");
    }
    std::fs::write(
        scratch("card_pressure_setup.json", &body)
            .parent()
            .unwrap()
            .join("material.fsmcdpk"),
        system_material_pack(2.0, false, true),
    )
    .unwrap();
    let (code, stderr) = refused("card_missing_pressure.json", &body);
    assert_eq!(code, 3, "{stderr}");
    assert!(
        stderr.contains("materials.metal") && stderr.contains("NoClaimInDomain"),
        "{stderr}"
    );
    let (code, _, stderr) = run_binary(&[
        "--json",
        "system",
        source_scene_file(
            "card_invalid_pack.json",
            &body,
            b"not a normalized material card pack",
        )
        .to_str()
        .unwrap(),
    ]);
    assert_eq!(code, 3, "{stderr}");
    assert!(stderr.contains("materials.metal.pack"), "{stderr}");
}

#[test]
fn material_domain_exit_preserves_native_scene_state_and_original_source_laws() {
    let bytes = system_material_pack(2.0, false, false);
    let text = source_body_scene(122_099.0, 2000.0, 1.0);
    let path = source_scene_file("card_runtime_domain.json", &text, &bytes);
    let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    let before = scene.render_json();
    // Admitted components own their prepared laws: the runtime does not reread
    // the external file, including on an unsuccessful trial.
    std::fs::write(
        path.parent().unwrap().join("material.fsmcdpk"),
        b"changed after admission",
    )
    .unwrap();
    scene_context(|cx, _| {
        let error = scene.advance(0.25, cx).unwrap_err();
        assert_eq!(error.exit_code, fs_cli::exit::REFUSED);
        assert!(error.message.contains("OutsideEnthalpyDomain"), "{error}");
        assert_eq!(scene.render_json(), before);
        assert_eq!(scene.steps(), 0);
        assert_eq!(scene.time_s(), 0.0);
    });
}

#[test]
fn coarse_prescribed_force_keeps_pulse_timing_in_fixed_and_adaptive_scenes() {
    let base = r#"{
  "schema": "frankensim.system.v0", "duration_s": 1.0, "step_s": 1.0,
  "sample_every_s": 0.125,
  "elements": {"body": {"kind": "mass", "mass_kg": 2.0}},
  "connections": [{"prescribe": "body.p", "force_n": {
    "table": [[0,0], [0.125,0], [0.25,4], [0.375,0], [1,0]]
  }}]
}"#;
    for adaptive in [false, true] {
        let scene = if adaptive {
            base.replace(
                "\"elements\"",
                r#""adaptive": {"tolerances": {
              "body.position_m": {"absolute": 1, "relative": 0},
              "body.velocity_m_s": {"absolute": 1, "relative": 0}
            }}, "elements""#,
            )
        } else {
            base.to_owned()
        };
        let path = scratch(&format!("scheduled_force_{adaptive}.json"), &scene);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        // The triangle has impulse 0.5 N s and time centroid 0.25 s.
        // At t=1, v=I/m=0.25 and x=I*(1-centroid)/m=0.1875.
        // Smearing its mean over the requested one-second step would give
        // the right impulse and kinetic energy but the wrong x=0.125.
        assert!((value(&json, &element("body"), "velocity_m_s") - 0.25).abs() < 1e-12);
        assert!((value(&json, &element("body"), "position_m") - 0.1875).abs() < 1e-12);
        assert!((value(&json, "\"totals\"", "stored_change_j") - 0.0625).abs() < 1e-12);
        assert!((value(&json, "\"totals\"", "boundary_energy_in_j") - 0.0625).abs() < 1e-12);
        assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-12);
        let samples: Vec<f64> = json
            .split("{\"time_s\":")
            .skip(1)
            .map(|sample| sample.split_once(',').unwrap().0.parse().unwrap())
            .collect();
        for knot in [0.125, 0.25, 0.375, 1.0] {
            assert!(samples.contains(&knot), "missing knot {knot}: {samples:?}");
        }
        if adaptive {
            assert!(value(&json, "\"adaptive\"", "max_error_norm") <= 1.0);
        } else {
            assert_eq!(value(&json, "\"schema\"", "steps"), 4.0);
        }
    }
}

#[test]
fn coarse_vector_boundary_uses_every_coordinate_knot_and_merges_signed_zero() {
    let scene = r#"{
  "schema": "frankensim.system.v0", "duration_s": 1.0, "step_s": 1.0,
  "elements": {"body": {
    "kind": "rigid_network", "gravity_m_s2": [0,0,0],
    "bodies": [{"mass_kg": 2.0, "principal_inertia_kg_m2": [1,1,1], "position_m": [0,0,0]}],
    "connectors": []
  }},
  "connections": [{"prescribe": "body.force0", "force_n": [
    {"table": [[0,0], [0.125,0], [0.25,4], [0.375,0], [1,0]]},
    {"table": [[-0.0,0], [0.5,0], [0.625,2], [0.75,0], [1,0]]},
    {"table": [[0,0]]}
  ]}]
}"#;
    let path = scratch("scheduled_vector_force.json", scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    // Independent impulse/first-moment integrals on two axes. The y pulse
    // has impulse 0.25 N s and centroid 0.625 s. Its knots are absent from x.
    for (field, expected) in [
        ("body0.velocity_x_m_s", 0.25),
        ("body0.position_x_m", 0.1875),
        ("body0.velocity_y_m_s", 0.125),
        ("body0.position_y_m", 0.046875),
        ("body0.velocity_z_m_s", 0.0),
        ("body0.position_z_m", 0.0),
    ] {
        let actual = value(&json, &element("body"), field);
        assert!(
            (actual - expected).abs() < 1e-11,
            "{field}: {actual} vs {expected}"
        );
    }
    assert_eq!(value(&json, "\"schema\"", "steps"), 7.0);
    assert!((value(&json, "\"totals\"", "stored_change_j") - 0.078125).abs() < 1e-11);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-11);
}

#[test]
fn a_coarse_prescribed_transmission_ratio_matches_its_per_knot_rotations() {
    let scene = std::fs::read_to_string(example("variable_transmission.json"))
        .unwrap()
        .replace("\"step_s\": 0.00390625", "\"step_s\": 2.0");
    let path = scratch("scheduled_transmission.json", &scene);
    let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
    assert_eq!(code, 0, "{stderr}");
    let (inertia, stiffness) = (0.02_f64, 400.0_f64);
    let frequency = (stiffness / inertia).sqrt();
    let (mut phase, mut angle, mut omega, mut extension) = (0.0_f64, 0.0, 1.0, 0.01);
    let (mut effort, mut flow) = (0.0, 0.0);
    // Integrate each affine ratio independently. Midpoint evolution of the
    // coupled conservative state is a rotation in energy coordinates.
    for (dt, ratio) in [(0.5, 0.01), (0.5, 0.0), (1.0, -0.01)] {
        phase += 2.0 * (0.5 * ratio * frequency * dt).atan();
        let next_w = phase.cos() + 0.01 * frequency * phase.sin();
        let next_q = 0.01 * phase.cos() - phase.sin() / frequency;
        angle += dt * 0.5 * (omega + next_w);
        effort = ratio * stiffness * 0.5 * (extension + next_q);
        flow = ratio * 0.5 * (omega + next_w);
        (omega, extension) = (next_w, next_q);
    }
    for (name, field, expected) in [
        ("rotor", "angle_rad", angle),
        ("rotor", "angular_velocity_rad_s", omega),
        ("spring", "extension_m", extension),
        ("transmission", "ratio", -0.01),
        ("transmission", "effort_a", effort),
        ("transmission", "flow_b", flow),
    ] {
        let actual = value(&json, &element(name), field);
        assert!(
            (actual - expected).abs() < 1e-10,
            "{name}.{field}: {actual} vs {expected}"
        );
    }
    assert_eq!(value(&json, "\"schema\"", "steps"), 3.0);
    let reached_w = value(&json, &element("rotor"), "angular_velocity_rad_s");
    let reached_q = value(&json, &element("spring"), "extension_m");
    let energy = 0.5 * (inertia * reached_w.powi(2) + stiffness * reached_q.powi(2));
    assert!((energy - 0.03).abs() < 1e-12);
    assert!(value(&json, "\"totals\"", "energy_residual_j").abs() < 1e-11);
}

#[test]
fn coarse_prescribed_slider_spans_schedule_their_own_knots_and_match_a_source() {
    let scene = std::fs::read_to_string(example("friction_belt.json")).unwrap();
    let original = r#"{"table": [[0.0, 0.0], [0.25, 1.0], [0.5, 3.0], [1.5, 3.0]]}"#;
    let span = r#"{"table": [[0,0], [0.125,0.25], [0.25,0.5], [0.5,0], [1,-0.5], [1.5,0]]}"#;
    // The power boundary describes exactly the same motion but omits a
    // collinear knot. Only the span can schedule t=0.125 in that scene.
    let power = r#"{"table": [[0,0], [0.25,0.5], [0.5,0], [1,-0.5], [1.5,0]]}"#;
    let source = scene
        .replace(original, span)
        .replace("\"step_s\": 0.005", "\"step_s\": 1.5")
        .replace("\"sample_every_s\": 0.25", "\"sample_every_s\": 0.125");
    let boundary = source
        .replace(
            &format!(
                "    \"belt\": {{\"kind\": \"velocity_source\", \"velocity_m_s\": {span}}},\n"
            ),
            "",
        )
        .replace(
            r#"{"join": ["belt.p", "friction.b"]}"#,
            &format!(r#"{{"prescribe": "friction.b", "velocity_m_s": {power}}}"#),
        )
        .replace(
            r#"{"signal": "belt.span", "to": ["friction.span_b"]}"#,
            &format!(r#"{{"prescribe": "friction.span_b", "velocity_m_s": {span}}}"#),
        );
    // Every ramp needs at most m*|a| = 2*2 = 4 N, below the 10 N static
    // limit. The block follows exactly, including reversal, without heat.
    // Positions are independent trapezoidal integrals of the belt speed.
    let expected = [
        (0.0, 0.0, 0.0),
        (0.125, 0.015625, 0.25),
        (0.25, 0.0625, 0.5),
        (0.5, 0.125, 0.0),
        (1.0, 0.0, -0.5),
        (1.5, -0.125, 0.0),
    ];
    let mut projections = Vec::new();
    for (name, scene) in [("source.json", source), ("boundary.json", boundary)] {
        let path = scratch(name, &scene);
        let (code, json, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{name}: {stderr}");
        assert_eq!(value(&json, "\"schema\"", "steps"), 5.0, "{name}");
        let samples: Vec<&str> = json.split("{\"time_s\":").skip(1).collect();
        assert_eq!(samples.len(), expected.len(), "{name}");
        let mut projection = Vec::new();
        for (sample, &(time, position, velocity)) in samples.iter().zip(&expected) {
            let reached: f64 = sample.split_once(',').unwrap().0.parse().unwrap();
            assert_eq!(reached, time, "{name}");
            for (field, wanted) in [
                ("block.position_m", position),
                ("block.velocity_m_s", velocity),
                ("friction.slip_m_s", 0.0),
                ("friction.stuck", 1.0),
                ("lining.temperature_k", 300.0),
            ] {
                let actual = value(sample, "\"values\"", field);
                assert!(
                    (actual - wanted).abs() < 1e-10,
                    "{name} at {time}: {field} = {actual}, expected {wanted}"
                );
                projection.push(actual);
            }
        }
        assert!(value(&json, &element("friction"), "dissipated_j").abs() < 1e-10);
        for field in [
            "stored_change_j",
            "boundary_energy_in_j",
            "energy_residual_j",
        ] {
            assert!(
                value(&json, "\"totals\"", field).abs() < 1e-10,
                "{name}: {field}"
            );
        }
        projections.push(projection);
    }
    for (source, boundary) in projections[0].iter().zip(&projections[1]) {
        assert!((source - boundary).abs() < 1e-10, "{source} vs {boundary}");
    }
}
