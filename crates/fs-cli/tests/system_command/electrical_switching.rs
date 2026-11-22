//! Switched electrical scenes use native Components and preserve every joule.

use super::{example, run_binary, scene_context};

const SIGNAL_CONNECTION: &str =
    r#"{"signal": "switch_command.output", "to": ["switch.closed"]}"#;
const COMMAND: &str = r#"{"held_table": [[0, 1], [0.2, 0], [0.35, 1]]}"#;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * (1.0 + expected.abs()),
        "{actual} != {expected}"
    );
}

fn probe(scene: &fs_cli::SystemScene, path: &str) -> f64 {
    scene.probe(path).unwrap().unwrap()
}

/// Solve C dv/dt = (supply - v) / R by the midpoint rule, independently
/// of the graph, its port conventions and its condensed load-line solve.
fn next_charge(charge: f64, supply: f64, dt: f64) -> f64 {
    let factor = dt / (2.0 * 10.0 * 0.01);
    ((1.0 - factor) * charge + dt * supply / 10.0) / (1.0 + factor)
}

#[test]
fn switching_rc_scene_resolves_events_and_matches_charge_work_heat_and_replay() {
    let path = example("electrical_switching.json");
    let authored = std::fs::read_to_string(&path).unwrap();
    // Both public command paths must use the same point value and event cap.
    for prescribed in [false, true] {
        let text = if prescribed {
            authored.replace(
                SIGNAL_CONNECTION,
                &format!(r#"{{"prescribe": "switch.closed", "value": {COMMAND}}}"#),
            )
        } else {
            authored.clone()
        };
        if prescribed {
            assert_ne!(text, authored);
        }
        let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        for field in ["diode.current_a", "diode.voltage_v", "switch.closed"] {
            assert_eq!(scene.probe(field).unwrap(), None);
        }
        let mut charge = [0.0; 2];
        let mut supplied = [0.0; 2];
        let mut heat = 0.0;
        let mut endpoints = Vec::new();
        let mut saved = None;
        scene_context(|cx, _| {
            while !scene.is_finished() {
                let start = scene.time_s();
                let diode_conducting = !(0.7..1.0).contains(&start);
                let switch_closed = !(0.2..0.35).contains(&start);
                let supply = [
                    if start < 0.7 { 12.0 } else if start < 1.0 { -6.0 } else { 18.0 },
                    12.0,
                ];
                assert!(scene.advance(f64::MAX, cx).unwrap());
                let end = scene.time_s();
                let dt = end - start;
                endpoints.push(end);
                for event in [0.2, 0.35, 0.7, 1.0] {
                    assert!(start >= event - 1e-14 || end <= event + 1e-14);
                }
                for (index, name) in ["diode", "switch"].into_iter().enumerate() {
                    let conducting = if index == 0 { diode_conducting } else { switch_closed };
                    let next = if conducting {
                        next_charge(charge[index], supply[index], dt)
                    } else {
                        charge[index]
                    };
                    let current = (next - charge[index]) / dt;
                    let mean_voltage = (next + charge[index]) / (2.0 * 0.01);
                    let device_voltage = if conducting { 0.0 } else { supply[index] - mean_voltage };
                    supplied[index] += supply[index] * (next - charge[index]);
                    heat += 10.0 * current * current * dt;
                    charge[index] = next;
                    close(probe(&scene, &format!("{name}_capacitor.charge_c")), next);
                    close(probe(&scene, &format!("{name}_capacitor.voltage_v")), next / 0.01);
                    close(probe(&scene, &format!("{name}.current_a")), current);
                    close(probe(&scene, &format!("{name}.voltage_v")), device_voltage);
                    close(probe(&scene, &format!("{name}_resistor.heat_w")), 10.0 * current * current);
                    close(probe(&scene, &format!("{name}_supply.energy_in_j")), -supplied[index]);
                }
                close(probe(&scene, "switch.closed"), if switch_closed { 1.0 } else { 0.0 });
                close(probe(&scene, "case.temperature_k"), 300.0 + heat / 10.0);
                let stored = (charge[0] * charge[0] + charge[1] * charge[1]) / (2.0 * 0.01);
                close(supplied[0] + supplied[1], stored + heat);
                assert!(scene.ledger().energy_residual_j().abs() < 2e-9);
                // Save after the switch has opened and both source clocks
                // have advanced, so continuation covers every later event.
                if saved.is_none() && end > 0.2 {
                    saved = Some(scene.checkpoint().unwrap());
                }
            }
            for event in [0.2, 0.35, 0.7, 1.0] {
                assert!(endpoints.iter().any(|time| (time - event).abs() < 1e-14));
            }
            let (bytes, root) = saved.unwrap();
            let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
            resumed.resume_checkpoint(&bytes, root).unwrap();
            resumed.run(cx).unwrap();
            assert_eq!(resumed.render_json(), scene.render_json());

            // The same knots with linear interpolation are a different
            // problem, even though their values agree at each knot.
            let linear = text.replace(COMMAND, r#"{"table": [[0, 1], [0.2, 0], [0.35, 1]]}"#);
            assert_ne!(linear, text);
            let mut changed = fs_cli::SystemScene::from_text(&linear, &path).unwrap();
            assert!(changed.resume_checkpoint(&bytes, root).is_err());
            assert_eq!(changed.time_s(), 0.0);
        });
        if !prescribed {
            let (code, json, stderr) = run_binary(&["--json", "system", &path]);
            assert_eq!(code, 0, "{stderr}");
            assert_eq!(json, scene.render_json());
        }
    }
}

#[test]
fn fractional_switch_commands_refuse_without_publishing_state_or_accounts() {
    let path = example("electrical_switching.json");
    let original = std::fs::read_to_string(&path).unwrap();
    for (command, successful_steps) in [
        ("0.5", 0),
        (r#"{"table": [[0, 1], [0.2, 0], [0.35, 1]]}"#, 1),
    ] {
        let text = original.replace(COMMAND, command);
        assert_ne!(text, original);
        let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        scene_context(|cx, _| {
            for _ in 0..successful_steps {
                scene.advance(f64::MAX, cx).unwrap();
            }
            let before = scene.render_json();
            let checkpoint = scene.checkpoint().unwrap();
            let error = scene.advance(f64::MAX, cx).unwrap_err();
            assert!(error.message.contains("closed command (0 or 1)"), "{error}");
            assert_eq!(scene.steps(), successful_steps);
            assert_eq!(scene.render_json(), before);
            assert_eq!(scene.checkpoint().unwrap(), checkpoint);
        });
    }
}

#[test]
fn prescribed_device_current_keeps_its_inward_sign_and_refuses_an_impossible_open_circuit() {
    for (kind, current, closed, accepted) in [
        ("ideal_diode", 2.0, None, true),
        ("ideal_diode", -2.0, None, false),
        ("ideal_diode", 0.0, None, true),
        ("ideal_switch", 2.0, Some(1.0), true),
        ("ideal_switch", -2.0, Some(1.0), true),
        ("ideal_switch", 2.0, Some(0.0), false),
        ("ideal_switch", -2.0, Some(0.0), false),
        ("ideal_switch", 0.0, Some(0.0), true),
    ] {
        let command = closed.map_or_else(String::new, |closed| {
            format!(r#", {{"prescribe":"device.closed","value":{closed}}}"#)
        });
        let text = format!(
            r#"{{"schema":"frankensim.system.v0","duration_s":0.1,"step_s":0.1,
                "elements":{{"device":{{"kind":"{kind}"}}}},
                "connections":[{{"prescribe":"device.p","current_a":{current}}}{command}]}}"#
        );
        let mut scene = fs_cli::SystemScene::from_text(&text, "prescribed-device.json").unwrap();
        let before = scene.render_json();
        let checkpoint = scene.checkpoint().unwrap();
        scene_context(|cx, _| {
            if accepted {
                assert!(scene.advance(f64::MAX, cx).unwrap());
                close(probe(&scene, "device.current_a"), current);
                close(probe(&scene, "device.voltage_v"), 0.0);
                close(scene.ledger().energy_residual_j(), 0.0);
                if let Some(closed) = closed {
                    close(probe(&scene, "device.closed"), closed);
                }
            } else {
                let error = scene.advance(f64::MAX, cx).unwrap_err();
                assert!(error.message.contains("prescribed nonzero current"), "{error}");
                assert_eq!(scene.steps(), 0);
                assert_eq!(scene.render_json(), before);
                assert_eq!(scene.checkpoint().unwrap(), checkpoint);
            }
        });
    }
}

#[test]
fn switching_scene_refuses_invalid_schedules_fields_wiring_and_unclaimed_error_control() {
    let path = example("electrical_switching.json");
    let text = std::fs::read_to_string(&path).unwrap();
    for (changed, expected) in [
        (
            text.replace(COMMAND, r#"{"held_table": [[0, 1], [0, 0]]}"#),
            "times do not strictly increase",
        ),
        (text.replace(COMMAND, r#"{"held_table": []}"#), "no knots"),
        (
            text.replace(r#"{"kind": "ideal_diode"}"#, r#"{"kind": "ideal_diode", "resistance_ohm": 1}"#),
            "resistance_ohm",
        ),
        (
            text.replace(r#"{"kind": "ideal_switch"}"#, r#"{"kind": "ideal_switch", "closed": 1}"#),
            "closed",
        ),
    ] {
        assert_ne!(changed, text);
        let error = fs_cli::SystemScene::from_text(&changed, &path).err().unwrap();
        assert!(error.message.contains(expected), "expected {expected}: {error}");
    }
    let unwired = text.replace(&format!("{SIGNAL_CONNECTION},"), "");
    assert_ne!(unwired, text);
    let error = fs_cli::SystemScene::from_text(&unwired, &path).err().unwrap();
    assert!(error.message.contains("switch") && error.message.contains("closed"), "{error}");

    let adaptive = text.replacen(
        "\"elements\":",
        r#""adaptive": {"tolerances": {}}, "elements":"#,
        1,
    );
    let error = fs_cli::SystemScene::from_text(&adaptive, &path).err().unwrap();
    assert!(error.message.contains("ideal_diode") && error.message.contains("fixed step_s"), "{error}");
}
