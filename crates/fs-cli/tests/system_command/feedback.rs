//! Authored feedback scenes use real plant, controller and actuator Components.

use super::{element, example, run_binary, scene_context, value};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * (1.0 + expected.abs()),
        "{actual} != {expected}"
    );
}

fn probe(scene: &fs_cli::SystemScene, path: &str) -> f64 {
    scene.probe(path).unwrap().unwrap()
}

#[test]
fn affine_and_sampled_feedback_follow_independent_mass_recurrences_and_resume() {
    for (name, sampled) in [
        ("feedback_position.json", false),
        ("sampled_pi_velocity.json", true),
    ] {
        let path = example(name);
        let authored = std::fs::read_to_string(&path).unwrap();
        // A second physics step does not divide the controller's sample period.
        // The kernel must still sample once per period and hold between events.
        for split_samples in [false, true] {
            if split_samples && !sampled {
                continue;
            }
            let text = if split_samples {
                authored.replace("\"step_s\": 0.03125", "\"step_s\": 0.046875")
            } else {
                authored.clone()
            };
            let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
            assert_eq!(scene.probe("controller.output0").unwrap(), None);
            assert_eq!(scene.probe("drive.command").unwrap(), None);
            let (mut position, mut velocity) = if sampled { (0.0, 0.0) } else { (0.25, -0.5) };
            let initial_energy = velocity * velocity; // m = 2 kg.
            let (mut integral, mut command, mut samples) = (0.0, 0.0, 0_u32);
            let mut saved = None;
            scene_context(|cx, _| {
                while !scene.is_finished() {
                    let start = scene.time_s();
                    let reference = if sampled {
                        1.0
                    } else if start < 0.5 {
                        1.0 - 3.0 * start
                    } else {
                        -0.5
                    };
                    if sampled {
                        if start >= f64::from(samples) * 0.125 {
                            let error = reference - velocity;
                            command = 2.0 * integral + 3.0 * error;
                            integral += 0.125 * error;
                            samples += 1;
                        }
                    } else {
                        command = 4.0 * (reference - position) - 3.0 * velocity;
                    }
                    assert!(scene.advance(f64::MAX, cx).unwrap());
                    let dt = scene.time_s() - start;
                    let next_velocity = velocity + dt * command / 2.0;
                    position += 0.5 * dt * (velocity + next_velocity);
                    velocity = next_velocity;
                    close(probe(&scene, "plant.position_m"), position);
                    close(probe(&scene, "plant.velocity_m_s"), velocity);
                    close(probe(&scene, "controller.output0"), command);
                    close(probe(&scene, "drive.command"), command);
                    close(probe(&scene, "reference.output0"), reference);
                    close(probe(&scene, "reference.time_s"), scene.time_s());
                    close(
                        probe(&scene, "drive.energy_in_j"),
                        initial_energy - velocity * velocity,
                    );
                    if sampled {
                        close(probe(&scene, "controller.x0"), integral);
                        assert_eq!(probe(&scene, "controller.sample_count"), f64::from(samples));
                        close(probe(&scene, "controller.time_s"), scene.time_s());
                    }
                    assert!(scene.ledger().energy_residual_j().abs() < 2e-9);
                    if saved.is_none() && scene.time_s() >= 0.5 {
                        saved = Some(scene.checkpoint().unwrap());
                    }
                }
                let (bytes, root) = saved.unwrap();
                let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
                resumed.resume_checkpoint(&bytes, root).unwrap();
                resumed.run(cx).unwrap();
                assert_eq!(resumed.render_json(), scene.render_json());
            });
            if sampled {
                assert_eq!(samples, 16);
            }
            if !split_samples {
                let (code, json, stderr) = run_binary(&["--json", "system", &path]);
                assert_eq!(code, 0, "{stderr}");
                assert_eq!(json, scene.render_json());
            }
        }
    }
}

#[test]
fn controlled_current_charges_with_the_declared_sign_and_point_sampled_commands() {
    let path = example("signal_current_charge.json");
    let authored = std::fs::read_to_string(&path).unwrap();
    for prescribed in [false, true] {
        let text = if prescribed {
            authored.replace(
                r#"{"signal": "reference.output", "to": ["drive.command"]}"#,
                r#"{"prescribe": "drive.command", "value": {"table": [[0, -2], [0.5, -1]]}}"#,
            )
        } else {
            authored.clone()
        };
        let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        let mut charge = 0.25;
        scene_context(|cx, _| {
            while !scene.is_finished() {
                let start = scene.time_s();
                let command = if prescribed { -2.0 + 2.0 * start } else { -2.0 };
                assert!(scene.advance(f64::MAX, cx).unwrap());
                charge -= command * (scene.time_s() - start);
                close(probe(&scene, "capacitor.charge_c"), charge);
                close(probe(&scene, "capacitor.voltage_v"), charge / 0.5);
                close(probe(&scene, "drive.command"), command);
                // Energy = q^2 / (2 C), with C = 0.5 F. Source work is signed
                // into the source, so charging the capacitor makes it negative.
                close(
                    probe(&scene, "drive.energy_in_j"),
                    0.25 * 0.25 - charge * charge,
                );
                assert!(scene.ledger().energy_residual_j().abs() < 2e-9);
            }
        });
        if !prescribed {
            let (code, json, stderr) = run_binary(&["--json", "system", &path]);
            assert_eq!(code, 0, "{stderr}");
            assert_eq!(json, scene.render_json());
            close(value(&json, &element("capacitor"), "charge_c"), 1.25);
            close(value(&json, &element("drive"), "energy_in_j"), -1.5);
        }
    }
}

#[test]
fn mapped_controller_outputs_drive_two_axes_with_independent_work_accounts() {
    let path = example("sampled_two_axis.json");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    let mut position = [0.2, -0.1];
    let mut velocity = [-0.25, 0.5];
    let mass = [2.0, 4.0];
    let initial_energy = [0.0625, 0.5];
    let mut integral = [0.1, -0.2];
    let mut output = [0.0; 2];
    let mut samples = 0_u32;
    let mut saved = None;
    scene_context(|cx, _| {
        while !scene.is_finished() {
            let start = scene.time_s();
            if start >= f64::from(samples) * 0.125 {
                let error = [1.0 - velocity[0], -0.5 - velocity[1]];
                // Independent sampled PI laws in controller coordinates,
                // followed by the authored two-axis actuator allocation.
                output = [
                    2.0 * integral[0] + 3.0 * error[0],
                    3.0 * integral[1] + 4.0 * error[1],
                ];
                for axis in 0..2 {
                    integral[axis] += 0.125 * error[axis];
                }
                samples += 1;
            }
            let forces = [output[0] + 0.5 * output[1], -0.25 * output[0] + output[1]];
            assert!(scene.advance(f64::MAX, cx).unwrap());
            let dt = scene.time_s() - start;
            for (axis, name) in ["x", "y"].into_iter().enumerate() {
                let next = velocity[axis] + dt * forces[axis] / mass[axis];
                position[axis] += 0.5 * dt * (velocity[axis] + next);
                velocity[axis] = next;
                close(
                    probe(&scene, &format!("plant_{name}.position_m")),
                    position[axis],
                );
                close(probe(&scene, &format!("plant_{name}.velocity_m_s")), next);
                close(
                    probe(&scene, &format!("drive_{name}.command")),
                    forces[axis],
                );
                close(
                    probe(&scene, &format!("drive_{name}.energy_in_j")),
                    initial_energy[axis] - 0.5 * mass[axis] * next * next,
                );
                close(
                    probe(&scene, &format!("controller.output{axis}")),
                    output[axis],
                );
                close(
                    probe(&scene, &format!("controller.x{axis}")),
                    integral[axis],
                );
            }
            assert_eq!(probe(&scene, "controller.sample_count"), f64::from(samples));
            assert!(scene.ledger().energy_residual_j().abs() < 2e-9);
            if saved.is_none() && scene.time_s() >= 0.5 {
                saved = Some(scene.checkpoint().unwrap());
            }
        }
        assert_eq!(samples, 16);
        let (bytes, root) = saved.unwrap();
        let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        resumed.resume_checkpoint(&bytes, root).unwrap();
        resumed.run(cx).unwrap();
        assert_eq!(resumed.render_json(), scene.render_json());

        let changed_map = text.replace("[[1, 0.5]]", "[[1, 0.25]]");
        let mut changed = fs_cli::SystemScene::from_text(&changed_map, &path).unwrap();
        assert!(changed.resume_checkpoint(&bytes, root).is_err());
        assert_eq!(changed.time_s(), 0.0);
    });
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(json, scene.render_json());
}

#[test]
fn signal_projection_admission_preserves_copy_connections_and_rejects_invalid_maps() {
    let path = example("sampled_pi_velocity.json");
    let original = std::fs::read_to_string(&path).unwrap();
    let projected = original.replace(
        r#"{"signal": "controller.output", "to": ["drive.command"]}"#,
        r#"{"signal": "controller.output", "to": ["drive.command"], "matrix": [[1]]}"#,
    );
    assert_ne!(projected, original);
    let mut copied = fs_cli::SystemScene::from_text(&original, &path).unwrap();
    let mut mapped = fs_cli::SystemScene::from_text(&projected, &path).unwrap();
    scene_context(|cx, _| {
        copied.run(cx).unwrap();
        mapped.run(cx).unwrap();
    });
    assert_eq!(mapped.render_json(), copied.render_json());

    let text = std::fs::read_to_string(example("sampled_two_axis.json")).unwrap();
    for (changed, expected) in [
        (
            text.replace("[[1, 0.5]]", "[[1, 0.5], [0, 1]]"),
            "row count",
        ),
        (text.replace("[[1, 0.5]]", "[[1]]"), "column count"),
        (text.replace("[[1, 0.5]]", "[[1, null]]"), "finite number"),
        (
            text.replace("\"drive_x.command\"", "\"plant_x.p\""),
            "signal input",
        ),
    ] {
        let error = fs_cli::SystemScene::from_text(&changed, "invalid-signal-map.json")
            .err()
            .unwrap();
        assert!(
            error.message.contains(expected),
            "expected {expected}: {error}"
        );
    }
    // Vector-to-scalar wiring must remain explicit, and the same command
    // input may never have two writers even when both use the same map.
    for changed in [
        text.replace(", \"matrix\": [[1, 0.5]]", ""),
        text.replace(
            "[\"drive_x.command\"]",
            "[\"drive_x.command\", \"drive_x.command\"]",
        ),
    ] {
        assert!(fs_cli::SystemScene::from_text(&changed, "invalid-signal-wiring.json").is_err());
    }
}

fn adaptive(text: &str, order: u32, source_work: bool) -> String {
    let work = if source_work {
        r#", "drive.energy_in_j": {"absolute": 0.001, "relative": 0}"#
    } else {
        ""
    };
    text.replacen(
        "\"elements\":",
        &format!(
            r#""adaptive": {{"order": {order}, "tolerances": {{
          "plant.position_m": {{"absolute": 0.001, "relative": 0}},
          "plant.velocity_m_s": {{"absolute": 0.001, "relative": 0}}{work}
        }}}}, "elements":"#
        ),
        1,
    )
}

#[test]
fn control_scene_validation_preserves_shapes_observation_opt_in_and_error_accounts() {
    let affine = std::fs::read_to_string(example("feedback_position.json")).unwrap();
    let sampled = std::fs::read_to_string(example("sampled_pi_velocity.json")).unwrap();
    for (text, expected) in [
        (
            affine.replace("\"publish\": true", "\"publish\": false"),
            "plant.position",
        ),
        (
            affine.replace("\"publish\": true", "\"publish\": 1"),
            "publish",
        ),
        (
            affine.replace(
                "\"publish\": true",
                "\"publish\": true, \"publish_state\": false",
            ),
            "publish",
        ),
        (
            affine.replace("[[4, -4, -3]]", "[[4, -4, -3], [1, 2]]"),
            "column count",
        ),
        (
            affine.replace("\"bias\": [0]", "\"bias\": [0, 1]"),
            "expected 1 coefficients",
        ),
        (affine.replace("\"translation\"", "\"thermal\""), "domain"),
        (
            sampled.replace("\"waveforms\": [1]", "\"waveforms\": []"),
            "waveforms",
        ),
        (
            sampled.replace("\"sample_period_s\": 0.125", "\"sample_period_s\": 0"),
            "sample_period_s",
        ),
        (
            sampled.replace("\"c\": [[2]]", "\"c\": [[2, 3]]"),
            "column count",
        ),
        (adaptive(&affine, 1, false), "drive.energy_in_j"),
        (adaptive(&affine, 2, true), "supported order is 1"),
        (adaptive(&sampled, 1, true), "fixed step_s"),
    ] {
        let error = fs_cli::SystemScene::from_text(&text, "control-invalid.json")
            .err()
            .unwrap();
        assert!(
            error.message.contains(expected),
            "expected {expected:?}: {error}"
        );
    }
    let mut scene =
        fs_cli::SystemScene::from_text(&adaptive(&affine, 1, true), "adaptive-control.json")
            .unwrap();
    scene_context(|cx, _| scene.run(cx).unwrap());
    let velocity = probe(&scene, "plant.velocity_m_s");
    close(
        probe(&scene, "drive.energy_in_j"),
        0.25 - velocity * velocity,
    );
    assert!(scene.ledger().energy_residual_j().abs() < 2e-9);

    // A generic point command needs the first-order guard even without an
    // affine map, SignalSource or published plant measurement in the graph.
    // SI-named mean commands keep their separate higher-order coverage in
    // the parent module's controlled-source tests.
    let direct = r#"{
      "schema": "frankensim.system.v0", "duration_s": 0.25, "step_s": 0.125,
      "elements": {
        "plant": {"kind": "mass", "mass_kg": 2},
        "drive": {"kind": "controlled_effort", "domain": "translation"}
      },
      "connections": [
        {"join": ["plant.p", "drive.p"]},
        {"prescribe": "drive.command", "value": -2}
      ]
    }"#;
    let error =
        fs_cli::SystemScene::from_text(&adaptive(direct, 2, true), "point-command-order.json")
            .err()
            .unwrap();
    assert!(error.message.contains("supported order is 1"), "{error}");
    let mut direct_scene =
        fs_cli::SystemScene::from_text(&adaptive(direct, 1, true), "point-command.json").unwrap();
    scene_context(|cx, _| direct_scene.run(cx).unwrap());
    close(probe(&direct_scene, "plant.velocity_m_s"), -0.25);
    close(probe(&direct_scene, "drive.energy_in_j"), -0.0625);

    // Existing accepted-state observations can also command an actuator;
    // requiring order one cannot depend on using a new signal element.
    let thermal_feedback = r#"{
      "schema": "frankensim.system.v0", "duration_s": 0.25, "step_s": 0.125,
      "adaptive": {"order": 2, "tolerances": {
        "plant.position_m": {"absolute": 0.001, "relative": 0},
        "plant.velocity_m_s": {"absolute": 0.001, "relative": 0},
        "sensor.temperature_k": {"absolute": 0.001, "relative": 0},
        "drive.energy_in_j": {"absolute": 0.001, "relative": 0}
      }},
      "elements": {
        "plant": {"kind": "mass", "mass_kg": 2},
        "sensor": {"kind": "heat_capacity", "capacity_j_k": 1, "temperature_k": 300},
        "drive": {"kind": "controlled_effort", "domain": "translation"}
      },
      "connections": [
        {"join": ["plant.p", "drive.p"]},
        {"prescribe": "sensor.heat", "heat_rate_w": 1},
        {"signal": "sensor.temperature", "to": ["drive.command"]}
      ]
    }"#;
    let error = fs_cli::SystemScene::from_text(thermal_feedback, "signal-feedback-order.json")
        .err()
        .unwrap();
    assert!(error.message.contains("supported order is 1"), "{error}");
}
