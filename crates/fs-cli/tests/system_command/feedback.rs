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

#[test]
fn controlled_hydraulics_close_volume_work_and_heat_and_resume() {
    let path = example("controlled_hydraulics.json");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    assert_eq!(scene.probe("pump.flow_m3_s").unwrap(), None);
    assert_eq!(scene.probe("reservoir.pressure_pa").unwrap(), None);
    let (capacity, resistance, initial_volume) = (1e-9, 1e9, 2e-4);
    let mut volume = initial_volume;
    let (mut pump_work, mut reservoir_work, mut heat) = (0.0, 0.0, 0.0);
    let mut flow_signs = [false; 2];
    let mut events = [false; 3];
    let mut saved = None;
    scene_context(|cx, _| {
        while !scene.is_finished() {
            let start = scene.time_s();
            let command = if start < 0.6 {
                -1e-4
            } else if start < 1.2 {
                1e-4
            } else {
                -2e-4
            };
            let reservoir_pressure = if start < 0.9 { 1e5 } else { 5e4 };
            assert!(scene.advance(f64::MAX, cx).unwrap());
            let dt = scene.time_s() - start;
            // C dp/dt = -q_pump - (p - p_reservoir) / R.
            // Integrate the linear balance independently at the midpoint.
            let alpha = dt / (2.0 * resistance * capacity);
            let next_volume = ((1.0 - alpha) * volume
                + dt * (reservoir_pressure / resistance - command))
                / (1.0 + alpha);
            let pressure = (volume + next_volume) / (2.0 * capacity);
            let flow = (pressure - reservoir_pressure) / resistance;
            let pump_increment = dt * pressure * command;
            pump_work += pump_increment;
            reservoir_work += dt * reservoir_pressure * flow;
            heat += dt * resistance * flow * flow;
            volume = next_volume;
            flow_signs[if command > 0.0 { 1 } else { 0 }] = true;
            assert!(pump_increment * command > 0.0);
            close(probe(&scene, "store.volume_m3"), volume);
            close(probe(&scene, "store.pressure_pa"), volume / capacity);
            close(probe(&scene, "pump.flow_m3_s"), command);
            close(probe(&scene, "pump.command"), command);
            close(probe(&scene, "reservoir.pressure_pa"), reservoir_pressure);
            close(probe(&scene, "line.flow_m3_s"), flow);
            close(probe(&scene, "line.heat_w"), resistance * flow * flow);
            close(probe(&scene, "pump.energy_in_j"), pump_work);
            close(probe(&scene, "reservoir.energy_in_j"), reservoir_work);
            close(probe(&scene, "oil.temperature_k"), 300.0 + heat / 50.0);
            close(
                (volume * volume - initial_volume * initial_volume) / (2.0 * capacity)
                    + pump_work
                    + reservoir_work
                    + heat,
                0.0,
            );
            assert!(scene.ledger().energy_residual_j().abs() < 2e-8);
            for (seen, event) in events.iter_mut().zip([0.6, 0.9, 1.2]) {
                assert!(!(start < event && scene.time_s() > event + 1e-12));
                *seen |= (scene.time_s() - event).abs() < 1e-12;
            }
            if saved.is_none() && scene.time_s() >= 0.9 {
                saved = Some(scene.checkpoint().unwrap());
            }
        }
        let (bytes, root) = saved.unwrap();
        let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        resumed.resume_checkpoint(&bytes, root).unwrap();
        resumed.run(cx).unwrap();
        assert_eq!(resumed.render_json(), scene.render_json());
        let mut changed =
            fs_cli::SystemScene::from_text(&text.replace("[0.9, 5e4]", "[0.9, 6e4]"), &path)
                .unwrap();
        assert!(changed.resume_checkpoint(&bytes, root).is_err());
        assert_eq!(changed.time_s(), 0.0);
    });
    assert!(flow_signs.into_iter().all(|seen| seen));
    assert!(events.into_iter().all(|seen| seen));
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(json, scene.render_json());
}

const HYDRAULIC_LAWS: &str = r#"{
  "schema": "frankensim.system.v0", "duration_s": 0.25, "step_s": 0.25,
  "elements": {
    "gyrator": {"kind": "gyrator", "ratio": 3, "a": "hydraulic", "b": "rotation"},
    "transformer": {"kind": "modulated_transformer", "a": "rotation", "b": "hydraulic"},
    "loss": {"kind": "resistance_matrix", "domain": "hydraulic", "matrix": [[3, 1], [1, 2]]},
    "storage": {"kind": "port_hamiltonian",
      "energy": [[4]], "structure": [[0]], "dissipation": [[0]], "input": [[1]],
      "state": [0.125], "ports": [{"name": "p", "domain": "hydraulic", "role": "flow_in"}]},
    "case": {"kind": "heat_capacity", "capacity_j_k": 10, "temperature_k": 300}
  },
  "connections": [
    {"prescribe": "gyrator.a", "flow_m3_s": 2},
    {"prescribe": "gyrator.b", "angular_velocity_rad_s": 4},
    {"prescribe": "transformer.a", "angular_velocity_rad_s": 2},
    {"prescribe": "transformer.b", "pressure_pa": 5},
    {"prescribe": "transformer.ratio", "ratio": -3},
    {"prescribe": "loss.u", "flow_m3_s": [2, -1]},
    {"prescribe": "storage.p", "flow_m3_s": -0.5},
    {"join": ["case.heat", "loss.heat", "storage.heat"]}
  ]
}"#;

#[test]
fn hydraulic_converter_matrix_and_phs_ports_keep_native_power_equations() {
    let mut scene =
        fs_cli::SystemScene::from_text(HYDRAULIC_LAWS, "hydraulic-laws.json").unwrap();
    scene_context(|cx, _| scene.run(cx).unwrap());
    // Prescribed attachment flows enter these generic ports with the
    // opposite sign, just as in the other admitted power domains.
    let a = probe(&scene, "gyrator.effort_a");
    let b = probe(&scene, "gyrator.effort_b");
    close(a, -12.0); // e_a = 3 f_b, f_b = -4 rad/s.
    close(b, 6.0); // e_b = -3 f_a, f_a = -2 m^3/s.
    close(-2.0 * a - 4.0 * b, 0.0);
    let effort = probe(&scene, "transformer.effort_a");
    let flow = probe(&scene, "transformer.flow_b");
    close(effort, -15.0);
    close(flow, -6.0);
    close(-2.0 * effort + 5.0 * flow, 0.0);
    // R [-2, 1]' = [-5, 0]'; u' R u = 10 W, all delivered as heat.
    close(probe(&scene, "loss.u0.pressure_pa"), -5.0);
    close(probe(&scene, "loss.u1.pressure_pa"), 0.0);
    close(probe(&scene, "loss.heat_w"), 10.0);
    close(probe(&scene, "case.temperature_k"), 300.25);
    // dx/dt = 0.5 m^3/s, H = 2 x^2 J: final x = 0.25 m^3.
    close(probe(&scene, "storage.x0"), 0.25);
    close(probe(&scene, "storage.p.pressure_pa"), 0.75);
    close(probe(&scene, "storage.heat_w"), 0.0);
    close(
        value(&scene.render_json(), &element("storage"), "stored_change_j"),
        0.09375,
    );
    assert!(scene.ledger().energy_residual_j().abs() < 2e-9);
}

#[test]
fn hydraulic_scene_commands_and_domains_refuse_wrong_units_and_unknown_names() {
    let text = std::fs::read_to_string(example("controlled_hydraulics.json")).unwrap();
    for (from, to, expected) in [
        (
            r#""flow_m3_s": {"held_table""#,
            r#""current_a": {"held_table""#,
            "flow_m3_s",
        ),
        (
            r#""pressure_pa": {"held_table""#,
            r#""force_n": {"held_table""#,
            "pressure_pa",
        ),
        (
            r#""domain": "hydraulic""#,
            r#""domain": "steam""#,
            "electrical, hydraulic",
        ),
    ] {
        let changed = text.replace(from, to);
        assert_ne!(changed, text);
        let error = fs_cli::SystemScene::from_text(&changed, "invalid-hydraulic-command.json")
            .err()
            .unwrap();
        assert!(error.message.contains(expected), "{error}");
    }
    // Every shared Domain::parse consumer lists the newly admitted domain
    // while continuing to reject thermal/unknown power-domain declarations.
    for (from, to) in [
        (r#""a": "hydraulic""#, r#""a": "thermal""#),
        (r#""b": "hydraulic""#, r#""b": "thermal""#),
        (
            r#""domain": "hydraulic", "matrix""#,
            r#""domain": "thermal", "matrix""#,
        ),
        (
            r#""domain": "hydraulic", "role""#,
            r#""domain": "thermal", "role""#,
        ),
    ] {
        let changed = HYDRAULIC_LAWS.replace(from, to);
        assert_ne!(changed, HYDRAULIC_LAWS);
        let error = fs_cli::SystemScene::from_text(&changed, "invalid-hydraulic-domain.json")
            .err()
            .unwrap();
        assert!(error.message.contains("electrical, hydraulic"), "{error}");
    }
}

const SIGNAL_LIMITS: &str = r#"{
  "schema": "frankensim.system.v0", "duration_s": 1, "step_s": 0.19,
  "sample_every_s": 1,
  "elements": {
    "reference": {"kind": "signal_source", "waveforms": [
      {"held_table": [[0, -1], [0.25, 0.5], [0.5, 2], [0.75, 0]]},
      {"table": [[0, -4], [1, 4]]}, 9
    ]},
    "limits": {"kind": "signal_limiter", "lower": [0, -2, 3], "upper": [1, 2, 3]}
  },
  "connections": [{"signal": "reference.output", "to": ["limits.input"]}]
}"#;

#[test]
fn signal_limiter_scenes_bound_vector_inputs_and_replay_prescribed_events() {
    for prescribed in [false, true] {
        let text = if prescribed {
            SIGNAL_LIMITS.replace(
                r#"{"signal": "reference.output", "to": ["limits.input"]}"#,
                r#"{"prescribe": "limits.input", "value": [
                    {"held_table": [[0, -1], [0.25, 0.5], [0.5, 2], [0.75, 0]]},
                    {"table": [[0, -4], [1, 4]]}, 9
                ]}"#,
            )
        } else {
            SIGNAL_LIMITS.to_owned()
        };
        let path = "signal-limits.json";
        let mut scene = fs_cli::SystemScene::from_text(&text, path).unwrap();
        for index in 0..3 {
            assert_eq!(scene.probe(&format!("limits.output{index}")).unwrap(), None);
        }
        let mut events = [false; 3];
        let mut saved = None;
        scene_context(|cx, _| {
            while !scene.is_finished() {
                let start = scene.time_s();
                let first = if start < 0.25 || start >= 0.75 {
                    0.0
                } else if start < 0.5 {
                    0.5
                } else {
                    1.0
                };
                // A varying prescribed signal is sampled at the left end,
                // exactly as a connected SignalSource is; it is not averaged.
                let second = (-4.0 + 8.0 * start).clamp(-2.0, 2.0);
                assert!(scene.advance(f64::MAX, cx).unwrap());
                close(probe(&scene, "limits.output0"), first);
                close(probe(&scene, "limits.output1"), second);
                close(probe(&scene, "limits.output2"), 3.0);
                for (seen, event) in events.iter_mut().zip([0.25, 0.5, 0.75]) {
                    assert!(!(start < event && scene.time_s() > event + 1e-13));
                    *seen |= (scene.time_s() - event).abs() < 1e-13;
                }
                if saved.is_none() && scene.time_s() >= 0.5 {
                    saved = Some(scene.checkpoint().unwrap());
                }
            }
            let (bytes, root) = saved.unwrap();
            let mut resumed = fs_cli::SystemScene::from_text(&text, path).unwrap();
            resumed.resume_checkpoint(&bytes, root).unwrap();
            resumed.run(cx).unwrap();
            assert_eq!(resumed.render_json(), scene.render_json());
            assert_eq!(resumed.checkpoint().unwrap(), scene.checkpoint().unwrap());
            for (from, to) in [
                (r#""lower": [0, -2, 3]"#, r#""lower": [0.1, -2, 3]"#),
                (r#""upper": [1, 2, 3]"#, r#""upper": [0.9, 2, 3]"#),
            ] {
                let changed_text = text.replace(from, to);
                assert_ne!(changed_text, text);
                let mut changed = fs_cli::SystemScene::from_text(&changed_text, path).unwrap();
                let before = changed.render_json();
                assert!(changed.resume_checkpoint(&bytes, root).is_err());
                assert_eq!(changed.time_s(), 0.0);
                assert_eq!(changed.render_json(), before);
            }
        });
        assert!(events.into_iter().all(|seen| seen), "{events:?}");
        assert_eq!(scene.ledger().energy_residual_j(), 0.0);
    }
}

#[test]
fn signal_limiter_scene_admission_rejects_bad_bounds_and_input_shapes() {
    for (from, to, expected) in [
        (r#""lower": [0, -2, 3]"#, r#""lower": []"#, "lower"),
        (r#""upper": [1, 2, 3]"#, r#""upper": [1, 2]"#, "upper"),
        (
            r#""lower": [0, -2, 3]"#,
            r#""lower": [2, -2, 3]"#,
            "lower bound",
        ),
        (
            r#""upper": [1, 2, 3]"#,
            r#""upper": [1, 2, 1e400]"#,
            "finite",
        ),
        (
            r#"{"signal": "reference.output", "to": ["limits.input"]}"#,
            r#"{"prescribe": "limits.input", "value": [0, 1]}"#,
            "3 waveforms",
        ),
        (
            r#"{"signal": "reference.output", "to": ["limits.input"]}"#,
            r#"{"prescribe": "limits.input", "temperature_k": [300, 300, 300]}"#,
            "value",
        ),
    ] {
        let text = SIGNAL_LIMITS.replace(from, to);
        assert_ne!(text, SIGNAL_LIMITS);
        let error = fs_cli::SystemScene::from_text(&text, "invalid-signal-limits.json")
            .err()
            .unwrap();
        assert!(
            error.message.contains(expected),
            "expected {expected}: {error}"
        );
    }
}
