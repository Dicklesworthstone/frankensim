//! A moving electric or magnetic gap exchanges energy with mechanics and heat.

use super::{example, run_binary, scene_context};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 5e-10,
        "{actual:e} vs {expected:e}"
    );
}

fn probe(scene: &fs_cli::SystemScene, path: &str) -> f64 {
    scene.probe(path).unwrap().unwrap()
}

fn coordinate(inductor: bool) -> &'static str {
    if inductor {
        "flux_linkage_wb"
    } else {
        "charge_c"
    }
}

fn held_scene(inductor: bool, adaptive: bool) -> String {
    let (kind, scale, input) = if inductor {
        ("gap_inductor", "mu_turns_area_h_m", "voltage_v")
    } else {
        ("gap_capacitor", "epsilon_area_f_m", "current_a")
    };
    let coordinate = coordinate(inductor);
    let adaptive = if adaptive {
        format!(
            r#""adaptive":{{"order":2,"tolerances":{{
          "field.{coordinate}":{{"absolute":1e-8,"relative":0}},
          "field.gap_m":{{"absolute":1e-8,"relative":0}},
          "body.position_m":{{"absolute":1e-8,"relative":0}},
          "body.velocity_m_s":{{"absolute":1e-8,"relative":0}},
          "push.energy_in_j":{{"absolute":1e-8,"relative":0}}
        }}}},"#
        )
    } else {
        String::new()
    };
    format!(
        r#"{{"schema":"frankensim.system.v0","duration_s":0.5,"step_s":0.0625,
      {adaptive}"elements":{{
        "field":{{"kind":"{kind}","{scale}":2.0,"gap_limits_m":[0.5,1.5],
          "gap_m":1.0,"{coordinate}":0.8}},
        "body":{{"kind":"mass","mass_kg":3.0,"position_m":1.0,"velocity_m_s":0.02}},
        "push":{{"kind":"force","force_n":0.1}}}},
      "connections":[{{"join":["body.p",{{"port":"field.gap","ratio":-1}},"push.p"]}},
        {{"prescribe":"field.field","{input}":0}}]}}"#
    )
}

#[test]
fn held_charge_and_flux_accelerate_the_mass_with_the_closed_form_force() {
    for inductor in [false, true] {
        for adaptive in [false, true] {
            let text = held_scene(inductor, adaptive);
            let mut scene = fs_cli::SystemScene::from_text(&text, "held-field.json").unwrap();
            for name in ["voltage_v", "current_a", "attraction_n", "gap_velocity_m_s"] {
                assert_eq!(scene.probe(&format!("field.{name}")).unwrap(), None);
            }
            scene_context(|cx, _| {
                scene.advance(f64::MAX, cx).unwrap();
                let (bytes, root) = scene.checkpoint().unwrap();
                let mut resumed = fs_cli::SystemScene::from_text(&text, "held-field.json").unwrap();
                resumed.resume_checkpoint(&bytes, root).unwrap();
                scene.run(cx).unwrap();
                resumed.run(cx).unwrap();
                assert_eq!(resumed.render_json(), scene.render_json());
            });
            // E = x^2 g / (2s), so a held x attracts with x^2 / (2s).
            let attraction = 0.8_f64.powi(2) / 4.0;
            let acceleration = (0.1 - attraction) / 3.0;
            let gap = 1.0 + 0.02 * 0.5 + 0.5 * acceleration * 0.5_f64.powi(2);
            let velocity = 0.02 + acceleration * 0.5;
            close(
                probe(&scene, &format!("field.{}", coordinate(inductor))),
                0.8,
            );
            close(probe(&scene, "field.gap_m"), gap);
            close(probe(&scene, "body.position_m"), gap);
            close(probe(&scene, "body.velocity_m_s"), velocity);
            close(probe(&scene, "field.attraction_n"), attraction);
            close(probe(&scene, "push.energy_in_j"), -0.1 * (gap - 1.0));
            close(
                attraction * (gap - 1.0) + 1.5 * (velocity.powi(2) - 0.02_f64.powi(2)),
                0.1 * (gap - 1.0),
            );
            close(scene.ledger().energy_residual_j(), 0.0);
            if adaptive {
                for name in [coordinate(inductor), "gap_m"] {
                    let missing = text.replace(
                        &format!(r#""field.{name}":{{"absolute":1e-8,"relative":0}},"#),
                        "",
                    );
                    assert_ne!(missing, text);
                    let error = fs_cli::SystemScene::from_text(&missing, "missing-tolerance.json")
                        .err()
                        .unwrap();
                    assert!(error.message.contains(&format!("field.{name}")), "{error}");
                }
            }
        }
    }
}

/// Eliminate motion from the two-coordinate energy gradient, then bisect the
/// remaining scalar electrical balance. No component or kernel solver is used.
fn reference_step([x0, g0, v0]: [f64; 3], loss: f64) -> ([f64; 3], [f64; 3], f64) {
    let (s, m, f, dt) = (1.5, 2.0, 0.05, 0.0625);
    let evaluate = |x1: f64| {
        let force = (x0 * x0 + x0 * x1 + x1 * x1) / (6.0 * s);
        let v1 = v0 + dt * (f - force) / m;
        let g1 = g0 + 0.5 * dt * (v0 + v1);
        let effort = ((2.0 * x0 + x1) * g0 + (x0 + 2.0 * x1) * g1) / (6.0 * s);
        (g1, v1, effort, force)
    };
    let (mut lower, mut upper) = (0.0, x0);
    for _ in 0..70 {
        let x1 = 0.5 * (lower + upper);
        if x1 - x0 + dt * loss * evaluate(x1).2 > 0.0 {
            upper = x1;
        } else {
            lower = x1;
        }
    }
    let x1 = 0.5 * (lower + upper);
    let (g1, v1, effort, force) = evaluate(x1);
    (
        [x1, g1, v1],
        [effort, force, 0.5 * (v0 + v1)],
        dt * loss * effort * effort,
    )
}

#[test]
fn authored_gaps_match_coupled_discharge_motion_heat_and_checkpoint_continuation() {
    for inductor in [false, true] {
        let path = example(if inductor {
            "gap_inductor.json"
        } else {
            "gap_capacitor.json"
        });
        let text = std::fs::read_to_string(&path).unwrap();
        let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        let loss = if inductor { 0.8 } else { 1.0 / 0.8 };
        let (mut state, mut heat) = ([1.2, 1.0, 0.0], 0.0);
        let mut saved = None;
        scene_context(|cx, _| {
            for step in 1..=16 {
                let (next, [effort, force, gap_velocity], joules) = reference_step(state, loss);
                state = next;
                heat += joules;
                assert!(scene.advance(f64::MAX, cx).unwrap());
                close(
                    probe(&scene, &format!("field.{}", coordinate(inductor))),
                    state[0],
                );
                close(probe(&scene, "field.gap_m"), state[1]);
                close(probe(&scene, "body.position_m"), state[1]);
                close(probe(&scene, "body.velocity_m_s"), state[2]);
                close(probe(&scene, "field.attraction_n"), force);
                close(probe(&scene, "field.gap_velocity_m_s"), gap_velocity);
                let (voltage, current) = if inductor {
                    (-loss * effort, effort)
                } else {
                    (effort, -loss * effort)
                };
                close(probe(&scene, "field.voltage_v"), voltage);
                close(probe(&scene, "field.current_a"), current);
                close(probe(&scene, "case.temperature_k"), 300.0 + heat / 1000.0);
                close(probe(&scene, "push.energy_in_j"), -0.05 * (state[1] - 1.0));
                // Check actual stored energy, independently of the ledger's receipts.
                let actual_x = probe(&scene, &format!("field.{}", coordinate(inductor)));
                let actual_g = probe(&scene, "field.gap_m");
                let actual_v = probe(&scene, "body.velocity_m_s");
                let actual_heat = 1000.0 * (probe(&scene, "case.temperature_k") - 300.0);
                close(
                    actual_x * actual_x * actual_g / 3.0 + actual_v * actual_v + actual_heat
                        - 1.2_f64.powi(2) / 3.0,
                    0.05 * (actual_g - 1.0),
                );
                close(scene.ledger().energy_residual_j(), 0.0);
                if step == 5 {
                    saved = Some(scene.checkpoint().unwrap());
                }
            }
            assert!(scene.is_finished());
            assert!(state[0] < 1.0 && state[1] < 0.95 && heat > 0.0);
            assert!(scene.ledger().entropy_generated_j_per_k() > 0.0);
            let (bytes, root) = saved.unwrap();
            let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
            resumed.resume_checkpoint(&bytes, root).unwrap();
            resumed.run(cx).unwrap();
            assert_eq!(resumed.render_json(), scene.render_json());
            for changed in [
                text.replacen(": 1.5,", ": 2.0,", 1),
                text.replace("[0.5, 1.5]", "[0.4, 1.5]"),
            ] {
                assert_ne!(changed, text);
                let mut foreign = fs_cli::SystemScene::from_text(&changed, &path).unwrap();
                let before = foreign.render_json();
                assert!(foreign.resume_checkpoint(&bytes, root).is_err());
                assert_eq!(foreign.render_json(), before);
            }
        });
        let (code, json, stderr) = run_binary(&["--json", "system", &path]);
        assert_eq!(code, 0, "{stderr}");
        assert_eq!(json, scene.render_json());
    }
}

#[test]
fn gap_support_is_open_and_a_refused_step_preserves_the_accepted_prefix() {
    for inductor in [false, true] {
        let text = held_scene(inductor, false);
        let (kind, scale, input) = if inductor {
            ("gap_inductor", "mu_turns_area_h_m", "voltage_v")
        } else {
            ("gap_capacitor", "epsilon_area_f_m", "current_a")
        };
        for changed in [
            text.replace("\"gap_m\":1.0", "\"gap_m\":0.5"),
            text.replace("\"gap_m\":1.0", "\"gap_m\":1.5"),
            text.replace("[0.5,1.5]", "[1.5,0.5]"),
            text.replace(
                &format!("\"{}\":0.8", coordinate(inductor)),
                &format!("\"{}\":1e140", coordinate(inductor)),
            )
            .replace(&format!("\"{scale}\":2.0"), &format!("\"{scale}\":1e-140")),
        ] {
            assert_ne!(changed, text);
            assert!(fs_cli::SystemScene::from_text(&changed, "bad-gap.json").is_err());
        }
        let text = format!(
            r#"{{"schema":"frankensim.system.v0","duration_s":0.25,"step_s":0.125,
          "elements":{{"field":{{"kind":"{kind}","{scale}":2,"gap_limits_m":[0.5,1.5],"gap_m":1}}}},
          "connections":[{{"prescribe":"field.field","{input}":0}},
            {{"prescribe":"field.gap","velocity_m_s":{{"held_table":[[0,0],[0.125,4]]}}}}]}}"#
        );
        let mut scene = fs_cli::SystemScene::from_text(&text, "gap-limit.json").unwrap();
        scene_context(|cx, _| {
            assert!(scene.advance(f64::MAX, cx).unwrap());
            let before = scene.render_json();
            let checkpoint = scene.checkpoint().unwrap();
            assert!(scene.advance(f64::MAX, cx).is_err());
            assert_eq!(scene.steps(), 1);
            assert_eq!(scene.render_json(), before);
            assert_eq!(scene.checkpoint().unwrap(), checkpoint);
        });
    }
}
