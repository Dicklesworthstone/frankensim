//! Live passive dissipation reaches ordinary scenes and checkpoint replay.

use super::{example, run_binary, scene_context};

const GAS_CONSTANT: f64 = 8.314_462_618_153_24;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5.0e-10 * (1.0 + expected.abs()),
        "{actual:e} against {expected:e}"
    );
}

fn probe(scene: &fs_cli::SystemScene, name: &str) -> f64 {
    scene.probe(name).unwrap().unwrap()
}

#[test]
fn authored_dissipation_runs_in_each_scalar_domain_at_the_live_temperature() {
    for (kind, coefficient, flow, effort) in [
        ("damper", "damping_n_s_m", "velocity_m_s", "force_n"),
        (
            "rotary_damper",
            "damping_n_m_s_rad",
            "angular_velocity_rad_s",
            "torque_n_m",
        ),
        ("resistor", "resistance_ohm", "current_a", "voltage_v"),
    ] {
        for wlf in [false, true] {
            let shift = if wlf {
                r#"{"wlf":{"c1":1.0,"c2_k":100.0,"reference_k":300.0,"valid_k":[290.0,400.0]}}"#
            } else {
                r#"{"arrhenius":{"activation_energy_j_mol":24000.0,"reference_k":300.0,"valid_k":[290.0,400.0]}}"#
            };
            for (rate, temperature) in [(-2.0_f64, 330.0_f64), (0.0, 300.0), (2.0, 300.0)] {
                let text = format!(
                    r#"{{"schema":"frankensim.system.v0","duration_s":0.125,"step_s":0.125,
                      "elements":{{"loss":{{"kind":"{kind}","{coefficient}":3.0,
                        "dissipation_law":{{"quadratic_flow_scale":2.0,"shift":{shift}}}}}}},
                      "connections":[{{"prescribe":"loss.a","{flow}":{rate}}},
                        {{"prescribe":"loss.b","{flow}":0.0}},
                        {{"prescribe":"loss.heat","temperature_k":{temperature}}}]}}"#
                );
                let mut scene = fs_cli::SystemScene::from_text(&text, "dissipation.json").unwrap();
                scene_context(|cx, _| scene.run(cx).unwrap());
                let delta = temperature - 300.0;
                // Independent printed laws, without the material shift implementation.
                let thermal = if wlf {
                    10.0_f64.powf(-delta / (100.0 + delta))
                } else {
                    ((24000.0 / GAS_CONSTANT) * (temperature.recip() - 300.0_f64.recip())).exp()
                };
                let expected = 3.0 * thermal * (1.0 + rate.abs() / 2.0) * rate;
                close(probe(&scene, &format!("loss.{effort}")), expected);
                close(probe(&scene, "loss.heat_w"), expected * rate);
                close(scene.ledger().energy_residual_j(), 0.0);
                assert!(scene.ledger().entropy_generated_j_per_k() >= 0.0);
                assert!(scene.render_json().contains("live heat-port temperature"));
            }
        }
    }
}

/// One coupled midpoint balance, eliminating temperature by total energy.
/// This scalar bisection is independent of the kernel's multiport Newton solve.
fn reference_step(velocity: f64, temperature: f64) -> (f64, f64) {
    let next_temperature = |next: f64| temperature + 0.5 * (velocity * velocity - next * next);
    let balance = |next: f64| {
        let mean_velocity = 0.5 * (velocity + next);
        let mean_temperature = 0.5 * (temperature + next_temperature(next));
        let thermal = (6000.0 * (mean_temperature.recip() - 300.0_f64.recip())).exp();
        let force = 2.0 * thermal * (1.0 + mean_velocity.abs() / 2.0) * mean_velocity;
        2.0 * (velocity - next) / 0.0625 - force
    };
    let (mut lower, mut upper) = (0.0, velocity);
    assert!(balance(lower) > 0.0 && balance(upper) < 0.0);
    for _ in 0..80 {
        let trial = lower + 0.5 * (upper - lower);
        if balance(trial) > 0.0 {
            lower = trial;
        } else {
            upper = trial;
        }
    }
    let next = lower + 0.5 * (upper - lower);
    (next, next_temperature(next))
}

#[test]
fn nonlinear_drag_heats_its_own_body_and_replays_the_same_scene() {
    let path = example("thermally_softening_drag.json");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    let (mut velocity, mut temperature) = (4.0_f64, 300.0_f64);
    let mut saved = None;
    scene_context(|cx, _| {
        for step in 1..=32 {
            let (next_velocity, next_temperature) = reference_step(velocity, temperature);
            let force = 2.0 * (velocity - next_velocity) / 0.0625;
            let heat = (velocity * velocity - next_velocity * next_velocity) / 0.0625;
            (velocity, temperature) = (next_velocity, next_temperature);
            assert!(scene.advance(f64::MAX, cx).unwrap());
            close(probe(&scene, "mass.velocity_m_s"), velocity);
            close(probe(&scene, "body.temperature_k"), temperature);
            // The flow junction gives the damper minus the mass velocity,
            // so its a-port effort is the signed braking force on the mass.
            close(probe(&scene, "drag.force_n"), -force);
            close(probe(&scene, "drag.heat_w"), heat);
            close(velocity * velocity + 2.0 * (temperature - 300.0), 16.0);
            assert!(scene.ledger().energy_residual_j().abs() < 1.0e-7);
            if step == 11 {
                saved = Some(scene.checkpoint().unwrap());
            }
        }
        assert!(scene.is_finished());
        assert!(velocity > 0.0 && velocity < 1.0);
        assert!(temperature > 307.0 && temperature < 308.0);
        assert!(scene.ledger().entropy_generated_j_per_k() > 0.0);
        let (bytes, root) = saved.unwrap();
        let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        resumed.resume_checkpoint(&bytes, root).unwrap();
        resumed.run(cx).unwrap();
        assert_eq!(resumed.render_json(), scene.render_json());

        let changed = text.replace("\"quadratic_flow_scale\": 2", "\"quadratic_flow_scale\": 4");
        assert_ne!(text, changed);
        let mut foreign = fs_cli::SystemScene::from_text(&changed, &path).unwrap();
        let before = foreign.render_json();
        assert!(foreign.resume_checkpoint(&bytes, root).is_err());
        assert_eq!(foreign.render_json(), before);
    });
    let (code, output, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(output, scene.render_json());
}

fn prescribed_scene(law: &str, temperature: f64) -> String {
    format!(
        r#"{{"schema":"frankensim.system.v0","duration_s":0.125,"step_s":0.125,
          "elements":{{"loss":{{"kind":"damper","damping_n_s_m":3.0,
            "dissipation_law":{law}}}}},
          "connections":[{{"prescribe":"loss.a","velocity_m_s":2.0}},
            {{"prescribe":"loss.b","velocity_m_s":0.0}},
            {{"prescribe":"loss.heat","temperature_k":{temperature}}}]}}"#
    )
}

#[test]
fn dissipation_declarations_and_temperatures_outside_support_refuse_without_advancing() {
    for law in [
        "{}",
        r#"{"quadratic_flow_scale":0.0}"#,
        r#"{"quadratic_flow_scale":-1.0}"#,
        r#"{"quadratic_flow_scale":1.0,"unknown":2}"#,
        r#"{"shift":{"arrhenius":{"activation_energy_j_mol":24000.0,"reference_k":300.0,"valid_k":[310.0,400.0]}}}"#,
    ] {
        let error = fs_cli::SystemScene::from_text(&prescribed_scene(law, 300.0), "refusal.json")
            .err()
            .expect("invalid law must refuse admission");
        assert!(
            error.message.contains("elements.loss.dissipation_law"),
            "{error}"
        );
    }
    let law = r#"{"shift":{"arrhenius":{"activation_energy_j_mol":24000.0,"reference_k":300.0,"valid_k":[290.0,340.0]}}}"#;
    let mut scene =
        fs_cli::SystemScene::from_text(&prescribed_scene(law, 400.0), "outside_support.json")
            .unwrap();
    let before = scene.render_json();
    scene_context(|cx, _| {
        assert!(scene.advance(f64::MAX, cx).is_err());
    });
    assert_eq!(scene.steps(), 0);
    assert_eq!(scene.render_json(), before);
}
