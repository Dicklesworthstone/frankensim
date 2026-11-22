//! G1/G3: authored multiport lossless laws compose through the ordinary scene path.

use super::{example, run_binary, scene_context};

const INITIAL: [f64; 3] = [1.0, -0.5, 2.0];
const PORTS: [&str; 3] = ["a", "b", "c"];

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * (1.0 + expected.abs()),
        "{actual} != {expected}"
    );
}

fn probe(scene: &fs_cli::SystemScene, path: &str) -> f64 {
    scene.probe(path).unwrap().unwrap()
}

/// With unit capacitances, q' = -J q is a rotation about (1, 1, 1).
/// This closed form uses no native component, matrix solver or time integrator.
fn rotated_charge(angle: f64) -> [f64; 3] {
    let mean = INITIAL.iter().sum::<f64>() / 3.0;
    let root_three = 3.0_f64.sqrt();
    let cross = [
        (INITIAL[2] - INITIAL[1]) / root_three,
        (INITIAL[0] - INITIAL[2]) / root_three,
        (INITIAL[1] - INITIAL[0]) / root_three,
    ];
    std::array::from_fn(|i| {
        mean + (INITIAL[i] - mean) * angle.cos() + cross[i] * angle.sin()
    })
}

#[test]
fn capacitor_scene_follows_midpoint_rotation_and_resumes_with_complete_reports() {
    let path = example("dirac_capacitors.json");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    for port in PORTS {
        for quantity in ["voltage_v", "current_a"] {
            assert_eq!(
                scene.probe(&format!("junction.{port}.{quantity}")).unwrap(),
                None
            );
        }
    }
    let (g, dt) = (0.75, 0.03125);
    let energy = 0.5 * INITIAL.iter().map(|q| q * q).sum::<f64>();
    let mut previous = INITIAL;
    let mut exchanged = [0.0_f64; 3];
    let mut saved = None;
    scene_context(|cx, _| {
        for step in 1..=128_u32 {
            assert!(scene.advance(f64::MAX, cx).unwrap());
            assert_eq!(scene.time_s(), f64::from(step) * dt);
            // The Cayley transform advances this rotation by a known angle.
            let angle = -2.0 * f64::from(step) * (3.0_f64.sqrt() * g * dt / 2.0).atan();
            let expected = rotated_charge(angle);
            let mean: [f64; 3] =
                std::array::from_fn(|i| 0.5 * (previous[i] + expected[i]));
            let flows = [
                g * (mean[2] - mean[1]),
                g * (mean[0] - mean[2]),
                g * (mean[1] - mean[0]),
            ];
            let mut charges = [0.0; 3];
            let mut power = 0.0;
            for (i, port) in PORTS.into_iter().enumerate() {
                charges[i] = probe(&scene, &format!("capacitor_{port}.charge_c"));
                close(charges[i], expected[i]);
                close(
                    probe(&scene, &format!("capacitor_{port}.voltage_v")),
                    expected[i],
                );
                let effort = probe(&scene, &format!("junction.{port}.voltage_v"));
                let flow = probe(&scene, &format!("junction.{port}.current_a"));
                close(effort, mean[i]);
                close(flow, flows[i]);
                close(
                    0.5 * (expected[i].powi(2) - previous[i].powi(2)) + dt * effort * flow,
                    0.0,
                );
                power += effort * flow;
                exchanged[i] = exchanged[i].max((charges[i] - INITIAL[i]).abs());
            }
            close(power, 0.0);
            close(charges.iter().sum(), INITIAL.iter().sum());
            close(0.5 * charges.iter().map(|q| q * q).sum::<f64>(), energy);
            close(scene.ledger().energy_residual_j(), 0.0);
            assert_eq!(scene.ledger().boundary_energy_in_j(), 0.0);
            previous = expected;
            if step == 37 {
                saved = Some(scene.checkpoint().unwrap());
            }
        }
        assert!(scene.is_finished());
        assert!(!scene.advance(f64::MAX, cx).unwrap());
        let (bytes, root) = saved.unwrap();
        let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        resumed.resume_checkpoint(&bytes, root).unwrap();
        resumed.run(cx).unwrap();
        assert_eq!(resumed.render_json(), scene.render_json());
        // Changing a still-valid skew law must not reuse the saved model.
        let changed = text.replace("0.75", "0.5");
        let mut foreign = fs_cli::SystemScene::from_text(&changed, &path).unwrap();
        let before = foreign.render_json();
        assert!(foreign.resume_checkpoint(&bytes, root).is_err());
        assert_eq!(foreign.render_json(), before);
    });
    assert!(exchanged.into_iter().all(|maximum| maximum > 0.5));
    for row in &scene.ledger().nodes {
        assert_eq!(row.dissipated_j.value(), 0.0);
        assert_eq!(row.numerical_loss_j.value(), 0.0);
    }
    for row in &scene.ledger().connections {
        close(row.energy_defect_j.value(), 0.0);
    }
    let (code, json, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(json, scene.render_json());
}

#[test]
fn named_ports_report_their_native_si_efforts_and_flows_in_every_domain() {
    for (domain, effort, flow) in [
        ("translation", "force_n", "velocity_m_s"),
        ("rotation", "torque_n_m", "angular_velocity_rad_s"),
        ("electrical", "voltage_v", "current_a"),
        ("hydraulic", "pressure_pa", "flow_m3_s"),
    ] {
        // This odd-dimensional matrix is singular; no inverse is needed.
        let text = format!(
            r#"{{"schema":"frankensim.system.v0","duration_s":0.125,"step_s":0.125,
              "elements":{{"junction":{{"kind":"dirac_junction","domain":"{domain}",
                "ports":["west","north","east"],"incidence":[[0,-2,2],[2,0,-2],[-2,2,0]]}}}},
              "connections":[{{"prescribe":"junction.west","{effort}":2}},
                {{"prescribe":"junction.north","{effort}":3}},
                {{"prescribe":"junction.east","{effort}":5}}]}}"#
        );
        for zero in [false, true] {
            let authored = if zero {
                text.replace("[[0,-2,2],[2,0,-2],[-2,2,0]]", "[[0,0,0],[0,0,0],[0,0,0]]")
            } else {
                text.clone()
            };
            let mut scene =
                fs_cli::SystemScene::from_text(&authored, "dirac-domains.json").unwrap();
            scene_context(|cx, _| scene.run(cx).unwrap());
            for (port, e, f) in [
                ("west", 2.0, 4.0),
                ("north", 3.0, -6.0),
                ("east", 5.0, 2.0),
            ] {
                close(probe(&scene, &format!("junction.{port}.{effort}")), e);
                close(
                    probe(&scene, &format!("junction.{port}.{flow}")),
                    if zero { 0.0 } else { f },
                );
            }
            close(scene.ledger().boundary_energy_in_j(), 0.0);
            close(scene.ledger().energy_residual_j(), 0.0);
        }
        let wrong_unit = text.replace(&format!("\"{effort}\":"), "\"temperature_k\":");
        let error = fs_cli::SystemScene::from_text(&wrong_unit, "dirac-wrong-unit.json")
            .err()
            .unwrap();
        assert!(error.message.contains(effort), "{error}");
    }
}

#[test]
fn adaptive_junction_needs_only_storage_accounts_and_keeps_checkpoint_history() {
    let path = example("dirac_capacitors.json");
    let text = std::fs::read_to_string(&path).unwrap().replace(
        "\"step_s\": 0.03125,",
        r#""step_s": 0.5,
          "adaptive": {"order": 2, "tolerances": {
            "capacitor_a.charge_c": {"absolute": 1e-7, "relative": 0},
            "capacitor_b.charge_c": {"absolute": 1e-7, "relative": 0},
            "capacitor_c.charge_c": {"absolute": 1e-7, "relative": 0}
          }},"#,
    );
    let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    scene_context(|cx, _| {
        assert!(scene.advance(f64::MAX, cx).unwrap());
        let (bytes, root) = scene.checkpoint().unwrap();
        let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        resumed.resume_checkpoint(&bytes, root).unwrap();
        scene.run(cx).unwrap();
        resumed.run(cx).unwrap();
        assert_eq!(resumed.render_json(), scene.render_json());
    });
    let expected = rotated_charge(-3.0_f64.sqrt() * 0.75 * scene.time_s());
    let mut energy = 0.0;
    for (port, exact) in PORTS.into_iter().zip(expected) {
        let charge = probe(&scene, &format!("capacitor_{port}.charge_c"));
        assert!((charge - exact).abs() < 2e-4, "{port}: {charge} vs {exact}");
        energy += 0.5 * charge * charge;
    }
    close(energy, 0.5 * INITIAL.iter().map(|q| q * q).sum::<f64>());
    close(scene.ledger().energy_residual_j(), 0.0);
}

#[test]
fn malformed_junctions_refuse_at_scene_admission() {
    let text = std::fs::read_to_string(example("dirac_capacitors.json")).unwrap();
    for (from, to, reason) in [
        ("[0, -0.75, 0.75]", "[1, -0.75, 0.75]", "skew"),
        ("[0, -0.75, 0.75]", "[0, 0.75, 0.75]", "skew"),
        ("[0, -0.75, 0.75]", "[0, -0.75]", "column count"),
        ("[0, -0.75, 0.75]", "[0, \"bad\", 0.75]", "finite number"),
        (
            "[\"a\", \"b\", \"c\"]",
            "[\"a\", \"a\", \"c\"]",
            "junction",
        ),
        ("[\"a\", \"b\", \"c\"]", "[\"a\"]", "ports"),
        ("[\"a\", \"b\", \"c\"]", "[]", "ports"),
        (
            "[\"a\", \"b\", \"c\"]",
            "[\"a\", \"b\", \"c\", \"d\"]",
            "row count",
        ),
        (
            "[\"a\", \"b\", \"c\"]",
            "[\"a\", \"b\", \"two words\"]",
            "ports",
        ),
        (
            "[\"a\", \"b\", \"c\"]",
            "[\"a\", \"b\", \"c.end\"]",
            "ports",
        ),
        (
            "\"domain\": \"electrical\"",
            "\"domain\": \"thermal\"",
            "domain",
        ),
        (
            "\"domain\": \"electrical\"",
            "\"domain\": \"electrical\", \"heat\": true",
            "unsupported field",
        ),
    ] {
        let changed = text.replace(from, to);
        assert_ne!(changed, text);
        let error = fs_cli::SystemScene::from_text(&changed, "invalid-dirac.json")
            .err()
            .unwrap();
        assert!(
            error.message.to_ascii_lowercase().contains(reason),
            "expected {reason:?}: {error}"
        );
    }
}
