//! G1/G3: radiative exchange runs through the common scene and checkpoint path.

use super::{element, example, run_binary, scene_context, value};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * (1.0 + expected.abs()),
        "{actual} != {expected}"
    );
}

fn probe(scene: &fs_cli::SystemScene, path: &str) -> f64 {
    scene.probe(path).unwrap().unwrap()
}

#[test]
fn gray_radiation_scenes_follow_the_surface_law_in_both_directions() {
    for (a, b) in [(600.0_f64, 300.0_f64), (300.0, 600.0), (400.0, 400.0)] {
        let text = format!(
            r#"{{"schema":"frankensim.system.v0","duration_s":0.5,"step_s":0.125,
              "elements":{{"exchange":{{"kind":"gray_radiation","area_m2":0.25,
                "emissivity_a":0.8,"emissivity_b":0.6}}}},
              "connections":[{{"prescribe":"exchange.a","temperature_k":{a}}},
                {{"prescribe":"exchange.b","temperature_k":{b}}}]}}"#
        );
        let mut scene = fs_cli::SystemScene::from_text(&text, "gray_radiation.json").unwrap();
        assert_eq!(scene.probe("exchange.heat_rate_w").unwrap(), None);
        scene_context(|cx, _| scene.run(cx).unwrap());
        let expected =
            5.670_374_419e-8 * 0.25 * (a.powi(4) - b.powi(4)) / (1.0 / 0.8 + 1.0 / 0.6 - 1.0);
        close(probe(&scene, "exchange.heat_rate_w"), expected);
        assert!(scene.ledger().energy_residual_j().abs() < 1e-10);
        let report = scene.render_json();
        assert!(report.contains("full mutual view"));
        assert!(value(&report, "\"totals\"", "entropy_generated_j_per_k") >= 0.0);
    }
}

/// Solve the independent scalar energy balance for one midpoint interval.
/// The oracle uses fourth powers, not the component's secant or its tangent.
fn exchanged_energy(a: f64, b: f64, dt: f64) -> f64 {
    let coefficient = 5.670_374_419e-8 * 0.25 / (1.0 / 0.8 + 1.0 / 0.6 - 1.0);
    let residual = |energy: f64| {
        energy - dt * coefficient * ((a - energy / 2000.0).powi(4) - (b + energy / 4000.0).powi(4))
    };
    let (mut low, mut high) = (0.0, (a - b) / (1.0 / 1000.0 + 1.0 / 2000.0));
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..80 {
        let middle = low + 0.5 * (high - low);
        if residual(middle) > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    low + 0.5 * (high - low)
}

#[test]
fn gray_radiation_couples_finite_bodies_and_resumes_the_same_scene() {
    let path = example("radiative_exchange.json");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut scene = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    let (mut a, mut b) = (600.0, 300.0);
    let mut saved = None;
    scene_context(|cx, _| {
        for step in 1..=32 {
            let heat = exchanged_energy(a, b, 0.125);
            a -= heat / 1000.0;
            b += heat / 2000.0;
            assert!(scene.advance(f64::MAX, cx).unwrap());
            close(probe(&scene, "hot.temperature_k"), a);
            close(probe(&scene, "cold.temperature_k"), b);
            close(probe(&scene, "exchange.heat_rate_w"), heat / 0.125);
            assert!(scene.ledger().energy_residual_j().abs() < 1e-7);
            if step == 13 {
                saved = Some(scene.checkpoint().unwrap());
            }
        }
        assert!(scene.is_finished());
        assert!(a < 599.0 && b > 300.5);
        close(1000.0 * a + 2000.0 * b, 1_200_000.0);
        let (bytes, root) = saved.unwrap();
        let mut resumed = fs_cli::SystemScene::from_text(&text, &path).unwrap();
        resumed.resume_checkpoint(&bytes, root).unwrap();
        resumed.run(cx).unwrap();
        assert_eq!(resumed.render_json(), scene.render_json());
        let changed = text.replace("\"emissivity_a\": 0.8", "\"emissivity_a\": 0.7");
        let mut foreign = fs_cli::SystemScene::from_text(&changed, &path).unwrap();
        let before = foreign.render_json();
        assert!(foreign.resume_checkpoint(&bytes, root).is_err());
        assert_eq!(foreign.render_json(), before);
    });
    let (code, report, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    assert_eq!(report, scene.render_json());
    close(value(&report, &element("hot"), "temperature_k"), a);
}

#[test]
fn gray_radiation_admission_names_invalid_surface_parameters() {
    let path = example("radiative_exchange.json");
    let text = std::fs::read_to_string(&path).unwrap();
    for (from, to) in [
        ("\"area_m2\": 0.25", "\"area_m2\": 0.0"),
        ("\"emissivity_a\": 0.8", "\"emissivity_a\": 0.0"),
        ("\"emissivity_b\": 0.6", "\"emissivity_b\": 1.1"),
        ("\"emissivity_a\": 0.8", "\"emissivity_a\": -0.1"),
        (
            "\"emissivity_b\": 0.6",
            "\"emissivity_b\": 0.6, \"view_factor\": 0.5",
        ),
    ] {
        let changed = text.replace(from, to);
        let refusal = fs_cli::SystemScene::from_text(&changed, &path)
            .err()
            .expect("invalid radiation description");
        assert!(refusal.message.contains("elements.exchange"), "{refusal}");
    }
}
