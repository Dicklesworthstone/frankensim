use super::*;

#[test]
fn adaptive_rigid_suspension_converges_to_its_independent_vector_motion() {
    let original = std::fs::read_to_string(example("adaptive_suspension.json")).unwrap();
    let mut errors = Vec::new();
    let mut steps = Vec::new();
    for (name, absolute) in [
        ("loose_suspension.json", "1e-5"),
        ("tight_suspension.json", "1e-7"),
    ] {
        let text = original.replace("\"absolute\": 1e-5", &format!("\"absolute\": {absolute}"));
        let path = scratch(name, &text);
        let (code, report, stderr) = run_binary(&["--json", "system", path.to_str().unwrap()]);
        assert_eq!(code, 0, "{stderr}");
        let time = value(&report, "\"schema\"", "time_s");
        assert!((time - 1.0).abs() < 1e-12);
        // Projection n=(0.6,-0.8,0) has q'' + (18/2)q=0,
        // q(0)=0.2, q'(0)=0.02. Both perpendicular velocities stay free.
        let q = 0.2 * (3.0 * time).cos() + (0.02 / 3.0) * (3.0 * time).sin();
        let rate = -0.6 * (3.0 * time).sin() + 0.02 * (3.0 * time).cos();
        let body = element("carriage");
        let positions =
            ["x", "y", "z"].map(|axis| value(&report, &body, &format!("body0.position_{axis}_m")));
        let velocities = ["x", "y", "z"]
            .map(|axis| value(&report, &body, &format!("body0.velocity_{axis}_m_s")));
        let projection = 0.6 * positions[0] - 0.8 * positions[1];
        let projected_rate = 0.6 * velocities[0] - 0.8 * velocities[1];
        errors.push((projection - q).hypot((projected_rate - rate) / 3.0));
        steps.push(value(&report, "\"schema\"", "steps"));
        assert!((0.8 * velocities[0] + 0.6 * velocities[1] - 0.36).abs() < 1e-10);
        assert!((positions[2] - (0.4 - 0.15 * time)).abs() < 1e-10);
        assert!((velocities[2] + 0.15).abs() < 1e-10);
        assert!((value(&report, &element("spring"), "extension_m") + projection).abs() < 1e-10);
        assert!(value(&report, "\"adaptive\"", "max_error_norm") <= 1.0);
        assert!(value(&report, "\"adaptive\"", "trial_steps") > steps[steps.len() - 1]);
        assert!(value(&report, "\"totals\"", "energy_residual_j").abs() < 1e-9);
    }
    assert!(
        errors[1] < 0.25 * errors[0],
        "loose/tight errors {errors:?}"
    );
    assert!(errors[1] < 5e-5, "tight error {}", errors[1]);
    assert!(steps[1] > steps[0], "loose/tight steps {steps:?}");
}

#[test]
fn adaptive_maxwell_network_preserves_independent_memory_motion_and_heat() {
    let path = example("adaptive_viscoelastic_mount.json");
    let (code, report, stderr) = run_binary(&["--json", "system", &path]);
    assert_eq!(code, 0, "{stderr}");
    let time = value(&report, "\"schema\"", "time_s");
    // u'=-c*u-k*x, x'=u-x/tau, from u=0 and x=0.5. This two-mode
    // closed form is independent of the native network's exponential solver.
    let (c, k, inverse_tau) = (0.4_f64, 1.69_f64, 5.0_f64);
    let separation = ((c + inverse_tau).powi(2) - 4.0 * (c * inverse_tau + k)).sqrt();
    let slow = (-(c + inverse_tau) + separation) * 0.5;
    let fast = (-(c + inverse_tau) - separation) * 0.5;
    let coefficient = -k * 0.5 / separation;
    let speed = coefficient * ((slow * time).exp() - (fast * time).exp());
    let acceleration = coefficient * (slow * (slow * time).exp() - fast * (fast * time).exp());
    let extension = -(acceleration + c * speed) / k;
    let memory_energy = 0.5 * k * extension * extension;
    let heat = 0.5 * k * 0.25 - 0.5 * speed * speed - memory_energy;
    let mount = element("mount");
    let actual_speed = value(&report, &mount, "body1.velocity_x_m_s")
        - value(&report, &mount, "body0.velocity_x_m_s");
    assert!((time - 0.25).abs() < 1e-12);
    assert!((actual_speed - speed).abs() < 1e-7);
    assert!((value(&report, &mount, "connector0.branch0.energy_j") - memory_energy).abs() < 1e-7);
    assert!((value(&report, &mount, "dissipated_j") - heat).abs() < 1e-7);
    assert!(
        (value(&report, &element("casing"), "temperature_k") - 300.0 - heat / 5.0).abs() < 1e-7
    );
    assert!(value(&report, "\"adaptive\"", "max_error_norm") <= 1.0);
    assert!(value(&report, "\"totals\"", "energy_residual_j").abs() < 1e-9);

    let text = std::fs::read_to_string(&path).unwrap();
    let mut native = fs_cli::SystemScene::from_text(&text, &path).unwrap();
    scene_context(|cx, _| native.run(cx).unwrap());
    assert_eq!(native.render_json(), report);
}

#[test]
fn network_accuracy_requires_rotation_and_full_memory_accounts() {
    let text = std::fs::read_to_string(example("adaptive_viscoelastic_mount.json")).unwrap();
    let missing_rotation = text.replace(
        "      \"mount.body1.orientation_rad\": {\"absolute\": 1e-7, \"relative\": 0},\n",
        "",
    );
    let missing_memory = text.replace(
        "      \"mount.connector0.branch0.amplitude_sqrt_j\": {\"absolute\": 1e-7, \"relative\": 0},\n",
        "",
    );
    let energy_only = text.replace("branch0.amplitude_sqrt_j", "branch0.energy_j");
    for (name, changed, expected) in [
        (
            "missing_rotation.json",
            missing_rotation,
            "mount.body1.orientation_rad",
        ),
        (
            "missing_memory.json",
            missing_memory,
            "mount.connector0.branch0.amplitude_sqrt_j",
        ),
        (
            "energy_only.json",
            energy_only,
            "mount.connector0.branch0.energy_j",
        ),
    ] {
        assert_ne!(changed, text);
        let (code, stderr) = refused(name, &changed);
        assert_eq!(code, 3, "{stderr}");
        assert!(stderr.contains(expected), "{stderr}");
    }
}

#[test]
fn exhausted_network_accuracy_trials_preserve_the_entire_accepted_scene() {
    let original = std::fs::read_to_string(example("adaptive_suspension.json")).unwrap();
    let text = original
        .replace("\"max_attempts\": 32", "\"max_attempts\": 1")
        .replace("\"absolute\": 1e-5", "\"absolute\": 1e-30");
    let mut scene = fs_cli::SystemScene::from_text(&text, "retained/suspension.json").unwrap();
    let before = scene.render_json();
    scene_context(|cx, _| {
        let refusal = scene.advance(0.25, cx).unwrap_err();
        assert_eq!(refusal.code, "system-adaptive");
        assert_eq!(scene.render_json(), before);
        assert_eq!(scene.time_s(), 0.0);
        assert_eq!(scene.steps(), 0);
    });
}
