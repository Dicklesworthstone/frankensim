use super::*;

fn saved_root(report: &str) -> String {
    let root = report
        .split_once("\"resume_root\":\"")
        .expect("checkpoint root in report")
        .1
        .split('"')
        .next()
        .unwrap();
    assert_eq!(root.len(), 64);
    assert!(root.bytes().all(|byte| byte.is_ascii_hexdigit()));
    root.to_owned()
}

fn checkpoint_continuation(scene: &std::path::Path, text: &str) -> String {
    let checkpoint = scene.with_extension("fssnap");
    let scene_path = scene.to_str().unwrap();
    let checkpoint_path = checkpoint.to_str().unwrap();
    let (code, prefix, stderr) = run_binary(&[
        "--json",
        "system",
        scene_path,
        "--until",
        "0.125",
        "--checkpoint",
        checkpoint_path,
    ]);
    assert_eq!(code, 0, "{}: {stderr}", scene.display());
    assert!(stderr.is_empty(), "{stderr}");
    assert_eq!(value(&prefix, "\"schema\"", "time_s"), 0.125);
    let root = saved_root(&prefix);
    assert!(!std::fs::read(&checkpoint).unwrap().is_empty());

    let (code, resumed, stderr) = run_binary(&[
        "--json",
        "system",
        scene_path,
        "--resume",
        checkpoint_path,
        "--resume-root",
        &root,
    ]);
    assert_eq!(code, 0, "{}: {stderr}", scene.display());
    assert!(stderr.is_empty(), "{stderr}");
    let mut reference = fs_cli::SystemScene::from_text(text, scene).unwrap();
    scene_context(|cx, _| {
        while reference.time_s() < 0.125 {
            reference.advance(0.125 - reference.time_s(), cx).unwrap();
        }
        reference.run(cx).unwrap();
    });
    assert_eq!(resumed, reference.render_json(), "{}", scene.display());
    root
}

#[test]
fn saved_fixed_adaptive_and_material_scenes_continue_with_identical_reports() {
    for name in [
        "variable_transmission.json",
        "adaptive_oscillator.json",
        "phase_change_store.json",
        "phase_change_conduction.json",
    ] {
        let text = std::fs::read_to_string(example(name)).unwrap();
        let scene = scratch(name, &text);
        checkpoint_continuation(&scene, &text);
    }
}

#[test]
fn saved_mechanical_networks_preserve_motion_memory_and_routed_heat() {
    for name in [
        "damped_pair.json",
        "projected_suspension.json",
        "anisotropic_damping.json",
        "viscoelastic_mount.json",
        "self_heating_damper.json",
    ] {
        let text = std::fs::read_to_string(example(name)).unwrap();
        let scene = scratch(name, &text);
        checkpoint_continuation(&scene, &text);
    }
}

#[test]
fn saved_thermal_maxwell_memory_continues_under_both_original_calibrations() {
    let mount = std::fs::read_to_string(example("viscoelastic_mount.json")).unwrap();
    for (name, calibration) in [
        (
            "arrhenius_mount.json",
            r#"{"arrhenius":{"activation_energy_j_mol":166289.2523630648,"reference_k":300,"valid_k":[280,360]}}"#,
        ),
        (
            "wlf_mount.json",
            r#"{"wlf":{"c1":4,"c2_k":100,"reference_k":300,"valid_k":[280,360]}}"#,
        ),
    ] {
        // The branch is already preloaded, and its released heat changes the
        // temperature law that acts on the same branch on subsequent steps.
        let text = mount
            .replace(
                "\"connectors\": [",
                &format!("\"shift\": {calibration}, \"connectors\": ["),
            )
            .replace(
                "{\"route\": \"mount.heat0\"",
                "{\"signal\": \"casing.temperature\", \"to\": [\"mount.temperature0\"]},\n    \
                 {\"route\": \"mount.heat0\"",
            );
        assert!(text.contains("\"shift\"") && text.contains("mount.temperature0"));
        let scene = scratch(name, &text);
        checkpoint_continuation(&scene, &text);
    }
}

#[test]
fn material_card_checkpoints_resume_original_laws_and_refuse_changed_source_packs() {
    for (name, text) in [
        ("card_body.json", source_body_scene(0.0, 2000.0, 2.0)),
        ("card_layer.json", source_layer_scene()),
    ] {
        let pack = system_material_pack(2.0, false, false);
        let scene = source_scene_file(name, &text, &pack);
        let root = checkpoint_continuation(&scene, &text);
        let checkpoint = scene.with_extension("fssnap");
        // The JSON stays byte-identical. Re-admission must bind the actual
        // prepared source, including material receipts, rather than its path.
        std::fs::write(
            scene.parent().unwrap().join("material.fsmcdpk"),
            system_material_pack(3.0, false, false),
        )
        .unwrap();
        let (code, stdout, stderr) = run_binary(&[
            "--json",
            "system",
            scene.to_str().unwrap(),
            "--resume",
            checkpoint.to_str().unwrap(),
            "--resume-root",
            &root,
        ]);
        assert_ne!(code, 0);
        assert!(stdout.is_empty());
        assert!(stderr.contains("system-checkpoint"), "{stderr}");
    }
}

#[test]
fn checkpoint_roots_corruption_and_changed_boundary_sources_are_refused() {
    let text = std::fs::read_to_string(example("heated_body.json")).unwrap();
    let scene = scratch("boundary.json", &text);
    let checkpoint = scene.with_extension("fssnap");
    let scene_path = scene.to_str().unwrap();
    let checkpoint_path = checkpoint.to_str().unwrap();
    let (code, report, stderr) = run_binary(&[
        "--json",
        "system",
        scene_path,
        "--until",
        "5",
        "--checkpoint",
        checkpoint_path,
    ]);
    assert_eq!(code, 0, "{stderr}");
    let root = saved_root(&report);
    let original = std::fs::read(&checkpoint).unwrap();
    let mut wrong_root = root.clone();
    wrong_root.replace_range(..1, if &root[..1] == "0" { "1" } else { "0" });
    let (code, stdout, stderr) = run_binary(&[
        "--json",
        "system",
        scene_path,
        "--resume",
        checkpoint_path,
        "--resume-root",
        &wrong_root,
    ]);
    assert_ne!(code, 0);
    assert!(stdout.is_empty());
    assert!(stderr.contains("system-checkpoint"), "{stderr}");

    let damaged = scene.with_file_name("damaged.fssnap");
    let mut bytes = original.clone();
    let last = bytes.len() - 1;
    bytes[last] ^= 1;
    std::fs::write(&damaged, bytes).unwrap();
    let (code, stdout, _) = run_binary(&[
        "--json",
        "system",
        scene_path,
        "--resume",
        damaged.to_str().unwrap(),
        "--resume-root",
        &root,
    ]);
    assert_ne!(code, 0);
    assert!(stdout.is_empty());

    // This changes only a prescribed closure, which the kernel's component
    // identity alone deliberately does not describe.
    let changed = text.replace(
        "\"prescribe\": \"wall.b\", \"temperature_k\": 293.15",
        "\"prescribe\": \"wall.b\", \"temperature_k\": 294.15",
    );
    assert_ne!(changed, text);
    std::fs::write(&scene, changed).unwrap();
    let (code, stdout, stderr) = run_binary(&[
        "--json",
        "system",
        scene_path,
        "--resume",
        checkpoint_path,
        "--resume-root",
        &root,
    ]);
    assert_ne!(code, 0);
    assert!(stdout.is_empty());
    assert!(stderr.contains("system-checkpoint"), "{stderr}");
    assert_eq!(std::fs::read(checkpoint).unwrap(), original);
}

#[test]
fn checkpoint_output_never_overwrites_an_existing_file() {
    let text = std::fs::read_to_string(example("adaptive_oscillator.json")).unwrap();
    let scene = scratch("scene.json", &text);
    let checkpoint = scratch("existing.fssnap", "retained checkpoint");
    let (code, stdout, stderr) = run_binary(&[
        "--json",
        "system",
        scene.to_str().unwrap(),
        "--until",
        "0",
        "--checkpoint",
        checkpoint.to_str().unwrap(),
    ]);
    assert_eq!(code, fs_cli::exit::INPUT as i32, "{stderr}");
    assert!(stdout.is_empty());
    assert!(stderr.contains("system-checkpoint-file"), "{stderr}");
    assert_eq!(
        std::fs::read_to_string(checkpoint).unwrap(),
        "retained checkpoint"
    );
}

#[test]
fn resume_requires_an_independent_root_and_valid_pause_horizon() {
    let path = example("adaptive_oscillator.json");
    let zero_root = "0".repeat(64);
    for args in [
        vec!["system", &path, "--resume", "untrusted.fssnap"],
        vec!["system", &path, "--resume-root", &zero_root],
        vec!["system", &path, "--until", "NaN"],
        vec!["system", &path, "--until", "-1"],
        vec!["system", &path, "--until", "1", "--until", "2"],
    ] {
        let (code, stdout, _) = run_binary(&args);
        assert_eq!(code, fs_cli::exit::USAGE as i32);
        assert!(stdout.is_empty());
    }
    let (code, stdout, stderr) = run_binary(&["system", &path, "--until", "5"]);
    assert_eq!(code, fs_cli::exit::INPUT as i32);
    assert!(stdout.is_empty());
    assert!(stderr.contains("system-until"), "{stderr}");
}
