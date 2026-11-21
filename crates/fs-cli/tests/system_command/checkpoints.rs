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
    checkpoint_continuation_at(scene, text, 0.125)
}

/// The scene saved at `until_s` and resumed by another process reports byte
/// for byte what one process reports that pauses at `until_s` and goes on.
/// Returns the saved root.
fn checkpoint_continuation_at(scene: &std::path::Path, text: &str, until_s: f64) -> String {
    let checkpoint = scene.with_extension("fssnap");
    let scene_path = scene.to_str().unwrap();
    let checkpoint_path = checkpoint.to_str().unwrap();
    let (code, prefix, stderr) = run_binary(&[
        "--json",
        "system",
        scene_path,
        "--until",
        &until_s.to_string(),
        "--checkpoint",
        checkpoint_path,
    ]);
    assert_eq!(code, 0, "{}: {stderr}", scene.display());
    assert!(stderr.is_empty(), "{stderr}");
    // The scene pauses at the time asked, to the bit.
    assert_eq!(
        value(&prefix, "\"schema\"", "time_s").to_bits(),
        until_s.to_bits()
    );
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
        while reference.time_s() < until_s {
            reference.advance(until_s - reference.time_s(), cx).unwrap();
        }
        reference.run(cx).unwrap();
    });
    assert_eq!(resumed, reference.render_json(), "{}", scene.display());
    root
}

#[test]
fn saved_fixed_adaptive_and_material_scenes_continue_with_identical_reports() {
    for name in [
        "controlled_sources.json",
        "closed_loop_speed.json",
        "variable_transmission.json",
        "adaptive_oscillator.json",
        "phase_change_store.json",
        "phase_change_conduction.json",
        "carreau_pipe.json",
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
        "adaptive_suspension.json",
        "adaptive_viscoelastic_mount.json",
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

/// The scenes of the whole-step claim, each with its `step_s`: a fixed step
/// that does not divide 0.125 s, an adaptive scene, a friction scene whose
/// owner shortens steps, and the bowed string.
const WHOLE_STEP_SCENES: [(&str, f64); 4] = [
    ("two_slabs.json", 0.675),
    ("adaptive_oscillator.json", 0.5),
    ("stribeck_block.json", 0.01),
    ("bowed_mode.json", 0.001),
];

/// Claims, declared before the first run: each of the four scenes, run with
/// `--until 0.125 --whole-steps --checkpoint`, stops at a time at or after
/// 0.125 s and less than one `step_s` past it, which its report gives;
/// `two_slabs.json` stops at 0.675 s, the end of its first step. The
/// checkpoint's time is that stop time: a resume asked to stop at 0 s is
/// refused by it. Resumed by another process, each scene reports byte for
/// byte what the binary reports for the scene with no pause.
#[test]
fn a_scene_saved_at_a_whole_step_and_resumed_reports_what_a_run_with_no_pause_reports() {
    for (name, step_s) in WHOLE_STEP_SCENES {
        let text = std::fs::read_to_string(example(name)).unwrap();
        assert!(text.contains(&format!("\"step_s\": {step_s}")), "{name}");
        let scene = scratch(name, &text);
        let scene_path = scene.to_str().unwrap();
        let checkpoint = scene.with_extension("fssnap");
        let checkpoint_path = checkpoint.to_str().unwrap();
        let (code, plain, stderr) = run_binary(&["--json", "system", scene_path]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let (code, prefix, stderr) = run_binary(&[
            "--json",
            "system",
            scene_path,
            "--until",
            "0.125",
            "--whole-steps",
            "--checkpoint",
            checkpoint_path,
        ]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let stopped_s = value(&prefix, "\"schema\"", "time_s");
        let steps = value(&prefix, "\"schema\"", "steps");
        eprintln!("{name}: stopped at {stopped_s} s after {steps} steps");
        assert!(
            (0.125..0.125 + step_s).contains(&stopped_s),
            "{name}: {stopped_s}"
        );
        if name == "two_slabs.json" {
            assert_eq!(
                (stopped_s.to_bits(), steps.to_bits()),
                (step_s.to_bits(), 1.0_f64.to_bits())
            );
        }
        let root = saved_root(&prefix);
        let resume = [
            "--json",
            "system",
            scene_path,
            "--resume",
            checkpoint_path,
            "--resume-root",
            &root,
        ];
        let (code, stdout, stderr) = run_binary(&[&resume[..], &["--until", "0"]].concat());
        assert_eq!(code, i32::from(fs_cli::exit::INPUT), "{name}: {stderr}");
        assert!(stdout.is_empty(), "{name}: {stdout}");
        assert!(
            stderr.contains(&format!("accepted time {stopped_s} s")),
            "{name}: {stderr}"
        );
        let (code, resumed, stderr) = run_binary(&resume);
        assert_eq!(code, 0, "{name}: {stderr}");
        assert_eq!(resumed, plain, "{name}");
    }
}

/// `--whole-steps` with no `--until`, or given twice, is a usage error with
/// no report.
#[test]
fn whole_steps_is_refused_without_until_and_when_repeated() {
    let path = example("adaptive_oscillator.json");
    let alone = "--whole-steps requires --until";
    let twice = "--whole-steps is given twice";
    for (options, said) in [
        ("--whole-steps", alone),
        ("--whole-steps --checkpoint never.fssnap", alone),
        ("--until 1 --whole-steps --whole-steps", twice),
    ] {
        let args: Vec<&str> = ["system", &path]
            .into_iter()
            .chain(options.split(' '))
            .collect();
        let (code, stdout, stderr) = run_binary(&args);
        assert_eq!(code, i32::from(fs_cli::exit::USAGE), "{args:?}: {stderr}");
        assert!(stdout.is_empty(), "{args:?}: {stdout}");
        assert!(stderr.contains(said), "{args:?}: {stderr}");
    }
}

/// The packs the conduction examples read, relative to the scene.
const STEEL_316_PACK: &str = "steel316/stainless-316-nasa.fsmcdpk";
const ATKINS_JOINT_PACK: &str = "atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk";

/// A conduction example copied into the test's scratch directory with the
/// pack it reads, if any: its path and its text.
fn conduction_scene(name: &str, pack: Option<&str>) -> (std::path::PathBuf, String) {
    let text = std::fs::read_to_string(example(name)).unwrap();
    let scene = scratch(name, &text);
    if let Some(pack) = pack {
        let beside = scene.parent().unwrap().join(pack);
        std::fs::create_dir_all(beside.parent().unwrap()).unwrap();
        std::fs::copy(example(pack), beside).unwrap();
    }
    (scene, text)
}

/// Claims, declared before the first run: every conduction example, saved
/// at 0.125 s and resumed by another process, reports byte for byte what
/// one process reports that pauses at 0.125 s and goes on
/// (`checkpoint_continuation`'s reference; the pause is itself a step
/// boundary, so this is not the report of a run with no pause, and was
/// wrongly called the uninterrupted run's here). The bodies are
/// fs-conduction's with a declared conductivity, with a card's k(T) on the
/// Newton path, and with a card's heat capacity (the enthalpy owner); the
/// contacts are a lumped conductance, trace contacts on matching and mapped
/// faces, a card-backed one and one that follows its card.
#[test]
fn saved_conduction_scenes_continue_with_identical_reports() {
    for (name, pack) in [
        ("two_slabs.json", None),
        ("two_slabs_trace.json", None),
        ("two_slabs_mapped.json", None),
        ("steel316_heated_bar.json", Some(STEEL_316_PACK)),
        ("steel316_blocks_trace.json", Some(STEEL_316_PACK)),
        ("steel316_enthalpy_bar.json", Some(STEEL_316_PACK)),
        ("two_blocks_atkins.json", Some(ATKINS_JOINT_PACK)),
        ("two_blocks_atkins_following.json", Some(ATKINS_JOINT_PACK)),
    ] {
        let (scene, text) = conduction_scene(name, pack);
        checkpoint_continuation(&scene, &text);
    }
}

/// The scene at `scene` refuses the checkpoint beside it: no report, and
/// the checkpoint diagnostic.
fn resume_is_refused(case: &str, scene: &std::path::Path, root: &str) {
    let checkpoint = scene.with_extension("fssnap");
    let (code, stdout, stderr) = run_binary(&[
        "--json",
        "system",
        scene.to_str().unwrap(),
        "--resume",
        checkpoint.to_str().unwrap(),
        "--resume-root",
        root,
    ]);
    eprintln!("{case}: exit {code}: {}", stderr.trim_end());
    assert_ne!(code, 0, "{case}");
    assert!(stdout.is_empty(), "{case}: {stdout}");
    assert!(stderr.contains("system-checkpoint"), "{case}: {stderr}");
}

/// Claims, declared before the first run: a conduction scene's checkpoint
/// is refused by the same scene with another conductivity, another grid,
/// another heat capacity, another Krylov tolerance of one body or another
/// contact resistance, and the read-once Atkins example's checkpoint is
/// refused by the example that follows its card.
#[test]
fn conduction_checkpoints_refuse_a_scene_whose_body_or_contact_changed() {
    let (scene, text) = conduction_scene("two_slabs_trace.json", None);
    let root = checkpoint_continuation(&scene, &text);
    for (case, from, to) in [
        (
            "another conductivity",
            "\"conductivity_w_m_k\": 15.0",
            "\"conductivity_w_m_k\": 15.5",
        ),
        (
            "another grid",
            "\"cells\": [4, 2, 2]",
            "\"cells\": [5, 2, 2]",
        ),
        (
            "another heat capacity",
            "\"volumetric_heat_capacity_j_m3_k\": 3600000.0",
            "\"volumetric_heat_capacity_j_m3_k\": 3700000.0",
        ),
        (
            "another Krylov tolerance",
            "\"tolerance\": 1e-13",
            "\"tolerance\": 1e-12",
        ),
        (
            "another contact resistance",
            "\"resistance_m2_k_w\": 0.0008",
            "\"resistance_m2_k_w\": 0.0009",
        ),
    ] {
        assert!(text.contains(from), "{case}: test fixture changed");
        std::fs::write(&scene, text.replacen(from, to, 1)).unwrap();
        resume_is_refused(case, &scene, &root);
    }
    let (scene, text) = conduction_scene("two_blocks_atkins.json", Some(ATKINS_JOINT_PACK));
    let root = checkpoint_continuation(&scene, &text);
    let following = std::fs::read_to_string(example("two_blocks_atkins_following.json")).unwrap();
    std::fs::write(&scene, following).unwrap();
    resume_is_refused("a contact that follows its card", &scene, &root);
}

/// Claims, declared before the first run: each friction example, saved at
/// 0.125 s and resumed by another process, reports byte for byte what one
/// process reports that pauses at 0.125 s and goes on. They are a sliding
/// interface under a prescribed slip, a Hertz patch under a spinning puck,
/// a friction interface that brings a block to rest, and a friction
/// interface that bows a modal string.
#[test]
fn saved_friction_scenes_continue_with_identical_reports() {
    for name in [
        "dragged_lining.json",
        "spun_puck.json",
        "stribeck_block.json",
        "bowed_mode.json",
    ] {
        let text = std::fs::read_to_string(example(name)).unwrap();
        let scene = scratch(name, &text);
        checkpoint_continuation(&scene, &text);
    }
}

/// The friction examples changed one passage at a time: the scene, the
/// passage and what replaces it.
const CHANGED_FRICTION_SCENES: &str = r#"
dragged_lining.json | "kinetic_mu": 0.2 | "kinetic_mu": 0.25
dragged_lining.json | "normal_n": 40.0 | "normal_n": 41.0
dragged_lining.json | [0.3, 0.4] | [0.3, 0.5]
dragged_lining.json | "normal_n": 40.0, | "normal_n": 40.0, "solve": {"work_tolerance_j": 2e-12},
stribeck_block.json | "static_mu": 0.5 | "static_mu": 0.55
stribeck_block.json | "normal_n": 20.0 | "normal_n": 21.0
stribeck_block.json | "normal_n": 20.0, | "normal_n": 20.0, "solve": {"max_iterations": 65},
spun_puck.json | "kinetic_mu": 0.3 | "kinetic_mu": 0.31
spun_puck.json | "kinetic_mu": 0.3, | "kinetic_mu": 0.3, "solve": {"maximum_evaluations": 513},
"#;

/// The bowed string changed one passage of its modal structure at a time.
const CHANGED_MODAL_SCENES: &str = r#"
bowed_mode.json | "frequency_hz": 2.0 | "frequency_hz": 2.5
bowed_mode.json | "damping_ratio": 0.0 | "damping_ratio": 0.01
bowed_mode.json | "weights": [[1.0]] | "weights": [[1.1]]
"#;

/// Claims, declared before the first run: the checkpoint of a friction
/// example is refused by the same scene with another coefficient of its
/// law, another normal force, another prescribed slip rate or another
/// solve budget ([`CHANGED_FRICTION_SCENES`]).
#[test]
fn friction_checkpoints_refuse_a_scene_whose_law_load_or_budget_changed() {
    each_changed_scene_is_refused(CHANGED_FRICTION_SCENES);
}

/// Claims, declared before the first run: the bowed string's checkpoint is
/// refused by the same scene with another frequency or damping of its mode
/// or another weight of its attachment ([`CHANGED_MODAL_SCENES`]).
#[test]
fn a_modal_checkpoint_refuses_a_scene_whose_mode_or_attachment_changed() {
    each_changed_scene_is_refused(CHANGED_MODAL_SCENES);
}

/// The struck mode changed one passage of its contact at a time.
const CHANGED_CONTACT_SCENES: &str = r#"
struck_mode.json | "exponent": 1.5 | "exponent": 1.6
struck_mode.json | "stiffness": 58423739.46721772 | "stiffness": 58423740.0
struck_mode.json | "exponent": 1.5, | "exponent": 1.5, "gap_m": 1e-07,
struck_mode.json | "maximum_penetration_m": 0.01 | "maximum_penetration_m": 0.02
"#;

/// Claims, declared before the first run: the struck mode, a stick on a
/// contact against a modal head, lasts 0.002 s. Saved at 0.001 s and
/// resumed by another process, it reports byte for byte what one process
/// reports that pauses at 0.001 s and goes on. Its checkpoint is refused by
/// the same scene with another exponent, stiffness, gap or penetration
/// budget of its contact ([`CHANGED_CONTACT_SCENES`]).
#[test]
fn a_saved_struck_mode_continues_alike_and_refuses_a_scene_whose_contact_changed() {
    each_changed_scene_is_refused_at(CHANGED_CONTACT_SCENES, 0.001);
}

/// The port-Hamiltonian oscillator changed one matrix entry at a time: its
/// energy, its structure, its dissipation and its input.
const CHANGED_PORT_HAMILTONIAN_SCENES: &str = r"
port_hamiltonian_oscillator.json | [[4.0, 0.0], [0.0, 400.0]] | [[4.0, 0.0], [0.0, 410.0]]
port_hamiltonian_oscillator.json | [[0.0, -1.0], [1.0, 0.0]] | [[0.0, -2.0], [2.0, 0.0]]
port_hamiltonian_oscillator.json | [[1.5, 0.0], [0.0, 0.0]] | [[1.6, 0.0], [0.0, 0.0]]
port_hamiltonian_oscillator.json | [[1.0], [0.0]] | [[1.1], [0.0]]
";

/// Claims, declared before the first run: the port-Hamiltonian oscillator,
/// saved at 0.125 s and resumed by another process, reports byte for byte
/// what one process reports that pauses at 0.125 s and goes on. Its
/// checkpoint is refused by the same scene with another entry of its
/// energy, its structure, its dissipation or its input matrix
/// ([`CHANGED_PORT_HAMILTONIAN_SCENES`]).
#[test]
fn a_saved_port_hamiltonian_oscillator_continues_alike_and_refuses_changed_matrices() {
    each_changed_scene_is_refused(CHANGED_PORT_HAMILTONIAN_SCENES);
}

/// Each line of `changes` (a scene, a passage and what replaces it) is a
/// scene that its example's checkpoint, saved at 0.125 s, refuses.
fn each_changed_scene_is_refused(changes: &str) {
    each_changed_scene_is_refused_at(changes, 0.125);
}

/// As [`each_changed_scene_is_refused`], with the examples saved at
/// `until_s`. An example is saved once, by [`checkpoint_continuation_at`],
/// when its first line is reached.
fn each_changed_scene_is_refused_at(changes: &str, until_s: f64) {
    let mut saved: Option<(&str, std::path::PathBuf, String, String)> = None;
    for line in changes.lines().filter(|line| !line.is_empty()) {
        let [name, from, to] = line.split(" | ").collect::<Vec<_>>()[..] else {
            panic!("{line}: three parts");
        };
        if saved.as_ref().is_none_or(|saved| saved.0 != name) {
            let text = std::fs::read_to_string(example(name)).unwrap();
            let scene = scratch(name, &text);
            let root = checkpoint_continuation_at(&scene, &text, until_s);
            saved = Some((name, scene, text, root));
        }
        let (_, scene, text, root) = saved.as_ref().unwrap();
        assert!(text.contains(from), "{line}: test fixture changed");
        std::fs::write(scene, text.replacen(from, to, 1)).unwrap();
        resume_is_refused(line, scene, root);
    }
}

/// The damped pair's lining changed one passage at a time: its law's chart,
/// a link's conductance, a cell's geometry and its step limit.
const CHANGED_TIMELINE_SCENES: &str = r#"
damped_pair_lining.json | [100000.0, 450.0, 0.0, 7850.0] | [100000.0, 460.0, 0.0, 7850.0]
damped_pair_lining.json | {"a": 1, "b": 2, "conductance_w_k": 1.0 | {"a": 1, "b": 2, "conductance_w_k": 1.5
damped_pair_lining.json | "end_m": [0.03, 0.0, 0.0] | "end_m": [0.035, 0.0, 0.0]
damped_pair_lining.json | "max_step_s": 0.125 | "max_step_s": 0.1
"#;

/// Claims, declared before the first run: the damped pair with a thermal
/// timeline for a lining, saved at 0.125 s and resumed by another process,
/// reports byte for byte what one process reports that pauses at 0.125 s and
/// goes on, and its checkpoint is refused by the same scene with another
/// law chart, link conductance, cell geometry or step limit of the lining.
#[test]
fn a_saved_thermal_timeline_continues_alike_and_refuses_a_changed_lining() {
    each_changed_scene_is_refused(CHANGED_TIMELINE_SCENES);
}

/// The packs the gas spring reads, relative to the scene.
const N2_PACKS: [&str; 3] = ["n2/species.pack", "n2/nasa9.pack", "n2/transport.fsmcdpk"];

/// The gas spring changed one passage at a time: its gas's temperature, its
/// volume and the lower end of its temperature bracket.
const CHANGED_GAS_SCENES: &str = r#"
n2_gas_spring.json | thermochemical 300 90000 | thermochemical 310 90000
n2_gas_spring.json | "volume_m3": 1.0e-4 | "volume_m3": 1.1e-4
n2_gas_spring.json | "temperature_min_k": 200 | "temperature_min_k": 210
"#;

/// Claims, declared before the first run: the piston on a `gas_cavity`,
/// saved at 0.125 s and resumed by another process, reports byte for byte
/// what one process reports that pauses at 0.125 s and goes on, and its
/// checkpoint is refused by the same scene with another gas temperature,
/// volume or temperature bracket. The scene reads its packs beside itself,
/// so they are put beside the scratch copy first.
#[test]
fn a_saved_gas_cavity_continues_alike_and_refuses_a_changed_gas() {
    let scene = scratch("n2_gas_spring.json", "");
    for pack in N2_PACKS {
        let beside = scene.parent().unwrap().join(pack);
        std::fs::create_dir_all(beside.parent().unwrap()).unwrap();
        std::fs::copy(example(pack), beside).unwrap();
    }
    each_changed_scene_is_refused(CHANGED_GAS_SCENES);
}

/// The examples the tests above save and resume by `checkpoint_continuation`
/// (28 scenes).
const SAVED_ABOVE: &str = "adaptive_oscillator adaptive_suspension \
    adaptive_viscoelastic_mount anisotropic_damping bowed_mode carreau_pipe damped_pair \
    damped_pair_lining dragged_lining phase_change_conduction phase_change_store \
    port_hamiltonian_oscillator projected_suspension self_heating_damper spun_puck \
    steel316_blocks_trace steel316_enthalpy_bar steel316_heated_bar stribeck_block struck_mode \
    two_blocks_atkins two_blocks_atkins_following two_slabs two_slabs_mapped two_slabs_trace \
    variable_transmission viscoelastic_mount n2_gas_spring";

/// Claims, declared before the first run: every other scene of
/// `examples/system`, run in place with `--until T --whole-steps
/// --checkpoint`, stops at or after T, and resumed by another process
/// reports byte for byte what the binary reports for the scene with no
/// pause. T is 0.125 s, or half the duration of a scene shorter than 0.25 s.
/// So every example is saved and resumed by a test here; a scene added to
/// the directory is saved here unless it is named above.
#[test]
fn every_other_example_saved_at_a_whole_step_and_resumed_reports_what_its_plain_run_reports() {
    let mut scenes: Vec<_> = std::fs::read_dir(example(""))
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|kind| kind == "json"))
        .collect();
    scenes.sort();
    let saves = scratch("saves", "");
    let named: Vec<&str> = SAVED_ABOVE.split_whitespace().collect();
    let mut saved = 0;
    for scene in &scenes {
        let name = scene.file_stem().unwrap().to_str().unwrap();
        if named.contains(&name) {
            continue;
        }
        let scene_path = scene.to_str().unwrap();
        let text = std::fs::read_to_string(scene).unwrap();
        let (_, duration) = text.split_once("\"duration_s\":").expect("a duration");
        let duration = duration.trim_start();
        let digits = duration.find(|c: char| c == ',' || c.is_whitespace());
        let duration_s: f64 = duration[..digits.unwrap()].parse().unwrap();
        let pause_s = if duration_s < 0.25 {
            0.5 * duration_s
        } else {
            0.125
        };
        let checkpoint = saves.with_file_name(name).with_extension("fssnap");
        let checkpoint_path = checkpoint.to_str().unwrap();
        let (code, plain, stderr) = run_binary(&["--json", "system", scene_path]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let (code, prefix, stderr) = run_binary(&[
            "--json",
            "system",
            scene_path,
            "--until",
            &pause_s.to_string(),
            "--whole-steps",
            "--checkpoint",
            checkpoint_path,
        ]);
        assert_eq!(code, 0, "{name}: {stderr}");
        let stopped_s = value(&prefix, "\"schema\"", "time_s");
        eprintln!("{name}: paused at {pause_s} s, stopped at {stopped_s} s");
        assert!(stopped_s >= pause_s, "{name}: {stopped_s}");
        let root = saved_root(&prefix);
        let (code, resumed, stderr) = run_binary(&[
            "--json",
            "system",
            scene_path,
            "--resume",
            checkpoint_path,
            "--resume-root",
            &root,
        ]);
        assert_eq!(code, 0, "{name}: {stderr}");
        assert_eq!(resumed, plain, "{name}");
        saved += 1;
    }
    assert!(saved >= 15, "{saved} scenes saved here");
}
