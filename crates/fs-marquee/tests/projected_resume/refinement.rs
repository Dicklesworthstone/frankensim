use super::*;
use fs_topols::robust_descent::MultiLoadProjectedOptimizer;

fn refine(source: &Path, output: &Path, solves: &str, pause: Option<&str>) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"));
    cmd.args(["--projected", "--refine"]).arg(source).arg(output)
        .args(["--updates", "2", "--max-solves", solves, "--recovery-solves", "4"]);
    if let Some(count) = pause { cmd.args(["--pause-after", count]); }
    cmd.output().unwrap()
}

fn restored(path: &Path) -> MultiLoadProjectedOptimizer {
    let bytes = std::fs::read(path).unwrap();
    assert_eq!(fs_ledger::hash_bytes(&bytes[65..]).to_string().as_bytes(), &bytes[..64]);
    MultiLoadProjectedOptimizer::restore_checkpoint(&bytes[65..], &mut 128).unwrap()
}

#[test]
fn refine_runs_real_fine_physics_and_preserves_stress_checkpoint_resume() {
    let root = workspace();
    let coarse = start(&root, "coarse", "64", Some("1"), true);
    assert_eq!(coarse.status.code(), Some(14), "{}", String::from_utf8_lossy(&coarse.stderr));
    let source = root.join("coarse/checkpoint.fscp");
    let source_bytes = std::fs::read(&source).unwrap();
    let paused = refine(&source, &root.join("fine-paused"), "64", Some("0"));
    assert_eq!(paused.status.code(), Some(14), "{}", String::from_utf8_lossy(&paused.stderr));
    let parent = restored(&source);
    let fine = restored(&root.join("fine-paused/checkpoint.fscp"));
    assert_eq!((parent.settings().level, fine.settings().level), (3, 4));
    assert_eq!(fine.next_iteration(), 0);
    assert_eq!(fine.solves_started(), 2);
    assert_eq!(fine.load_cases(), parent.load_cases());
    assert_eq!(fine.stress_limit(), parent.stress_limit());
    assert_eq!(fine.current(), *fine.baseline());
    let independently_solved = fs_topols::evaluate_robust_sampled_stress(
        fine.geometry(), fine.load_cases(), fine.settings(), fine.aggregate()).unwrap();
    assert_eq!(fine.current_stress(), Some(&independently_solved));
    let handoff = text(&root.join("fine-paused/refinement.json"));
    assert!(handoff.contains("\"coarse_level\":3,\"fine_level\":4"));
    assert!(handoff.contains("\"source_updates\":1"));
    assert!(handoff.contains("\"cross_grid_improvement_claimed\":false"));
    let full = refine(&source, &root.join("fine-full"), "64", None);
    assert!(matches!(full.status.code(), Some(0 | 11)), "{}", String::from_utf8_lossy(&full.stderr));
    assert!(restored(&root.join("fine-full/checkpoint.fscp")).next_iteration() > 0,
        "the fixture must execute a real accepted fine-grid update");
    let continued = resume(&root.join("fine-paused/checkpoint.fscp"), &root.join("fine-resumed"), "4");
    assert_eq!(continued.status.code(), full.status.code());
    for name in ["baseline-level-set.csv", "level-set.csv", "load-cases.csv", "checkpoint.fscp", "trajectory.jsonl"] {
        assert_eq!(std::fs::read(root.join("fine-full").join(name)).unwrap(),
            std::fs::read(root.join("fine-resumed").join(name)).unwrap(), "{name}");
    }
    assert_eq!(std::fs::read(source).unwrap(), source_bytes);
}

#[test]
fn new_fine_budget_is_explicit_and_source_budget_remains_exhausted() {
    let root = workspace();
    assert_eq!(start(&root, "coarse", "2", Some("0"), false).status.code(), Some(14));
    let source = root.join("coarse/checkpoint.fscp");
    let bytes = std::fs::read(&source).unwrap();
    let result = refine(&source, &root.join("fine"), "2", None);
    assert_eq!(result.status.code(), Some(13), "{}", String::from_utf8_lossy(&result.stderr));
    let fine = restored(&root.join("fine/checkpoint.fscp"));
    assert_eq!(fine.settings().level, 4);
    assert_eq!((fine.next_iteration(), fine.solves_started(), fine.max_solves()), (0, 2, 2));
    assert_eq!(std::fs::read(&source).unwrap(), bytes);
    assert_eq!(resume(&source, &root.join("source-resumed"), "4").status.code(), Some(13));
}

#[test]
fn refinement_refuses_corruption_underfunding_and_policy_replacement_without_outputs() {
    let root = workspace();
    assert_eq!(start(&root, "coarse", "2", Some("0"), false).status.code(), Some(14));
    let source = root.join("coarse/checkpoint.fscp");
    let mut bytes = std::fs::read(&source).unwrap();
    *bytes.last_mut().unwrap() ^= 1;
    std::fs::write(root.join("bad.fscp"), bytes).unwrap();
    let bad = root.join("bad.fscp");
    for (input, target, solves) in [(&bad, "bad", "2"), (&source, "underfunded", "1")] {
        let result = refine(input, &root.join(target), solves, None);
        assert_eq!(result.status.code(), Some(1));
        assert!(result.stdout.is_empty());
        assert!(!root.join(target).exists());
    }
    let result = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"))
        .args(["--projected", "--refine"]).arg(root.join("missing.fscp")).arg(root.join("out"))
        .args(["--updates", "2", "--stress-limit", "100", "--recovery-solves", "4"])
        .output().unwrap();
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("physical policies are inherited"));
    assert!(!root.join("out").exists());
    std::fs::create_dir(root.join("occupied")).unwrap();
    std::fs::write(root.join("occupied/keep"), b"original").unwrap();
    assert_eq!(refine(&source, &root.join("occupied"), "2", None).status.code(), Some(1));
    assert_eq!(std::fs::read(root.join("occupied/keep")).unwrap(), b"original");
}
