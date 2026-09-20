use super::*;

fn start_restoration(root: &Path, output: &str, budget: &str, bound: &str) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"))
        .arg("--projected").arg(root.join(output)).arg(root.join("loads.csv"))
        .args(["3", "2", "0.6", "8", "sum", budget, "--checkpoint", "--stress-limit", bound,
            "--restore-stress", "--restoration-reduction", "0.01"])
        .output().unwrap()
}

#[test]
fn infeasible_budget_stop_is_exported_and_resumes_without_a_budget_refund() {
    let root = workspace();
    let run = start_restoration(&root, "repair", "2", "1e-30");
    assert_eq!(run.status.code(), Some(13), "{}", String::from_utf8_lossy(&run.stderr));
    let summary = text(&root.join("repair/summary.json"));
    assert!(summary.contains("projected-multiload-restoration-v1"));
    assert!(summary.contains("\"phase\":\"restoring\""));
    assert!(summary.contains("\"restoration_updates\":0"));
    assert!(summary.contains("sampled_limit_exceeded"));
    assert!(!summary.contains("sampled_feasible"));
    let source = root.join("repair/checkpoint.fscp");
    let bytes = std::fs::read(&source).unwrap();
    let resumed = resume(&source, &root.join("continued"), "4");
    assert_eq!(resumed.status.code(), Some(13), "{}", String::from_utf8_lossy(&resumed.stderr));
    assert_eq!(bytes, std::fs::read(root.join("continued/checkpoint.fscp")).unwrap());
    assert_eq!(bytes, std::fs::read(&source).unwrap());
    assert!(text(&root.join("continued/summary.json")).contains("\"solves_started\":2"));
    assert!(text(&root.join("continued/summary.json")).contains("\"recovery_solves_started\":4"));
}

#[test]
fn already_feasible_opt_in_preserves_the_original_numerical_path() {
    let root = workspace();
    let strict = start(&root, "strict", "64", None, true);
    let restoration = start_restoration(&root, "restoration", "64", "100000000");
    assert!(matches!(strict.status.code(), Some(0 | 11)), "{}", String::from_utf8_lossy(&strict.stderr));
    assert_eq!(strict.status.code(), restoration.status.code());
    for name in ["baseline-level-set.csv", "level-set.csv", "load-cases.csv"] {
        assert_eq!(std::fs::read(root.join("strict").join(name)).unwrap(),
            std::fs::read(root.join("restoration").join(name)).unwrap());
    }
    let summary = text(&root.join("restoration/summary.json"));
    assert!(summary.contains("\"phase\":\"sampled_feasible\""));
    assert!(summary.contains("\"restoration_updates\":0"));
    let strict_summary = text(&root.join("strict/summary.json"));
    assert!(!strict_summary.contains("stress_restoration"));
}

#[test]
fn actual_restoration_step_is_labeled_and_matches_recovered_physics() {
    use fs_cutfem::DesignBoxEdge;
    use fs_topols::{GridSdf, OptimizeSettings, RobustAggregate, RobustLoadCase};
    use fs_topols::robust_descent::{MultiLoadProjectedOptimizer, MultiLoadProjectedSettings};
    use fs_topols::volume::VolumeProjectionSettings;
    let root = workspace();
    let field = GridSdf::from_fn(8, &|x, y| (y - 0.5).abs() - (0.15 + 0.3*(2.0*x - 1.0).powi(2)));
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.7).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.5, 0.0], 0.3).unwrap(),
    ];
    let fixed = field.nodes().iter().copied().enumerate().filter(|(i, _)|
        i % 9 == 0 || i % 9 == 8 || i / 9 == 0 || i / 9 == 8).collect();
    let settings = OptimizeSettings { level: 3, iterations: 1, volfrac: 0.5, move_cells: 0.35,
        nucleation_period: 4, hole_radius_cells: 1.5, ..OptimizeSettings::default() };
    let baseline = MultiLoadProjectedOptimizer::new(field.clone(), &cases, settings, RobustAggregate::WeightedSum,
        fixed, VolumeProjectionSettings { target: 0.5, tolerance: 1e-4, max_shift: 2.0, max_evaluations: 64 },
        MultiLoadProjectedSettings { max_candidates: 16, max_solves: 64,
            ..MultiLoadProjectedSettings::default() }).unwrap();
    let measured = fs_topols::evaluate_robust_sampled_stress(baseline.geometry(), &cases, settings,
        RobustAggregate::WeightedSum).unwrap();
    let limit = format!("{:.17e}", measured.worst_sampled_von_mises * (1.0 - 1e-6));
    let mut csv = "x_normalized,y_normalized,phi_normalized\n".to_string();
    for j in 0..=8 { for i in 0..=8 {
        let [x, y] = field.pos(i, j);
        csv.push_str(&format!("{x:.17e},{y:.17e},{:.17e}\n", field.node(i, j)));
    }}
    let input = root.join("neck.csv");
    std::fs::write(&input, csv.as_bytes()).unwrap();
    let run = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"))
        .arg("--projected").arg(root.join("repair")).arg(root.join("loads.csv"))
        .args(["3", "1", "0.5", "16", "sum", "64"]).arg(&input)
        .args(["--checkpoint", "--stress-limit", &limit, "--restore-stress", "--restoration-reduction", "0"])
        .output().unwrap();
    assert!(matches!(run.status.code(), Some(0 | 15)), "{}", String::from_utf8_lossy(&run.stderr));
    let trace = text(&root.join("repair/trajectory.jsonl"));
    assert_eq!(trace.lines().count(), 1, "must execute a real restoration update");
    assert!(trace.contains("\"acceptance_phase\":\"stress_restoration\""));
    let envelope = std::fs::read(root.join("repair/checkpoint.fscp")).unwrap();
    let mut recovery = 4;
    let recovered = MultiLoadProjectedOptimizer::restore_checkpoint(&envelope[65..], &mut recovery).unwrap();
    assert_eq!(recovered.restoration_updates(), 1);
    assert!(recovered.current_stress().unwrap().worst_sampled_von_mises < measured.worst_sampled_von_mises);
    assert_eq!(std::fs::read(&input).unwrap(), csv.as_bytes());
}

#[test]
fn missing_allowable_and_resume_policy_replacement_refuse_before_input_reads() {
    let root = workspace();
    let run = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"))
        .arg("--projected").arg(root.join("out")).arg(root.join("missing.csv"))
        .arg("--restore-stress").output().unwrap();
    assert_eq!(run.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&run.stderr).contains("requires --stress-limit"));
    assert!(run.stdout.is_empty());
    assert!(!root.join("out").exists());
    let changed = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"))
        .args(["--projected", "--resume"]).arg(root.join("missing.fscp")).arg(root.join("out"))
        .args(["--restoration-reduction", "0"]).output().unwrap();
    assert_eq!(changed.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&changed.stderr).contains("immutable"));
    assert!(!root.join("out").exists());
}
