//! Real-binary coverage for the authored-load projected study, not a mock runner.
use fs_cutfem::DesignBoxEdge;
use fs_topols::{GridSdf, OptimizeSettings, RobustAggregate, RobustLoadCase, evaluate_robust_design};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

static NEXT: AtomicUsize = AtomicUsize::new(0);

fn workspace() -> PathBuf {
    let root = std::env::temp_dir().join(format!("frankensim-projected-multiload-{}-{}",
        std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&root).unwrap();
    std::fs::write(root.join("loads.csv"),
        "right,0.375,0.625,0,-1,0.7\nright,0.375,0.625,0.5,0,0.3\n").unwrap();
    root
}

fn run(root: &Path, output: &str, max_solves: &str, initial: Option<&Path>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"));
    command.arg("--projected").arg(root.join(output)).arg(root.join("loads.csv"))
        .args(["3", "1", "0.6", "8", "sum", max_solves]);
    if let Some(path) = initial { command.arg(path); }
    command.output().unwrap()
}

fn field(path: &Path) -> GridSdf {
    let text = std::fs::read_to_string(path).unwrap();
    let values = text.lines().skip(1).map(|row| {
        row.split(',').nth(2).unwrap().parse::<f64>().unwrap()
    }).collect::<Vec<_>>();
    let mut field = GridSdf::from_fn(8, &|_, _| 0.0);
    assert_eq!(values.len(), field.nodes().len());
    field.nodes_mut().copy_from_slice(&values);
    field
}

fn check_final(output: &Path) {
    let field = field(&output.join("level-set.csv"));
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.7).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.5, 0.0], 0.3).unwrap(),
    ];
    let evaluation = evaluate_robust_design(&field, &cases,
        OptimizeSettings { level: 3, ..OptimizeSettings::default() },
        RobustAggregate::WeightedSum).unwrap();
    let text = std::fs::read_to_string(output.join("summary.json")).unwrap();
    let final_state = text.split("\"final\":").nth(1).unwrap();
    assert!((evaluation.volume - 0.6).abs() <= 1e-4);
    assert!(final_state.contains(&format!("\"objective\":{:.17e}", evaluation.objective)));
    assert!(final_state.contains(&format!("\"snapshot\":\"{:#018x}\"", evaluation.snapshot)));
    for compliance in evaluation.case_compliances {
        assert!(compliance > 0.0);
        assert!(final_state.contains(&format!("{compliance:.17e}")));
    }
    assert!(text.contains("\"authority\":\"estimated\""));
}

#[test]
fn actual_binary_publishes_an_accepted_same_area_multiload_update() {
    let root = workspace();
    let result = run(&root, "study", "64", None);
    assert_eq!(result.status.code(), Some(0), "{}", String::from_utf8_lossy(&result.stderr));
    let text = std::fs::read_to_string(root.join("study/summary.json")).unwrap();
    assert!(text.contains("\"accepted_updates\":1"));
    assert!(text.contains("\"status\":\"iteration_limit\""));
    assert_eq!(std::fs::read_to_string(root.join("study/trajectory.jsonl")).unwrap().lines().count(), 1);
    check_final(&root.join("study"));
}

#[test]
fn spent_budget_exports_the_real_baseline_without_claiming_an_update() {
    let root = workspace();
    let result = run(&root, "bounded", "2", None);
    assert_eq!(result.status.code(), Some(13), "{}", String::from_utf8_lossy(&result.stderr));
    let text = std::fs::read_to_string(root.join("bounded/summary.json")).unwrap();
    assert!(text.contains("\"status\":\"solve_budget\""));
    assert!(text.contains("\"solves_started\":2"));
    assert!(text.contains("\"accepted_updates\":0"));
    assert!(std::fs::read(root.join("bounded/trajectory.jsonl")).unwrap().is_empty());
    assert_eq!(std::fs::read(root.join("bounded/baseline-level-set.csv")).unwrap(),
        std::fs::read(root.join("bounded/level-set.csv")).unwrap());
    check_final(&root.join("bounded"));
}

#[test]
fn exported_geometry_can_start_a_new_study_without_altering_the_input() {
    let root = workspace();
    assert_eq!(run(&root, "first", "2", None).status.code(), Some(13));
    let initial = root.join("first/level-set.csv");
    let original = std::fs::read(&initial).unwrap();
    let result = run(&root, "warm", "2", Some(&initial));
    assert_eq!(result.status.code(), Some(13), "{}", String::from_utf8_lossy(&result.stderr));
    assert_eq!(std::fs::read(&initial).unwrap(), original);
    assert_eq!(std::fs::read(root.join("warm/baseline-level-set.csv")).unwrap(), original);
    check_final(&root.join("warm"));
}

#[test]
fn malformed_authored_geometry_refuses_before_creating_outputs() {
    let root = workspace();
    let path = root.join("bad-field.csv");
    std::fs::write(&path, "x_normalized,y_normalized,phi_normalized\n0.125,0,-1\n").unwrap();
    let result = run(&root, "invalid", "2", Some(&path));
    assert_eq!(result.status.code(), Some(1));
    assert!(!root.join("invalid").exists());
    assert!(result.stdout.is_empty());
}

#[test]
fn existing_output_is_not_overwritten_or_reported_as_success() {
    let root = workspace();
    std::fs::create_dir(root.join("occupied")).unwrap();
    std::fs::write(root.join("occupied/keep.txt"), "original").unwrap();
    let result = run(&root, "occupied", "2", None);
    assert_eq!(result.status.code(), Some(1));
    assert!(result.stdout.is_empty());
    assert!(!root.join("occupied/summary.json").exists());
    assert_eq!(std::fs::read_to_string(root.join("occupied/keep.txt")).unwrap(), "original");
}

fn run_stress(root: &Path, output: &str, max_solves: &str, limit: f64) -> Output {
    Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"))
        .arg("--projected").arg(root.join(output)).arg(root.join("loads.csv"))
        .args(["3", "1", "0.6", "8", "sum", max_solves, "--stress-limit"])
        .arg(format!("{limit:.17e}"))
        .output().unwrap()
}

#[test]
fn actual_stress_limited_update_replays_the_exported_geometry_and_samples() {
    let root = workspace();
    let result = run_stress(&root, "stress", "64", 1e8);
    assert_eq!(result.status.code(), Some(0), "{}", String::from_utf8_lossy(&result.stderr));
    let output = root.join("stress");
    check_final(&output);
    let summary = std::fs::read_to_string(output.join("summary.json")).unwrap();
    assert!(summary.contains("\"schema\":\"projected-multiload-stress-v1\""));
    assert!(summary.contains("\"accepted_updates\":1"));
    let final_stress = summary.split("\"final_stress\":").nth(1).unwrap()
        .split(",\"refusal\":").next().unwrap();
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 0.7).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.5, 0.0], 0.3).unwrap(),
    ];
    let reference = fs_topols::evaluate_robust_sampled_stress(
        &field(&output.join("level-set.csv")), &cases,
        OptimizeSettings { level: 3, ..OptimizeSettings::default() }, RobustAggregate::WeightedSum,
    ).unwrap();
    assert!(final_stress.contains("\"status\":\"sampled_feasible\""));
    assert!(final_stress.contains("\"continuous_maximum_certificate\":false"));
    assert!(final_stress.contains(&format!("\"snapshot\":\"{:#018x}\"", reference.snapshot)));
    assert!(final_stress.contains(&format!("\"worst_von_mises\":{:.17e}", reference.worst_sampled_von_mises)));
    assert!(final_stress.contains(&format!("\"case_maxima\":[{:.17e},{:.17e}]",
        reference.case_sampled_max_von_mises[0], reference.case_sampled_max_von_mises[1])));
    let trace = std::fs::read_to_string(output.join("trajectory.jsonl")).unwrap();
    assert!(trace.contains("\"sampled_stress\":{\"scope\":\"q1-positive-volume-material-cell-probes-v1\""));
    assert_eq!(trace.lines().count(), 1);
}

#[test]
fn stress_budget_stop_retains_feasible_baseline_without_extra_case_solves() {
    let root = workspace();
    let result = run_stress(&root, "stress-budget", "2", 1e8);
    assert_eq!(result.status.code(), Some(13), "{}", String::from_utf8_lossy(&result.stderr));
    let summary = std::fs::read_to_string(root.join("stress-budget/summary.json")).unwrap();
    assert!(summary.contains("\"solves_started\":2"));
    assert!(summary.contains("\"accepted_updates\":0"));
    assert_eq!(summary.matches("\"status\":\"sampled_feasible\"").count(), 2);
    assert_eq!(std::fs::read(root.join("stress-budget/baseline-level-set.csv")).unwrap(),
        std::fs::read(root.join("stress-budget/level-set.csv")).unwrap());
    assert!(std::fs::read(root.join("stress-budget/trajectory.jsonl")).unwrap().is_empty());
    assert_eq!(run(&root, "unconstrained", "2", None).status.code(), Some(13));
    let original = std::fs::read_to_string(root.join("unconstrained/summary.json")).unwrap();
    assert!(original.contains("\"schema\":\"projected-multiload-v1\""));
    assert!(!original.contains("sampled_stress") && !original.contains("sampled_feasible"));
}

#[test]
fn zero_weight_overload_refuses_without_outputs_or_input_changes() {
    let root = workspace();
    let loads = "right,0.375,0.625,0,-1,1\nright,0.375,0.625,0,4,0\n";
    std::fs::write(root.join("loads.csv"), loads).unwrap();
    assert_eq!(run(&root, "reference", "2", None).status.code(), Some(13));
    let cases = [
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, -1.0], 1.0).unwrap(),
        RobustLoadCase::new(DesignBoxEdge::Right, 0.375, 0.625, [0.0, 4.0], 0.0).unwrap(),
    ];
    let reference = fs_topols::evaluate_robust_sampled_stress(
        &field(&root.join("reference/level-set.csv")), &cases,
        OptimizeSettings { level: 3, ..OptimizeSettings::default() }, RobustAggregate::WeightedSum,
    ).unwrap();
    assert_eq!(reference.worst_stress_case, 1);
    let limit = 2.0 * reference.case_sampled_max_von_mises[0];
    assert!(reference.worst_sampled_von_mises > limit);
    let result = run_stress(&root, "refused-stress", "64", limit);
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("sampled stress limit exceeded in case 1"));
    assert!(result.stdout.is_empty());
    assert!(!root.join("refused-stress").exists());
    assert_eq!(std::fs::read_to_string(root.join("loads.csv")).unwrap(), loads);
}

#[test]
fn malformed_stress_policy_refuses_before_reading_a_missing_input() {
    let root = workspace();
    let result = Command::new(env!("CARGO_BIN_EXE_fs-marquee-elasticity-robust"))
        .arg("--projected").arg(root.join("malformed-stress")).arg(root.join("missing.csv"))
        .args(["--stress-tolerance", "0.1"]).output().unwrap();
    assert_eq!(result.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&result.stderr).contains("--stress-tolerance requires --stress-limit"));
    assert!(result.stdout.is_empty());
    assert!(!root.join("malformed-stress").exists());
}
