//! Real command -> physical adjoints -> whole-experiment information selection.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
fn examples() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/equilibrium-uncertainty")
}
fn run_with(design: &Path, tail: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_equilibrium_sensitivity"))
        .arg(examples().join("sensitivity-independent.model")).arg(design)
        .args(["--rank-relative-tolerance", "0.00001", "--max-observations", "16",
            "--max-adjoints", "16", "--point-x", "left", "0", "--point-x", "right", "0"])
        .args(tail).output().unwrap()
}
fn run(tail: &[&str]) -> Output { run_with(&examples().join("experiment-candidates.fit"), tail) }
fn text(output: &Output) -> String {
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    String::from_utf8(output.stdout.clone()).unwrap()
}
fn value(text: &str, key: &str) -> f64 {
    text.split(&format!("\"{key}\":")).nth(1).unwrap().split([',', '}']).next().unwrap().parse().unwrap()
}
fn selected(text: &str) -> &str { text.split("\"experiment_selection\":").nth(1).unwrap() }
fn copy_design(contents: &str) -> PathBuf {
    let nonce = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let directory = std::env::temp_dir().join(format!("experiment-selection-{}-{nonce}", std::process::id()));
    std::fs::create_dir(&directory).unwrap();
    let path = directory.join("candidates.fit"); std::fs::write(&path, contents).unwrap(); path
}
const SELECT: [&str; 6] = ["--select-cases", "2", "--design-ridge", "1", "--max-design-factorizations", "16"];

#[test]
fn chooses_a_complementary_case_at_an_exact_fit_without_extra_physics() {
    let json = text(&run(&SELECT)); let selection = selected(&json);
    assert!(selection.contains("\"stop_reason\":\"selection-limit\""));
    assert_eq!(value(selection, "selected_count"), 2.0);
    assert_eq!(value(selection, "factorizations"), 6.0);
    assert!(selection.contains("\"case\":\"strong-left\""));
    assert!(selection.contains("\"case\":\"independent-right\""));
    assert!(!selection.contains("\"case\":\"redundant-left\""));
    assert!(selection.contains("\"observation_indices\":[0]"));
    assert!(selection.contains("\"observation_indices\":[2]"));
    assert!((value(selection, "log_determinant_gain") - 1717.0_f64.ln()).abs() < 1e-10);
    assert!(value(&json, "objective").abs() < 1e-20);
    assert_eq!(value(&json, "evaluations"), 1.0);
    assert_eq!(value(&json, "case_solves"), 3.0);
    assert_eq!(value(&json, "observation_adjoints"), 3.0);
    assert_eq!(json, text(&run(&SELECT)));
}
#[test]
fn default_output_is_byte_identical_to_the_unextended_portion_and_inputs_are_unchanged() {
    let path = examples().join("experiment-candidates.fit"); let before = std::fs::read(&path).unwrap();
    let baseline = text(&run(&[])); let extended = text(&run(&SELECT));
    let start = extended.rfind(",\"experiment_selection\":").unwrap();
    assert_eq!(baseline, format!("{}}}\n", &extended[..start]));
    assert_eq!(before, std::fs::read(&path).unwrap());
}
#[test]
fn limited_scoring_keeps_last_complete_selection_and_does_not_claim_cardinality() {
    let mut tail = SELECT; tail[5] = "5";
    let json = text(&run(&tail)); let selection = selected(&json);
    assert!(selection.contains("\"stop_reason\":\"factorization-budget\""));
    assert_eq!(value(selection, "requested_cases"), 2.0);
    assert_eq!(value(selection, "selected_count"), 1.0);
    assert_eq!(value(selection, "factorizations"), 4.0);
    assert_eq!(value(&json, "case_solves"), 3.0);
}
#[test]
fn zero_weight_measurements_keep_their_physical_derivatives_but_do_not_drive_selection() {
    let original = std::fs::read_to_string(examples().join("experiment-candidates.fit")).unwrap();
    let path = copy_design(&original.replace("target 1 0 0.0625 0.015625 1", "target 1 0 0.0625 0.015625 0"));
    let json = text(&run_with(&path, &SELECT));
    let selection = selected(&json);
    assert!(selection.contains("\"case\":\"redundant-left\""));
    assert!(!selection.contains("\"case\":\"independent-right\""));
    assert_eq!(value(&json, "numerical_rank"), 1.0); // full-family diagnostic has no ridge
    assert!(json.contains("\"d_displacement_d_x\":[0.00000000000000000e0,-6.25000000000000000e-2]"));
}
#[test]
fn malformed_controls_and_failed_candidates_never_publish_a_selection() {
    for tail in [vec!["--select-cases", "2"], vec!["--design-ridge", "1"],
        vec!["--select-cases", "0", "--design-ridge", "1", "--max-design-factorizations", "16"],
        vec!["--select-cases", "2", "--design-ridge", "0", "--max-design-factorizations", "16"],
        vec!["--select-cases", "2", "--design-ridge", "1", "--max-design-factorizations", "0"],
        vec!["--select-cases", "4", "--design-ridge", "1", "--max-design-factorizations", "16"]] {
        let out = run(&tail); assert!(!out.status.success()); assert!(out.stdout.is_empty());
    }
    let original = std::fs::read_to_string(examples().join("experiment-candidates.fit")).unwrap();
    let path = copy_design(&original.replace("load 1 0 4", "load 1 0 1000000000"));
    let out = run_with(&path, &SELECT);
    assert!(!out.status.success()); assert!(out.stdout.is_empty());
    assert!(String::from_utf8_lossy(&out.stderr).contains("attempted 3 cases"));
}
