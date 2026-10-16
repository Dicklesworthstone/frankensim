//! Exercise the actual command, authored constraints, physical solves and UQ.
use std::path::PathBuf;
use std::process::{Command, Output};
fn root() -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..") }
fn run(method: &str, samples: &str, tail: &[&str]) -> Output {
    let root = root();
    let mut c = Command::new(env!("CARGO_BIN_EXE_equilibrium_uq"));
    c.arg(root.join("examples/equilibrium-uncertainty/joint-reliability.model"))
        .arg(root.join("examples/equilibrium-uncertainty/joint-reliability.fit"))
        .args(["--method", method, "--samples", samples, "--seed", "73", "--independent"]);
    c.args(tail).output().unwrap()
}
fn text(out: &Output) -> String {
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout.clone()).unwrap()
}
fn scalar(s: &str, key: &str) -> f64 {
    s.split(&format!("\"{key}\":")).nth(1).unwrap().split([',','}']).next().unwrap().parse().unwrap()
}
fn row<'a>(s: &'a str, name: &str) -> &'a str {
    s.split(&format!("\"name\":\"{name}\"")).nth(1).unwrap().split('}').next().unwrap()
}

#[test]
fn actual_qmc_joint_probability_retains_dependent_responses_and_every_case() {
    let tail = ["--replicates", "4", "--all-constraints", "--uniform-x", "force-N", "-0.5", "0.5",
        "--equality-tolerance", "settled-band", "0.00390625"];
    let first = run("rqmc", "512", &tail); let json = text(&first);
    assert!(json.contains("frankensim-equilibrium-reliability-v1"));
    assert_eq!(scalar(&json, "samples"), 512.0);
    assert_eq!(scalar(&json, "case_solves"), 1024.0);
    assert_eq!(scalar(&json, "unassessed_response_constraints"), 0.0);
    let joint = scalar(&json, "compliance_probability");
    assert!((joint-0.5).abs() < 0.01, "{json}");
    let a = scalar(row(&json, "force-ceiling"), "empirical_compliance_probability");
    let b = scalar(row(&json, "minimum-displacement"), "empirical_compliance_probability");
    assert!((a-0.75).abs() < 0.01 && (b-0.75).abs() < 0.01);
    assert!((joint-a*b).abs() > 0.04);
    assert!(json.contains("complete-independent-scramble-means"));
    assert!(!json.contains("\"mean_m\"") && !json.contains("\"limit_m\""));
    assert_eq!(first.stdout, run("rqmc", "512", &tail).stdout);
}

#[test]
fn actual_mc_joint_confidence_can_decide_both_sides_without_exhausting_the_budget() {
    for (x, decision, probability) in [("0", "satisfied", 1.0), ("0.5", "violated", 0.0)] {
        let json = text(&run("mc", "4096", &["--all-constraints", "--fixed-x", "force-N", x,
            "--equality-tolerance", "settled-band", "0.00390625", "--require-probability", "0.5", "--confidence-alpha", "0.05"]));
        assert!(json.contains(&format!("\"decision\":\"{decision}\"")), "{json}");
        assert_eq!(scalar(&json, "compliance_probability"), probability);
        assert!(scalar(&json, "samples") < 4096.0);
        assert_eq!(scalar(&json, "case_solves"), 2.0*scalar(&json, "samples"));
        assert!(scalar(&json, "upper") > scalar(&json, "lower"));
    }
}

#[test]
fn missing_equality_bands_wrong_selection_and_qmc_confidence_never_emit_results() {
    for tail in [vec!["--all-constraints", "--fixed-x", "force-N", "0"],
        vec!["--all-constraints", "--case", "load-a", "--fixed-x", "force-N", "0"],
        vec!["--all-constraints", "--fixed-x", "force-N", "0", "--equality-tolerance", "force-ceiling", "1"]] {
        let out = run("mc", "16", &tail); assert!(!out.status.success()); assert!(out.stdout.is_empty());
    }
    let out = run("rqmc", "32", &["--replicates", "4", "--all-constraints", "--fixed-x", "force-N", "0",
        "--equality-tolerance", "settled-band", "0.01", "--require-probability", "0.9", "--confidence-alpha", "0.05"]);
    assert!(!out.status.success()); assert!(out.stdout.is_empty());
    assert!(String::from_utf8_lossy(&out.stderr).contains("MC-only"));
}

#[test]
fn legacy_displacement_output_and_source_files_remain_unchanged() {
    let root = root(); let dir = root.join("examples/equilibrium-uncertainty");
    let before_model = std::fs::read(dir.join("joint-reliability.model")).unwrap();
    let before_design = std::fs::read(dir.join("joint-reliability.fit")).unwrap();
    let json = text(&run("mc", "16", &["--case", "load-a", "--target", "0", "--limit-m", "0.02", "--fixed-x", "force-N", "0"]));
    assert!(json.contains("frankensim-equilibrium-uq-v2"));
    assert_eq!(scalar(&json, "unassessed_response_constraints"), 5.0);
    assert!((scalar(&json, "mean_m")-1.0/64.0).abs() < 1e-11);
    assert_eq!(before_model, std::fs::read(dir.join("joint-reliability.model")).unwrap());
    assert_eq!(before_design, std::fs::read(dir.join("joint-reliability.fit")).unwrap());
}

#[test]
fn a_late_physical_failure_cannot_publish_partial_joint_or_confidence_statistics() {
    let root = root(); let examples = root.join("examples/equilibrium-uncertainty");
    let original = std::fs::read_to_string(examples.join("joint-reliability.fit")).unwrap();
    let bad = original.replace("case load-b 1 1\nload 0 0 1\n", "case load-b 1 1\nload 0 0 1000000000\n")
        .replace("variable force-N 1 1 0.5 1.5 2\n", "variable force-N 1 1 0.5 1.5 1\n")
        .replace("bind actuator-force 1 0\n", "");
    let nonce = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("frankensim-joint-refusal-{}-{nonce}", std::process::id()));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("late-failure.fit"); std::fs::write(&path, bad).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_equilibrium_uq"))
        .arg(examples.join("joint-reliability.model")).arg(&path)
        .args(["--method", "mc", "--samples", "16", "--seed", "73", "--independent", "--all-constraints",
            "--fixed-x", "force-N", "0", "--equality-tolerance", "settled-band", "0.00390625",
            "--require-probability", "0.5", "--confidence-alpha", "0.05"]).output().unwrap();
    assert!(!out.status.success()); assert!(out.stdout.is_empty());
    assert!(String::from_utf8_lossy(&out.stderr).contains("2 case solves"));
    assert_eq!(original, std::fs::read_to_string(examples.join("joint-reliability.fit")).unwrap());
}
