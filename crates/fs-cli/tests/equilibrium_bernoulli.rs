//! Actual mechanical solves and explicit likelihood-mixture stopping selection.
use std::path::PathBuf;
use std::process::{Command, Output};

fn root() -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..") }
fn run(joint: bool, method: Option<&str>, requirement: &str, law: &[&str], samples: &str) -> Output {
    let examples = root().join("examples/equilibrium-uncertainty");
    let mut c = Command::new(env!("CARGO_BIN_EXE_equilibrium_uq"));
    c.arg(examples.join("joint-reliability.model")).arg(examples.join("joint-reliability.fit"))
        .args(["--method", "mc", "--samples", samples, "--seed", "73", "--independent",
            "--require-probability", requirement, "--confidence-alpha", "0.05"]);
    if joint { c.args(["--all-constraints", "--equality-tolerance", "settled-band", "0.00390625"]); }
    else { c.args(["--case", "load-a", "--target", "0", "--limit-m", "0.02"]); }
    if let Some(method) = method { c.args(["--confidence-method", method]); }
    c.args(law).output().unwrap()
}
fn text(out: &Output) -> String {
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout.clone()).unwrap()
}
fn scalar(text: &str, key: &str) -> f64 {
    text.split(&format!("\"{key}\":")).nth(1).unwrap()
        .split([',', '}']).next().unwrap().parse().unwrap()
}

#[test]
fn selected_displacement_reaches_high_reliability_without_changing_its_physics() {
    let law = ["--uniform-x", "force-N", "-0.1", "0.1"];
    let out = run(false, Some("bernoulli-mixture"), "0.99", &law, "2048");
    let json = text(&out);
    assert!(json.contains("frankensim-equilibrium-uq-v2"));
    assert!(json.contains("\"decision\":\"satisfied\""));
    assert!(json.contains("beta-half-bernoulli-mixture-confidence-sequence"));
    assert_eq!(scalar(&json, "samples"), 1024.0);
    assert_eq!(scalar(&json, "case_solves"), 2048.0);
    assert!(scalar(&json, "lower") > 0.99 && scalar(&json, "lower") < 1.0);
    assert_eq!(scalar(&json, "upper"), 1.0);
    assert!((scalar(&json, "mean_m") - 1.0/64.0).abs() < 0.002);
    assert_eq!(out.stdout, run(false, Some("bernoulli-mixture"), "0.99", &law, "2048").stdout);
    let generic = text(&run(false, None, "0.99", &law, "2048"));
    assert!(generic.contains("\"decision\":\"inconclusive\""));
    assert_eq!(scalar(&generic, "samples"), 2048.0);
}

#[test]
fn joint_event_preserves_all_response_requirements_in_both_decision_directions() {
    let dir = root().join("examples/equilibrium-uncertainty");
    let model = std::fs::read(dir.join("joint-reliability.model")).unwrap();
    let design = std::fs::read(dir.join("joint-reliability.fit")).unwrap();
    for (x, required, decision, observed) in [("0", "0.99", "satisfied", 1.0), ("0.5", "0.01", "violated", 0.0)] {
        let json = text(&run(true, Some("bernoulli-mixture"), required,
            &["--fixed-x", "force-N", x], "2048"));
        assert!(json.contains("frankensim-equilibrium-reliability-v1"));
        assert!(json.contains(&format!("\"decision\":\"{decision}\"")));
        assert_eq!(scalar(&json, "samples"), 1024.0);
        assert_eq!(scalar(&json, "case_solves"), 2048.0);
        assert_eq!(scalar(&json, "compliance_probability"), observed);
        assert_eq!(scalar(&json, "unassessed_response_constraints"), 0.0);
        assert_eq!(json.matches("\"failure_count\"").count(), 5);
        assert!(scalar(&json, "upper") > scalar(&json, "lower"));
    }
    assert_eq!(model, std::fs::read(dir.join("joint-reliability.model")).unwrap());
    assert_eq!(design, std::fs::read(dir.join("joint-reliability.fit")).unwrap());
}

#[test]
fn explicit_gaussian_retains_the_default_output_and_short_bernoulli_runs_are_inconclusive() {
    let law = ["--fixed-x", "force-N", "0"];
    let default = run(true, None, "0.99", &law, "19");
    let explicit = run(true, Some("gaussian-mixture"), "0.99", &law, "19");
    text(&default); text(&explicit); assert_eq!(default.stdout, explicit.stdout);
    let json = text(&run(true, Some("bernoulli-mixture"), "0.99", &law, "19"));
    assert!(json.contains("\"decision\":\"inconclusive\""));
    assert_eq!(scalar(&json, "samples"), 19.0);
    assert!(scalar(&json, "lower") < 0.99);
}

#[test]
fn unpaired_unknown_duplicate_and_qmc_method_choices_refuse_before_files_are_opened() {
    let base = ["does-not-exist.model", "does-not-exist.fit", "--method", "mc", "--samples", "64",
        "--seed", "73", "--independent", "--all-constraints"];
    for (tail, diagnostic) in [
        (vec!["--confidence-method", "bernoulli-mixture"], "requires --require-probability"),
        (vec!["--confidence-method", "best-after-peeking"], "must be gaussian-mixture or bernoulli-mixture"),
        (vec!["--confidence-method", "gaussian-mixture", "--confidence-method", "bernoulli-mixture"], "duplicate option"),
    ] {
        let out = Command::new(env!("CARGO_BIN_EXE_equilibrium_uq")).args(base).args(tail).output().unwrap();
        assert!(!out.status.success()); assert!(out.stdout.is_empty());
        assert!(String::from_utf8_lossy(&out.stderr).contains(diagnostic));
    }
    let mut qmc = base; qmc[3] = "rqmc";
    let out = Command::new(env!("CARGO_BIN_EXE_equilibrium_uq")).args(qmc)
        .args(["--replicates", "4", "--confidence-method", "bernoulli-mixture",
            "--require-probability", "0.99", "--confidence-alpha", "0.05"]).output().unwrap();
    assert!(!out.status.success()); assert!(out.stdout.is_empty());
    assert!(String::from_utf8_lossy(&out.stderr).contains("MC-only"));
}
