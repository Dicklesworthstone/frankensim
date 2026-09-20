//! Actual command -> authored model -> existing equilibrium -> UQ observations.
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

const MODEL: &str = include_str!("../../../examples/equilibrium-uncertainty/linear.model");
const DESIGN: &str = include_str!("../../../examples/equilibrium-uncertainty/linear.fit");
static SERIAL: AtomicUsize = AtomicUsize::new(0);
fn inputs(model: &str, design: &str) -> (PathBuf, PathBuf) {
    let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let root = std::env::temp_dir().join(format!("equilibrium-uq-{}-{stamp}-{}",
        std::process::id(), SERIAL.fetch_add(1, Ordering::Relaxed)));
    std::fs::create_dir(&root).unwrap();
    let m = root.join("source.model"); let d = root.join("source.fit");
    std::fs::write(&m, model).unwrap(); std::fs::write(&d, design).unwrap(); (m, d)
}
fn run(model: &Path, design: &Path, flags: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_equilibrium_uq")).arg(model).arg(design).args(flags).output().unwrap()
}
fn valid_flags() -> Vec<&'static str> {
    vec!["--method", "mc", "--samples", "1024", "--seed", "73", "--case", "applied-force",
        "--target", "0", "--limit-m", "0.0078125", "--independent", "--uniform-x", "force-N", "-1", "1"]
}
fn successful(result: Output) -> String {
    assert!(result.status.success(), "{}", String::from_utf8_lossy(&result.stderr));
    String::from_utf8(result.stdout).unwrap()
}
fn number(text: &str, key: &str) -> f64 {
    text.split_once(key).unwrap().1.split([',', '}', ']']).next().unwrap().parse().unwrap()
}

#[test]
fn physical_force_uncertainty_reaches_the_selected_displacement_and_replays() {
    // Modal K=256 and physical port shape=2 imply u=4F/256=F/64.
    // Uniform x in [-1,1] means F in [0,1] N, not [-1,1] N.
    let (m, d) = inputs(MODEL, DESIGN);
    let output = successful(run(&m, &d, &valid_flags()));
    assert_eq!(output.lines().count(), 1);
    assert!(output.contains("\"evidence\":\"Estimated\""));
    assert!(output.contains("\"unit\":\"m\""));
    assert!((number(&output, "\"mean_m\":") - 1.0 / 128.0).abs() < 0.001);
    assert!((number(&output, "\"displacement_std_dev_m\":") - 1.0 / (64.0 * 12.0_f64.sqrt())).abs() < 0.0006);
    assert!((number(&output, "\"compliance_probability\":") - 0.5).abs() < 0.08);
    assert_eq!(number(&output, "\"samples\":"), 1024.0);
    assert_eq!(number(&output, "\"case_solves\":"), 1024.0);
    let (other_m, other_d) = inputs(MODEL, DESIGN);
    assert_eq!(successful(run(&other_m, &other_d, &valid_flags())), output);
    assert_eq!(std::fs::read_to_string(m).unwrap(), MODEL);
    assert_eq!(std::fs::read_to_string(d).unwrap(), DESIGN);
}

#[test]
fn fixed_shared_contact_parameters_execute_all_cases_without_a_replacement_solver() {
    let model = include_str!("../../fs-couple/examples/equilibrium-design.model");
    let design = include_str!("../../fs-couple/examples/equilibrium-design.fit");
    let (m, d) = inputs(model, design);
    let output = successful(run(&m, &d, &["--method", "mc", "--samples", "4", "--seed", "93",
        "--case", "load-1.6N", "--target", "1", "--limit-m", "0.00014", "--independent",
        "--fixed-x", "gap-m", "0.5", "--fixed-x", "support-N-per-m", "0.5",
        "--fixed-x", "contact-N-per-m2", "0.8"]));
    // This observation is the RECEIVER under the SECOND load case, not the
    // support displacement, sum of targets, or an optimizer objective.
    assert!((number(&output, "\"mean_m\":") - 0.00013472559545424769).abs() < 1e-10);
    assert_eq!(number(&output, "\"displacement_std_dev_m\":"), 0.0);
    assert_eq!(number(&output, "\"compliance_probability\":"), 1.0);
    assert_eq!(number(&output, "\"case_solves\":"), 12.0);
}

#[test]
fn invalid_laws_targets_and_physics_refuse_without_publishing_partial_statistics() {
    let (m, d) = inputs(MODEL, DESIGN);
    let base = valid_flags();
    for (from, to) in [("applied-force", "missing"), ("force-N", "unknown"), ("-1", "-3")] {
        let bad: Vec<_> = base.iter().map(|&v| if v == from { to } else { v }).collect();
        let result = run(&m, &d, &bad);
        assert!(!result.status.success()); assert!(result.stdout.is_empty());
    }
    let mut missing = valid_flags(); missing.retain(|v| *v != "--independent");
    assert!(!run(&m, &d, &missing).status.success());
    let mut missing = valid_flags(); missing.truncate(missing.len() - 4);
    assert!(!run(&m, &d, &missing).status.success());
    // Tight source energy budget stays enforced: no clipped displacement or
    // skipped draw may turn this physical refusal into a compliance report.
    let tight = MODEL.replace("1000 1000 1000 1e-10", "1e-20 1000 1000 1e-10");
    let (m, d) = inputs(&tight, DESIGN);
    let failed = run(&m, &d, &base);
    assert!(!failed.status.success()); assert!(failed.stdout.is_empty());
    assert!(!failed.stderr.is_empty());
}

#[test]
fn randomized_qmc_runs_complete_nets_and_uses_net_level_errors() {
    let (m, d) = inputs(MODEL, DESIGN);
    let mut flags = valid_flags();
    flags[1] = "rqmc";
    flags.extend(["--replicates", "4"]);
    let output = successful(run(&m, &d, &flags));
    assert!(output.contains("\"method\":\"rqmc\""));
    assert!(output.contains("\"standard_error_basis\":\"complete-independent-scramble-means\""));
    assert!(output.contains("\"displacement_std_dev_m\":null"));
    assert_eq!(number(&output, "\"completed_replicates\":"), 4.0);
    assert_eq!(number(&output, "\"case_solves\":"), 1024.0);
    assert!((number(&output, "\"mean_m\":") - 1.0 / 128.0).abs() < 1.0 / (64.0 * 256.0));
    assert_eq!(number(&output, "\"compliance_probability\":"), 0.5);
    assert_eq!(number(&output, "\"compliance_standard_error\":"), 0.0);
    assert!(number(&output, "\"mean_standard_error_m\":") < 0.0001);
    // Zero between-net compliance variation is explicitly not a certainty claim.
    assert!(output.contains("no physical-validation"));
    assert_eq!(successful(run(&m, &d, &flags)), output);
    for count in ["1", "3", "257"] {
        let mut bad = flags.clone(); *bad.last_mut().unwrap() = count;
        let result = run(&m, &d, &bad);
        assert!(!result.status.success()); assert!(result.stdout.is_empty());
    }
}

#[test]
fn independent_material_and_gap_variation_reaches_nonlinear_contact_observations() {
    let model = include_str!("../../fs-couple/examples/equilibrium-design.model");
    let design = include_str!("../../fs-couple/examples/equilibrium-design.fit");
    let (m, d) = inputs(model, design);
    let output = successful(run(&m, &d, &["--method", "rqmc", "--replicates", "4", "--samples", "256",
        "--seed", "73", "--case", "load-1.6N", "--target", "1", "--limit-m", "0.00014", "--independent",
        "--uniform-x", "support-N-per-m", "0", "0.5", "--uniform-x", "contact-N-per-m2", "0", "0.8",
        "--uniform-x", "gap-m", "0", "0.5"]));
    // Independent closed-form static balance, not the production modal inverse.
    let receiver = |support: f64, stiffness: f64, gap: f64| {
        let excess = 1.6 / support - gap;
        let compliance = 1.0 / support + 1.0 / 10000.0;
        let penetration = 2.0 * excess / (1.0 + (1.0 + 4.0 * stiffness * compliance * excess).sqrt());
        stiffness * penetration * penetration / 10000.0
    };
    let minimum = receiver(600.0, 1e8, 0.0002);
    let maximum = receiver(400.0, 1.8e8, 0.0001);
    let mean = number(&output, "\"mean_m\":");
    assert!(mean > minimum && mean < maximum);
    // Tensor Gauss quadrature of the independent physical-coordinate formula.
    assert!((mean - 0.00014043493391165647).abs() < 5e-7);
    assert_eq!(number(&output, "\"case_solves\":"), 768.0);
    assert_eq!(number(&output, "\"completed_replicates\":"), 4.0);
    assert!(number(&output, "\"mean_standard_error_m\":").is_finite());
}
