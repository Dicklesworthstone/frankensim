use std::path::PathBuf;
use std::process::Command;

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}

#[test]
fn actual_binary_uq_runs_real_cooling_samples() {
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .arg("--json")
        .arg("cooling-network-uq")
        .arg(example("fan-correlated-hotspot.json"))
        .arg(example("uq-fan-hotspot.json"))
        .output()
        .expect("run UQ binary");
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("frankensim.cooling-network-uq.result.v1"));
    assert!(stdout.contains("\"status\":\"complete\""));
    assert!(stdout.contains("\"samples_evaluated\":8"));
    assert!(stdout.contains("empirical_probability_of_compliance"));
}

#[test]
fn unknown_multivariate_dependence_refuses_before_sampling() {
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .arg("--json")
        .arg("cooling-network-uq")
        .arg(example("fan-correlated-hotspot.json"))
        .arg(example("uq-unknown-dependence.json"))
        .output()
        .expect("run UQ negative case");
    assert!(!output.status.success());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("explicit dependence") || stderr.contains("joint probability"));
}
