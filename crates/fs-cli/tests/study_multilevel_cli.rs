//! The authored solver declaration must reach the actual study executable.
#![cfg(feature = "sdf3-study")]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use fs_blake3::ContentHash;
use fs_ledger::Ledger;
use json::JsonValue as J;
use std::fs;
use std::process::Command;

const FIXTURE: &str = include_str!("../../../examples/marquee/bracket-3d-multilevel.fsim");

#[test]
fn executable_uses_sparse_multilevel_work_for_an_actual_accepted_design() {
    let nonce = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let root = std::env::temp_dir().join(format!("fs-multilevel-cli-{}-{nonce}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let source = root.join("study.fsim");
    let database = root.join("study.db");
    // The full adaptive/recovery fixture is exercised by native library tests.
    let text = FIXTURE.replace(":schedule ((1.0 1.0) (2.0 2.0))", ":schedule ((1.0 1.0))");
    fs::write(&source, text).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "study"]).arg(&source).arg(&database).output().unwrap();
    assert_eq!(output.status.code(), Some(0), "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
    let result = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(result.str_field("status"), Some("completed"));
    assert!(result.path(&["receipt", "work", "preconditioner_galerkin_products"])
        .and_then(J::as_f64).unwrap() > 0.0, "declared multilevel cannot run as Jacobi");
    assert_eq!(result.path(&["receipt", "work", "preconditioner_operator_applications"])
        .and_then(J::as_f64), Some(0.0));
    let ledger = Ledger::open(database.to_str().unwrap()).unwrap();
    let read = |key| {
        let hash = result.path(&["receipt", key]).and_then(J::as_str)
            .and_then(ContentHash::from_hex).unwrap();
        let bytes = ledger.get_artifact_bounded(&hash, 16 * 1024 * 1024).unwrap().unwrap();
        J::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
    };
    let trace = read("iterations");
    let stages = trace.get("stages").and_then(J::as_array).unwrap();
    assert_eq!(stages.len(), 1);
    assert_eq!(stages[0].get("gradient_check_passed"), Some(&J::Bool(true)));
    let history = stages[0].get("history").and_then(J::as_array).unwrap();
    assert!(history.len() > 1, "require a real accepted update after the baseline");
    let number = |v: &J, key: &str| v.get(key).and_then(J::as_f64).unwrap();
    for pair in history.windows(2) {
        assert!(number(&pair[1], "compliance_j") <= number(&pair[0], "compliance_j"));
        assert!(number(&pair[1], "volume_fraction") <= 0.5 + 1e-8);
    }
    let design = read("design");
    let fields = design.get("displacements_m").and_then(J::as_array).unwrap();
    assert_eq!(fields.len(), 2);
    assert_ne!(fields[0], fields[1]);
    assert!(design.get("cells").and_then(J::as_array).unwrap().len() > 8);
}
