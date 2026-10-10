//! Exercise the actual executable, not only the native library entry point.
#![cfg(feature = "sdf3-study")]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use fs_blake3::ContentHash;
use fs_ledger::Ledger;
use json::JsonValue as J;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const FIXTURE: &str = include_str!("../../../examples/marquee/bracket-3d-embedded.fsim");

fn scratch() -> PathBuf {
    let nonce = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let root = std::env::temp_dir().join(format!("fs-embedded-cli-{}-{nonce}", std::process::id()));
    fs::create_dir(&root).unwrap();
    root
}
fn invoke(source: &Path, ledger: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_frankensim")).args(["--json", "study"])
        .arg(source).arg(ledger).output().unwrap()
}
fn number(value: &J, key: &str) -> f64 { value.get(key).and_then(J::as_f64).unwrap() }
fn run(root: &Path, name: &str) -> (Vec<u8>, Vec<u8>) {
    let source = root.join(format!("{name}.fsim"));
    let ledger_path = root.join(format!("{name}.db"));
    // One stage bounds this subprocess test; adaptive/recovery coverage remains
    // in the native library's full two-stage fixture.
    let text = FIXTURE.replace(":schedule ((1.0 1.0) (2.0 2.0))", ":schedule ((1.0 1.0))");
    fs::write(&source, text).unwrap();
    let output = invoke(&source, &ledger_path);
    assert_eq!(output.status.code(), Some(0), "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
    let result = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(result.str_field("status"), Some("completed"));
    let ledger = Ledger::open(ledger_path.to_str().unwrap()).unwrap();
    let read = |key| {
        let hash = result.path(&["receipt", key]).and_then(J::as_str)
            .and_then(ContentHash::from_hex).unwrap();
        ledger.get_artifact_bounded(&hash, 16 * 1024 * 1024).unwrap().unwrap()
    };
    (read("design"), read("iterations"))
}

#[test]
fn executable_optimizes_on_interior_supports_and_replays_the_solved_fields() {
    let root = scratch();
    let (design, iterations) = run(&root, "first");
    let data = J::parse(std::str::from_utf8(&design).unwrap()).unwrap();
    let cells = data.get("cells").and_then(J::as_array).unwrap();
    let volume: f64 = cells.iter().map(|cell| number(cell, "cut_volume_m3")).sum();
    assert!((volume - 0.66).abs() < 1e-8);
    let material: f64 = cells.iter().map(|cell|
        number(cell, "cut_volume_m3") * number(cell, "projected_density")).sum();
    assert!(material / volume <= 0.50000001);
    let nodes = data.get("physical_nodes_m").and_then(J::as_array).unwrap();
    let fields = data.get("displacements_m").and_then(J::as_array).unwrap();
    assert_eq!(fields.len(), 2);
    assert_ne!(fields[0], fields[1]);
    let tension = fields[0].as_array().unwrap();
    assert_eq!(tension.len(), 3 * nodes.len());
    assert!(nodes.iter().enumerate().any(|(i, p)| {
        p.as_array().unwrap()[0].as_f64() == Some(0.0)
            && tension[3 * i].as_f64().unwrap().abs() > 1e-8
    }), "no hidden strong clamp may replace the support at x=0.17");
    let trace = J::parse(std::str::from_utf8(&iterations).unwrap()).unwrap();
    let stages = trace.get("stages").and_then(J::as_array).unwrap();
    assert_eq!(stages.len(), 1);
    assert_eq!(stages[0].get("gradient_check_passed"), Some(&J::Bool(true)));
    let history = stages[0].get("history").and_then(J::as_array).unwrap();
    assert!(history.len() >= 2, "require an accepted update after the solved baseline");
    for pair in history.windows(2) {
        assert!(number(&pair[1], "compliance_j") <= number(&pair[0], "compliance_j"));
    }
    let repeated = run(&root, "repeated");
    assert_eq!(design, repeated.0);
    assert_eq!(iterations, repeated.1);
}

#[test]
fn malformed_embedded_support_refuses_before_creating_a_native_ledger() {
    let root = scratch();
    let source = root.join("invalid.fsim");
    let ledger = root.join("invalid.db");
    fs::write(&source, FIXTURE.replace(":fraction (0.0 0.5)", ":fraction (0.3 0.5)")).unwrap();
    let output = invoke(&source, &ledger);
    assert!(!output.status.success());
    assert!(!ledger.exists(), "a malformed fixture must not start ledger or physics work");
    let diagnostic = format!("{}{}", String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr));
    assert!(diagnostic.contains("embedded support band must align"), "{diagnostic}");
}
