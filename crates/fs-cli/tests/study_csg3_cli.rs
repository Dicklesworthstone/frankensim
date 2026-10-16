//! Actual native study/ledger round-trip on an authored bored 3-D domain.
#![cfg(feature = "sdf3-study")]
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use fs_blake3::ContentHash;
use fs_ledger::Ledger;
use json::JsonValue as J;
use std::fs;
use std::path::Path;
use std::process::Command;

const FIXTURE: &str = include_str!("../../../examples/marquee/bracket-3d-bored.fsim");

fn run(root: &Path, name: &str) -> (J, Vec<u8>, Vec<u8>) {
    let source = root.join(format!("{name}.fsim"));
    let db = root.join(format!("{name}.db"));
    fs::write(&source, FIXTURE).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "study"]).arg(&source).arg(&db).output().unwrap();
    assert!(output.status.success(), "stdout={} stderr={}",
        String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
    let result = J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap();
    assert_eq!(result.str_field("status"), Some("completed"));
    let ledger = Ledger::open(db.to_str().unwrap()).unwrap();
    let retained = |name| {
        let hash = result.path(&["receipt", name]).and_then(J::as_str)
            .and_then(ContentHash::from_hex).unwrap();
        ledger.get_artifact_bounded(&hash, 16 * 1024 * 1024).unwrap().unwrap()
    };
    let design = retained("design");
    let iterations = retained("iterations");
    (result, design, iterations)
}

fn number(value: &J, name: &str) -> f64 { value.get(name).and_then(J::as_f64).unwrap() }

#[test]
fn authored_bore_drives_native_optimization_and_retained_fields() {
    let nonce = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
    let root = std::env::temp_dir().join(format!("fs-csg3-cli-{}-{nonce}", std::process::id()));
    fs::create_dir(&root).unwrap();
    let (_, design, iterations) = run(&root, "first");
    let data = J::parse(std::str::from_utf8(&design).unwrap()).unwrap();
    let cells = data.get("cells").and_then(J::as_array).unwrap();
    let volume: f64 = cells.iter().map(|c| number(c, "cut_volume_m3")).sum();
    let exact = 0.63 - std::f64::consts::PI * 0.12 * 0.12;
    assert!((volume - exact).abs() < 0.003, "real bore volume {volume} vs {exact}");
    let material: f64 = cells.iter().map(|c| number(c, "cut_volume_m3") * number(c, "projected_density")).sum();
    assert!(material / volume <= 0.50000001);
    let fields = data.get("displacements_m").and_then(J::as_array).unwrap();
    assert_eq!(fields.len(), 2);
    assert_ne!(fields[0], fields[1], "independent load displacements cannot be conflated");
    let trace = J::parse(std::str::from_utf8(&iterations).unwrap()).unwrap();
    let stages = trace.get("stages").and_then(J::as_array).unwrap();
    assert_eq!(stages.len(), 1);
    assert_eq!(stages[0].get("gradient_check_passed"), Some(&J::Bool(true)));
    let history = stages[0].get("history").and_then(J::as_array).unwrap();
    assert!(history.len() >= 2, "must perform a real design update after its baseline");
    for pair in history.windows(2) {
        assert!(number(&pair[1], "compliance_j") <= number(&pair[0], "compliance_j") * (1.0 + 1e-12));
    }
    let (_, repeated_design, repeated_iterations) = run(&root, "repeated");
    assert_eq!(design, repeated_design);
    assert_eq!(iterations, repeated_iterations);
}
