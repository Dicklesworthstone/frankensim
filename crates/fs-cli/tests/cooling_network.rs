//! Execute the real binary dispatch, not only the library producers.
use std::process::Command;

#[test]
fn binary_admits_a_file_driven_mixed_solid_network() {
    let request = concat!(env!("CARGO_MANIFEST_DIR"), "/../../examples/cooling-network/mixed-slab.json");
    let output = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network", request]).output().expect("start binary");
    assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
    let result = String::from_utf8(output.stdout).unwrap();
    assert!(result.contains("\"schema\":\"frankensim.cooling-network.result.v1\""));
    assert!(result.contains("\"solid_temperatures_k\":["));
    assert!(result.contains("\"dmean_dinlet_k\":["));
    assert!(output.stderr.is_empty());
}

#[test]
fn binary_help_and_missing_request_have_stable_exit_classes() {
    let help = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["cooling-network", "--help"]).output().unwrap();
    assert!(help.status.success());
    assert!(String::from_utf8(help.stdout).unwrap().contains("tetrahedral solid"));
    let bad = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(["--json", "cooling-network"]).output().unwrap();
    assert_eq!(bad.status.code(), Some(2));
    assert!(bad.stdout.is_empty());
}
