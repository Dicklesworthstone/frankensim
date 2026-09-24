//! Admission must finish before reserving output or launching a physical solve.
//! Keep these zero-observation cases separate from the costly solver replay tests.
use std::{fs, path::{Path, PathBuf}, process::{Command, Output}, time::{SystemTime, UNIX_EPOCH}};

fn example(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/cooling-network").join(name)
}

fn command(base: &Path, replicates: Option<&str>) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.arg("--json").arg("cooling-network-uq").arg(base).arg(example("uq-fan-hotspot.json"));
    if let Some(replicates) = replicates { command.args(["--qmc-replicates", replicates]); }
    command
}

fn empty_checkpoint(name: &str) -> (PathBuf, PathBuf, Vec<u8>) {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let dir = std::env::temp_dir().join(format!("fs-qmc-admission-{name}-{}-{nonce}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    let checkpoint = dir.join("original.bin");
    let output = command(&example("fan-correlated-hotspot.json"), Some("2"))
        .arg("--checkpoint").arg(&checkpoint).args(["--max-new-samples", "0"]).output().unwrap();
    assert_eq!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)),
        "{}", String::from_utf8_lossy(&output.stderr));
    let bytes = fs::read(&checkpoint).unwrap();
    assert!(bytes.starts_with(b"FSQMC001"));
    assert_eq!(bytes.len(), 81, "zero allowance must retain no physical observations");
    (dir, checkpoint, bytes)
}

fn refused_before_output(output: &Output, destination: &Path) {
    assert!(!output.status.success());
    assert_ne!(output.status.code(), Some(i32::from(fs_cli::exit::BUDGET)));
    assert!(output.stdout.is_empty());
    assert!(!destination.exists(), "invalid resume reserved a destination");
}

#[test]
fn changed_model_layout_and_mc_interpretation_refuse_before_file_creation() {
    let (dir, original, bytes) = empty_checkpoint("identity");
    let base = example("fan-correlated-hotspot.json");
    let fresh = dir.join("must-not-exist.bin");
    let changed_base = dir.join("changed.json");
    fs::write(&changed_base, format!("{}\n", fs::read_to_string(&base).unwrap())).unwrap();
    // Even semantically neutral base edits change the exact fixed-input identity.
    let changed = command(&changed_base, Some("2")).arg("--resume").arg(&original)
        .arg("--checkpoint").arg(&fresh).output().unwrap();
    refused_before_output(&changed, &fresh);
    assert!(String::from_utf8_lossy(&changed.stderr).contains("identity"));
    for replicates in [Some("4"), None] {
        let output = command(&base, replicates).arg("--resume").arg(&original)
            .arg("--checkpoint").arg(&fresh).output().unwrap();
        refused_before_output(&output, &fresh);
    }
    assert_eq!(fs::read(&original).unwrap(), bytes);
}

#[test]
fn corrupt_oversized_and_colliding_checkpoints_cannot_destroy_the_retained_prefix() {
    let (dir, original, bytes) = empty_checkpoint("bytes");
    let base = example("fan-correlated-hotspot.json");
    let fresh = dir.join("must-not-exist.bin");
    let damaged = dir.join("damaged.bin");
    let mut corrupt = bytes.clone();
    let last = corrupt.len() - 1; corrupt[last] ^= 1;
    fs::write(&damaged, &corrupt).unwrap();
    let output = command(&base, Some("2")).arg("--resume").arg(&damaged)
        .arg("--checkpoint").arg(&fresh).output().unwrap();
    refused_before_output(&output, &fresh);
    assert!(String::from_utf8_lossy(&output.stderr).contains("checksum"));
    let oversized = dir.join("oversized.bin");
    fs::write(&oversized, vec![0_u8; 1_000_000]).unwrap();
    let output = command(&base, Some("2")).arg("--resume").arg(&oversized)
        .arg("--checkpoint").arg(&fresh).output().unwrap();
    refused_before_output(&output, &fresh);
    let collision = command(&base, Some("2")).arg("--resume").arg(&original)
        .arg("--checkpoint").arg(&original).output().unwrap();
    assert!(!collision.status.success()); assert!(collision.stdout.is_empty());
    assert_eq!(fs::read(&original).unwrap(), bytes);
    assert_eq!(fs::read(&damaged).unwrap(), corrupt);
}
