//! Real binary and on-disk ledger regressions, not a surrogate optimizer.
#[allow(dead_code)]
#[path = "../src/json_read.rs"]
mod json;
use json::JsonValue as J;
use fs_blake3::ContentHash;
use fs_ledger::Ledger;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

const FIXTURE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-2d.fsim"));

fn scratch(label: &str) -> PathBuf {
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
    let path = std::env::temp_dir().join(format!("fs-native-study-{label}-{}-{nonce}", std::process::id()));
    fs::create_dir(&path).unwrap();
    path
}
fn source(dir: &Path) -> PathBuf {
    let path = dir.join("study.fsim");
    fs::write(&path, FIXTURE.replace(":mesh-level 4", ":mesh-level 2")
        .replace(":max-iterations 8", ":max-iterations 3")
        .replace(":steps 8", ":steps 3")
        .replace(":move-cells 0.35", ":move-cells 0.05")
        .replace(":nucleation-period 4", ":nucleation-period 0")).unwrap();
    path
}
fn command(verb: &str) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_frankensim"));
    command.args(["--json", verb]);
    command
}
fn document(output: &Output, exit: u8) -> J {
    assert_eq!(output.status.code(), Some(i32::from(exit)), "{}", String::from_utf8_lossy(&output.stderr));
    J::parse(std::str::from_utf8(&output.stdout).unwrap()).unwrap()
}
fn run_id(result: &J) -> &str { result.str_field("run_id").unwrap() }
fn updates(result: &J) -> f64 {
    result.path(&["receipt", "continuation", "updates_this_invocation"]).and_then(J::as_f64).unwrap()
}
fn retained(ledger: &Path, result: &J, key: &str) -> Vec<u8> {
    let hash = result.path(&["receipt", key]).and_then(J::as_str).and_then(ContentHash::from_hex).unwrap();
    Ledger::open(ledger.to_str().unwrap()).unwrap()
        .get_artifact_bounded(&hash, 16 * 1024 * 1024).unwrap().unwrap()
}

#[test]
fn native_study_chunks_keep_identical_geometry_and_rows_then_export_without_solving() {
    let dir = scratch("chunks");
    let source = source(&dir);
    let whole_db = dir.join("whole.db");
    let chunk_db = dir.join("chunk.db");
    let whole = document(&command("study").arg(&source).arg(&whole_db).output().unwrap(), fs_cli::exit::SUCCESS);
    assert_eq!(whole.str_field("status"), Some("completed"));
    assert_eq!(updates(&whole), 3.0);
    let first = document(&command("study").arg(&source).arg(&chunk_db)
        .args(["--budget", "1"]).output().unwrap(), fs_cli::exit::BUDGET);
    assert_eq!(first.str_field("status"), Some("budget-exhausted"));
    assert_eq!(updates(&first), 1.0);
    let first_design = retained(&chunk_db, &first, "design");
    let resumed = document(&command("study").arg("--resume").arg(run_id(&first)).arg(&chunk_db)
        .args(["--budget", "2"]).output().unwrap(), fs_cli::exit::SUCCESS);
    assert_eq!(updates(&resumed), 2.0);
    assert_eq!(resumed.path(&["receipt", "continuation", "legacy_prefix_updates_replayed"]).and_then(J::as_f64), Some(0.0));
    for key in ["design", "iterations"] {
        assert_eq!(retained(&whole_db, &whole, key), retained(&chunk_db, &resumed, key));
    }
    assert_eq!(retained(&chunk_db, &first, "design"), first_design);
    assert_eq!(whole.path(&["receipt", "trace_hash"]), resumed.path(&["receipt", "trace_hash"]));
    for verb in ["report", "package"] {
        let exported = document(&command(verb).arg(run_id(&resumed)).arg(&chunk_db)
            .output().unwrap(), fs_cli::exit::SUCCESS);
        assert_eq!(exported.str_field("status"), Some("ok"));
        assert_eq!(exported.str_field("study_status"), Some("completed"));
    }
    let again = document(&command("study").arg("--resume").arg(run_id(&resumed)).arg(&chunk_db)
        .output().unwrap(), fs_cli::exit::SUCCESS);
    assert_eq!(again, resumed, "a completed receipt is returned unchanged");
}

#[test]
fn public_resume_finalizes_a_committed_final_design_without_repeating_an_update() {
    let dir = scratch("finalize");
    let source = source(&dir);
    let database = dir.join("study.db");
    let complete = document(&command("study").arg(&source).arg(&database)
        .output().unwrap(), fs_cli::exit::SUCCESS);
    let predecessor = complete.path(&["receipt", "predecessor"]).and_then(J::as_str).unwrap();
    let pointer = format!("study-{predecessor}");
    let recovered = document(&command("study").arg("--resume").arg(&pointer).arg(&database)
        .output().unwrap(), fs_cli::exit::SUCCESS);
    assert_eq!(updates(&recovered), 0.0);
    assert_eq!(recovered.str_field("status"), Some("completed"));
    for key in ["design", "iterations"] {
        assert_eq!(retained(&database, &complete, key), retained(&database, &recovered, key));
    }
}

#[path = "study_checkpoint_cli/projected.rs"]
mod projected;

#[path = "study_checkpoint_cli/multi_load.rs"]
mod multi_load;
