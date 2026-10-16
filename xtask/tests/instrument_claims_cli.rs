//! E2E for `check-instrument-claims` through the REAL xtask binary (bead
//! `frankensim-music-claims-registry-mc31g`): the full row lifecycle —
//! seed ungated -> green with evidence -> live-default refusals ->
//! D21 deletion refusal -> demotion to refused — driven against a scratch
//! git repository, asserting the JSON-lines verdicts at every step
//! (agents parse verdicts, not prose). Plus one run against the live
//! repository registry, which must be clean.
//!
//! The scratch repo carries one minimal crate manifest because
//! `load_workspace` refuses an empty `crates/` directory; the check itself
//! never reads manifests, but the binary's shared preamble does.

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_xtask")
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("xtask has a parent")
        .to_path_buf()
}

/// Run `xtask check-instrument-claims` in `dir`, returning (success, stdout).
/// Stdout carries the JSON-lines verdicts; stderr is surfaced on panic so a
/// failing step can be diagnosed from the test log alone.
fn run_check(dir: &Path) -> (bool, String) {
    let output = Command::new(bin())
        .arg("check-instrument-claims")
        .current_dir(dir)
        .output()
        .expect("spawn xtask");
    let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    if !stderr.trim().is_empty() && !output.status.success() {
        eprintln!("xtask stderr:\n{stderr}");
    }
    (output.status.success(), stdout)
}

fn git(dir: &Path, args: &[&str]) {
    let output = Command::new("git")
        .args(args)
        .current_dir(dir)
        .output()
        .expect("spawn git");
    assert!(
        output.status.success(),
        "git {args:?} failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn write(dir: &Path, rel: &str, content: &str) {
    let path = dir.join(rel);
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).expect("mkdir");
    }
    std::fs::write(path, content).expect("write");
}

fn registry(rows: &str) -> String {
    format!("{{\n  \"schema\": \"frankensim-instrument-claims-v1\",\n  \"rows\": [{rows}]\n}}\n")
}

fn row(gate: &str, live_default: &str, evidence: &str, budget: &str, notes: &str) -> String {
    format!(
        "{{\"filling\":\"scratch\",\"image\":\"probe\",\"qoi\":\"lifecycle\",\
         \"owner_crates\":[\"fs-scratch\"],\"exactness\":[\"X-Exact\"],\
         \"gate\":\"{gate}\",\"live_default\":\"{live_default}\",\
         \"determinism\":\"one-host\",\"evidence\":{evidence},\
         \"budget_row\":{budget},\"corpus_refs\":[],\"notes\":\"{notes}\"}}"
    )
}

#[test]
fn live_repository_registry_is_clean() {
    let (ok, stdout) = run_check(&repo_root());
    assert!(ok, "live registry check failed:\n{stdout}");
    assert!(
        stdout.contains("\"check\":\"instrument-claims\""),
        "missing check verdicts:\n{stdout}"
    );
    assert!(stdout.contains("rows="), "missing summary note:\n{stdout}");
    // The human-facing "policy OK" line goes to stderr; the machine verdict
    // is the JSON summary — agents parse verdicts, not prose.
    assert!(
        stdout.contains("\"checks\":\"instrument-claims\"") && stdout.contains("\"violations\":0"),
        "missing zero-violation summary verdict:\n{stdout}"
    );
}

#[test]
fn scratch_lifecycle_through_the_real_binary() {
    let scratch = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join(format!("instrument-claims-e2e-{}", std::process::id()));
    if scratch.exists() {
        std::fs::remove_dir_all(&scratch).expect("clean stale scratch");
    }
    std::fs::create_dir_all(&scratch).expect("mkdir scratch");

    // Minimal workspace so the binary's manifest preamble loads.
    write(
        &scratch,
        "crates/fs-scratch/Cargo.toml",
        "[package]\nname = \"fs-scratch\"\nversion = \"0.0.1\"\nedition = \"2024\"\n\n\
         [package.metadata.frankensim]\nlayer = \"UTIL\"\n",
    );
    git(&scratch, &["init", "-q"]);
    git(&scratch, &["config", "user.email", "e2e@frankensim.test"]);
    git(&scratch, &["config", "user.name", "instrument-claims e2e"]);

    // Step 1: seed an ungated row; commit it as the predecessor.
    write(
        &scratch,
        "instrument-claims.json",
        &registry(&row("ungated", "no", "[]", "null", "seed")),
    );
    git(&scratch, &["add", "-A"]);
    git(&scratch, &["commit", "-qm", "seed"]);
    let (ok, stdout) = run_check(&scratch);
    assert!(ok, "seed row must pass:\n{stdout}");
    assert!(
        stdout.contains("rows=1 ungated=1 green=0 refused=0"),
        "summary wrong:\n{stdout}"
    );

    // Step 2: green without evidence refuses by name.
    write(
        &scratch,
        "instrument-claims.json",
        &registry(&row("green", "no", "[]", "null", "premature")),
    );
    let (ok, stdout) = run_check(&scratch);
    assert!(!ok, "green-without-evidence must fail:\n{stdout}");
    assert!(
        stdout.contains("gate=green with empty evidence"),
        "wrong refusal:\n{stdout}"
    );
    assert!(
        stdout.contains("\"crate\":\"scratch/probe/lifecycle\""),
        "refusal must name the row:\n{stdout}"
    );

    // Step 3: green with resolvable evidence passes and surfaces the
    // transition against the committed predecessor.
    let evidence = "[{\"kind\":\"test\",\"ref\":\"crates/fs-scratch/Cargo.toml\"}]";
    write(
        &scratch,
        "instrument-claims.json",
        &registry(&row("green", "no", evidence, "null", "earned")),
    );
    let (ok, stdout) = run_check(&scratch);
    assert!(ok, "evidenced green must pass:\n{stdout}");
    assert!(
        stdout.contains("ungated -> green"),
        "missing transition note:\n{stdout}"
    );

    // Step 4: live_default without a budget row refuses (D25), even green.
    write(
        &scratch,
        "instrument-claims.json",
        &registry(&row("green", "yes", evidence, "null", "no budget yet")),
    );
    let (ok, stdout) = run_check(&scratch);
    assert!(!ok, "live-default without budget must fail:\n{stdout}");
    assert!(stdout.contains("no budget_row"), "wrong refusal:\n{stdout}");

    // Step 5: with a budget row it passes; commit as the new predecessor.
    write(
        &scratch,
        "instrument-claims.json",
        &registry(&row(
            "green",
            "yes",
            evidence,
            "\"blake3:budget-fixture\"",
            "measured",
        )),
    );
    let (ok, stdout) = run_check(&scratch);
    assert!(ok, "budgeted live-default must pass:\n{stdout}");
    git(&scratch, &["add", "-A"]);
    git(&scratch, &["commit", "-qm", "green live default"]);

    // Step 6: deleting the row refuses (D21) — rows never vanish.
    write(&scratch, "instrument-claims.json", &registry(""));
    let (ok, stdout) = run_check(&scratch);
    assert!(!ok, "deletion must fail:\n{stdout}");
    assert!(stdout.contains("row deleted"), "wrong refusal:\n{stdout}");
    assert!(stdout.contains("D21"), "refusal cites doctrine:\n{stdout}");

    // Step 6b: the determinism composition lint through the binary — a
    // cross-isa claim over a default-ceiling owner refuses by name
    // (weakest-operand law; zero cross-ISA goldens exist).
    write(
        &scratch,
        "instrument-claims.json",
        &registry(
            &row("green", "no", evidence, "null", "over-claimed replay")
                .replace("\"one-host\"", "\"cross-isa\""),
        ),
    );
    let (ok, stdout) = run_check(&scratch);
    assert!(!ok, "cross-isa over-claim must fail:\n{stdout}");
    assert!(
        stdout.contains("exceeds owner crate"),
        "wrong refusal:\n{stdout}"
    );

    // Step 7: demotion to refused (with its reason) is the legal path, and
    // the transition is surfaced. live_default must drop with it.
    write(
        &scratch,
        "instrument-claims.json",
        &registry(&row(
            "refused",
            "no",
            evidence,
            "null",
            "lost the bake-off; kept per D21",
        )),
    );
    let (ok, stdout) = run_check(&scratch);
    assert!(ok, "reasoned refusal must pass:\n{stdout}");
    assert!(
        stdout.contains("green -> refused"),
        "missing transition note:\n{stdout}"
    );
    assert!(
        stdout.contains("rows=1 ungated=0 green=0 refused=1"),
        "summary wrong:\n{stdout}"
    );

    std::fs::remove_dir_all(&scratch).expect("cleanup scratch");
}
