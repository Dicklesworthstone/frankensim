//! E10.0 battery (bead wf-root-guzez.11.1): the registry-conformance
//! audit. Real tracked registry audits CONFORMANT with per-row
//! oracles; the POST-HOC-WIDENING hostile twin is EXECUTED (a band
//! widened in a frozen artifact is enumerated with its identity
//! history — frozen vs observed hash); the doctored-registry twin
//! refuses before judging rows; file-count caps at cap AND cap+1;
//! the committed audit trail matches a regeneration byte for byte.
//! Repro: cargo test -p fs-flyer --test registryaudit_battery

use fs_flyer::registryaudit::{MAX_FROZEN_FILES, audit_registry};
use std::fs;
use std::path::PathBuf;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-registryaudit\",\"case\":\"{case}\",{payload}}}");
}

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/wright-flyer")
}

fn scratch(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("wf-registryaudit-{tag}-{}", std::process::id()));
    let _ = fs::create_dir_all(&d);
    d
}

const REGISTRY: &str = "evidence-registry-freeze-v1.json";

#[test]
fn real_registry_audits_conformant_per_row() {
    let receipt = audit_registry(&data_dir(), REGISTRY, None).unwrap();
    assert_eq!(receipt.verdict, "CONFORMANT");
    assert_eq!(receipt.violations, 0);
    assert_eq!(
        receipt.rows.len(),
        8,
        "the v1 freeze covers eight artifacts"
    );
    assert_eq!(receipt.registry_status, "FROZEN");
    for row in &receipt.rows {
        assert!(
            row.conformant,
            "{}: {} vs {}",
            row.file, row.frozen_hex, row.current_hex
        );
        assert_eq!(row.frozen_hex, row.current_hex, "{}", row.file);
    }
    // Citing the CORRECT id passes; determinism.
    let again = audit_registry(&data_dir(), REGISTRY, Some(&receipt.evidence_registry_id)).unwrap();
    assert_eq!(
        again.receipt_digest, receipt.receipt_digest,
        "bit-identical twice"
    );
    jlog(
        "conformant",
        &format!(
            "\"evidence_registry_id\":\"{}\",\"receipt_digest\":\"{}\"",
            receipt.evidence_registry_id, receipt.receipt_digest
        ),
    );
}

#[test]
fn post_hoc_widening_hostile_twin_is_enumerated_with_identity_history() {
    // Copy the tracked data dir, then WIDEN a frozen band in the
    // canard dossier (the exact attack the freeze exists to stop).
    let dir = scratch("widen");
    for entry in fs::read_dir(data_dir()).unwrap() {
        let e = entry.unwrap();
        if e.path().is_file() {
            fs::copy(e.path(), dir.join(e.file_name())).unwrap();
        }
    }
    let victim = dir.join("canard-mechanics-v1.json");
    let text = fs::read_to_string(&victim).unwrap();
    let widened = text.replacen("0.", "9.", 1); // any numeric widening
    assert_ne!(text, widened, "the attack must actually edit bytes");
    fs::write(&victim, widened).unwrap();
    let receipt = audit_registry(&dir, REGISTRY, None).unwrap();
    assert_eq!(receipt.verdict, "VIOLATED");
    assert_eq!(receipt.violations, 1, "exactly the widened artifact");
    let bad = receipt
        .rows
        .iter()
        .find(|r| !r.conformant)
        .expect("the violation is enumerated");
    assert_eq!(bad.file, "canard-mechanics-v1.json");
    // Identity history: BOTH hashes present and different.
    assert_ne!(bad.frozen_hex, bad.current_hex);
    assert_eq!(bad.frozen_hex.len(), 64);
    assert_eq!(bad.current_hex.len(), 64);
    // Every OTHER row still conformant (per-item, not totals-only).
    for row in receipt.rows.iter().filter(|r| r.file != bad.file) {
        assert!(row.conformant, "{}", row.file);
    }
    jlog(
        "hostile-widening",
        &format!(
            "\"violated\":\"{}\",\"frozen\":\"{}\",\"observed\":\"{}\"",
            bad.file, bad.frozen_hex, bad.current_hex
        ),
    );
}

#[test]
fn doctored_registry_refuses_before_judging_rows() {
    let dir = scratch("doctor");
    for entry in fs::read_dir(data_dir()).unwrap() {
        let e = entry.unwrap();
        if e.path().is_file() {
            fs::copy(e.path(), dir.join(e.file_name())).unwrap();
        }
    }
    // The honest id of the TRACKED registry.
    let honest = audit_registry(&data_dir(), REGISTRY, None).unwrap();
    // Doctor the copied registry (swap a pinned hash) and cite the
    // honest id: the audit must refuse — a doctored registry never
    // gets to judge rows.
    let reg = dir.join(REGISTRY);
    let text = fs::read_to_string(&reg).unwrap();
    let doctored = text.replacen("0a75f54d", "deadbeef", 1);
    assert_ne!(text, doctored);
    fs::write(&reg, doctored).unwrap();
    let err = audit_registry(&dir, REGISTRY, Some(&honest.evidence_registry_id)).unwrap_err();
    assert_eq!(err.code, "registry-identity-mismatch");
    jlog("doctored-registry", &format!("\"code\":\"{}\"", err.code));
}

#[test]
fn caps_and_refusals() {
    let dir = scratch("caps");
    // AT the cap admits: a synthetic registry with exactly 64 files.
    let mk = |n: usize| {
        let rows: Vec<String> = (0..n)
            .map(|i| format!("  \"f{i}.json\": \"{:064x}\",", i))
            .collect();
        format!(
            "{{\n \"status\": \"FROZEN\",\n \"files\": {{\n{}\n }}\n}}\n",
            rows.join("\n").trim_end_matches(',')
        )
    };
    fs::write(dir.join("at-cap.json"), mk(MAX_FROZEN_FILES)).unwrap();
    let at = audit_registry(&dir, "at-cap.json", None).unwrap();
    assert_eq!(at.rows.len(), MAX_FROZEN_FILES);
    assert_eq!(
        at.verdict, "VIOLATED",
        "synthetic files are absent, honestly violated"
    );
    // One past the cap refuses.
    fs::write(dir.join("over-cap.json"), mk(MAX_FROZEN_FILES + 1)).unwrap();
    assert_eq!(
        audit_registry(&dir, "over-cap.json", None)
            .unwrap_err()
            .code,
        "registry-files-invalid"
    );
    // Empty files block refuses.
    fs::write(
        dir.join("empty.json"),
        "{\n \"status\": \"X\",\n \"files\": {\n }\n}\n",
    )
    .unwrap();
    assert_eq!(
        audit_registry(&dir, "empty.json", None).unwrap_err().code,
        "registry-files-invalid"
    );
    // Missing registry refuses.
    assert_eq!(
        audit_registry(&dir, "absent.json", None).unwrap_err().code,
        "registry-unreadable"
    );
    // Shapeless registry refuses.
    fs::write(dir.join("shapeless.json"), "{}").unwrap();
    assert_eq!(
        audit_registry(&dir, "shapeless.json", None)
            .unwrap_err()
            .code,
        "registry-shape-invalid"
    );
    jlog("caps", &format!("\"max_files\":{MAX_FROZEN_FILES}"));
}

#[test]
fn committed_audit_trail_matches_regeneration() {
    let tracked = data_dir().join("registry-audit-v1.json");
    let receipt = audit_registry(&data_dir(), REGISTRY, None).unwrap();
    let regenerated = receipt.to_json();
    let committed = fs::read_to_string(&tracked).unwrap_or_default();
    assert_eq!(
        committed, regenerated,
        "the committed audit trail must be the regeneration, byte for byte \
         (regenerate data/wright-flyer/registry-audit-v1.json via the \
         registry_audit bin and commit it with the change that moved it)"
    );
    jlog(
        "audit-trail",
        &format!("\"digest\":\"{}\"", receipt.receipt_digest),
    );
}
