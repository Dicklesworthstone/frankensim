//! E1.7 registry-freeze guard (bead wf-root-guzez.2.7). The hostile twin
//! against post-hoc widening: every frozen evidence artifact's fs-blake3
//! content hash is pinned in the registry; this battery recomputes each
//! from the exact bytes on disk. ANY in-place edit to a frozen artifact —
//! band widening, partition flip, claim relaxation — turns CI red here.
//! Legitimate evolution mints a v2 beside the frozen v1.
//! Repro: cargo test -p fs-flyer --test registry_freeze

use std::fs;
use std::path::PathBuf;

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/wright-flyer")
}

fn registry() -> serde_free::Registry {
    serde_free::parse(&fs::read_to_string(data_dir().join("evidence-registry-freeze-v1.json")).unwrap())
}

/// Minimal dependency-free JSON field extraction for the two structures
/// this guard needs (workspace law: no serde in production paths; a test
/// helper stays honest by parsing only what it asserts).
mod serde_free {
    pub struct Registry {
        pub files: Vec<(String, String)>,
        pub status: String,
        pub decision: String,
    }

    pub fn parse(text: &str) -> Registry {
        // The files block: "name.json": "hex" pairs inside "files": { ... }.
        let files_block = text
            .split("\"files\": {")
            .nth(1)
            .and_then(|rest| rest.split('}').next())
            .expect("registry must carry a files block");
        let mut files = Vec::new();
        for line in files_block.lines() {
            let line = line.trim().trim_end_matches(',');
            if let Some((k, v)) = line.split_once("\": \"") {
                files.push((
                    k.trim_start_matches('"').to_string(),
                    v.trim_end_matches('"').to_string(),
                ));
            }
        }
        let field = |key: &str| -> String {
            text.split(&format!("\"{key}\": \""))
                .nth(1)
                .and_then(|rest| rest.split('"').next())
                .unwrap_or_default()
                .to_string()
        };
        Registry { files, status: field("status"), decision: field("decision") }
    }
}

#[test]
fn every_frozen_artifact_hash_matches_disk() {
    let reg = registry();
    assert_eq!(reg.status, "FROZEN");
    assert_eq!(reg.files.len(), 8, "the freeze covers exactly the eight evidence artifacts");
    let mut receipts = Vec::new();
    let mut mismatches = Vec::new();
    for (name, pinned) in &reg.files {
        let bytes = fs::read(data_dir().join(name))
            .unwrap_or_else(|e| panic!("frozen artifact {name} unreadable: {e}"));
        let got = fs_blake3::hash_bytes(&bytes).to_hex();
        receipts.push(format!("{{\"file\":\"{name}\",\"blake3\":\"{got}\"}}"));
        if &got != pinned {
            mismatches.push(format!("{name}: disk {got} vs pinned {pinned}"));
        }
    }
    println!(
        "{{\"suite\":\"wf-registry-freeze\",\"case\":\"hashes\",\"receipts\":[{}]}}",
        receipts.join(",")
    );
    assert!(
        mismatches.is_empty(),
        "FROZEN artifacts changed on disk — post-hoc edits to frozen evidence are forbidden; \
         mint a v2 beside the artifact (or pin at the initial measure-then-pin):\n{}",
        mismatches.join("\n")
    );
}

#[test]
fn a7_decision_is_recorded_and_ceiling_stands() {
    let reg = registry();
    assert_eq!(reg.decision, "A7NoGo-v1", "the go/no-go must be explicit");
    // The Estimated ceiling must still stand in the canard dossier (the
    // NoGo consequence): its evidence class and ceiling strings survive.
    let canard = fs::read_to_string(data_dir().join("canard-mechanics-v1.json")).unwrap();
    assert!(canard.contains("\"evidence_class\": \"sign-tendency-only\""));
    assert!(canard.contains("\"quantitative_ceiling\": \"Estimated\""));
    println!("{{\"suite\":\"wf-registry-freeze\",\"case\":\"a7\",\"decision\":\"A7NoGo-v1\"}}");
}
