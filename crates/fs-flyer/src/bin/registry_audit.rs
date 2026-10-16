//! E10.0 audit-trail generator (bead wf-root-guzez.11.1): runs the
//! registry-conformance audit against the tracked data directory and
//! writes data/wright-flyer/registry-audit-v1.json (the committed
//! trail the battery verifies byte-for-byte).
//! Repro: cargo run -p fs-flyer --bin registry_audit

use fs_flyer::registryaudit::audit_registry;
use std::path::PathBuf;

fn main() {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../data/wright-flyer");
    match audit_registry(&dir, "evidence-registry-freeze-v1.json", None) {
        Ok(receipt) => {
            let out = dir.join("registry-audit-v1.json");
            std::fs::write(&out, receipt.to_json()).expect("write audit trail");
            println!(
                "{{\"suite\":\"wf-registry-audit\",\"case\":\"emit\",\"verdict\":\"{}\",\"violations\":{},\"receipt_digest\":\"{}\",\"path\":\"{}\"}}",
                receipt.verdict,
                receipt.violations,
                receipt.receipt_digest,
                out.display()
            );
        }
        Err(e) => {
            println!(
                "{{\"suite\":\"wf-registry-audit\",\"case\":\"refusal\",\"code\":\"{}\",\"message\":\"{}\"}}",
                e.code, e.message
            );
            std::process::exit(40);
        }
    }
}
