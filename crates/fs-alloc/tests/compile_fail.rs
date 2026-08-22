//! Compile-fail battery for the affine lease authority surface
//! (frankensim-epic-bedrock-6ys.21.1.3.1, acceptance criterion 4): every
//! misuse clause must be rejected by the compiler itself — no runtime
//! fallback that could silently accept a cloned owner, a reissued charge,
//! a use-after-close, an escaped parent lifetime, or a fabricated terminal
//! receipt.
//!
//! In-house harness (no trybuild — Franken-only law, mirroring the
//! `fs-soa` battery): a scratch cargo project with a path dependency back
//! to `fs-alloc` is `cargo check`ed offline per fixture and stderr is
//! asserted to contain the expected diagnostic text.
//!
//! Typed N/A rationale: the AC's "unbounded-to-bounded evidence" clause is
//! a runtime VALUE axis, not a type axis — bounded and unbounded roots are
//! both `OperationMemoryLease`, so there is deliberately no compile-time
//! type to reject. Its nearest strict equivalent is the existing typed
//! refusal coverage in `lease_delegation.rs`
//! (`delegation_configuration_is_fallible_bounded_and_pristine_only`
//! refuses unbounded-root configuration; `delegate_from` refuses
//! unbounded parents with `unbounded_parent`), asserted there by
//! deterministic refusal codes rather than here by rustc diagnostics.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

const CASES: &[(&str, &str, &str)] = &[
    (
        "clone_delegated_owner",
        "use fs_alloc::{LeaseIdentity, OperationMemoryLease};\n\
         let root = OperationMemoryLease::bounded(8);\n\
         let root_id = LeaseIdentity::root(*b\"cfclon01\", [1; 32]);\n\
         let child_id = root_id.child([2; 32], 0).unwrap();\n\
         root.enable_delegation(root_id, \"run\", 1).unwrap();\n\
         let child = root.delegate_capacity(child_id, \"run/child\", 8).unwrap();\n\
         let _duplicate_owner = child.clone();\n",
        "no method named `clone`",
    ),
    (
        "clone_charge",
        "use fs_alloc::{LeaseIdentity, OperationMemoryLease};\n\
         let root = OperationMemoryLease::bounded(8);\n\
         let root_id = LeaseIdentity::root(*b\"cfcharge\", [1; 32]);\n\
         let child_id = root_id.child([2; 32], 0).unwrap();\n\
         root.enable_delegation(root_id, \"run\", 1).unwrap();\n\
         let child = root.delegate_capacity(child_id, \"run/child\", 8).unwrap();\n\
         let charge = child.reserve(\"payload\", 4).unwrap();\n\
         let _reissued_charge = charge.clone();\n",
        "no method named `clone`",
    ),
    (
        "reissue_after_ownership_transfer",
        "use fs_alloc::{LeaseIdentity, OperationMemoryLease};\n\
         let root = OperationMemoryLease::bounded(8);\n\
         let root_id = LeaseIdentity::root(*b\"cfxfer01\", [1; 32]);\n\
         let child_id = root_id.child([2; 32], 0).unwrap();\n\
         root.enable_delegation(root_id, \"run\", 1).unwrap();\n\
         let child = root.delegate_capacity(child_id, \"run/child\", 8).unwrap();\n\
         let outcome = child.close();\n\
         let _after_move = child.capacity_bytes();\n\
         let _ = outcome;\n",
        "borrow of moved value",
    ),
    (
        "use_after_explicit_close",
        "use fs_alloc::{LeaseIdentity, OperationMemoryLease};\n\
         let root = OperationMemoryLease::bounded(8);\n\
         let root_id = LeaseIdentity::root(*b\"cfuac001\", [1; 32]);\n\
         let child_id = root_id.child([2; 32], 0).unwrap();\n\
         root.enable_delegation(root_id, \"run\", 1).unwrap();\n\
         let child = root.delegate_capacity(child_id, \"run/child\", 8).unwrap();\n\
         let _receipt = child.close().unwrap();\n\
         let _after_close = child.capacity_bytes();\n",
        "borrow of moved value",
    ),
    (
        "parent_escape",
        "use fs_alloc::{LeaseIdentity, OperationMemoryLease};\n\
         let child = {\n\
         let root = OperationMemoryLease::bounded(8);\n\
         let root_id = LeaseIdentity::root(*b\"cfescape\", [1; 32]);\n\
         let child_id = root_id.child([2; 32], 0).unwrap();\n\
         root.enable_delegation(root_id, \"run\", 1).unwrap();\n\
         root.delegate_capacity(child_id, \"run/child\", 8).unwrap()\n\
         };\n\
         let _escaped = child.capacity_bytes();\n",
        "does not live long enough",
    ),
    (
        "fabricated_sealed_receipt",
        "use fs_alloc::{LeaseIdentity, SealedLeaseReceipt};\n\
         let root_id = LeaseIdentity::root(*b\"cfforg01\", [1; 32]);\n\
         let forged = SealedLeaseReceipt {\n\
         schema_version: 1,\n\
         root_identity: root_id,\n\
         root_id: \"run\",\n\
         limit_bytes: 8,\n\
         metadata_limit: 1,\n\
         delegation_count: 0,\n\
         direct_granted_bytes: 0,\n\
         direct_returned_bytes: 0,\n\
         delegated_bytes: 0,\n\
         returned_delegated_bytes: 0,\n\
         child_granted_bytes: 0,\n\
         child_returned_bytes: 0,\n\
         child_published_bytes: 0,\n\
         publication_record_count: 0,\n\
         published_transfer_count: 0,\n\
         rolled_back_transfer_count: 0,\n\
         refused_requests: 0,\n\
         refused_bytes: 0,\n\
         peak_used_bytes: 0,\n\
         final_used_bytes: 0,\n\
         active_delegations: 0,\n\
         release_invariant_violations: 0,\n\
         counter_overflowed: false,\n\
         seal_sequence: 0,\n\
         close_sequence: 0,\n\
         refusal_root: [0; 32],\n\
         delegation_root: [0; 32],\n\
         publication_root: [0; 32],\n\
         receipt_root: [0; 32],\n\
         };\n\
         let _ = forged;\n",
        "private field",
    ),
];

#[test]
fn lease_misuse_clauses_are_rejected_at_compile_time() {
    let alloc_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let root = std::env::temp_dir().join(format!("fs_alloc_compile_fail_{}", std::process::id()));
    let src = root.join("src");
    fs::create_dir_all(&src).expect("scratch dirs");
    fs::write(
        root.join("Cargo.toml"),
        format!(
            "[package]\nname = \"alloc-fixture\"\nversion = \"0.0.0\"\nedition = \"2024\"\n\n\
             [dependencies]\nfs-alloc = {{ path = {alloc_dir:?} }}\n\n[workspace]\n"
        ),
    )
    .expect("scratch manifest");

    for (case, source, expected) in CASES {
        // Fixtures are statement-level misuse snippets; they must live in a
        // function body for rustc to reach the intended borrow/move/private
        // diagnostics instead of item-grammar errors.
        let wrapped = format!("fn main() {{\n{source}}}");
        fs::write(src.join("lib.rs"), wrapped).expect("fixture source");
        let out = Command::new("cargo")
            .args(["check", "--offline", "--quiet"])
            .current_dir(&root)
            .env("RCH_DISABLE", "1")
            .env("CARGO_TARGET_DIR", root.join("target"))
            .output()
            .expect("cargo check runs");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(
            !out.status.success(),
            "fixture `{case}` compiled but must fail:\n{source}"
        );
        assert!(
            stderr.contains(expected),
            "fixture `{case}`: expected diagnostic containing {expected:?}, got:\n{stderr}"
        );
        println!(
            "{{\"suite\":\"fs-alloc\",\"case\":\"compile-fail-{case}\",\"verdict\":\"pass\",\
             \"detail\":\"rejected with expected diagnostic\"}}"
        );
    }
}
