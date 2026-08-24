//! Formal proofs test suite (bead `frankensim-extreal-program-f85xj.3.8.2`).
//!
//! Validates:
//! - Verification status of all 4 minimum core primitives
//! - Assumption inventory and toolchain lock
//! - Deterministic receipt fingerprinting
//! - Mutant refusal tests (unverified theorem, bad toolchain, empty assumptions)

use fs_ivl::formal_manifest::FROZEN_FORMAL_MANIFEST;
use fs_ivl::formal_proofs::{
    ProofArtifactReceipt, TheoremStatus, ToolchainLock, FROZEN_TOOLCHAIN_LOCK,
    FROZEN_VERIFICATION_RECORDS,
};

#[test]
fn frozen_proof_receipt_validates_all_core_theorems() {
    let receipt = ProofArtifactReceipt::frozen_receipt();
    assert_eq!(receipt.validate(), Ok(()));
    assert_eq!(receipt.records.len(), 4);
    for r in receipt.records {
        assert_eq!(r.status, TheoremStatus::Verified);
        assert!(!r.assumptions_used.is_empty());
        assert!(r.proof_lines > 0);
    }
}

#[test]
fn receipt_fingerprint_is_deterministic_and_non_empty() {
    let receipt = ProofArtifactReceipt::frozen_receipt();
    let fp1 = receipt.receipt_fingerprint();
    let fp2 = receipt.receipt_fingerprint();
    assert_eq!(fp1, fp2);
    assert!(!fp1.to_hex().is_empty());
}

#[test]
fn unverified_theorem_refuses_validation() {
    let mut bad_records = FROZEN_VERIFICATION_RECORDS.clone();
    bad_records[0].status = TheoremStatus::Refused;

    let receipt = ProofArtifactReceipt {
        manifest: &FROZEN_FORMAL_MANIFEST,
        toolchain: FROZEN_TOOLCHAIN_LOCK,
        records: &bad_records,
    };
    assert_eq!(
        receipt.validate(),
        Err("required theorem is not machine-checked verified")
    );
}

#[test]
fn bad_toolchain_version_refuses_validation() {
    let bad_toolchain = ToolchainLock {
        proof_system: "Coq",
        version: "8.14", // old version
        library_version: "Flocq 3.4.0",
    };

    let receipt = ProofArtifactReceipt {
        manifest: &FROZEN_FORMAL_MANIFEST,
        toolchain: bad_toolchain,
        records: &FROZEN_VERIFICATION_RECORDS,
    };
    assert_eq!(receipt.validate(), Err("toolchain version mismatch"));
}

#[test]
fn empty_assumptions_refuse_validation() {
    let mut bad_records = FROZEN_VERIFICATION_RECORDS.clone();
    bad_records[0].assumptions_used = &[];

    let receipt = ProofArtifactReceipt {
        manifest: &FROZEN_FORMAL_MANIFEST,
        toolchain: FROZEN_TOOLCHAIN_LOCK,
        records: &bad_records,
    };
    assert_eq!(
        receipt.validate(),
        Err("verification record has empty assumptions list")
    );
}
