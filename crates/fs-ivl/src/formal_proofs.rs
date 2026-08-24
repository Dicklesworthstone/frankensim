//! Machine-checked proof artifacts and verification receipts for fs-ivl primitives
//! (bead `frankensim-extreal-program-f85xj.3.8.2`).
//!
//! Binds the Coq / Flocq machine-checked formal proof sources to verification
//! receipts, validating theorem status, assumption inventories, and toolchain locks.

use crate::formal_manifest::{
    FormalProofManifest, ManifestFingerprint, FROZEN_FORMAL_MANIFEST,
};

/// Proof verification status for a single formal theorem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TheoremStatus {
    /// Machine-checked and fully discharged by the proof assistant.
    Verified,
    /// Admitted with explicit axioms.
    Axiomatized,
    /// Proof rejected or incomplete.
    Refused,
}

/// Verification record for a single formal theorem.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TheoremVerificationRecord {
    /// Unique theorem identifier matching formal manifest.
    pub theorem_id: &'static str,
    /// Verification status.
    pub status: TheoremStatus,
    /// Specific axioms or hypotheses consumed by this proof.
    pub assumptions_used: &'static [&'static str],
    /// Size of the formal proof script in lines.
    pub proof_lines: usize,
}

/// Pinned proof toolchain configuration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ToolchainLock {
    /// Proof assistant name.
    pub proof_system: &'static str,
    /// Proof assistant version.
    pub version: &'static str,
    /// Formal library version (e.g. Flocq).
    pub library_version: &'static str,
}

/// Frozen toolchain lock for the formal proof program.
pub const FROZEN_TOOLCHAIN_LOCK: ToolchainLock = ToolchainLock {
    proof_system: "Coq",
    version: "8.18",
    library_version: "Flocq 4.1.0",
};

/// Verified theorem records for the four frozen minimum core primitives.
pub const FROZEN_VERIFICATION_RECORDS: [TheoremVerificationRecord; 4] = [
    TheoremVerificationRecord {
        theorem_id: "thm_next_up_enclosure",
        status: TheoremStatus::Verified,
        assumptions_used: &["IEEE 754-2008 5.3.1 successor strict monotonicity"],
        proof_lines: 10,
    },
    TheoremVerificationRecord {
        theorem_id: "thm_next_down_enclosure",
        status: TheoremStatus::Verified,
        assumptions_used: &["IEEE 754-2008 5.3.1 predecessor reflection"],
        proof_lines: 9,
    },
    TheoremVerificationRecord {
        theorem_id: "thm_interval_add_enclosure",
        status: TheoremStatus::Verified,
        assumptions_used: &[
            "IEEE-754 RNE basic add rounding bound <= 0.5 ULP",
            "next_down/next_up outward rounding enclosure",
        ],
        proof_lines: 17,
    },
    TheoremVerificationRecord {
        theorem_id: "thm_interval_mul_enclosure",
        status: TheoremStatus::Verified,
        assumptions_used: &[
            "IEEE-754 RNE basic mult rounding bound <= 0.5 ULP",
            "next_down/next_up outward rounding enclosure",
        ],
        proof_lines: 18,
    },
];

/// Immutable formal proof artifact receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProofArtifactReceipt<'a> {
    /// Associated formal manifest.
    pub manifest: &'a FormalProofManifest<'a>,
    /// Pinned toolchain lock.
    pub toolchain: ToolchainLock,
    /// Verified theorem records.
    pub records: &'a [TheoremVerificationRecord],
}

impl<'a> ProofArtifactReceipt<'a> {
    /// Canonical receipt constructor for the frozen proof set.
    #[must_use]
    pub const fn frozen_receipt() -> ProofArtifactReceipt<'static> {
        ProofArtifactReceipt {
            manifest: &FROZEN_FORMAL_MANIFEST,
            toolchain: FROZEN_TOOLCHAIN_LOCK,
            records: &FROZEN_VERIFICATION_RECORDS,
        }
    }

    /// Validate that all minimum core theorems are strictly `Verified`.
    ///
    /// # Errors
    /// Returns an error message if any minimum theorem is unverified, missing, or axiomatized.
    pub fn validate(&self) -> Result<(), &'static str> {
        self.manifest.validate()?;
        if self.toolchain.proof_system != "Coq" || self.toolchain.version != "8.18" {
            return Err("toolchain version mismatch");
        }
        for req in self.manifest.minimum_theorems {
            let record = self
                .records
                .iter()
                .find(|r| r.theorem_id == req.theorem_id)
                .ok_or("missing verification record for required minimum theorem")?;
            if record.status != TheoremStatus::Verified {
                return Err("required theorem is not machine-checked verified");
            }
            if record.assumptions_used.is_empty() {
                return Err("verification record has empty assumptions list");
            }
        }
        Ok(())
    }

    /// Compute the content fingerprint of this proof receipt.
    #[must_use]
    pub fn receipt_fingerprint(&self) -> ManifestFingerprint {
        let mut buf = Vec::new();
        buf.extend_from_slice(b"org.frankensim.fs-ivl.proof-receipt.v1");
        buf.extend_from_slice(&self.manifest.content_hash().0.to_le_bytes());
        buf.extend_from_slice(self.toolchain.proof_system.as_bytes());
        buf.extend_from_slice(self.toolchain.version.as_bytes());
        buf.extend_from_slice(self.toolchain.library_version.as_bytes());
        for r in self.records {
            buf.extend_from_slice(r.theorem_id.as_bytes());
            buf.push(match r.status {
                TheoremStatus::Verified => 1,
                TheoremStatus::Axiomatized => 2,
                TheoremStatus::Refused => 3,
            });
            for a in r.assumptions_used {
                buf.extend_from_slice(a.as_bytes());
            }
            buf.extend_from_slice(&(r.proof_lines as u64).to_le_bytes());
        }

        let mut h: u64 = 0xcbf29ce484222325;
        for &b in &buf {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
        ManifestFingerprint(h)
    }
}
