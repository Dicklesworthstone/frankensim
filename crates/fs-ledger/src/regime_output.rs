//! Exact artifact retention for final operating-envelope demotion receipts.
//!
//! `fs-regime` owns the receipt schema and semantic identity. This module only
//! gives those already-typed receipts a bounded, content-addressed ledger
//! representation and an exact caller-pinned read path. It does not parse
//! arbitrary artifact bytes into scientific authority.

use fs_regime::OutputClaimReceipt;

use crate::{ContentHash, Ledger, LedgerError, PutReceipt};

/// Artifact kind reserved for canonical `fs-regime` demotion receipts.
pub const REGIME_DEMOTION_RECEIPT_ARTIFACT_KIND: &str = "regime-output-demotion-receipt-v1";

/// Maximum canonical receipt bytes accepted by the ledger adapter.
pub const MAX_REGIME_DEMOTION_RECEIPT_BYTES: usize = 1024 * 1024;
const MAX_REGIME_DEMOTION_RECEIPT_BYTES_U64: u64 = 1024 * 1024;

/// The two distinct identities returned after retaining one demotion receipt.
///
/// `artifact` addresses the exact stored bytes under the ledger's artifact
/// identity. `receipt_id` is the `fs-regime` semantic receipt identity. Neither
/// is an authentication result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RegimeDemotionArtifactReceipt {
    /// Ledger artifact write result for the exact canonical JSON bytes.
    pub artifact: PutReceipt,
    /// Domain-separated `fs-regime` identity of the typed receipt.
    pub receipt_id: fs_blake3::ContentHash,
}

impl Ledger {
    /// Retain one already-demoted output receipt as exact canonical JSON.
    ///
    /// The write is content-addressed and idempotent under the ledger's normal
    /// agreeing-envelope rule. Fully in-domain receipts are refused because
    /// this artifact kind is specifically the no-claim/demotion trail.
    ///
    /// # Errors
    ///
    /// Returns [`LedgerError::Invalid`] for an in-domain receipt or a canonical
    /// payload above [`MAX_REGIME_DEMOTION_RECEIPT_BYTES`]. Storage errors and
    /// envelope conflicts are forwarded from [`Ledger::put_artifact`].
    pub fn put_regime_demotion_receipt(
        &self,
        receipt: &OutputClaimReceipt,
    ) -> Result<RegimeDemotionArtifactReceipt, LedgerError> {
        if !receipt.demoted() {
            return Err(LedgerError::Invalid {
                field: "regime_output_receipt.coverage".to_string(),
                problem: "receipt is fully in-domain; only demotion receipts belong in the no-claim artifact stream"
                    .to_string(),
            });
        }

        let canonical = receipt.to_canonical_json();
        if canonical.len() > MAX_REGIME_DEMOTION_RECEIPT_BYTES {
            return Err(LedgerError::Invalid {
                field: "regime_output_receipt.canonical_json".to_string(),
                problem: format!(
                    "{} bytes exceeds the {}-byte demotion-receipt limit",
                    canonical.len(),
                    MAX_REGIME_DEMOTION_RECEIPT_BYTES
                ),
            });
        }

        let artifact = self.put_artifact(
            REGIME_DEMOTION_RECEIPT_ARTIFACT_KIND,
            canonical.as_bytes(),
            None,
        )?;
        Ok(RegimeDemotionArtifactReceipt {
            artifact,
            receipt_id: receipt.content_id(),
        })
    }

    /// Read a retained demotion receipt only under an exact typed expectation.
    ///
    /// The caller supplies the `OutputClaimReceipt` it expects. This method
    /// preflights the artifact kind and byte limit, materializes within the
    /// fixed cap, and returns bytes only if they equal that receipt's canonical
    /// encoding exactly. It intentionally offers no unpinned parse path.
    ///
    /// # Errors
    ///
    /// Returns [`LedgerError::NotFound`] when `artifact` is absent,
    /// [`LedgerError::Invalid`] for a wrong artifact kind or exact-receipt
    /// mismatch, and the normal bounded artifact-read errors for corrupt or
    /// oversized storage.
    pub fn read_regime_demotion_receipt(
        &self,
        artifact: &ContentHash,
        expected: &OutputClaimReceipt,
    ) -> Result<Vec<u8>, LedgerError> {
        if !expected.demoted() {
            return Err(LedgerError::Invalid {
                field: "regime_output_receipt.coverage".to_string(),
                problem: "caller-pinned receipt is fully in-domain; the demotion artifact stream cannot admit it"
                    .to_string(),
            });
        }
        let info = self
            .artifact_info(artifact)?
            .ok_or_else(|| LedgerError::NotFound {
                what: format!("regime demotion receipt artifact {}", artifact.to_hex()),
            })?;
        if info.kind != REGIME_DEMOTION_RECEIPT_ARTIFACT_KIND {
            return Err(LedgerError::Invalid {
                field: "artifact.kind".to_string(),
                problem: format!(
                    "expected {REGIME_DEMOTION_RECEIPT_ARTIFACT_KIND:?}, found {:?}",
                    info.kind
                ),
            });
        }

        let bytes = self
            .get_artifact_bounded(artifact, MAX_REGIME_DEMOTION_RECEIPT_BYTES_U64)?
            .ok_or_else(|| LedgerError::NotFound {
                what: format!("regime demotion receipt artifact {}", artifact.to_hex()),
            })?;
        let expected_bytes = expected.to_canonical_json();
        if bytes != expected_bytes.as_bytes() {
            return Err(LedgerError::Invalid {
                field: "regime_output_receipt.expected".to_string(),
                problem: format!(
                    "artifact {} does not equal caller-pinned fs-regime receipt {}",
                    artifact.to_hex(),
                    expected.content_id().to_hex()
                ),
            });
        }
        Ok(bytes)
    }
}
