//! Package projection for a typed `NoUsefulBound` refusal.
//!
//! This record deliberately sits beside, rather than inside, the v8
//! `EvidencePackage` claim list. A failed usefulness criterion is not a
//! scientific color and must never be converted into a certificate claim.

use core::fmt;

use fs_evidence::NoUsefulBound;

/// Maximum UTF-8 bytes in a package-local refusal claim id.
pub const MAX_NO_USEFUL_BOUND_CLAIM_ID_BYTES: usize = 512;

/// Machine-readable package projection of one useful-bound refusal.
#[derive(Debug, Clone, PartialEq)]
pub struct NoUsefulBoundRecord {
    claim_id: String,
    refusal: NoUsefulBound,
}

impl NoUsefulBoundRecord {
    /// Construct a package-side refusal record.
    pub fn try_new(
        claim_id: impl Into<String>,
        refusal: NoUsefulBound,
    ) -> Result<Self, NoUsefulBoundRecordError> {
        let claim_id = claim_id.into();
        if claim_id.trim().is_empty()
            || claim_id.len() > MAX_NO_USEFUL_BOUND_CLAIM_ID_BYTES
            || claim_id.chars().any(char::is_control)
        {
            return Err(NoUsefulBoundRecordError::InvalidClaimId {
                bytes: claim_id.len(),
            });
        }
        Ok(Self { claim_id, refusal })
    }

    /// Stable claim id of the decision that could not be bounded usefully.
    #[must_use]
    pub fn claim_id(&self) -> &str {
        &self.claim_id
    }

    /// Exact typed refusal.
    #[must_use]
    pub const fn refusal(&self) -> &NoUsefulBound {
        &self.refusal
    }

    /// This projection never supplies an `EvidencePackage` certificate claim.
    #[must_use]
    pub const fn has_certificate_claim(&self) -> bool {
        false
    }

    /// Deterministic package-inspection rendering with an explicit no-claim
    /// boundary.
    #[must_use]
    pub fn render_manifest(&self) -> String {
        format!(
            "package-record=no-useful-bound\nclaim-id={}\n{}certificate-claim=none\nno-claim=valid enclosure retained; engineering compliance and scientific color not established\n",
            self.claim_id,
            self.refusal.render_report(),
        )
    }
}

/// Refusal while constructing a package-side `NoUsefulBound` record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoUsefulBoundRecordError {
    /// Claim id was empty, oversized, or contained a control character.
    InvalidClaimId {
        /// Presented UTF-8 byte count.
        bytes: usize,
    },
}

impl fmt::Display for NoUsefulBoundRecordError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidClaimId { bytes } => write!(
                f,
                "NoUsefulBound claim id must contain 1..={MAX_NO_USEFUL_BOUND_CLAIM_ID_BYTES} non-control UTF-8 bytes; found {bytes}"
            ),
        }
    }
}

impl std::error::Error for NoUsefulBoundRecordError {}
