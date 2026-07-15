//! Opaque admitted scientific color (bead 6pf9, stage S1).
//!
//! [`crate::color::Color`] is a plain DECLARATION enum: any crate can write
//! `Color::Verified { .. }` and structural payload validation cannot tell an
//! admitted certificate from a fabricated literal. This module splits the
//! epistemic-rank REPRESENTATION from admission AUTHORITY: an
//! [`AdmittedColor`] has private fields and exactly one constructor,
//! [`AdmittedColor::from_receipt`], which demands three things —
//!
//! 1. a structurally valid, POSITIVE (Verified/Validated) candidate color;
//! 2. an [`AdmissionReceipt`] carrying the color-write row identity that
//!    admitted it (node provenance hash, row schema version, color-algebra
//!    version, admitting policy fingerprint);
//! 3. an [`AdmissionVerifier`] capability that authenticates the pair.
//!
//! Receipts are plain data — anyone can build one — because authority lives
//! in the VERIFIER, exactly like the waiver and source-origin capabilities:
//! the default [`NoAdmissionVerifier`] refuses everything, so at this layer
//! NOTHING admits. The real verifier is injected by the admission authority
//! (HELM-side `ColorGraph`, which re-derives the receipt from its replay-
//! audited node state); a lying verifier is visible at the composition root,
//! the same trust model as a lying `WaiverVerifier`. This crate stays at
//! UTIL layer: it never depends on the ledger, it only names the shape the
//! ledger must authenticate.
//!
//! Defense in depth: even a verifier that accepts everything cannot mint an
//! [`AdmittedColor`] from a malformed payload, a non-positive rank, a
//! stale-algebra receipt, or an accepting decision whose policy identity does
//! not match the receipt — those refusals fire HERE around the capability
//! boundary.

use fs_blake3::ContentHash;

use crate::color::{
    COLOR_ALGEBRA_VERSION, Color, ColorPayloadError, ColorRank, validate_color_payload,
};

/// The color-write row identity under which a color was admitted. Plain
/// data: constructible by anyone, authoritative only once an
/// [`AdmissionVerifier`] authenticates it against the admission ledger.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdmissionReceipt {
    node_hash: ContentHash,
    row_schema_version: u32,
    color_algebra_version: u32,
    policy_fingerprint: ContentHash,
}

impl AdmissionReceipt {
    /// Assemble a receipt from its row identity parts.
    #[must_use]
    pub fn from_parts(
        node_hash: ContentHash,
        row_schema_version: u32,
        color_algebra_version: u32,
        policy_fingerprint: ContentHash,
    ) -> Self {
        AdmissionReceipt {
            node_hash,
            row_schema_version,
            color_algebra_version,
            policy_fingerprint,
        }
    }

    /// Provenance hash of the admitted color-graph node.
    #[must_use]
    pub fn node_hash(&self) -> ContentHash {
        self.node_hash
    }

    /// Color-write row schema version at admission time.
    #[must_use]
    pub fn row_schema_version(&self) -> u32 {
        self.row_schema_version
    }

    /// Color-algebra version bound at admission time.
    #[must_use]
    pub fn color_algebra_version(&self) -> u32 {
        self.color_algebra_version
    }

    /// Fingerprint of the policy that admitted the node.
    #[must_use]
    pub fn policy_fingerprint(&self) -> ContentHash {
        self.policy_fingerprint
    }
}

/// A verifier's answer, binding the policy identity that decided.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AdmissionDecision {
    accepted: bool,
    policy: ContentHash,
}

impl AdmissionDecision {
    /// Accept under the given policy identity.
    #[must_use]
    pub fn accept(policy: ContentHash) -> Self {
        AdmissionDecision {
            accepted: true,
            policy,
        }
    }

    /// Refuse under the given policy identity.
    #[must_use]
    pub fn reject(policy: ContentHash) -> Self {
        AdmissionDecision {
            accepted: false,
            policy,
        }
    }

    /// Whether the verifier accepted.
    #[must_use]
    pub fn accepted(&self) -> bool {
        self.accepted
    }

    /// The deciding policy identity.
    #[must_use]
    pub fn policy(&self) -> ContentHash {
        self.policy
    }
}

/// Injected capability that authenticates a (candidate, receipt) pair
/// against the admission ledger. Implemented by the HELM-side authority;
/// everything below it uses [`NoAdmissionVerifier`].
pub trait AdmissionVerifier {
    /// Authenticate one candidate against one receipt.
    fn verify(&self, candidate: &Color, receipt: &AdmissionReceipt) -> AdmissionDecision;
}

/// Deny-all default: at UTIL layer nothing admits scientific evidence.
#[derive(Debug, Clone, Copy, Default)]
pub struct NoAdmissionVerifier;

/// Stable policy identity for the deny-all default.
#[must_use]
pub fn no_admission_policy() -> ContentHash {
    fs_blake3::hash_bytes(b"fs-evidence/no-admission-verifier/deny-all/v1")
}

impl AdmissionVerifier for NoAdmissionVerifier {
    fn verify(&self, _candidate: &Color, _receipt: &AdmissionReceipt) -> AdmissionDecision {
        AdmissionDecision::reject(no_admission_policy())
    }
}

/// Why a candidate could not become admitted scientific evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AdmissionRejection {
    /// The candidate payload is structurally malformed. No verifier can
    /// override this.
    MalformedPayload(ColorPayloadError),
    /// Only positive ranks (Verified/Validated) carry scientific admission;
    /// Estimated evidence stays a declared candidate.
    NotPositive {
        /// The candidate's rank.
        rank: ColorRank,
    },
    /// The receipt was minted under a different color algebra than this
    /// build composes with. Stale-algebra evidence cannot convert.
    StaleAlgebra {
        /// Algebra version bound in the receipt.
        receipt: u32,
        /// Algebra version of this build.
        current: u32,
    },
    /// The verifier refused the (candidate, receipt) pair.
    Refused {
        /// The deciding policy identity.
        policy: ContentHash,
    },
    /// An accepting verifier named a policy other than the policy bound into
    /// the receipt. Acceptance cannot silently rewrite or detach lineage.
    PolicyMismatch {
        /// Policy identity committed by the receipt.
        receipt: ContentHash,
        /// Policy identity returned by the accepting decision.
        decision: ContentHash,
    },
}

impl core::fmt::Display for AdmissionRejection {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MalformedPayload(error) => {
                write!(f, "admission refused: malformed candidate payload: {error}")
            }
            Self::NotPositive { rank } => write!(
                f,
                "admission refused: {rank:?} is not positive scientific evidence"
            ),
            Self::StaleAlgebra { receipt, current } => write!(
                f,
                "admission refused: receipt bound color-algebra v{receipt}, this build composes v{current}"
            ),
            Self::Refused { policy } => {
                write!(f, "admission refused by policy {policy}")
            }
            Self::PolicyMismatch { receipt, decision } => write!(
                f,
                "admission refused: receipt policy {receipt} differs from accepting decision policy {decision}"
            ),
        }
    }
}

impl std::error::Error for AdmissionRejection {}

/// Positive scientific evidence with an authenticated admission lineage.
/// Fields are private and the only constructor is [`Self::from_receipt`]:
/// holding a value of this type MEANS an [`AdmissionVerifier`] accepted the
/// (candidate, receipt) pair after the local structural gates passed.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedColor {
    color: Color,
    receipt: AdmissionReceipt,
}

impl AdmittedColor {
    /// Admit a candidate color under a receipt, authenticated by `verifier`.
    ///
    /// Local gates fire before the capability is consulted, so even an
    /// accept-everything verifier cannot admit a malformed payload, a
    /// non-positive rank, or a stale-algebra receipt.
    ///
    /// # Errors
    /// [`AdmissionRejection`] naming the exact refusing gate.
    pub fn from_receipt(
        color: Color,
        receipt: AdmissionReceipt,
        verifier: &dyn AdmissionVerifier,
    ) -> Result<Self, AdmissionRejection> {
        validate_color_payload(&color).map_err(AdmissionRejection::MalformedPayload)?;
        let rank = color.rank();
        if rank == ColorRank::Estimated {
            return Err(AdmissionRejection::NotPositive { rank });
        }
        if receipt.color_algebra_version != COLOR_ALGEBRA_VERSION {
            return Err(AdmissionRejection::StaleAlgebra {
                receipt: receipt.color_algebra_version,
                current: COLOR_ALGEBRA_VERSION,
            });
        }
        let decision = verifier.verify(&color, &receipt);
        if !decision.accepted() {
            return Err(AdmissionRejection::Refused {
                policy: decision.policy(),
            });
        }
        if decision.policy() != receipt.policy_fingerprint {
            return Err(AdmissionRejection::PolicyMismatch {
                receipt: receipt.policy_fingerprint,
                decision: decision.policy(),
            });
        }
        Ok(AdmittedColor { color, receipt })
    }

    /// The admitted color. Named to read as evidence, not declaration.
    #[must_use]
    pub fn admitted_color(&self) -> &Color {
        &self.color
    }

    /// The admitted rank (always Verified or Validated).
    #[must_use]
    pub fn rank(&self) -> ColorRank {
        self.color.rank()
    }

    /// The receipt this value was admitted under.
    #[must_use]
    pub fn receipt(&self) -> &AdmissionReceipt {
        &self.receipt
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn receipt(algebra: u32) -> AdmissionReceipt {
        AdmissionReceipt::from_parts(
            fs_blake3::hash_bytes(b"test-node"),
            7,
            algebra,
            fs_blake3::hash_bytes(b"test-policy"),
        )
    }

    struct LyingVerifier;
    impl AdmissionVerifier for LyingVerifier {
        fn verify(&self, _c: &Color, _r: &AdmissionReceipt) -> AdmissionDecision {
            AdmissionDecision::accept(fs_blake3::hash_bytes(b"lying-policy"))
        }
    }

    struct ReceiptPolicyVerifier;
    impl AdmissionVerifier for ReceiptPolicyVerifier {
        fn verify(&self, _c: &Color, receipt: &AdmissionReceipt) -> AdmissionDecision {
            AdmissionDecision::accept(receipt.policy_fingerprint())
        }
    }

    #[test]
    fn deny_all_default_refuses_a_well_formed_candidate() {
        let error = AdmittedColor::from_receipt(
            Color::Verified { lo: 0.0, hi: 1.0 },
            receipt(COLOR_ALGEBRA_VERSION),
            &NoAdmissionVerifier,
        )
        .expect_err("deny-all must refuse");
        assert_eq!(
            error,
            AdmissionRejection::Refused {
                policy: no_admission_policy()
            }
        );
    }

    #[test]
    fn local_gates_fire_before_the_capability_even_for_a_lying_verifier() {
        // Malformed payload: inverted interval.
        let malformed = AdmittedColor::from_receipt(
            Color::Verified { lo: 1.0, hi: 0.0 },
            receipt(COLOR_ALGEBRA_VERSION),
            &LyingVerifier,
        )
        .expect_err("inverted interval must refuse");
        assert!(matches!(malformed, AdmissionRejection::MalformedPayload(_)));

        // Non-positive rank.
        let estimated = AdmittedColor::from_receipt(
            Color::Estimated {
                estimator: "probe-a".to_string(),
                dispersion: 0.5,
            },
            receipt(COLOR_ALGEBRA_VERSION),
            &LyingVerifier,
        )
        .expect_err("estimated must refuse");
        assert_eq!(
            estimated,
            AdmissionRejection::NotPositive {
                rank: ColorRank::Estimated
            }
        );

        // Stale algebra.
        let stale = AdmittedColor::from_receipt(
            Color::Verified { lo: 0.0, hi: 1.0 },
            receipt(COLOR_ALGEBRA_VERSION - 1),
            &LyingVerifier,
        )
        .expect_err("stale algebra must refuse");
        assert_eq!(
            stale,
            AdmissionRejection::StaleAlgebra {
                receipt: COLOR_ALGEBRA_VERSION - 1,
                current: COLOR_ALGEBRA_VERSION
            }
        );
    }

    #[test]
    fn an_accepting_verifier_mints_and_the_value_reports_its_lineage() {
        let admitted = AdmittedColor::from_receipt(
            Color::Verified { lo: -2.0, hi: 3.0 },
            receipt(COLOR_ALGEBRA_VERSION),
            &ReceiptPolicyVerifier,
        )
        .expect("accepting verifier mints");
        assert_eq!(admitted.rank(), ColorRank::Verified);
        assert_eq!(
            admitted.admitted_color(),
            &Color::Verified { lo: -2.0, hi: 3.0 }
        );
        assert_eq!(
            admitted.receipt().node_hash(),
            fs_blake3::hash_bytes(b"test-node")
        );
    }

    #[test]
    fn accepting_decision_must_match_the_receipt_policy() {
        let receipt = receipt(COLOR_ALGEBRA_VERSION);
        let error = AdmittedColor::from_receipt(
            Color::Verified { lo: -2.0, hi: 3.0 },
            receipt,
            &LyingVerifier,
        )
        .expect_err("an accepting decision cannot substitute another policy identity");
        assert_eq!(
            error,
            AdmissionRejection::PolicyMismatch {
                receipt: receipt.policy_fingerprint(),
                decision: fs_blake3::hash_bytes(b"lying-policy"),
            }
        );
    }
}
