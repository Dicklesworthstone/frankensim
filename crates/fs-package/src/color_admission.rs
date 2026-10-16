//! Package-side scientific-color admission (bead 6pf9, stage S2).
//!
//! [`crate::AdmittedClaim::scientific_color`] exposes a raw
//! [`Color`] that any downstream consumer can copy and re-launder as if it
//! carried admission authority. This module is the package-side twin of the
//! ledger's `ColorGraph` admission oracle: the authority is a
//! [`VerifiedPackage`] — a package inseparably paired with the checker
//! report/receipt that admitted it — and admission converts a
//! package-scientific claim into an [`fs_evidence::AdmittedColor`] whose
//! lineage names the exact claim.
//!
//! The receipt's node identity is the claim's domain-separated declaration
//! hash ([`crate::Claim::declared_content_hash_unverified`], the
//! `fs-package:v8:claim` surface), its row schema version is
//! [`crate::CLAIM_DECLARATION_IDENTITY_VERSION`], and its policy fingerprint
//! is [`package_color_admission_policy_fingerprint`] — deliberately DISTINCT
//! from the ledger authority's fingerprint so package-admitted and
//! ledger-admitted colors stay separately auditable and a receipt minted by
//! one authority can never satisfy the other's verifier.
//!
//! Minting refuses waiver-tainted claims (transitively, via the receipt's
//! admission class), non-positive ranks, and any retained package/report
//! pair whose binding no longer re-derives; the verifier re-derives every
//! acceptance from the retained pair, so a receipt outlives neither
//! tampering nor substitution.

use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::{
    AdmissionDecision, AdmissionReceipt, AdmissionVerifier, COLOR_ALGEBRA_VERSION, Color, ColorRank,
};

use crate::{CLAIM_DECLARATION_IDENTITY_VERSION, VerifiedPackage};

/// Fingerprint of the package claim-admission policy. Distinct from the
/// ledger color-admission policy by construction: consumers that require
/// ledger-grade admission can filter on the policy identity.
#[must_use]
pub fn package_color_admission_policy_fingerprint() -> ContentHash {
    hash_bytes(
        format!(
            "fs-package/claim-admission-policy/v1/claim-schema={CLAIM_DECLARATION_IDENTITY_VERSION}/algebra={COLOR_ALGEBRA_VERSION}"
        )
        .as_bytes(),
    )
}

/// Why a package claim could not mint an admission receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PackageColorAdmissionRefusal {
    /// No claim in the verified package carries this identity.
    UnknownClaim {
        /// The requested claim id.
        id: String,
    },
    /// The claim is directly waived or derived from a waiver-dependent
    /// claim; waiver-tainted evidence never acquires scientific admission.
    WaiverTainted {
        /// The refusing claim id.
        id: String,
    },
    /// Only positive ranks (Verified/Validated) carry scientific admission.
    NotPositive {
        /// The refusing claim id.
        id: String,
        /// The claim's declared rank.
        rank: ColorRank,
    },
    /// The retained package/report pair no longer re-derives: admissions,
    /// roots, or summaries diverge, so nothing about it can be admitted.
    BindingDivergence,
}

impl core::fmt::Display for PackageColorAdmissionRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::UnknownClaim { id } => {
                write!(f, "admission refused: claim {id:?} not in the package")
            }
            Self::WaiverTainted { id } => write!(
                f,
                "admission refused: claim {id:?} depends on an authenticated waiver"
            ),
            Self::NotPositive { id, rank } => write!(
                f,
                "admission refused: claim {id:?} declares {rank:?}, not positive evidence"
            ),
            Self::BindingDivergence => write!(
                f,
                "admission refused: the retained package/report binding does not re-derive"
            ),
        }
    }
}

impl std::error::Error for PackageColorAdmissionRefusal {}

impl VerifiedPackage {
    /// Mint an admission receipt for one scientific claim (bead 6pf9).
    /// Receipts convert declared colors into [`fs_evidence::AdmittedColor`]
    /// through [`PackageColorAdmissionVerifier`]; minting refuses unknown,
    /// waiver-tainted, and non-positive claims, and refuses entirely when
    /// the retained package/report binding no longer re-derives.
    ///
    /// # Errors
    /// [`PackageColorAdmissionRefusal`] naming the refusing gate.
    pub fn claim_admission_receipt(
        &self,
        claim_id: &str,
    ) -> Result<AdmissionReceipt, PackageColorAdmissionRefusal> {
        if !self.validate_binding() {
            return Err(PackageColorAdmissionRefusal::BindingDivergence);
        }
        let admitted = self
            .admitted_claims()
            .find(|claim| claim.id() == claim_id)
            .ok_or_else(|| PackageColorAdmissionRefusal::UnknownClaim {
                id: claim_id.to_string(),
            })?;
        let color = admitted.scientific_color().ok_or_else(|| {
            PackageColorAdmissionRefusal::WaiverTainted {
                id: claim_id.to_string(),
            }
        })?;
        let rank = color.rank();
        if rank == ColorRank::Estimated {
            return Err(PackageColorAdmissionRefusal::NotPositive {
                id: claim_id.to_string(),
                rank,
            });
        }
        Ok(AdmissionReceipt::from_parts(
            admitted.claim_declaration_hash(),
            CLAIM_DECLARATION_IDENTITY_VERSION,
            COLOR_ALGEBRA_VERSION,
            package_color_admission_policy_fingerprint(),
        ))
    }
}

/// The package-side admission oracle (bead 6pf9): authenticates a
/// (candidate, receipt) pair by re-deriving it from the retained
/// package/report pair. Acceptance requires the receipt's node hash to name
/// a claim in this package whose admission class is scientific, whose
/// declared color is bit-exactly the candidate (canonical bytes, not
/// display JSON), whose receipt versions match this build, whose policy
/// fingerprint is this authority's own, and whose package/report binding
/// still re-derives.
#[derive(Debug, Clone, Copy)]
pub struct PackageColorAdmissionVerifier<'p> {
    verified: &'p VerifiedPackage,
}

impl<'p> PackageColorAdmissionVerifier<'p> {
    /// Wrap the verified package that will re-derive admissions.
    #[must_use]
    pub fn new(verified: &'p VerifiedPackage) -> Self {
        PackageColorAdmissionVerifier { verified }
    }
}

impl AdmissionVerifier for PackageColorAdmissionVerifier<'_> {
    fn verify(&self, candidate: &Color, receipt: &AdmissionReceipt) -> AdmissionDecision {
        let policy = package_color_admission_policy_fingerprint();
        if receipt.row_schema_version() != CLAIM_DECLARATION_IDENTITY_VERSION
            || receipt.color_algebra_version() != COLOR_ALGEBRA_VERSION
            || receipt.policy_fingerprint() != policy
        {
            return AdmissionDecision::reject(policy);
        }
        let Some(admitted) = self
            .verified
            .admitted_claims()
            .find(|claim| claim.claim_declaration_hash() == receipt.node_hash())
        else {
            return AdmissionDecision::reject(policy);
        };
        let Some(scientific) = admitted.scientific_color() else {
            return AdmissionDecision::reject(policy);
        };
        if scientific.canonical_bytes() != candidate.canonical_bytes() {
            return AdmissionDecision::reject(policy);
        }
        if !self.verified.validate_binding() {
            return AdmissionDecision::reject(policy);
        }
        AdmissionDecision::accept(policy)
    }
}
