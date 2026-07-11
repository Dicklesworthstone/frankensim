//! CLAIM ORIGINS (schema v5, bead krym): where a claim's certificate
//! CAME FROM — the missing half of "machine-checkable".
//!
//! Schema v4 made the content address collision-resistant, but content
//! consistency is not evidence origin: `Color::Verified { lo, hi }` is
//! public algebra (fs-evidence composition needs it), so any producer
//! could mint a finite interval and the standalone checker would pass
//! it. v5 seals claims AT THE PACKAGE BOUNDARY: every claim must carry
//! a [`ClaimOrigin`] consistent with its color, bound into the content
//! address, and re-derivable by the checker — while the Color algebra
//! itself stays public and untouched.
//!
//! The five origins and their re-derivation obligations:
//! - [`ClaimOrigin::SourceCertificate`] — a named producer plus the
//!   64-hex content hash of its certificate artifact (solver
//!   certificate, proof object). The checker verifies shape and the
//!   color class (Verified); the artifact hash makes the certificate
//!   subpoenable without shipping it.
//! - [`ClaimOrigin::AnchoredSource`] — a validated claim's reference
//!   dataset by id + content hash; must MATCH the color's named
//!   dataset exactly (an unrelated anchor is refused).
//! - [`ClaimOrigin::EstimatedSource`] — the estimator identity; must
//!   match the color's estimator string exactly.
//! - [`ClaimOrigin::Derived`] — a composition receipt; the checker
//!   re-runs `compose` over the parents and the result must equal the
//!   claimed color bit-exactly (the v3 receipt machinery, now the
//!   origin itself).
//! - [`ClaimOrigin::AuthenticatedWaiver`] — an explicit, expiring,
//!   MAC'd grant. NEVER self-authorizing: verification requires an
//!   INJECTED [`WaiverVerifier`] capability plus a date context; the
//!   in-tree default refuses everything (the fs-ledger fail-closed
//!   pattern). The MAC binds the claim's canonical bytes, so a waiver
//!   replayed onto a different claim fails.

use core::fmt;

/// An explicit waiver grant that travels WITH its claim (schema v5).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WaiverGrant {
    /// Stable, non-blank waiver identity (audit key).
    pub waiver_id: String,
    /// Last day (days since the Unix epoch) this waiver is valid.
    pub expiry_day: u64,
    /// Authenticator over the waiver id, expiry, and the CLAIM'S
    /// canonical bytes (replay onto another claim changes the message).
    /// Opaque here: only an injected [`WaiverVerifier`] can accept it.
    pub mac: String,
}

/// The waiver-verification CAPABILITY (injected; fs-package ships no
/// cryptography — the same fail-closed pattern as the checker's
/// [`SignatureVerifier`] and fs-ledger's waivers).
pub trait WaiverVerifier {
    /// True iff `grant.mac` authenticates `message` (the claim's
    /// canonical bytes plus the grant's own id/expiry) for this grant.
    fn verify(&self, grant: &WaiverGrant, message: &[u8]) -> bool;
}

/// The in-tree default: nothing authenticates. A package whose claims
/// carry waiver origins can NEVER verify without an explicitly
/// injected capability and date context.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoWaiverVerifier;

impl WaiverVerifier for NoWaiverVerifier {
    fn verify(&self, _grant: &WaiverGrant, _message: &[u8]) -> bool {
        false
    }
}

/// Where a claim's certificate came from (schema v5). Bound into the
/// content address and re-derived by the standalone checker.
#[derive(Debug, Clone, PartialEq)]
pub enum ClaimOrigin {
    /// A named producer's certificate artifact (64-hex content hash).
    SourceCertificate {
        /// Non-blank producer identity (e.g. "fs-solver/ivp-cert").
        producer: String,
        /// Canonical 64-hex lowercase content hash of the certificate.
        certificate_hash: String,
    },
    /// The validated color's reference dataset, by id + content hash.
    AnchoredSource {
        /// Must equal the color's named dataset exactly.
        dataset_id: String,
        /// Canonical 64-hex lowercase content hash of the dataset.
        content_hash: String,
    },
    /// The estimated color's estimator identity.
    EstimatedSource {
        /// Must equal the color's estimator string exactly.
        estimator: String,
    },
    /// Derived from earlier claims: the composition receipt IS the
    /// origin (parents by index, fold op) — re-run by the checker.
    Derived,
    /// An explicit, expiring, MAC'd waiver (see [`WaiverGrant`]).
    AuthenticatedWaiver(WaiverGrant),
}

impl ClaimOrigin {
    /// Stable kind tag for hashing, JSON, and refusal messages.
    #[must_use]
    pub const fn kind(&self) -> &'static str {
        match self {
            ClaimOrigin::SourceCertificate { .. } => "source-certificate",
            ClaimOrigin::AnchoredSource { .. } => "anchored-source",
            ClaimOrigin::EstimatedSource { .. } => "estimated-source",
            ClaimOrigin::Derived => "derived",
            ClaimOrigin::AuthenticatedWaiver(_) => "authenticated-waiver",
        }
    }

    /// The canonical atom sequence bound into the claim's content
    /// hash (length-prefixed strings via the caller's `push_atom`
    /// discipline; this returns the ordered raw parts).
    #[must_use]
    pub fn canonical_parts(&self) -> Vec<String> {
        match self {
            ClaimOrigin::SourceCertificate {
                producer,
                certificate_hash,
            } => vec![
                self.kind().to_string(),
                producer.clone(),
                certificate_hash.clone(),
            ],
            ClaimOrigin::AnchoredSource {
                dataset_id,
                content_hash,
            } => vec![
                self.kind().to_string(),
                dataset_id.clone(),
                content_hash.clone(),
            ],
            ClaimOrigin::EstimatedSource { estimator } => {
                vec![self.kind().to_string(), estimator.clone()]
            }
            ClaimOrigin::Derived => vec![self.kind().to_string()],
            ClaimOrigin::AuthenticatedWaiver(grant) => vec![
                self.kind().to_string(),
                grant.waiver_id.clone(),
                grant.expiry_day.to_string(),
                grant.mac.clone(),
            ],
        }
    }
}

/// A structured origin-validation refusal (field-level, teaching).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OriginError {
    /// Which claim.
    pub claim: String,
    /// The refusal.
    pub why: String,
}

impl fmt::Display for OriginError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "claim '{}': {}", self.claim, self.why)
    }
}

impl core::error::Error for OriginError {}

/// Shape-level validation shared by construction and parsing: non-blank
/// identities, canonical 64-hex hashes where required. Color-class
/// consistency and re-derivation live with the package verifier (they
/// need the claim's color, siblings, and the injected capabilities).
///
/// # Errors
/// [`OriginError`] naming the field.
pub fn validate_origin_shape(
    claim_id: &str,
    origin: &ClaimOrigin,
    is_canonical_hash: &dyn Fn(&str) -> bool,
) -> Result<(), OriginError> {
    let refuse = |why: String| {
        Err(OriginError {
            claim: claim_id.to_string(),
            why,
        })
    };
    match origin {
        ClaimOrigin::SourceCertificate {
            producer,
            certificate_hash,
        } => {
            if producer.trim().is_empty() {
                return refuse("source-certificate origin has a blank producer".to_string());
            }
            if !is_canonical_hash(certificate_hash) {
                return refuse(
                    "source-certificate origin needs a canonical 64-hex certificate hash"
                        .to_string(),
                );
            }
            Ok(())
        }
        ClaimOrigin::AnchoredSource {
            dataset_id,
            content_hash,
        } => {
            if dataset_id.trim().is_empty() {
                return refuse("anchored-source origin has a blank dataset id".to_string());
            }
            if !is_canonical_hash(content_hash) {
                return refuse(
                    "anchored-source origin needs a canonical 64-hex dataset hash".to_string(),
                );
            }
            Ok(())
        }
        ClaimOrigin::EstimatedSource { estimator } => {
            if estimator.trim().is_empty() {
                return refuse("estimated-source origin has a blank estimator".to_string());
            }
            Ok(())
        }
        ClaimOrigin::Derived => Ok(()),
        ClaimOrigin::AuthenticatedWaiver(grant) => {
            if grant.waiver_id.trim().is_empty() {
                return refuse("waiver origin has a blank waiver id".to_string());
            }
            if grant.mac.trim().is_empty() {
                return refuse("waiver origin has a blank authenticator".to_string());
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex64() -> String {
        "0123456789abcdef".repeat(4)
    }

    fn canonical(h: &str) -> bool {
        h.len() == 64
            && h.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    }

    #[test]
    fn shape_validation_fails_closed_per_field() {
        let ok = |o: &ClaimOrigin| validate_origin_shape("c", o, &canonical).is_ok();
        assert!(ok(&ClaimOrigin::SourceCertificate {
            producer: "fs-solver/ivp".to_string(),
            certificate_hash: hex64(),
        }));
        assert!(!ok(&ClaimOrigin::SourceCertificate {
            producer: "  ".to_string(),
            certificate_hash: hex64(),
        }));
        assert!(!ok(&ClaimOrigin::SourceCertificate {
            producer: "p".to_string(),
            certificate_hash: "deadbeef".to_string(),
        }));
        assert!(!ok(&ClaimOrigin::AnchoredSource {
            dataset_id: String::new(),
            content_hash: hex64(),
        }));
        assert!(!ok(&ClaimOrigin::EstimatedSource {
            estimator: " ".to_string(),
        }));
        assert!(ok(&ClaimOrigin::Derived));
        assert!(!ok(&ClaimOrigin::AuthenticatedWaiver(WaiverGrant {
            waiver_id: "w1".to_string(),
            expiry_day: 20_000,
            mac: "  ".to_string(),
        })));
    }

    #[test]
    fn canonical_parts_are_kind_prefixed_and_distinct() {
        let a = ClaimOrigin::EstimatedSource {
            estimator: "surrogate-v2".to_string(),
        };
        let b = ClaimOrigin::SourceCertificate {
            producer: "surrogate-v2".to_string(),
            certificate_hash: hex64(),
        };
        assert_eq!(a.canonical_parts()[0], "estimated-source");
        assert_ne!(a.canonical_parts(), b.canonical_parts());
        // The waiver's expiry and mac are bound (tamper moves the parts).
        let w1 = ClaimOrigin::AuthenticatedWaiver(WaiverGrant {
            waiver_id: "w".to_string(),
            expiry_day: 1,
            mac: "m".to_string(),
        });
        let w2 = ClaimOrigin::AuthenticatedWaiver(WaiverGrant {
            waiver_id: "w".to_string(),
            expiry_day: 2,
            mac: "m".to_string(),
        });
        assert_ne!(w1.canonical_parts(), w2.canonical_parts());
    }

    #[test]
    fn the_default_waiver_verifier_refuses_everything() {
        let grant = WaiverGrant {
            waiver_id: "w1".to_string(),
            expiry_day: u64::MAX,
            mac: "anything".to_string(),
        };
        assert!(!NoWaiverVerifier.verify(&grant, b"message"));
    }
}
