//! fs-checker — the standalone evidence-package checker (plan addendum,
//! Proposal 12). Layer: L6.
//!
//! "Don't trust us; here is the checker." A third party — an auditor or a
//! regulator — runs THIS to re-verify a FrankenSim [`EvidencePackage`] without
//! trusting the vendor and, crucially, WITHOUT running any solver. The whole
//! proposition rides the check/produce asymmetry (P6): re-verification is cheap
//! and needs neither the solver stack nor a license.
//!
//! # Hard distribution constraint
//! This crate depends only on `fs-package`; that package's production cone is
//! dependency-free `fs-blake3`, the static `fs-crosswalk` vocabulary, and
//! `fs-evidence` plus its observability utility. There is NO solver, geometry
//! kernel, or license gate anywhere in the cone, so by construction the checker
//! CANNOT run a solve. It carries its own protocol version
//! ([`CHECKER_PROTOCOL_VERSION`]) because it is distributed independently.
//!
//! What it re-verifies: format support + per-claim completeness (delegated to
//! [`EvidencePackage::verify`]), the content address (Merkle root, optionally
//! against an expected value — tamper detection), and signature validity when
//! an external capability is supplied. It renders a by-color budget pie only
//! after package verification succeeds. Everything is deterministic.

pub use fs_package::{
    ColorBreakdown, ContentHash, EvidencePackage, MagnitudeBudget, PackageError, ParseError,
};

/// The checker's own protocol version (it is distributed independently).
pub const CHECKER_PROTOCOL_VERSION: u32 = 2;

/// The one evidence-package format understood by this checker protocol.
///
/// Keep this as an explicit protocol literal rather than deriving it from
/// `fs-package`: a package-format change must make this crate fail to compile
/// until the independently distributed checker ABI is reviewed and versioned.
pub const CHECKER_SUPPORTED_PACKAGE_FORMAT: u32 = 4;
const _: () = assert!(CHECKER_SUPPORTED_PACKAGE_FORMAT == fs_package::FORMAT_VERSION);

/// Whether the package carried a detached signature and how far it was
/// verified (bead qmao.6.1): a present signature is only ASSERTED valid
/// when a [`SignatureVerifier`] capability accepts it over the
/// recomputed content root — presence alone is recorded, never
/// promoted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SignatureStatus {
    /// No signature — the bundle stands on its content address alone.
    Unsigned,
    /// A detached signature is present but NOT cryptographically
    /// verified (no capability supplied, or none exists in-tree).
    Unverified(String),
    /// The supplied verifier accepted the signature over the
    /// recomputed content root.
    Valid(String),
}

/// The signature-verification CAPABILITY (injected; this crate ships no
/// cryptography — the same fail-closed pattern as fs-ledger waivers).
pub trait SignatureVerifier {
    /// True iff `signature` authenticates the package's recomputed
    /// 32-byte BLAKE3 content root.
    fn verify(&self, merkle_root: &ContentHash, signature: &str) -> bool;
}

/// The in-tree default: nothing authenticates (no-crypto no-claim).
#[derive(Debug, Default, Clone, Copy)]
pub struct NoSignatureVerifier;

impl SignatureVerifier for NoSignatureVerifier {
    fn verify(&self, _merkle_root: &ContentHash, _signature: &str) -> bool {
        false
    }
}

/// The checker's overall verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// The package re-verified.
    Pass,
    /// The package failed re-verification (see the findings).
    Fail,
}

/// One reason a check failed.
#[derive(Debug, Clone, PartialEq)]
pub struct Finding {
    /// A short kind slug.
    pub kind: &'static str,
    /// Human detail.
    pub detail: String,
}

/// The result of running the checker over a package.
#[derive(Debug, Clone, PartialEq)]
pub struct CheckReport {
    /// Pass/Fail.
    pub verdict: Verdict,
    /// The recomputed content address (domain-separated BLAKE3 root).
    pub merkle_root: ContentHash,
    /// The by-color budget pie.
    pub breakdown: ColorBreakdown,
    /// Signature presence.
    pub signature: SignatureStatus,
    /// Reasons for failure (empty on Pass).
    pub findings: Vec<Finding>,
}

impl CheckReport {
    /// Did the package pass?
    #[must_use]
    pub fn passed(&self) -> bool {
        matches!(self.verdict, Verdict::Pass)
    }

    /// Render the by-color budget pie as a deterministic text chart.
    #[must_use]
    pub fn render_pie(&self) -> String {
        use core::fmt::Write as _;
        let b = &self.breakdown;
        let total = b.verified + b.validated + b.estimated;
        if total == 0 {
            return "budget pie: no claims".to_string();
        }
        let mut out = format!("budget pie ({total} claims):\n");
        for (label, count) in [
            ("verified ", b.verified),
            ("validated", b.validated),
            ("estimated", b.estimated),
        ] {
            // ten-cell bar, deterministic integer rounding.
            let filled = (count * 10 + total / 2) / total;
            let pct = (count * 100 + total / 2) / total;
            let bar: String = (0..10)
                .map(|i| if i < filled { '#' } else { '.' })
                .collect();
            writeln!(out, "  {label} {bar} {count} ({pct}%)").expect("write to String");
        }
        out
    }
}

/// Re-verify a package (no expected content address, no signature
/// capability — presence recorded, never asserted).
#[must_use]
pub fn check(pkg: &EvidencePackage) -> CheckReport {
    build_report(pkg, None, None)
}

/// Re-verify a package AND confirm its content address matches `expected_root`
/// — a mismatch (tamper, or the wrong package) fails the check.
#[must_use]
pub fn check_against_root(pkg: &EvidencePackage, expected_root: ContentHash) -> CheckReport {
    build_report(pkg, Some(expected_root), None)
}

/// The full third-party entry point (bead qmao.6.1): parse the
/// serialized package STRICTLY (schema v4 — the parser itself
/// recomputes the content root and re-derives the magnitude budget
/// from the parsed claims), then re-verify semantics, optionally
/// against an expected root and a signature capability. A package that
/// fails parsing never produces a Pass.
#[must_use]
pub fn check_json(
    text: &str,
    expected_root: Option<ContentHash>,
    verifier: Option<&dyn SignatureVerifier>,
) -> CheckReport {
    match EvidencePackage::from_json(text) {
        Ok(pkg) => build_report(&pkg, expected_root, verifier),
        Err(e) => CheckReport {
            verdict: Verdict::Fail,
            // Fail-closed sentinel: parsing refused, so there is no
            // recomputed root. The Fail verdict is authoritative; the
            // zero bytes are only a deterministic placeholder.
            merkle_root: ContentHash([0u8; 32]),
            breakdown: ColorBreakdown::default(),
            signature: SignatureStatus::Unsigned,
            findings: vec![Finding {
                kind: "parse-refused",
                detail: e.to_string(),
            }],
        },
    }
}

/// [`check`] with a signature-verification capability.
#[must_use]
pub fn check_with(
    pkg: &EvidencePackage,
    expected_root: Option<ContentHash>,
    verifier: &dyn SignatureVerifier,
) -> CheckReport {
    build_report(pkg, expected_root, Some(verifier))
}

fn build_report(
    pkg: &EvidencePackage,
    expected_root: Option<ContentHash>,
    verifier: Option<&dyn SignatureVerifier>,
) -> CheckReport {
    let mut findings = Vec::new();

    // 1. delegate format + per-claim completeness to the package format.
    let breakdown = match pkg.verify() {
        Ok(report) => report.breakdown,
        Err(e) => {
            findings.push(describe(&e));
            // Invalid claims must not retain a normal-looking positive
            // evidence summary. The finding still identifies the exact
            // refusal; the pie fails closed to no admitted claims.
            ColorBreakdown::default()
        }
    };

    // 2. content address (recomputed here, independently).
    let merkle_root = pkg.merkle_root();
    if let Some(expected) = expected_root
        && merkle_root != expected
    {
        findings.push(Finding {
            kind: "content-address-mismatch",
            detail: format!("recomputed root {merkle_root} != expected {expected}"),
        });
    }

    // 3. the magnitude budget must reconcile with its parts (the pie
    // is over error magnitudes, not claim counts — and it must not be
    // able to drift from the claims it summarizes).
    let mb = pkg.magnitude_budget();
    if mb.quantified_total.to_bits() != (mb.verified_width + mb.estimated_dispersion).to_bits() {
        findings.push(Finding {
            kind: "magnitude-budget-drift",
            detail: "quantified total does not reconcile with its parts".to_string(),
        });
    }

    // 4. signature: presence recorded; VALIDITY only through the
    // supplied capability, over the recomputed root (fail closed).
    let signature = match (&pkg.signature, verifier) {
        (None, _) => SignatureStatus::Unsigned,
        (Some(s), None) => SignatureStatus::Unverified(s.clone()),
        (Some(s), Some(v)) => {
            if v.verify(&merkle_root, s) {
                SignatureStatus::Valid(s.clone())
            } else {
                findings.push(Finding {
                    kind: "signature-invalid",
                    detail: "the supplied verifier rejected the detached signature over the \
                             recomputed content root"
                        .to_string(),
                });
                SignatureStatus::Unverified(s.clone())
            }
        }
    };

    let verdict = if findings.is_empty() {
        Verdict::Pass
    } else {
        Verdict::Fail
    };
    CheckReport {
        verdict,
        merkle_root,
        breakdown,
        signature,
        findings,
    }
}

/// Translate a package error into a checker finding.
fn describe(e: &PackageError) -> Finding {
    match e {
        PackageError::IncompleteProvenance { missing } => Finding {
            kind: "incomplete-provenance",
            detail: format!("package provenance is missing {missing}"),
        },
        PackageError::InvalidClaimId { index, id, reason } => Finding {
            kind: "invalid-claim-id",
            detail: format!("claim at index {index} has {reason} id {id:?}"),
        },
        PackageError::IncompleteValidatedClaim { claim, missing } => Finding {
            kind: "incomplete-validated-claim",
            detail: format!("claim '{claim}' is missing its {missing}"),
        },
        PackageError::IncompleteVerifiedClaim { claim } => Finding {
            kind: "incomplete-verified-claim",
            detail: format!("claim '{claim}' has no valid certificate interval"),
        },
        PackageError::InvalidValidatedRegime { claim, axis } => Finding {
            kind: "invalid-validated-regime",
            detail: format!("claim '{claim}' has an invalid validity axis {axis:?}"),
        },
        PackageError::IncompleteEstimatedClaim { claim, missing } => Finding {
            kind: "incomplete-estimated-claim",
            detail: format!("claim '{claim}' is missing its {missing}"),
        },
        PackageError::InvalidEstimatedDispersion { claim } => Finding {
            kind: "invalid-estimated-dispersion",
            detail: format!("claim '{claim}' has a NaN or negative dispersion"),
        },
        PackageError::MagnitudeOverflow { claim, component } => Finding {
            kind: "magnitude-overflow",
            detail: format!(
                "claim '{claim}' made finite {component} evidence overflow; explicit +infinity \
                 estimated dispersion is the only unbounded sentinel"
            ),
        },
        PackageError::TransportLimit { what, limit } => Finding {
            kind: "transport-limit",
            detail: format!("{what} exceeds the standalone checker limit {limit}"),
        },
        PackageError::UnsupportedFormat { found } => Finding {
            kind: "unsupported-format",
            detail: format!("package format version {found} is not supported"),
        },
        PackageError::ReceiptMismatch { claim } => Finding {
            kind: "receipt-mismatch",
            detail: format!(
                "claim '{claim}': re-running its composition receipt does not reproduce the \
                 claimed color — forged or stale derivation"
            ),
        },
        PackageError::BadReceiptParent { claim, parent } => Finding {
            kind: "bad-receipt-parent",
            detail: format!(
                "claim '{claim}': receipt parent {parent} is out of range or not strictly \
                 earlier in the package"
            ),
        },
        PackageError::RefutedClaim { claim, falsifier } => Finding {
            kind: "refuted-claim",
            detail: format!("claim '{claim}' was REFUTED by falsifier '{falsifier}'"),
        },
        PackageError::InvalidFalsifierRecord {
            claim,
            falsifier,
            field,
        } => Finding {
            kind: "invalid-falsifier-record",
            detail: format!(
                "claim '{claim}' has invalid falsifier record {falsifier}: {field} is missing or \
                 invalid"
            ),
        },
        PackageError::InvalidAnchorRecord {
            claim,
            anchor,
            field,
        } => Finding {
            kind: "invalid-anchor-record",
            detail: format!(
                "claim '{claim}' has invalid anchor record {anchor}: {field} is missing or \
                 non-canonical"
            ),
        },
    }
}
