//! Package-grounded crosswalk coverage (bead qmao.6.1): which
//! regulatory concepts are EVIDENCED by the fields actually present in
//! a package — the static fs-crosswalk mapping alone can never produce
//! "covered".

// ---------------------------------------------------------------------------
// Coverage from ACTUAL package fields (bead qmao.6.1; moved here from
// fs-crosswalk under the layer rule — L6 may know UTIL, never the
// reverse): the static fs-crosswalk table maps concepts to standards; whether a concept is EVIDENCED must
// come from a parsed package, never from the mapping itself — a mapped
// concept with no evidence in hand is "mapped but absent", not covered.
// ---------------------------------------------------------------------------

use crate::{
    AdmissionClass, EvidencePackage, PackageReport, SignaturePurpose, SignatureStatus,
    VerificationCapabilities, VerificationReceipt, VerifiedPackage,
};
use fs_crosswalk::{
    CROSSWALK_VERSION, Counterpart, PackageConcept, SUPPORTED_PACKAGE_FORMAT, Standard, crosswalk,
};
use fs_evidence::{Color, ColorRank};

/// Whether one concept is actually evidenced by a package, with the
/// reason either way.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConceptPresence {
    /// The concept.
    concept: PackageConcept,
    /// Is it evidenced by fields ACTUALLY PRESENT in this package?
    present: bool,
    /// Why / why not (teaching, deterministic).
    why: String,
}

impl ConceptPresence {
    /// Concept decided by this row.
    #[must_use]
    pub const fn concept(&self) -> PackageConcept {
        self.concept
    }

    /// Whether the verified package evidence establishes the concept. This bare
    /// value is non-authoritative when detached from a report whose decision
    /// hash and receipt have been validated.
    #[must_use]
    pub const fn present(&self) -> bool {
        self.present
    }

    /// Deterministic rationale for the decision.
    #[must_use]
    pub fn why(&self) -> &str {
        &self.why
    }
}

/// Receipt-bearing package-presence decision. Positive rows are evidentiary
/// only when `receipt` is present.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackagePresenceReport {
    /// One row per package concept.
    rows: Vec<ConceptPresence>,
    /// Exact package/policy decision that produced the rows.
    receipt: Option<VerificationReceipt>,
    /// Domain-separated integrity digest over the rows and package receipt.
    decision_hash: crate::ContentHash,
}

impl PackagePresenceReport {
    /// One sealed decision per package concept.
    #[must_use]
    pub fn rows(&self) -> &[ConceptPresence] {
        &self.rows
    }

    /// Borrow rows without discarding the report's receipt/hash context.
    pub fn iter(&self) -> core::slice::Iter<'_, ConceptPresence> {
        self.rows.iter()
    }

    /// Exact verification receipt used to derive positive rows, when package
    /// admission succeeded.
    #[must_use]
    pub fn receipt(&self) -> Option<&VerificationReceipt> {
        self.receipt.as_ref()
    }

    /// Stored integrity digest for this presence decision.
    #[must_use]
    pub const fn decision_hash(&self) -> crate::ContentHash {
        self.decision_hash
    }

    /// Recompute the digest over every authority-bearing field.
    #[must_use]
    pub fn validate_decision_hash(&self) -> bool {
        self.decision_hash == presence_report_hash(self)
    }
}

impl<'a> IntoIterator for &'a PackagePresenceReport {
    type Item = &'a ConceptPresence;
    type IntoIter = core::slice::Iter<'a, ConceptPresence>;

    fn into_iter(self) -> Self::IntoIter {
        self.rows.iter()
    }
}

fn is_scientific(report: &PackageReport, index: usize) -> bool {
    report
        .receipt()
        .admissions()
        .get(index)
        .is_some_and(|decision| decision.class() == AdmissionClass::Scientific)
}

fn anchoring_dataset_presence(pkg: &EvidencePackage, report: &PackageReport) -> (bool, String) {
    let named = pkg
        .claims
        .iter()
        .enumerate()
        .filter(|claim| {
            is_scientific(report, claim.0)
                && matches!(&claim.1.color, Color::Validated { dataset, .. } if !dataset.trim().is_empty())
        })
        .count();
    let matched = pkg
        .claims
        .iter()
        .enumerate()
        .filter(|(index, claim)| {
            is_scientific(report, *index)
                && claim.has_declared_matching_validated_anchor_unverified()
        })
        .count();
    (
        matched > 0,
        format!(
            "{matched} admitted validated claim(s) have an exact dataset-matching authenticated \
             anchor; {named} scientific dataset declaration(s); unrelated extra anchors are \
             excluded"
        ),
    )
}

fn signature_presence(report: &PackageReport) -> (bool, String) {
    match report.receipt().signature() {
        SignatureStatus::Authenticated(authenticated)
            if matches!(
                authenticated.purpose(),
                SignaturePurpose::ReleaseApproval { .. }
            ) =>
        {
            (
                true,
                "policy-authenticated signature carrying a release-approval purpose bound to an \
                 explicit checker protocol, expected package root, and scientific admission \
                 context; signature coverage does not establish checker release admission"
                    .to_string(),
            )
        }
        SignatureStatus::Authenticated(_) => (
            false,
            "generic package-root attestation is integrity evidence, not a release-purpose \
             signature and not checker release admission"
                .to_string(),
        ),
        SignatureStatus::Unverified(_) => (
            false,
            "detached signature present but not authenticated; raw presence cannot establish \
             a policy-authenticated root attestation"
                .to_string(),
        ),
        SignatureStatus::Refused { reason } => (
            false,
            format!("signature decision refused before bounded authentication: {reason}"),
        ),
        SignatureStatus::Unsigned => (false, "no detached signature".to_string()),
    }
}

fn rank_presence(
    pkg: &EvidencePackage,
    report: &PackageReport,
    rank: ColorRank,
    description: &str,
) -> (bool, String) {
    let count = pkg
        .claims
        .iter()
        .enumerate()
        .filter(|(index, claim)| is_scientific(report, *index) && claim.color.rank() == rank)
        .count();
    (
        count > 0,
        format!("{count} admitted scientific {description} claim(s)"),
    )
}

fn certificate_presence(pkg: &EvidencePackage, report: &PackageReport) -> (bool, String) {
    let count = pkg
        .claims
        .iter()
        .enumerate()
        .filter(|(index, claim)| {
            is_scientific(report, *index)
                && matches!(
                    &claim.color,
                    Color::Verified { lo, hi } if lo.is_finite() && hi.is_finite()
                )
        })
        .count();
    (
        count > 0,
        if count > 0 {
            format!("{count} admitted scientific Verified certificate claim(s)")
        } else {
            "no admitted scientific informative Verified certificate claim; vacuous infinite enclosures do not establish certificate coverage".to_string()
        },
    )
}

fn falsifier_presence(pkg: &EvidencePackage, report: &PackageReport) -> (bool, String) {
    let count: usize = pkg
        .claims
        .iter()
        .enumerate()
        .filter(|(index, _)| is_scientific(report, *index))
        .map(|(_, claim)| claim.falsifiers.len())
        .sum();
    (
        count > 0,
        if count > 0 {
            format!("{count} falsifier record(s) attached to claims")
        } else {
            "no falsifier records attached; absent evidence cannot claim coverage".to_string()
        },
    )
}

fn regime_presence(pkg: &EvidencePackage, report: &PackageReport) -> (bool, String) {
    let count = pkg
        .claims
        .iter()
        .enumerate()
        .filter(|(index, claim)| {
            is_scientific(report, *index)
                && matches!(&claim.color, Color::Validated { regime, .. } if !regime.bounds().is_empty())
        })
        .count();
    (
        count > 0,
        format!("{count} validated claim(s) with regime bounds"),
    )
}

fn provenance_presence(pkg: &EvidencePackage) -> (bool, String) {
    (
        false,
        format!(
            "declared code version {:?} and constellation lock {:?} are root-bound but no \
             provenance-artifact verifier authenticates their referenced content",
            pkg.provenance.code_version, pkg.provenance.constellation_lock
        ),
    )
}

fn claim_origin_presence(report: &PackageReport) -> (bool, String) {
    let policies = report.receipt().policy_fingerprints();
    let count = report
        .receipt()
        .admissions()
        .iter()
        .filter(|admission| {
            admission.class() == AdmissionClass::Scientific
                && match admission.origin_kind() {
                    crate::AdmissionOriginKind::SourceCertificate => {
                        policies.source_certificates().is_some()
                    }
                    crate::AdmissionOriginKind::AnchoredSource => {
                        policies.anchored_sources().is_some()
                    }
                    crate::AdmissionOriginKind::Derived => policies.derivations().is_some(),
                    crate::AdmissionOriginKind::AuthenticatedWaiver
                    | crate::AdmissionOriginKind::EstimatedSource => false,
                }
        })
        .count();
    (
        count > 0,
        if count > 0 {
            format!("{count} claim origin(s) authenticated by an external/derivation policy")
        } else {
            "no externally or derivationally authenticated claim origins; estimated-source \
             identities are declarations only"
                .to_string()
        },
    )
}

fn waiver_authorization_presence(pkg: &EvidencePackage, report: &PackageReport) -> (bool, String) {
    let count = pkg.waiver_claims();
    let dependent = report.breakdown().waived;
    (
        count > 0,
        if count > 0 {
            format!(
                "{count} waiver authorization(s) authenticated and unexpired; {dependent} \
                 direct/derived claim(s) remain waiver-dependent"
            )
        } else {
            "no authenticated waiver origins".to_string()
        },
    )
}

fn concept_presence(
    pkg: &EvidencePackage,
    report: &PackageReport,
    concept: PackageConcept,
) -> ConceptPresence {
    let (present, why) = match concept {
        PackageConcept::VerifiedColor => rank_presence(
            pkg,
            report,
            ColorRank::Verified,
            "verified interval-certificate",
        ),
        PackageConcept::ValidatedColor => {
            rank_presence(pkg, report, ColorRank::Validated, "validated")
        }
        PackageConcept::EstimatedColor => {
            rank_presence(pkg, report, ColorRank::Estimated, "estimated")
        }
        PackageConcept::Certificate => certificate_presence(pkg, report),
        PackageConcept::FalsifierLog => falsifier_presence(pkg, report),
        PackageConcept::RegimeTag => regime_presence(pkg, report),
        PackageConcept::AnchoringDataset => anchoring_dataset_presence(pkg, report),
        PackageConcept::Provenance => provenance_presence(pkg),
        PackageConcept::MerkleRoot => (
            true,
            format!(
                "content root {} recomputed and verified",
                report.merkle_root()
            ),
        ),
        PackageConcept::Signature => signature_presence(report),
        PackageConcept::ClaimOrigin => claim_origin_presence(report),
        PackageConcept::WaiverAuthorization => waiver_authorization_presence(pkg, report),
    };
    ConceptPresence {
        concept,
        present,
        why,
    }
}

/// Judge every concept against the fields actually present in `pkg`.
#[must_use]
pub fn package_presence(pkg: &EvidencePackage) -> PackagePresenceReport {
    package_presence_with(pkg, &VerificationCapabilities::deny_all())
}

/// Judge every concept only after package verification with the supplied
/// source, anchor, falsifier, derivation, waiver, and signature capabilities.
#[must_use]
pub fn package_presence_with(
    pkg: &EvidencePackage,
    capabilities: &VerificationCapabilities<'_>,
) -> PackagePresenceReport {
    let Ok(report) = pkg.verify_with(capabilities) else {
        let why = "package verification failed; fail-closed coverage suppresses every concept"
            .to_string();
        let mut report = PackagePresenceReport {
            rows: PackageConcept::ALL
                .iter()
                .map(|&concept| ConceptPresence {
                    concept,
                    present: false,
                    why: why.clone(),
                })
                .collect(),
            receipt: None,
            decision_hash: crate::ContentHash([0; 32]),
        };
        report.decision_hash = presence_report_hash(&report);
        return report;
    };
    package_presence_from_report(pkg, &report)
}

fn package_presence_from_report(
    pkg: &EvidencePackage,
    report: &PackageReport,
) -> PackagePresenceReport {
    let mut presence = PackagePresenceReport {
        rows: PackageConcept::ALL
            .iter()
            .map(|&concept| concept_presence(pkg, report, concept))
            .collect(),
        receipt: Some(report.receipt().clone()),
        decision_hash: crate::ContentHash([0; 32]),
    };
    presence.decision_hash = presence_report_hash(&presence);
    presence
}

fn presence_report_hash(report: &PackagePresenceReport) -> crate::ContentHash {
    fn atom(bytes: &mut Vec<u8>, value: &[u8]) {
        bytes.extend_from_slice(&(value.len() as u64).to_le_bytes());
        bytes.extend_from_slice(value);
    }

    let mut canonical = Vec::new();
    atom(&mut canonical, &(report.rows.len() as u64).to_le_bytes());
    match &report.receipt {
        Some(receipt) => {
            atom(&mut canonical, b"receipt");
            atom(&mut canonical, receipt.receipt_hash().as_bytes());
        }
        None => atom(&mut canonical, b"no-package-receipt"),
    }
    for row in &report.rows {
        atom(&mut canonical, row.concept.label().as_bytes());
        atom(&mut canonical, &[u8::from(row.present)]);
        atom(&mut canonical, row.why.as_bytes());
    }
    fs_blake3::hash_domain("fs-package:v7:presence-decision", &canonical)
}

/// Derive coverage rows from an already verified package without invoking any
/// external verifier a second time.
#[must_use]
pub fn verified_package_presence(verified: &VerifiedPackage) -> PackagePresenceReport {
    package_presence_from_report(verified.package(), verified.report())
}

/// One descriptive row status. A detached value is not an evidence decision;
/// authority resides only in a validated [`PackageCoverageReport`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CoverageStatus {
    /// Mapped to the standard AND evidenced by this package.
    Covered,
    /// Mapped to the standard but the evidence is ABSENT from this
    /// package — never reported as covered.
    MappedButAbsent,
    /// Deliberately unmapped for this standard (see the table's reason).
    NoClaim,
}

/// Receipt-bearing standards coverage decision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackageCoverageReport {
    /// Crosswalk rows for one standard.
    rows: Vec<(PackageConcept, CoverageStatus, String)>,
    /// Exact package/policy decision that produced the rows.
    receipt: Option<VerificationReceipt>,
    /// Regulatory standard whose mapping produced the rows.
    standard: Standard,
    /// Crosswalk vocabulary version.
    crosswalk_version: u32,
    /// Package format interpreted by this mapping.
    package_format: u32,
    /// Domain-separated integrity digest over mapping context, rows, and the
    /// package verification receipt.
    decision_hash: crate::ContentHash,
}

impl PackageCoverageReport {
    /// Standards mapping rows in deterministic crosswalk order.
    #[must_use]
    pub fn rows(&self) -> &[(PackageConcept, CoverageStatus, String)] {
        &self.rows
    }

    /// Borrow mapping rows without consuming the authority-bearing report.
    pub fn iter(&self) -> core::slice::Iter<'_, (PackageConcept, CoverageStatus, String)> {
        self.rows.iter()
    }
    /// Package verification receipt used for positive evidence, when present.
    #[must_use]
    pub fn receipt(&self) -> Option<&VerificationReceipt> {
        self.receipt.as_ref()
    }
    /// Regulatory standard evaluated by this report.
    #[must_use]
    pub const fn standard(&self) -> Standard {
        self.standard
    }
    /// Static crosswalk vocabulary version.
    #[must_use]
    pub const fn crosswalk_version(&self) -> u32 {
        self.crosswalk_version
    }
    /// Package format version interpreted by the crosswalk.
    #[must_use]
    pub const fn package_format(&self) -> u32 {
        self.package_format
    }
    /// Stored domain-separated coverage decision digest.
    #[must_use]
    pub const fn decision_hash(&self) -> crate::ContentHash {
        self.decision_hash
    }
    /// Whether the stored digest binds every report field.
    #[must_use]
    pub fn validate_decision_hash(&self) -> bool {
        self.decision_hash == coverage_report_hash(self)
    }
}

impl<'a> IntoIterator for &'a PackageCoverageReport {
    type Item = &'a (PackageConcept, CoverageStatus, String);
    type IntoIter = core::slice::Iter<'a, (PackageConcept, CoverageStatus, String)>;

    fn into_iter(self) -> Self::IntoIter {
        self.rows.iter()
    }
}

/// Coverage of `standard` derived from the fields actually present in
/// `pkg` (bead qmao.6.1): the static mapping alone can never produce
/// `Covered` — only the intersection with parsed package evidence can.
#[must_use]
pub fn package_coverage(pkg: &EvidencePackage, standard: Standard) -> PackageCoverageReport {
    package_coverage_with(pkg, standard, &VerificationCapabilities::deny_all())
}

/// Package-grounded standards coverage after verification with explicit
/// origin capabilities.
#[must_use]
pub fn package_coverage_with(
    pkg: &EvidencePackage,
    standard: Standard,
    capabilities: &VerificationCapabilities<'_>,
) -> PackageCoverageReport {
    let presence = package_presence_with(pkg, capabilities);
    coverage_from_presence(standard, presence)
}

fn coverage_from_presence(
    standard: Standard,
    presence: PackagePresenceReport,
) -> PackageCoverageReport {
    let PackagePresenceReport {
        rows: presence_rows,
        receipt,
        ..
    } = presence;
    let rows = crosswalk()
        .iter()
        .filter(|e| e.standard == standard)
        .map(|e| {
            let p = presence_rows
                .iter()
                .find(|pr| pr.concept() == e.concept)
                .expect("every concept judged");
            let status = if !e.is_mapped() {
                CoverageStatus::NoClaim
            } else if p.present() {
                CoverageStatus::Covered
            } else {
                CoverageStatus::MappedButAbsent
            };
            let mapping = match &e.counterpart {
                Counterpart::Mapped { clause, note } => {
                    format!("mapped to {clause}: {note}")
                }
                Counterpart::NoCounterpart { reason } => {
                    format!("no standard counterpart: {reason}")
                }
            };
            (
                e.concept,
                status,
                format!("{mapping}; package evidence: {}", p.why()),
            )
        })
        .collect();
    let mut report = PackageCoverageReport {
        rows,
        receipt,
        standard,
        crosswalk_version: CROSSWALK_VERSION,
        package_format: SUPPORTED_PACKAGE_FORMAT,
        decision_hash: crate::ContentHash([0; 32]),
    };
    report.decision_hash = coverage_report_hash(&report);
    report
}

fn coverage_report_hash(report: &PackageCoverageReport) -> crate::ContentHash {
    use core::fmt::Write as _;

    let mut canonical = String::new();
    let _ = write!(
        canonical,
        "standard:{}|crosswalk:{}|package:{}|rows:{}|",
        report.standard.label(),
        report.crosswalk_version,
        report.package_format,
        report.rows.len()
    );
    match &report.receipt {
        Some(receipt) => canonical.push_str(&receipt.receipt_hash().to_hex()),
        None => canonical.push_str("no-package-receipt"),
    }
    canonical.push('|');
    for (concept, status, why) in &report.rows {
        let _ = write!(
            canonical,
            "{}:{}:{}:",
            concept.label(),
            match status {
                CoverageStatus::Covered => "covered",
                CoverageStatus::MappedButAbsent => "mapped-but-absent",
                CoverageStatus::NoClaim => "no-claim",
            },
            why.len()
        );
        canonical.push_str(why);
        canonical.push('|');
    }
    fs_blake3::hash_domain("fs-package:v7:coverage-decision", canonical.as_bytes())
}

/// Derive standards coverage from an already verified package without
/// re-invoking external policy callbacks.
#[must_use]
pub fn verified_package_coverage(
    verified: &VerifiedPackage,
    standard: Standard,
) -> PackageCoverageReport {
    coverage_from_presence(standard, verified_package_presence(verified))
}
