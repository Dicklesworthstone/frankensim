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

/// Semantic version of the receipt-bearing package-presence decision.
pub const PRESENCE_DECISION_IDENTITY_VERSION: u32 = 8;
/// Exact BLAKE3 domain for the package-presence decision.
pub const PRESENCE_DECISION_IDENTITY_DOMAIN: &str = "fs-package:v8:presence-decision";

/// Semantic version of the receipt-bearing standards-coverage decision.
pub const COVERAGE_DECISION_IDENTITY_VERSION: u32 = 8;
/// Exact BLAKE3 domain for the standards-coverage decision.
pub const COVERAGE_DECISION_IDENTITY_DOMAIN: &str = "fs-package:v8:coverage-decision";
const _: () = assert!(PRESENCE_DECISION_IDENTITY_VERSION == crate::FORMAT_VERSION);
const _: () = assert!(COVERAGE_DECISION_IDENTITY_VERSION == crate::FORMAT_VERSION);

/// Owner-local declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const PRESENCE_DECISION_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-package:presence-decision",
    "version_const=PRESENCE_DECISION_IDENTITY_VERSION",
    "version=8",
    "domain=fs-package:v8:presence-decision",
    "domain_const=PRESENCE_DECISION_IDENTITY_DOMAIN",
    "encoder=presence_report_hash",
    "encoder_helpers=presence_report_hash_with_domain,append_presence_atom,admit_decision_hash",
    "schema_constants=PRESENCE_DECISION_IDENTITY_VERSION,PRESENCE_DECISION_IDENTITY_DOMAIN,crates/fs-package/src/lib.rs#FORMAT_VERSION,crates/fs-crosswalk/src/lib.rs#CROSSWALK_VERSION,crates/fs-crosswalk/src/lib.rs#SUPPORTED_PACKAGE_FORMAT",
    "schema_functions=admit_decision_hash,crates/fs-crosswalk/src/lib.rs#PackageConcept::label,package_presence_with,package_presence_from_report,concept_presence,is_scientific,anchoring_dataset_presence,signature_presence,rank_presence,certificate_presence,falsifier_presence,regime_presence,provenance_presence,claim_origin_presence,waiver_authorization_presence",
    "schema_dependencies=fs-package:verification-receipt",
    "digest=blake3-derive-key",
    "encoding=typed-binary",
    "sources=PackagePresenceReport,ConceptPresence",
    "source_fields=PackagePresenceReport.rows:derived:expanded-into-concept-presence-fields,PackagePresenceReport.receipt:semantic,PackagePresenceReport.decision_hash:derived:recomputed-from-semantic-fields,ConceptPresence.concept:semantic,ConceptPresence.present:semantic,ConceptPresence.why:semantic",
    "source_bindings=PackagePresenceReport.receipt>receipt-presence-and-hash,ConceptPresence.concept>ordered-concept-labels,ConceptPresence.present>ordered-presence-bits,ConceptPresence.why>ordered-rationale-utf8",
    "external_semantic_fields=identity-version,digest-domain,row-count",
    "semantic_fields=identity-version,digest-domain,row-count,receipt-presence-and-hash,ordered-concept-labels,ordered-presence-bits,ordered-rationale-utf8",
    "excluded_fields=none",
    "consumers=PackagePresenceReport::decision_hash,PackagePresenceReport::validate_decision_hash,package_presence,package_presence_with,verified_package_presence,coverage_from_presence",
    "mutations=identity-version:crates/fs-package/tests/package.rs#package_identity_versions_and_transports_fail_closed,digest-domain:crates/fs-package/src/coverage.rs#presence_decision_identity_fields_move_independently,row-count:crates/fs-package/src/coverage.rs#presence_decision_identity_fields_move_independently,receipt-presence-and-hash:crates/fs-package/src/coverage.rs#presence_decision_identity_fields_move_independently,ordered-concept-labels:crates/fs-package/src/coverage.rs#presence_decision_identity_fields_move_independently,ordered-presence-bits:crates/fs-package/src/coverage.rs#presence_decision_identity_fields_move_independently,ordered-rationale-utf8:crates/fs-package/src/coverage.rs#presence_decision_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_presence_decision_identity_fields",
    "transport_guard=PackagePresenceReport::admit_retained_decision_hash",
    "version_guard=crates/fs-package/tests/package.rs#package_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-package:presence-decision",
];

/// Owner-local declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const COVERAGE_DECISION_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-package:coverage-decision",
    "version_const=COVERAGE_DECISION_IDENTITY_VERSION",
    "version=8",
    "domain=fs-package:v8:coverage-decision",
    "domain_const=COVERAGE_DECISION_IDENTITY_DOMAIN",
    "encoder=coverage_report_hash",
    "encoder_helpers=coverage_report_hash_with_domain,CoverageDecisionRow::from_tuple",
    "schema_constants=COVERAGE_DECISION_IDENTITY_VERSION,COVERAGE_DECISION_IDENTITY_DOMAIN,crates/fs-package/src/lib.rs#FORMAT_VERSION,crates/fs-crosswalk/src/lib.rs#CROSSWALK_VERSION,crates/fs-crosswalk/src/lib.rs#SUPPORTED_PACKAGE_FORMAT",
    "schema_functions=admit_decision_hash,crates/fs-crosswalk/src/lib.rs#PackageConcept::label,crates/fs-crosswalk/src/lib.rs#Standard::label,crates/fs-blake3/src/lib.rs#ContentHash::to_hex,coverage_from_presence,crates/fs-crosswalk/src/lib.rs#crosswalk",
    "schema_dependencies=fs-package:verification-receipt",
    "digest=blake3-derive-key",
    "encoding=typed-binary",
    "sources=PackageCoverageReport,CoverageDecisionRow",
    "source_fields=PackageCoverageReport.rows:derived:expanded-into-coverage-row-fields,PackageCoverageReport.receipt:semantic,PackageCoverageReport.standard:semantic,PackageCoverageReport.crosswalk_version:semantic,PackageCoverageReport.package_format:semantic,PackageCoverageReport.decision_hash:derived:recomputed-from-semantic-fields,CoverageDecisionRow.concept:semantic,CoverageDecisionRow.status:semantic,CoverageDecisionRow.why:semantic",
    "source_bindings=PackageCoverageReport.receipt>receipt-presence-and-hash,PackageCoverageReport.standard>standard-label,PackageCoverageReport.crosswalk_version>crosswalk-version,PackageCoverageReport.package_format>package-format,CoverageDecisionRow.concept>ordered-concept-labels,CoverageDecisionRow.status>ordered-status-tags,CoverageDecisionRow.why>ordered-rationale-byte-counts-and-utf8",
    "external_semantic_fields=identity-version,digest-domain,row-count",
    "semantic_fields=identity-version,digest-domain,row-count,receipt-presence-and-hash,standard-label,crosswalk-version,package-format,ordered-concept-labels,ordered-status-tags,ordered-rationale-byte-counts-and-utf8",
    "excluded_fields=none",
    "consumers=PackageCoverageReport::decision_hash,PackageCoverageReport::validate_decision_hash,package_coverage,package_coverage_with,verified_package_coverage,standards-release-reports",
    "mutations=identity-version:crates/fs-package/tests/package.rs#package_identity_versions_and_transports_fail_closed,digest-domain:crates/fs-package/src/coverage.rs#coverage_decision_identity_fields_move_independently,row-count:crates/fs-package/src/coverage.rs#coverage_decision_identity_fields_move_independently,receipt-presence-and-hash:crates/fs-package/src/coverage.rs#coverage_decision_identity_fields_move_independently,standard-label:crates/fs-package/src/coverage.rs#coverage_decision_identity_fields_move_independently,crosswalk-version:crates/fs-package/src/coverage.rs#coverage_decision_identity_fields_move_independently,package-format:crates/fs-package/src/coverage.rs#coverage_decision_identity_fields_move_independently,ordered-concept-labels:crates/fs-package/src/coverage.rs#coverage_decision_identity_fields_move_independently,ordered-status-tags:crates/fs-package/src/coverage.rs#coverage_decision_identity_fields_move_independently,ordered-rationale-byte-counts-and-utf8:crates/fs-package/src/coverage.rs#coverage_decision_identity_fields_move_independently",
    "nonsemantic_mutations=none",
    "field_guard=classify_coverage_decision_identity_fields",
    "transport_guard=PackageCoverageReport::admit_retained_decision_hash",
    "version_guard=crates/fs-package/tests/package.rs#package_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-package:coverage-decision",
];

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

#[allow(dead_code)]
fn classify_presence_decision_identity_fields(
    report: &PackagePresenceReport,
    row: &ConceptPresence,
) {
    let PackagePresenceReport {
        rows,
        receipt,
        decision_hash,
    } = report;
    let ConceptPresence {
        concept,
        present,
        why,
    } = row;
    let _ = (rows, receipt, decision_hash, concept, present, why);
}

fn admit_decision_hash(
    found_version: u32,
    expected_version: u32,
    bytes: &[u8],
) -> Option<crate::ContentHash> {
    if found_version != expected_version || bytes.len() != 32 {
        return None;
    }
    let mut exact = [0_u8; 32];
    exact.copy_from_slice(bytes);
    Some(crate::ContentHash(exact))
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

    /// Admit a retained presence-decision digest only under the exact schema
    /// version and fixed-width binary transport.
    #[must_use]
    pub fn admit_retained_decision_hash(version: u32, bytes: &[u8]) -> Option<crate::ContentHash> {
        admit_decision_hash(version, PRESENCE_DECISION_IDENTITY_VERSION, bytes)
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

fn append_presence_atom(bytes: &mut Vec<u8>, value: &[u8]) {
    bytes.extend_from_slice(&(value.len() as u64).to_le_bytes());
    bytes.extend_from_slice(value);
}

fn presence_report_hash(report: &PackagePresenceReport) -> crate::ContentHash {
    presence_report_hash_with_domain(report, PRESENCE_DECISION_IDENTITY_DOMAIN)
}

fn presence_report_hash_with_domain(
    report: &PackagePresenceReport,
    domain: &str,
) -> crate::ContentHash {
    let mut canonical = Vec::new();
    append_presence_atom(&mut canonical, &(report.rows.len() as u64).to_le_bytes());
    match &report.receipt {
        Some(receipt) => {
            append_presence_atom(&mut canonical, b"receipt");
            append_presence_atom(&mut canonical, receipt.receipt_hash().as_bytes());
        }
        None => append_presence_atom(&mut canonical, b"no-package-receipt"),
    }
    for row in &report.rows {
        append_presence_atom(&mut canonical, row.concept.label().as_bytes());
        append_presence_atom(&mut canonical, &[u8::from(row.present)]);
        append_presence_atom(&mut canonical, row.why.as_bytes());
    }
    fs_blake3::hash_domain(domain, &canonical)
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

#[derive(Debug, Clone, Copy)]
struct CoverageDecisionRow<'a> {
    concept: PackageConcept,
    status: &'static str,
    why: &'a str,
}

impl<'a> CoverageDecisionRow<'a> {
    fn from_tuple(row: &'a (PackageConcept, CoverageStatus, String)) -> Self {
        Self {
            concept: row.0,
            status: match &row.1 {
                CoverageStatus::Covered => "covered",
                CoverageStatus::MappedButAbsent => "mapped-but-absent",
                CoverageStatus::NoClaim => "no-claim",
            },
            why: &row.2,
        }
    }
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

#[allow(dead_code)]
fn classify_coverage_decision_identity_fields(
    report: &PackageCoverageReport,
    row: &CoverageDecisionRow<'_>,
) {
    let PackageCoverageReport {
        rows,
        receipt,
        standard,
        crosswalk_version,
        package_format,
        decision_hash,
    } = report;
    let CoverageDecisionRow {
        concept,
        status,
        why,
    } = row;
    let _ = (
        rows,
        receipt,
        standard,
        crosswalk_version,
        package_format,
        decision_hash,
        concept,
        status,
        why,
    );
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
    /// Admit a retained coverage-decision digest only under the exact schema
    /// version and fixed-width binary transport.
    #[must_use]
    pub fn admit_retained_decision_hash(version: u32, bytes: &[u8]) -> Option<crate::ContentHash> {
        admit_decision_hash(version, COVERAGE_DECISION_IDENTITY_VERSION, bytes)
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
    coverage_report_hash_with_domain(report, COVERAGE_DECISION_IDENTITY_DOMAIN)
}

fn coverage_report_hash_with_domain(
    report: &PackageCoverageReport,
    domain: &str,
) -> crate::ContentHash {
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
    for tuple in &report.rows {
        let row = CoverageDecisionRow::from_tuple(tuple);
        let _ = write!(
            canonical,
            "{}:{}:{}:",
            row.concept.label(),
            row.status,
            row.why.len()
        );
        canonical.push_str(row.why);
        canonical.push('|');
    }
    fs_blake3::hash_domain(domain, canonical.as_bytes())
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

#[cfg(test)]
mod identity_tests {
    use super::*;

    fn empty_package_receipt() -> VerificationReceipt {
        EvidencePackage::new(crate::Provenance::new("commit-a", "lock-a"))
            .verify()
            .expect("an empty bounded package admits under deny-all policy")
            .receipt()
            .clone()
    }

    fn presence_fixture() -> PackagePresenceReport {
        let mut report = PackagePresenceReport {
            rows: vec![ConceptPresence {
                concept: PackageConcept::MerkleRoot,
                present: true,
                why: "bounded package root recomputed".to_string(),
            }],
            receipt: None,
            decision_hash: crate::ContentHash([0; 32]),
        };
        report.decision_hash = presence_report_hash(&report);
        report
    }

    fn coverage_fixture() -> PackageCoverageReport {
        let mut report = PackageCoverageReport {
            rows: vec![(
                PackageConcept::MerkleRoot,
                CoverageStatus::Covered,
                "mapped root and package evidence agree".to_string(),
            )],
            receipt: None,
            standard: Standard::AsmeVvV10,
            crosswalk_version: CROSSWALK_VERSION,
            package_format: SUPPORTED_PACKAGE_FORMAT,
            decision_hash: crate::ContentHash([0; 32]),
        };
        report.decision_hash = coverage_report_hash(&report);
        report
    }

    fn assert_hash_moves(baseline: crate::ContentHash, changed: crate::ContentHash, field: &str) {
        assert_ne!(
            baseline, changed,
            "semantic identity field did not move: {field}"
        );
    }

    #[test]
    fn presence_decision_identity_fields_move_independently() {
        let report = presence_fixture();
        let baseline = presence_report_hash(&report);
        assert_hash_moves(
            baseline,
            presence_report_hash_with_domain(&report, "fs-package:v8:alternate-presence-decision"),
            "digest-domain",
        );

        let mut changed = report.clone();
        changed.rows.push(ConceptPresence {
            concept: PackageConcept::Provenance,
            present: false,
            why: "missing reproducibility identity".to_string(),
        });
        assert_hash_moves(baseline, presence_report_hash(&changed), "row-count");

        let mut changed = report.clone();
        changed.receipt = Some(empty_package_receipt());
        assert_hash_moves(
            baseline,
            presence_report_hash(&changed),
            "receipt-presence-and-hash",
        );

        let mut changed = report.clone();
        changed.rows[0].concept = PackageConcept::Provenance;
        assert_hash_moves(
            baseline,
            presence_report_hash(&changed),
            "ordered-concept-labels",
        );

        let mut changed = report.clone();
        changed.rows[0].present = false;
        assert_hash_moves(
            baseline,
            presence_report_hash(&changed),
            "ordered-presence-bits",
        );

        let mut changed = report;
        changed.rows[0].why.push_str(" under schema v8");
        assert_hash_moves(
            baseline,
            presence_report_hash(&changed),
            "ordered-rationale-utf8",
        );
    }

    #[test]
    fn coverage_decision_identity_fields_move_independently() {
        let report = coverage_fixture();
        let baseline = coverage_report_hash(&report);
        assert_hash_moves(
            baseline,
            coverage_report_hash_with_domain(&report, "fs-package:v8:alternate-coverage-decision"),
            "digest-domain",
        );

        let mut changed = report.clone();
        changed.rows.push((
            PackageConcept::Provenance,
            CoverageStatus::MappedButAbsent,
            "mapped but missing".to_string(),
        ));
        assert_hash_moves(baseline, coverage_report_hash(&changed), "row-count");

        let mut changed = report.clone();
        changed.receipt = Some(empty_package_receipt());
        assert_hash_moves(
            baseline,
            coverage_report_hash(&changed),
            "receipt-presence-and-hash",
        );

        let mut changed = report.clone();
        changed.standard = Standard::AsmeVvV20;
        assert_hash_moves(baseline, coverage_report_hash(&changed), "standard-label");

        let mut changed = report.clone();
        changed.crosswalk_version += 1;
        assert_hash_moves(
            baseline,
            coverage_report_hash(&changed),
            "crosswalk-version",
        );

        let mut changed = report.clone();
        changed.package_format += 1;
        assert_hash_moves(baseline, coverage_report_hash(&changed), "package-format");

        let mut changed = report.clone();
        changed.rows[0].0 = PackageConcept::Provenance;
        assert_hash_moves(
            baseline,
            coverage_report_hash(&changed),
            "ordered-concept-labels",
        );

        let mut changed = report.clone();
        changed.rows[0].1 = CoverageStatus::NoClaim;
        assert_hash_moves(
            baseline,
            coverage_report_hash(&changed),
            "ordered-status-tags",
        );

        let mut changed = report;
        changed.rows[0].2.push_str(" under schema v8");
        assert_hash_moves(
            baseline,
            coverage_report_hash(&changed),
            "ordered-rationale-byte-counts-and-utf8",
        );
    }
}
