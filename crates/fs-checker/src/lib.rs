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
//! What it re-verifies, in authority order: the content address (optionally
//! against an expected value), callback-free package structure, every attached
//! portable semantic witness through a closed built-in registry, then external
//! origins and signatures through explicit capabilities. It renders a by-color
//! budget pie only after package verification succeeds. Everything is
//! deterministic for deterministic injected capabilities.

mod semantic;

pub use semantic::{
    BOUNDED_LINF_RESIDUAL_FAMILY, EXACT_INTERVAL_FAMILY, INITIAL_SEMANTIC_SCHEMA_VERSION,
    MAX_INTERVAL_NODES, MAX_RESIDUAL_DIMENSION, MAX_RESIDUAL_MATRIX_ENTRIES,
    MAX_SEMANTIC_OPERATIONS, MAX_SEMANTIC_PAYLOAD_BYTES, MAX_SEMANTIC_WITNESS_BYTES,
    MAX_SEMANTIC_WITNESSES, SEMANTIC_IMPLEMENTATION_VERSION, SEMANTIC_PLUGIN_IDENTITY_DOMAIN,
    SEMANTIC_PLUGIN_IDENTITY_VERSION, SEMANTIC_REGISTRY_IDENTITY_DOMAIN,
    SEMANTIC_REGISTRY_IDENTITY_VERSION, SEMANTIC_REPORT_IDENTITY_DOMAIN,
    SEMANTIC_REPORT_IDENTITY_VERSION, SemanticClaimReceipt, SemanticClaimStatus, SemanticFailure,
    SemanticFailureKind, SemanticPluginDescriptor, SemanticReport, SemanticStatus,
    admit_retained_semantic_registry_fingerprint, semantic_plugin_registry,
    semantic_registry_fingerprint, verify_portable_semantics,
};

pub use fs_package::{
    AdmissionClass, AdmissionOriginKind, AnchoredSourceRequest, AnchoredSourceVerifier,
    ClaimAdmission, ColorBreakdown, ContentHash, DerivationRequest, DerivationVerifier,
    EvidencePackage, FalsifierRequest, FalsifierVerifier, MagnitudeBudget,
    NoAnchoredSourceVerifier, NoDerivationVerifier, NoFalsifierVerifier, NoSignatureVerifier,
    NoSourceCertificateVerifier, NoWaiverVerifier, PackageError, ParseError, PolicyFingerprint,
    SignatureIntent, SignaturePurpose, SignatureRequest, SignatureStatus, SignatureVerification,
    SignatureVerifier, SourceCertificateRequest, SourceCertificateVerifier,
    VerificationCapabilities, VerificationDecision, VerificationPolicyFingerprints,
    VerificationReceipt, VerifiedPackage, WaiverGrant, WaiverVerification, WaiverVerifier,
    hash_checker_decision, signature_subject_hash,
};

/// The checker's own protocol version (it is distributed independently).
pub const CHECKER_PROTOCOL_VERSION: u32 = 6;

/// The one evidence-package format understood by this checker protocol.
///
/// Keep this as an explicit protocol literal rather than deriving it from
/// `fs-package`: a package-format change must make this crate fail to compile
/// until the independently distributed checker ABI is reviewed and versioned.
pub const CHECKER_SUPPORTED_PACKAGE_FORMAT: u32 = 8;
const _: () = assert!(CHECKER_SUPPORTED_PACKAGE_FORMAT == fs_package::FORMAT_VERSION);

/// Semantic version of the retained checker-decision digest.
pub const CHECKER_DECISION_IDENTITY_VERSION: u32 = 8;
/// Exact derive-key domain used by `fs-package` for checker decisions.
pub use fs_package::CHECKER_DECISION_IDENTITY_DOMAIN;
const _: () = assert!(CHECKER_DECISION_IDENTITY_VERSION == CHECKER_SUPPORTED_PACKAGE_FORMAT);

/// Owner-local checker-decision declaration consumed by `xtask check-identities`.
#[allow(dead_code)]
pub const CHECKER_DECISION_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-checker:decision-report",
    "version_const=CHECKER_DECISION_IDENTITY_VERSION",
    "version=8",
    "domain=fs-package:v8:checker-decision",
    "domain_const=crates/fs-package/src/lib.rs#CHECKER_DECISION_IDENTITY_DOMAIN",
    "encoder=checker_report_hash",
    "encoder_helpers=checker_report_hash_with_protocol,checker_decision_atom,append_signature_identity,append_authenticated_signature_identity",
    "schema_constants=CHECKER_DECISION_IDENTITY_VERSION,CHECKER_PROTOCOL_VERSION,CHECKER_SUPPORTED_PACKAGE_FORMAT,crates/fs-package/src/lib.rs#CHECKER_DECISION_IDENTITY_DOMAIN,crates/fs-package/src/lib.rs#FORMAT_VERSION",
    "schema_functions=CheckReport::admit_retained_decision_hash,CheckReport::validate_decision_hash,crates/fs-checker/src/semantic.rs#SemanticReport::validate_context_hash,crates/fs-package/src/lib.rs#hash_checker_decision",
    "schema_dependencies=fs-checker:semantic-report,fs-package:package-root,fs-package:signature-subject,fs-package:verification-receipt,fs-package:release-admission-context",
    "digest=blake3-derive-key",
    "encoding=typed-binary",
    "sources=CheckReport,Finding",
    "source_fields=CheckReport.verdict:semantic,CheckReport.merkle_root:semantic,CheckReport.breakdown:semantic,CheckReport.integrity_status:semantic,CheckReport.semantic_report:semantic,CheckReport.origin_status:semantic,CheckReport.signature:semantic,CheckReport.receipt:semantic,CheckReport.findings:semantic,CheckReport.policy:semantic,CheckReport.expected_root:semantic,CheckReport.decision_hash:derived:recomputed-from-semantic-fields,Finding.kind:semantic,Finding.detail:semantic",
    "source_bindings=CheckReport.verdict>verdict,CheckReport.merkle_root>package-root,CheckReport.breakdown>verified-count+validated-count+estimated-count+waived-count,CheckReport.integrity_status>integrity-status,CheckReport.semantic_report>semantic-status+semantic-context-hash,CheckReport.origin_status>origin-status,CheckReport.signature>signature-status+signature-payload+signature-purpose,CheckReport.receipt>verification-receipt-presence+verification-receipt-hash,CheckReport.findings>finding-count+finding-order,CheckReport.policy>policy,CheckReport.expected_root>expected-root-presence+expected-root,Finding.kind>finding-kind,Finding.detail>finding-detail",
    "external_semantic_fields=identity-version,digest-domain,checker-protocol-version",
    "semantic_fields=identity-version,digest-domain,checker-protocol-version,verdict,package-root,verified-count,validated-count,estimated-count,waived-count,integrity-status,semantic-status,semantic-context-hash,origin-status,signature-status,signature-payload,signature-purpose,verification-receipt-presence,verification-receipt-hash,finding-count,finding-order,policy,expected-root-presence,expected-root,finding-kind,finding-detail",
    "excluded_fields=none",
    "consumers=CheckReport::decision_hash,CheckReport::validate_decision_hash,CheckReport::release_admitted,release-approval-auditors",
    "mutations=identity-version:crates/fs-checker/src/lib.rs#checker_decision_identity_versions_and_transports_fail_closed,digest-domain:crates/fs-checker/src/lib.rs#checker_decision_identity_versions_and_transports_fail_closed,checker-protocol-version:crates/fs-checker/src/lib.rs#checker_decision_identity_versions_and_transports_fail_closed,verdict:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,package-root:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,verified-count:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,validated-count:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,estimated-count:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,waived-count:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,integrity-status:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,semantic-status:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,semantic-context-hash:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,origin-status:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,signature-status:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,signature-payload:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,signature-purpose:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,verification-receipt-presence:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,verification-receipt-hash:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,finding-count:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,finding-order:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,policy:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,expected-root-presence:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,expected-root:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,finding-kind:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field,finding-detail:crates/fs-checker/src/lib.rs#decision_hash_binds_every_checker_authority_field",
    "nonsemantic_mutations=none",
    "field_guard=classify_checker_decision_identity_fields",
    "transport_guard=CheckReport::admit_retained_decision_hash",
    "version_guard=crates/fs-checker/src/lib.rs#checker_decision_identity_versions_and_transports_fail_closed",
    "coupling_surface=fs-checker:decision-report",
];

/// The checker's overall verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// The package re-verified.
    Pass,
    /// The package failed re-verification (see the findings).
    Fail,
}

/// Checker gate whose policy produced a report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckPolicy {
    /// Integrity/origin admission; not release approval.
    Integrity,
    /// Non-admitting release readiness inventory.
    ReleasePreflight,
    /// Strong release admission under protocol v6.
    ReleaseAdmission,
}

/// Callback-free package-integrity stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntegrityStatus {
    /// Transport, format, content binding, and claim structure were rechecked.
    Verified,
    /// Integrity refused before semantic or origin authority could be granted.
    Refused,
}

/// External-origin authentication stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OriginStatus {
    /// Package admission, including every required external origin, succeeded.
    Authenticated,
    /// An invoked package/origin capability refused admission.
    Refused,
    /// An earlier integrity or semantic stage stopped callback dispatch.
    NotRun,
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
    verdict: Verdict,
    /// The recomputed content address (domain-separated BLAKE3 root).
    merkle_root: ContentHash,
    /// The by-color budget pie.
    breakdown: ColorBreakdown,
    /// Callback-free package-integrity outcome.
    integrity_status: IntegrityStatus,
    /// Independent portable-witness transcript.
    semantic_report: SemanticReport,
    /// External-origin authentication outcome.
    origin_status: OriginStatus,
    /// Signature presence.
    signature: SignatureStatus,
    /// Policy-bound package receipt. Present only after successful package
    /// verification; parse and capability refusals expose no partial receipt.
    receipt: Option<VerificationReceipt>,
    /// Reasons for failure (empty on Pass).
    findings: Vec<Finding>,
    /// Gate intent bound into `decision_hash`.
    policy: CheckPolicy,
    /// External expected root, when the gate required one.
    expected_root: Option<ContentHash>,
    /// Domain-separated integrity digest over checker protocol, gate context,
    /// package receipt, verdict, and findings.
    decision_hash: ContentHash,
}

#[allow(dead_code)]
fn classify_checker_decision_identity_fields(report: &CheckReport, finding: &Finding) {
    let CheckReport {
        verdict,
        merkle_root,
        breakdown,
        integrity_status,
        semantic_report,
        origin_status,
        signature,
        receipt,
        findings,
        policy,
        expected_root,
        decision_hash,
    } = report;
    let ColorBreakdown {
        verified,
        validated,
        estimated,
        waived,
    } = breakdown;
    let Finding { kind, detail } = finding;
    let _ = (
        verdict,
        merkle_root,
        verified,
        validated,
        estimated,
        waived,
        integrity_status,
        semantic_report,
        origin_status,
        signature,
        receipt,
        findings,
        policy,
        expected_root,
        decision_hash,
        kind,
        detail,
    );
}

impl CheckReport {
    /// Overall gate verdict.
    #[must_use]
    pub const fn verdict(&self) -> Verdict {
        self.verdict
    }
    /// Bounded recomputed package root, or the zero refusal sentinel.
    #[must_use]
    pub const fn merkle_root(&self) -> ContentHash {
        self.merkle_root
    }
    /// Admitted evidence-color breakdown.
    #[must_use]
    pub const fn breakdown(&self) -> &ColorBreakdown {
        &self.breakdown
    }
    /// Callback-free package-integrity outcome.
    #[must_use]
    pub const fn integrity_status(&self) -> IntegrityStatus {
        self.integrity_status
    }
    /// Independent portable-witness package outcome.
    #[must_use]
    pub const fn semantic_status(&self) -> SemanticStatus {
        self.semantic_report.status()
    }
    /// Full sealed portable-witness transcript.
    #[must_use]
    pub const fn semantic_report(&self) -> &SemanticReport {
        &self.semantic_report
    }
    /// External-origin authentication outcome.
    #[must_use]
    pub const fn origin_status(&self) -> OriginStatus {
        self.origin_status
    }
    /// Detached-signature decision.
    #[must_use]
    pub fn signature(&self) -> &SignatureStatus {
        &self.signature
    }
    /// Package verification receipt, only after successful package admission.
    #[must_use]
    pub fn receipt(&self) -> Option<&VerificationReceipt> {
        self.receipt.as_ref()
    }
    /// Ordered deterministic refusal findings.
    #[must_use]
    pub fn findings(&self) -> &[Finding] {
        &self.findings
    }
    /// Gate policy bound into this report.
    #[must_use]
    pub const fn policy(&self) -> CheckPolicy {
        self.policy
    }
    /// Expected package root supplied by the caller, when any.
    #[must_use]
    pub const fn expected_root(&self) -> Option<ContentHash> {
        self.expected_root
    }
    /// Stored domain-separated checker decision digest.
    #[must_use]
    pub const fn decision_hash(&self) -> ContentHash {
        self.decision_hash
    }

    /// Admit retained checker-decision bytes only under the exact v8 schema
    /// and fixed-width digest transport. Stale and future versions refuse.
    #[must_use]
    pub fn admit_retained_decision_hash(version: u32, bytes: &[u8]) -> Option<ContentHash> {
        if version != CHECKER_DECISION_IDENTITY_VERSION || bytes.len() != 32 {
            return None;
        }
        let mut exact = [0_u8; 32];
        exact.copy_from_slice(bytes);
        Some(ContentHash(exact))
    }
    /// Did the package pass?
    #[must_use]
    pub fn passed(&self) -> bool {
        matches!(self.verdict, Verdict::Pass)
    }

    /// Whether this is a hash-valid, receipt-bearing release-admission pass.
    /// [`Self::passed`] is policy-local and must not be substituted for this at
    /// a release boundary.
    #[must_use]
    pub fn release_admitted(&self) -> bool {
        self.policy == CheckPolicy::ReleaseAdmission
            && self.passed()
            && self.receipt.is_some()
            && self.integrity_status == IntegrityStatus::Verified
            && matches!(
                self.semantic_report.status(),
                SemanticStatus::NotProvided | SemanticStatus::Verified
            )
            && self.semantic_report.validate_context_hash()
            && self.origin_status == OriginStatus::Authenticated
            && self.validate_decision_hash()
    }

    /// Whether release admission also independently recomputed at least one
    /// attached portable semantic witness.
    #[must_use]
    pub fn release_independently_verified(&self) -> bool {
        self.release_admitted() && self.semantic_report.status() == SemanticStatus::Verified
    }

    /// Recompute the checker gate-context integrity digest.
    #[must_use]
    pub fn recomputed_decision_hash(&self) -> ContentHash {
        checker_report_hash(self)
    }

    /// Whether the stored checker decision digest matches the report fields.
    #[must_use]
    pub fn validate_decision_hash(&self) -> bool {
        self.semantic_report.validate_context_hash()
            && self.decision_hash == self.recomputed_decision_hash()
    }

    /// Render the by-color budget pie as a deterministic text chart.
    #[must_use]
    pub fn render_pie(&self) -> String {
        use core::fmt::Write as _;
        let b = &self.breakdown;
        let total =
            b.verified as u128 + b.validated as u128 + b.estimated as u128 + b.waived as u128;
        if total == 0 {
            return "budget pie: no claims".to_string();
        }
        let mut out = format!("budget pie ({total} claims):\n");
        for (label, count) in [
            ("verified ", b.verified),
            ("validated", b.validated),
            ("estimated", b.estimated),
            ("waived   ", b.waived),
        ] {
            // ten-cell bar, deterministic integer rounding.
            let count_u128 = count as u128;
            let filled = ((count_u128 * 10 + total / 2) / total) as usize;
            let pct = (count_u128 * 100 + total / 2) / total;
            let bar: String = (0..10)
                .map(|i| if i < filled { '#' } else { '.' })
                .collect();
            writeln!(out, "  {label} {bar} {count} ({pct}%)").expect("write to String");
        }
        out
    }
}

fn checker_report_hash(report: &CheckReport) -> ContentHash {
    checker_report_hash_with_protocol(report, CHECKER_PROTOCOL_VERSION)
}

fn checker_decision_atom(out: &mut Vec<u8>, bytes: &[u8]) {
    out.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
    out.extend_from_slice(bytes);
}

fn append_authenticated_signature_identity(
    canonical: &mut Vec<u8>,
    signature: &str,
    purpose: SignaturePurpose,
) {
    checker_decision_atom(canonical, signature.as_bytes());
    match purpose {
        SignaturePurpose::PackageRootAttestation => {
            checker_decision_atom(canonical, b"package-root-attestation");
        }
        SignaturePurpose::ReleaseApproval {
            checker_protocol,
            expected_root,
            admission_context,
            semantic_context,
        } => {
            checker_decision_atom(canonical, b"release-approval");
            checker_decision_atom(canonical, &checker_protocol.to_le_bytes());
            checker_decision_atom(canonical, expected_root.as_bytes());
            checker_decision_atom(canonical, admission_context.as_bytes());
            checker_decision_atom(canonical, semantic_context.as_bytes());
        }
    }
}

fn append_signature_identity(canonical: &mut Vec<u8>, signature: &SignatureStatus) {
    match signature {
        SignatureStatus::Unsigned => checker_decision_atom(canonical, b"signature:unsigned"),
        SignatureStatus::Refused { reason } => {
            checker_decision_atom(canonical, b"signature:refused");
            checker_decision_atom(canonical, reason.as_bytes());
        }
        SignatureStatus::Unverified(signature) => {
            checker_decision_atom(canonical, b"signature:unverified");
            checker_decision_atom(canonical, signature.as_bytes());
        }
        SignatureStatus::Authenticated(authenticated) => {
            checker_decision_atom(canonical, b"signature:authenticated");
            append_authenticated_signature_identity(
                canonical,
                authenticated.signature(),
                authenticated.purpose(),
            );
        }
    }
}

fn checker_report_hash_with_protocol(
    report: &CheckReport,
    checker_protocol_version: u32,
) -> ContentHash {
    let mut canonical = Vec::new();
    checker_decision_atom(&mut canonical, &checker_protocol_version.to_le_bytes());
    checker_decision_atom(
        &mut canonical,
        match report.policy {
            CheckPolicy::Integrity => b"integrity",
            CheckPolicy::ReleasePreflight => b"release-preflight",
            CheckPolicy::ReleaseAdmission => b"release-admission",
        },
    );
    match report.expected_root {
        Some(root) => checker_decision_atom(&mut canonical, root.as_bytes()),
        None => checker_decision_atom(&mut canonical, b"no-expected-root"),
    }
    checker_decision_atom(&mut canonical, report.merkle_root.as_bytes());
    checker_decision_atom(
        &mut canonical,
        match report.integrity_status {
            IntegrityStatus::Verified => b"integrity:verified",
            IntegrityStatus::Refused => b"integrity:refused",
        },
    );
    checker_decision_atom(
        &mut canonical,
        match report.semantic_report.status() {
            SemanticStatus::NotProvided => b"semantics:not-provided",
            SemanticStatus::Verified => b"semantics:verified",
            SemanticStatus::Refused => b"semantics:refused",
            SemanticStatus::NotRun => b"semantics:not-run",
        },
    );
    checker_decision_atom(
        &mut canonical,
        report.semantic_report.context_hash().as_bytes(),
    );
    checker_decision_atom(
        &mut canonical,
        match report.origin_status {
            OriginStatus::Authenticated => b"origin:authenticated",
            OriginStatus::Refused => b"origin:refused",
            OriginStatus::NotRun => b"origin:not-run",
        },
    );
    checker_decision_atom(
        &mut canonical,
        match report.verdict {
            Verdict::Pass => b"pass",
            Verdict::Fail => b"fail",
        },
    );
    match &report.receipt {
        Some(receipt) => checker_decision_atom(&mut canonical, receipt.receipt_hash().as_bytes()),
        None => checker_decision_atom(&mut canonical, b"no-package-receipt"),
    }
    append_signature_identity(&mut canonical, &report.signature);
    let ColorBreakdown {
        verified,
        validated,
        estimated,
        waived,
    } = report.breakdown;
    for count in [verified, validated, estimated, waived] {
        checker_decision_atom(&mut canonical, &(count as u64).to_le_bytes());
    }
    checker_decision_atom(
        &mut canonical,
        &(report.findings.len() as u64).to_le_bytes(),
    );
    for finding in &report.findings {
        checker_decision_atom(&mut canonical, finding.kind.as_bytes());
        checker_decision_atom(&mut canonical, finding.detail.as_bytes());
    }
    hash_checker_decision(&canonical)
}

/// Re-verify a package (no expected content address, no signature
/// capability — presence recorded, never asserted).
#[must_use]
pub fn check(pkg: &EvidencePackage) -> CheckReport {
    check_with_capabilities(pkg, None, None, &VerificationCapabilities::deny_all())
}

/// Re-verify a package AND confirm its content address matches `expected_root`
/// — a mismatch (tamper, or the wrong package) fails the check.
#[must_use]
pub fn check_against_root(pkg: &EvidencePackage, expected_root: ContentHash) -> CheckReport {
    check_with_capabilities(
        pkg,
        Some(expected_root),
        None,
        &VerificationCapabilities::deny_all(),
    )
}

/// The full third-party entry point (bead qmao.6.1): parse the
/// serialized package under the deterministic schema-v8 JSON profile (the parser itself
/// recomputes the content root and re-derives the magnitude budget
/// from the parsed claims), then re-verify semantics, optionally
/// against an expected root and a signature capability. A package that
/// fails parsing never produces a Pass. Every external artifact/authorization
/// capability is denied by default; use [`check_json_with_capabilities`] for
/// explicit typed admission.
#[must_use]
pub fn check_json(
    text: &str,
    expected_root: Option<ContentHash>,
    verifier: Option<&dyn SignatureVerifier>,
) -> CheckReport {
    check_json_with_capabilities(
        text,
        expected_root,
        verifier,
        &VerificationCapabilities::deny_all(),
    )
}

/// Deterministic-profile JSON checking with explicit source, anchor, falsifier,
/// derivation, waiver, and signature capabilities. The separate signature
/// verifier remains optional for integrity decisions.
#[must_use]
pub fn check_json_with_capabilities(
    text: &str,
    expected_root: Option<ContentHash>,
    signature_verifier: Option<&dyn SignatureVerifier>,
    capabilities: &VerificationCapabilities<'_>,
) -> CheckReport {
    match EvidencePackage::from_json(text) {
        Ok(pkg) => build_report(
            &pkg,
            expected_root,
            signature_verifier,
            capabilities,
            CheckPolicy::Integrity,
        ),
        Err(e) => parse_refusal(&e, CheckPolicy::Integrity, expected_root),
    }
}

/// [`check`] with an independent signature-verification capability. Package
/// origins remain deny-all; use [`check_with_capabilities`] when source
/// certificates, anchors, falsifiers, derivations, or waivers are part of the
/// decision.
#[must_use]
pub fn check_with(
    pkg: &EvidencePackage,
    expected_root: Option<ContentHash>,
    verifier: &dyn SignatureVerifier,
) -> CheckReport {
    check_with_capabilities(
        pkg,
        expected_root,
        Some(verifier),
        &VerificationCapabilities::deny_all(),
    )
}

/// Re-verify an in-memory package with explicit origin-verification
/// capabilities. Detached-signature verification is a separate, optional
/// capability and does not become required merely because an origin verifier
/// was supplied.
#[must_use]
pub fn check_with_capabilities(
    pkg: &EvidencePackage,
    expected_root: Option<ContentHash>,
    signature_verifier: Option<&dyn SignatureVerifier>,
    capabilities: &VerificationCapabilities<'_>,
) -> CheckReport {
    build_report(
        pkg,
        expected_root,
        signature_verifier,
        capabilities,
        CheckPolicy::Integrity,
    )
}

/// Fail-closed release preflight without scientific-origin capabilities.
///
/// This helper intentionally supplies deny-all scientific capabilities, so it
/// can inventory release blockers but cannot admit certificate-class evidence.
/// Use [`check_for_release_with_capabilities`] for the actual release gate.
#[must_use]
pub fn check_release_preflight(
    pkg: &EvidencePackage,
    expected_root: ContentHash,
    verifier: &dyn SignatureVerifier,
) -> CheckReport {
    let preflight =
        preflight_content_address(pkg, Some(expected_root), CheckPolicy::ReleasePreflight);
    let mut report = match preflight {
        Ok(preflight) => {
            let mut report = build_report_from_preflight(
                pkg,
                Some(expected_root),
                Some(verifier),
                &VerificationCapabilities::deny_all(),
                CheckPolicy::ReleasePreflight,
                preflight,
            );
            append_release_findings(pkg, &mut report);
            report
        }
        Err(report) => *report,
    };
    report.findings.push(Finding {
        kind: "release-preflight-only",
        detail: "preflight inventories blockers but never grants release admission; invoke the \
                 capability-bearing release gate"
            .to_string(),
    });
    report.verdict = Verdict::Fail;
    report.decision_hash = checker_report_hash(&report);
    report
}

/// Actual release-admission gate with explicit source, anchor, falsifier,
/// derivation, waiver, and signature capabilities. It requires a non-empty
/// package, at least one scientifically admitted finite Verified or
/// authenticated Validated claim,
/// purpose-bound release signature, authenticated falsifiers for every
/// certificate-class claim, and an exact authenticated anchor for every
/// Validated claim.
#[must_use]
pub fn check_for_release_with_capabilities(
    pkg: &EvidencePackage,
    expected_root: ContentHash,
    verifier: &dyn SignatureVerifier,
    capabilities: &VerificationCapabilities<'_>,
) -> CheckReport {
    let preflight =
        match preflight_content_address(pkg, Some(expected_root), CheckPolicy::ReleaseAdmission) {
            Ok(preflight) => preflight,
            Err(report) => return *report,
        };
    if let Some(report) = release_shape_refusal(
        pkg,
        expected_root,
        preflight.merkle_root,
        &preflight.semantic_report,
    ) {
        return report;
    }
    let mut report = build_report_from_preflight(
        pkg,
        Some(expected_root),
        Some(verifier),
        capabilities,
        CheckPolicy::ReleaseAdmission,
        preflight,
    );
    append_release_findings(pkg, &mut report);
    report.decision_hash = checker_report_hash(&report);
    report
}

/// Strict JSON counterpart of [`check_release_preflight`]. This inventories
/// blockers under deny-all scientific capabilities and cannot grant release.
#[must_use]
pub fn check_json_release_preflight(
    text: &str,
    expected_root: ContentHash,
    verifier: &dyn SignatureVerifier,
) -> CheckReport {
    match EvidencePackage::from_json(text) {
        Ok(pkg) => check_release_preflight(&pkg, expected_root, verifier),
        Err(error) => parse_refusal(&error, CheckPolicy::ReleasePreflight, Some(expected_root)),
    }
}

/// Strict JSON counterpart of [`check_for_release_with_capabilities`].
/// Structural parse refusal and capability refusal both fail closed.
#[must_use]
pub fn check_json_for_release_with_capabilities(
    text: &str,
    expected_root: ContentHash,
    verifier: &dyn SignatureVerifier,
    capabilities: &VerificationCapabilities<'_>,
) -> CheckReport {
    match EvidencePackage::from_json(text) {
        Ok(pkg) => check_for_release_with_capabilities(&pkg, expected_root, verifier, capabilities),
        Err(e) => parse_refusal(&e, CheckPolicy::ReleaseAdmission, Some(expected_root)),
    }
}

fn parse_refusal(
    error: &ParseError,
    policy: CheckPolicy,
    expected_root: Option<ContentHash>,
) -> CheckReport {
    let mut report = CheckReport {
        verdict: Verdict::Fail,
        // Fail-closed sentinel: parsing refused, so there is no recomputed
        // root. The Fail verdict is authoritative; the zero bytes are only a
        // deterministic placeholder.
        merkle_root: ContentHash([0u8; 32]),
        breakdown: ColorBreakdown::default(),
        integrity_status: IntegrityStatus::Refused,
        semantic_report: SemanticReport::not_run(error.to_string()),
        origin_status: OriginStatus::NotRun,
        signature: SignatureStatus::Unsigned,
        receipt: None,
        findings: vec![Finding {
            kind: "parse-refused",
            detail: error.to_string(),
        }],
        policy,
        expected_root,
        decision_hash: ContentHash([0; 32]),
    };
    report.decision_hash = checker_report_hash(&report);
    report
}

fn scientific_evidence_required_finding() -> Finding {
    Finding {
        kind: "release-scientific-evidence-required",
        detail: "release admission requires at least one scientifically admitted finite \
                 Verified interval or authenticated Validated claim; vacuous infinite \
                 enclosures, preflight capability refusals, all-estimated packages, and \
                 all-waived packages remain publication/no-claim artifacts"
            .to_string(),
    }
}

fn release_signature_required_finding() -> Finding {
    Finding {
        kind: "release-signature-required",
        detail: "release admission requires a policy-authenticated detached release-approval \
                 signature bound to this checker protocol, expected content root, and exact \
                 scientific admission context and exact semantic-checker transcript"
            .to_string(),
    }
}

fn semantic_failure_finding(failure: &SemanticFailure) -> Finding {
    let kind = match failure.kind() {
        SemanticFailureKind::StructuralIntegrity => "semantic-structural-integrity",
        SemanticFailureKind::UnknownFamily => "semantic-unknown-family",
        SemanticFailureKind::UnsupportedVersion => "semantic-unsupported-version",
        SemanticFailureKind::MalformedPayload => "semantic-malformed-witness",
        SemanticFailureKind::ResourceLimit => "semantic-resource-limit",
        SemanticFailureKind::ClaimMismatch => "semantic-claim-mismatch",
        SemanticFailureKind::VerifierPanic => "semantic-verifier-panic",
    };
    let subject = match (failure.claim_index(), failure.claim_id()) {
        (Some(index), Some(id)) => format!("claim {index} ('{id}')"),
        (Some(index), None) => format!("claim {index}"),
        (None, Some(id)) => format!("claim '{id}'"),
        (None, None) => "package".to_string(),
    };
    let dispatch = match (failure.family(), failure.schema_version()) {
        (Some(family), Some(version)) => format!(" family '{family}' schema {version}"),
        (Some(family), None) => format!(" family '{family}'"),
        (None, Some(version)) => format!(" schema {version}"),
        (None, None) => String::new(),
    };
    Finding {
        kind,
        detail: format!("{subject}{dispatch}: {}", failure.detail()),
    }
}

fn release_shape_findings(pkg: &EvidencePackage) -> Option<Vec<Finding>> {
    if !pkg.is_structurally_inspectable_unverified() {
        return None;
    }
    let claims = pkg.declared_claims_unverified();
    let mut waiver_dependent = Vec::with_capacity(claims.len());
    for claim in claims {
        let depends = match claim.declared_origin_unverified() {
            fs_package::ClaimOrigin::AuthenticatedWaiver(_) => true,
            fs_package::ClaimOrigin::Derived => {
                claim.declared_receipt_unverified().is_none_or(|receipt| {
                    receipt
                        .parents
                        .iter()
                        .any(|parent| waiver_dependent.get(*parent).copied().unwrap_or(true))
                })
            }
            fs_package::ClaimOrigin::SourceCertificate { .. }
            | fs_package::ClaimOrigin::AnchoredSource { .. }
            | fs_package::ClaimOrigin::EstimatedSource { .. } => false,
        };
        waiver_dependent.push(depends);
    }

    let mut findings = Vec::new();
    if pkg.declared_claims_unverified().is_empty() {
        findings.push(Finding {
            kind: "release-empty-package",
            detail: "release admission requires at least one claim".to_string(),
        });
    }
    if !claims.iter().zip(&waiver_dependent).any(|(claim, waived)| {
        !waived && claim.declared_is_release_scientific_evidence_unverified()
    }) {
        findings.push(scientific_evidence_required_finding());
    }
    if pkg.signature.is_none() {
        findings.push(release_signature_required_finding());
    }
    for claim in claims {
        if claim.declared_requires_release_falsifier_unverified()
            && claim.declared_falsifiers_unverified().is_empty()
        {
            findings.push(Finding {
                kind: "release-falsifier-required",
                detail: format!(
                    "certificate-class claim '{}' cannot ship without an attached falsifier \
                     record",
                    claim.id()
                ),
            });
        }
        if claim.declared_requires_validated_anchor_unverified()
            && !claim.has_declared_matching_validated_anchor_unverified()
        {
            findings.push(Finding {
                kind: "release-anchor-required",
                detail: format!(
                    "validated claim '{}' cannot ship without a canonical content-hash anchor \
                     for its named dataset",
                    claim.id()
                ),
            });
        }
    }
    Some(findings)
}

fn push_unique_finding(report: &mut CheckReport, finding: Finding) {
    if !report
        .findings
        .iter()
        .any(|existing| existing.kind == finding.kind && existing.detail == finding.detail)
    {
        report.findings.push(finding);
    }
}

fn unverified_signature_status(pkg: &EvidencePackage) -> SignatureStatus {
    pkg.signature
        .as_ref()
        .map_or(SignatureStatus::Unsigned, |signature| {
            SignatureStatus::Unverified(signature.clone())
        })
}

fn release_shape_refusal(
    pkg: &EvidencePackage,
    expected_root: ContentHash,
    merkle_root: ContentHash,
    semantic_report: &SemanticReport,
) -> Option<CheckReport> {
    let findings = release_shape_findings(pkg)?;
    if findings.is_empty() {
        return None;
    }
    Some(callback_free_refusal(
        CheckPolicy::ReleaseAdmission,
        Some(expected_root),
        merkle_root,
        unverified_signature_status(pkg),
        findings,
        IntegrityStatus::Verified,
        semantic_report.clone(),
        OriginStatus::NotRun,
    ))
}

fn append_release_findings(pkg: &EvidencePackage, report: &mut CheckReport) {
    // Complete structural inspection is the bounded precondition for every raw
    // declaration scan. Oversized or malformed builders remain single-refusal
    // inputs and are never amplified.
    let Some(shape_findings) = release_shape_findings(pkg) else {
        report.verdict = Verdict::Fail;
        return;
    };
    for finding in shape_findings {
        push_unique_finding(report, finding);
    }
    let has_scientific_certificate = report.receipt.as_ref().is_some_and(|receipt| {
        pkg.declared_claims_unverified()
            .iter()
            .zip(receipt.admissions())
            .any(|(claim, admission)| {
                admission.class() == AdmissionClass::Scientific
                    && claim.declared_is_release_scientific_evidence_unverified()
            })
    });
    if !has_scientific_certificate {
        push_unique_finding(report, scientific_evidence_required_finding());
    }
    let receipt_admission_context = report
        .receipt
        .as_ref()
        .map(VerificationReceipt::release_admission_context);
    let release_signature = matches!(
        &report.signature,
        SignatureStatus::Authenticated(authenticated)
            if matches!(authenticated.purpose(), SignaturePurpose::ReleaseApproval {
                checker_protocol,
                expected_root,
                admission_context,
                semantic_context,
            } if checker_protocol == CHECKER_PROTOCOL_VERSION
                && Some(expected_root) == report.expected_root
                && Some(admission_context) == receipt_admission_context
                && semantic_context == report.semantic_report.context_hash())
    );
    if !release_signature {
        push_unique_finding(report, release_signature_required_finding());
    }
    if !report.findings.is_empty() {
        report.verdict = Verdict::Fail;
    }
}

fn callback_free_refusal(
    policy: CheckPolicy,
    expected_root: Option<ContentHash>,
    merkle_root: ContentHash,
    signature: SignatureStatus,
    findings: Vec<Finding>,
    integrity_status: IntegrityStatus,
    semantic_report: SemanticReport,
    origin_status: OriginStatus,
) -> CheckReport {
    let mut report = CheckReport {
        verdict: Verdict::Fail,
        merkle_root,
        breakdown: ColorBreakdown::default(),
        integrity_status,
        semantic_report,
        origin_status,
        signature,
        receipt: None,
        findings,
        policy,
        expected_root,
        decision_hash: ContentHash([0; 32]),
    };
    report.decision_hash = checker_report_hash(&report);
    report
}

struct CheckerPreflight {
    merkle_root: ContentHash,
    semantic_report: SemanticReport,
}

fn preflight_content_address(
    pkg: &EvidencePackage,
    expected_root: Option<ContentHash>,
    policy: CheckPolicy,
) -> Result<CheckerPreflight, Box<CheckReport>> {
    let merkle_root = pkg.try_merkle_root().map_err(|error| {
        Box::new(callback_free_refusal(
            policy,
            expected_root,
            ContentHash([0; 32]),
            SignatureStatus::Refused {
                reason: "package transport envelope refused",
            },
            vec![describe(&error)],
            IntegrityStatus::Refused,
            SemanticReport::not_run("content root unavailable before semantic verification"),
            OriginStatus::NotRun,
        ))
    })?;
    if let Some(expected) = expected_root
        && merkle_root != expected
    {
        return Err(Box::new(callback_free_refusal(
            policy,
            expected_root,
            merkle_root,
            unverified_signature_status(pkg),
            vec![Finding {
                kind: "content-address-mismatch",
                detail: format!("recomputed root {merkle_root} != expected {expected}"),
            }],
            IntegrityStatus::Refused,
            SemanticReport::not_run(
                "expected-root mismatch stopped structural and semantic verification",
            ),
            OriginStatus::NotRun,
        )));
    }

    let structural_root = match pkg.verify_structural_integrity() {
        Ok(root) => root,
        Err(error) => {
            return Err(Box::new(callback_free_refusal(
                policy,
                expected_root,
                merkle_root,
                unverified_signature_status(pkg),
                vec![describe(&error)],
                IntegrityStatus::Refused,
                SemanticReport::not_run(error.to_string()),
                OriginStatus::NotRun,
            )));
        }
    };
    if structural_root != merkle_root {
        return Err(Box::new(callback_free_refusal(
            policy,
            expected_root,
            merkle_root,
            unverified_signature_status(pkg),
            vec![Finding {
                kind: "structural-root-drift",
                detail: format!(
                    "transport root {merkle_root} differs from structural root {structural_root}"
                ),
            }],
            IntegrityStatus::Refused,
            SemanticReport::not_run("structural verification recomputed a different root"),
            OriginStatus::NotRun,
        )));
    }

    let semantic_report = semantic::verify_portable_semantics_after_integrity(pkg, structural_root);
    if matches!(
        semantic_report.status(),
        SemanticStatus::Refused | SemanticStatus::NotRun
    ) {
        let findings = semantic_report
            .failures()
            .iter()
            .map(semantic_failure_finding)
            .collect();
        return Err(Box::new(callback_free_refusal(
            policy,
            expected_root,
            merkle_root,
            unverified_signature_status(pkg),
            findings,
            IntegrityStatus::Verified,
            semantic_report,
            OriginStatus::NotRun,
        )));
    }

    Ok(CheckerPreflight {
        merkle_root,
        semantic_report,
    })
}

/// The full report builder. Missing origin capabilities fail closed only for
/// the origin kinds that require them. Any package-verification refusal yields
/// a zeroed breakdown, so unauthenticated bytes never retain a positive pie.
fn build_report(
    pkg: &EvidencePackage,
    expected_root: Option<ContentHash>,
    signature_verifier: Option<&dyn SignatureVerifier>,
    capabilities: &VerificationCapabilities<'_>,
    policy: CheckPolicy,
) -> CheckReport {
    // The caller-supplied root is a trust boundary, not merely another
    // finding. Reject an oversized or substituted package before dispatching
    // any injected capability: those callbacks may fetch attacker-selected
    // artifact addresses or perform other trusted side effects.
    let preflight = match preflight_content_address(pkg, expected_root, policy) {
        Ok(preflight) => preflight,
        Err(report) => return *report,
    };

    build_report_from_preflight(
        pkg,
        expected_root,
        signature_verifier,
        capabilities,
        policy,
        preflight,
    )
}

fn build_report_from_preflight(
    pkg: &EvidencePackage,
    expected_root: Option<ContentHash>,
    signature_verifier: Option<&dyn SignatureVerifier>,
    capabilities: &VerificationCapabilities<'_>,
    policy: CheckPolicy,
    preflight: CheckerPreflight,
) -> CheckReport {
    let CheckerPreflight {
        merkle_root,
        semantic_report,
    } = preflight;
    let mut findings = Vec::new();

    // Signature verification is now part of the package capability ledger.
    // The legacy separate argument remains as a compatible checker entry
    // point and overrides an optional signature capability in `capabilities`.
    let signatures = signature_verifier
        .map(|verifier| SignatureVerification {
            verifier,
            intent: match policy {
                CheckPolicy::Integrity => SignatureIntent::PackageRootAttestation,
                CheckPolicy::ReleasePreflight | CheckPolicy::ReleaseAdmission => {
                    SignatureIntent::ReleaseApproval {
                        checker_protocol: CHECKER_PROTOCOL_VERSION,
                        expected_root: expected_root.unwrap_or(ContentHash([0; 32])),
                        semantic_context: semantic_report.context_hash(),
                    }
                }
            },
        })
        .or(capabilities.signatures);
    let effective_capabilities = VerificationCapabilities {
        source_certificates: capabilities.source_certificates,
        anchored_sources: capabilities.anchored_sources,
        falsifiers: capabilities.falsifiers,
        derivations: capabilities.derivations,
        waivers: capabilities.waivers,
        signatures,
    };

    // Package-format structural and portable semantic checks already completed
    // callback-free. Only now may origin/signature capabilities run.
    let verified = pkg.verify_with(&effective_capabilities);
    let breakdown = match &verified {
        Ok(report) => *report.breakdown(),
        Err(e) => {
            findings.push(describe(e));
            // Invalid claims must not retain a normal-looking positive
            // evidence summary. The finding still identifies the exact
            // refusal; the pie fails closed to no admitted claims.
            ColorBreakdown::default()
        }
    };

    // 2. the magnitude budget must reconcile with its parts (the pie
    // is over error magnitudes, not claim counts — and it must not be
    // able to drift from the claims it summarizes).
    if let Ok(report) = &verified {
        let mb = report.magnitude_budget();
        if mb.quantified_total.to_bits() != (mb.verified_width + mb.estimated_dispersion).to_bits()
        {
            findings.push(Finding {
                kind: "magnitude-budget-drift",
                detail: "quantified total does not reconcile with its parts".to_string(),
            });
        }
    }

    // 3. The package verifier owns signature authentication so coverage and
    // checker reports consume the exact same decision.
    let signature = match &verified {
        Ok(report) => report.receipt().signature().clone(),
        Err(_) => match &pkg.signature {
            Some(signature) => SignatureStatus::Unverified(signature.clone()),
            None => SignatureStatus::Unsigned,
        },
    };
    let receipt = verified
        .as_ref()
        .ok()
        .map(|report| report.receipt().clone());
    let origin_status = match &verified {
        Ok(_)
        | Err(PackageError::SignatureRefused { .. } | PackageError::InvalidSignature { .. }) => {
            OriginStatus::Authenticated
        }
        Err(_) => OriginStatus::Refused,
    };

    let verdict = if findings.is_empty() {
        Verdict::Pass
    } else {
        Verdict::Fail
    };
    let mut report = CheckReport {
        verdict,
        merkle_root,
        breakdown,
        integrity_status: IntegrityStatus::Verified,
        semantic_report,
        origin_status,
        signature,
        receipt,
        findings,
        policy,
        expected_root,
        decision_hash: ContentHash([0; 32]),
    };
    report.decision_hash = checker_report_hash(&report);
    report
}

fn rejection_policy_suffix(policy_fingerprint: Option<ContentHash>) -> String {
    match policy_fingerprint {
        Some(fingerprint) => format!("; rejecting policy fingerprint {fingerprint}"),
        None => "; no policy decision was produced".to_string(),
    }
}

/// Translate a package error into a checker finding.
#[allow(clippy::too_many_lines)] // exhaustive one-to-one diagnostic mapping for every package refusal
fn describe(e: &PackageError) -> Finding {
    match e {
        PackageError::IncompleteProvenance { missing } => Finding {
            kind: "incomplete-provenance",
            detail: format!("package provenance is missing {missing}"),
        },
        PackageError::InvalidIdentity {
            claim,
            field,
            reason,
        } => Finding {
            kind: "invalid-identity",
            detail: match claim {
                Some(claim) => format!("claim '{claim}' has {reason} identity field {field}"),
                None => format!("package has {reason} identity field {field}"),
            },
        },
        PackageError::InvalidClaimId { index, id, reason } => Finding {
            kind: "invalid-claim-id",
            detail: format!("claim at index {index} has {reason} id {id:?}"),
        },
        PackageError::InvalidClaimStatement { claim, reason } => Finding {
            kind: "invalid-claim-statement",
            detail: format!("claim '{claim}' has a {reason} statement"),
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
        PackageError::UnsupportedColorAlgebra {
            claim,
            found,
            supported,
        } => Finding {
            kind: "unsupported-color-algebra",
            detail: format!(
                "claim '{claim}' uses color algebra version {found}; this checker supports \
                 version {supported}"
            ),
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
        PackageError::InvalidDerivationArtifact { claim } => Finding {
            kind: "invalid-derivation-artifact",
            detail: format!("claim '{claim}' has a non-canonical derivation artifact address"),
        },
        PackageError::DerivationRefused {
            claim,
            why,
            policy_fingerprint,
        } => Finding {
            kind: "derivation-refused",
            detail: format!(
                "claim '{claim}': derivation artifact refused — {why}{}",
                rejection_policy_suffix(*policy_fingerprint)
            ),
        },
        PackageError::InvalidOrigin { claim, why } => Finding {
            kind: "invalid-origin",
            detail: format!("claim '{claim}' has a malformed origin: {why}"),
        },
        PackageError::OriginMismatch { claim, origin } => Finding {
            kind: "origin-mismatch",
            detail: format!(
                "claim '{claim}': its {origin} origin cannot justify its color class — a raw \
                 color without a consistent origin is not evidence"
            ),
        },
        PackageError::SourceCertificateRefused {
            claim,
            producer,
            why,
            policy_fingerprint,
        } => Finding {
            kind: "source-certificate-refused",
            detail: format!(
                "claim '{claim}': source certificate from '{producer}' refused — {why}{}",
                rejection_policy_suffix(*policy_fingerprint)
            ),
        },
        PackageError::AnchoredSourceRefused {
            claim,
            dataset,
            why,
            policy_fingerprint,
        } => Finding {
            kind: "anchored-source-refused",
            detail: format!(
                "claim '{claim}': anchoring dataset '{dataset}' refused — {why}{}",
                rejection_policy_suffix(*policy_fingerprint)
            ),
        },
        PackageError::FalsifierRefused {
            claim,
            falsifier,
            why,
            policy_fingerprint,
        } => Finding {
            kind: "falsifier-refused",
            detail: format!(
                "claim '{claim}': falsifier '{falsifier}' refused — {why}{}",
                rejection_policy_suffix(*policy_fingerprint)
            ),
        },
        PackageError::WaiverRefused {
            claim,
            waiver,
            why,
            policy_fingerprint,
        } => Finding {
            kind: "waiver-refused",
            detail: format!(
                "claim '{claim}': waiver '{waiver}' refused — {why}{}",
                rejection_policy_suffix(*policy_fingerprint)
            ),
        },
        PackageError::PolicyFingerprintRefused {
            capability,
            why,
            previous,
            observed,
        } => Finding {
            kind: "policy-fingerprint-refused",
            detail: format!(
                "{capability} verification policy refused — {why}; previous {previous}, observed \
                 {observed}"
            ),
        },
        PackageError::SignatureRefused {
            why,
            policy_fingerprint,
        } => Finding {
            kind: "signature-invalid",
            detail: format!(
                "detached signature refused — {why}{}",
                rejection_policy_suffix(*policy_fingerprint)
            ),
        },
        PackageError::InvalidSignature { why } => Finding {
            kind: "invalid-signature",
            detail: format!("detached signature has invalid transport shape — {why}"),
        },
        PackageError::DuplicateWaiverId {
            waiver,
            first_claim,
            duplicate_claim,
        } => Finding {
            kind: "duplicate-waiver-id",
            detail: format!(
                "waiver '{waiver}' is reused by claims '{first_claim}' and '{duplicate_claim}'"
            ),
        },
        PackageError::InvalidWaiverTarget { index } => Finding {
            kind: "invalid-waiver-target",
            detail: format!("claim index {index} is absent or does not carry a waiver origin"),
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
        PackageError::InvalidSemanticWitness {
            claim,
            field,
            reason,
        } => Finding {
            kind: "invalid-semantic-witness",
            detail: format!("claim '{claim}' has invalid {field}: {reason}"),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_package::{Claim, Provenance};

    fn report_fixture() -> CheckReport {
        let package =
            EvidencePackage::new(Provenance::new("checker-test", "lock-test")).with_claim(
                Claim::estimated("estimate", "bounded estimate", "test-estimator", 1.0),
            );
        let mut report = check(&package);
        report.findings = vec![
            Finding {
                kind: "fixture-one",
                detail: "first fixture finding".to_string(),
            },
            Finding {
                kind: "fixture-two",
                detail: "second fixture finding".to_string(),
            },
        ];
        report.decision_hash = checker_report_hash(&report);
        report
    }

    #[test]
    fn checker_decision_identity_versions_and_transports_fail_closed() {
        assert_eq!(CHECKER_DECISION_IDENTITY_VERSION, 8);
        assert_eq!(
            CHECKER_DECISION_IDENTITY_VERSION,
            CHECKER_SUPPORTED_PACKAGE_FORMAT
        );
        assert_eq!(CHECKER_SUPPORTED_PACKAGE_FORMAT, fs_package::FORMAT_VERSION);
        assert_eq!(
            CHECKER_DECISION_IDENTITY_DOMAIN,
            "fs-package:v8:checker-decision"
        );
        assert_eq!(CHECKER_PROTOCOL_VERSION, 6);
        assert_ne!(CHECKER_PROTOCOL_VERSION, CHECKER_DECISION_IDENTITY_VERSION);

        let report = report_fixture();
        assert_eq!(
            CheckReport::admit_retained_decision_hash(
                CHECKER_DECISION_IDENTITY_VERSION,
                report.decision_hash().as_bytes(),
            ),
            Some(report.decision_hash())
        );
        for stale in [0, 7, 9, u32::MAX] {
            assert_eq!(
                CheckReport::admit_retained_decision_hash(stale, report.decision_hash().as_bytes(),),
                None
            );
        }
        for malformed in [&[0_u8; 31][..], &[0_u8; 33][..]] {
            assert_eq!(
                CheckReport::admit_retained_decision_hash(
                    CHECKER_DECISION_IDENTITY_VERSION,
                    malformed,
                ),
                None
            );
        }
        assert_ne!(
            checker_report_hash_with_protocol(&report, CHECKER_PROTOCOL_VERSION + 1),
            report.decision_hash()
        );
    }

    #[test]
    fn decision_hash_binds_every_checker_authority_field() {
        let report = report_fixture();
        assert!(report.validate_decision_hash());

        let authenticated_signature_hash = |signature: &str, purpose: SignaturePurpose| {
            let mut canonical = Vec::new();
            checker_decision_atom(&mut canonical, b"signature:authenticated");
            append_authenticated_signature_identity(&mut canonical, signature, purpose);
            hash_checker_decision(&canonical)
        };
        let package_root_purpose = SignaturePurpose::PackageRootAttestation;
        let release_purpose = SignaturePurpose::ReleaseApproval {
            checker_protocol: CHECKER_PROTOCOL_VERSION,
            expected_root: ContentHash([3; 32]),
            admission_context: ContentHash([4; 32]),
            semantic_context: ContentHash([5; 32]),
        };
        let authenticated_baseline =
            authenticated_signature_hash("authenticated-payload-a", package_root_purpose);
        assert_ne!(
            authenticated_baseline,
            authenticated_signature_hash("authenticated-payload-b", package_root_purpose),
            "signature payload must move while purpose stays fixed"
        );
        assert_ne!(
            authenticated_baseline,
            authenticated_signature_hash("authenticated-payload-a", release_purpose),
            "signature purpose must move while payload stays fixed"
        );

        macro_rules! assert_mutation_refused {
            ($mutation:expr) => {{
                let mut changed = report.clone();
                $mutation(&mut changed);
                assert!(!changed.validate_decision_hash());
            }};
        }

        assert_mutation_refused!(|changed: &mut CheckReport| changed.verdict = Verdict::Fail);
        assert_mutation_refused!(
            |changed: &mut CheckReport| changed.merkle_root = ContentHash([0; 32])
        );
        assert_mutation_refused!(|changed: &mut CheckReport| changed.breakdown.verified = 1);
        assert_mutation_refused!(|changed: &mut CheckReport| changed.breakdown.validated = 1);
        assert_mutation_refused!(|changed: &mut CheckReport| changed.breakdown.estimated += 1);
        assert_mutation_refused!(|changed: &mut CheckReport| changed.breakdown.waived = 1);
        assert_mutation_refused!(
            |changed: &mut CheckReport| changed.integrity_status = IntegrityStatus::Refused
        );
        assert_mutation_refused!(|changed: &mut CheckReport| changed.semantic_report =
            SemanticReport::not_run("mutated semantic transcript"));
        assert_mutation_refused!(
            |changed: &mut CheckReport| changed.origin_status = OriginStatus::NotRun
        );
        assert_mutation_refused!(|changed: &mut CheckReport| changed.signature =
            SignatureStatus::Refused {
                reason: "mutated test status"
            });
        assert_mutation_refused!(|changed: &mut CheckReport| changed.receipt = None);
        let replacement_receipt = {
            let replacement_package =
                EvidencePackage::new(Provenance::new("checker-test", "replacement-lock"))
                    .with_claim(Claim::estimated(
                        "replacement-estimate",
                        "replacement bounded estimate",
                        "replacement-estimator",
                        2.0,
                    ));
            check(&replacement_package)
                .receipt
                .expect("replacement fixture must retain a receipt")
        };
        assert_ne!(
            report
                .receipt
                .as_ref()
                .expect("fixture must retain a receipt")
                .receipt_hash(),
            replacement_receipt.receipt_hash()
        );
        assert_mutation_refused!(
            |changed: &mut CheckReport| changed.receipt = Some(replacement_receipt.clone())
        );
        assert_mutation_refused!(|changed: &mut CheckReport| changed.findings.push(Finding {
            kind: "mutated",
            detail: "mutated finding".to_string(),
        }));
        assert_mutation_refused!(|changed: &mut CheckReport| changed.findings.swap(0, 1));
        assert_mutation_refused!(|changed: &mut CheckReport| changed.findings[0].kind = "mutated");
        assert_mutation_refused!(|changed: &mut CheckReport| changed.findings[0].detail.push('x'));
        assert_mutation_refused!(
            |changed: &mut CheckReport| changed.policy = CheckPolicy::ReleasePreflight
        );
        assert_mutation_refused!(
            |changed: &mut CheckReport| changed.expected_root = Some(ContentHash([7; 32]))
        );
        let mut present_expected_root = report.clone();
        present_expected_root.expected_root = Some(ContentHash([7; 32]));
        present_expected_root.decision_hash = checker_report_hash(&present_expected_root);
        assert!(present_expected_root.validate_decision_hash());
        present_expected_root.expected_root = Some(ContentHash([9; 32]));
        assert!(
            !present_expected_root.validate_decision_hash(),
            "expected-root value must move while option presence stays fixed"
        );
        assert_mutation_refused!(
            |changed: &mut CheckReport| changed.decision_hash = ContentHash([8; 32])
        );
    }

    #[test]
    fn budget_pie_uses_wide_arithmetic_for_extreme_counts() {
        let mut report = report_fixture();
        report.breakdown = ColorBreakdown {
            verified: usize::MAX,
            validated: usize::MAX,
            estimated: usize::MAX,
            waived: usize::MAX,
        };
        let rendered = report.render_pie();
        assert!(rendered.contains("verified"));
        assert!(rendered.contains("waived"));
    }
}
