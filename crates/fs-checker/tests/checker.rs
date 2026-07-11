//! Battery for the standalone evidence-package checker (addendum Proposal 12).
//! Covers a clean pass, completeness-failure findings, content-address
//! (Merkle) tamper detection, signature-presence reporting, budget-pie
//! rendering (including the empty case), the protocol version, and
//! determinism. The checker uses only the package format — no solver.

use fs_checker::{
    AnchoredSourceRequest, AnchoredSourceVerifier, CHECKER_PROTOCOL_VERSION, ColorBreakdown,
    ContentHash, FalsifierRequest, FalsifierVerifier, SignaturePurpose, SignatureRequest,
    SignatureStatus, SourceCertificateRequest, SourceCertificateVerifier, Verdict,
    VerificationCapabilities, VerificationDecision, WaiverGrant, WaiverVerifier, check,
    check_against_root, check_for_release_with_capabilities, check_json,
    check_json_for_release_with_capabilities, check_json_release_preflight,
    check_json_with_capabilities, check_release_preflight, check_with_capabilities,
};
use fs_evidence::{Color, ValidityDomain};
use fs_package::{Claim, EvidencePackage, FalsifierRecord, Provenance};

const ARTIFACT_HASH: &str = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

/// A deliberately corrupted content root (one byte flipped): the v5
/// 32-byte replacement for the old `root ^ 0xdead` tamper idiom.
fn flip(root: ContentHash) -> ContentHash {
    let mut bytes = *root.as_bytes();
    bytes[0] ^= 0xde;
    ContentHash(bytes)
}

fn prov() -> Provenance {
    Provenance::new("commit-abc", "lock-def")
}

fn package_root(package: &EvidencePackage) -> ContentHash {
    package
        .try_merkle_root()
        .expect("bounded checker fixture has a content root")
}

fn package_json(package: &EvidencePackage) -> String {
    package
        .to_json()
        .expect("bounded checker fixture serializes")
}

fn verified(id: &str) -> Claim {
    Claim::from_certificate(id, "ok", -1.0, 1.0, "test-solver/cert", ARTIFACT_HASH)
}
fn estimated(id: &str) -> Claim {
    Claim::estimated(id, "maybe", "surrogate", 2.0)
}
fn validated(id: &str, regime: ValidityDomain) -> Claim {
    Claim::anchored(id, "matches", regime, "wt-2026", ARTIFACT_HASH)
}
fn good_regime() -> ValidityDomain {
    ValidityDomain::unconstrained().with("Re", 1e5, 3e5)
}

struct ExactSourceVerifier<'a> {
    claim_id: &'a str,
}

struct AlternatePolicySourceVerifier<'a> {
    claim_id: &'a str,
}

struct ExactAnchorVerifier;
struct ExactFalsifierVerifier;

fn policy_fingerprint(label: &str) -> ContentHash {
    let byte = match label {
        "exact-source-verifier" => 0x11,
        "exact-anchor-verifier" => 0x22,
        "release-verifier" => 0x33,
        "exact-waiver-verifier" => 0x44,
        "mac-verifier" => 0x55,
        "exact-falsifier-verifier" => 0x66,
        _ => 0xff,
    };
    ContentHash([byte; 32])
}

impl SourceCertificateVerifier for ExactSourceVerifier<'_> {
    fn verify(&self, request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let accepted = request.package_provenance == &prov()
            && request.claim_index == 0
            && request.claim_id == self.claim_id
            && request.statement == "ok"
            && request.lo.to_bits() == (-1.0f64).to_bits()
            && request.hi.to_bits() == 1.0f64.to_bits()
            && request.producer == "test-solver/cert"
            && request.certificate_hash.to_hex() == ARTIFACT_HASH;
        VerificationDecision {
            accepted,
            policy_fingerprint: policy_fingerprint("exact-source-verifier"),
        }
    }
}

impl SourceCertificateVerifier for AlternatePolicySourceVerifier<'_> {
    fn verify(&self, request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let accepted = ExactSourceVerifier {
            claim_id: self.claim_id,
        }
        .verify(request)
        .accepted;
        VerificationDecision {
            accepted,
            policy_fingerprint: ContentHash([0x77; 32]),
        }
    }
}

impl AnchoredSourceVerifier for ExactAnchorVerifier {
    fn verify(&self, request: &AnchoredSourceRequest<'_>) -> VerificationDecision {
        let accepted = request.package_provenance == &prov()
            && request.statement == "matches"
            && request.dataset_id == "wt-2026"
            && request.content_hash.to_hex() == ARTIFACT_HASH
            && request.regime == &good_regime()
            && matches!(
                (request.claim_index, request.claim_id),
                (1, "c2" | "validated") | (0, "v")
            );
        VerificationDecision {
            accepted,
            policy_fingerprint: policy_fingerprint("exact-anchor-verifier"),
        }
    }
}

impl FalsifierVerifier for ExactFalsifierVerifier {
    fn verify(&self, request: &FalsifierRequest<'_>) -> VerificationDecision {
        VerificationDecision {
            accepted: request.artifact_hash.to_hex() == ARTIFACT_HASH,
            policy_fingerprint: policy_fingerprint("exact-falsifier-verifier"),
        }
    }
}

static EXACT_ANCHOR_VERIFIER: ExactAnchorVerifier = ExactAnchorVerifier;
static EXACT_FALSIFIER_VERIFIER: ExactFalsifierVerifier = ExactFalsifierVerifier;

fn source_capabilities(verifier: &dyn SourceCertificateVerifier) -> VerificationCapabilities<'_> {
    VerificationCapabilities::deny_all()
        .with_source_certificates(verifier)
        .with_anchored_sources(&EXACT_ANCHOR_VERIFIER)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER)
}

#[test]
fn a_valid_package_passes_with_no_findings() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(validated("c2", good_regime()))
        .with_claim(estimated("c3"));
    let source_verifier = ExactSourceVerifier { claim_id: "c1" };
    let capabilities = source_capabilities(&source_verifier);
    let report = check_with_capabilities(&pkg, None, None, &capabilities);
    assert!(report.passed());
    assert!(!report.release_admitted());
    assert!(report.validate_decision_hash());
    assert_eq!(report.verdict(), Verdict::Pass);
    assert!(report.findings().is_empty());
    assert_eq!(report.merkle_root(), package_root(&pkg));
    assert_eq!(report.breakdown().verified, 1);
    assert_eq!(report.breakdown().validated, 1);
    assert_eq!(report.breakdown().estimated, 1);
    let receipt = report.receipt().expect("successful check retains receipt");
    assert_eq!(receipt.package_root(), package_root(&pkg));
    assert_eq!(receipt.admissions().len(), 3);
    assert_eq!(
        receipt.policy_fingerprints().source_certificates(),
        Some(policy_fingerprint("exact-source-verifier"))
    );
    assert_eq!(
        receipt.policy_fingerprints().anchored_sources(),
        Some(policy_fingerprint("exact-anchor-verifier"))
    );
}

#[test]
fn an_incomplete_validated_claim_fails_the_check() {
    // unconstrained regime = missing regime tag.
    let pkg =
        EvidencePackage::new(prov()).with_claim(validated("v", ValidityDomain::unconstrained()));
    let report = check(&pkg);
    assert!(!report.passed());
    assert_eq!(report.verdict(), Verdict::Fail);
    assert_eq!(report.findings().len(), 1);
    assert_eq!(report.findings()[0].kind, "incomplete-validated-claim");
}

#[test]
fn a_semantically_empty_falsifier_record_fails_the_check() {
    let pkg =
        EvidencePackage::new(prov()).with_claim(verified("v").with_falsifier(FalsifierRecord {
            name: " ".to_string(),
            attempts: 0,
            refuted: false,
            detail: " ".to_string(),
            artifact_hash: ARTIFACT_HASH.to_string(),
        }));
    let report = check(&pkg);
    assert!(!report.passed());
    assert_eq!(report.findings()[0].kind, "invalid-falsifier-record");
    assert_eq!(report.breakdown(), &ColorBreakdown::default());
    assert_eq!(report.render_pie(), "budget pie: no claims");
}

#[test]
fn placeholder_claim_and_falsifier_text_fail_the_check() {
    let placeholder_statement = EvidencePackage::new(prov()).with_claim(Claim::from_certificate(
        "claim",
        "TODO",
        0.0,
        1.0,
        "test-solver/cert",
        ARTIFACT_HASH,
    ));
    let report = check(&placeholder_statement);
    assert!(!report.passed());
    assert_eq!(report.findings()[0].kind, "invalid-claim-statement");

    let placeholder_falsifier = EvidencePackage::new(prov()).with_claim(
        verified("claim").with_falsifier(FalsifierRecord {
            name: "independent-probe".to_string(),
            attempts: 1,
            refuted: false,
            detail: "placeholder".to_string(),
            artifact_hash: ARTIFACT_HASH.to_string(),
        }),
    );
    let report = check(&placeholder_falsifier);
    assert!(!report.passed());
    assert_eq!(report.findings()[0].kind, "invalid-falsifier-record");
}

#[test]
fn content_address_mismatch_is_caught() {
    let pkg = EvidencePackage::new(prov()).with_claim(estimated("c1"));
    let real_root = package_root(&pkg);
    // the right root passes.
    assert!(check_against_root(&pkg, real_root).passed());
    // a wrong expected root (a tampered/substituted package) fails.
    let report = check_against_root(&pkg, flip(real_root));
    assert!(!report.passed());
    assert!(
        report
            .findings()
            .iter()
            .any(|f| f.kind == "content-address-mismatch")
    );
}

#[test]
fn content_address_mismatch_catches_provenance_tamper() {
    let pkg = EvidencePackage::new(prov()).with_claim(estimated("c1"));
    let root = package_root(&pkg);
    let tampered = EvidencePackage::new(Provenance::new("commit-evil", "lock-def"))
        .with_claim(estimated("c1"));

    let report = check_against_root(&tampered, root);
    assert!(!report.passed());
    assert!(
        report
            .findings()
            .iter()
            .any(|f| f.kind == "content-address-mismatch")
    );
}

#[test]
fn signature_presence_is_reported() {
    let unsigned = EvidencePackage::new(prov()).with_claim(estimated("e1"));
    assert_eq!(check(&unsigned).signature(), &SignatureStatus::Unsigned);
    let signed = unsigned.signed("ed25519:cafe");
    assert_eq!(
        check(&signed).signature(),
        &SignatureStatus::Unverified("ed25519:cafe".to_string())
    );
}

#[test]
fn the_budget_pie_renders_deterministically() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(estimated("c2"))
        .with_claim(estimated("c3"));
    let source_verifier = ExactSourceVerifier { claim_id: "c1" };
    let capabilities = source_capabilities(&source_verifier);
    let pie = check_with_capabilities(&pkg, None, None, &capabilities).render_pie();
    assert_eq!(
        pie,
        check_with_capabilities(&pkg, None, None, &capabilities).render_pie()
    );
    assert!(pie.contains("budget pie (3 claims)"));
    assert!(pie.contains("verified") && pie.contains("estimated"));
    assert!(pie.contains('#') && pie.contains('.'));
}

#[test]
fn the_budget_pie_handles_an_empty_package() {
    let pkg = EvidencePackage::new(prov());
    // an empty package still verifies (vacuously) and renders a no-claims pie.
    let report = check(&pkg);
    assert!(report.passed());
    assert_eq!(report.render_pie(), "budget pie: no claims");
}

struct ReleaseVerifier;

impl fs_checker::SignatureVerifier for ReleaseVerifier {
    fn verify(&self, request: &SignatureRequest<'_>) -> VerificationDecision {
        VerificationDecision {
            accepted: request.signature == format!("release-test:{}", request.subject_hash())
                && matches!(
                    request.purpose,
                    SignaturePurpose::ReleaseApproval {
                        checker_protocol: 4,
                        expected_root,
                        admission_context,
                    } if expected_root == request.package_root
                        && admission_context != ContentHash([0; 32])
                ),
            policy_fingerprint: policy_fingerprint("release-verifier"),
        }
    }
}

fn signed_for_release(
    package: EvidencePackage,
    capabilities: &VerificationCapabilities<'_>,
) -> EvidencePackage {
    let root = package_root(&package);
    let unsigned = package
        .verify_with(capabilities)
        .expect("unsigned release subject verifies");
    let purpose = SignaturePurpose::ReleaseApproval {
        checker_protocol: CHECKER_PROTOCOL_VERSION,
        expected_root: root,
        admission_context: unsigned.receipt().release_admission_context(),
    };
    package.signed(format!(
        "release-test:{}",
        fs_checker::signature_subject_hash(root, purpose)
    ))
}

fn passed_falsifier() -> FalsifierRecord {
    FalsifierRecord {
        name: "independent-interval-probe".to_string(),
        attempts: 64,
        refuted: false,
        detail: "64 boundary-biased probes found no violation".to_string(),
        artifact_hash: ARTIFACT_HASH.to_string(),
    }
}

fn assert_capability_refusal(report: &fs_checker::CheckReport, kind: &str) {
    assert!(!report.passed(), "capability refusal unexpectedly passed");
    assert_eq!(
        report.breakdown(),
        &ColorBreakdown::default(),
        "refused origin retained a positive evidence breakdown"
    );
    assert!(
        report.findings().iter().any(|finding| finding.kind == kind),
        "missing {kind} finding: {:?}",
        report.findings()
    );
    assert!(
        report.receipt().is_none(),
        "capability-refused checks must not expose a partial verification receipt"
    );
}

fn fixture_waiver_mac(message: &[u8]) -> String {
    let mut state = 0xcbf2_9ce4_8422_2325u64;
    for byte in message {
        state ^= u64::from(*byte);
        state = state.wrapping_mul(0x0000_0100_0000_01b3);
    }
    format!("checker-fixture:{state:016x}")
}

struct ExactWaiverVerifier;

impl WaiverVerifier for ExactWaiverVerifier {
    fn verify(&self, mac: &str, message: &[u8]) -> VerificationDecision {
        VerificationDecision {
            accepted: mac == fixture_waiver_mac(message),
            policy_fingerprint: policy_fingerprint("exact-waiver-verifier"),
        }
    }
}

fn waived_fixture() -> EvidencePackage {
    let pending = EvidencePackage::new(prov()).with_claim(
        Claim::waived(
            "waived",
            "authorized interval",
            Color::Verified { lo: -1.0, hi: 1.0 },
            WaiverGrant {
                waiver_id: "checker-waiver-2026".to_string(),
                expiry_day: 300,
                mac: "pending-authenticator".to_string(),
            },
        )
        .with_falsifier(passed_falsifier()),
    );
    let message = pending.waiver_message(0).expect("waiver target");
    pending
        .with_waiver_mac(0, fixture_waiver_mac(&message))
        .expect("install exact fixture authenticator")
}

#[test]
fn source_certificates_are_capability_gated_across_every_entry_path() {
    let unsigned = EvidencePackage::new(prov())
        .with_claim(verified("source").with_falsifier(passed_falsifier()));
    let root = package_root(&unsigned);
    let exact = ExactSourceVerifier { claim_id: "source" };
    let capabilities = source_capabilities(&exact);
    let signed = signed_for_release(unsigned.clone(), &capabilities);

    for report in [
        check(&unsigned),
        check_json(&package_json(&unsigned), Some(root), None),
        check_release_preflight(&signed, root, &ReleaseVerifier),
        check_json_release_preflight(&package_json(&signed), root, &ReleaseVerifier),
    ] {
        assert_capability_refusal(&report, "source-certificate-refused");
    }

    for report in [
        check_with_capabilities(&unsigned, Some(root), None, &capabilities),
        check_json_with_capabilities(&package_json(&unsigned), Some(root), None, &capabilities),
    ] {
        assert!(report.passed(), "{:?}", report.findings());
        assert_eq!(report.signature(), &SignatureStatus::Unsigned);
        assert_eq!(report.breakdown().verified, 1);
    }
    for report in [
        check_for_release_with_capabilities(&signed, root, &ReleaseVerifier, &capabilities),
        check_json_for_release_with_capabilities(
            &package_json(&signed),
            root,
            &ReleaseVerifier,
            &capabilities,
        ),
    ] {
        assert!(report.passed(), "{:?}", report.findings());
        assert!(matches!(
            report.signature(),
            SignatureStatus::Authenticated(authenticated)
                if matches!(authenticated.purpose(), SignaturePurpose::ReleaseApproval { .. })
        ));
    }

    let wrong_subject = ExactSourceVerifier {
        claim_id: "different-source",
    };
    let wrong_capabilities = source_capabilities(&wrong_subject);
    let rejected = check_with_capabilities(&unsigned, None, None, &wrong_capabilities);
    assert_capability_refusal(&rejected, "source-certificate-refused");
    assert!(
        rejected.findings()[0]
            .detail
            .contains(&policy_fingerprint("exact-source-verifier").to_string()),
        "rejected atomic policy identity must survive into the checker decision"
    );
}

#[test]
fn waivers_are_capability_gated_across_every_entry_path() {
    let unsigned = waived_fixture();
    let root = package_root(&unsigned);
    let waiver_verifier = ExactWaiverVerifier;
    let capabilities = VerificationCapabilities::deny_all()
        .with_waivers(&waiver_verifier, 250)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER);
    let signed = signed_for_release(unsigned.clone(), &capabilities);

    for report in [
        check(&unsigned),
        check_json(&package_json(&unsigned), Some(root), None),
        check_release_preflight(&signed, root, &ReleaseVerifier),
        check_json_release_preflight(&package_json(&signed), root, &ReleaseVerifier),
    ] {
        assert_capability_refusal(&report, "waiver-refused");
    }

    for report in [
        check_with_capabilities(&unsigned, Some(root), None, &capabilities),
        check_json_with_capabilities(&package_json(&unsigned), Some(root), None, &capabilities),
    ] {
        assert!(report.passed(), "{:?}", report.findings());
        assert_eq!(report.signature(), &SignatureStatus::Unsigned);
        assert_eq!(report.breakdown().verified, 0);
        assert_eq!(report.breakdown().waived, 1);
    }
    for report in [
        check_for_release_with_capabilities(&signed, root, &ReleaseVerifier, &capabilities),
        check_json_for_release_with_capabilities(
            &package_json(&signed),
            root,
            &ReleaseVerifier,
            &capabilities,
        ),
    ] {
        assert!(
            !report.passed(),
            "all-waived packages cannot be release evidence"
        );
        assert!(
            report
                .findings()
                .iter()
                .any(|finding| finding.kind == "release-scientific-evidence-required"),
            "missing scientific-evidence refusal: {:?}",
            report.findings()
        );
        assert!(matches!(
            report.signature(),
            SignatureStatus::Authenticated(authenticated)
                if matches!(authenticated.purpose(), SignaturePurpose::ReleaseApproval { .. })
        ));
    }

    let expired = VerificationCapabilities::deny_all()
        .with_waivers(&waiver_verifier, 301)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER);
    assert_capability_refusal(
        &check_with_capabilities(&unsigned, None, None, &expired),
        "waiver-refused",
    );
}

#[test]
fn release_gate_requires_certificate_obligations() {
    let unsigned = EvidencePackage::new(prov())
        .with_claim(verified("verified").with_falsifier(passed_falsifier()))
        .with_claim(validated("validated", good_regime()).with_falsifier(passed_falsifier()))
        .with_claim(estimated("honest-estimate"));
    let source_verifier = ExactSourceVerifier {
        claim_id: "verified",
    };
    let capabilities = source_capabilities(&source_verifier);
    let pkg = signed_for_release(unsigned, &capabilities);
    let root = package_root(&pkg);
    let preflight = check_release_preflight(&pkg, root, &ReleaseVerifier);
    assert!(!preflight.release_admitted());
    assert_eq!(
        preflight.policy(),
        fs_checker::CheckPolicy::ReleasePreflight
    );
    assert!(
        preflight
            .findings()
            .iter()
            .any(|finding| finding.kind == "release-preflight-only")
    );
    let report = check_for_release_with_capabilities(&pkg, root, &ReleaseVerifier, &capabilities);
    assert!(report.passed(), "{:?}", report.findings());
    assert!(report.release_admitted());
    assert!(report.validate_decision_hash());
    assert_eq!(report.policy(), fs_checker::CheckPolicy::ReleaseAdmission);
    assert_ne!(preflight.decision_hash(), report.decision_hash());
    assert!(matches!(
        report.signature(),
        SignatureStatus::Authenticated(authenticated)
            if matches!(authenticated.purpose(), SignaturePurpose::ReleaseApproval { .. })
    ));
    assert!(
        check_json_for_release_with_capabilities(
            &package_json(&pkg),
            root,
            &ReleaseVerifier,
            &capabilities,
        )
        .passed()
    );
}

#[test]
fn release_gate_rejects_all_estimated_even_with_valid_signature_and_root() {
    let package = signed_for_release(
        EvidencePackage::new(prov()).with_claim(estimated("estimate-only")),
        &VerificationCapabilities::deny_all(),
    );
    let root = package_root(&package);
    let report = check_for_release_with_capabilities(
        &package,
        root,
        &ReleaseVerifier,
        &VerificationCapabilities::deny_all(),
    );
    assert!(!report.passed());
    assert!(report.validate_decision_hash());
    assert!(matches!(
        report.signature(),
        SignatureStatus::Authenticated(authenticated)
            if matches!(authenticated.purpose(), SignaturePurpose::ReleaseApproval { .. })
    ));
    assert!(
        report
            .findings()
            .iter()
            .any(|finding| finding.kind == "release-scientific-evidence-required")
    );
}

#[test]
fn release_gate_refuses_package_root_attestation_substitution() {
    struct RootAttestationVerifier;

    impl fs_checker::SignatureVerifier for RootAttestationVerifier {
        fn verify(&self, request: &SignatureRequest<'_>) -> VerificationDecision {
            VerificationDecision {
                accepted: request.signature == format!("integrity-test:{}", request.subject_hash())
                    && request.purpose == SignaturePurpose::PackageRootAttestation,
                policy_fingerprint: ContentHash([0x88; 32]),
            }
        }
    }

    let unsigned = EvidencePackage::new(prov()).with_claim(estimated("integrity-only"));
    let root = package_root(&unsigned);
    let attestation_subject =
        fs_checker::signature_subject_hash(root, SignaturePurpose::PackageRootAttestation);
    let package = unsigned.signed(format!("integrity-test:{attestation_subject}"));
    let report = check_for_release_with_capabilities(
        &package,
        root,
        &RootAttestationVerifier,
        &VerificationCapabilities::deny_all(),
    );
    assert_capability_refusal(&report, "signature-invalid");
}

#[test]
fn release_approval_refuses_policy_or_waiver_clock_replay() {
    let source_package = EvidencePackage::new(prov())
        .with_claim(verified("source").with_falsifier(passed_falsifier()));
    let primary_source = ExactSourceVerifier { claim_id: "source" };
    let primary_capabilities = source_capabilities(&primary_source);
    let signed_source = signed_for_release(source_package.clone(), &primary_capabilities);

    let alternate_source = AlternatePolicySourceVerifier { claim_id: "source" };
    let alternate_capabilities = source_capabilities(&alternate_source);
    let primary_context = source_package
        .verify_with(&primary_capabilities)
        .expect("primary source policy admits the unsigned subject")
        .receipt()
        .release_admission_context();
    let alternate_context = source_package
        .verify_with(&alternate_capabilities)
        .expect("alternate source policy admits the same unsigned subject")
        .receipt()
        .release_admission_context();
    assert_ne!(
        primary_context, alternate_context,
        "the release subject must bind the scientific policy fingerprint"
    );
    let replayed_source = check_for_release_with_capabilities(
        &signed_source,
        package_root(&signed_source),
        &ReleaseVerifier,
        &alternate_capabilities,
    );
    assert_capability_refusal(&replayed_source, "signature-invalid");

    let waiver_package = waived_fixture();
    let waiver_verifier = ExactWaiverVerifier;
    let day_250 = VerificationCapabilities::deny_all()
        .with_waivers(&waiver_verifier, 250)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER);
    let day_249 = VerificationCapabilities::deny_all()
        .with_waivers(&waiver_verifier, 249)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER);
    let signed_waiver = signed_for_release(waiver_package.clone(), &day_250);
    let context_250 = waiver_package
        .verify_with(&day_250)
        .expect("waiver is valid on signing day")
        .receipt()
        .release_admission_context();
    let context_249 = waiver_package
        .verify_with(&day_249)
        .expect("waiver is also valid on replay day")
        .receipt()
        .release_admission_context();
    assert_ne!(
        context_250, context_249,
        "the release subject must bind the explicit waiver clock"
    );
    let replayed_waiver = check_for_release_with_capabilities(
        &signed_waiver,
        package_root(&signed_waiver),
        &ReleaseVerifier,
        &day_249,
    );
    assert_capability_refusal(&replayed_waiver, "signature-invalid");
}

#[test]
fn release_gate_refuses_vacuous_or_unpaired_packages() {
    let empty = signed_for_release(
        EvidencePackage::new(prov()),
        &VerificationCapabilities::deny_all(),
    );
    assert!(
        check(&empty).passed(),
        "ordinary integrity check stays vacuous"
    );
    let report = check_release_preflight(&empty, package_root(&empty), &ReleaseVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings()
            .iter()
            .any(|finding| finding.kind == "release-empty-package")
    );

    let source_verifier = ExactSourceVerifier { claim_id: "v" };
    let capabilities = source_capabilities(&source_verifier);
    let unpaired = signed_for_release(
        EvidencePackage::new(prov()).with_claim(verified("v")),
        &capabilities,
    );
    let report = check_for_release_with_capabilities(
        &unpaired,
        package_root(&unpaired),
        &ReleaseVerifier,
        &capabilities,
    );
    assert!(!report.passed());
    assert!(
        report
            .findings()
            .iter()
            .any(|finding| finding.kind == "release-falsifier-required")
    );
    let report = check_json_for_release_with_capabilities(
        &package_json(&unpaired),
        package_root(&unpaired),
        &ReleaseVerifier,
        &capabilities,
    );
    assert!(!report.passed(), "JSON must not bypass release policy");
    assert!(
        report
            .findings()
            .iter()
            .any(|finding| finding.kind == "release-falsifier-required")
    );
}

#[test]
fn release_preflight_does_not_amplify_rejected_oversized_builders() {
    let oversized = EvidencePackage::new(prov())
        .with_claim(estimated("bounded"))
        .signed("s".repeat(fs_package::MAX_JSON_STRING_BYTES + 1));
    let report = check_release_preflight(&oversized, ContentHash([0; 32]), &ReleaseVerifier);
    assert!(!report.passed());
    assert!(report.validate_decision_hash());
    assert_eq!(report.findings().len(), 2);
    assert_eq!(report.findings()[0].kind, "transport-limit");
    assert_eq!(report.findings()[1].kind, "release-preflight-only");
    assert!(report.receipt().is_none());
    assert!(matches!(
        report.signature(),
        SignatureStatus::Refused {
            reason: "package transport envelope refused"
        }
    ));
}

#[test]
fn release_gate_requires_matching_anchor_signature_and_root() {
    // Schema v6: the sealed `anchored` constructor attaches the matching
    // anchor, so an in-memory validated-without-anchor package is
    // unconstructible. The release anchor gate is now exercised through
    // the PARSE path: strip the matching anchor from the transported
    // JSON and the recomputed root refuses before the gate is even
    // reached — the transported form cannot lose its anchor silently.
    let anchor_capabilities = VerificationCapabilities::deny_all()
        .with_anchored_sources(&EXACT_ANCHOR_VERIFIER)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER);
    let anchored_pkg = signed_for_release(
        EvidencePackage::new(prov())
            .with_claim(validated("v", good_regime()).with_falsifier(passed_falsifier())),
        &anchor_capabilities,
    );
    let json = package_json(&anchored_pkg);
    let stripped = json.replacen(
        "{\"dataset_id\":\"wt-2026\",\"content_hash\"",
        "{\"dataset_id\":\"different-dataset\",\"content_hash\"",
        1,
    );
    assert!(
        fs_checker::EvidencePackage::from_json(&stripped).is_err(),
        "anchor tamper breaks the content address at parse"
    );
    assert_capability_refusal(
        &check_release_preflight(&anchored_pkg, package_root(&anchored_pkg), &ReleaseVerifier),
        "anchored-source-refused",
    );
    let report = check_for_release_with_capabilities(
        &anchored_pkg,
        package_root(&anchored_pkg),
        &ReleaseVerifier,
        &anchor_capabilities,
    );
    assert!(report.passed(), "{:?}", report.findings());
    assert!(
        !report
            .findings()
            .iter()
            .any(|finding| finding.kind == "release-anchor-required")
    );

    let unsigned = EvidencePackage::new(prov()).with_claim(estimated("v"));
    let report = check_release_preflight(&unsigned, package_root(&unsigned), &ReleaseVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings()
            .iter()
            .any(|finding| finding.kind == "release-signature-required")
    );

    let signed = signed_for_release(unsigned, &VerificationCapabilities::deny_all());
    let report = check_release_preflight(&signed, flip(package_root(&signed)), &ReleaseVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings()
            .iter()
            .any(|finding| finding.kind == "content-address-mismatch")
    );
}

#[test]
fn the_checker_advertises_its_protocol_version() {
    assert_eq!(CHECKER_PROTOCOL_VERSION, 4);
    assert_eq!(fs_checker::CHECKER_SUPPORTED_PACKAGE_FORMAT, 6);
    assert_eq!(
        fs_checker::CHECKER_SUPPORTED_PACKAGE_FORMAT,
        fs_package::FORMAT_VERSION
    );
}

#[test]
fn checking_is_deterministic() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(estimated("c1"))
        .with_claim(estimated("c2"));
    assert_eq!(check(&pkg), check(&pkg));
}

/// qmao.6.1 — the third-party JSON path: parse-refused inputs never
/// pass; signature validity is asserted only through a capability over
/// the canonical root-attestation subject; tamper anywhere fails.
#[test]
fn checker_json_path_and_signature_capability() {
    use fs_checker::{NoSignatureVerifier, SignatureVerifier, check_json, check_with};
    struct MacVerifier;
    fn mac(subject: ContentHash) -> String {
        format!("test-key/{subject}")
    }
    impl SignatureVerifier for MacVerifier {
        fn verify(&self, request: &SignatureRequest<'_>) -> VerificationDecision {
            VerificationDecision {
                accepted: request.signature == mac(request.subject_hash())
                    && request.purpose == SignaturePurpose::PackageRootAttestation,
                policy_fingerprint: policy_fingerprint("mac-verifier"),
            }
        }
    }
    let base = EvidencePackage::new(Provenance::new("v1.0", "lock:abc"))
        .with_claim(Claim::estimated("c1", "bounded", "surrogate", 1.0));
    let root = package_root(&base);
    let subject =
        fs_checker::signature_subject_hash(root, SignaturePurpose::PackageRootAttestation);
    let pkg = base.signed(mac(subject));
    // Valid signature via the capability.
    let report = check_with(&pkg, Some(root), &MacVerifier);
    assert!(report.passed(), "{:?}", report.findings());
    assert!(matches!(
        report.signature(),
        fs_checker::SignatureStatus::Authenticated(authenticated)
            if authenticated.purpose() == SignaturePurpose::PackageRootAttestation
    ));
    // The no-crypto default cannot assert validity (fail closed: the
    // signature stays Unverified and a finding is raised).
    let report = check_with(&pkg, Some(root), &NoSignatureVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings()
            .iter()
            .any(|f| f.kind == "signature-invalid")
    );
    // Full JSON path: round trip passes; tampered JSON is parse-refused
    // (never a Pass with quietly wrong content).
    let json = package_json(&pkg);
    assert!(check_json(&json, Some(root), Some(&MacVerifier)).passed());
    let tampered = json.replace("bounded", "PROVEN");
    let report = check_json(&tampered, Some(root), Some(&MacVerifier));
    assert!(!report.passed());
    assert!(report.validate_decision_hash());
    assert_eq!(report.findings()[0].kind, "parse-refused");
}
