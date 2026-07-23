//! Battery for the standalone evidence-package checker (addendum Proposal 12).
//! Covers a clean pass, completeness-failure findings, content-address
//! (Merkle) tamper detection, signature-presence reporting, budget-pie
//! rendering (including the empty case), the protocol version, and
//! determinism. The checker uses only the package format — no solver.

use fs_checker::{
    AnchoredSourceRequest, AnchoredSourceVerifier, CHECKER_PROTOCOL_VERSION, ColorBreakdown,
    ContentHash, DerivationRequest, DerivationVerifier, FalsifierRequest, FalsifierVerifier,
    IntegrityStatus, OriginStatus, SemanticFailureKind, SemanticStatus, SignaturePurpose,
    SignatureRequest, SignatureStatus, SignatureVerifier, SourceCertificateRequest,
    SourceCertificateVerifier, Verdict, VerificationCapabilities, VerificationDecision,
    WaiverGrant, WaiverVerifier, check, check_against_root, check_for_release_with_capabilities,
    check_json, check_json_for_release_with_capabilities, check_json_release_preflight,
    check_json_with_capabilities, check_release_preflight, check_with_capabilities,
};
use fs_evidence::{Color, IntervalOp, ValidityDomain};
use fs_package::{Claim, EvidencePackage, FalsifierRecord, Provenance, SemanticWitness};
use std::sync::atomic::{AtomicUsize, Ordering};

const ARTIFACT_HASH: &str = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

/// A deliberately corrupted 32-byte content root (one byte flipped).
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
        let policy = policy_fingerprint("exact-source-verifier");
        if accepted {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
        }
    }
}

impl SourceCertificateVerifier for AlternatePolicySourceVerifier<'_> {
    fn verify(&self, request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let accepted = ExactSourceVerifier {
            claim_id: self.claim_id,
        }
        .verify(request)
        .accepted();
        let policy = ContentHash([0x77; 32]);
        if accepted {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
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
        let policy = policy_fingerprint("exact-anchor-verifier");
        if accepted {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
        }
    }
}

impl FalsifierVerifier for ExactFalsifierVerifier {
    fn verify(&self, request: &FalsifierRequest<'_>) -> VerificationDecision {
        let policy = policy_fingerprint("exact-falsifier-verifier");
        if request.artifact_hash.to_hex() == ARTIFACT_HASH {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
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
    assert_eq!(report.integrity_status(), IntegrityStatus::Verified);
    assert_eq!(report.semantic_status(), SemanticStatus::NotProvided);
    assert_eq!(report.origin_status(), OriginStatus::Authenticated);
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
fn numerical_only_validation_promotion_has_a_specific_checker_finding() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("mesh-bound"))
        .with_claim(
            Claim::derived(
                "forged-validation",
                "a small residual proves agreement with experiment",
                Color::Validated {
                    regime: good_regime(),
                    dataset: "invented-validation-dataset".to_string(),
                },
                vec![0],
                IntervalOp::Hull,
                ARTIFACT_HASH,
            )
            .with_anchor("invented-validation-dataset", ARTIFACT_HASH),
        );
    let report = check(&pkg);
    assert!(!report.passed());
    assert_eq!(report.verdict(), Verdict::Fail);
    assert_eq!(report.findings().len(), 1);
    assert_eq!(
        report.findings()[0].kind,
        "validation-authority-promotion-refused"
    );
    assert!(report.findings()[0].detail.contains("cannot mint"));
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
    assert_eq!(report.integrity_status(), IntegrityStatus::Refused);
    assert_eq!(report.semantic_status(), SemanticStatus::NotRun);
    assert_eq!(report.origin_status(), OriginStatus::NotRun);
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
        let accepted = request.signature == format!("release-test:{}", request.subject_hash())
            && matches!(
                request.purpose,
                SignaturePurpose::ReleaseApproval {
                    checker_protocol,
                    expected_root,
                    admission_context,
                    semantic_context,
                } if expected_root == request.package_root
                    && checker_protocol == CHECKER_PROTOCOL_VERSION
                    && admission_context != ContentHash([0; 32])
                    && semantic_context != ContentHash([0; 32])
            );
        let policy = policy_fingerprint("release-verifier");
        if accepted {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
        }
    }
}

fn signed_for_release(
    package: EvidencePackage,
    capabilities: &VerificationCapabilities<'_>,
) -> EvidencePackage {
    let semantic_context = fs_checker::verify_portable_semantics(&package).context_hash();
    signed_for_release_with_semantic_context(package, capabilities, semantic_context)
}

fn signed_for_release_with_semantic_context(
    package: EvidencePackage,
    capabilities: &VerificationCapabilities<'_>,
    semantic_context: ContentHash,
) -> EvidencePackage {
    let root = package_root(&package);
    let unsigned = package
        .verify_with(capabilities)
        .expect("unsigned release subject verifies");
    let purpose = SignaturePurpose::ReleaseApproval {
        checker_protocol: CHECKER_PROTOCOL_VERSION,
        expected_root: root,
        admission_context: unsigned.receipt().release_admission_context(),
        semantic_context,
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
        let policy = policy_fingerprint("exact-waiver-verifier");
        if mac == fixture_waiver_mac(message) {
            VerificationDecision::accept(policy)
        } else {
            VerificationDecision::reject(policy)
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

#[derive(Default)]
struct UnexpectedCallbacks {
    source: AtomicUsize,
    anchor: AtomicUsize,
    falsifier: AtomicUsize,
    derivation: AtomicUsize,
    waiver: AtomicUsize,
    signature: AtomicUsize,
}

impl UnexpectedCallbacks {
    fn accept(counter: &AtomicUsize) -> VerificationDecision {
        counter.fetch_add(1, Ordering::SeqCst);
        VerificationDecision::accept(ContentHash([0xa5; 32]))
    }

    fn assert_never_called(&self) {
        for (name, counter) in [
            ("source", &self.source),
            ("anchor", &self.anchor),
            ("falsifier", &self.falsifier),
            ("derivation", &self.derivation),
            ("waiver", &self.waiver),
            ("signature", &self.signature),
        ] {
            assert_eq!(
                counter.load(Ordering::SeqCst),
                0,
                "{name} callback ran before expected-root admission"
            );
        }
    }

    fn counts(&self) -> [usize; 6] {
        [
            self.source.load(Ordering::SeqCst),
            self.anchor.load(Ordering::SeqCst),
            self.falsifier.load(Ordering::SeqCst),
            self.derivation.load(Ordering::SeqCst),
            self.waiver.load(Ordering::SeqCst),
            self.signature.load(Ordering::SeqCst),
        ]
    }
}

impl SourceCertificateVerifier for UnexpectedCallbacks {
    fn verify(&self, _request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        Self::accept(&self.source)
    }
}

impl AnchoredSourceVerifier for UnexpectedCallbacks {
    fn verify(&self, _request: &AnchoredSourceRequest<'_>) -> VerificationDecision {
        Self::accept(&self.anchor)
    }
}

impl FalsifierVerifier for UnexpectedCallbacks {
    fn verify(&self, _request: &FalsifierRequest<'_>) -> VerificationDecision {
        Self::accept(&self.falsifier)
    }
}

impl DerivationVerifier for UnexpectedCallbacks {
    fn verify(&self, _request: &DerivationRequest<'_>) -> VerificationDecision {
        Self::accept(&self.derivation)
    }
}

impl WaiverVerifier for UnexpectedCallbacks {
    fn verify(&self, _mac: &str, _message: &[u8]) -> VerificationDecision {
        Self::accept(&self.waiver)
    }
}

impl SignatureVerifier for UnexpectedCallbacks {
    fn verify(&self, _request: &SignatureRequest<'_>) -> VerificationDecision {
        Self::accept(&self.signature)
    }
}

fn all_unexpected_capabilities(callbacks: &UnexpectedCallbacks) -> VerificationCapabilities<'_> {
    VerificationCapabilities::deny_all()
        .with_source_certificates(callbacks)
        .with_anchored_sources(callbacks)
        .with_falsifiers(callbacks)
        .with_derivations(callbacks)
        .with_waivers(callbacks, 100)
        .with_signatures(callbacks)
}

struct PortableSourceVerifier;

impl SourceCertificateVerifier for PortableSourceVerifier {
    fn verify(&self, request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let accepted = request.package_provenance == &prov()
            && request.package_root != ContentHash([0; 32])
            && request.claim_index == 0
            && request.claim_id == "portable"
            && request.statement == "portable result"
            && request.claim_subject_hash != ContentHash([0; 32])
            && request.producer == "portable-test/cert"
            && request.semantic_witness.is_some_and(|witness| {
                let family = witness.family();
                witness.content_hash() == request.certificate_hash
                    && (family == fs_checker::EXACT_INTERVAL_FAMILY
                        || family == fs_checker::BOUNDED_LINF_RESIDUAL_FAMILY)
            });
        if accepted {
            VerificationDecision::accept(ContentHash([0xb4; 32]))
        } else {
            VerificationDecision::reject(ContentHash([0xb4; 32]))
        }
    }
}

fn exact_sum_payload(left: i64, right: i64) -> Vec<u8> {
    let mut payload = Vec::new();
    payload.extend_from_slice(&3_u32.to_le_bytes());
    payload.push(0);
    payload.extend_from_slice(&left.to_le_bytes());
    payload.push(0);
    payload.extend_from_slice(&right.to_le_bytes());
    payload.push(2);
    payload.extend_from_slice(&0_u32.to_le_bytes());
    payload.extend_from_slice(&1_u32.to_le_bytes());
    payload.extend_from_slice(&2_u32.to_le_bytes());
    payload
}

fn exact_division_payload(numerator: i64, denominator: i64) -> Vec<u8> {
    let mut payload = Vec::new();
    payload.extend_from_slice(&3_u32.to_le_bytes());
    payload.push(0);
    payload.extend_from_slice(&numerator.to_le_bytes());
    payload.push(0);
    payload.extend_from_slice(&denominator.to_le_bytes());
    payload.push(5);
    payload.extend_from_slice(&0_u32.to_le_bytes());
    payload.extend_from_slice(&1_u32.to_le_bytes());
    payload.extend_from_slice(&2_u32.to_le_bytes());
    payload
}

fn interval_binary_payload(tag: u8, left: (f64, f64), right: (f64, f64)) -> Vec<u8> {
    let mut payload = Vec::new();
    payload.extend_from_slice(&3_u32.to_le_bytes());
    for (lo, hi) in [left, right] {
        payload.push(1);
        payload.extend_from_slice(&lo.to_bits().to_le_bytes());
        payload.extend_from_slice(&hi.to_bits().to_le_bytes());
    }
    payload.push(tag);
    payload.extend_from_slice(&0_u32.to_le_bytes());
    payload.extend_from_slice(&1_u32.to_le_bytes());
    payload.extend_from_slice(&2_u32.to_le_bytes());
    payload
}

fn residual_1x1_payload(matrix: f64, candidate: f64, right_hand_side: f64) -> Vec<u8> {
    let mut payload = vec![0];
    payload.extend_from_slice(&1_u32.to_le_bytes());
    payload.extend_from_slice(&1_u32.to_le_bytes());
    payload.extend_from_slice(&matrix.to_bits().to_le_bytes());
    payload.extend_from_slice(&candidate.to_bits().to_le_bytes());
    payload.extend_from_slice(&right_hand_side.to_bits().to_le_bytes());
    payload
}

fn residual_2x1_payload() -> Vec<u8> {
    let mut payload = vec![0];
    payload.extend_from_slice(&2_u32.to_le_bytes());
    payload.extend_from_slice(&1_u32.to_le_bytes());
    for value in [1.0_f64, 2.0, 1.0, 2.0, 4.0] {
        payload.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    payload
}

fn exact_one_work_payload(multiplications: usize, hulls: usize) -> Vec<u8> {
    let nodes = 1 + multiplications + hulls;
    let mut payload = Vec::new();
    payload.extend_from_slice(&(nodes as u32).to_le_bytes());
    payload.push(0);
    payload.extend_from_slice(&1_i64.to_le_bytes());
    for current in 1..=multiplications {
        payload.push(4);
        payload.extend_from_slice(&((current - 1) as u32).to_le_bytes());
        payload.extend_from_slice(&0_u32.to_le_bytes());
    }
    for current in (multiplications + 1)..nodes {
        payload.push(7);
        payload.extend_from_slice(&((current - 1) as u32).to_le_bytes());
        payload.extend_from_slice(&0_u32.to_le_bytes());
    }
    payload.extend_from_slice(&((nodes - 1) as u32).to_le_bytes());
    payload
}

fn portable_package(witness: SemanticWitness, lo: f64, hi: f64) -> EvidencePackage {
    EvidencePackage::new(prov()).with_claim(Claim::from_portable_certificate(
        "portable",
        "portable result",
        lo,
        hi,
        "portable-test/cert",
        witness,
    ))
}

#[test]
fn built_in_semantic_families_verify_in_memory_and_json() {
    let fixtures = [
        (
            SemanticWitness::new(
                fs_checker::EXACT_INTERVAL_FAMILY,
                fs_checker::INITIAL_SEMANTIC_SCHEMA_VERSION,
                exact_sum_payload(1, 2),
            ),
            3.0,
            3.0,
        ),
        // Interval leaves deliberately bypass the exact-i53 shortcut. Raw
        // [1, 2] + [3, 4] is [4, 6], expanded by exactly one adjacent value.
        (
            SemanticWitness::new(
                fs_checker::EXACT_INTERVAL_FAMILY,
                fs_checker::INITIAL_SEMANTIC_SCHEMA_VERSION,
                interval_binary_payload(2, (1.0, 2.0), (3.0, 4.0)),
            ),
            f64::from_bits(0x400f_ffff_ffff_ffff),
            f64::from_bits(0x4018_0000_0000_0001),
        ),
        // Raw positive products span [2*4, 3*5] = [8, 15].
        (
            SemanticWitness::new(
                fs_checker::EXACT_INTERVAL_FAMILY,
                fs_checker::INITIAL_SEMANTIC_SCHEMA_VERSION,
                interval_binary_payload(4, (2.0, 3.0), (4.0, 5.0)),
            ),
            f64::from_bits(0x401f_ffff_ffff_ffff),
            f64::from_bits(0x402e_0000_0000_0001),
        ),
        // Raw positive quotients span [6/4, 8/2] = [1.5, 4].
        (
            SemanticWitness::new(
                fs_checker::EXACT_INTERVAL_FAMILY,
                fs_checker::INITIAL_SEMANTIC_SCHEMA_VERSION,
                interval_binary_payload(5, (6.0, 8.0), (2.0, 4.0)),
            ),
            f64::from_bits(0x3ff7_ffff_ffff_ffff),
            f64::from_bits(0x4010_0000_0000_0001),
        ),
        // A=[[1],[2]], x=[1], b=[2,4] has residual [1,2]. The
        // independently hand-derived local enclosure ends two ulps above 2.
        (
            SemanticWitness::new(
                fs_checker::BOUNDED_LINF_RESIDUAL_FAMILY,
                fs_checker::INITIAL_SEMANTIC_SCHEMA_VERSION,
                residual_2x1_payload(),
            ),
            0.0,
            f64::from_bits(0x4000_0000_0000_0002),
        ),
    ];
    let source_verifier = PortableSourceVerifier;
    let capabilities =
        VerificationCapabilities::deny_all().with_source_certificates(&source_verifier);

    for (witness, lo, hi) in fixtures {
        let package = portable_package(witness, lo, hi);
        let missing_origin = check(&package);
        assert_eq!(missing_origin.integrity_status(), IntegrityStatus::Verified);
        assert_eq!(missing_origin.semantic_status(), SemanticStatus::Verified);
        assert_eq!(missing_origin.origin_status(), OriginStatus::Refused);
        assert_eq!(missing_origin.breakdown(), &ColorBreakdown::default());

        let in_memory = check_with_capabilities(&package, None, None, &capabilities);
        assert!(in_memory.passed(), "{:?}", in_memory.findings());
        assert_eq!(in_memory.integrity_status(), IntegrityStatus::Verified);
        assert_eq!(in_memory.semantic_status(), SemanticStatus::Verified);
        assert_eq!(in_memory.origin_status(), OriginStatus::Authenticated);
        assert!(in_memory.semantic_report().validate_context_hash());
        assert_eq!(
            in_memory.semantic_report().claims()[0].status(),
            fs_checker::SemanticClaimStatus::Verified
        );

        let from_json =
            check_json_with_capabilities(&package_json(&package), None, None, &capabilities);
        assert_eq!(from_json, in_memory);
    }
}

#[test]
fn semantic_refusals_precede_every_external_callback() {
    let mut trailing = exact_sum_payload(1, 2);
    trailing.push(0xff);
    let mut nan_leaf = Vec::new();
    nan_leaf.extend_from_slice(&1_u32.to_le_bytes());
    nan_leaf.push(1);
    nan_leaf.extend_from_slice(&f64::NAN.to_bits().to_le_bytes());
    nan_leaf.extend_from_slice(&1.0_f64.to_bits().to_le_bytes());
    nan_leaf.extend_from_slice(&0_u32.to_le_bytes());
    let oversized_node_count = (fs_checker::MAX_INTERVAL_NODES as u32 + 1)
        .to_le_bytes()
        .to_vec();
    let mut forward_reference = Vec::new();
    forward_reference.extend_from_slice(&3_u32.to_le_bytes());
    forward_reference.push(0);
    forward_reference.extend_from_slice(&1_i64.to_le_bytes());
    forward_reference.push(2);
    forward_reference.extend_from_slice(&2_u32.to_le_bytes());
    forward_reference.extend_from_slice(&0_u32.to_le_bytes());
    let mut self_reference = Vec::new();
    self_reference.extend_from_slice(&2_u32.to_le_bytes());
    self_reference.push(0);
    self_reference.extend_from_slice(&1_i64.to_le_bytes());
    self_reference.push(2);
    self_reference.extend_from_slice(&1_u32.to_le_bytes());
    self_reference.extend_from_slice(&0_u32.to_le_bytes());
    let mut truncated = Vec::new();
    truncated.extend_from_slice(&1_u32.to_le_bytes());
    truncated.push(0);
    truncated.extend_from_slice(&[1, 2, 3]);
    let mut oversized_dimension = vec![0];
    oversized_dimension
        .extend_from_slice(&(fs_checker::MAX_RESIDUAL_DIMENSION as u32 + 1).to_le_bytes());
    oversized_dimension.extend_from_slice(&1_u32.to_le_bytes());
    let mut residual_extra_element = residual_1x1_payload(0.0, 0.0, 0.0);
    residual_extra_element.extend_from_slice(&0.0_f64.to_bits().to_le_bytes());

    let fixtures = [
        (
            SemanticWitness::new(
                fs_checker::EXACT_INTERVAL_FAMILY,
                1,
                exact_sum_payload(1, 2),
            ),
            4.0,
            4.0,
            SemanticFailureKind::ClaimMismatch,
            "semantic-claim-mismatch",
        ),
        (
            SemanticWitness::new(
                fs_checker::BOUNDED_LINF_RESIDUAL_FAMILY,
                1,
                residual_1x1_payload(2.0, 3.0, 5.0),
            ),
            0.0,
            0.0,
            SemanticFailureKind::ClaimMismatch,
            "semantic-claim-mismatch",
        ),
        (
            SemanticWitness::new("frankensim/unknown-proof", 1, exact_sum_payload(1, 2)),
            3.0,
            3.0,
            SemanticFailureKind::UnknownFamily,
            "semantic-unknown-family",
        ),
        (
            SemanticWitness::new(
                fs_checker::EXACT_INTERVAL_FAMILY,
                2,
                exact_sum_payload(1, 2),
            ),
            3.0,
            3.0,
            SemanticFailureKind::UnsupportedVersion,
            "semantic-unsupported-version",
        ),
        (
            SemanticWitness::new(fs_checker::EXACT_INTERVAL_FAMILY, 1, trailing),
            3.0,
            3.0,
            SemanticFailureKind::MalformedPayload,
            "semantic-malformed-witness",
        ),
        (
            SemanticWitness::new(
                fs_checker::BOUNDED_LINF_RESIDUAL_FAMILY,
                1,
                residual_1x1_payload(f64::NAN, 0.0, 0.0),
            ),
            0.0,
            0.0,
            SemanticFailureKind::MalformedPayload,
            "semantic-malformed-witness",
        ),
        (
            SemanticWitness::new(fs_checker::EXACT_INTERVAL_FAMILY, 1, nan_leaf),
            0.0,
            1.0,
            SemanticFailureKind::MalformedPayload,
            "semantic-malformed-witness",
        ),
        (
            SemanticWitness::new(
                fs_checker::EXACT_INTERVAL_FAMILY,
                1,
                exact_division_payload(1, 0),
            ),
            0.0,
            0.0,
            SemanticFailureKind::MalformedPayload,
            "semantic-malformed-witness",
        ),
        (
            SemanticWitness::new(fs_checker::EXACT_INTERVAL_FAMILY, 1, oversized_node_count),
            0.0,
            0.0,
            SemanticFailureKind::ResourceLimit,
            "semantic-resource-limit",
        ),
        (
            SemanticWitness::new(fs_checker::EXACT_INTERVAL_FAMILY, 1, forward_reference),
            1.0,
            1.0,
            SemanticFailureKind::MalformedPayload,
            "semantic-malformed-witness",
        ),
        (
            SemanticWitness::new(fs_checker::EXACT_INTERVAL_FAMILY, 1, self_reference),
            1.0,
            1.0,
            SemanticFailureKind::MalformedPayload,
            "semantic-malformed-witness",
        ),
        (
            SemanticWitness::new(fs_checker::EXACT_INTERVAL_FAMILY, 1, truncated),
            1.0,
            1.0,
            SemanticFailureKind::MalformedPayload,
            "semantic-malformed-witness",
        ),
        (
            SemanticWitness::new(
                fs_checker::BOUNDED_LINF_RESIDUAL_FAMILY,
                1,
                oversized_dimension,
            ),
            0.0,
            0.0,
            SemanticFailureKind::ResourceLimit,
            "semantic-resource-limit",
        ),
        (
            SemanticWitness::new(
                fs_checker::BOUNDED_LINF_RESIDUAL_FAMILY,
                1,
                residual_extra_element,
            ),
            0.0,
            f64::from_bits(3),
            SemanticFailureKind::MalformedPayload,
            "semantic-malformed-witness",
        ),
    ];

    for (witness, lo, hi, failure_kind, finding_kind) in fixtures {
        let package = portable_package(witness, lo, hi);
        let callbacks = UnexpectedCallbacks::default();
        let capabilities = all_unexpected_capabilities(&callbacks);
        let report =
            check_with_capabilities(&package, Some(package_root(&package)), None, &capabilities);
        callbacks.assert_never_called();
        assert_eq!(report.integrity_status(), IntegrityStatus::Verified);
        assert_eq!(report.semantic_status(), SemanticStatus::Refused);
        assert_eq!(report.origin_status(), OriginStatus::NotRun);
        assert_eq!(report.breakdown(), &ColorBreakdown::default());
        assert!(report.receipt().is_none());
        assert!(report.validate_decision_hash());
        assert_eq!(report.semantic_report().failures()[0].kind(), failure_kind);
        assert!(
            report
                .findings()
                .iter()
                .any(|finding| finding.kind == finding_kind),
            "missing {finding_kind}: {:?}",
            report.findings()
        );
    }
}

#[test]
fn semantic_payload_tamper_with_a_recomputed_root_still_refuses() {
    let original = portable_package(
        SemanticWitness::new(
            fs_checker::EXACT_INTERVAL_FAMILY,
            1,
            exact_sum_payload(1, 2),
        ),
        3.0,
        3.0,
    );
    let tampered = portable_package(
        SemanticWitness::new(
            fs_checker::EXACT_INTERVAL_FAMILY,
            1,
            exact_sum_payload(1, 3),
        ),
        3.0,
        3.0,
    );
    assert_ne!(package_root(&original), package_root(&tampered));

    let callbacks = UnexpectedCallbacks::default();
    let capabilities = all_unexpected_capabilities(&callbacks);
    let expected = package_root(&tampered);
    let in_memory = check_with_capabilities(&tampered, Some(expected), None, &capabilities);
    callbacks.assert_never_called();
    assert_eq!(in_memory.merkle_root(), expected);
    assert_eq!(in_memory.semantic_status(), SemanticStatus::Refused);
    assert!(
        in_memory
            .findings()
            .iter()
            .any(|finding| finding.kind == "semantic-claim-mismatch")
    );

    let from_json = check_json_with_capabilities(
        &package_json(&tampered),
        Some(expected),
        None,
        &capabilities,
    );
    callbacks.assert_never_called();
    assert_eq!(from_json, in_memory);
}

#[test]
fn aggregate_semantic_operation_budget_refuses_at_limit_plus_one() {
    // Thirty maximum-size multiplication programs charge 982,830 operations.
    // The final program charges exactly 17,171 more, so its last hull would
    // cross the package budget from 1,000,000 to 1,000,001.
    let maximum_work = exact_one_work_payload(fs_checker::MAX_INTERVAL_NODES - 1, 0);
    let limit_plus_one_work = exact_one_work_payload(2_146, 2);
    let mut package = EvidencePackage::new(prov());
    for index in 0..30 {
        package = package.with_claim(Claim::from_portable_certificate(
            format!("operation-budget-{index}"),
            "bounded exact multiplication work",
            1.0,
            1.0,
            "portable-test/cert",
            SemanticWitness::new(fs_checker::EXACT_INTERVAL_FAMILY, 1, maximum_work.clone()),
        ));
    }
    package = package.with_claim(Claim::from_portable_certificate(
        "operation-budget-limit-plus-one",
        "bounded exact multiplication work",
        1.0,
        1.0,
        "portable-test/cert",
        SemanticWitness::new(fs_checker::EXACT_INTERVAL_FAMILY, 1, limit_plus_one_work),
    ));

    let callbacks = UnexpectedCallbacks::default();
    let report = check_with_capabilities(
        &package,
        Some(package_root(&package)),
        None,
        &all_unexpected_capabilities(&callbacks),
    );
    callbacks.assert_never_called();
    assert_eq!(report.semantic_status(), SemanticStatus::Refused);
    assert_eq!(
        report.semantic_report().operations(),
        fs_checker::MAX_SEMANTIC_OPERATIONS
    );
    assert_eq!(
        report
            .semantic_report()
            .failures()
            .last()
            .map(|failure| failure.kind()),
        Some(SemanticFailureKind::ResourceLimit)
    );
    assert!(
        report
            .semantic_report()
            .failures()
            .last()
            .is_some_and(|failure| failure.detail().contains("operation limit"))
    );
}

#[test]
fn release_signature_cannot_replay_a_different_semantic_context() {
    let witness = SemanticWitness::new(
        fs_checker::EXACT_INTERVAL_FAMILY,
        1,
        exact_sum_payload(1, 2),
    );
    let package = EvidencePackage::new(prov())
        .with_claim(
            Claim::from_portable_certificate(
                "portable",
                "portable result",
                3.0,
                3.0,
                "portable-test/cert",
                witness,
            )
            .with_falsifier(passed_falsifier()),
        )
        .with_claim(estimated("context-companion"));
    let source_verifier = PortableSourceVerifier;
    let capabilities = VerificationCapabilities::deny_all()
        .with_source_certificates(&source_verifier)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER);

    let replayed = signed_for_release_with_semantic_context(
        package.clone(),
        &capabilities,
        ContentHash([0x42; 32]),
    );
    let refused = check_for_release_with_capabilities(
        &replayed,
        package_root(&replayed),
        &ReleaseVerifier,
        &capabilities,
    );
    assert_eq!(refused.semantic_status(), SemanticStatus::Verified);
    assert_eq!(refused.origin_status(), OriginStatus::Authenticated);
    assert_capability_refusal(&refused, "signature-invalid");

    let correctly_signed = signed_for_release(package, &capabilities);
    let admitted = check_for_release_with_capabilities(
        &correctly_signed,
        package_root(&correctly_signed),
        &ReleaseVerifier,
        &capabilities,
    );
    assert!(admitted.release_admitted(), "{:?}", admitted.findings());
    assert!(admitted.release_independently_verified());
}

fn mixed_origin_signed_fixture() -> EvidencePackage {
    let estimate = estimated("estimate");
    let derived_color = fs_evidence::compose(
        estimate.declared_color_unverified(),
        estimate.declared_color_unverified(),
        IntervalOp::Hull,
    );
    EvidencePackage::new(prov())
        .with_claim(verified("certificate").with_falsifier(passed_falsifier()))
        .with_claim(validated("anchor", good_regime()).with_falsifier(passed_falsifier()))
        .with_claim(estimate)
        .with_claim(
            Claim::derived(
                "derived",
                "derived estimate",
                derived_color,
                vec![2, 2],
                IntervalOp::Hull,
                ARTIFACT_HASH,
            )
            .with_falsifier(passed_falsifier()),
        )
        .with_claim(
            Claim::waived(
                "waived",
                "approved exception",
                Color::Verified { lo: -1.0, hi: 1.0 },
                WaiverGrant {
                    waiver_id: "wrong-root-waiver".to_string(),
                    expiry_day: 200,
                    mac: "wrong-root-mac".to_string(),
                },
            )
            .with_falsifier(passed_falsifier()),
        )
        .signed("wrong-root-signature")
}

fn assert_root_mismatch_refusal(
    report: &fs_checker::CheckReport,
    actual: ContentHash,
    expected: ContentHash,
) {
    assert_eq!(report.verdict(), Verdict::Fail);
    assert_eq!(report.merkle_root(), actual);
    assert_eq!(report.expected_root(), Some(expected));
    assert_eq!(report.breakdown(), &ColorBreakdown::default());
    assert!(report.receipt().is_none());
    assert!(report.validate_decision_hash());
    assert!(report.findings().iter().any(|finding| {
        finding.kind == "content-address-mismatch"
            && finding.detail.contains(&actual.to_string())
            && finding.detail.contains(&expected.to_string())
    }));
}

fn assert_release_shape_refuses_without_callbacks(package: &EvidencePackage, kind: &str) {
    let root = package_root(package);
    let callbacks = UnexpectedCallbacks::default();
    let capabilities = all_unexpected_capabilities(&callbacks);

    let in_memory = check_for_release_with_capabilities(package, root, &callbacks, &capabilities);
    callbacks.assert_never_called();
    assert_eq!(in_memory.verdict(), Verdict::Fail);
    assert_eq!(in_memory.merkle_root(), root);
    assert_eq!(in_memory.expected_root(), Some(root));
    assert_eq!(in_memory.breakdown(), &ColorBreakdown::default());
    assert!(in_memory.receipt().is_none());
    assert!(in_memory.validate_decision_hash());
    assert!(!matches!(
        in_memory.signature(),
        SignatureStatus::Authenticated(_)
    ));
    assert!(
        in_memory
            .findings()
            .iter()
            .any(|finding| finding.kind == kind),
        "missing {kind} finding: {:?}",
        in_memory.findings()
    );

    let json = check_json_for_release_with_capabilities(
        &package_json(package),
        root,
        &callbacks,
        &capabilities,
    );
    callbacks.assert_never_called();
    assert_eq!(json, in_memory);
}

#[test]
fn expected_root_mismatch_precedes_all_integrity_callbacks_in_memory_and_json() {
    let package = mixed_origin_signed_fixture();
    let actual = package_root(&package);
    let expected = flip(actual);
    let callbacks = UnexpectedCallbacks::default();
    let capabilities = all_unexpected_capabilities(&callbacks);

    let in_memory =
        check_with_capabilities(&package, Some(expected), Some(&callbacks), &capabilities);
    callbacks.assert_never_called();
    assert_root_mismatch_refusal(&in_memory, actual, expected);
    assert!(matches!(
        in_memory.signature(),
        SignatureStatus::Unverified(_)
    ));

    let json = check_json_with_capabilities(
        &package_json(&package),
        Some(expected),
        Some(&callbacks),
        &capabilities,
    );
    callbacks.assert_never_called();
    assert_eq!(json, in_memory);
}

#[test]
fn expected_root_mismatch_precedes_all_release_callbacks_in_memory_and_json() {
    let package = mixed_origin_signed_fixture();
    let actual = package_root(&package);
    let expected = flip(actual);
    let callbacks = UnexpectedCallbacks::default();
    let capabilities = all_unexpected_capabilities(&callbacks);

    let in_memory =
        check_for_release_with_capabilities(&package, expected, &callbacks, &capabilities);
    callbacks.assert_never_called();
    assert_root_mismatch_refusal(&in_memory, actual, expected);
    assert!(!in_memory.release_admitted());

    let json = check_json_for_release_with_capabilities(
        &package_json(&package),
        expected,
        &callbacks,
        &capabilities,
    );
    callbacks.assert_never_called();
    assert_eq!(json, in_memory);
}

#[test]
fn release_empty_package_refuses_before_all_callbacks_in_memory_and_json() {
    let package = EvidencePackage::new(prov()).signed("shape-test-signature");
    assert_release_shape_refuses_without_callbacks(&package, "release-empty-package");
}

#[test]
fn release_unsigned_package_refuses_before_all_callbacks_in_memory_and_json() {
    let package = EvidencePackage::new(prov())
        .with_claim(verified("unsigned").with_falsifier(passed_falsifier()));
    assert_release_shape_refuses_without_callbacks(&package, "release-signature-required");
}

#[test]
fn release_all_waived_package_refuses_before_all_callbacks_in_memory_and_json() {
    let package = EvidencePackage::new(prov())
        .with_claim(
            Claim::waived(
                "waived-only",
                "authorized exception",
                Color::Verified { lo: -1.0, hi: 1.0 },
                WaiverGrant {
                    waiver_id: "shape-waiver".to_string(),
                    expiry_day: 200,
                    mac: "shape-waiver-mac".to_string(),
                },
            )
            .with_falsifier(passed_falsifier()),
        )
        .signed("shape-test-signature");
    assert_release_shape_refuses_without_callbacks(
        &package,
        "release-scientific-evidence-required",
    );
}

#[test]
fn release_missing_falsifier_refuses_before_all_callbacks_in_memory_and_json() {
    let package = EvidencePackage::new(prov())
        .with_claim(verified("missing-falsifier"))
        .signed("shape-test-signature");
    assert_release_shape_refuses_without_callbacks(&package, "release-falsifier-required");
}

#[test]
fn release_missing_anchor_refuses_before_all_callbacks_in_memory_and_json() {
    let parent = validated("anchored-parent", good_regime()).with_falsifier(passed_falsifier());
    let derived_color = parent.declared_color_unverified().clone();
    let package = EvidencePackage::new(prov())
        .with_claim(parent)
        .with_claim(
            Claim::derived(
                "missing-anchor",
                "derived validation",
                derived_color,
                vec![0],
                IntervalOp::Hull,
                ARTIFACT_HASH,
            )
            .with_falsifier(passed_falsifier()),
        )
        .signed("shape-test-signature");
    assert_release_shape_refuses_without_callbacks(&package, "release-anchor-required");
}

#[test]
fn release_complete_shape_dispatches_every_callback_in_memory_and_json() {
    let package = mixed_origin_signed_fixture();
    let root = package_root(&package);
    let callbacks = UnexpectedCallbacks::default();
    let capabilities = all_unexpected_capabilities(&callbacks);

    let in_memory = check_for_release_with_capabilities(&package, root, &callbacks, &capabilities);
    assert!(in_memory.release_admitted(), "{:?}", in_memory.findings());
    assert_eq!(callbacks.counts(), [1, 1, 4, 1, 1, 1]);

    let json = check_json_for_release_with_capabilities(
        &package_json(&package),
        root,
        &callbacks,
        &capabilities,
    );
    assert_eq!(json, in_memory);
    assert_eq!(callbacks.counts(), [2, 2, 8, 2, 2, 2]);
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
    let release =
        check_for_release_with_capabilities(&signed, root, &ReleaseVerifier, &capabilities);
    let release_json = check_json_for_release_with_capabilities(
        &package_json(&signed),
        root,
        &ReleaseVerifier,
        &capabilities,
    );
    assert_eq!(release_json, release);
    for report in [&release, &release_json] {
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
        assert_eq!(report.breakdown(), &ColorBreakdown::default());
        assert!(report.receipt().is_none());
        assert!(matches!(report.signature(), SignatureStatus::Unverified(_)));
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
    let json = check_json_for_release_with_capabilities(
        &package_json(&package),
        root,
        &ReleaseVerifier,
        &VerificationCapabilities::deny_all(),
    );
    assert_eq!(json, report);
    assert!(!report.passed());
    assert!(report.validate_decision_hash());
    assert_eq!(report.breakdown(), &ColorBreakdown::default());
    assert!(report.receipt().is_none());
    assert!(matches!(report.signature(), SignatureStatus::Unverified(_)));
    assert!(
        report
            .findings()
            .iter()
            .any(|finding| finding.kind == "release-scientific-evidence-required")
    );
}

#[test]
fn release_gate_rejects_sole_vacuous_verified_enclosure_under_valid_policies() {
    struct VacuousSourceVerifier;

    impl SourceCertificateVerifier for VacuousSourceVerifier {
        fn verify(&self, request: &SourceCertificateRequest<'_>) -> VerificationDecision {
            let accepted = request.package_provenance == &prov()
                && request.claim_index == 0
                && request.claim_id == "vacuous"
                && request.statement == "sound but vacuous enclosure"
                && request.lo == f64::NEG_INFINITY
                && request.hi == f64::INFINITY
                && request.producer == "test-solver/cert"
                && request.certificate_hash.to_hex() == ARTIFACT_HASH;
            let policy = ContentHash([0x7a; 32]);
            if accepted {
                VerificationDecision::accept(policy)
            } else {
                VerificationDecision::reject(policy)
            }
        }
    }

    let source_verifier = VacuousSourceVerifier;
    let capabilities = VerificationCapabilities::deny_all()
        .with_source_certificates(&source_verifier)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER);
    let unsigned = EvidencePackage::new(prov()).with_claim(
        Claim::from_certificate(
            "vacuous",
            "sound but vacuous enclosure",
            f64::NEG_INFINITY,
            f64::INFINITY,
            "test-solver/cert",
            ARTIFACT_HASH,
        )
        .with_falsifier(passed_falsifier()),
    );

    let integrity = check_with_capabilities(&unsigned, None, None, &capabilities);
    assert!(
        integrity.passed(),
        "vacuous enclosure remains an honest integrity artifact: {:?}",
        integrity.findings()
    );
    assert_eq!(integrity.breakdown().verified, 1);

    let package = signed_for_release(unsigned, &capabilities);
    let root = package_root(&package);
    let report =
        check_for_release_with_capabilities(&package, root, &ReleaseVerifier, &capabilities);
    let json = check_json_for_release_with_capabilities(
        &package_json(&package),
        root,
        &ReleaseVerifier,
        &capabilities,
    );
    assert_eq!(json, report);
    for report in [&report, &json] {
        assert!(!report.passed());
        assert!(!report.release_admitted());
        assert!(report.validate_decision_hash());
        assert_eq!(report.breakdown(), &ColorBreakdown::default());
        assert!(report.receipt().is_none());
        assert!(matches!(report.signature(), SignatureStatus::Unverified(_)));
        assert!(
            report
                .findings()
                .iter()
                .any(|finding| finding.kind == "release-scientific-evidence-required"),
            "vacuous evidence must be named as the release blocker: {:?}",
            report.findings()
        );
    }
}

#[test]
fn release_gate_refuses_package_root_attestation_substitution() {
    struct RootAttestationVerifier;

    impl fs_checker::SignatureVerifier for RootAttestationVerifier {
        fn verify(&self, request: &SignatureRequest<'_>) -> VerificationDecision {
            let accepted = request.signature
                == format!("integrity-test:{}", request.subject_hash())
                && request.purpose == SignaturePurpose::PackageRootAttestation;
            if accepted {
                VerificationDecision::accept(ContentHash([0x88; 32]))
            } else {
                VerificationDecision::reject(ContentHash([0x88; 32]))
            }
        }
    }

    let unsigned = EvidencePackage::new(prov())
        .with_claim(verified("integrity-only").with_falsifier(passed_falsifier()));
    let root = package_root(&unsigned);
    let attestation_subject =
        fs_checker::signature_subject_hash(root, SignaturePurpose::PackageRootAttestation);
    let package = unsigned.signed(format!("integrity-test:{attestation_subject}"));
    let report = check_for_release_with_capabilities(
        &package,
        root,
        &RootAttestationVerifier,
        &source_capabilities(&ExactSourceVerifier {
            claim_id: "integrity-only",
        }),
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

    let pending_waiver = EvidencePackage::new(prov())
        .with_claim(verified("source").with_falsifier(passed_falsifier()))
        .with_claim(
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
    let waiver_message = pending_waiver
        .waiver_message(1)
        .expect("combined waiver target");
    let waiver_package = pending_waiver
        .with_waiver_mac(1, fixture_waiver_mac(&waiver_message))
        .expect("install combined fixture authenticator");
    let waiver_verifier = ExactWaiverVerifier;
    let waiver_source = ExactSourceVerifier { claim_id: "source" };
    let day_250 = source_capabilities(&waiver_source)
        .with_waivers(&waiver_verifier, 250)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER);
    let day_249 = source_capabilities(&waiver_source)
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
fn release_preflight_inventories_declaration_blockers_after_capability_refusal() {
    let validated_color = Color::Validated {
        regime: good_regime(),
        dataset: "wt-2026".to_string(),
    };
    let package = EvidencePackage::new(prov())
        .with_claim(validated("parent", good_regime()))
        .with_claim(Claim::derived(
            "derived",
            "matches",
            validated_color,
            vec![0],
            IntervalOp::Hull,
            ARTIFACT_HASH,
        ));
    assert!(package.is_structurally_inspectable_unverified());

    let report = check_release_preflight(&package, package_root(&package), &ReleaseVerifier);
    assert!(
        report.receipt().is_none(),
        "deny-all cannot admit the anchor"
    );
    for expected in [
        "anchored-source-refused",
        "release-scientific-evidence-required",
        "release-signature-required",
        "release-falsifier-required",
        "release-anchor-required",
        "release-preflight-only",
    ] {
        assert!(
            report
                .findings()
                .iter()
                .any(|finding| finding.kind == expected),
            "missing {expected}: {:?}",
            report.findings()
        );
    }
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
fn release_gate_rejects_derived_validated_anchor_substitution() {
    struct DerivedAnchorVerifier;
    struct ExactDerivationVerifier;

    impl AnchoredSourceVerifier for DerivedAnchorVerifier {
        fn verify(&self, request: &AnchoredSourceRequest<'_>) -> VerificationDecision {
            let accepted = request.package_provenance == &prov()
                && request.statement == "matches"
                && request.dataset_id == "wt-2026"
                && request.content_hash.to_hex() == ARTIFACT_HASH
                && request.regime == &good_regime()
                && matches!(
                    (request.claim_index, request.claim_id),
                    (0, "parent") | (1, "derived")
                );
            let policy = policy_fingerprint("exact-anchor-verifier");
            if accepted {
                VerificationDecision::accept(policy)
            } else {
                VerificationDecision::reject(policy)
            }
        }
    }

    impl DerivationVerifier for ExactDerivationVerifier {
        fn verify(&self, request: &DerivationRequest<'_>) -> VerificationDecision {
            if request.artifact_hash.to_hex() == ARTIFACT_HASH {
                VerificationDecision::accept(ContentHash([0x99; 32]))
            } else {
                VerificationDecision::reject(ContentHash([0x99; 32]))
            }
        }
    }

    let substituted_hash = "1123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";
    let package = EvidencePackage::new(prov())
        .with_claim(validated("parent", good_regime()).with_falsifier(passed_falsifier()))
        .with_claim(
            Claim::derived(
                "derived",
                "matches",
                Color::Validated {
                    regime: good_regime(),
                    dataset: "wt-2026".to_string(),
                },
                vec![0],
                IntervalOp::Hull,
                ARTIFACT_HASH,
            )
            .with_anchor("wt-2026", ARTIFACT_HASH)
            .with_anchor("wt-2026", substituted_hash)
            .with_falsifier(passed_falsifier()),
        )
        .signed("release-test:forged");
    let capabilities = VerificationCapabilities::deny_all()
        .with_anchored_sources(&DerivedAnchorVerifier)
        .with_derivations(&ExactDerivationVerifier)
        .with_falsifiers(&EXACT_FALSIFIER_VERIFIER);
    let report = check_for_release_with_capabilities(
        &package,
        package_root(&package),
        &ReleaseVerifier,
        &capabilities,
    );
    assert_capability_refusal(&report, "anchored-source-refused");
    assert!(
        report.findings().iter().any(|finding| {
            finding.kind == "anchored-source-refused" && finding.detail.contains("claim 'derived'")
        }),
        "release must reject the substituted derived anchor specifically: {:?}",
        report.findings()
    );
}

#[test]
fn the_checker_advertises_its_protocol_version() {
    assert_eq!(CHECKER_PROTOCOL_VERSION, 7);
    assert_eq!(fs_checker::CHECKER_SUPPORTED_PACKAGE_FORMAT, 9);
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
            let accepted = request.signature == mac(request.subject_hash())
                && request.purpose == SignaturePurpose::PackageRootAttestation;
            let policy = policy_fingerprint("mac-verifier");
            if accepted {
                VerificationDecision::accept(policy)
            } else {
                VerificationDecision::reject(policy)
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
