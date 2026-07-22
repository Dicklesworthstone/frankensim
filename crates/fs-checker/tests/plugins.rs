//! Conformance battery for the public semantic-plugin facade (bead 9e8n).
//!
//! Fixtures are hand-authored canonical bytes with independently pinned
//! witness, plugin, and registry identities. No checker-side producer or
//! arithmetic helper constructs the expected evidence. Every positive case
//! travels through `SemanticWitness -> Claim -> EvidencePackage -> strict JSON`
//! and the same package-bound semantic report consumed by release admission.

use fs_checker::plugins::{SemanticPluginDescriptor, SemanticReport};
use fs_checker::{
    CHECKER_PROTOCOL_VERSION, ContentHash, FalsifierRequest, FalsifierVerifier, IntegrityStatus,
    OriginStatus, SemanticClaimStatus, SemanticFailureKind, SemanticStatus, SignaturePurpose,
    SignatureRequest, SignatureVerifier, SourceCertificateRequest, SourceCertificateVerifier,
    VerificationCapabilities, VerificationDecision, check_for_release_with_capabilities,
    check_json_for_release_with_capabilities, check_json_with_capabilities,
    check_with_capabilities, signature_subject_hash,
};
use fs_package::{Claim, EvidencePackage, FalsifierRecord, Provenance, SemanticWitness};

const ARTIFACT_HASH: &str = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";
const EXACT_WITNESS_HASH: &str = "a0e653f9901ee63e4f03c147715362cff0ae6f002e786e992f3a008138cc3345";
const RESIDUAL_WITNESS_HASH: &str =
    "3744095132bc420c754081385509648b71c88c5de7a0d1b5dddf22eea37e3730";
const EXACT_PLUGIN_HASH: &str = "14af87a06c1aa93da84bcc5cb020301e9d4a0aa1b75d14e1fbcb0056c1b66276";
const RESIDUAL_PLUGIN_HASH: &str =
    "2b18acdd01d48caecadc6602d90f7d35f98b2e68275d6b07a56782659feebcda";
const REGISTRY_HASH: &str = "c9bef0915f642dcded1ce599a4705d7884c84a45ba4adc4c39b8740b676ee402";

// Canonical exact-interval v1 program: integer leaves 1 and 2, exact addition,
// result node 2. All integer/index fields are little-endian by schema.
const EXACT_SUM_PAYLOAD: &[u8] = &[
    0x03, 0x00, 0x00, 0x00, // node count
    0x00, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // exact 1
    0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, // exact 2
    0x02, 0x00, 0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, // add(0, 1)
    0x02, 0x00, 0x00, 0x00, // result index
];

// Canonical residual v1 program: norm tag 0, A=[[1],[2]], x=[1], b=[2,4].
// Its independently hand-derived L-infinity enclosure is [0, next_up(next_up(2))].
const RESIDUAL_2X1_PAYLOAD: &[u8] = &[
    0x00, // norm tag
    0x02, 0x00, 0x00, 0x00, // rows
    0x01, 0x00, 0x00, 0x00, // columns
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f, // A[0] = 1
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, // A[1] = 2
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xf0, 0x3f, // x[0] = 1
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, // b[0] = 2
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10, 0x40, // b[1] = 4
];

#[derive(Clone, Copy)]
struct Fixture {
    claim_id: &'static str,
    family: &'static str,
    payload: &'static [u8],
    lo: f64,
    hi: f64,
    witness_hash: &'static str,
    plugin_hash: &'static str,
}

const FIXTURES: [Fixture; 2] = [
    Fixture {
        claim_id: "exact-sum",
        family: fs_checker::plugins::EXACT_INTERVAL_FAMILY,
        payload: EXACT_SUM_PAYLOAD,
        lo: 3.0,
        hi: 3.0,
        witness_hash: EXACT_WITNESS_HASH,
        plugin_hash: EXACT_PLUGIN_HASH,
    },
    Fixture {
        claim_id: "bounded-residual",
        family: fs_checker::plugins::BOUNDED_LINF_RESIDUAL_FAMILY,
        payload: RESIDUAL_2X1_PAYLOAD,
        lo: 0.0,
        hi: f64::from_bits(0x4000_0000_0000_0002),
        witness_hash: RESIDUAL_WITNESS_HASH,
        plugin_hash: RESIDUAL_PLUGIN_HASH,
    },
];

fn provenance() -> Provenance {
    Provenance::new("plugin-fixture-commit", "plugin-fixture-lock")
}

fn witness(fixture: Fixture) -> SemanticWitness {
    SemanticWitness::new(
        fixture.family,
        fs_checker::plugins::INITIAL_SEMANTIC_SCHEMA_VERSION,
        fixture.payload.to_vec(),
    )
}

fn package(fixture: Fixture, release_ready: bool) -> EvidencePackage {
    let claim = Claim::from_portable_certificate(
        fixture.claim_id,
        "portable semantic fixture",
        fixture.lo,
        fixture.hi,
        "fixture/portable-certificate",
        witness(fixture),
    );
    let claim = if release_ready {
        claim.with_falsifier(FalsifierRecord {
            name: "fixture-independent-probe".to_string(),
            attempts: 64,
            refuted: false,
            detail: "64 retained boundary probes found no violation".to_string(),
            artifact_hash: ARTIFACT_HASH.to_string(),
        })
    } else {
        claim
    };
    EvidencePackage::new(provenance()).with_claim(claim)
}

fn package_root(package: &EvidencePackage) -> ContentHash {
    package
        .try_merkle_root()
        .expect("bounded package fixture has a root")
}

fn package_json(package: &EvidencePackage) -> String {
    package
        .to_json()
        .expect("bounded package fixture has strict JSON")
}

fn flip(hash: ContentHash) -> ContentHash {
    let mut bytes = *hash.as_bytes();
    bytes[0] ^= 0x80;
    ContentHash(bytes)
}

struct FixtureSourceVerifier;

impl SourceCertificateVerifier for FixtureSourceVerifier {
    fn verify(&self, request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let accepted = request.package_provenance == &provenance()
            && request.package_root != ContentHash([0; 32])
            && request.claim_index == 0
            && request.statement == "portable semantic fixture"
            && request.claim_subject_hash != ContentHash([0; 32])
            && request.producer == "fixture/portable-certificate"
            && request.semantic_witness.is_some_and(|witness| {
                let exact_identity = match witness.family() {
                    fs_checker::plugins::EXACT_INTERVAL_FAMILY => {
                        request.claim_id == "exact-sum"
                            && request.lo.to_bits() == 3.0_f64.to_bits()
                            && request.hi.to_bits() == 3.0_f64.to_bits()
                            && request.certificate_hash.to_hex() == EXACT_WITNESS_HASH
                            && witness.canonical_payload() == EXACT_SUM_PAYLOAD
                    }
                    fs_checker::plugins::BOUNDED_LINF_RESIDUAL_FAMILY => {
                        request.claim_id == "bounded-residual"
                            && request.lo.to_bits() == 0.0_f64.to_bits()
                            && request.hi.to_bits() == 0x4000_0000_0000_0002
                            && request.certificate_hash.to_hex() == RESIDUAL_WITNESS_HASH
                            && witness.canonical_payload() == RESIDUAL_2X1_PAYLOAD
                    }
                    _ => false,
                };
                exact_identity
                    && witness.schema_version()
                        == fs_checker::plugins::INITIAL_SEMANTIC_SCHEMA_VERSION
                    && witness.content_hash() == request.certificate_hash
            });
        if accepted {
            VerificationDecision::accept(ContentHash([0x31; 32]))
        } else {
            VerificationDecision::reject(ContentHash([0x31; 32]))
        }
    }
}

struct FixtureFalsifierVerifier;

impl FalsifierVerifier for FixtureFalsifierVerifier {
    fn verify(&self, request: &FalsifierRequest<'_>) -> VerificationDecision {
        let accepted = request.package_provenance == &provenance()
            && request.package_root != ContentHash([0; 32])
            && request.claim_index == 0
            && matches!(request.claim_id, "exact-sum" | "bounded-residual")
            && request.statement == "portable semantic fixture"
            && request.claim_subject_hash != ContentHash([0; 32])
            && request.falsifier_index == 0
            && request.name == "fixture-independent-probe"
            && request.attempts == 64
            && !request.refuted
            && request.detail == "64 retained boundary probes found no violation"
            && request.artifact_hash.to_hex() == ARTIFACT_HASH;
        if accepted {
            VerificationDecision::accept(ContentHash([0x32; 32]))
        } else {
            VerificationDecision::reject(ContentHash([0x32; 32]))
        }
    }
}

struct ReleaseVerifier;

impl SignatureVerifier for ReleaseVerifier {
    fn verify(&self, request: &SignatureRequest<'_>) -> VerificationDecision {
        let accepted = request.signature == format!("release-fixture:{}", request.subject_hash())
            && matches!(
                request.purpose,
                SignaturePurpose::ReleaseApproval {
                    checker_protocol,
                    expected_root,
                    admission_context,
                    semantic_context,
                } if checker_protocol == CHECKER_PROTOCOL_VERSION
                    && expected_root == request.package_root
                    && admission_context != ContentHash([0; 32])
                    && semantic_context != ContentHash([0; 32])
            );
        if accepted {
            VerificationDecision::accept(ContentHash([0x33; 32]))
        } else {
            VerificationDecision::reject(ContentHash([0x33; 32]))
        }
    }
}

static SOURCE_VERIFIER: FixtureSourceVerifier = FixtureSourceVerifier;
static FALSIFIER_VERIFIER: FixtureFalsifierVerifier = FixtureFalsifierVerifier;

fn capabilities() -> VerificationCapabilities<'static> {
    VerificationCapabilities::deny_all()
        .with_source_certificates(&SOURCE_VERIFIER)
        .with_falsifiers(&FALSIFIER_VERIFIER)
}

fn signed_for_release(package: EvidencePackage, semantic_context: ContentHash) -> EvidencePackage {
    let root = package_root(&package);
    let verified = package
        .verify_with(&capabilities())
        .expect("origin and falsifier fixture authenticate");
    let purpose = SignaturePurpose::ReleaseApproval {
        checker_protocol: CHECKER_PROTOCOL_VERSION,
        expected_root: root,
        admission_context: verified.receipt().release_admission_context(),
        semantic_context,
    };
    package.signed(format!(
        "release-fixture:{}",
        signature_subject_hash(root, purpose)
    ))
}

#[test]
fn facade_exposes_the_one_closed_registry_and_pinned_identities() {
    let registry = fs_checker::plugins::semantic_plugin_registry();
    assert_eq!(registry, fs_checker::semantic_plugin_registry());
    assert_eq!(registry.len(), FIXTURES.len());
    assert_eq!(registry[0].family(), FIXTURES[0].family);
    assert_eq!(registry[1].family(), FIXTURES[1].family);
    assert!(registry.iter().all(|descriptor| {
        descriptor.schema_version() == fs_checker::plugins::INITIAL_SEMANTIC_SCHEMA_VERSION
            && descriptor.maximum_payload_bytes() == fs_checker::plugins::MAX_SEMANTIC_WITNESS_BYTES
    }));
    assert_eq!(registry[0].fingerprint().to_hex(), EXACT_PLUGIN_HASH);
    assert_eq!(registry[1].fingerprint().to_hex(), RESIDUAL_PLUGIN_HASH);

    let registry_fingerprint = fs_checker::plugins::semantic_registry_fingerprint();
    assert_eq!(registry_fingerprint.to_hex(), REGISTRY_HASH);
    assert_eq!(
        registry_fingerprint,
        fs_checker::semantic_registry_fingerprint()
    );

    for descriptor in registry {
        let fingerprint = descriptor.fingerprint();
        assert_eq!(
            SemanticPluginDescriptor::admit_retained_fingerprint(
                fs_checker::plugins::SEMANTIC_PLUGIN_IDENTITY_VERSION,
                fingerprint.as_bytes(),
            ),
            Some(fingerprint)
        );
        assert!(
            SemanticPluginDescriptor::admit_retained_fingerprint(
                fs_checker::plugins::SEMANTIC_PLUGIN_IDENTITY_VERSION + 1,
                fingerprint.as_bytes(),
            )
            .is_none()
        );
        assert!(
            SemanticPluginDescriptor::admit_retained_fingerprint(
                fs_checker::plugins::SEMANTIC_PLUGIN_IDENTITY_VERSION,
                &fingerprint.as_bytes()[..31],
            )
            .is_none()
        );
    }
}

#[test]
fn both_families_reach_package_json_standalone_and_release_gates() {
    for fixture in FIXTURES {
        let portable_witness = witness(fixture);
        assert_eq!(portable_witness.family(), fixture.family);
        assert_eq!(portable_witness.canonical_payload(), fixture.payload);
        assert_eq!(
            portable_witness.content_hash().to_hex(),
            fixture.witness_hash
        );

        let base_package = package(fixture, false);
        let root = package_root(&base_package);
        let semantic = fs_checker::plugins::verify_portable_semantics(&base_package);
        assert_eq!(semantic.status(), SemanticStatus::Verified);
        assert_eq!(semantic.package_root(), root);
        assert_eq!(semantic.registry_fingerprint().to_hex(), REGISTRY_HASH);
        assert_eq!(semantic.witnesses(), 1);
        assert_eq!(semantic.payload_bytes(), fixture.payload.len());
        assert!(semantic.operations() > 0);
        assert!(semantic.failures().is_empty());
        assert!(semantic.validate_context_hash());
        let receipt = &semantic.claims()[0];
        assert_eq!(receipt.claim_index(), 0);
        assert_eq!(receipt.claim_id(), fixture.claim_id);
        assert_eq!(receipt.family(), Some(fixture.family));
        assert_eq!(
            receipt.schema_version(),
            Some(fs_checker::plugins::INITIAL_SEMANTIC_SCHEMA_VERSION)
        );
        assert_eq!(receipt.status(), SemanticClaimStatus::Verified);
        assert_eq!(
            receipt.witness_hash().map(|hash| hash.to_hex()),
            Some(fixture.witness_hash.to_string())
        );
        assert_eq!(
            receipt.plugin_fingerprint().map(|hash| hash.to_hex()),
            Some(fixture.plugin_hash.to_string())
        );

        let in_memory = check_with_capabilities(&base_package, Some(root), None, &capabilities());
        assert!(in_memory.passed(), "{:?}", in_memory.findings());
        assert_eq!(in_memory.integrity_status(), IntegrityStatus::Verified);
        assert_eq!(in_memory.semantic_status(), SemanticStatus::Verified);
        assert_eq!(in_memory.origin_status(), OriginStatus::Authenticated);
        assert_eq!(in_memory.semantic_report(), &semantic);
        assert_eq!(
            check_json_with_capabilities(
                &package_json(&base_package),
                Some(root),
                None,
                &capabilities(),
            ),
            in_memory
        );

        let release_package = package(fixture, true);
        let release_semantics = fs_checker::plugins::verify_portable_semantics(&release_package);
        assert_eq!(release_semantics.status(), SemanticStatus::Verified);
        let signed = signed_for_release(release_package, release_semantics.context_hash());
        let signed_root = package_root(&signed);
        let release = check_for_release_with_capabilities(
            &signed,
            signed_root,
            &ReleaseVerifier,
            &capabilities(),
        );
        assert!(release.release_admitted(), "{:?}", release.findings());
        assert!(release.release_independently_verified());
        assert_eq!(release.semantic_report(), &release_semantics);
        assert_eq!(
            check_json_for_release_with_capabilities(
                &package_json(&signed),
                signed_root,
                &ReleaseVerifier,
                &capabilities(),
            ),
            release
        );
    }
}

#[test]
fn package_root_family_schema_and_payload_substitution_fail_closed() {
    let build = |family: &str, schema_version: u32, payload: Vec<u8>| {
        EvidencePackage::new(provenance()).with_claim(Claim::from_portable_certificate(
            "tamper-target",
            "portable semantic fixture",
            3.0,
            3.0,
            "fixture/portable-certificate",
            SemanticWitness::new(family, schema_version, payload),
        ))
    };
    let baseline = build(
        fs_checker::plugins::EXACT_INTERVAL_FAMILY,
        fs_checker::plugins::INITIAL_SEMANTIC_SCHEMA_VERSION,
        EXACT_SUM_PAYLOAD.to_vec(),
    );
    let baseline_root = package_root(&baseline);

    let unknown_family = build(
        "frankensim/unknown-portable-proof",
        fs_checker::plugins::INITIAL_SEMANTIC_SCHEMA_VERSION,
        EXACT_SUM_PAYLOAD.to_vec(),
    );
    let unsupported_schema = build(
        fs_checker::plugins::EXACT_INTERVAL_FAMILY,
        fs_checker::plugins::INITIAL_SEMANTIC_SCHEMA_VERSION + 1,
        EXACT_SUM_PAYLOAD.to_vec(),
    );
    let mut changed_math = EXACT_SUM_PAYLOAD.to_vec();
    changed_math[14] = 3; // the second exact leaf: 1 + 2 becomes 1 + 3
    let changed_payload = build(
        fs_checker::plugins::EXACT_INTERVAL_FAMILY,
        fs_checker::plugins::INITIAL_SEMANTIC_SCHEMA_VERSION,
        changed_math,
    );

    let changed_roots = [
        package_root(&unknown_family),
        package_root(&unsupported_schema),
        package_root(&changed_payload),
    ];
    assert!(changed_roots.iter().all(|root| *root != baseline_root));
    assert_ne!(changed_roots[0], changed_roots[1]);
    assert_ne!(changed_roots[0], changed_roots[2]);
    assert_ne!(changed_roots[1], changed_roots[2]);

    // An expected-root mismatch stops before semantic bytes or callbacks are
    // trusted, even when the substituted payload is otherwise well-formed.
    let root_refusal =
        check_with_capabilities(&changed_payload, Some(baseline_root), None, &capabilities());
    assert_eq!(root_refusal.integrity_status(), IntegrityStatus::Refused);
    assert_eq!(root_refusal.semantic_status(), SemanticStatus::NotRun);
    assert_eq!(root_refusal.origin_status(), OriginStatus::NotRun);
    assert!(
        root_refusal
            .findings()
            .iter()
            .any(|finding| finding.kind == "content-address-mismatch")
    );
    assert_eq!(
        check_json_with_capabilities(
            &package_json(&changed_payload),
            Some(baseline_root),
            None,
            &capabilities(),
        ),
        root_refusal
    );

    for (package, expected_kind) in [
        (&unknown_family, SemanticFailureKind::UnknownFamily),
        (&unsupported_schema, SemanticFailureKind::UnsupportedVersion),
        (&changed_payload, SemanticFailureKind::ClaimMismatch),
    ] {
        let semantic = fs_checker::plugins::verify_portable_semantics(package);
        assert_eq!(semantic.status(), SemanticStatus::Refused);
        assert_eq!(semantic.registry_fingerprint().to_hex(), REGISTRY_HASH);
        assert!(semantic.validate_context_hash());
        assert_eq!(semantic.failures()[0].kind(), expected_kind);
        assert_eq!(semantic.claims()[0].status(), SemanticClaimStatus::Refused);

        let root = package_root(package);
        let report = check_with_capabilities(package, Some(root), None, &capabilities());
        assert_eq!(report.integrity_status(), IntegrityStatus::Verified);
        assert_eq!(report.semantic_status(), SemanticStatus::Refused);
        assert_eq!(report.origin_status(), OriginStatus::NotRun);
        assert_eq!(report.semantic_report(), &semantic);
        assert_eq!(
            check_json_with_capabilities(&package_json(package), Some(root), None, &capabilities(),),
            report
        );
    }
}

#[test]
fn registry_and_context_identity_transports_are_exact_version_and_width() {
    let registry = fs_checker::plugins::semantic_registry_fingerprint();
    assert_eq!(
        fs_checker::plugins::admit_retained_semantic_registry_fingerprint(
            fs_checker::plugins::SEMANTIC_REGISTRY_IDENTITY_VERSION,
            registry.as_bytes(),
        ),
        Some(registry)
    );
    assert!(
        fs_checker::plugins::admit_retained_semantic_registry_fingerprint(
            fs_checker::plugins::SEMANTIC_REGISTRY_IDENTITY_VERSION + 1,
            registry.as_bytes(),
        )
        .is_none()
    );
    assert!(
        fs_checker::plugins::admit_retained_semantic_registry_fingerprint(
            fs_checker::plugins::SEMANTIC_REGISTRY_IDENTITY_VERSION,
            &registry.as_bytes()[..31],
        )
        .is_none()
    );
    let mut extended_registry = registry.as_bytes().to_vec();
    extended_registry.push(0);
    assert!(
        fs_checker::plugins::admit_retained_semantic_registry_fingerprint(
            fs_checker::plugins::SEMANTIC_REGISTRY_IDENTITY_VERSION,
            &extended_registry,
        )
        .is_none()
    );

    // The transport guard validates version and shape, not authority. An exact
    // width mutation remains an opaque retained hash and must still be compared
    // against the compiled registry identity before it can be trusted.
    let foreign_registry = flip(registry);
    assert_eq!(
        fs_checker::plugins::admit_retained_semantic_registry_fingerprint(
            fs_checker::plugins::SEMANTIC_REGISTRY_IDENTITY_VERSION,
            foreign_registry.as_bytes(),
        ),
        Some(foreign_registry)
    );
    assert_ne!(foreign_registry, registry);

    let report = fs_checker::plugins::verify_portable_semantics(&package(FIXTURES[0], false));
    let context = report.context_hash();
    assert!(report.validate_context_hash());
    assert_eq!(report.registry_fingerprint(), registry);
    assert_eq!(
        SemanticReport::admit_retained_context_hash(
            fs_checker::plugins::SEMANTIC_REPORT_IDENTITY_VERSION,
            context.as_bytes(),
        ),
        Some(context)
    );
    assert!(
        SemanticReport::admit_retained_context_hash(
            fs_checker::plugins::SEMANTIC_REPORT_IDENTITY_VERSION + 1,
            context.as_bytes(),
        )
        .is_none()
    );
    let mut extended_context = context.as_bytes().to_vec();
    extended_context.push(0);
    assert!(
        SemanticReport::admit_retained_context_hash(
            fs_checker::plugins::SEMANTIC_REPORT_IDENTITY_VERSION,
            &extended_context,
        )
        .is_none()
    );
    assert!(
        SemanticReport::admit_retained_context_hash(
            fs_checker::plugins::SEMANTIC_REPORT_IDENTITY_VERSION,
            &context.as_bytes()[..31],
        )
        .is_none()
    );
    assert_ne!(
        context,
        fs_checker::plugins::verify_portable_semantics(&package(FIXTURES[1], false)).context_hash()
    );
}

#[test]
fn release_signature_cannot_replay_a_different_semantic_context() {
    let unsigned = package(FIXTURES[0], true);
    let root = package_root(&unsigned);
    let semantics = fs_checker::plugins::verify_portable_semantics(&unsigned);
    assert_eq!(semantics.status(), SemanticStatus::Verified);

    let wrong_context = flip(semantics.context_hash());
    let replayed = signed_for_release(unsigned.clone(), wrong_context);
    assert_eq!(package_root(&replayed), root);
    let refused =
        check_for_release_with_capabilities(&replayed, root, &ReleaseVerifier, &capabilities());
    assert_eq!(refused.integrity_status(), IntegrityStatus::Verified);
    assert_eq!(refused.semantic_status(), SemanticStatus::Verified);
    assert_eq!(refused.origin_status(), OriginStatus::Authenticated);
    assert!(!refused.release_admitted());
    assert!(
        refused
            .findings()
            .iter()
            .any(|finding| finding.kind == "signature-invalid")
    );
    assert_eq!(
        check_json_for_release_with_capabilities(
            &package_json(&replayed),
            root,
            &ReleaseVerifier,
            &capabilities(),
        ),
        refused
    );

    let correctly_signed = signed_for_release(unsigned, semantics.context_hash());
    let admitted = check_for_release_with_capabilities(
        &correctly_signed,
        root,
        &ReleaseVerifier,
        &capabilities(),
    );
    assert!(admitted.release_admitted(), "{:?}", admitted.findings());
    assert!(admitted.release_independently_verified());
}

#[test]
fn structurally_valid_but_mathematically_false_claims_refuse_for_both_families() {
    let false_packages = [
        EvidencePackage::new(provenance()).with_claim(Claim::from_portable_certificate(
            "false-exact-sum",
            "portable semantic fixture",
            4.0,
            4.0,
            "fixture/portable-certificate",
            witness(FIXTURES[0]),
        )),
        EvidencePackage::new(provenance()).with_claim(Claim::from_portable_certificate(
            "false-residual-bound",
            "portable semantic fixture",
            0.0,
            2.0,
            "fixture/portable-certificate",
            witness(FIXTURES[1]),
        )),
    ];

    for package in false_packages {
        let root = package_root(&package);
        let semantic = fs_checker::plugins::verify_portable_semantics(&package);
        assert_eq!(semantic.status(), SemanticStatus::Refused);
        assert_eq!(semantic.failures().len(), 1);
        assert_eq!(
            semantic.failures()[0].kind(),
            SemanticFailureKind::ClaimMismatch
        );
        assert_eq!(semantic.failures()[0].claim_index(), Some(0));
        assert!(semantic.failures()[0].family().is_some());
        assert!(semantic.validate_context_hash());

        let report = check_with_capabilities(&package, Some(root), None, &capabilities());
        assert_eq!(report.integrity_status(), IntegrityStatus::Verified);
        assert_eq!(report.semantic_status(), SemanticStatus::Refused);
        assert_eq!(report.origin_status(), OriginStatus::NotRun);
        assert!(!report.passed());
        assert!(report.receipt().is_none());
        assert_eq!(
            check_json_with_capabilities(
                &package_json(&package),
                Some(root),
                None,
                &capabilities(),
            ),
            report
        );
    }
}
