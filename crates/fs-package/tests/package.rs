//! Battery for evidence packages (addendum Proposal 12). Covers a complete
//! mixed-color package, the all-estimated boundary (still valid, round-trips),
//! completeness failures (validated claim missing regime / dataset, verified
//! claim with a bad interval), Merkle content-addressing (determinism + tamper
//! detection), the format-version gate, optional signature, the color
//! breakdown, and deterministic JSON.

use fs_evidence::{Color, IntervalOp, ValidityDomain};
use fs_package::{
    AdmissionClass, AnchoredSourceRequest, AnchoredSourceVerifier, Claim, ContentHash,
    DerivationRequest, DerivationVerifier, EvidencePackage, FalsifierRecord, FalsifierRequest,
    FalsifierVerifier, MAX_JSON_CONTAINER_ITEMS, MAX_JSON_DEPTH, MAX_JSON_NUMBER_BYTES,
    MAX_JSON_STRING_BYTES, PackageError, PackageReport, Provenance, SignaturePurpose,
    SignatureRequest, SignatureStatus, SignatureVerifier, SourceCertificateRequest,
    SourceCertificateVerifier, VerificationCapabilities, VerificationDecision, WaiverGrant,
    WaiverVerifier, signature_subject_hash,
};

const CANONICAL_DATASET_HASH: &str =
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

fn package_root(package: &EvidencePackage) -> ContentHash {
    package
        .try_merkle_root()
        .expect("bounded test package has a content root")
}

fn package_json(package: &EvidencePackage) -> String {
    package.to_json().expect("bounded test package serializes")
}

struct FixtureSourceVerifier;
struct FixtureAnchorVerifier;
struct FixtureFalsifierVerifier;
struct FixtureDerivationVerifier;

fn fixture_policy(label: &str) -> ContentHash {
    fs_blake3::hash_domain("fs-package:test:v6:policy", label.as_bytes())
}

fn fixture_source_hash(
    provenance: &Provenance,
    claim_index: usize,
    claim_id: &str,
    statement: &str,
    lo: f64,
    hi: f64,
    producer: &str,
) -> String {
    use core::fmt::Write as _;

    fn push_atom(out: &mut String, atom: &str) {
        use core::fmt::Write as _;
        let _ = write!(out, "{}:{atom}|", atom.len());
    }

    let mut subject = String::from("fs-package:test:source-certificate-subject|");
    push_atom(&mut subject, &provenance.code_version);
    push_atom(&mut subject, &provenance.constellation_lock);
    let _ = write!(subject, "index:{claim_index}|");
    push_atom(&mut subject, claim_id);
    push_atom(&mut subject, statement);
    let _ = write!(subject, "lo:{}|hi:{}|", lo.to_bits(), hi.to_bits());
    push_atom(&mut subject, producer);
    fs_blake3::hash_domain("fs-package:test:source-certificate", subject.as_bytes()).to_hex()
}

impl SourceCertificateVerifier for FixtureSourceVerifier {
    fn verify(&self, request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let accepted = request.certificate_hash.to_hex()
            == fixture_source_hash(
                request.package_provenance,
                request.claim_index,
                request.claim_id,
                request.statement,
                request.lo,
                request.hi,
                request.producer,
            );
        VerificationDecision {
            accepted,
            policy_fingerprint: fixture_policy("fixture-source-verifier"),
        }
    }
}

impl AnchoredSourceVerifier for FixtureAnchorVerifier {
    fn verify(&self, request: &AnchoredSourceRequest<'_>) -> VerificationDecision {
        let accepted = request.content_hash.to_hex() == CANONICAL_DATASET_HASH
            && !request.dataset_id.is_empty()
            && !request.regime.bounds().is_empty();
        VerificationDecision {
            accepted,
            policy_fingerprint: fixture_policy("fixture-anchor-verifier"),
        }
    }
}

impl FalsifierVerifier for FixtureFalsifierVerifier {
    fn verify(&self, request: &FalsifierRequest<'_>) -> VerificationDecision {
        VerificationDecision {
            accepted: request.artifact_hash.to_hex() == CANONICAL_DATASET_HASH,
            policy_fingerprint: fixture_policy("fixture-falsifier-verifier"),
        }
    }
}

impl DerivationVerifier for FixtureDerivationVerifier {
    fn verify(&self, request: &DerivationRequest<'_>) -> VerificationDecision {
        VerificationDecision {
            accepted: request.artifact_hash.to_hex() == CANONICAL_DATASET_HASH,
            policy_fingerprint: fixture_policy("fixture-derivation-verifier"),
        }
    }
}

static FIXTURE_SOURCE_VERIFIER: FixtureSourceVerifier = FixtureSourceVerifier;
static FIXTURE_ANCHOR_VERIFIER: FixtureAnchorVerifier = FixtureAnchorVerifier;
static FIXTURE_FALSIFIER_VERIFIER: FixtureFalsifierVerifier = FixtureFalsifierVerifier;
static FIXTURE_DERIVATION_VERIFIER: FixtureDerivationVerifier = FixtureDerivationVerifier;

fn source_capabilities() -> VerificationCapabilities<'static> {
    VerificationCapabilities::deny_all()
        .with_source_certificates(&FIXTURE_SOURCE_VERIFIER)
        .with_anchored_sources(&FIXTURE_ANCHOR_VERIFIER)
        .with_falsifiers(&FIXTURE_FALSIFIER_VERIFIER)
        .with_derivations(&FIXTURE_DERIVATION_VERIFIER)
}

fn verify_package(pkg: &EvidencePackage) -> Result<PackageReport, PackageError> {
    pkg.verify_with(&source_capabilities())
}

fn prov() -> Provenance {
    Provenance::new("commit-abc123", "lock-deadbeef")
}

fn source_claim(
    provenance: &Provenance,
    claim_index: usize,
    id: &str,
    statement: &str,
    lo: f64,
    hi: f64,
) -> Claim {
    let producer = "test-solver/cert";
    Claim::from_certificate(
        id,
        statement,
        lo,
        hi,
        producer,
        fixture_source_hash(provenance, claim_index, id, statement, lo, hi, producer),
    )
}

fn verified(id: &str) -> Claim {
    source_claim(
        &prov(),
        0,
        id,
        &format!("{id}: stress <= sigma*"),
        -1.0,
        1.0,
    )
}
fn estimated(id: &str) -> Claim {
    Claim::estimated(id, format!("{id}: surrogate says ok"), "surrogate", 2.0)
}
fn validated(id: &str, regime: ValidityDomain, dataset: &str) -> Claim {
    Claim::anchored(
        id,
        format!("{id}: matches data"),
        regime,
        dataset,
        CANONICAL_DATASET_HASH,
    )
}
fn good_regime() -> ValidityDomain {
    ValidityDomain::unconstrained().with("Re", 1e5, 3e5)
}

fn derived_package() -> EvidencePackage {
    let verified_color = |lo: f64, hi: f64| Color::Verified { lo, hi };
    let provenance = Provenance::new("v", "l");
    EvidencePackage::new(provenance.clone())
        .with_claim(source_claim(&provenance, 0, "a", "left", 1.0, 2.0))
        .with_claim(source_claim(&provenance, 1, "b", "right", 10.0, 20.0))
        .with_claim(
            Claim::derived(
                "c",
                "sum",
                fs_evidence::compose(
                    &verified_color(1.0, 2.0),
                    &verified_color(10.0, 20.0),
                    IntervalOp::Add,
                ),
                vec![0, 1],
                IntervalOp::Add,
                CANONICAL_DATASET_HASH,
            )
            .with_falsifier(FalsifierRecord {
                name: "interval-probe".to_string(),
                attempts: 512,
                refuted: false,
                detail: "no violation found".to_string(),
                artifact_hash: CANONICAL_DATASET_HASH.to_string(),
            })
            .with_anchor("wt-2026-run9", CANONICAL_DATASET_HASH),
        )
}

fn assert_serialized_refuses(pkg: &EvidencePackage) {
    if let Ok(json) = pkg.to_json() {
        assert!(
            EvidencePackage::from_json(&json).is_err(),
            "serialized package must refuse the same semantics as verify(): {pkg:?}"
        );
    }
}

#[test]
fn fixture_source_authority_binds_every_typed_request_field() {
    let provenance = prov();
    let other_provenance = Provenance::new("different-commit", "lock-deadbeef");
    let hash = fixture_source_hash(
        &provenance,
        0,
        "claim",
        "certified subject",
        1.0,
        2.0,
        "test-solver/cert",
    );
    let request = SourceCertificateRequest {
        package_provenance: &provenance,
        claim_index: 0,
        claim_id: "claim",
        statement: "certified subject",
        lo: 1.0,
        hi: 2.0,
        producer: "test-solver/cert",
        certificate_hash: fs_package::ContentHash::from_hex(&hash).expect("fixture hash"),
    };
    assert!(FIXTURE_SOURCE_VERIFIER.verify(&request).accepted);
    for altered in [
        SourceCertificateRequest {
            package_provenance: &other_provenance,
            ..request
        },
        SourceCertificateRequest {
            claim_index: 1,
            ..request
        },
        SourceCertificateRequest {
            claim_id: "other-claim",
            ..request
        },
        SourceCertificateRequest {
            statement: "different subject",
            ..request
        },
        SourceCertificateRequest { lo: 0.0, ..request },
        SourceCertificateRequest { hi: 3.0, ..request },
        SourceCertificateRequest {
            producer: "other-solver/cert",
            ..request
        },
    ] {
        assert!(
            !FIXTURE_SOURCE_VERIFIER.verify(&altered).accepted,
            "fixture authority accepted a modified typed request: {altered:?}"
        );
    }
}

#[test]
fn a_complete_mixed_color_package_verifies() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(validated("c2", good_regime(), "wind-tunnel-2026"))
        .with_claim(estimated("c3"));
    let report = verify_package(&pkg).expect("complete package verifies");
    assert_eq!(report.claims(), 3);
    assert_eq!(report.breakdown().verified, 1);
    assert_eq!(report.breakdown().validated, 1);
    assert_eq!(report.breakdown().estimated, 1);
    assert_eq!(report.merkle_root(), package_root(&pkg));
    assert_ne!(report.merkle_root().as_bytes(), &[0u8; 32]);
}

#[test]
fn verified_package_binding_recomputes_admissions_and_summaries() {
    let package = EvidencePackage::new(prov())
        .with_claim(estimated("estimate-a"))
        .with_claim(estimated("estimate-b"));
    let verified = package
        .into_verified()
        .expect("deny-all estimated package verifies");
    assert!(verified.validate_binding());
    assert_eq!(verified.admitted_claims().len(), 2);
    assert!(
        verified
            .admitted_claims()
            .all(|claim| claim.scientific_color().is_some())
    );
}

#[test]
fn an_all_estimated_package_is_still_valid_and_round_trips() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(estimated("e1"))
        .with_claim(estimated("e2"));
    let report = pkg.verify().expect("all-estimated is honest, not invalid");
    assert_eq!(report.breakdown().estimated, 2);
    assert_eq!(report.breakdown().verified, 0);
    let json = package_json(&pkg);
    assert!(json.contains("\"estimated\""));
    assert!(json.contains("e1") && json.contains("e2"));
}

#[test]
fn a_validated_claim_missing_its_regime_fails_completeness() {
    // an unconstrained (empty) regime = no regime tag.
    let pkg = EvidencePackage::new(prov()).with_claim(validated(
        "v",
        ValidityDomain::unconstrained(),
        "some-data",
    ));
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::IncompleteValidatedClaim {
            missing: "regime",
            ..
        })
    ));
}

#[test]
fn a_validated_claim_missing_its_dataset_fails_completeness() {
    let pkg = EvidencePackage::new(prov()).with_claim(validated("v", good_regime(), "   "));
    // Schema v5: the sealed constructor anchors the claim to its own
    // dataset, so a blank dataset now refuses at the anchor-record gate
    // (earlier, and just as closed) rather than at color completeness.
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::InvalidAnchorRecord {
            field: "dataset_id",
            ..
        })
    ));
}

#[test]
fn a_verified_claim_with_a_bad_interval_fails() {
    let pkg = EvidencePackage::new(prov()).with_claim(Claim::from_certificate(
        "v",
        "backwards",
        5.0,
        1.0,
        "test-solver/cert",
        CANONICAL_DATASET_HASH,
    ));
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::IncompleteVerifiedClaim { .. })
    ));
}

#[test]
fn in_memory_and_serialized_provenance_and_identity_gates_are_identical() {
    let bad_code = EvidencePackage::new(Provenance::new(" ", "lock")).with_claim(verified("c"));
    assert!(matches!(
        bad_code.verify(),
        Err(PackageError::IncompleteProvenance {
            missing: "code_version"
        })
    ));
    assert_serialized_refuses(&bad_code);

    let bad_lock = EvidencePackage::new(Provenance::new("commit", "\t")).with_claim(verified("c"));
    assert!(matches!(
        bad_lock.verify(),
        Err(PackageError::IncompleteProvenance {
            missing: "constellation_lock"
        })
    ));
    assert_serialized_refuses(&bad_lock);

    let blank_id = EvidencePackage::new(prov()).with_claim(verified("  "));
    assert!(matches!(
        blank_id.verify(),
        Err(PackageError::InvalidClaimId {
            index: 0,
            reason: "blank",
            ..
        })
    ));
    assert_serialized_refuses(&blank_id);

    let duplicate_id = EvidencePackage::new(prov())
        .with_claim(verified("same"))
        .with_claim(estimated("same"));
    assert!(matches!(
        duplicate_id.verify(),
        Err(PackageError::InvalidClaimId {
            index: 1,
            reason: "duplicate",
            ..
        })
    ));
    assert_serialized_refuses(&duplicate_id);
}

#[test]
fn in_memory_and_serialized_claim_statements_must_be_meaningful() {
    for (statement, reason) in [
        (" \t", "blank"),
        ("TODO", "placeholder"),
        (" n/A ", "placeholder"),
    ] {
        let pkg = EvidencePackage::new(prov()).with_claim(Claim::from_certificate(
            "claim",
            statement,
            0.0,
            1.0,
            "test-solver/cert",
            CANONICAL_DATASET_HASH,
        ));
        assert!(matches!(
            pkg.verify(),
            Err(PackageError::InvalidClaimStatement {
                claim,
                reason: found,
            }) if claim == "claim" && found == reason
        ));
        assert_serialized_refuses(&pkg);
    }
}

#[test]
fn in_memory_and_serialized_estimate_gates_are_identical() {
    let blank_estimator = EvidencePackage::new(prov()).with_claim(Claim::estimated(
        "e",
        "missing estimator identity",
        " ",
        1.0,
    ));
    assert!(matches!(
        blank_estimator.verify(),
        Err(PackageError::IncompleteEstimatedClaim {
            missing: "estimator",
            ..
        })
    ));
    assert_serialized_refuses(&blank_estimator);

    for dispersion in [-1.0, f64::NEG_INFINITY, f64::NAN] {
        let pkg = EvidencePackage::new(prov()).with_claim(Claim::estimated(
            "e",
            "invalid dispersion",
            "probe",
            dispersion,
        ));
        assert!(matches!(
            pkg.verify(),
            Err(PackageError::InvalidEstimatedDispersion { .. })
        ));
        assert_serialized_refuses(&pkg);
    }

    let explicitly_unbounded = EvidencePackage::new(prov()).with_claim(Claim::estimated(
        "unbounded",
        "honest no-spread-claim sentinel",
        "regime-exit",
        f64::INFINITY,
    ));
    assert!(explicitly_unbounded.verify().is_ok());
    assert!(
        explicitly_unbounded
            .magnitude_budget()
            .expect("estimated-only package verifies")
            .estimated_dispersion
            .is_infinite()
    );
    assert_eq!(
        EvidencePackage::from_json(&package_json(&explicitly_unbounded)).unwrap(),
        explicitly_unbounded
    );
}

#[test]
fn in_memory_and_serialized_magnitude_gates_are_identical() {
    let width_overflow = EvidencePackage::new(prov()).with_claim(Claim::from_certificate(
        "wide",
        "finite endpoints whose width overflows",
        -f64::MAX,
        f64::MAX,
        "test-solver/cert",
        CANONICAL_DATASET_HASH,
    ));
    assert!(matches!(
        width_overflow.verify(),
        Err(PackageError::MagnitudeOverflow {
            component: "verified_width",
            ..
        })
    ));
    assert_serialized_refuses(&width_overflow);

    let dispersion_overflow = EvidencePackage::new(prov())
        .with_claim(Claim::estimated(
            "d1",
            "large finite dispersion",
            "probe-1",
            f64::MAX,
        ))
        .with_claim(Claim::estimated(
            "d2",
            "second large finite dispersion",
            "probe-2",
            f64::MAX,
        ));
    assert!(matches!(
        dispersion_overflow.verify(),
        Err(PackageError::MagnitudeOverflow {
            component: "estimated_dispersion",
            ..
        })
    ));
    assert_serialized_refuses(&dispersion_overflow);

    let large = f64::MAX * 0.75;
    let cross_component_overflow = EvidencePackage::new(prov())
        .with_claim(Claim::from_certificate(
            "wide-finite",
            "large but finite interval width",
            0.0,
            large,
            "test-solver/cert",
            CANONICAL_DATASET_HASH,
        ))
        .with_claim(Claim::estimated(
            "spread-finite",
            "large but finite estimated spread",
            "probe",
            large,
        ));
    assert!(matches!(
        cross_component_overflow.verify(),
        Err(PackageError::MagnitudeOverflow {
            component: "quantified_total",
            ..
        })
    ));
    assert_serialized_refuses(&cross_component_overflow);
}

#[test]
fn in_memory_and_serialized_regime_gates_are_identical() {
    for regime in [
        ValidityDomain::unconstrained().with(" ", 1.0, 2.0),
        ValidityDomain::unconstrained().with("Re", 1.0, f64::INFINITY),
    ] {
        let pkg = EvidencePackage::new(prov()).with_claim(validated("v", regime, "dataset"));
        assert!(matches!(
            pkg.verify(),
            Err(PackageError::InvalidValidatedRegime { .. })
        ));
        assert_serialized_refuses(&pkg);
    }

    let ordered = EvidencePackage::new(prov()).with_claim(validated(
        "ordered",
        ValidityDomain::unconstrained().with("Re", 1.0, 2.0),
        "dataset",
    ));
    let ordered_pair = format!(
        "[\"{:016x}\",\"{:016x}\"]",
        1.0_f64.to_bits(),
        2.0_f64.to_bits()
    );
    let inverted_pair = format!(
        "[\"{:016x}\",\"{:016x}\"]",
        2.0_f64.to_bits(),
        1.0_f64.to_bits()
    );
    let inverted_json = package_json(&ordered).replace(&ordered_pair, &inverted_pair);
    let err = EvidencePackage::from_json(&inverted_json).expect_err("inverted regime refuses");
    assert!(err.why.contains("inverted bounds"), "{err}");
}

#[test]
fn falsifier_attempt_counts_round_trip_without_f64_precision_loss() {
    let first_unrepresentable_integer = (1_u64 << 53) + 1;
    let pkg = EvidencePackage::new(prov()).with_claim(
        verified("count")
            .with_falsifier(FalsifierRecord {
                name: "precision-probe".to_string(),
                attempts: first_unrepresentable_integer,
                refuted: false,
                detail: "one above f64's exact-integer range".to_string(),
                artifact_hash: CANONICAL_DATASET_HASH.to_string(),
            })
            .with_falsifier(FalsifierRecord {
                name: "exhaustive-probe".to_string(),
                attempts: u64::MAX,
                refuted: false,
                detail: "full-width counter".to_string(),
                artifact_hash: CANONICAL_DATASET_HASH.to_string(),
            }),
    );
    let json = package_json(&pkg);
    assert!(json.contains("\"attempts\":9007199254740993"));
    assert!(json.contains("\"attempts\":18446744073709551615"));
    let back = EvidencePackage::from_json(&json).expect("full-width u64 values parse exactly");
    assert_eq!(back, pkg);
    assert_eq!(
        back.declared_claims_unverified()[0].declared_falsifiers_unverified()[0].attempts,
        first_unrepresentable_integer
    );
    assert_eq!(
        back.declared_claims_unverified()[0].declared_falsifiers_unverified()[1].attempts,
        u64::MAX
    );

    let overflow = json.replace("18446744073709551615", "18446744073709551616");
    let err = EvidencePackage::from_json(&overflow).expect_err("u64 overflow must refuse");
    assert!(err.why.contains("out of range"), "{err}");
}

#[test]
fn falsifier_records_require_identity_work_and_outcome() {
    let assert_invalid = |record: FalsifierRecord, expected_field: &'static str| {
        let pkg = EvidencePackage::new(prov()).with_claim(verified("f").with_falsifier(record));
        assert!(matches!(
            pkg.verify(),
            Err(PackageError::InvalidFalsifierRecord {
                falsifier: 0,
                field,
                ..
            }) if field == expected_field
        ));
        assert_serialized_refuses(&pkg);
    };

    assert_invalid(
        FalsifierRecord {
            name: " \t".to_string(),
            attempts: 1,
            refuted: false,
            detail: "no violation".to_string(),
            artifact_hash: CANONICAL_DATASET_HASH.to_string(),
        },
        "name",
    );
    assert_invalid(
        FalsifierRecord {
            name: "interval-probe".to_string(),
            attempts: 0,
            refuted: false,
            detail: "no work ran".to_string(),
            artifact_hash: CANONICAL_DATASET_HASH.to_string(),
        },
        "attempts",
    );
    assert_invalid(
        FalsifierRecord {
            name: "interval-probe".to_string(),
            attempts: 1,
            refuted: false,
            detail: "  ".to_string(),
            artifact_hash: CANONICAL_DATASET_HASH.to_string(),
        },
        "detail",
    );
    assert_invalid(
        FalsifierRecord {
            name: "interval-probe".to_string(),
            attempts: 1,
            refuted: false,
            detail: "no violation".to_string(),
            artifact_hash: "deadbeef".to_string(),
        },
        "artifact_hash",
    );
    for placeholder in ["TODO", " placeholder ", "N/A", "not run", "unknown"] {
        assert_invalid(
            FalsifierRecord {
                name: placeholder.to_string(),
                attempts: 1,
                refuted: false,
                detail: "no violation observed".to_string(),
                artifact_hash: CANONICAL_DATASET_HASH.to_string(),
            },
            "name",
        );
        assert_invalid(
            FalsifierRecord {
                name: "interval-probe".to_string(),
                attempts: 1,
                refuted: false,
                detail: placeholder.to_string(),
                artifact_hash: CANONICAL_DATASET_HASH.to_string(),
            },
            "detail",
        );
    }
}

#[test]
fn anchor_records_require_identity_and_canonical_content_hash() {
    let cases = [
        (" ", CANONICAL_DATASET_HASH, "dataset_id"),
        ("dataset", "deadbeef", "content_hash"),
        (
            "dataset",
            "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdeg",
            "content_hash",
        ),
        (
            "dataset",
            "0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF0123456789ABCDEF",
            "content_hash",
        ),
    ];
    for (dataset_id, content_hash, expected_field) in cases {
        let claim = verified("a").with_anchor(dataset_id.to_string(), content_hash.to_string());
        let pkg = EvidencePackage::new(prov()).with_claim(claim);
        assert!(matches!(
            pkg.verify(),
            Err(PackageError::InvalidAnchorRecord {
                anchor: 0,
                field,
                ..
            }) if field == expected_field
        ));
        assert_serialized_refuses(&pkg);
    }
}

#[test]
fn untrusted_json_resource_limits_fail_before_schema_mapping() {
    let nested = format!(
        "{}null{}",
        "[".repeat(MAX_JSON_DEPTH + 1),
        "]".repeat(MAX_JSON_DEPTH + 1)
    );
    let err = EvidencePackage::from_json(&nested).expect_err("deep nesting must refuse");
    assert!(err.why.contains("nesting depth"), "{err}");

    let oversized_string = format!("\"{}\"", "x".repeat(MAX_JSON_STRING_BYTES + 1));
    let err =
        EvidencePackage::from_json(&oversized_string).expect_err("oversized string must refuse");
    assert!(err.why.contains("decoded string"), "{err}");

    let oversized_number = "9".repeat(MAX_JSON_NUMBER_BYTES + 1);
    let err =
        EvidencePackage::from_json(&oversized_number).expect_err("long number token must refuse");
    assert!(err.why.contains("number token"), "{err}");

    let oversized_array = format!(
        "[{}]",
        std::iter::repeat_n("null", MAX_JSON_CONTAINER_ITEMS + 1)
            .collect::<Vec<_>>()
            .join(",")
    );
    let err =
        EvidencePackage::from_json(&oversized_array).expect_err("oversized array must refuse");
    assert!(err.why.contains("array element count"), "{err}");

    let oversized_in_memory = EvidencePackage::new(prov()).with_claim(Claim::estimated(
        "large",
        "x".repeat(MAX_JSON_STRING_BYTES + 1),
        "probe",
        1.0,
    ));
    assert!(matches!(
        oversized_in_memory.verify(),
        Err(PackageError::TransportLimit { .. })
    ));
    assert!(matches!(
        oversized_in_memory.try_merkle_root(),
        Err(PackageError::TransportLimit { .. })
    ));
    assert!(matches!(
        oversized_in_memory.to_json(),
        Err(PackageError::TransportLimit { .. })
    ));
    assert_serialized_refuses(&oversized_in_memory);
}

#[test]
fn the_merkle_root_is_deterministic_and_tamper_evident() {
    let build = || {
        EvidencePackage::new(prov())
            .with_claim(verified("c1"))
            .with_claim(estimated("c2"))
    };
    // identical packages -> identical content address.
    assert_eq!(package_root(&build()), package_root(&build()));
    assert_eq!(
        package_root(&build()).to_hex(),
        "1a917f759d541819f863a02787aceb0408cacae9bf55c8a50230d8ca9db89465",
        "schema-v6 package-root fixture (re-pinned for the v6 domain separation)"
    );
    // tampering with a claim changes the root.
    let tampered = EvidencePackage::new(prov())
        .with_claim(Claim::from_certificate(
            "c1",
            "TAMPERED",
            -1.0,
            1.0,
            "test-solver/cert",
            CANONICAL_DATASET_HASH,
        ))
        .with_claim(estimated("c2"));
    assert_ne!(package_root(&build()), package_root(&tampered));
}

#[test]
fn the_merkle_root_covers_reproducibility_provenance() {
    let pkg = EvidencePackage::new(prov()).with_claim(verified("c1"));
    let changed_code = EvidencePackage::new(Provenance::new("commit-other", "lock-deadbeef"))
        .with_claim(verified("c1"));
    let changed_lock = EvidencePackage::new(Provenance::new("commit-abc123", "lock-other"))
        .with_claim(verified("c1"));

    assert_ne!(package_root(&pkg), package_root(&changed_code));
    assert_ne!(package_root(&pkg), package_root(&changed_lock));
}

#[test]
fn an_unsupported_format_version_is_rejected() {
    let mut pkg = EvidencePackage::new(prov()).with_claim(estimated("e1"));
    pkg.format_version = 999;
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::UnsupportedFormat { found: 999 })
    ));
}

#[test]
fn a_signature_is_optional_and_detached() {
    let unsigned = EvidencePackage::new(prov()).with_claim(estimated("e1"));
    assert!(
        unsigned.verify().is_ok(),
        "no signature is fine (content-addressed)"
    );
    assert!(package_json(&unsigned).contains("\"signature\":null"));
    let signed = unsigned.clone().signed("ed25519:deadbeef");
    // signing does not change the content address (detached).
    assert_eq!(package_root(&unsigned), package_root(&signed));
    assert!(package_json(&signed).contains("ed25519:deadbeef"));
}

struct ExactSignatureVerifier;
struct PermissiveSignatureVerifier;

impl SignatureVerifier for ExactSignatureVerifier {
    fn verify(&self, request: &SignatureRequest<'_>) -> VerificationDecision {
        VerificationDecision {
            accepted: request.signature == format!("test-signature:{}", request.subject_hash())
                && match request.purpose {
                    SignaturePurpose::PackageRootAttestation => true,
                    SignaturePurpose::ReleaseApproval {
                        checker_protocol: 4,
                        expected_root,
                        admission_context,
                    } => {
                        expected_root == request.package_root
                            && admission_context != ContentHash([0; 32])
                    }
                    SignaturePurpose::ReleaseApproval { .. } => false,
                },
            policy_fingerprint: fixture_policy("exact-signature-verifier"),
        }
    }
}

impl SignatureVerifier for PermissiveSignatureVerifier {
    fn verify(&self, _request: &SignatureRequest<'_>) -> VerificationDecision {
        VerificationDecision::accept(fixture_policy("permissive-signature-verifier"))
    }
}

#[test]
fn signature_subject_hash_separates_every_purpose_axis() {
    let root = ContentHash([0x11; 32]);
    let other_root = ContentHash([0x22; 32]);
    let context = ContentHash([0x33; 32]);
    let other_context = ContentHash([0x44; 32]);
    let release =
        |checker_protocol, expected_root, admission_context| SignaturePurpose::ReleaseApproval {
            checker_protocol,
            expected_root,
            admission_context,
        };

    let subjects = [
        signature_subject_hash(root, SignaturePurpose::PackageRootAttestation),
        signature_subject_hash(other_root, SignaturePurpose::PackageRootAttestation),
        signature_subject_hash(root, release(4, root, context)),
        signature_subject_hash(root, release(3, root, context)),
        signature_subject_hash(root, release(4, other_root, context)),
        signature_subject_hash(root, release(4, root, other_context)),
    ];
    for (left_index, left) in subjects.iter().enumerate() {
        for right in &subjects[left_index + 1..] {
            assert_ne!(left, right, "distinct signature axes must not alias");
        }
    }
    assert_eq!(
        subjects[2],
        signature_subject_hash(root, release(4, root, context)),
        "the exact typed subject must replay deterministically"
    );
}

#[test]
fn signature_coverage_requires_authentication_and_rejection_fails_closed() {
    use fs_crosswalk::PackageConcept;

    let unsigned = EvidencePackage::new(prov()).with_claim(estimated("signed-evidence"));
    let root = package_root(&unsigned);
    let integrity_purpose = SignaturePurpose::PackageRootAttestation;
    let signed = unsigned.clone().signed(format!(
        "test-signature:{}",
        signature_subject_hash(root, integrity_purpose)
    ));
    let raw_presence = fs_package::package_presence(&signed);
    assert!(
        !raw_presence
            .iter()
            .find(|row| row.concept() == PackageConcept::Signature)
            .expect("signature concept")
            .present()
    );

    let capabilities =
        VerificationCapabilities::deny_all().with_signatures(&ExactSignatureVerifier);
    let report = signed
        .verify_with(&capabilities)
        .expect("exact signature authenticates");
    assert!(matches!(
        report.receipt().signature(),
        SignatureStatus::Authenticated(authenticated)
            if authenticated.purpose() == SignaturePurpose::PackageRootAttestation
    ));
    assert_eq!(
        report.receipt().policy_fingerprints().signatures(),
        Some(fixture_policy("exact-signature-verifier"))
    );
    let integrity_presence = fs_package::package_presence_with(&signed, &capabilities);
    assert!(integrity_presence.validate_decision_hash());
    assert!(
        !integrity_presence
            .iter()
            .find(|row| row.concept() == PackageConcept::Signature)
            .expect("signature concept")
            .present(),
        "a generic integrity attestation is not regulatory release approval"
    );

    let release_capabilities = VerificationCapabilities::deny_all().with_release_signatures(
        &ExactSignatureVerifier,
        4,
        root,
    );
    let unsigned_report = unsigned
        .verify()
        .expect("unsigned release subject verifies");
    let release_purpose = SignaturePurpose::ReleaseApproval {
        checker_protocol: 4,
        expected_root: root,
        admission_context: unsigned_report.receipt().release_admission_context(),
    };
    let release_signed = unsigned.signed(format!(
        "test-signature:{}",
        signature_subject_hash(root, release_purpose)
    ));
    assert!(
        fs_package::package_presence_with(&release_signed, &release_capabilities)
            .iter()
            .find(|row| row.concept() == PackageConcept::Signature)
            .expect("signature concept")
            .present(),
        "an authenticated, purpose-bound release approval establishes signature coverage"
    );

    let forged = EvidencePackage::new(prov())
        .with_claim(estimated("signed-evidence"))
        .signed("test-signature:wrong-root");
    assert!(matches!(
        forged.verify_with(&capabilities),
        Err(PackageError::SignatureRefused {
            why: "rejected by the injected verifier",
            ..
        })
    ));
    assert!(
        fs_package::package_presence_with(&forged, &capabilities)
            .iter()
            .all(|row| !row.present()),
        "an invalid supplied signature suppresses all coverage"
    );
}

#[test]
fn release_signature_purpose_must_name_the_recomputed_package_root() {
    let unsigned = EvidencePackage::new(prov()).with_claim(estimated("release-subject"));
    let root = package_root(&unsigned);
    let signed = unsigned.signed(format!("test-signature:{root}"));
    let mut wrong_root_bytes = *root.as_bytes();
    wrong_root_bytes[0] ^= 1;
    let capabilities = VerificationCapabilities::deny_all().with_release_signatures(
        &PermissiveSignatureVerifier,
        4,
        ContentHash(wrong_root_bytes),
    );
    assert!(matches!(
        signed.verify_with(&capabilities),
        Err(PackageError::SignatureRefused {
            why: "release-approval purpose names a different package root",
            policy_fingerprint: None,
        })
    ));
    assert!(
        fs_package::package_presence_with(&signed, &capabilities)
            .iter()
            .all(|row| !row.present()),
        "a permissive verifier cannot authenticate release approval for another root"
    );
}

#[test]
fn all_estimated_packages_do_not_claim_certificate_coverage() {
    use fs_crosswalk::PackageConcept;

    let pkg = EvidencePackage::new(prov()).with_claim(estimated("estimate-only"));
    let presence = fs_package::package_presence(&pkg);
    assert!(
        presence
            .iter()
            .find(|row| row.concept() == PackageConcept::EstimatedColor)
            .expect("estimated concept")
            .present()
    );
    assert!(
        !presence
            .iter()
            .find(|row| row.concept() == PackageConcept::Certificate)
            .expect("certificate concept")
            .present()
    );
    for concept in [PackageConcept::ClaimOrigin, PackageConcept::Provenance] {
        let row = presence
            .iter()
            .find(|row| row.concept() == concept)
            .expect("concept is always judged");
        assert!(!row.present(), "{}", row.why());
    }
}

struct ExactFalsifierArtifactVerifier {
    package_root: ContentHash,
    claim_subject_hash: ContentHash,
}

impl FalsifierVerifier for ExactFalsifierArtifactVerifier {
    fn verify(&self, request: &FalsifierRequest<'_>) -> VerificationDecision {
        let accepted = request.package_provenance == &prov()
            && request.package_root == self.package_root
            && request.claim_index == 0
            && request.claim_id == "falsified"
            && request.statement == "bounded estimate"
            && matches!(
                request.color,
                Color::Estimated {
                    estimator,
                    dispersion,
                } if estimator == "probe" && dispersion.to_bits() == 1.0f64.to_bits()
            )
            && matches!(request.origin, fs_package::ClaimOrigin::EstimatedSource { estimator } if estimator == "probe")
            && request.claim_subject_hash == self.claim_subject_hash
            && request.falsifier_index == 0
            && request.name == "boundary-probe"
            && request.attempts == 32
            && !request.refuted
            && request.detail == "no counterexample"
            && request.artifact_hash.to_hex() == CANONICAL_DATASET_HASH;
        VerificationDecision {
            accepted,
            policy_fingerprint: fixture_policy("exact-falsifier-artifact-verifier"),
        }
    }
}

fn falsifier_subject(detail: &str, artifact_hash: &str) -> EvidencePackage {
    EvidencePackage::new(prov()).with_claim(
        Claim::estimated("falsified", "bounded estimate", "probe", 1.0).with_falsifier(
            FalsifierRecord {
                name: "boundary-probe".to_string(),
                attempts: 32,
                refuted: false,
                detail: detail.to_string(),
                artifact_hash: artifact_hash.to_string(),
            },
        ),
    )
}

#[test]
fn falsifier_artifacts_bind_the_exact_claim_and_record() {
    use fs_crosswalk::PackageConcept;

    let pkg = falsifier_subject("no counterexample", CANONICAL_DATASET_HASH);
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::FalsifierRefused {
            why: "falsifier capability missing",
            ..
        })
    ));
    let verifier = ExactFalsifierArtifactVerifier {
        package_root: package_root(&pkg),
        claim_subject_hash: pkg.declared_claims_unverified()[0]
            .declared_verification_subject_hash_unverified(),
    };
    let capabilities = VerificationCapabilities::deny_all().with_falsifiers(&verifier);
    let report = pkg
        .verify_with(&capabilities)
        .expect("exact falsifier subject authenticates");
    assert_eq!(
        report.receipt().policy_fingerprints().falsifiers(),
        Some(fixture_policy("exact-falsifier-artifact-verifier"))
    );
    let coverage = fs_package::package_presence_with(&pkg, &capabilities);
    assert!(coverage.receipt().is_some());
    assert!(
        coverage
            .iter()
            .find(|row| row.concept() == PackageConcept::FalsifierLog)
            .expect("falsifier concept")
            .present()
    );

    for forged in [
        falsifier_subject("different outcome", CANONICAL_DATASET_HASH),
        falsifier_subject(
            "no counterexample",
            "1123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
        ),
    ] {
        assert!(matches!(
            forged.verify_with(&capabilities),
            Err(PackageError::FalsifierRefused {
                why: "rejected by the injected verifier",
                ..
            })
        ));
    }
}

#[test]
fn json_is_deterministic_and_carries_the_root() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(validated("c2", good_regime(), "wt-2026"));
    let j = package_json(&pkg);
    assert_eq!(j, package_json(&pkg));
    assert!(j.starts_with('{') && j.ends_with('}'));
    assert!(j.contains(&package_root(&pkg).to_hex()));
    assert!(j.contains("\"format_version\":6"), "schema v6");
    // v3 carries COMPLETE payloads, not just rank labels.
    assert!(j.contains("\"lo_bits\":") && j.contains("\"dataset\":"));
}

#[test]
fn schema_v6_magnitude_shape_is_closed_and_waiver_aware() {
    let pkg = EvidencePackage::new(prov()).with_claim(estimated("shape"));
    let json = package_json(&pkg);
    assert!(json.contains("\"waived_unquantified\":0"));
    assert_eq!(EvidencePackage::from_json(&json).unwrap(), pkg);

    let missing_v6_field = json.replace(",\"waived_unquantified\":0", "");
    let error = EvidencePackage::from_json(&missing_v6_field)
        .expect_err("a v5-shaped magnitude budget must not parse as v6");
    assert!(error.why.contains("waived_unquantified"), "{error}");

    let stale_v5 = missing_v6_field.replacen("\"format_version\":6", "\"format_version\":5", 1);
    let error = EvidencePackage::from_json(&stale_v5).expect_err("schema v5 is not v6");
    assert!(error.why.contains("unsupported version 5"), "{error}");
}

/// qmao.6.1 — the schema-v3 round trip and its fail-closed walls: a
/// package decode-encodes stably (golden), floats travel bit-exactly,
/// and hostile/missing/unknown/non-finite/forged/tampered inputs each
/// refuse with a structured parse error.
#[test]
fn v3_round_trip_and_fail_closed_walls() {
    let pkg = EvidencePackage::new(Provenance::new("frankensim@abc123", "lock:deadbeef"))
        .with_claim(Claim::from_certificate(
            "c-verified",
            "tip deflection within bound",
            0.1875,
            0.25,
            "test-solver/cert",
            CANONICAL_DATASET_HASH,
        ))
        .with_claim(Claim::anchored(
            "c-validated",
            "k-epsilon matched within regime",
            fs_evidence::ValidityDomain::unconstrained().with("reynolds", 1e3, 1e5),
            "tunnel-run-9",
            CANONICAL_DATASET_HASH,
        ))
        .with_claim(Claim::estimated(
            "c-estimated",
            "surrogate prediction",
            "pod-deim",
            0.02,
        ))
        .signed("test-key/1234abcd");
    // Golden decode-encode stability: parse(to_json) == pkg, and the
    // re-emission is byte-identical (semantic AND textual round trip).
    let json = package_json(&pkg);
    let back = EvidencePackage::from_json(&json).expect("canonical JSON parses");
    assert_eq!(back, pkg, "semantic round trip");
    assert_eq!(package_json(&back), json, "textual round trip");
    let leading_zero = json.replacen("\"format_version\":6", "\"format_version\":06", 1);
    assert!(
        EvidencePackage::from_json(&leading_zero).is_err(),
        "non-JSON leading-zero integer refuses"
    );
    // Forged Verified claim: structurally shaped, but the embedded root
    // no longer recomputes from the tampered content — refused at parse.
    let forged = json.replace("tip deflection within bound", "tip deflection PROVEN SAFE");
    let err = EvidencePackage::from_json(&forged).expect_err("forgery refused");
    assert!(err.why.contains("does not recompute"), "{err}");
    // Widening the certified interval (payload tamper) also refused.
    let widened = json.replace(
        &format!("{:016x}", 0.25f64.to_bits()),
        &format!("{:016x}", 2.5f64.to_bits()),
    );
    assert!(
        EvidencePackage::from_json(&widened).is_err(),
        "payload tamper refused"
    );
    // Non-finite certificate bits: fail closed.
    let nan = json.replace(
        &format!("{:016x}", 0.1875f64.to_bits()),
        &format!("{:016x}", f64::NAN.to_bits()),
    );
    let err = EvidencePackage::from_json(&nan).expect_err("NaN certificate refused");
    assert!(err.why.contains("non-finite"), "{err}");
    // Unknown fields: closed schema.
    let unknown = json.replacen(
        "{\"format_version\"",
        "{\"vendor_extra\":1,\"format_version\"",
        1,
    );
    let err = EvidencePackage::from_json(&unknown).expect_err("unknown field refused");
    assert!(err.why.contains("unknown field"), "{err}");
    // Missing fields: fail closed.
    let missing = json.replacen("\"signature\":\"test-key/1234abcd\",", "", 1);
    assert!(
        EvidencePackage::from_json(&missing).is_err(),
        "missing field refused"
    );
    // Unknown color kind: fail closed.
    let bad_kind = json.replacen("\"kind\":\"verified\"", "\"kind\":\"blessed\"", 1);
    let err = EvidencePackage::from_json(&bad_kind).expect_err("unknown kind refused");
    assert!(err.why.contains("unknown color kind"), "{err}");
    // Drifted magnitude budget: refused (it must re-derive).
    let mb_tag = "\"validated_unquantified\":1";
    assert!(json.contains(mb_tag), "fixture sanity");
    let drifted = json.replacen(mb_tag, "\"validated_unquantified\":7", 1);
    let err = EvidencePackage::from_json(&drifted).expect_err("budget drift refused");
    assert!(err.why.contains("re-derive"), "{err}");
    // Hostile garbage and truncation: structured refusals, no panic.
    assert!(EvidencePackage::from_json("{\"format_version\":2").is_err());
    assert!(EvidencePackage::from_json(&json[..json.len() / 2]).is_err());
    assert!(EvidencePackage::from_json("").is_err());
    // Canonical fixed-width roots and JSON's control-character rule are
    // enforced at the parser boundary, before semantic verification.
    let root = package_root(&pkg).to_hex();
    let short_root = json.replacen(&format!("\"{root}\""), &format!("\"{}\"", &root[1..]), 1);
    let err = EvidencePackage::from_json(&short_root).expect_err("short root refused");
    assert!(err.why.contains("64 lowercase hex chars"), "{err}");
    let raw_control = json.replacen("surrogate prediction", "surrogate\nprediction", 1);
    let err = EvidencePackage::from_json(&raw_control).expect_err("raw control refused");
    assert!(err.why.contains("control character"), "{err}");
}

/// qmao.6.1 — the magnitude budget attributes ERROR MAGNITUDES, not
/// claim counts, and reconciles with an independent recomputation.
#[test]
fn magnitude_budget_reconciles() {
    let provenance = Provenance::new("v", "l");
    let pkg = EvidencePackage::new(provenance.clone())
        .with_claim(source_claim(&provenance, 0, "a", "s", 0.0, 0.5))
        .with_claim(source_claim(&provenance, 1, "b", "s", 1.0, 1.25))
        .with_claim(Claim::estimated("c", "s", "e", 0.125))
        .with_claim(Claim::anchored(
            "d",
            "s",
            fs_evidence::ValidityDomain::unconstrained().with("re", 1.0, 2.0),
            "ds",
            CANONICAL_DATASET_HASH,
        ));
    let mb = pkg
        .magnitude_budget_with(&source_capabilities())
        .expect("scientific magnitude requires admitted origins");
    assert_eq!(mb.verified_width.to_bits(), 0.75f64.to_bits());
    assert_eq!(mb.estimated_dispersion.to_bits(), 0.125f64.to_bits());
    assert_eq!(
        mb.validated_unquantified, 1,
        "regional trust counted, never numerified"
    );
    assert_eq!(
        mb.quantified_total.to_bits(),
        (mb.verified_width + mb.estimated_dispersion).to_bits(),
        "total reconciles with its parts"
    );
}

/// qmao.6.1 — crosswalk coverage derives from fields ACTUALLY PRESENT
/// in a parsed package: a mapped concept with absent evidence is
/// "mapped but absent", never covered; no static mapping can report
/// evidence that the package does not carry.
#[test]
fn coverage_cannot_claim_absent_evidence() {
    use fs_crosswalk::{PackageConcept, Standard};
    use fs_package::{CoverageStatus, package_coverage_with, package_presence_with};
    // Unsigned, verified-only package: no validated claims, no regime
    // tags, no datasets, no signature.
    let provenance = Provenance::new("v", "lock");
    let pkg = EvidencePackage::new(provenance.clone()).with_claim(source_claim(
        &provenance,
        0,
        "c",
        "bounded",
        0.0,
        1.0,
    ));
    let capabilities = source_capabilities();
    let presence = package_presence_with(&pkg, &capabilities);
    assert!(presence.validate_decision_hash());
    assert!(presence.receipt().is_some());
    let get = |c: PackageConcept| presence.iter().find(|p| p.concept() == c).unwrap();
    assert!(get(PackageConcept::VerifiedColor).present());
    assert!(!get(PackageConcept::ValidatedColor).present());
    assert!(!get(PackageConcept::RegimeTag).present());
    assert!(!get(PackageConcept::AnchoringDataset).present());
    assert!(!get(PackageConcept::Signature).present());
    assert!(
        !get(PackageConcept::FalsifierLog).present(),
        "absent falsifier records can never read as present"
    );
    for standard in Standard::ALL {
        let coverage = package_coverage_with(&pkg, standard, &capabilities);
        assert!(coverage.validate_decision_hash());
        assert_eq!(coverage.standard(), standard);
        assert_eq!(
            coverage.crosswalk_version(),
            fs_crosswalk::CROSSWALK_VERSION
        );
        assert_eq!(coverage.package_format(), fs_package::FORMAT_VERSION);
        assert_eq!(
            coverage
                .receipt()
                .map(fs_package::VerificationReceipt::receipt_hash),
            presence
                .receipt()
                .map(fs_package::VerificationReceipt::receipt_hash)
        );
        assert!(
            coverage
                .iter()
                .all(|(_, _, why)| why.contains("package evidence:"))
        );
        for (concept, status, _why) in &coverage {
            if matches!(status, CoverageStatus::Covered) {
                let p = presence.iter().find(|p| p.concept() == *concept).unwrap();
                assert!(
                    p.present(),
                    "{concept:?} covered without evidence for {standard:?}"
                );
            }
            if *concept == PackageConcept::FalsifierLog {
                assert!(
                    !matches!(status, CoverageStatus::Covered),
                    "falsifier logs without records cannot be covered ({standard:?})"
                );
            }
        }
    }
}

#[test]
fn coverage_requires_a_valid_package_and_authenticated_evidence() {
    use fs_crosswalk::{PackageConcept, Standard};
    use fs_package::{CoverageStatus, package_coverage, package_presence};

    let refuted = EvidencePackage::new(prov())
        .with_claim(verified("refuted").with_falsifier(FalsifierRecord {
            name: "interval-probe".to_string(),
            attempts: 1,
            refuted: true,
            detail: "counterexample found".to_string(),
            artifact_hash: CANONICAL_DATASET_HASH.to_string(),
        }))
        .signed("raw-detached-signature");
    assert!(package_presence(&refuted).iter().all(|row| !row.present()));
    for standard in Standard::ALL {
        assert!(
            package_coverage(&refuted, standard)
                .iter()
                .all(|(_, status, _)| !matches!(status, CoverageStatus::Covered)),
            "invalid package covered a concept for {standard:?}"
        );
    }

    let merely_signed = EvidencePackage::new(prov())
        .with_claim(verified("signed"))
        .signed("unauthenticated-bytes");
    let presence = package_presence(&merely_signed);
    let signature = presence
        .iter()
        .find(|row| row.concept() == PackageConcept::Signature)
        .expect("signature concept judged");
    assert!(!signature.present(), "{}", signature.why());
    for standard in Standard::ALL {
        let coverage = package_coverage(&merely_signed, standard);
        let status = coverage
            .iter()
            .find(|(concept, _, _)| *concept == PackageConcept::Signature)
            .expect("signature is mapped for every standard")
            .1
            .clone();
        assert!(!matches!(status, CoverageStatus::Covered));
    }
}

#[test]
fn dataset_coverage_requires_a_matching_valid_anchor() {
    use fs_crosswalk::PackageConcept;
    use fs_package::package_presence_with;

    // Schema v6: a validated claim can only be CONSTRUCTED anchored to
    // its own dataset (the AnchoredSource origin attaches the matching
    // record), so "validated but unanchored" is unrepresentable in
    // memory; an EXTRA unrelated anchor neither helps nor hides it.
    let unrelated = EvidencePackage::new(prov()).with_claim(
        validated("v", good_regime(), "wind-tunnel-2026")
            .with_anchor("different-dataset", CANONICAL_DATASET_HASH),
    );
    let capabilities = source_capabilities();
    unrelated
        .verify_with(&capabilities)
        .expect("an extra unrelated anchor is structurally valid");
    let unrelated_presence = package_presence_with(&unrelated, &capabilities);
    let anchor = unrelated_presence
        .iter()
        .find(|row| row.concept() == PackageConcept::AnchoringDataset)
        .expect("anchor concept judged");
    assert!(anchor.present(), "{}", anchor.why());

    let matching = EvidencePackage::new(prov()).with_claim(
        validated("v", good_regime(), "wind-tunnel-2026")
            .with_anchor("wind-tunnel-2026", CANONICAL_DATASET_HASH),
    );
    matching
        .verify_with(&capabilities)
        .expect("matching anchor verifies");
    let matching_presence = package_presence_with(&matching, &capabilities);
    let anchor = matching_presence
        .iter()
        .find(|row| row.concept() == PackageConcept::AnchoringDataset)
        .expect("anchor concept judged");
    assert!(anchor.present(), "{}", anchor.why());
}

/// xfxq (schema v3) — composition receipts re-run solver-free, refuted
/// falsifiers fail, anchors and falsifier logs travel with claims, and
/// every new field round-trips through the strict parser bound into the
/// content address.
#[test]
fn v3_receipts_falsifiers_anchors() {
    use fs_package::PackageError;
    let ve = |lo: f64, hi: f64| Color::Verified { lo, hi };
    // A well-formed derived package: c = a + b with a valid receipt.
    let good = derived_package();
    verify_package(&good).expect("receipt re-derives");
    // Round trip: the v3 fields survive the strict parser bit-for-bit.
    let back = EvidencePackage::from_json(&package_json(&good)).expect("v3 parses");
    assert_eq!(back, good);
    assert_eq!(
        back.declared_claims_unverified()[2].declared_falsifiers_unverified()[0].attempts,
        512
    );
    assert_eq!(
        back.declared_claims_unverified()[2].declared_anchors_unverified()[0].dataset_id,
        "wt-2026-run9"
    );
    // FORGED receipt: claiming Verified while a parent is Estimated —
    // the re-run composition cannot reproduce it (semantic catch, not
    // just the content-address catch).
    let forged = EvidencePackage::new(Provenance::new("v", "l"))
        .with_claim(Claim::estimated("a", "shaky", "guess", 0.5))
        .with_claim(Claim::from_certificate(
            "b",
            "solid",
            1.0,
            2.0,
            "test-solver/cert",
            CANONICAL_DATASET_HASH,
        ))
        .with_claim(Claim::derived(
            "c",
            "laundered",
            ve(2.0, 4.0),
            vec![0, 1],
            IntervalOp::Add,
            CANONICAL_DATASET_HASH,
        ));
    assert!(matches!(
        forged.verify(),
        Err(PackageError::ReceiptMismatch { claim }) if claim == "c"
    ));
    // Forward/self parent references refuse.
    let cyclic = EvidencePackage::new(Provenance::new("v", "l")).with_claim(Claim::derived(
        "a",
        "s",
        ve(0.0, 1.0),
        vec![0],
        IntervalOp::Hull,
        CANONICAL_DATASET_HASH,
    ));
    assert!(matches!(
        cyclic.verify(),
        Err(PackageError::BadReceiptParent { parent: 0, .. })
    ));
    // A refuted falsifier fails the whole claim.
    let refuted = EvidencePackage::new(Provenance::new("v", "l")).with_claim(
        Claim::from_certificate(
            "a",
            "wrong",
            0.0,
            1.0,
            "test-solver/cert",
            CANONICAL_DATASET_HASH,
        )
        .with_falsifier(FalsifierRecord {
            name: "adversary".to_string(),
            attempts: 3,
            refuted: true,
            detail: "counterexample at x=0.7".to_string(),
            artifact_hash: CANONICAL_DATASET_HASH.to_string(),
        }),
    );
    assert!(matches!(
        refuted.verify(),
        Err(PackageError::RefutedClaim { falsifier, .. }) if falsifier == "adversary"
    ));
}

struct ExactDerivationArtifactVerifier {
    package_root: ContentHash,
    child_subject_hash: ContentHash,
    parent_claim_hash: ContentHash,
}

impl DerivationVerifier for ExactDerivationArtifactVerifier {
    fn verify(&self, request: &DerivationRequest<'_>) -> VerificationDecision {
        let accepted = request.package_provenance == &prov()
            && request.package_root == self.package_root
            && request.claim_index == 1
            && request.claim_id == "derived"
            && request.statement == "same bounded estimate"
            && matches!(
                request.color,
                Color::Estimated {
                    estimator,
                    dispersion,
                } if estimator == "probe" && dispersion.to_bits() == 1.0f64.to_bits()
            )
            && request.child_subject_hash == self.child_subject_hash
            && request.anchors.is_empty()
            && request.op == IntervalOp::Hull
            && request.parent_indices == [0]
            && request.parent_claim_hashes == [self.parent_claim_hash]
            && request.artifact_hash.to_hex() == CANONICAL_DATASET_HASH;
        VerificationDecision {
            accepted,
            policy_fingerprint: fixture_policy("exact-derivation-artifact-verifier"),
        }
    }
}

fn derivation_subject(statement: &str, artifact_hash: &str) -> EvidencePackage {
    let color = Color::Estimated {
        estimator: "probe".to_string(),
        dispersion: 1.0,
    };
    EvidencePackage::new(prov())
        .with_claim(Claim::estimated("parent", "bounded estimate", "probe", 1.0))
        .with_claim(Claim::derived(
            "derived",
            statement,
            color,
            vec![0],
            IntervalOp::Hull,
            artifact_hash,
        ))
}

#[test]
fn derivation_artifacts_bind_child_statement_color_and_ordered_parents() {
    let pkg = derivation_subject("same bounded estimate", CANONICAL_DATASET_HASH);
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::DerivationRefused {
            why: "derivation capability missing",
            ..
        })
    ));
    let claims = pkg.declared_claims_unverified();
    let verifier = ExactDerivationArtifactVerifier {
        package_root: package_root(&pkg),
        child_subject_hash: claims[1].declared_verification_subject_hash_unverified(),
        parent_claim_hash: claims[0].declared_content_hash_unverified(),
    };
    let capabilities = VerificationCapabilities::deny_all().with_derivations(&verifier);
    let report = pkg
        .verify_with(&capabilities)
        .expect("exact derivation subject authenticates");
    assert_eq!(
        report.receipt().policy_fingerprints().derivations(),
        Some(fixture_policy("exact-derivation-artifact-verifier"))
    );
    assert!(report.receipt().validate_hash());

    for forged in [
        derivation_subject("unrelated safety assertion", CANONICAL_DATASET_HASH),
        derivation_subject(
            "same bounded estimate",
            "1123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef",
        ),
    ] {
        assert!(matches!(
            forged.verify_with(&capabilities),
            Err(PackageError::DerivationRefused {
                why: "rejected by the injected verifier",
                ..
            })
        ));
    }
}

#[test]
fn derived_receipts_compare_signed_zero_by_exact_bits() {
    let estimated_parent = Color::Estimated {
        estimator: "signed-zero-probe".to_string(),
        dispersion: -0.0,
    };
    let verified_parent = Color::Verified { lo: 0.0, hi: 0.0 };
    let recomputed = fs_evidence::compose(&estimated_parent, &verified_parent, IntervalOp::Add);
    assert!(matches!(
        recomputed,
        Color::Estimated { dispersion, .. } if dispersion.to_bits() == 0.0f64.to_bits()
    ));

    let package = EvidencePackage::new(prov())
        .with_claim(Claim::estimated(
            "estimated-parent",
            "signed zero estimate",
            "signed-zero-probe",
            -0.0,
        ))
        .with_claim(verified("verified-parent"))
        .with_claim(Claim::derived(
            "derived",
            "bit-exact derived identity",
            Color::Estimated {
                estimator: "signed-zero-probe".to_string(),
                dispersion: -0.0,
            },
            vec![0, 1],
            IntervalOp::Add,
            CANONICAL_DATASET_HASH,
        ));
    assert!(matches!(
        package.verify(),
        Err(PackageError::ReceiptMismatch { claim }) if claim == "derived"
    ));
}

struct DivergentDerivationPolicy;

impl DerivationVerifier for DivergentDerivationPolicy {
    fn verify(&self, request: &DerivationRequest<'_>) -> VerificationDecision {
        VerificationDecision::accept(ContentHash([request.claim_index as u8; 32]))
    }
}

#[test]
fn one_capability_cannot_rotate_policy_fingerprints_mid_package() {
    let color = Color::Estimated {
        estimator: "probe".to_string(),
        dispersion: 1.0,
    };
    let pkg = EvidencePackage::new(prov())
        .with_claim(Claim::estimated("parent", "bounded", "probe", 1.0))
        .with_claim(Claim::derived(
            "d1",
            "first derivation",
            color.clone(),
            vec![0],
            IntervalOp::Hull,
            CANONICAL_DATASET_HASH,
        ))
        .with_claim(Claim::derived(
            "d2",
            "second derivation",
            color,
            vec![1],
            IntervalOp::Hull,
            CANONICAL_DATASET_HASH,
        ));
    let capabilities =
        VerificationCapabilities::deny_all().with_derivations(&DivergentDerivationPolicy);
    assert!(matches!(
        pkg.verify_with(&capabilities),
        Err(PackageError::PolicyFingerprintRefused {
            capability: "derivations",
            why: "policy fingerprint changed during package verification",
            previous,
            observed,
        }) if previous == ContentHash([1; 32]) && observed == ContentHash([2; 32])
    ));

    let unused = EvidencePackage::new(prov()).with_claim(estimated("unused-policy"));
    let report = unused
        .verify_with(&capabilities)
        .expect("unused capability changes no proof state");
    assert_eq!(report.receipt().policy_fingerprints().derivations(), None);
}

#[test]
fn v3_falsifier_coverage_and_content_binding() {
    let good = derived_package();
    let capabilities = source_capabilities();
    let presence = fs_package::package_presence_with(&good, &capabilities);
    let falsifier = presence
        .iter()
        .find(|row| row.concept() == fs_crosswalk::PackageConcept::FalsifierLog)
        .expect("falsifier concept judged");
    assert!(falsifier.present(), "{}", falsifier.why());

    let tampered = package_json(&good).replace("\"refuted\":false", "\"refuted\":true");
    assert!(
        EvidencePackage::from_json(&tampered).is_err(),
        "content root binds falsifier fields"
    );
}

/// 7uq9 (schema v4) — the content address is a domain-separated 32-byte
/// BLAKE3 root; legacy v3 transports and legacy 16-hex FNV roots are
/// refused with messages that NAME the incompatibility.
#[test]
fn v4_blake3_root_refuses_legacy_transports() {
    let pkg = EvidencePackage::new(prov()).with_claim(verified("c1"));
    let json = package_json(&pkg);
    let root = package_root(&pkg).to_hex();
    assert_eq!(root.len(), 64, "32-byte BLAKE3 root renders as 64 hex");
    assert!(json.contains(&format!("\"merkle_root\":\"{root}\"")));

    // A v3 transport is refused BY VERSION before any field is read.
    let v4 = json.replacen("\"format_version\":6", "\"format_version\":4", 1);
    let err = EvidencePackage::from_json(&v4).expect_err("v4 refused");
    assert!(err.why.contains("unsupported version 4"), "{err}");

    // A legacy 16-hex FNV root inside a v4 envelope is named as such.
    let legacy = json.replacen(&format!("\"{root}\""), "\"deadbeefcafe0123\"", 1);
    let err = EvidencePackage::from_json(&legacy).expect_err("legacy root refused");
    assert!(err.why.contains("legacy v3 FNV"), "{err}");

    // Right length, non-hex content: refused at the boundary.
    let garbled = json.replacen(
        &format!("\"{root}\""),
        &format!("\"{}\"", "z".repeat(64)),
        1,
    );
    let err = EvidencePackage::from_json(&garbled).expect_err("non-hex root refused");
    assert!(err.why.contains("lowercase hexadecimal"), "{err}");

    // A validly-formatted but WRONG root: recomputation refuses it.
    let mut flipped = root.clone().into_bytes();
    flipped[0] = if flipped[0] == b'0' { b'1' } else { b'0' };
    let wrong = json.replacen(&root, core::str::from_utf8(&flipped).unwrap(), 1);
    let err = EvidencePackage::from_json(&wrong).expect_err("wrong root refused");
    assert!(err.why.contains("does not recompute"), "{err}");
}

/// 7uq9 — domain separation: the root is not the bare BLAKE3 of any
/// undomained serialization. The package header binds claim count, so
/// zero/one and duplicate-tail shapes have distinct roots; the BLAKE3
/// owner separately locks cross-mode and cross-domain behavior.
#[test]
fn v4_root_is_domain_separated() {
    let pkg = EvidencePackage::new(prov()).with_claim(verified("c1"));
    let root = package_root(&pkg);
    // Not the undomained hash of the canonical JSON or of the claim.
    assert_ne!(root, fs_blake3::hash_bytes(package_json(&pkg).as_bytes()));
    assert_ne!(
        root,
        package_root(&EvidencePackage::new(prov())),
        "header claim count distinguishes zero and one claim"
    );
    assert_ne!(
        root,
        package_root(
            &EvidencePackage::new(prov())
                .with_claim(verified("c1"))
                .with_claim(verified("c1"))
        ),
        "duplicating the final claim changes the tree identity"
    );
    // Stable across recomputation and bound to claim ORDER.
    let swapped = EvidencePackage::new(prov())
        .with_claim(verified("c2"))
        .with_claim(verified("c1"));
    let ordered = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(verified("c2"));
    assert_ne!(package_root(&swapped), package_root(&ordered));
}

struct ExactSourceVerifier;

impl SourceCertificateVerifier for ExactSourceVerifier {
    fn verify(&self, request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        let accepted = request.package_provenance == &prov()
            && request.claim_index == 0
            && request.claim_id == "source"
            && request.statement == "certified interval"
            && request.lo.to_bits() == 1.0f64.to_bits()
            && request.hi.to_bits() == 2.0f64.to_bits()
            && request.producer == "test-solver/cert"
            && request.certificate_hash.to_hex() == CANONICAL_DATASET_HASH;
        VerificationDecision {
            accepted,
            policy_fingerprint: fixture_policy("exact-source-verifier"),
        }
    }
}

static EXACT_SOURCE_VERIFIER: ExactSourceVerifier = ExactSourceVerifier;

struct ExactAnchorSubjectVerifier;

impl AnchoredSourceVerifier for ExactAnchorSubjectVerifier {
    fn verify(&self, request: &AnchoredSourceRequest<'_>) -> VerificationDecision {
        let accepted = request.package_provenance == &prov()
            && request.claim_index == 0
            && request.claim_id == "anchor"
            && request.statement == "matches exact reference data"
            && request.regime == &good_regime()
            && request.dataset_id == "wind-tunnel-2026"
            && request.content_hash.to_hex() == CANONICAL_DATASET_HASH;
        VerificationDecision {
            accepted,
            policy_fingerprint: fixture_policy("exact-anchor-subject-verifier"),
        }
    }
}

static EXACT_ANCHOR_SUBJECT_VERIFIER: ExactAnchorSubjectVerifier = ExactAnchorSubjectVerifier;

fn exact_anchor_claim(statement: &str, regime: ValidityDomain, dataset: &str, hash: &str) -> Claim {
    Claim::anchored("anchor", statement, regime, dataset, hash)
}

#[test]
fn anchoring_hashes_require_exact_typed_external_verification() {
    let pkg = EvidencePackage::new(prov()).with_claim(exact_anchor_claim(
        "matches exact reference data",
        good_regime(),
        "wind-tunnel-2026",
        CANONICAL_DATASET_HASH,
    ));
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::AnchoredSourceRefused {
            why: "anchored-source capability missing",
            ..
        })
    ));
    let capabilities =
        VerificationCapabilities::deny_all().with_anchored_sources(&EXACT_ANCHOR_SUBJECT_VERIFIER);
    let report = pkg
        .verify_with(&capabilities)
        .expect("exact anchor subject verifies");
    assert_eq!(report.breakdown().validated, 1);
    assert_eq!(
        report.receipt().policy_fingerprints().anchored_sources(),
        Some(fixture_policy("exact-anchor-subject-verifier"))
    );

    let different_hash = "1123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";
    let forgeries = [
        EvidencePackage::new(Provenance::new("other-commit", "lock-deadbeef")).with_claim(
            exact_anchor_claim(
                "matches exact reference data",
                good_regime(),
                "wind-tunnel-2026",
                CANONICAL_DATASET_HASH,
            ),
        ),
        EvidencePackage::new(prov()).with_claim(exact_anchor_claim(
            "different assertion",
            good_regime(),
            "wind-tunnel-2026",
            CANONICAL_DATASET_HASH,
        )),
        EvidencePackage::new(prov()).with_claim(exact_anchor_claim(
            "matches exact reference data",
            ValidityDomain::unconstrained().with("Re", 1e5, 4e5),
            "wind-tunnel-2026",
            CANONICAL_DATASET_HASH,
        )),
        EvidencePackage::new(prov()).with_claim(exact_anchor_claim(
            "matches exact reference data",
            good_regime(),
            "other-dataset",
            CANONICAL_DATASET_HASH,
        )),
        EvidencePackage::new(prov()).with_claim(exact_anchor_claim(
            "matches exact reference data",
            good_regime(),
            "wind-tunnel-2026",
            different_hash,
        )),
    ];
    for forged in forgeries {
        assert!(matches!(
            forged.verify_with(&capabilities),
            Err(PackageError::AnchoredSourceRefused {
                why: "rejected by the injected verifier",
                ..
            })
        ));
    }
}

#[test]
fn source_certificate_hashes_require_typed_external_verification() {
    use fs_crosswalk::{PackageConcept, Standard};
    use fs_package::CoverageStatus;

    let pkg = EvidencePackage::new(prov()).with_claim(Claim::from_certificate(
        "source",
        "certified interval",
        1.0,
        2.0,
        "test-solver/cert",
        CANONICAL_DATASET_HASH,
    ));
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::SourceCertificateRefused {
            why: "source-certificate capability missing",
            ..
        })
    ));
    assert!(pkg.color_breakdown().is_err());
    assert!(
        pkg.magnitude_budget().is_err(),
        "an unverified source cannot expose a scientific magnitude"
    );
    assert!(
        fs_package::package_presence(&pkg)
            .iter()
            .all(|row| !row.present()),
        "unverified source bytes must not produce positive coverage"
    );

    let capabilities =
        VerificationCapabilities::deny_all().with_source_certificates(&EXACT_SOURCE_VERIFIER);
    let report = pkg
        .verify_with(&capabilities)
        .expect("typed verifier recognizes the exact certificate subject");
    assert_eq!(report.breakdown().verified, 1);
    assert_eq!(
        pkg.color_breakdown_with(&capabilities)
            .expect("authenticated breakdown")
            .verified,
        1
    );
    let presence = fs_package::package_presence_with(&pkg, &capabilities);
    assert!(
        presence
            .iter()
            .find(|row| row.concept() == PackageConcept::VerifiedColor)
            .expect("verified concept")
            .present()
    );
    assert!(
        presence
            .iter()
            .find(|row| row.concept() == PackageConcept::ClaimOrigin)
            .expect("origin concept")
            .present(),
        "claim-origin presence requires the successful typed verifier above"
    );
    assert!(
        !presence
            .iter()
            .find(|row| row.concept() == PackageConcept::WaiverAuthorization)
            .expect("waiver concept")
            .present()
    );
    let origin_coverage =
        fs_package::package_coverage_with(&pkg, Standard::AsmeVvV40, &capabilities);
    let origin_status = origin_coverage
        .iter()
        .find(|(concept, _, _)| *concept == PackageConcept::ClaimOrigin)
        .expect("claim-origin crosswalk row")
        .1
        .clone();
    assert_eq!(origin_status, CoverageStatus::Covered);
    let decoded = EvidencePackage::from_json_with(&package_json(&pkg), &capabilities)
        .expect("strict transport plus capability authentication");
    assert_eq!(decoded.package(), &pkg);
    assert_eq!(
        decoded.report().receipt().package_root(),
        package_root(&pkg)
    );

    let widened = EvidencePackage::new(prov()).with_claim(Claim::from_certificate(
        "source",
        "certified interval",
        1.0,
        3.0,
        "test-solver/cert",
        CANONICAL_DATASET_HASH,
    ));
    assert!(matches!(
        widened.verify_with(&capabilities),
        Err(PackageError::SourceCertificateRefused {
            why: "rejected by the injected verifier",
            ..
        })
    ));
}

struct PanickingSourceVerifier;

impl SourceCertificateVerifier for PanickingSourceVerifier {
    fn verify(&self, _request: &SourceCertificateRequest<'_>) -> VerificationDecision {
        panic!("hostile external verifier")
    }
}

#[test]
fn verifier_panics_become_structured_refusals() {
    let pkg = EvidencePackage::new(prov()).with_claim(Claim::from_certificate(
        "source",
        "certified interval",
        1.0,
        2.0,
        "test-solver/cert",
        CANONICAL_DATASET_HASH,
    ));
    let capabilities =
        VerificationCapabilities::deny_all().with_source_certificates(&PanickingSourceVerifier);
    assert!(matches!(
        pkg.verify_with(&capabilities),
        Err(PackageError::SourceCertificateRefused {
            why: "verifier callback panicked",
            ..
        })
    ));
}

struct HashWaiverVerifier;

fn waiver_mac(message: &[u8]) -> String {
    fs_blake3::hash_domain("fs-package:test:waiver-mac", message).to_hex()
}

impl WaiverVerifier for HashWaiverVerifier {
    fn verify(&self, mac: &str, message: &[u8]) -> VerificationDecision {
        VerificationDecision {
            accepted: mac == waiver_mac(message),
            policy_fingerprint: fixture_policy("hash-waiver-verifier"),
        }
    }
}

static HASH_WAIVER_VERIFIER: HashWaiverVerifier = HashWaiverVerifier;

fn waived_package(provenance: Provenance, waiver_id: &str, expiry_day: u64) -> EvidencePackage {
    let pending = EvidencePackage::new(provenance).with_claim(Claim::waived(
        "waived",
        "authorized interval",
        Color::Verified { lo: 0.0, hi: 1.0 },
        WaiverGrant {
            waiver_id: waiver_id.to_string(),
            expiry_day,
            mac: "pending-authenticator".to_string(),
        },
    ));
    authorize_waiver(pending, 0)
}

fn authorize_waiver(pending: EvidencePackage, claim_index: usize) -> EvidencePackage {
    let message = pending.waiver_message(claim_index).expect("waiver target");
    let mac = waiver_mac(&message);
    let authorized = pending
        .with_waiver_mac(claim_index, mac)
        .expect("install waiver authenticator");
    assert_eq!(
        authorized.waiver_message(claim_index).as_deref(),
        Ok(message.as_slice())
    );
    authorized
}

#[test]
fn waiver_authentication_enables_capability_aware_coverage() {
    use fs_crosswalk::{PackageConcept, Standard};
    use fs_package::CoverageStatus;

    let pkg = waived_package(prov(), "waiver-2026-01", 200);
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::WaiverRefused {
            why: "waiver capability missing",
            ..
        })
    ));
    assert!(pkg.color_breakdown().is_err());
    assert!(
        fs_package::package_presence(&pkg)
            .iter()
            .all(|row| !row.present())
    );

    let capabilities =
        VerificationCapabilities::deny_all().with_waivers(&HASH_WAIVER_VERIFIER, 199);
    let report = pkg
        .verify_with(&capabilities)
        .expect("unexpired, context-bound waiver verifies");
    assert_eq!(report.breakdown().verified, 0);
    assert_eq!(report.breakdown().waived, 1);
    assert_eq!(report.magnitude_budget().verified_width.to_bits(), 0);
    assert_eq!(report.magnitude_budget().waived_unquantified, 1);
    assert_eq!(report.receipt().package_root(), report.merkle_root());
    assert_eq!(report.receipt().waiver_day(), Some(199));
    assert_eq!(
        report.receipt().policy_fingerprints().waivers(),
        Some(fixture_policy("hash-waiver-verifier"))
    );
    assert_eq!(
        report.receipt().admissions()[0].class(),
        AdmissionClass::WaiverDependent
    );
    assert_eq!(
        report.receipt().admissions()[0].origin_kind(),
        fs_package::AdmissionOriginKind::AuthenticatedWaiver
    );
    assert_eq!(report.receipt().admissions()[0].direct_waiver(), Some(0));
    assert!(report.receipt().admissions()[0].waiver_parents().is_empty());
    assert_eq!(report.receipt().waiver_registry().len(), 1);
    assert_eq!(
        report.receipt().waiver_registry()[0].waiver_id(),
        "waiver-2026-01"
    );
    let different_day =
        VerificationCapabilities::deny_all().with_waivers(&HASH_WAIVER_VERIFIER, 198);
    let different_day_report = pkg
        .verify_with(&different_day)
        .expect("earlier day is also inside the grant");
    assert_ne!(
        report.receipt().receipt_hash(),
        different_day_report.receipt().receipt_hash(),
        "the decision receipt binds the explicit waiver clock"
    );
    let presence = fs_package::package_presence_with(&pkg, &capabilities);
    assert!(
        !presence
            .iter()
            .find(|row| row.concept() == PackageConcept::ClaimOrigin)
            .expect("claim-origin concept")
            .present(),
        "administrative waiver authority must not become scientific source traceability"
    );
    assert!(
        presence
            .iter()
            .find(|row| row.concept() == PackageConcept::WaiverAuthorization)
            .expect("waiver-authorization concept")
            .present(),
        "successful waiver authentication belongs only in waiver coverage"
    );
    let waiver_coverage =
        fs_package::package_coverage_with(&pkg, Standard::FaaEasaCbA, &capabilities);
    let waiver_status = waiver_coverage
        .iter()
        .find(|(concept, _, _)| *concept == PackageConcept::WaiverAuthorization)
        .expect("waiver crosswalk row")
        .1
        .clone();
    assert_eq!(waiver_status, CoverageStatus::Covered);
    let asme_coverage = fs_package::package_coverage_with(&pkg, Standard::AsmeVvV40, &capabilities);
    let asme_status = asme_coverage
        .iter()
        .find(|(concept, _, _)| *concept == PackageConcept::WaiverAuthorization)
        .expect("waiver no-counterpart row")
        .1
        .clone();
    assert_eq!(asme_status, CoverageStatus::NoClaim);
    let decoded = EvidencePackage::from_json_with(&package_json(&pkg), &capabilities)
        .expect("waiver capability is available at the JSON boundary");
    assert_eq!(decoded.package(), &pkg);
    assert_eq!(
        decoded.report().receipt().package_root(),
        package_root(&pkg)
    );
}

fn pending_waiver_claim(id: &str, color: Color) -> Claim {
    Claim::waived(
        id,
        format!("{id} administrative exception"),
        color,
        WaiverGrant {
            waiver_id: "waiver-taint-root".to_string(),
            expiry_day: 400,
            mac: "pending-authenticator".to_string(),
        },
    )
}

#[test]
fn one_parent_derivation_cannot_launder_a_waiver() {
    let color = Color::Verified { lo: 0.0, hi: 1.0 };
    let pending = EvidencePackage::new(prov())
        .with_claim(pending_waiver_claim("waived-root", color.clone()))
        .with_claim(Claim::derived(
            "child",
            "identity derivation",
            color,
            vec![0],
            IntervalOp::Hull,
            CANONICAL_DATASET_HASH,
        ));
    let pkg = authorize_waiver(pending, 0);
    let capabilities = VerificationCapabilities::deny_all()
        .with_waivers(&HASH_WAIVER_VERIFIER, 300)
        .with_derivations(&FIXTURE_DERIVATION_VERIFIER);
    let report = pkg
        .verify_with(&capabilities)
        .expect("waiver authenticates");
    assert_eq!(report.breakdown().verified, 0);
    assert_eq!(report.breakdown().waived, 2);
    assert_eq!(report.magnitude_budget().verified_width.to_bits(), 0);
    assert_eq!(report.magnitude_budget().waived_unquantified, 2);
    assert_eq!(report.receipt().admissions()[1].waiver_parents(), [0]);
    assert_eq!(report.receipt().waiver_registry().len(), 1);
}

#[test]
fn multihop_derivation_preserves_the_original_waiver_identity() {
    let color = Color::Verified { lo: 0.0, hi: 1.0 };
    let pending = EvidencePackage::new(prov())
        .with_claim(pending_waiver_claim("waived-root", color.clone()))
        .with_claim(Claim::derived(
            "child",
            "first hop",
            color.clone(),
            vec![0],
            IntervalOp::Hull,
            CANONICAL_DATASET_HASH,
        ))
        .with_claim(Claim::derived(
            "grandchild",
            "second hop",
            color,
            vec![1],
            IntervalOp::Hull,
            CANONICAL_DATASET_HASH,
        ));
    let pkg = authorize_waiver(pending, 0);
    let capabilities = VerificationCapabilities::deny_all()
        .with_waivers(&HASH_WAIVER_VERIFIER, 300)
        .with_derivations(&FIXTURE_DERIVATION_VERIFIER);
    let report = pkg
        .verify_with(&capabilities)
        .expect("waiver authenticates");
    assert_eq!(report.breakdown().waived, 3);
    assert!(
        report
            .receipt()
            .admissions()
            .iter()
            .all(|decision| decision.class() == AdmissionClass::WaiverDependent)
    );
    assert_eq!(report.receipt().admissions()[1].waiver_parents(), [0]);
    assert_eq!(report.receipt().admissions()[2].waiver_parents(), [1]);
    assert_eq!(report.receipt().waiver_registry().len(), 1);
}

#[test]
fn long_waiver_ids_are_interned_once_across_deep_taint_chains() {
    let color = Color::Verified { lo: 0.0, hi: 1.0 };
    let waiver_id = format!("waiver-{}", "w".repeat(8 * 1024));
    let mut pending = EvidencePackage::new(prov()).with_claim(Claim::waived(
        "waived-root",
        "long waiver identity",
        color.clone(),
        WaiverGrant {
            waiver_id: waiver_id.clone(),
            expiry_day: 400,
            mac: "pending-authenticator".to_string(),
        },
    ));
    for index in 1..32 {
        pending = pending.with_claim(Claim::derived(
            format!("derived-{index}"),
            format!("taint hop {index}"),
            color.clone(),
            vec![index - 1],
            IntervalOp::Hull,
            CANONICAL_DATASET_HASH,
        ));
    }
    let package = authorize_waiver(pending, 0);
    let capabilities = VerificationCapabilities::deny_all()
        .with_waivers(&HASH_WAIVER_VERIFIER, 300)
        .with_derivations(&FIXTURE_DERIVATION_VERIFIER);
    let report = package
        .verify_with(&capabilities)
        .expect("deep waiver chain verifies");
    assert_eq!(report.receipt().waiver_registry().len(), 1);
    assert_eq!(report.receipt().waiver_registry()[0].waiver_id(), waiver_id);
    for (index, admission) in report.receipt().admissions().iter().enumerate().skip(1) {
        assert_eq!(admission.waiver_parents(), [index - 1]);
        assert_eq!(admission.direct_waiver(), None);
    }
    assert!(report.receipt().validate_hash());
}

#[test]
fn multiparent_derivation_is_tainted_if_any_parent_depends_on_a_waiver() {
    let waived_color = Color::Verified { lo: 0.0, hi: 1.0 };
    let clean_color = Color::Estimated {
        estimator: "clean-estimator".to_string(),
        dispersion: 2.0,
    };
    let derived_color = fs_evidence::compose(&waived_color, &clean_color, IntervalOp::Add);
    let pending = EvidencePackage::new(prov())
        .with_claim(pending_waiver_claim("waived-root", waived_color))
        .with_claim(Claim::estimated(
            "clean",
            "independent estimate",
            "clean-estimator",
            2.0,
        ))
        .with_claim(Claim::derived(
            "mixed",
            "waived plus clean",
            derived_color,
            vec![0, 1],
            IntervalOp::Add,
            CANONICAL_DATASET_HASH,
        ));
    let pkg = authorize_waiver(pending, 0);
    let capabilities = VerificationCapabilities::deny_all()
        .with_waivers(&HASH_WAIVER_VERIFIER, 300)
        .with_derivations(&FIXTURE_DERIVATION_VERIFIER);
    let report = pkg
        .verify_with(&capabilities)
        .expect("waiver authenticates");
    assert_eq!(report.breakdown().estimated, 1);
    assert_eq!(report.breakdown().waived, 2);
    assert_eq!(
        report.magnitude_budget().estimated_dispersion.to_bits(),
        2.0f64.to_bits()
    );
    assert_eq!(report.magnitude_budget().waived_unquantified, 2);
    assert_eq!(
        report.receipt().admissions()[1].class(),
        AdmissionClass::Scientific
    );
    assert_eq!(
        report.receipt().admissions()[2].class(),
        AdmissionClass::WaiverDependent
    );
    assert_eq!(report.receipt().admissions()[2].waiver_parents(), [0]);
}

#[test]
fn waiver_authentication_binds_package_context_id_expiry_and_claim() {
    let pkg = waived_package(prov(), "waiver-2026-01", 200);
    let capabilities =
        VerificationCapabilities::deny_all().with_waivers(&HASH_WAIVER_VERIFIER, 199);

    let old_mac = pkg
        .declared_claims_unverified()
        .iter()
        .find_map(|claim| match claim.declared_origin_unverified() {
            fs_package::ClaimOrigin::AuthenticatedWaiver(grant) => Some(grant.mac.clone()),
            _ => None,
        })
        .unwrap_or_default();
    assert!(!old_mac.is_empty(), "waiver fixture changed origin");
    for replay in [
        EvidencePackage::new(Provenance::new("different-commit", "lock-deadbeef")).with_claim(
            Claim::waived(
                "waived",
                "authorized interval",
                Color::Verified { lo: 0.0, hi: 1.0 },
                WaiverGrant {
                    waiver_id: "waiver-2026-01".to_string(),
                    expiry_day: 200,
                    mac: old_mac.clone(),
                },
            ),
        ),
        EvidencePackage::new(prov()).with_claim(Claim::waived(
            "waived",
            "different assertion",
            Color::Verified { lo: 0.0, hi: 1.0 },
            WaiverGrant {
                waiver_id: "waiver-2026-01".to_string(),
                expiry_day: 200,
                mac: old_mac.clone(),
            },
        )),
        EvidencePackage::new(prov()).with_claim(Claim::waived(
            "waived",
            "authorized interval",
            Color::Verified { lo: 0.0, hi: 1.0 },
            WaiverGrant {
                waiver_id: "different-waiver".to_string(),
                expiry_day: 200,
                mac: old_mac.clone(),
            },
        )),
        EvidencePackage::new(prov()).with_claim(Claim::waived(
            "waived",
            "authorized interval",
            Color::Verified { lo: 0.0, hi: 1.0 },
            WaiverGrant {
                waiver_id: "waiver-2026-01".to_string(),
                expiry_day: 201,
                mac: old_mac.clone(),
            },
        )),
    ] {
        assert!(matches!(
            replay.verify_with(&capabilities),
            Err(PackageError::WaiverRefused {
                why: "rejected by the injected verifier",
                ..
            })
        ));
    }

    let expired = VerificationCapabilities::deny_all().with_waivers(&HASH_WAIVER_VERIFIER, 201);
    assert!(matches!(
        pkg.verify_with(&expired),
        Err(PackageError::WaiverRefused { why: "expired", .. })
    ));
}

#[test]
fn waiver_ids_are_unique_and_builder_targets_are_typed() {
    let duplicate = EvidencePackage::new(prov())
        .with_claim(Claim::waived(
            "first",
            "first waiver",
            Color::Estimated {
                estimator: "probe".to_string(),
                dispersion: 1.0,
            },
            WaiverGrant {
                waiver_id: "same-waiver".to_string(),
                expiry_day: 200,
                mac: "opaque-one".to_string(),
            },
        ))
        .with_claim(Claim::waived(
            "second",
            "second waiver",
            Color::Estimated {
                estimator: "probe".to_string(),
                dispersion: 1.0,
            },
            WaiverGrant {
                waiver_id: "same-waiver".to_string(),
                expiry_day: 200,
                mac: "opaque-two".to_string(),
            },
        ));
    assert!(matches!(
        duplicate.verify(),
        Err(PackageError::DuplicateWaiverId {
            first_claim,
            duplicate_claim,
            ..
        }) if first_claim == "first" && duplicate_claim == "second"
    ));

    let ordinary = EvidencePackage::new(prov()).with_claim(estimated("ordinary"));
    assert!(matches!(
        ordinary.with_waiver_mac(0, "mac"),
        Err(PackageError::InvalidWaiverTarget { index: 0 })
    ));
}

#[test]
fn multiple_waivers_share_one_mac_independent_authorization_context() {
    let pending = EvidencePackage::new(prov())
        .with_claim(Claim::waived(
            "first",
            "first authorization",
            Color::Verified { lo: 0.0, hi: 1.0 },
            WaiverGrant {
                waiver_id: "waiver-one".to_string(),
                expiry_day: 200,
                mac: "pending-one".to_string(),
            },
        ))
        .with_claim(Claim::waived(
            "second",
            "second authorization",
            Color::Verified { lo: 2.0, hi: 3.0 },
            WaiverGrant {
                waiver_id: "waiver-two".to_string(),
                expiry_day: 200,
                mac: "pending-two".to_string(),
            },
        ));
    let first_message = pending.waiver_message(0).expect("first waiver target");
    let second_message = pending.waiver_message(1).expect("second waiver target");
    let package = pending
        .with_waiver_mac(0, waiver_mac(&first_message))
        .expect("first authenticator")
        .with_waiver_mac(1, waiver_mac(&second_message))
        .expect("second authenticator");
    assert_eq!(
        package.waiver_message(0).as_deref(),
        Ok(first_message.as_slice())
    );
    assert_eq!(
        package.waiver_message(1).as_deref(),
        Ok(second_message.as_slice())
    );
    let capabilities =
        VerificationCapabilities::deny_all().with_waivers(&HASH_WAIVER_VERIFIER, 199);
    package
        .verify_with(&capabilities)
        .expect("both waivers authenticate against the shared context");
}

#[test]
fn origin_transport_and_machine_identity_boundaries_fail_closed() {
    let oversized_source = EvidencePackage::new(prov()).with_claim(Claim::from_certificate(
        "source",
        "bounded",
        0.0,
        1.0,
        "p".repeat(MAX_JSON_STRING_BYTES + 1),
        CANONICAL_DATASET_HASH,
    ));
    assert!(matches!(
        oversized_source.verify(),
        Err(PackageError::TransportLimit { what, .. }) if what == "origin.producer"
    ));

    let oversized_waiver = EvidencePackage::new(prov()).with_claim(Claim::waived(
        "waived",
        "bounded",
        Color::Verified { lo: 0.0, hi: 1.0 },
        WaiverGrant {
            waiver_id: "waiver-id".to_string(),
            expiry_day: 200,
            mac: "m".repeat(MAX_JSON_STRING_BYTES + 1),
        },
    ));
    assert!(matches!(
        oversized_waiver.verify(),
        Err(PackageError::TransportLimit { what, .. }) if what == "origin.mac"
    ));

    let placeholder_id = EvidencePackage::new(prov()).with_claim(estimated("TODO"));
    assert!(matches!(
        placeholder_id.verify(),
        Err(PackageError::InvalidClaimId {
            reason: "placeholder",
            ..
        })
    ));
    let padded_provenance =
        EvidencePackage::new(Provenance::new(" commit", "lock")).with_claim(estimated("e"));
    assert!(matches!(
        padded_provenance.verify(),
        Err(PackageError::InvalidIdentity {
            field: "provenance.code_version",
            reason: "surrounding-whitespace",
            ..
        })
    ));
    let placeholder_estimator =
        EvidencePackage::new(prov()).with_claim(Claim::estimated("e", "estimate", "unknown", 1.0));
    assert!(matches!(
        placeholder_estimator.verify(),
        Err(PackageError::InvalidIdentity {
            field: "color.estimator",
            reason: "placeholder",
            ..
        })
    ));
    let placeholder_producer = EvidencePackage::new(prov()).with_claim(Claim::from_certificate(
        "source",
        "bounded",
        0.0,
        1.0,
        "pending",
        CANONICAL_DATASET_HASH,
    ));
    assert!(matches!(
        placeholder_producer.verify(),
        Err(PackageError::InvalidOrigin { .. })
    ));
}
