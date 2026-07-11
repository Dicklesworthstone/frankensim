//! Battery for evidence packages (addendum Proposal 12). Covers a complete
//! mixed-color package, the all-estimated boundary (still valid, round-trips),
//! completeness failures (validated claim missing regime / dataset, verified
//! claim with a bad interval), Merkle content-addressing (determinism + tamper
//! detection), the format-version gate, optional signature, the color
//! breakdown, and deterministic JSON.

use fs_evidence::{Color, ValidityDomain};
use fs_package::{
    AnchorRecord, Claim, EvidencePackage, FalsifierRecord, MAX_JSON_CONTAINER_ITEMS,
    MAX_JSON_DEPTH, MAX_JSON_NUMBER_BYTES, MAX_JSON_STRING_BYTES, PackageError, Provenance,
};

const CANONICAL_DATASET_HASH: &str =
    "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

fn prov() -> Provenance {
    Provenance::new("commit-abc123", "lock-deadbeef")
}
fn verified(id: &str) -> Claim {
    Claim::from_certificate(
        id,
        format!("{id}: stress <= sigma*"),
        -1.0,
        1.0,
        "test-solver/cert",
        CANONICAL_DATASET_HASH,
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

fn assert_serialized_refuses(pkg: &EvidencePackage) {
    assert!(
        EvidencePackage::from_json(&pkg.to_json()).is_err(),
        "serialized package must refuse the same semantics as verify(): {pkg:?}"
    );
}

#[test]
fn a_complete_mixed_color_package_verifies() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(validated("c2", good_regime(), "wind-tunnel-2026"))
        .with_claim(estimated("c3"));
    let report = pkg.verify().expect("complete package verifies");
    assert_eq!(report.claims, 3);
    assert_eq!(report.breakdown.verified, 1);
    assert_eq!(report.breakdown.validated, 1);
    assert_eq!(report.breakdown.estimated, 1);
    assert_eq!(report.merkle_root, pkg.merkle_root());
    assert_ne!(report.merkle_root.as_bytes(), &[0u8; 32]);
}

#[test]
fn an_all_estimated_package_is_still_valid_and_round_trips() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(estimated("e1"))
        .with_claim(estimated("e2"));
    let report = pkg.verify().expect("all-estimated is honest, not invalid");
    assert_eq!(report.breakdown.estimated, 2);
    assert_eq!(report.breakdown.verified, 0);
    let json = pkg.to_json();
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
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::IncompleteValidatedClaim {
            missing: "dataset",
            ..
        })
    ));
}

#[test]
fn a_verified_claim_with_a_bad_interval_fails() {
    let pkg = EvidencePackage::new(prov()).with_claim(Claim::from_certificate("v", "backwards", 5.0, 1.0, "test-solver/cert", CANONICAL_DATASET_HASH));
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
        let pkg = EvidencePackage::new(prov()).with_claim(Claim::from_certificate("claim", statement, 0.0, 1.0, "test-solver/cert", CANONICAL_DATASET_HASH));
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
    let blank_estimator = EvidencePackage::new(prov()).with_claim(Claim::estimated("e", "missing estimator identity", " ".to_string(), 1.0));
    assert!(matches!(
        blank_estimator.verify(),
        Err(PackageError::IncompleteEstimatedClaim {
            missing: "estimator",
            ..
        })
    ));
    assert_serialized_refuses(&blank_estimator);

    for dispersion in [-1.0, f64::NEG_INFINITY, f64::NAN] {
        let pkg = EvidencePackage::new(prov())
            .with_claim(Claim::estimated("e", "invalid dispersion", "probe", dispersion));
        assert!(matches!(
            pkg.verify(),
            Err(PackageError::InvalidEstimatedDispersion { .. })
        ));
        assert_serialized_refuses(&pkg);
    }

    let explicitly_unbounded = EvidencePackage::new(prov()).with_claim(Claim::estimated("unbounded", "honest no-spread-claim sentinel", "regime-exit".to_string(), f64::INFINITY));
    assert!(explicitly_unbounded.verify().is_ok());
    assert!(
        explicitly_unbounded
            .magnitude_budget()
            .estimated_dispersion
            .is_infinite()
    );
    assert_eq!(
        EvidencePackage::from_json(&explicitly_unbounded.to_json()).unwrap(),
        explicitly_unbounded
    );
}

#[test]
fn in_memory_and_serialized_magnitude_gates_are_identical() {
    let width_overflow = EvidencePackage::new(prov()).with_claim(Claim::from_certificate("wide", "finite endpoints whose width overflows", -f64::MAX, f64::MAX, "test-solver/cert", CANONICAL_DATASET_HASH));
    assert!(matches!(
        width_overflow.verify(),
        Err(PackageError::MagnitudeOverflow {
            component: "verified_width",
            ..
        })
    ));
    assert_serialized_refuses(&width_overflow);

    let dispersion_overflow = EvidencePackage::new(prov())
        .with_claim(Claim::estimated("d1", "large finite dispersion", "probe-1".to_string(), f64::MAX))
        .with_claim(Claim::estimated("d2", "second large finite dispersion", "probe-2".to_string(), f64::MAX));
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
        .with_claim(Claim::from_certificate("wide-finite", "large but finite interval width", 0.0, large, "test-solver/cert", CANONICAL_DATASET_HASH))
        .with_claim(Claim::estimated("spread-finite", "large but finite estimated spread", "probe".to_string(), large));
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
    let inverted_json = ordered.to_json().replace(&ordered_pair, &inverted_pair);
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
            })
            .with_falsifier(FalsifierRecord {
                name: "exhaustive-probe".to_string(),
                attempts: u64::MAX,
                refuted: false,
                detail: "full-width counter".to_string(),
            }),
    );
    let json = pkg.to_json();
    assert!(json.contains("\"attempts\":9007199254740993"));
    assert!(json.contains("\"attempts\":18446744073709551615"));
    let back = EvidencePackage::from_json(&json).expect("full-width u64 values parse exactly");
    assert_eq!(back, pkg);
    assert_eq!(
        back.claims[0].falsifiers()[0].attempts,
        first_unrepresentable_integer
    );
    assert_eq!(back.claims[0].falsifiers()[1].attempts, u64::MAX);

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
        },
        "name",
    );
    assert_invalid(
        FalsifierRecord {
            name: "interval-probe".to_string(),
            attempts: 0,
            refuted: false,
            detail: "no work ran".to_string(),
        },
        "attempts",
    );
    assert_invalid(
        FalsifierRecord {
            name: "interval-probe".to_string(),
            attempts: 1,
            refuted: false,
            detail: "  ".to_string(),
        },
        "detail",
    );
    for placeholder in ["TODO", " placeholder ", "N/A", "not run", "unknown"] {
        assert_invalid(
            FalsifierRecord {
                name: placeholder.to_string(),
                attempts: 1,
                refuted: false,
                detail: "no violation observed".to_string(),
            },
            "name",
        );
        assert_invalid(
            FalsifierRecord {
                name: "interval-probe".to_string(),
                attempts: 1,
                refuted: false,
                detail: placeholder.to_string(),
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

    let oversized_in_memory = EvidencePackage::new(prov()).with_claim(Claim::estimated("large", "x".repeat(MAX_JSON_STRING_BYTES + 1), "probe".to_string(), 1.0));
    assert!(matches!(
        oversized_in_memory.verify(),
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
    assert_eq!(build().merkle_root(), build().merkle_root());
    assert_eq!(
        build().merkle_root().to_hex(),
        "61a09cc114970a8400b77c52192e4909cb85c98e245cd8b8aa056f6d8f1e134a",
        "schema-v4 package-root ABI; changes require a format-version migration"
    );
    // tampering with a claim changes the root.
    let tampered = EvidencePackage::new(prov())
        .with_claim(Claim::from_certificate("c1", "TAMPERED", -1.0, 1.0, "test-solver/cert", CANONICAL_DATASET_HASH))
        .with_claim(estimated("c2"));
    assert_ne!(build().merkle_root(), tampered.merkle_root());
}

#[test]
fn the_merkle_root_covers_reproducibility_provenance() {
    let pkg = EvidencePackage::new(prov()).with_claim(verified("c1"));
    let changed_code = EvidencePackage::new(Provenance::new("commit-other", "lock-deadbeef"))
        .with_claim(verified("c1"));
    let changed_lock = EvidencePackage::new(Provenance::new("commit-abc123", "lock-other"))
        .with_claim(verified("c1"));

    assert_ne!(pkg.merkle_root(), changed_code.merkle_root());
    assert_ne!(pkg.merkle_root(), changed_lock.merkle_root());
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
    assert!(unsigned.to_json().contains("\"signature\":null"));
    let signed = unsigned.clone().signed("ed25519:deadbeef");
    // signing does not change the content address (detached).
    assert_eq!(unsigned.merkle_root(), signed.merkle_root());
    assert!(signed.to_json().contains("ed25519:deadbeef"));
}

#[test]
fn json_is_deterministic_and_carries_the_root() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(validated("c2", good_regime(), "wt-2026"));
    let j = pkg.to_json();
    assert_eq!(j, pkg.to_json());
    assert!(j.starts_with('{') && j.ends_with('}'));
    assert!(j.contains(&pkg.merkle_root().to_hex()));
    assert!(j.contains("\"format_version\":4"), "schema v4 (7uq9)");
    // v3 carries COMPLETE payloads, not just rank labels.
    assert!(j.contains("\"lo_bits\":") && j.contains("\"dataset\":"));
}

/// qmao.6.1 — the schema-v3 round trip and its fail-closed walls: a
/// package decode-encodes stably (golden), floats travel bit-exactly,
/// and hostile/missing/unknown/non-finite/forged/tampered inputs each
/// refuse with a structured parse error.
#[test]
fn v3_round_trip_and_fail_closed_walls() {
    let pkg = EvidencePackage::new(Provenance::new("frankensim@abc123", "lock:deadbeef"))
        .with_claim(Claim::from_certificate("c-verified", "tip deflection within bound", 0.1875, 0.25, "test-solver/cert", CANONICAL_DATASET_HASH))
        .with_claim(Claim::anchored(
            "c-validated",
            "k-epsilon matched within regime",
            fs_evidence::ValidityDomain::unconstrained().with("reynolds", 1e3, 1e5),
            "tunnel-run-9",
            CANONICAL_DATASET_HASH,
        ))
        .with_claim(Claim::estimated("c-estimated", "surrogate prediction", "pod-deim".to_string(), 0.02))
        .signed("test-key/1234abcd");
    // Golden decode-encode stability: parse(to_json) == pkg, and the
    // re-emission is byte-identical (semantic AND textual round trip).
    let json = pkg.to_json();
    let back = EvidencePackage::from_json(&json).expect("canonical JSON parses");
    assert_eq!(back, pkg, "semantic round trip");
    assert_eq!(back.to_json(), json, "textual round trip");
    let leading_zero = json.replacen("\"format_version\":4", "\"format_version\":04", 1);
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
    let root = pkg.merkle_root().to_hex();
    let short_root = json.replacen(&format!("\"{root}\""), &format!("\"{}\"", &root[1..]), 1);
    let err = EvidencePackage::from_json(&short_root).expect_err("short root refused");
    assert!(err.why.contains("64 hex chars"), "{err}");
    let raw_control = json.replacen("surrogate prediction", "surrogate\nprediction", 1);
    let err = EvidencePackage::from_json(&raw_control).expect_err("raw control refused");
    assert!(err.why.contains("control character"), "{err}");
}

/// qmao.6.1 — the magnitude budget attributes ERROR MAGNITUDES, not
/// claim counts, and reconciles with an independent recomputation.
#[test]
fn magnitude_budget_reconciles() {
    let pkg = EvidencePackage::new(Provenance::new("v", "l"))
        .with_claim(Claim::from_certificate("a", "s", 0.0, 0.5, "test-solver/cert", CANONICAL_DATASET_HASH))
        .with_claim(Claim::from_certificate("b", "s", 1.0, 1.25, "test-solver/cert", CANONICAL_DATASET_HASH))
        .with_claim(Claim::estimated("c", "s", "e".to_string(), 0.125))
        .with_claim(Claim::anchored(
            "d",
            "s",
            fs_evidence::ValidityDomain::unconstrained().with("re", 1.0, 2.0),
            "ds",
            CANONICAL_DATASET_HASH,
        ));
    let mb = pkg.magnitude_budget();
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
    use fs_evidence::Color;
    use fs_package::{Claim, EvidencePackage, Provenance};
    use fs_package::{CoverageStatus, package_coverage, package_presence};
    // Unsigned, verified-only package: no validated claims, no regime
    // tags, no datasets, no signature.
    let pkg = EvidencePackage::new(Provenance::new("v", "lock")).with_claim(Claim::from_certificate("c", "bounded", 0.0, 1.0, "test-solver/cert", CANONICAL_DATASET_HASH));
    let presence = package_presence(&pkg);
    let get = |c: PackageConcept| presence.iter().find(|p| p.concept == c).unwrap();
    assert!(get(PackageConcept::VerifiedColor).present);
    assert!(!get(PackageConcept::ValidatedColor).present);
    assert!(!get(PackageConcept::RegimeTag).present);
    assert!(!get(PackageConcept::AnchoringDataset).present);
    assert!(!get(PackageConcept::Signature).present);
    assert!(
        !get(PackageConcept::FalsifierLog).present,
        "absent falsifier records can never read as present"
    );
    for standard in Standard::ALL {
        for (concept, status, _why) in package_coverage(&pkg, standard) {
            if matches!(status, CoverageStatus::Covered) {
                let p = presence.iter().find(|p| p.concept == concept).unwrap();
                assert!(
                    p.present,
                    "{concept:?} covered without evidence for {standard:?}"
                );
            }
            if concept == PackageConcept::FalsifierLog {
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
        }))
        .signed("raw-detached-signature");
    assert!(package_presence(&refuted).iter().all(|row| !row.present));
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
    let signature = package_presence(&merely_signed)
        .into_iter()
        .find(|row| row.concept == PackageConcept::Signature)
        .expect("signature concept judged");
    assert!(!signature.present, "{}", signature.why);
    for standard in Standard::ALL {
        let status = package_coverage(&merely_signed, standard)
            .into_iter()
            .find(|(concept, _, _)| *concept == PackageConcept::Signature)
            .expect("signature is mapped for every standard")
            .1;
        assert!(!matches!(status, CoverageStatus::Covered));
    }
}

#[test]
fn dataset_coverage_requires_a_matching_valid_anchor() {
    use fs_crosswalk::PackageConcept;
    use fs_package::package_presence;

    let unrelated = EvidencePackage::new(prov()).with_claim(
        validated("v", good_regime(), "wind-tunnel-2026")
            .with_anchor("different-dataset", CANONICAL_DATASET_HASH),
    );
    unrelated
        .verify()
        .expect("unrelated anchor is structurally valid");
    let unrelated_presence = package_presence(&unrelated);
    let anchor = unrelated_presence
        .iter()
        .find(|row| row.concept == PackageConcept::AnchoringDataset)
        .expect("anchor concept judged");
    assert!(!anchor.present, "{}", anchor.why);

    let matching = EvidencePackage::new(prov()).with_claim(
        validated("v", good_regime(), "wind-tunnel-2026")
            .with_anchor("wind-tunnel-2026", CANONICAL_DATASET_HASH),
    );
    matching.verify().expect("matching anchor verifies");
    let matching_presence = package_presence(&matching);
    let anchor = matching_presence
        .iter()
        .find(|row| row.concept == PackageConcept::AnchoringDataset)
        .expect("anchor concept judged");
    assert!(anchor.present, "{}", anchor.why);
}

/// xfxq (schema v3) — composition receipts re-run solver-free, refuted
/// falsifiers fail, anchors and falsifier logs travel with claims, and
/// every new field round-trips through the strict parser bound into the
/// content address.
#[test]
fn v3_receipts_falsifiers_anchors() {
    use fs_evidence::IntervalOp;
    use fs_package::PackageError;
    let ve = |lo: f64, hi: f64| Color::Verified { lo, hi };
    // A well-formed derived package: c = a + b with a valid receipt.
    let good = EvidencePackage::new(Provenance::new("v", "l"))
        .with_claim(Claim::from_certificate("a", "left", 1.0, 2.0, "test-solver/cert", CANONICAL_DATASET_HASH))
        .with_claim(Claim::from_certificate("b", "right", 10.0, 20.0, "test-solver/cert", CANONICAL_DATASET_HASH))
        .with_claim(
            Claim::derived(
                "c",
                "sum",
                fs_evidence::compose(&ve(1.0, 2.0), &ve(10.0, 20.0), IntervalOp::Add),
                vec![0, 1],
                IntervalOp::Add,
            )
            .with_falsifier(FalsifierRecord {
                name: "interval-probe".to_string(),
                attempts: 512,
                refuted: false,
                detail: "no violation found".to_string(),
            })
            .with_anchor("wt-2026-run9", CANONICAL_DATASET_HASH),
        );
    good.verify().expect("receipt re-derives");
    // Round trip: the v3 fields survive the strict parser bit-for-bit.
    let back = EvidencePackage::from_json(&good.to_json()).expect("v3 parses");
    assert_eq!(back, good);
    assert_eq!(back.claims[2].falsifiers()[0].attempts, 512);
    assert_eq!(back.claims[2].anchors()[0].dataset_id, "wt-2026-run9");
    // FORGED receipt: claiming Verified while a parent is Estimated —
    // the re-run composition cannot reproduce it (semantic catch, not
    // just the content-address catch).
    let forged = EvidencePackage::new(Provenance::new("v", "l"))
        .with_claim(Claim::estimated("a", "shaky", "guess".to_string(), 0.5))
        .with_claim(Claim::from_certificate("b", "solid", 1.0, 2.0, "test-solver/cert", CANONICAL_DATASET_HASH))
        .with_claim(
            Claim::derived("c", "laundered", ve(2.0, 4.0), vec![0, 1], IntervalOp::Add),
        );
    assert!(matches!(
        forged.verify(),
        Err(PackageError::ReceiptMismatch { claim }) if claim == "c"
    ));
    // Forward/self parent references refuse.
    let cyclic = EvidencePackage::new(Provenance::new("v", "l"))
        .with_claim(Claim::derived("a", "s", ve(0.0, 1.0), vec![0], IntervalOp::Hull));
    assert!(matches!(
        cyclic.verify(),
        Err(PackageError::BadReceiptParent { parent: 0, .. })
    ));
    // A refuted falsifier fails the whole claim.
    let refuted = EvidencePackage::new(Provenance::new("v", "l")).with_claim(
        Claim::from_certificate("a", "wrong", 0.0, 1.0, "test-solver/cert", CANONICAL_DATASET_HASH).with_falsifier(FalsifierRecord {
            name: "adversary".to_string(),
            attempts: 3,
            refuted: true,
            detail: "counterexample at x=0.7".to_string(),
        }),
    );
    assert!(matches!(
        refuted.verify(),
        Err(PackageError::RefutedClaim { falsifier, .. }) if falsifier == "adversary"
    ));
    // Crosswalk: falsifier logs are now REPRESENTABLE — records present
    // flips the concept to present.
    let presence = fs_package::package_presence(&good);
    let fal = presence
        .iter()
        .find(|p| p.concept == fs_crosswalk::PackageConcept::FalsifierLog)
        .unwrap();
    assert!(fal.present, "{}", fal.why);
    // Tampering with a falsifier flag flips the content address (bound).
    let json = good.to_json();
    let tampered = json.replace("\"refuted\":false", "\"refuted\":true");
    assert!(
        EvidencePackage::from_json(&tampered).is_err(),
        "root binds v3 fields"
    );
}

/// 7uq9 (schema v4) — the content address is a domain-separated 32-byte
/// BLAKE3 root; legacy v3 transports and legacy 16-hex FNV roots are
/// refused with messages that NAME the incompatibility.
#[test]
fn v4_blake3_root_refuses_legacy_transports() {
    let pkg = EvidencePackage::new(prov()).with_claim(verified("c1"));
    let json = pkg.to_json();
    let root = pkg.merkle_root().to_hex();
    assert_eq!(root.len(), 64, "32-byte BLAKE3 root renders as 64 hex");
    assert!(json.contains(&format!("\"merkle_root\":\"{root}\"")));

    // A v3 transport is refused BY VERSION before any field is read.
    let v3 = json.replacen("\"format_version\":4", "\"format_version\":3", 1);
    let err = EvidencePackage::from_json(&v3).expect_err("v3 refused");
    assert!(err.why.contains("unsupported version 3"), "{err}");

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
    assert!(err.why.contains("non-hex"), "{err}");

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
    let root = pkg.merkle_root();
    // Not the undomained hash of the canonical JSON or of the claim.
    assert_ne!(root, fs_blake3::hash_bytes(pkg.to_json().as_bytes()));
    assert_ne!(
        root,
        EvidencePackage::new(prov()).merkle_root(),
        "header claim count distinguishes zero and one claim"
    );
    assert_ne!(
        root,
        EvidencePackage::new(prov())
            .with_claim(verified("c1"))
            .with_claim(verified("c1"))
            .merkle_root(),
        "duplicating the final claim changes the tree identity"
    );
    // Stable across recomputation and bound to claim ORDER.
    let swapped = EvidencePackage::new(prov())
        .with_claim(verified("c2"))
        .with_claim(verified("c1"));
    let ordered = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(verified("c2"));
    assert_ne!(swapped.merkle_root(), ordered.merkle_root());
}
