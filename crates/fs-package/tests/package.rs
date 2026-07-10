//! Battery for evidence packages (addendum Proposal 12). Covers a complete
//! mixed-color package, the all-estimated boundary (still valid, round-trips),
//! completeness failures (validated claim missing regime / dataset, verified
//! claim with a bad interval), Merkle content-addressing (determinism + tamper
//! detection), the format-version gate, optional signature, the color
//! breakdown, and deterministic JSON.

use fs_evidence::{Color, ValidityDomain};
use fs_package::{Claim, EvidencePackage, PackageError, Provenance};

fn prov() -> Provenance {
    Provenance::new("commit-abc123", "lock-deadbeef")
}
fn verified(id: &str) -> Claim {
    Claim::new(
        id,
        format!("{id}: stress <= sigma*"),
        Color::Verified { lo: -1.0, hi: 1.0 },
    )
}
fn estimated(id: &str) -> Claim {
    Claim::new(
        id,
        format!("{id}: surrogate says ok"),
        Color::Estimated {
            estimator: "surrogate".into(),
            dispersion: 2.0,
        },
    )
}
fn validated(id: &str, regime: ValidityDomain, dataset: &str) -> Claim {
    Claim::new(
        id,
        format!("{id}: matches data"),
        Color::Validated {
            regime,
            dataset: dataset.into(),
        },
    )
}
fn good_regime() -> ValidityDomain {
    ValidityDomain::unconstrained().with("Re", 1e5, 3e5)
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
    assert_ne!(report.merkle_root, 0);
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
    let pkg = EvidencePackage::new(prov()).with_claim(Claim::new(
        "v",
        "backwards",
        Color::Verified { lo: 5.0, hi: 1.0 },
    ));
    assert!(matches!(
        pkg.verify(),
        Err(PackageError::IncompleteVerifiedClaim { .. })
    ));
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
    // tampering with a claim changes the root.
    let tampered = EvidencePackage::new(prov())
        .with_claim(Claim::new(
            "c1",
            "TAMPERED",
            Color::Verified { lo: -1.0, hi: 1.0 },
        ))
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
    assert!(j.contains(&format!("{:016x}", pkg.merkle_root())));
    assert!(j.contains("\"format_version\":2"), "schema v2 (qmao.6.1)");
    // v2 carries COMPLETE payloads, not just rank labels.
    assert!(j.contains("\"lo_bits\":") && j.contains("\"dataset\":"));
}

/// qmao.6.1 — the schema-v2 round trip and its fail-closed walls: a
/// package decode-encodes stably (golden), floats travel bit-exactly,
/// and hostile/missing/unknown/non-finite/forged/tampered inputs each
/// refuse with a structured parse error.
#[test]
fn v2_round_trip_and_fail_closed_walls() {
    let pkg = EvidencePackage::new(Provenance::new("frankensim@abc123", "lock:deadbeef"))
        .with_claim(Claim::new(
            "c-verified",
            "tip deflection within bound",
            Color::Verified {
                lo: 0.1875,
                hi: 0.25,
            },
        ))
        .with_claim(Claim::new(
            "c-validated",
            "k-epsilon matched within regime",
            Color::Validated {
                regime: fs_evidence::ValidityDomain::unconstrained().with("reynolds", 1e3, 1e5),
                dataset: "tunnel-run-9".to_string(),
            },
        ))
        .with_claim(Claim::new(
            "c-estimated",
            "surrogate prediction",
            Color::Estimated {
                estimator: "pod-deim".to_string(),
                dispersion: 0.02,
            },
        ))
        .signed("test-key/1234abcd");
    // Golden decode-encode stability: parse(to_json) == pkg, and the
    // re-emission is byte-identical (semantic AND textual round trip).
    let json = pkg.to_json();
    let back = EvidencePackage::from_json(&json).expect("canonical JSON parses");
    assert_eq!(back, pkg, "semantic round trip");
    assert_eq!(back.to_json(), json, "textual round trip");
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
}

/// qmao.6.1 — the magnitude budget attributes ERROR MAGNITUDES, not
/// claim counts, and reconciles with an independent recomputation.
#[test]
fn magnitude_budget_reconciles() {
    let pkg = EvidencePackage::new(Provenance::new("v", "l"))
        .with_claim(Claim::new("a", "s", Color::Verified { lo: 0.0, hi: 0.5 }))
        .with_claim(Claim::new("b", "s", Color::Verified { lo: 1.0, hi: 1.25 }))
        .with_claim(Claim::new(
            "c",
            "s",
            Color::Estimated {
                estimator: "e".to_string(),
                dispersion: 0.125,
            },
        ))
        .with_claim(Claim::new(
            "d",
            "s",
            Color::Validated {
                regime: fs_evidence::ValidityDomain::unconstrained().with("re", 1.0, 2.0),
                dataset: "ds".to_string(),
            },
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
/// "mapped but absent", never covered; unrepresentable concepts
/// (falsifier logs in schema v2) can never report covered.
#[test]
fn coverage_cannot_claim_absent_evidence() {
    use fs_crosswalk::{PackageConcept, Standard};
    use fs_evidence::Color;
    use fs_package::{Claim, EvidencePackage, Provenance};
    use fs_package::{CoverageStatus, package_coverage, package_presence};
    // Unsigned, verified-only package: no validated claims, no regime
    // tags, no datasets, no signature.
    let pkg = EvidencePackage::new(Provenance::new("v", "lock")).with_claim(Claim::new(
        "c",
        "bounded",
        Color::Verified { lo: 0.0, hi: 1.0 },
    ));
    let presence = package_presence(&pkg);
    let get = |c: PackageConcept| presence.iter().find(|p| p.concept == c).unwrap();
    assert!(get(PackageConcept::VerifiedColor).present);
    assert!(!get(PackageConcept::ValidatedColor).present);
    assert!(!get(PackageConcept::RegimeTag).present);
    assert!(!get(PackageConcept::AnchoringDataset).present);
    assert!(!get(PackageConcept::Signature).present);
    assert!(
        !get(PackageConcept::FalsifierLog).present,
        "unrepresentable evidence can never read as present"
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
                    "falsifier logs cannot be covered in v2 ({standard:?})"
                );
            }
        }
    }
}
