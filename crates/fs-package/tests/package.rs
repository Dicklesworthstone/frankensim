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
    assert!(j.contains("\"format_version\":1"));
}
