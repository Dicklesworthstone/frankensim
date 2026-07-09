//! Battery for the standalone evidence-package checker (addendum Proposal 12).
//! Covers a clean pass, completeness-failure findings, content-address
//! (Merkle) tamper detection, signature-presence reporting, budget-pie
//! rendering (including the empty case), the protocol version, and
//! determinism. The checker uses only the package format — no solver.

use fs_checker::{CHECKER_PROTOCOL_VERSION, SignatureStatus, Verdict, check, check_against_root};
use fs_evidence::{Color, ValidityDomain};
use fs_package::{Claim, EvidencePackage, Provenance};

fn prov() -> Provenance {
    Provenance::new("commit-abc", "lock-def")
}
fn verified(id: &str) -> Claim {
    Claim::new(id, "ok", Color::Verified { lo: -1.0, hi: 1.0 })
}
fn estimated(id: &str) -> Claim {
    Claim::new(
        id,
        "maybe",
        Color::Estimated {
            estimator: "surrogate".into(),
            dispersion: 2.0,
        },
    )
}
fn validated(id: &str, regime: ValidityDomain) -> Claim {
    Claim::new(
        id,
        "matches",
        Color::Validated {
            regime,
            dataset: "wt-2026".into(),
        },
    )
}
fn good_regime() -> ValidityDomain {
    ValidityDomain::unconstrained().with("Re", 1e5, 3e5)
}

#[test]
fn a_valid_package_passes_with_no_findings() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(validated("c2", good_regime()))
        .with_claim(estimated("c3"));
    let report = check(&pkg);
    assert!(report.passed());
    assert_eq!(report.verdict, Verdict::Pass);
    assert!(report.findings.is_empty());
    assert_eq!(report.merkle_root, pkg.merkle_root());
    assert_eq!(report.breakdown.verified, 1);
    assert_eq!(report.breakdown.validated, 1);
    assert_eq!(report.breakdown.estimated, 1);
}

#[test]
fn an_incomplete_validated_claim_fails_the_check() {
    // unconstrained regime = missing regime tag.
    let pkg =
        EvidencePackage::new(prov()).with_claim(validated("v", ValidityDomain::unconstrained()));
    let report = check(&pkg);
    assert!(!report.passed());
    assert_eq!(report.verdict, Verdict::Fail);
    assert_eq!(report.findings.len(), 1);
    assert_eq!(report.findings[0].kind, "incomplete-validated-claim");
}

#[test]
fn content_address_mismatch_is_caught() {
    let pkg = EvidencePackage::new(prov()).with_claim(verified("c1"));
    let real_root = pkg.merkle_root();
    // the right root passes.
    assert!(check_against_root(&pkg, real_root).passed());
    // a wrong expected root (a tampered/substituted package) fails.
    let report = check_against_root(&pkg, real_root ^ 0xdead_beef);
    assert!(!report.passed());
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.kind == "content-address-mismatch")
    );
}

#[test]
fn content_address_mismatch_catches_provenance_tamper() {
    let pkg = EvidencePackage::new(prov()).with_claim(verified("c1"));
    let root = pkg.merkle_root();
    let tampered =
        EvidencePackage::new(Provenance::new("commit-evil", "lock-def")).with_claim(verified("c1"));

    let report = check_against_root(&tampered, root);
    assert!(!report.passed());
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.kind == "content-address-mismatch")
    );
}

#[test]
fn signature_presence_is_reported() {
    let unsigned = EvidencePackage::new(prov()).with_claim(estimated("e1"));
    assert_eq!(check(&unsigned).signature, SignatureStatus::Unsigned);
    let signed = unsigned.signed("ed25519:cafe");
    assert_eq!(
        check(&signed).signature,
        SignatureStatus::Present("ed25519:cafe".to_string())
    );
}

#[test]
fn the_budget_pie_renders_deterministically() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(verified("c2"))
        .with_claim(estimated("c3"));
    let pie = check(&pkg).render_pie();
    assert_eq!(pie, check(&pkg).render_pie());
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

#[test]
fn the_checker_advertises_its_protocol_version() {
    assert_eq!(CHECKER_PROTOCOL_VERSION, 1);
}

#[test]
fn checking_is_deterministic() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(estimated("c2"));
    assert_eq!(check(&pkg), check(&pkg));
}
