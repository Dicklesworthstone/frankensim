//! Battery for the standalone evidence-package checker (addendum Proposal 12).
//! Covers a clean pass, completeness-failure findings, content-address
//! (Merkle) tamper detection, signature-presence reporting, budget-pie
//! rendering (including the empty case), the protocol version, and
//! determinism. The checker uses only the package format — no solver.

use fs_checker::{
    CHECKER_PROTOCOL_VERSION, ColorBreakdown, ContentHash, SignatureStatus, Verdict, check,
    check_against_root, check_for_release, check_json_for_release,
};
use fs_evidence::{Color, ValidityDomain};
use fs_package::{Claim, EvidencePackage, FalsifierRecord, Provenance};

/// A deliberately corrupted content root (one byte flipped): the v4
/// 32-byte replacement for the old `root ^ 0xdead` tamper idiom.
fn flip(root: ContentHash) -> ContentHash {
    let mut bytes = *root.as_bytes();
    bytes[0] ^= 0xde;
    ContentHash(bytes)
}

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
fn a_semantically_empty_falsifier_record_fails_the_check() {
    let pkg =
        EvidencePackage::new(prov()).with_claim(verified("v").with_falsifier(FalsifierRecord {
            name: " ".to_string(),
            attempts: 0,
            refuted: false,
            detail: " ".to_string(),
        }));
    let report = check(&pkg);
    assert!(!report.passed());
    assert_eq!(report.findings[0].kind, "invalid-falsifier-record");
    assert_eq!(report.breakdown, ColorBreakdown::default());
    assert_eq!(report.render_pie(), "budget pie: no claims");
}

#[test]
fn placeholder_claim_and_falsifier_text_fail_the_check() {
    let placeholder_statement = EvidencePackage::new(prov()).with_claim(Claim::new(
        "claim",
        "TODO",
        Color::Verified { lo: 0.0, hi: 1.0 },
    ));
    let report = check(&placeholder_statement);
    assert!(!report.passed());
    assert_eq!(report.findings[0].kind, "invalid-claim-statement");

    let placeholder_falsifier = EvidencePackage::new(prov()).with_claim(
        verified("claim").with_falsifier(FalsifierRecord {
            name: "independent-probe".to_string(),
            attempts: 1,
            refuted: false,
            detail: "placeholder".to_string(),
        }),
    );
    let report = check(&placeholder_falsifier);
    assert!(!report.passed());
    assert_eq!(report.findings[0].kind, "invalid-falsifier-record");
}

#[test]
fn content_address_mismatch_is_caught() {
    let pkg = EvidencePackage::new(prov()).with_claim(verified("c1"));
    let real_root = pkg.merkle_root();
    // the right root passes.
    assert!(check_against_root(&pkg, real_root).passed());
    // a wrong expected root (a tampered/substituted package) fails.
    let report = check_against_root(&pkg, flip(real_root));
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
        SignatureStatus::Unverified("ed25519:cafe".to_string())
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

struct ReleaseVerifier;

impl fs_checker::SignatureVerifier for ReleaseVerifier {
    fn verify(&self, merkle_root: &ContentHash, signature: &str) -> bool {
        signature == format!("release-test:{merkle_root}")
    }
}

fn signed_for_release(pkg: EvidencePackage) -> EvidencePackage {
    let root = pkg.merkle_root();
    pkg.signed(format!("release-test:{root}"))
}

fn passed_falsifier() -> FalsifierRecord {
    FalsifierRecord {
        name: "independent-interval-probe".to_string(),
        attempts: 64,
        refuted: false,
        detail: "64 boundary-biased probes found no violation".to_string(),
    }
}

const DATASET_HASH: &str = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";

#[test]
fn release_gate_requires_certificate_obligations() {
    let pkg = signed_for_release(
        EvidencePackage::new(prov())
            .with_claim(verified("verified").with_falsifier(passed_falsifier()))
            .with_claim(
                validated("validated", good_regime())
                    .with_falsifier(passed_falsifier())
                    .with_anchor("wt-2026", DATASET_HASH),
            )
            .with_claim(estimated("honest-estimate")),
    );
    let root = pkg.merkle_root();
    let report = check_for_release(&pkg, root, &ReleaseVerifier);
    assert!(report.passed(), "{:?}", report.findings);
    assert!(matches!(report.signature, SignatureStatus::Valid(_)));
    assert!(check_json_for_release(&pkg.to_json(), root, &ReleaseVerifier).passed());
}

#[test]
fn release_gate_refuses_vacuous_or_unpaired_packages() {
    let empty = signed_for_release(EvidencePackage::new(prov()));
    assert!(
        check(&empty).passed(),
        "ordinary integrity check stays vacuous"
    );
    let report = check_for_release(&empty, empty.merkle_root(), &ReleaseVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.kind == "release-empty-package")
    );

    let unpaired = signed_for_release(EvidencePackage::new(prov()).with_claim(verified("v")));
    let report = check_for_release(&unpaired, unpaired.merkle_root(), &ReleaseVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.kind == "release-falsifier-required")
    );
    let report = check_json_for_release(
        &unpaired.to_json(),
        unpaired.merkle_root(),
        &ReleaseVerifier,
    );
    assert!(!report.passed(), "JSON must not bypass release policy");
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.kind == "release-falsifier-required")
    );
}

#[test]
fn release_gate_requires_matching_anchor_signature_and_root() {
    let unanchored = signed_for_release(
        EvidencePackage::new(prov()).with_claim(
            validated("v", good_regime())
                .with_falsifier(passed_falsifier())
                .with_anchor("different-dataset", DATASET_HASH),
        ),
    );
    let report = check_for_release(&unanchored, unanchored.merkle_root(), &ReleaseVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.kind == "release-anchor-required")
    );

    let unsigned =
        EvidencePackage::new(prov()).with_claim(verified("v").with_falsifier(passed_falsifier()));
    let report = check_for_release(&unsigned, unsigned.merkle_root(), &ReleaseVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.kind == "release-signature-required")
    );

    let signed = signed_for_release(unsigned);
    let report = check_for_release(&signed, flip(signed.merkle_root()), &ReleaseVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings
            .iter()
            .any(|finding| finding.kind == "content-address-mismatch")
    );
}

#[test]
fn the_checker_advertises_its_protocol_version() {
    assert_eq!(CHECKER_PROTOCOL_VERSION, 2);
    assert_eq!(fs_checker::CHECKER_SUPPORTED_PACKAGE_FORMAT, 4);
    assert_eq!(
        fs_checker::CHECKER_SUPPORTED_PACKAGE_FORMAT,
        fs_package::FORMAT_VERSION
    );
}

#[test]
fn checking_is_deterministic() {
    let pkg = EvidencePackage::new(prov())
        .with_claim(verified("c1"))
        .with_claim(estimated("c2"));
    assert_eq!(check(&pkg), check(&pkg));
}

/// qmao.6.1 — the third-party JSON path: parse-refused inputs never
/// pass; signature validity is asserted only through a capability over
/// the recomputed root; tamper anywhere fails.
#[test]
fn checker_json_path_and_signature_capability() {
    use fs_checker::{NoSignatureVerifier, SignatureVerifier, check_json, check_with};
    use fs_evidence::Color;
    struct MacVerifier;
    fn mac(root: &ContentHash) -> String {
        format!("test-key/{root}")
    }
    impl SignatureVerifier for MacVerifier {
        fn verify(&self, merkle_root: &ContentHash, signature: &str) -> bool {
            signature == mac(merkle_root)
        }
    }
    let base = EvidencePackage::new(Provenance::new("v1.0", "lock:abc")).with_claim(Claim::new(
        "c1",
        "bounded",
        Color::Verified { lo: 0.0, hi: 1.0 },
    ));
    let root = base.merkle_root();
    let pkg = base.signed(mac(&root));
    // Valid signature via the capability.
    let report = check_with(&pkg, Some(root), &MacVerifier);
    assert!(report.passed(), "{:?}", report.findings);
    assert!(matches!(
        report.signature,
        fs_checker::SignatureStatus::Valid(_)
    ));
    // The no-crypto default cannot assert validity (fail closed: the
    // signature stays Unverified and a finding is raised).
    let report = check_with(&pkg, Some(root), &NoSignatureVerifier);
    assert!(!report.passed());
    assert!(
        report
            .findings
            .iter()
            .any(|f| f.kind == "signature-invalid")
    );
    // Full JSON path: round trip passes; tampered JSON is parse-refused
    // (never a Pass with quietly wrong content).
    let json = pkg.to_json();
    assert!(check_json(&json, Some(root), Some(&MacVerifier)).passed());
    let tampered = json.replace("bounded", "PROVEN");
    let report = check_json(&tampered, Some(root), Some(&MacVerifier));
    assert!(!report.passed());
    assert_eq!(report.findings[0].kind, "parse-refused");
}
