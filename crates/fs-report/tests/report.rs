//! Battery for automatic lab notebooks + semantic diffs (fs-report). Covers
//! deterministic Markdown rendering (units on every metric), the reproducibility
//! loop (content-addressed + the reproducing IR), and semantic design-diff
//! attribution recovering known edits ranked by significance.

use std::collections::BTreeMap;

use fs_evidence::{Ambition, Color, ModelCard, ValidityDomain};
use fs_package::{Claim, EvidencePackage, Provenance};
use fs_regime::{OperatingPoint, OverrideAcknowledgement, QoiClaim, audit_product_output};
use fs_report::{
    LabNotebook, Quantity, REGIME_DEMOTION_PACKAGE_ESTIMATOR, RegimePackageError, ReproStep,
    project_regime_audit_outputs, regime_demotion_package_claim_id, regime_no_claims_markdown,
    retain_regime_demotions_in_package, semantic_diff,
};

fn study() -> LabNotebook {
    let mut nb = LabNotebook::new("Bracket study", 42, "0.1.0");
    nb.prose("Optimized the bracket for mass under a stiffness floor.")
        .metric("mass", 1.4, "kg")
        .metric("max_stress", 180.0, "MPa")
        .step("optimize", vec!["lbfgs".into(), "50".into()])
        .step("verify", vec!["stiffness".into()]);
    nb
}

#[test]
fn the_notebook_renders_all_sections_with_units() {
    let md = study().render_markdown();
    assert!(md.contains("# Bracket study"));
    assert!(md.contains("seed: 42") && md.contains("version: 0.1.0")); // provenance
    assert!(md.contains("Optimized the bracket"));
    // units on every value (P10).
    assert!(
        md.contains("**mass**: 1.4 kg"),
        "missing unit-labelled metric:\n{md}"
    );
    assert!(md.contains("**max_stress**: 180 MPa"));
    assert!(md.contains("repro: `optimize(lbfgs, 50)`"));
}

#[test]
fn metrics_carry_their_units() {
    let nb = study();
    let metrics = nb.metrics();
    assert_eq!(metrics.len(), 2);
    assert_eq!(metrics[0], ("mass", &Quantity::new(1.4, "kg")));
}

#[test]
fn the_notebook_carries_the_exact_reproducing_ir() {
    let ir = study().repro_ir();
    assert_eq!(
        ir,
        vec![
            ReproStep {
                op: "optimize".into(),
                args: vec!["lbfgs".into(), "50".into()]
            },
            ReproStep {
                op: "verify".into(),
                args: vec!["stiffness".into()]
            },
        ]
    );
}

#[test]
fn the_reproducibility_loop_closes_by_content_hash() {
    // rebuilding the study from the same inputs reproduces the exact artifact.
    let h1 = study().content_hash();
    let h2 = study().content_hash();
    assert_eq!(h1, h2);
    // a changed metric changes the content hash (no silent drift).
    let mut altered = LabNotebook::new("Bracket study", 42, "0.1.0");
    altered
        .prose("Optimized the bracket for mass under a stiffness floor.")
        .metric("mass", 1.5, "kg"); // 1.4 -> 1.5
    assert_ne!(altered.content_hash(), h1);
    // gp3.14 regression: a Prose block whose body embeds a rendered
    // step line is not the same notebook as an actual Step + Prose —
    // but both RENDER byte-identically. The former markdown-render
    // hash collided here; the structural hash refuses.
    let mut real = LabNotebook::new("x", 1, "0.1.0");
    real.step("x", Vec::new()).prose("q");
    let mut forged = LabNotebook::new("x", 1, "0.1.0");
    forged.prose("- repro: `x()`\nq");
    assert_eq!(
        forged.render_markdown(),
        real.render_markdown(),
        "the adversarial pair must render identically for the gate to mean anything"
    );
    assert_ne!(
        forged.content_hash(),
        real.content_hash(),
        "prose imitating a step must not share the content address"
    );
}

#[test]
fn semantic_diff_recovers_known_edits() {
    let before = BTreeMap::from([
        ("wall_thickness".to_string(), Quantity::new(2.0, "mm")),
        ("lip_curvature".to_string(), Quantity::new(1.0, "1/mm")),
        ("mass".to_string(), Quantity::new(1.4, "kg")),
    ]);
    let after = BTreeMap::from([
        ("wall_thickness".to_string(), Quantity::new(1.6, "mm")), // thinned 0.4 mm (-20%)
        ("lip_curvature".to_string(), Quantity::new(0.82, "1/mm")), // -18%
        ("mass".to_string(), Quantity::new(1.4, "kg")),           // unchanged
    ]);
    let d = semantic_diff(&before, &after);
    assert_eq!(d.len(), 3);
    // ranked by significance: wall_thickness (-20%) before lip_curvature (-18%).
    assert_eq!(d[0].name, "wall_thickness");
    assert!((d[0].abs_change - (-0.4)).abs() < 1e-12);
    assert!((d[0].rel_change - (-0.2)).abs() < 1e-12);
    assert_eq!(d[1].name, "lip_curvature");
    assert!((d[1].rel_change - (-0.18)).abs() < 1e-12);
    // the unchanged feature sorts last.
    assert_eq!(d[2].name, "mass");
    assert!(d[2].abs_change.abs() < 1e-12);
    // the attribution string carries units + the percentage.
    assert!(d[0].describe().contains("mm") && d[0].describe().contains("-20.0%"));
}

#[test]
fn reporting_is_deterministic() {
    assert_eq!(study().content_hash(), study().content_hash());
    assert_eq!(study().render_markdown(), study().render_markdown());
}

fn regime_card() -> ModelCard {
    ModelCard::new(
        "forced-convection",
        "2.1.0",
        Ambition::Solid,
        vec![],
        ValidityDomain::unconstrained().with("Re", 10.0, 100.0),
        vec![],
        0.05,
    )
}

fn regime_claim(qoi: &str, acknowledgement: bool) -> QoiClaim {
    QoiClaim {
        qoi: qoi.to_string(),
        color: Color::Validated {
            regime: ValidityDomain::unconstrained().with("Re", 10.0, 100.0),
            dataset: "forced-convection-reference-v1".to_string(),
        },
        model_cards: vec!["forced-convection".to_string()],
        override_acknowledgement: acknowledgement.then(|| OverrideAcknowledgement {
            actor: "reviewer-7".to_string(),
            reason: "exploratory-only".to_string(),
        }),
    }
}

#[test]
fn report_renders_exact_regime_receipts_in_the_no_claim_section() {
    let audit = audit_product_output(
        &[regime_card()],
        &[
            OperatingPoint {
                id: "inside".to_string(),
                groups: BTreeMap::from([("Re".to_string(), 50.0)]),
            },
            OperatingPoint {
                id: "outside".to_string(),
                groups: BTreeMap::from([("Re".to_string(), 1_000.0)]),
            },
        ],
        &[
            regime_claim("temperature:max", true),
            regime_claim("temperature:mean", false),
        ],
    )
    .expect("valid final-envelope audit");
    let markdown = regime_no_claims_markdown(&audit).expect("demotions render");

    for receipt in &audit.receipts {
        assert!(markdown.contains(&receipt.content_id().to_string()));
        assert!(markdown.contains(&receipt.to_canonical_json()));
    }
    for expected in [
        "## Operating-envelope no-claim boundaries",
        "estimated / no dispersion claim",
        "coverage `partial`",
        "`outside` / `forced-convection` / `Re`",
        "Override acknowledged by `reviewer-7`",
        "acknowledgement does not restore color",
        "cannot authenticate model-card or calibration authorities",
    ] {
        assert!(
            markdown.contains(expected),
            "missing {expected:?}:\n{markdown}"
        );
    }
    assert_eq!(markdown, regime_no_claims_markdown(&audit).unwrap());
}

#[test]
fn fully_in_domain_audit_does_not_invent_a_no_claim_section() {
    let audit = audit_product_output(
        &[regime_card()],
        &[OperatingPoint {
            id: "inside".to_string(),
            groups: BTreeMap::from([("Re".to_string(), 50.0)]),
        }],
        &[regime_claim("temperature:max", false)],
    )
    .expect("valid final-envelope audit");

    assert_eq!(regime_no_claims_markdown(&audit), None);
}

#[test]
fn demotion_receipts_round_trip_in_packages_without_certificate_laundering() {
    let audit = audit_product_output(
        &[regime_card()],
        &[
            OperatingPoint {
                id: "inside".to_string(),
                groups: BTreeMap::from([("Re".to_string(), 50.0)]),
            },
            OperatingPoint {
                id: "outside".to_string(),
                groups: BTreeMap::from([("Re".to_string(), 1_000.0)]),
            },
        ],
        &[
            regime_claim("temperature:max", true),
            regime_claim("temperature:mean", false),
        ],
    )
    .expect("valid final-envelope audit");
    let package = retain_regime_demotions_in_package(
        EvidencePackage::new(Provenance::new("regime-package-test", "Cargo.lock:test")),
        &audit,
    )
    .expect("demotion receipts package");

    let demoted = audit
        .receipts
        .iter()
        .filter(|receipt| receipt.demoted())
        .collect::<Vec<_>>();
    assert_eq!(package.declared_claims_unverified().len(), demoted.len());
    for receipt in demoted {
        let claim = package
            .declared_claims_unverified()
            .iter()
            .find(|claim| claim.id() == regime_demotion_package_claim_id(&audit, receipt))
            .expect("receipt claim is present");
        assert!(matches!(
            claim.declared_color_unverified(),
            Color::Estimated {
                estimator,
                dispersion,
            } if estimator == REGIME_DEMOTION_PACKAGE_ESTIMATOR && dispersion.is_infinite()
        ));
        assert_eq!(claim.declared_semantic_witness_unverified(), None);
        assert_eq!(
            claim.statement(),
            format!(
                "{{\"schema\":\"fs-regime-output-demotion-package-v1\",\"audit_provenance\":\"{:016x}\",\"receipt_content_id\":\"{}\",\"receipt\":{}}}",
                audit.provenance.0,
                receipt.content_id(),
                receipt.to_canonical_json()
            )
        );
    }

    let json = package.to_json().expect("bounded package serializes");
    let decoded = EvidencePackage::from_json(&json).expect("package round trips");
    assert_eq!(decoded, package);
    let idempotent = retain_regime_demotions_in_package(package.clone(), &audit)
        .expect("exact retry is idempotent");
    assert_eq!(idempotent, package);
}

#[test]
fn package_projection_omits_in_domain_receipts_and_refuses_id_conflicts() {
    let in_domain = audit_product_output(
        &[regime_card()],
        &[OperatingPoint {
            id: "inside".to_string(),
            groups: BTreeMap::from([("Re".to_string(), 50.0)]),
        }],
        &[regime_claim("temperature:max", false)],
    )
    .expect("valid final-envelope audit");
    let base = EvidencePackage::new(Provenance::new(
        "regime-package-in-domain-test",
        "Cargo.lock:test",
    ));
    let unchanged = retain_regime_demotions_in_package(base.clone(), &in_domain)
        .expect("fully in-domain audit is a no-op");
    assert_eq!(unchanged, base);

    let mut partial = audit_product_output(
        &[regime_card()],
        &[OperatingPoint {
            id: "outside".to_string(),
            groups: BTreeMap::from([("Re".to_string(), 1_000.0)]),
        }],
        &[regime_claim("temperature:max", false)],
    )
    .expect("valid demotion audit");
    let receipt = partial.receipts.pop().expect("one receipt");
    let conflict = EvidencePackage::new(Provenance::new(
        "regime-package-conflict-test",
        "Cargo.lock:test",
    ))
    .with_claim(Claim::estimated(
        regime_demotion_package_claim_id(&partial, &receipt),
        "different retained bytes",
        REGIME_DEMOTION_PACKAGE_ESTIMATOR,
        f64::INFINITY,
    ));
    partial.receipts.push(receipt);
    assert!(matches!(
        retain_regime_demotions_in_package(conflict, &partial),
        Err(RegimePackageError::ClaimIdConflict { .. })
    ));
}

#[test]
fn coupled_projection_cannot_drop_either_side_of_a_demoted_audit() {
    let audit = audit_product_output(
        &[regime_card()],
        &[
            OperatingPoint {
                id: "inside".to_string(),
                groups: BTreeMap::from([("Re".to_string(), 50.0)]),
            },
            OperatingPoint {
                id: "outside".to_string(),
                groups: BTreeMap::from([("Re".to_string(), 1_000.0)]),
            },
        ],
        &[
            regime_claim("temperature:max", true),
            regime_claim("temperature:mean", false),
        ],
    )
    .expect("valid final-envelope audit");
    let base = EvidencePackage::new(Provenance::new(
        "coupled-regime-projection-test",
        "Cargo.lock:test",
    ));
    let outputs = project_regime_audit_outputs(base, &audit).expect("coupled projection");
    let markdown = outputs
        .no_claims_markdown
        .expect("partial envelope must render no-claim boundaries");

    for receipt in audit.receipts.iter().filter(|receipt| receipt.demoted()) {
        assert!(markdown.contains(&receipt.content_id().to_string()));
        assert!(markdown.contains(&receipt.to_canonical_json()));
        let claim_id = regime_demotion_package_claim_id(&audit, receipt);
        let claim = outputs
            .package
            .declared_claims_unverified()
            .iter()
            .find(|claim| claim.id() == claim_id)
            .expect("the same receipt must be retained in the package");
        assert!(
            claim
                .statement()
                .contains(&receipt.content_id().to_string())
        );
        assert!(claim.statement().contains(&receipt.to_canonical_json()));
    }

    let in_domain = audit_product_output(
        &[regime_card()],
        &[OperatingPoint {
            id: "inside".to_string(),
            groups: BTreeMap::from([("Re".to_string(), 50.0)]),
        }],
        &[regime_claim("temperature:max", false)],
    )
    .expect("valid in-domain audit");
    let base = EvidencePackage::new(Provenance::new(
        "coupled-regime-in-domain-test",
        "Cargo.lock:test",
    ));
    let unchanged = project_regime_audit_outputs(base.clone(), &in_domain)
        .expect("in-domain coupled projection");
    assert_eq!(unchanged.no_claims_markdown, None);
    assert_eq!(unchanged.package, base);
}
