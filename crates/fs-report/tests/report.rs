//! Battery for automatic lab notebooks + semantic diffs (fs-report). Covers
//! deterministic Markdown rendering (units on every metric), the reproducibility
//! loop (content-addressed + the reproducing IR), and semantic design-diff
//! attribution recovering known edits ranked by significance.

use std::collections::BTreeMap;

use fs_evidence::{Ambition, Color, ModelCard, ValidityDomain};
use fs_regime::{OperatingPoint, OverrideAcknowledgement, QoiClaim, audit_product_output};
use fs_report::{LabNotebook, Quantity, ReproStep, regime_no_claims_markdown, semantic_diff};

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
