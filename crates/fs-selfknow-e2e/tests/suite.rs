//! The STRUCTURE & SELF-KNOWLEDGE battery (bead knh1.7; runs under
//! `selfknow-e2e`): six stages, one script, REAL fs-ledger events per
//! stage, and every stage's FAIL-SAFE assertion — unknown pairings
//! rejected, no false speedups, no crying wolf, no infinite descent,
//! no confabulated explanations, no useless VoI purchases.
//! Gauntlet: G0 (type/VoI laws), G2 (drag + VoI audit), G3
//! (explanation reconciliation, symmetry falsifier), G4 (degraded-gap
//! refusal).
#![cfg(feature = "selfknow-e2e")]

use fs_evidence::Color;
use fs_ledger::{EventRow, Ledger};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-selfknow-e2e\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

struct Harness {
    ledger: Ledger,
    t: std::cell::Cell<i64>,
}

impl Harness {
    fn new(tag: &str) -> Harness {
        let dir = std::env::temp_dir().join(format!("selfknow-e2e-{}-{tag}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("tempdir");
        Harness {
            ledger: Ledger::open(dir.join("suite.led").to_str().expect("utf8")).expect("ledger"),
            t: std::cell::Cell::new(0),
        }
    }

    fn log(&self, kind: &str, payload: &str) {
        self.t.set(self.t.get() + 1);
        self.ledger
            .append_event(&EventRow {
                session: None,
                t: self.t.get(),
                kind,
                payload: Some(payload),
            })
            .expect("event");
    }
}

#[test]
fn stage_1_interface_types() {
    use fs_iface::{CouplingGraph, CouplingRole, PairingRegistry, SpaceType, check};
    let h = Harness::new("iface");
    let registry = PairingRegistry::standard();
    // LEGAL: the Stokes-class certified pairing passes.
    let legal = CouplingGraph::new()
        .field("velocity", SpaceType::HDiv)
        .field("pressure", SpaceType::L2)
        .couple("stokes", "velocity", "pressure", CouplingRole::Saddle);
    let report = check(&legal, &registry);
    h.log(
        "iface.legal",
        &format!("{{\"admitted\":{}}}", report.admitted),
    );
    assert!(report.admitted, "the certified pairing is admitted");
    // ILLEGAL: an unstable pairing is rejected PRE-RUN with a located
    // diagnostic.
    let illegal = CouplingGraph::new()
        .field("velocity", SpaceType::HGrad)
        .field("pressure", SpaceType::HGrad)
        .couple("equal-order", "velocity", "pressure", CouplingRole::Saddle);
    let report = check(&illegal, &registry);
    h.log(
        "iface.illegal",
        &format!(
            "{{\"admitted\":{},\"findings\":{}}}",
            report.admitted,
            report.findings.len()
        ),
    );
    assert!(!report.admitted, "the unstable pairing is rejected");
    assert!(
        report.findings.iter().any(|f| f.coupling == "equal-order"),
        "the diagnostic is LOCALIZED to the offending coupling"
    );
    // UNKNOWN: a pairing absent from the registry is rejected
    // CONSERVATIVELY (illegal-until-certified — the fail-safe).
    let unknown = CouplingGraph::new()
        .field("a", SpaceType::HCurl)
        .field("b", SpaceType::HCurl)
        .couple("mystery", "a", "b", CouplingRole::Saddle);
    let report = check(&unknown, &registry);
    h.log(
        "iface.unknown",
        &format!("{{\"admitted\":{}}}", report.admitted),
    );
    assert!(
        !report.admitted,
        "unknown pairings are illegal until certified"
    );
    verdict(
        "stage-1",
        "certified pairing admitted; unstable pairing rejected with a localized \
         diagnostic; unknown pairing rejected conservatively (G0)",
    );
}

#[test]
fn stage_2_symmetry_harvest() {
    use fs_symmetry::{cyclic_residual, solve_circulant, symmetrized_solve};
    let h = Harness::new("symmetry");
    // A 4-fold circulant system: the harvested solve must agree with
    // the full solve BITWISE-CLASS (the full-solve falsifier).
    let first_row = vec![4.0, -1.0, 0.5, -1.0];
    let rhs_sym = vec![1.0, 2.0, 1.0, 2.0]; // 2-fold symmetric
    let direct = solve_circulant(&first_row, &rhs_sym).expect("direct");
    let harvested = symmetrized_solve(&first_row, &rhs_sym, 2).expect("harvest");
    let gap: f64 = direct
        .iter()
        .zip(&harvested.symmetric_solution)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    h.log(
        "symmetry.exact",
        &format!(
            "{{\"gap\":{gap:.3e},\"residual\":{:.3e}}}",
            harvested.asymmetry_residual
        ),
    );
    assert!(gap < 1e-12, "exact symmetry: harvested == full solve");
    assert!(
        harvested.asymmetry_residual < 1e-12,
        "no asymmetry detected"
    );
    // APPROXIMATE symmetry: the certified perturbation correction
    // CONTAINS the true asymmetry effect.
    let rhs_approx = vec![1.0, 2.0, 1.1, 2.0]; // slightly broken
    let bound = symmetrized_solve(&first_row, &rhs_approx, 2).expect("approx");
    let truth = solve_circulant(&first_row, &rhs_approx).expect("truth");
    let true_gap: f64 = truth
        .iter()
        .zip(&bound.symmetric_solution)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    h.log(
        "symmetry.approx",
        &format!(
            "{{\"true_gap\":{true_gap:.3e},\"bound\":{:.3e}}}",
            bound.correction_bound
        ),
    );
    assert!(
        true_gap <= bound.correction_bound + 1e-12,
        "the certified correction bound contains the true asymmetry: {true_gap:.2e} \
         <= {:.2e}",
        bound.correction_bound
    );
    // FULLY ASYMMETRIC input: the residual diagnostic reports a large
    // asymmetric fraction — no false speedup is on offer.
    let res = cyclic_residual(&[3.0, -1.0, 0.2, 5.0], 2).expect("residual");
    h.log(
        "symmetry.asym",
        &format!("{{\"fraction\":{:.3}}}", res.relative),
    );
    assert!(
        res.relative > 0.3,
        "asymmetric input is flagged, not falsely harvested: {}",
        res.relative
    );
    verdict(
        "stage-2",
        "harvested solve matches the full-solve falsifier to 1e-12; the perturbation \
         bound contains the true asymmetry; a broken-symmetry input is flagged (G3)",
    );
}

#[test]
fn stage_3_spectral_health() {
    use fs_spectral::{GapHealthMonitor, Health, propagate, spectral_gap};
    let h = Harness::new("spectral");
    // DEGRADED gap: the triage must REFUSE confidence (refusal is the
    // pass) and the flag must survive color propagation into a merge.
    let degraded = [0.001, 0.0012, 0.9, 1.1]; // collapsed leading gap
    let gap = spectral_gap(&degraded).expect("gap");
    let mut monitor = GapHealthMonitor::new(0.2, 0.5);
    let health = monitor.update(gap.ratio);
    h.log(
        "spectral.degraded",
        &format!("{{\"ratio\":{:.4},\"health\":\"{health:?}\"}}", gap.ratio),
    );
    assert_eq!(health, Health::Degraded, "degraded gap refuses confidence");
    let merged = propagate(Color::Verified { lo: 0.0, hi: 1.0 }, health);
    assert!(
        !matches!(merged, Color::Verified { .. }),
        "the low-confidence flag DEMOTES the merge color: {merged:?}"
    );
    // WELL-CONDITIONED assembly: no crying wolf.
    let healthy = [1.0, 5.0, 5.5, 6.0];
    let gap = spectral_gap(&healthy).expect("gap");
    let health = monitor.update(gap.ratio);
    h.log(
        "spectral.healthy",
        &format!("{{\"ratio\":{:.4},\"health\":\"{health:?}\"}}", gap.ratio),
    );
    assert_eq!(health, Health::Healthy, "a clean gap is not flagged");
    let kept = propagate(Color::Verified { lo: 0.0, hi: 1.0 }, health);
    assert!(
        matches!(kept, Color::Verified { .. }),
        "healthy assemblies keep their color"
    );
    verdict(
        "stage-3",
        "collapsed lambda-gap refuses confidence and demotes the merge color; a clean \
         assembly is not flagged — refusal without crying wolf (G4)",
    );
}

#[test]
fn stage_4_abstraction_ladder() {
    use fs_surrogate::ladder::Ladder;
    let h = Harness::new("ladder");
    let ladder = Ladder::build(150, (0.0, 4.0), &[5, 2], true).expect("bounded ladder plan");
    // The leak alarm fires and auto-drills; termination at full order
    // is guaranteed (no infinite descent — the fail-safe).
    let ans = ladder
        .at_level(ladder.top())
        .expect("declared top level")
        .query(1.7, 1e-32)
        .expect("in-range ladder query");
    h.log(
        "ladder.drill",
        &format!(
            "{{\"level_used\":{},\"leaks\":{:?}}}",
            ans.level_used(),
            ans.leaks()
        ),
    );
    assert_eq!(
        ans.level_used(),
        0,
        "a query that leaks everywhere ends at truth"
    );
    assert_eq!(
        ans.leaks().len(),
        3,
        "every rung above it is a recorded leak"
    );
    // A satisfiable query stays high; an estimate-tolerant query gets
    // the concept rung with ESTIMATED color (never Verified).
    let quiet = ladder
        .at_level(1)
        .expect("declared RB level")
        .query(1.7, 1e-2)
        .expect("in-range ladder query");
    assert_eq!(quiet.level_used(), 1);
    assert!(quiet.leaks().is_empty());
    let concept = ladder
        .at_level(ladder.top())
        .expect("declared top level")
        .query(1.7, 0.5)
        .expect("in-range ladder query");
    assert_eq!(
        concept.level_used(),
        ladder.top(),
        "the estimate-tolerant control must actually remain at the concept rung"
    );
    h.log(
        "ladder.concept",
        &format!(
            "{{\"color_estimated\":{}}}",
            matches!(concept.color(), Color::Estimated { .. })
        ),
    );
    assert!(matches!(
        concept.color(),
        Color::Estimated { estimator, .. }
            if estimator == "fs-surrogate.concept-cross-rung-v1"
    ));
    verdict(
        "stage-4",
        "the leak alarm drills to full order with a complete trail (no infinite \
         descent); satisfiable queries stay at their rung; the concept rung answers \
         with honest estimated color",
    );
}

#[test]
fn stage_5_explanation_objects() {
    use fs_adjoint::explain::{Elliptic1d, Explanation, adjoint_attribution, finalize};
    let h = Harness::new("explain");
    let fixture = Elliptic1d::new(100).expect("bounded explanation fixture");
    let a0 = vec![1.0f64; 101];
    let mut a1 = a0.clone();
    for (e, ae) in a1.iter_mut().enumerate() {
        *ae = if e < 50 { 1.3 } else { 0.8 };
    }
    let observed = fixture
        .compliance(&fixture.solve(&a1).expect("positive conductivity solves"))
        .expect("finite compliance")
        - fixture
            .compliance(&fixture.solve(&a0).expect("positive conductivity solves"))
            .expect("finite compliance");
    // Full channel set: reconciles within bounds (G3 — the permanent
    // invariant).
    let full = [
        ("left-half", (0..50).collect::<Vec<_>>()),
        ("right-half", (50..=100).collect::<Vec<_>>()),
    ];
    let ok = finalize(
        adjoint_attribution(&fixture, &a0, &a1, &full).expect("valid full attribution"),
        observed,
        1e-8,
    )
    .expect("valid explanation finalization");
    h.log(
        "explain.reconcile",
        &format!("{{\"ok\":{}}}", ok.reconciles()),
    );
    assert!(matches!(ok, Explanation::Explained { .. }) && ok.reconciles());
    // Hidden channel: the honesty gate refuses AND the rendering layer
    // refuses too (no confabulated explanation).
    let partial = [("left-half", (0..50).collect::<Vec<_>>())];
    let refused = finalize(
        adjoint_attribution(&fixture, &a0, &a1, &partial).expect("valid partial attribution"),
        observed,
        1e-8,
    )
    .expect("valid explanation finalization");
    let narrative = refused.render_narrative();
    h.log(
        "explain.refuse",
        &format!(
            "{{\"refused\":{}}}",
            matches!(refused, Explanation::Refused { .. })
        ),
    );
    assert!(matches!(refused, Explanation::Refused { .. }));
    assert!(
        narrative.contains("REFUSED") && narrative.contains("NON-AUTHORITATIVE"),
        "the NL layer refuses alongside the gate"
    );
    verdict(
        "stage-5",
        "channels + residual reconcile to the observed change; a hidden channel makes \
         BOTH the gate and the narrative refuse — no confabulation (G3)",
    );
}

#[test]
fn stage_6_value_of_information() {
    use fs_plan::voi::{
        LiveDecision, Probe, ProbeKind, UncertaintyNode, hint_for_query, rank_purchases,
    };
    let h = Harness::new("voi");
    let margin = |v: &[f64]| 2.0 * v[0] - 1.0;
    let decision = LiveDecision {
        margin: &margin,
        arity: 1,
    };
    // A decision that CAN flip: a priced ranking comes back.
    let live = vec![UncertaintyNode {
        name: "drag".to_string(),
        lo: 0.0,
        hi: 1.0,
        nominal: 0.6,
    }];
    let menu = vec![Probe {
        name: "anchor".to_string(),
        target: "drag".to_string(),
        cost: 20.0,
        shrink: 0.1,
        kind: ProbeKind::Physical,
    }];
    let ranked = rank_purchases(&decision, &live, &menu, 64).expect("valid live VoI request");
    h.log(
        "voi.live",
        &format!("{{\"top_score\":{:.5}}}", ranked[0].score()),
    );
    assert!(
        ranked[0].score() > 0.0,
        "a flippable decision prices its evidence"
    );
    assert!(
        hint_for_query(&ranked)
            .expect("valid live hint")
            .contains("anchor")
    );
    // A decision that CANNOT flip: the recommendation is EMPTY (no
    // useless purchase — the fail-safe).
    let settled = vec![UncertaintyNode {
        name: "drag".to_string(),
        lo: 0.9,
        hi: 1.0,
        nominal: 0.95,
    }];
    let ranked = rank_purchases(&decision, &settled, &menu, 64).expect("valid settled VoI request");
    h.log(
        "voi.settled",
        &format!("{{\"top_score\":{:.5}}}", ranked[0].score()),
    );
    assert_eq!(
        ranked[0].score().to_bits(),
        0.0f64.to_bits(),
        "an unflippable decision buys nothing"
    );
    assert!(
        hint_for_query(&ranked)
            .expect("valid settled hint")
            .contains("spend nothing"),
        "the hint says so explicitly"
    );
    verdict(
        "stage-6",
        "a flippable decision returns priced top-k purchases; an unflippable one \
         returns 'spend nothing' — no useless VoI purchase (G0/G2)",
    );
}
