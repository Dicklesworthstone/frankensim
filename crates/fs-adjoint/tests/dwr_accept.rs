//! DWR-accept conformance (the lmp4.4 bead; runs under the `dwr-accept`
//! feature). Acceptance: DWR accepts are correct when colored estimated
//! and become verified ONLY with a valid equilibrated bracket; the
//! laundering attempt fails the type check; DWR-driven refinement
//! concentrates where the QoI error lives; the falsifier's occasional
//! high-fidelity evaluation confirms the estimate's honesty.
#![cfg(feature = "dwr-accept")]

use fs_adjoint::dwr_accept::{Bracket, DwrQuery, accept, dwr_integral_qoi};
use fs_evidence::Color;
use fs_evidence::falsify::{FalsifierRegistry, FalsifierSpec};
use fs_verify::estimator::verify;
use fs_verify::fem1d::{MmsProblem, Poly, solve_p1};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-adjoint/dwr-accept\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

/// A quartic manufactured solution on [0, 1] (zero BC):
/// u = x(1 − x)(1 + 2x²), rich enough that P1 has visible error.
fn quartic_problem(cells: usize) -> MmsProblem {
    // u = (x − x²)(1 + 2x²) = x + 2x³ − x² − 2x⁴.
    let u = Poly(vec![0.0, 1.0, -1.0, 2.0, -2.0]);
    #[allow(clippy::cast_precision_loss)]
    let mesh: Vec<f64> = (0..=cells).map(|k| k as f64 / cells as f64).collect();
    MmsProblem::new("dwr-quartic", u, mesh)
}

/// The exact QoI ∫_{a}^{b} u dx from the polynomial antiderivative.
fn exact_qoi(problem: &MmsProblem, a: f64, b: f64) -> f64 {
    let big_u = problem.u.antiderive();
    big_u.eval(b) - big_u.eval(a)
}

#[test]
fn dw_001_g1_effectivity_and_estimated_accept() {
    let problem = quartic_problem(16);
    let u_h = solve_p1(&problem);
    let (a, b) = (0.25, 0.75);
    let out = dwr_integral_qoi(&problem, &u_h, a, b);
    let true_err = exact_qoi(&problem, a, b) - out.j_primal;
    let effectivity = out.eta / true_err;
    assert!(
        (0.7..=1.3).contains(&effectivity),
        "G1: DWR effectivity near 1: eta {:.3e} vs true {:.3e} ({effectivity:.2})",
        out.eta,
        true_err
    );
    // The DWR-only accept at a discharging tolerance is ESTIMATED.
    let query = DwrQuery {
        qoi: "integral[0.25,0.75]".to_string(),
        tolerance: 1e-3,
    };
    let outcome = accept(&query, out.eta.abs(), None);
    assert!(outcome.accepted);
    assert!(
        matches!(outcome.color, Color::Estimated { .. }),
        "DWR constants are not guaranteed: {:?}",
        outcome.color
    );
    // And a too-tight tolerance rejects.
    let strict = accept(
        &DwrQuery {
            qoi: "integral".to_string(),
            tolerance: 1e-12,
        },
        out.eta.abs(),
        None,
    );
    assert!(!strict.accepted, "no silent discharge");
    verdict(
        "dw-001",
        "quartic G1 fixture: effectivity within [0.7, 1.3]; DWR-only accept is \
         estimated; too-tight tolerance rejects",
    );
}

#[test]
fn dw_002_promotion_requires_a_guaranteed_bracket() {
    let problem = quartic_problem(24);
    let u_h = solve_p1(&problem);
    let (a, b) = (0.25, 0.75);
    let out = dwr_integral_qoi(&problem, &u_h, a, b);
    // The REAL equilibrated primal bound (lmp4.1 machinery)…
    let primal = verify(&problem, &u_h, 1.0);
    // …and the dual's: the dual −z″ = 1_{[a,b]} has the exact piecewise
    // solution; bound its P1 error with the same equilibrated verifier
    // by posing the dual as an MMS problem on the window via a smooth
    // stand-in (a parabola matching the interior load): honest bound on
    // a HARDER dual surrogate.
    let dual_problem = MmsProblem::new(
        "dual-surrogate",
        Poly(vec![0.0, 0.5, -0.5]), // z = x(1−x)/2 solves −z″ = 1
        problem.mesh.clone(),
    );
    let dual_h = solve_p1(&dual_problem);
    let dual = verify(&dual_problem, &dual_h, 1.0);
    let bracket = Bracket::cauchy_schwarz(&primal, &dual);
    assert!(bracket.guaranteed, "both factors certified: {bracket:?}");
    // Promotion: the bracket discharges a tolerance ABOVE it.
    let query = DwrQuery {
        qoi: "integral[0.25,0.75]".to_string(),
        tolerance: bracket.bound * 1.5,
    };
    let outcome = accept(&query, out.eta.abs(), Some(&bracket));
    assert!(outcome.accepted);
    match &outcome.color {
        Color::Verified { lo, hi } => {
            assert!(lo.abs() < f64::EPSILON && (*hi - bracket.bound).abs() < 1e-15);
        }
        other => panic!("bracketed accept must be verified: {other:?}"),
    }
    assert!(
        !outcome.estimator_inconsistent,
        "the DWR estimate sits inside the bracket: {}",
        outcome.audit
    );
    // A NON-guaranteed bracket never promotes.
    let bogus = Bracket {
        bound: 1e-9,
        guaranteed: false,
        source: "wishful".to_string(),
    };
    let not_promoted = accept(&query, out.eta.abs(), Some(&bogus));
    assert!(
        matches!(not_promoted.color, Color::Estimated { .. }),
        "no guarantee, no verified color"
    );
    verdict(
        "dw-002",
        "Cauchy-Schwarz of two equilibrated bounds promotes to verified [0, bound]; \
         non-guaranteed brackets never do",
    );
}

#[test]
fn dw_003_laundering_fails_the_type_check() {
    let problem = quartic_problem(16);
    let u_h = solve_p1(&problem);
    let out = dwr_integral_qoi(&problem, &u_h, 0.25, 0.75);
    let outcome = accept(
        &DwrQuery {
            qoi: "integral".to_string(),
            tolerance: 1e-3,
        },
        out.eta.abs(),
        None,
    );
    // The unbracketed accept is estimated; writing it into the ledger
    // claiming VERIFIED must fail the type check.
    let mut graph = fs_ledger::ColorGraph::new();
    let node = graph.source("dwr-accept", outcome.color.clone());
    let laundered = graph.derive(
        "query-report",
        &[node],
        fs_evidence::IntervalOp::Hull,
        Some(Color::Verified { lo: 0.0, hi: 1e-3 }),
        &std::collections::BTreeMap::new(),
        None,
    );
    assert!(
        laundered.is_err(),
        "an unbracketed DWR accept can never be written as verified"
    );
    verdict(
        "dw-003",
        "the adversarial upgrade of an estimated DWR accept is refused at ledger-write \
         time",
    );
}

#[test]
fn dw_004_refinement_concentrates_where_the_qoi_lives() {
    let problem = quartic_problem(20);
    let u_h = solve_p1(&problem);
    // QoI window on the RIGHT fifth of the domain.
    let out = dwr_integral_qoi(&problem, &u_h, 0.8, 1.0);
    let total: f64 = out.indicators.iter().sum();
    let right: f64 = out.indicators[12..].iter().sum();
    assert!(
        right > 0.6 * total,
        "goal-oriented indicators concentrate near the QoI window: right-40% share {:.2}",
        right / total
    );
    // Control: a CENTERED QoI does not pile mass on the right.
    let centered = dwr_integral_qoi(&problem, &u_h, 0.4, 0.6);
    let ctotal: f64 = centered.indicators.iter().sum();
    let cright: f64 = centered.indicators[12..].iter().sum();
    assert!(
        cright < 0.5 * ctotal,
        "centered QoI spreads differently: {:.2}",
        cright / ctotal
    );
    verdict(
        "dw-004",
        "right-window QoI puts >60% of indicator mass on the right; the centered \
         control does not",
    );
}

#[test]
fn dw_005_falsifier_spot_check_and_pairing() {
    // The high-fidelity spot check: a much finer solve's QoI is the
    // reference; the DWR-corrected coarse QoI must land near it.
    let coarse = quartic_problem(12);
    let u_c = solve_p1(&coarse);
    let out = dwr_integral_qoi(&coarse, &u_c, 0.25, 0.75);
    let reference = exact_qoi(&coarse, 0.25, 0.75);
    let corrected = out.j_primal + out.eta;
    let raw_err = (out.j_primal - reference).abs();
    let corr_err = (corrected - reference).abs();
    assert!(
        corr_err < 0.35 * raw_err,
        "the DWR correction moves TOWARD the high-fidelity truth: {corr_err:.3e} vs \
         raw {raw_err:.3e}"
    );
    // The falsifier pairing (Proposal 6): dwr-accept ships only once
    // paired with the high-fidelity spot check.
    let mut registry = FalsifierRegistry::standard();
    let blocked = registry.ship_gate(&["dwr-accept"]);
    assert_eq!(blocked.len(), 1, "unpaired: blocked by name");
    registry
        .register(
            "dwr-accept",
            vec![FalsifierSpec {
                name: "high-fidelity-qoi-spot-check".to_string(),
                method: "occasional full fine-grid QoI evaluation, independent of \
                         the dual machinery"
                    .to_string(),
            }],
        )
        .expect("pairing registers");
    assert!(
        registry.ship_gate(&["dwr-accept"]).is_empty(),
        "paired: ships"
    );
    verdict(
        "dw-005",
        "DWR correction converges toward the reference (falsifier honesty); the class \
         ships only once paired",
    );
}

#[test]
fn guaranteed_over_tolerance_bracket_refuses_acceptance() {
    // Bead 9sf6 F4 regression: a GUARANTEED bracket whose bound
    // exceeds tolerance is rigorous evidence the tolerance cannot be
    // confirmed — the old path silently dropped it and accepted on
    // the bare estimate. Fail closed.
    let query = DwrQuery {
        qoi: "test-qoi".to_string(),
        tolerance: 1e-3,
    };
    let bracket = Bracket {
        bound: 5e-2,
        guaranteed: true,
        source: "test-enclosure".to_string(),
    };
    // The estimate LOOKS fine (below tolerance) — that is exactly the
    // trap: an uncertified estimate must not override certified
    // inconclusiveness.
    let outcome = accept(&query, 5e-4, Some(&bracket));
    assert!(
        !outcome.accepted,
        "must refuse when rigorous evidence cannot confirm tolerance: {}",
        outcome.audit
    );
    assert!(
        outcome.audit.contains("REFUSED"),
        "audit must state the refusal: {}",
        outcome.audit
    );
    // A NON-guaranteed over-tolerance bracket still falls through to
    // the estimate path (nothing rigorous was dropped).
    let soft = Bracket {
        bound: 5e-2,
        guaranteed: false,
        source: "heuristic".to_string(),
    };
    let outcome2 = accept(&query, 5e-4, Some(&soft));
    assert!(outcome2.accepted, "heuristic brackets do not veto");
}
