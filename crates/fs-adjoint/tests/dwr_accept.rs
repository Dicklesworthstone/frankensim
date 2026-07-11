//! DWR-accept conformance (the lmp4.4 bead; runs under the `dwr-accept`
//! feature). Acceptance: DWR accepts are correct when colored estimated
//! and remain estimated until an exact QoI-dual relation is certified; the
//! laundering attempt fails the type check; DWR-driven refinement
//! concentrates where the QoI error lives; the falsifier's occasional
//! high-fidelity evaluation confirms the estimate's honesty.
#![cfg(feature = "dwr-accept")]

use fs_adjoint::dwr_accept::{
    Bracket, DwrError, DwrQuery, MAX_DWR_MESH_NODES, MAX_DWR_POLY_COEFFICIENTS, MAX_DWR_WORK_UNITS,
    accept, dwr_integral_qoi,
};
use fs_evidence::Color;
use fs_evidence::falsify::{FalsifierRegistry, FalsifierSpec};
use fs_verify::estimator::{EstimatorFamily, VerifierReport};
use fs_verify::fem1d::{MmsProblem, Poly, solve_p1};
use fs_verify::interval::Iv;

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
    let u_h = solve_p1(&problem).expect("quartic primal fixture must solve");
    let (a, b) = (0.25, 0.75);
    let out = dwr_integral_qoi(&problem, &u_h, a, b).expect("valid DWR inputs");
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
    assert!(!outcome.refused);
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
    assert!(
        !strict.refused,
        "over-tolerance is a decision, not a refusal"
    );
    verdict(
        "dw-001",
        "quartic G1 fixture: effectivity within [0.7, 1.3]; DWR-only accept is \
         estimated; too-tight tolerance rejects",
    );
}

#[test]
fn dw_002_reverified_energy_product_does_not_promote_without_a_typed_dual_relation() {
    let problem = quartic_problem(24);
    let u_h = solve_p1(&problem).expect("quartic primal fixture must solve");
    let (a, b) = (0.25, 0.75);
    let out = dwr_integral_qoi(&problem, &u_h, a, b).expect("valid DWR inputs");
    // The dual −z″ = 1_{[a,b]} has the exact piecewise
    // solution; bound its P1 error with the same equilibrated verifier
    // by posing the dual as an MMS problem on the window via a smooth
    // stand-in (a parabola matching the interior load): honest bound on
    // a HARDER dual surrogate.
    let dual_problem = MmsProblem::new(
        "dual-surrogate",
        Poly(vec![0.0, 0.5, -0.5]), // z = x(1−x)/2 solves −z″ = 1
        problem.mesh.clone(),
    );
    let dual_h = solve_p1(&dual_problem).expect("dual surrogate fixture must solve");
    let bracket = Bracket::cauchy_schwarz(&problem, &u_h, &dual_problem, &dual_h)
        .expect("both factors independently reverify");
    assert!(
        bracket.bound().is_finite(),
        "bounded diagnostic: {bracket:?}"
    );
    let query = DwrQuery {
        qoi: "integral[0.25,0.75]".to_string(),
        tolerance: bracket.bound() * 1.5,
    };
    let outcome = accept(&query, out.eta.abs(), Some(&bracket));
    assert!(outcome.accepted);
    assert!(
        matches!(outcome.color, Color::Estimated { .. }),
        "an unbound dual relation cannot mint Verified: {:?}",
        outcome.color
    );
    assert!(
        outcome.audit.contains("QoI-dual relation unverified"),
        "audit retains the exact no-claim: {}",
        outcome.audit
    );
    verdict(
        "dw-002",
        "both energy factors are independently reverified, but the product remains an \
         Estimated diagnostic until a typed QoI-dual relation exists",
    );
}

#[test]
fn forged_public_verifier_report_is_not_bracket_authority() {
    let forged = VerifierReport {
        bound: Iv { lo: 0.0, hi: 0.0 },
        accept: true,
        color: Some(Color::Verified { lo: 0.0, hi: 0.0 }),
        tolerance: f64::MAX,
        family: EstimatorFamily::EquilibratedFlux.id(),
        flux_hash: 0,
        refusal: None,
    };
    assert!(
        forged.accept,
        "raw reports are demonstrably caller-forgeable"
    );
    let outcome = accept(
        &DwrQuery {
            qoi: "forgery-probe".to_string(),
            tolerance: 1.0,
        },
        0.0,
        None,
    );
    assert!(matches!(outcome.color, Color::Estimated { .. }));
    // There is intentionally no API that consumes `forged`: Bracket fields are
    // private and its only constructor reruns verification from exact inputs.
}

#[test]
fn dw_003_laundering_fails_the_type_check() {
    let problem = quartic_problem(16);
    let u_h = solve_p1(&problem).expect("quartic primal fixture must solve");
    let out = dwr_integral_qoi(&problem, &u_h, 0.25, 0.75).expect("valid DWR inputs");
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
    let node = graph
        .source("dwr-accept", outcome.color.clone())
        .expect("unbracketed DWR accept is Estimated");
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
    let u_h = solve_p1(&problem).expect("quartic primal fixture must solve");
    // QoI window on the RIGHT fifth of the domain.
    let out = dwr_integral_qoi(&problem, &u_h, 0.8, 1.0).expect("valid DWR inputs");
    let total: f64 = out.indicators.iter().sum();
    let right: f64 = out.indicators[12..].iter().sum();
    assert!(
        right > 0.6 * total,
        "goal-oriented indicators concentrate near the QoI window: right-40% share {:.2}",
        right / total
    );
    // Control: a CENTERED QoI does not pile mass on the right.
    let centered = dwr_integral_qoi(&problem, &u_h, 0.4, 0.6).expect("valid DWR inputs");
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
    let u_c = solve_p1(&coarse).expect("coarse falsifier fixture must solve");
    let out = dwr_integral_qoi(&coarse, &u_c, 0.25, 0.75).expect("valid DWR inputs");
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
fn unrelated_over_tolerance_energy_product_cannot_veto_the_dwr_decision() {
    let query = DwrQuery {
        qoi: "test-qoi".to_string(),
        tolerance: 1e-3,
    };
    let primal = quartic_problem(4);
    let primal_h = solve_p1(&primal).expect("primal fixture must solve");
    let unrelated_dual = MmsProblem::new(
        "unrelated-dual",
        Poly(vec![0.0, 4.0, -4.0]),
        primal.mesh.clone(),
    );
    let unrelated_h = solve_p1(&unrelated_dual).expect("unrelated dual fixture must solve");
    let bracket = Bracket::cauchy_schwarz(&primal, &primal_h, &unrelated_dual, &unrelated_h)
        .expect("both unrelated energy factors still reverify");
    let outcome = accept(&query, 5e-4, Some(&bracket));
    assert!(
        outcome.accepted,
        "an unrelated dual product cannot veto an Estimated DWR decision: {}",
        outcome.audit
    );
    assert!(
        matches!(outcome.color, Color::Estimated { .. }),
        "unrelated genuine reports cannot promote"
    );
    assert!(
        outcome.audit.contains("QoI-dual relation unverified"),
        "audit must state the no-claim: {}",
        outcome.audit
    );
}

#[test]
fn malformed_accept_inputs_fail_closed_without_minting_invalid_colors() {
    use fs_evidence::validate_color_payload;

    let base = DwrQuery {
        qoi: "hostile display label (not a machine id)".to_string(),
        tolerance: 1e-3,
    };
    for estimate in [f64::NAN, f64::NEG_INFINITY, -1.0] {
        let outcome = accept(&base, estimate, None);
        assert!(!outcome.accepted);
        assert!(outcome.refused);
        validate_color_payload(&outcome.color).expect("refusal color remains structurally valid");
        assert!(matches!(
            outcome.color,
            Color::Estimated { dispersion, .. } if dispersion.is_infinite()
        ));
    }

    for tolerance in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        let outcome = accept(
            &DwrQuery {
                tolerance,
                ..base.clone()
            },
            1e-4,
            None,
        );
        assert!(!outcome.accepted);
        assert!(outcome.refused);
        validate_color_payload(&outcome.color).expect("invalid tolerance refusal is valid");
    }

    let problem = quartic_problem(4);
    let candidate = solve_p1(&problem).expect("malformed-input control fixture must solve");
    assert!(
        Bracket::cauchy_schwarz(
            &problem,
            &candidate[..candidate.len() - 1],
            &problem,
            &candidate,
        )
        .is_err(),
        "malformed candidates fail before verifier execution"
    );
    let bracket = Bracket::cauchy_schwarz(&problem, &candidate, &problem, &candidate)
        .expect("valid diagnostic");
    let independent = accept(&base, f64::NAN, Some(&bracket));
    assert!(
        !independent.accepted,
        "an unbound energy product cannot discharge malformed DWR"
    );
    assert!(independent.refused);
    validate_color_payload(&independent.color).expect("refusal color is valid");
}

#[test]
fn dwr_two_node_mesh_is_a_finite_supported_boundary_case() {
    let problem = MmsProblem::new("two-node", Poly(vec![0.0]), vec![0.0, 1.0]);
    let output = dwr_integral_qoi(&problem, &[0.0, 0.0], 0.0, 1.0)
        .expect("two boundary nodes refine to one interior dual degree of freedom");
    assert_eq!(output.indicators.len(), 1);
    assert!(output.j_primal.is_finite());
    assert!(output.eta.is_finite());
    assert!(output.indicators[0].is_finite());
}

#[test]
fn dwr_refuses_malformed_meshes_and_candidates_before_execution() {
    for mesh in [Vec::new(), vec![0.0]] {
        let problem = MmsProblem::new("too-small", Poly(vec![0.0]), mesh.clone());
        assert!(matches!(
            dwr_integral_qoi(&problem, &vec![0.0; mesh.len()], 0.0, 1.0),
            Err(DwrError::MeshNodeCount { .. })
        ));
    }

    let base = MmsProblem::new("base", Poly(vec![0.0]), vec![0.0, 1.0]);
    assert!(matches!(
        dwr_integral_qoi(&base, &[0.0], 0.0, 1.0),
        Err(DwrError::CandidateLengthMismatch { .. })
    ));
    assert!(matches!(
        dwr_integral_qoi(&base, &[0.0, f64::NAN], 0.0, 1.0),
        Err(DwrError::NonFiniteCandidate { index: 1 })
    ));

    for mesh in [vec![0.0, 0.0], vec![1.0, 0.0]] {
        let problem = MmsProblem::new("unordered", Poly(vec![0.0]), mesh);
        assert!(matches!(
            dwr_integral_qoi(&problem, &[0.0, 0.0], 0.0, 1.0),
            Err(DwrError::NonIncreasingMeshCell { cell: 0 })
        ));
    }
    let nonfinite = MmsProblem::new("nonfinite", Poly(vec![0.0]), vec![0.0, f64::INFINITY]);
    assert!(matches!(
        dwr_integral_qoi(&nonfinite, &[0.0, 0.0], 0.0, 1.0),
        Err(DwrError::NonFiniteMeshNode { index: 1 })
    ));
    let overflowing_width = MmsProblem::new(
        "overflowing-width",
        Poly(vec![0.0]),
        vec![-f64::MAX, f64::MAX],
    );
    assert!(matches!(
        dwr_integral_qoi(&overflowing_width, &[0.0, 0.0], -1.0, 1.0),
        Err(DwrError::NonFiniteCellWidth { cell: 0 })
    ));
    let adjacent = MmsProblem::new("adjacent", Poly(vec![0.0]), vec![1.0, 1.0_f64.next_up()]);
    assert!(matches!(
        dwr_integral_qoi(&adjacent, &[0.0, 0.0], 1.0, 2.0),
        Err(DwrError::NonInteriorMidpoint { cell: 0 })
    ));
    let tiny = MmsProblem::new("tiny", Poly(vec![0.0]), vec![0.0, f64::from_bits(2)]);
    assert!(matches!(
        dwr_integral_qoi(&tiny, &[0.0, 0.0], 0.0, 1.0),
        Err(DwrError::NonFiniteReciprocal {
            cell: 0,
            refined_half: None
        })
    ));
    let refined_tiny = MmsProblem::new("refined-tiny", Poly(vec![0.0]), vec![0.0, 1.0e-308]);
    assert!((1.0_f64 / 1.0e-308_f64).is_finite());
    assert!(matches!(
        dwr_integral_qoi(&refined_tiny, &[0.0, 0.0], 0.0, 1.0),
        Err(DwrError::NonFiniteReciprocal {
            cell: 0,
            refined_half: Some(0 | 1)
        })
    ));
}

#[test]
fn dwr_refuses_invalid_windows_polynomials_and_resource_counts() {
    let base = MmsProblem::new("base", Poly(vec![0.0]), vec![0.0, 1.0]);
    for (lo, hi) in [
        (f64::NAN, 1.0),
        (0.0, f64::INFINITY),
        (1.0, 0.0),
        (0.5, 0.5),
    ] {
        assert!(matches!(
            dwr_integral_qoi(&base, &[0.0, 0.0], lo, hi),
            Err(DwrError::InvalidQoiWindow { .. })
        ));
    }

    let empty_poly = MmsProblem::new("empty-poly", Poly(Vec::new()), vec![0.0, 1.0]);
    assert!(matches!(
        dwr_integral_qoi(&empty_poly, &[0.0, 0.0], 0.0, 1.0),
        Err(DwrError::PolynomialCoefficientCount { count: 0, .. })
    ));
    let nonfinite_poly = MmsProblem::new("nan-poly", Poly(vec![0.0, f64::NAN]), vec![0.0, 1.0]);
    assert!(matches!(
        dwr_integral_qoi(&nonfinite_poly, &[0.0, 0.0], 0.0, 1.0),
        Err(DwrError::NonFinitePolynomialCoefficient { index: 1 })
    ));
    let oversized_poly = MmsProblem::new(
        "oversized-poly",
        Poly(vec![0.0; MAX_DWR_POLY_COEFFICIENTS + 1]),
        vec![0.0, 1.0],
    );
    assert!(matches!(
        dwr_integral_qoi(&oversized_poly, &[0.0, 0.0], 0.0, 1.0),
        Err(DwrError::PolynomialCoefficientCount { .. })
    ));

    let too_many_candidates = vec![0.0; MAX_DWR_MESH_NODES + 1];
    assert!(matches!(
        dwr_integral_qoi(&base, &too_many_candidates, 0.0, 1.0),
        Err(DwrError::CandidateNodeCount { .. })
    ));
    let oversized_mesh = MmsProblem::new(
        "oversized-mesh",
        Poly(vec![0.0]),
        vec![0.0; MAX_DWR_MESH_NODES + 1],
    );
    assert!(matches!(
        dwr_integral_qoi(&oversized_mesh, &[], 0.0, 1.0),
        Err(DwrError::MeshNodeCount { .. })
    ));

    let mesh_nodes = 100_001;
    let coefficients = 100;
    assert!(mesh_nodes <= MAX_DWR_MESH_NODES);
    assert!(coefficients <= MAX_DWR_POLY_COEFFICIENTS);
    let cross_product = MmsProblem::new(
        "work-cross-product",
        Poly(vec![0.0; coefficients]),
        vec![0.0; mesh_nodes],
    );
    let cross_product_candidate = vec![0.0; mesh_nodes];
    assert!(matches!(
        dwr_integral_qoi(&cross_product, &cross_product_candidate, 0.0, 1.0),
        Err(DwrError::WorkBudgetExceeded {
            estimated_work: Some(work),
            maximum: MAX_DWR_WORK_UNITS,
            ..
        }) if work > MAX_DWR_WORK_UNITS
    ));
}

#[test]
fn dwr_refuses_finite_inputs_that_overflow_derived_arithmetic() {
    let overflowing_forcing = MmsProblem::new(
        "overflowing-forcing",
        Poly(vec![0.0, 0.0, f64::MAX]),
        vec![0.0, 1.0],
    );
    assert!(matches!(
        dwr_integral_qoi(&overflowing_forcing, &[0.0, 0.0], 0.0, 1.0),
        Err(DwrError::NonFiniteDerived {
            quantity: "forcing coefficient",
            ..
        })
    ));

    let base = MmsProblem::new("base", Poly(vec![0.0]), vec![0.0, 1.0]);
    assert!(matches!(
        dwr_integral_qoi(&base, &[-f64::MAX, f64::MAX], 2.0, 3.0),
        Err(DwrError::NonFiniteDerived {
            quantity: "primal slope",
            index: Some(0)
        })
    ));
}
