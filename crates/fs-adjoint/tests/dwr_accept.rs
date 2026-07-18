//! DWR-accept conformance (the lmp4.4 bead; runs under the `dwr-accept`
//! feature). Acceptance: DWR accepts are correct when colored estimated
//! and remain estimated until an exact QoI-dual relation is certified; the
//! laundering attempt fails the type check; DWR-driven refinement
//! concentrates where the QoI error lives; the falsifier's occasional
//! high-fidelity evaluation confirms the estimate's honesty.
#![cfg(feature = "dwr-accept")]

use fs_adjoint::dwr_accept::{
    AcceptOutcome, Bracket, BracketError, DWR_EVIDENCE_IDENTITY_VERSION, DWR_POLL_POLICY_VERSION,
    DWR_POLL_STRIDE_ITEMS, DWR_QUERY_SCHEMA_VERSION, DWR_WORK_PLAN_VERSION, DwrError, DwrOutput,
    DwrQuery, EndpointConvention, MAX_BRACKET_MESH_NODES, MAX_DWR_MESH_NODES,
    MAX_DWR_POLY_COEFFICIENTS, MAX_DWR_QOI_BYTES, MAX_DWR_WORK_UNITS, QoiSemantics,
    accept as accept_with_cx, dwr_integral_qoi as dwr_integral_qoi_with_cx,
};
use fs_evidence::Color;
use fs_evidence::falsify::{FalsifierRegistry, FalsifierSpec};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_verify::estimator::{EstimatorFamily, VerifierReport};
use fs_verify::fem1d::{
    Fem1dError, MAX_FEM1D_CLASS_NAME_BYTES, MAX_FEM1D_MESH_NODES, MAX_FEM1D_POLY_COEFFICIENTS,
    MmsClass, MmsProblem, Poly, solve_p1,
};
use fs_verify::interval::Iv;

fn with_cx<R>(
    cancelled: bool,
    mode: ExecMode,
    budget: Budget,
    stream: StreamKey,
    f: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let result = pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, stream, budget, mode);
        f(&cx)
    });
    let stats = pool.stats();
    assert!(
        stats.quiescent(),
        "Cx arena must be quiescent after scope: {}",
        stats.to_json()
    );
    result
}

fn default_stream() -> StreamKey {
    StreamKey {
        seed: 0xD0_00_00_01,
        kernel_id: 0xAD_10_17,
        tile: 3,
        iteration: 5,
    }
}

fn with_default_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_cx(
        false,
        ExecMode::Deterministic,
        Budget::INFINITE,
        default_stream(),
        f,
    )
}

fn dwr_integral_qoi(
    problem: &MmsProblem,
    candidate: &[f64],
    w_lo: f64,
    w_hi: f64,
) -> Result<DwrOutput, DwrError> {
    with_default_cx(|cx| {
        dwr_integral_qoi_with_cx(problem, candidate, w_lo, w_hi, cx, &VirtualClock::new())
    })
}

fn accept(query: &DwrQuery, estimate: &DwrOutput, bracket: Option<&Bracket>) -> AcceptOutcome {
    with_default_cx(|cx| accept_with_cx(query, estimate, bracket, cx, &VirtualClock::new()))
        .expect("healthy DWR acceptance context")
}

/// A small genuinely solved estimate for acceptance tests whose concern is
/// not the estimate itself (sealed-accept has no scalar path — sj31i.1).
fn solved_estimate() -> DwrOutput {
    let problem = quartic_problem(4);
    let u_h = solve_p1(&problem).expect("standard acceptance fixture must solve");
    dwr_integral_qoi(&problem, &u_h, 0.25, 0.75).expect("standard acceptance estimate")
}

fn cauchy_schwarz(
    primal_problem: &MmsProblem,
    primal_candidate: &[f64],
    dual_problem: &MmsProblem,
    dual_candidate: &[f64],
) -> Result<Bracket, BracketError> {
    with_default_cx(|cx| {
        Bracket::cauchy_schwarz(
            primal_problem,
            primal_candidate,
            dual_problem,
            dual_candidate,
            cx,
            &VirtualClock::new(),
        )
    })
}

fn assert_same_dwr_semantics(left: &DwrOutput, right: &DwrOutput) {
    assert_eq!(left.j_primal().to_bits(), right.j_primal().to_bits());
    assert_eq!(left.eta().to_bits(), right.eta().to_bits());
    assert_eq!(left.indicators().len(), right.indicators().len());
    for (left, right) in left.indicators().iter().zip(right.indicators()) {
        assert_eq!(left.to_bits(), right.to_bits());
    }
    let (left_lo, left_hi) = left.window();
    let (right_lo, right_hi) = right.window();
    assert_eq!(left_lo.to_bits(), right_lo.to_bits());
    assert_eq!(left_hi.to_bits(), right_hi.to_bits());
}

fn assert_same_accept_semantics(left: &AcceptOutcome, right: &AcceptOutcome) {
    assert_eq!(left.accepted(), right.accepted());
    assert_eq!(left.color(), right.color());
    assert_eq!(left.refused(), right.refused());
    assert_eq!(left.audit(), right.audit());
}

fn assert_same_bracket_semantics(left: &Bracket, right: &Bracket) {
    assert_eq!(left.bound().to_bits(), right.bound().to_bits());
    assert_eq!(left.source(), right.source());
}

fn poly(coefficients: Vec<f64>) -> Poly {
    Poly::new(coefficients).expect("valid DWR fixture polynomial")
}

fn admitted_problem(name: &str, coefficients: Vec<f64>, mesh: Vec<f64>) -> MmsProblem {
    MmsProblem::new(name, poly(coefficients), mesh).expect("valid DWR fixture problem")
}

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
    #[allow(clippy::cast_precision_loss)]
    let mesh: Vec<f64> = (0..=cells).map(|k| k as f64 / cells as f64).collect();
    admitted_problem("dwr-quartic", vec![0.0, 1.0, -1.0, 2.0, -2.0], mesh)
}

/// The exact QoI ∫_{a}^{b} u dx from the polynomial antiderivative.
fn exact_qoi(problem: &MmsProblem, a: f64, b: f64) -> f64 {
    let big_u = problem
        .exact_solution()
        .antiderive()
        .expect("quartic fixture antiderivative stays inside the shared cap");
    big_u.eval(b) - big_u.eval(a)
}

#[test]
fn dw_001_g1_effectivity_and_estimated_accept() {
    let problem = quartic_problem(16);
    let u_h = solve_p1(&problem).expect("quartic primal fixture must solve");
    let (a, b) = (0.25, 0.75);
    let out = dwr_integral_qoi(&problem, &u_h, a, b).expect("valid DWR inputs");
    let true_err = exact_qoi(&problem, a, b) - out.j_primal();
    let effectivity = out.eta() / true_err;
    assert!(
        (0.7..=1.3).contains(&effectivity),
        "G1: DWR effectivity near 1: eta {:.3e} vs true {:.3e} ({effectivity:.2})",
        out.eta(),
        true_err
    );
    // The DWR-only accept at a discharging tolerance is ESTIMATED.
    let query = DwrQuery {
        qoi: "integral[0.25,0.75]".to_string(),
        schema_version: DWR_QUERY_SCHEMA_VERSION,
        window: (0.25, 0.75),
        endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
        semantics: QoiSemantics::WindowedIntegral,
        dims: [0; 6],
        tolerance: 1e-3,
    };
    let outcome = accept(&query, &out, None);
    assert!(outcome.accepted());
    assert!(!outcome.refused());
    assert!(
        matches!(outcome.color(), Color::Estimated { .. }),
        "DWR constants are not guaranteed: {:?}",
        outcome.color()
    );
    // And a too-tight tolerance rejects.
    let strict = accept(
        &DwrQuery {
            qoi: "integral".to_string(),
            schema_version: DWR_QUERY_SCHEMA_VERSION,
            window: (0.25, 0.75),
            endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
            semantics: QoiSemantics::WindowedIntegral,
            dims: [0; 6],
            tolerance: 1e-12,
        },
        &out,
        None,
    );
    assert!(!strict.accepted(), "no silent discharge");
    assert!(
        !strict.refused(),
        "over-tolerance is a decision, not a refusal"
    );
    verdict(
        "dw-001",
        "quartic G1 fixture: effectivity within [0.7, 1.3]; DWR-only accept is \
         estimated; too-tight tolerance rejects",
    );
}

#[test]
fn dw_001a_narrow_nonaligned_qoi_window_cannot_alias_to_zero() {
    let problem = admitted_problem("narrow-window", vec![0.0, 0.5, -0.5], vec![0.0, 0.5, 1.0]);
    let candidate = [0.0, 1.0, 0.0];
    let output = dwr_integral_qoi(&problem, &candidate, 0.01, 0.011)
        .expect("a representable nonaligned window is admitted");
    let expected = 0.011_f64.powi(2) - 0.01_f64.powi(2);
    assert!(
        (output.j_primal() - expected).abs() <= 4.0 * f64::EPSILON * expected,
        "exact clipped P1 integral {} differs from {expected}",
        output.j_primal()
    );
    assert!(
        output.eta() != 0.0,
        "the clipped dual load must not disappear with the whole-cell Gauss nodes"
    );

    let refined_problem = admitted_problem(
        "narrow-window-refined",
        vec![0.0, 0.5, -0.5],
        vec![0.0, 0.25, 0.5, 0.75, 1.0],
    );
    let refined = dwr_integral_qoi(&refined_problem, &[0.0, 0.5, 1.0, 0.5, 0.0], 0.01, 0.011)
        .expect("mesh refinement preserves the same P1 field/window integral");
    assert!(
        (refined.j_primal() - output.j_primal()).abs() <= 32.0 * f64::EPSILON * expected,
        "G3: inserting nodes outside the clipped window preserves the QoI"
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
    let dual_problem = admitted_problem(
        "dual-surrogate",
        vec![0.0, 0.5, -0.5], // z = x(1−x)/2 solves −z″ = 1
        problem.mesh().to_vec(),
    );
    let dual_h = solve_p1(&dual_problem).expect("dual surrogate fixture must solve");
    let bracket = cauchy_schwarz(&problem, &u_h, &dual_problem, &dual_h)
        .expect("both factors independently reverify");
    assert!(
        bracket.bound().is_finite(),
        "bounded diagnostic: {bracket:?}"
    );
    let query = DwrQuery {
        qoi: "integral[0.25,0.75]".to_string(),
        schema_version: DWR_QUERY_SCHEMA_VERSION,
        window: (0.25, 0.75),
        endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
        semantics: QoiSemantics::WindowedIntegral,
        dims: [0; 6],
        tolerance: bracket.bound() * 1.5,
    };
    let outcome = accept(&query, &out, Some(&bracket));
    assert!(outcome.accepted());
    assert!(
        matches!(outcome.color(), Color::Estimated { .. }),
        "an unbound dual relation cannot mint Verified: {:?}",
        outcome.color()
    );
    assert!(
        outcome.audit().contains("QoI-dual relation unverified"),
        "audit retains the exact no-claim: {}",
        outcome.audit()
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
            schema_version: DWR_QUERY_SCHEMA_VERSION,
            window: (0.25, 0.75),
            endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
            semantics: QoiSemantics::WindowedIntegral,
            dims: [0; 6],
            tolerance: 1.0,
        },
        &solved_estimate(),
        None,
    );
    assert!(matches!(outcome.color(), Color::Estimated { .. }));
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
            schema_version: DWR_QUERY_SCHEMA_VERSION,
            window: (0.25, 0.75),
            endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
            semantics: QoiSemantics::WindowedIntegral,
            dims: [0; 6],
            tolerance: 1e-3,
        },
        &out,
        None,
    );
    // The unbracketed accept is estimated; writing it into the ledger
    // claiming VERIFIED must fail the type check.
    let mut graph = fs_ledger::ColorGraph::new();
    let node = graph
        .source("dwr-accept", outcome.color().clone())
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
    let total: f64 = out.indicators().iter().sum();
    let right: f64 = out.indicators()[12..].iter().sum();
    assert!(
        right > 0.6 * total,
        "goal-oriented indicators concentrate near the QoI window: right-40% share {:.2}",
        right / total
    );
    // Control: a CENTERED QoI does not pile mass on the right.
    let centered = dwr_integral_qoi(&problem, &u_h, 0.4, 0.6).expect("valid DWR inputs");
    let ctotal: f64 = centered.indicators().iter().sum();
    let cright: f64 = centered.indicators()[12..].iter().sum();
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
    // Manufactured-solution oracle spot check: the fixture's analytic QoI is
    // independent of the DWR dual calculation, and the corrected coarse QoI
    // must land nearer it.
    let coarse = quartic_problem(12);
    let u_c = solve_p1(&coarse).expect("coarse falsifier fixture must solve");
    let out = dwr_integral_qoi(&coarse, &u_c, 0.25, 0.75).expect("valid DWR inputs");
    let reference = exact_qoi(&coarse, 0.25, 0.75);
    let corrected = out.j_primal() + out.eta();
    let raw_err = (out.j_primal() - reference).abs();
    let corr_err = (corrected - reference).abs();
    assert!(
        corr_err < 0.35 * raw_err,
        "the DWR correction moves TOWARD the manufactured-solution oracle: {corr_err:.3e} vs \
         raw {raw_err:.3e}"
    );
    // The declaration catalog (Proposal 6) reports dwr-accept until its
    // manufactured-solution oracle is declared. This is not release authority.
    let mut registry = FalsifierRegistry::standard();
    let blocked = registry
        .catalog_gate(&["dwr-accept"])
        .expect("bounded valid catalog query");
    assert_eq!(blocked.len(), 1, "unpaired: reported by name");
    registry
        .register(
            "dwr-accept",
            vec![FalsifierSpec {
                name: "manufactured-solution-qoi-oracle".to_string(),
                method: "analytic fixture QoI evaluation independently derived from \
                         the manufactured solution, outside the dual machinery"
                    .to_string(),
            }],
        )
        .expect("pairing registers");
    assert!(
        registry
            .catalog_gate(&["dwr-accept"])
            .expect("bounded valid catalog query")
            .is_empty(),
        "paired: catalog metadata is complete"
    );
    verdict(
        "dw-005",
        "DWR correction converges toward the manufactured-solution oracle; catalog \
         completeness is distinct from exact-instance release admission",
    );
}

#[test]
fn unrelated_over_tolerance_energy_product_cannot_veto_the_dwr_decision() {
    let estimate = solved_estimate();
    // Discharging tolerance derived from the sealed estimate itself: the
    // concern here is bracket-veto behavior, not the estimate magnitude.
    let query = DwrQuery {
        qoi: "test-qoi".to_string(),
        schema_version: DWR_QUERY_SCHEMA_VERSION,
        window: (0.25, 0.75),
        endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
        semantics: QoiSemantics::WindowedIntegral,
        dims: [0; 6],
        tolerance: estimate.eta().abs().max(1e-12) * 2.0,
    };
    let primal = quartic_problem(4);
    let primal_h = solve_p1(&primal).expect("primal fixture must solve");
    let unrelated_dual = admitted_problem(
        "unrelated-dual",
        vec![0.0, 4.0, -4.0],
        primal.mesh().to_vec(),
    );
    let unrelated_h = solve_p1(&unrelated_dual).expect("unrelated dual fixture must solve");
    let bracket = cauchy_schwarz(&primal, &primal_h, &unrelated_dual, &unrelated_h)
        .expect("both unrelated energy factors still reverify");
    let outcome = accept(&query, &estimate, Some(&bracket));
    assert!(
        outcome.accepted(),
        "an unrelated dual product cannot veto an Estimated DWR decision: {}",
        outcome.audit()
    );
    assert!(
        matches!(outcome.color(), Color::Estimated { .. }),
        "unrelated genuine reports cannot promote"
    );
    assert!(
        outcome.audit().contains("QoI-dual relation unverified"),
        "audit must state the no-claim: {}",
        outcome.audit()
    );
}

#[test]
fn malformed_accept_inputs_fail_closed_without_minting_invalid_colors() {
    use fs_evidence::validate_color_payload;

    let base = DwrQuery {
        qoi: "hostile display label (not a machine id)".to_string(),
        schema_version: DWR_QUERY_SCHEMA_VERSION,
        window: (0.25, 0.75),
        endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
        semantics: QoiSemantics::WindowedIntegral,
        dims: [0; 6],
        tolerance: 1e-3,
    };
    // Scalar injection (NaN/-inf/-1.0 estimates) is now UNREPRESENTABLE:
    // acceptance consumes only a constructor-sealed DwrOutput whose eta is
    // finite by construction (bead sj31i.1) — the malformed legs below
    // exercise the surviving refusal surface (tolerance) instead.

    for tolerance in [f64::NAN, f64::INFINITY, 0.0, -1.0] {
        let outcome = accept(
            &DwrQuery {
                tolerance,
                ..base.clone()
            },
            &solved_estimate(),
            None,
        );
        assert!(!outcome.accepted());
        assert!(outcome.refused());
        validate_color_payload(outcome.color()).expect("invalid tolerance refusal is valid");
    }

    let problem = quartic_problem(4);
    let candidate = solve_p1(&problem).expect("malformed-input control fixture must solve");
    assert!(
        cauchy_schwarz(
            &problem,
            &candidate[..candidate.len() - 1],
            &problem,
            &candidate,
        )
        .is_err(),
        "malformed candidates fail before verifier execution"
    );
    let bracket =
        cauchy_schwarz(&problem, &candidate, &problem, &candidate).expect("valid diagnostic");
    let independent = accept(
        &DwrQuery {
            tolerance: f64::NAN,
            ..base.clone()
        },
        &solved_estimate(),
        Some(&bracket),
    );
    assert!(
        !independent.accepted(),
        "an unbound energy product cannot discharge malformed DWR"
    );
    assert!(independent.refused());
    validate_color_payload(independent.color()).expect("refusal color is valid");
}

#[test]
fn typed_functional_binding_rejects_invalid_and_cross_semantic_replay() {
    let estimate = solved_estimate();
    let base = DwrQuery {
        qoi: "typed-functional-binding".to_string(),
        schema_version: DWR_QUERY_SCHEMA_VERSION,
        tolerance: 1.0,
        window: (0.25, 0.75),
        endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
        semantics: QoiSemantics::WindowedIntegral,
        dims: [0; 6],
    };

    let cases = [
        (
            "non-finite window",
            DwrQuery {
                window: (f64::NAN, 0.75),
                ..base.clone()
            },
            "functional window must be finite and increasing",
        ),
        (
            "empty window",
            DwrQuery {
                window: (0.5, 0.5),
                ..base.clone()
            },
            "functional window must be finite and increasing",
        ),
        (
            "reversed window",
            DwrQuery {
                window: (0.8, 0.2),
                ..base.clone()
            },
            "functional window must be finite and increasing",
        ),
        (
            "cross-window replay",
            DwrQuery {
                window: (0.2, 0.8),
                ..base.clone()
            },
            "does not match the sealed estimate's retained window",
        ),
        (
            "dimension relabel",
            DwrQuery {
                dims: [1, 0, 0, 0, 0, 0],
                ..base.clone()
            },
            "does not match the sealed estimate dimensions",
        ),
        (
            "unknown schema",
            DwrQuery {
                schema_version: DWR_QUERY_SCHEMA_VERSION + 1,
                ..base.clone()
            },
            "functional schema version",
        ),
    ];

    for (case, query, diagnostic) in cases {
        let outcome = accept(&query, &estimate, None);
        assert!(!outcome.accepted(), "{case} must not discharge the query");
        assert!(outcome.refused(), "{case} must fail closed as malformed");
        assert!(
            outcome.audit().contains(diagnostic),
            "{case} emitted the wrong bounded diagnostic: {}",
            outcome.audit()
        );
    }

    verdict(
        "typed-functional-binding",
        "invalid windows plus cross-window, dimensional, and schema relabeling refuse distinctly",
    );
}

#[test]
fn dwr_two_node_mesh_is_a_finite_supported_boundary_case() {
    let problem = admitted_problem("two-node", vec![0.0], vec![0.0, 1.0]);
    let output = dwr_integral_qoi(&problem, &[0.0, 0.0], 0.0, 1.0)
        .expect("two boundary nodes refine to one interior dual degree of freedom");
    assert_eq!(output.indicators().len(), 1);
    assert!(output.j_primal().is_finite());
    assert!(output.eta().is_finite());
    assert!(output.indicators()[0].is_finite());
}

#[test]
fn dwr_relies_on_problem_admission_and_refuses_bad_candidates() {
    for mesh in [Vec::new(), vec![0.0]] {
        assert!(matches!(
            MmsProblem::new("too-small", poly(vec![0.0]), mesh),
            Err(Fem1dError::ResourceLimit {
                resource: "mesh nodes",
                ..
            })
        ));
    }

    let base = admitted_problem("base", vec![0.0], vec![0.0, 1.0]);
    assert!(matches!(
        dwr_integral_qoi(&base, &[0.0], 0.0, 1.0),
        Err(DwrError::CandidateLengthMismatch { .. })
    ));
    assert!(matches!(
        dwr_integral_qoi(&base, &[0.0, f64::NAN], 0.0, 1.0),
        Err(DwrError::NonFiniteCandidate { index: 1 })
    ));

    for mesh in [vec![0.0, 0.5, 0.5, 1.0], vec![0.0, 0.75, 0.5, 1.0]] {
        assert!(matches!(
            MmsProblem::new("unordered", poly(vec![0.0]), mesh),
            Err(Fem1dError::NonIncreasingMeshCell { .. })
        ));
    }
    assert!(matches!(
        MmsProblem::new("nonfinite", poly(vec![0.0]), vec![0.0, f64::NAN, 1.0],),
        Err(Fem1dError::NonFiniteMeshNode { index: 1 })
    ));
    assert!(matches!(
        MmsProblem::new("tiny", poly(vec![0.0]), vec![0.0, f64::from_bits(2), 1.0],),
        Err(Fem1dError::NonFiniteReciprocal { cell: 0 })
    ));

    let adjacent = admitted_problem("adjacent", vec![0.0], vec![0.0, 1.0_f64.next_down(), 1.0]);
    assert!(matches!(
        dwr_integral_qoi(&adjacent, &[0.0, 0.0, 0.0], 0.0, 1.0),
        Err(DwrError::NonInteriorMidpoint { cell: 1 })
    ));

    let refined_tiny = admitted_problem("refined-tiny", vec![0.0], vec![0.0, 1.0e-308, 1.0]);
    assert!((1.0_f64 / 1.0e-308_f64).is_finite());
    assert!(matches!(
        dwr_integral_qoi(&refined_tiny, &[0.0, 0.0, 0.0], 0.0, 1.0),
        Err(DwrError::NonFiniteReciprocal {
            cell: 0,
            refined_half: Some(0 | 1)
        })
    ));
}

#[test]
fn dwr_refuses_invalid_windows_and_resource_counts_at_owner_boundaries() {
    let base = admitted_problem("base", vec![0.0], vec![0.0, 1.0]);
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

    assert!(matches!(
        Poly::new(Vec::new()),
        Err(Fem1dError::PolynomialCoefficientCount { count: 0, .. })
    ));
    assert!(matches!(
        Poly::new(vec![0.0, f64::NAN]),
        Err(Fem1dError::NonFinitePolynomialCoefficient { index: 1, .. })
    ));
    let mut too_many_semantic_coefficients = vec![0.0; MAX_DWR_POLY_COEFFICIENTS + 1];
    *too_many_semantic_coefficients.last_mut().unwrap() = 1.0;
    assert!(matches!(
        Poly::new(too_many_semantic_coefficients),
        Err(Fem1dError::PolynomialCoefficientCount { .. })
    ));
    assert_eq!(
        MAX_DWR_POLY_COEFFICIENTS, MAX_FEM1D_POLY_COEFFICIENTS,
        "DWR must not advertise a larger polynomial class than fs-verify admits"
    );
    assert_eq!(MAX_DWR_MESH_NODES, MAX_FEM1D_MESH_NODES);
    assert_eq!(MAX_BRACKET_MESH_NODES, MAX_FEM1D_MESH_NODES);

    assert!(
        matches!(
            dwr_integral_qoi(&base, &[1.0, 1.0], 0.0, 1.0),
            Err(DwrError::CandidateBoundary)
        ),
        "a constant non-homogeneous candidate must not produce a zero-error accept"
    );
    assert!(
        matches!(
            dwr_integral_qoi(&base, &[-0.0, 0.0], 0.0, 1.0),
            Err(DwrError::CandidateBoundary)
        ),
        "DWR shares fs-verify's bit-canonical +0.0 endpoint rule"
    );
    assert!(matches!(
        cauchy_schwarz(&base, &[1.0, 1.0], &base, &[0.0, 0.0]),
        Err(BracketError::InvalidInput {
            factor: "primal",
            reason: "candidate endpoints must be canonical homogeneous +0.0",
        })
    ));

    let too_many_candidates = vec![0.0; MAX_DWR_MESH_NODES + 1];
    assert!(matches!(
        dwr_integral_qoi(&base, &too_many_candidates, 0.0, 1.0),
        Err(DwrError::CandidateNodeCount { .. })
    ));
    assert!(matches!(
        MmsProblem::new(
            "oversized-mesh",
            poly(vec![0.0]),
            vec![0.0; MAX_DWR_MESH_NODES + 1],
        ),
        Err(Fem1dError::ResourceLimit {
            resource: "mesh nodes",
            ..
        })
    ));

    let maximum_admitted_work = (MAX_DWR_MESH_NODES - 1)
        .checked_mul(MAX_DWR_POLY_COEFFICIENTS * 10 + 15)
        .expect("shared finite caps have a representable work product");
    assert!(maximum_admitted_work <= MAX_DWR_WORK_UNITS);

    let maximum_qoi = DwrQuery {
        qoi: "q".repeat(MAX_DWR_QOI_BYTES),
        schema_version: DWR_QUERY_SCHEMA_VERSION,
        window: (0.25, 0.75),
        endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
        semantics: QoiSemantics::WindowedIntegral,
        dims: [0; 6],
        tolerance: 1.0,
    };
    let maximum_qoi_outcome = accept(&maximum_qoi, &solved_estimate(), None);
    assert!(maximum_qoi_outcome.accepted());
    let plus_one_qoi = DwrQuery {
        qoi: "q".repeat(MAX_DWR_QOI_BYTES + 1),
        schema_version: DWR_QUERY_SCHEMA_VERSION,
        window: (0.25, 0.75),
        endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
        semantics: QoiSemantics::WindowedIntegral,
        dims: [0; 6],
        tolerance: 1.0,
    };
    let plus_one = with_default_cx(|cx| {
        accept_with_cx(
            &plus_one_qoi,
            &solved_estimate(),
            None,
            cx,
            &VirtualClock::new(),
        )
    });
    assert!(matches!(
        plus_one,
        Err(DwrError::QoiLabelTooLong {
            bytes,
            maximum: MAX_DWR_QOI_BYTES,
        }) if bytes == MAX_DWR_QOI_BYTES + 1
    ));
}

#[test]
fn dwr_refuses_finite_inputs_with_unrepresentable_or_unresolved_derived_arithmetic() {
    let overflowing_forcing =
        MmsClass::new("overflowing-forcing", poly(vec![0.0, f64::MAX, -f64::MAX]));
    assert!(matches!(
        overflowing_forcing,
        Err(Fem1dError::NonFiniteIntermediate {
            stage: "polynomial derivative",
            index: Some(1),
        })
    ));

    let base = admitted_problem("base", vec![0.0], vec![0.0, 0.5, 1.0]);
    assert!(matches!(
        dwr_integral_qoi(&base, &[0.0, f64::MAX, 0.0], 2.0, 3.0),
        Err(DwrError::NonFiniteDerived {
            quantity: "primal slope",
            index: Some(0)
        })
    ));

    let subnormal = admitted_problem("subnormal-integral", vec![0.0], vec![0.0, 0.5, 1.0]);
    assert!(matches!(
        dwr_integral_qoi(&subnormal, &[0.0, f64::from_bits(1), 0.0], 0.0, 0.25,),
        Err(DwrError::UnresolvedZeroIntegral {
            quantity: "primal QoI",
            cell: 0,
        })
    ));

    // fs-verify 665d499 canonicalized fem1d meshes to [+0.0, 1.0]; the
    // fixture keeps its structure — an antisymmetric middle cell with the
    // window symmetric about the crossing — on dyadic in-domain nodes.
    let cancellation_problem =
        admitted_problem("zero-cancellation", vec![0.0], vec![0.0, 0.25, 0.75, 1.0]);
    let exact_cancellation =
        dwr_integral_qoi(&cancellation_problem, &[0.0, 1.0, -1.0, 0.0], 0.375, 0.625)
            .expect("an exactly symmetric clipped P1 integral is admitted");
    assert_eq!(exact_cancellation.j_primal().to_bits(), 0.0_f64.to_bits());
    // The unresolved-zero fixture needs BIT-exact cancellation whose zero
    // the prover cannot certify: the window offsets from the cell edges
    // are decimal-equal but float-unequal (0.2505 - 0.25 differs from
    // 0.75 - 0.7495 at the ULP), so `two_diff` symmetry fails while the
    // subnormal-scale halves still cancel exactly — the same structure the
    // pre-canonicalization fixture used on [-1, 2].
    assert!(matches!(
        dwr_integral_qoi(
            &cancellation_problem,
            &[0.0, 1.0e-300, -1.0e-300, 0.0],
            0.2505,
            0.7495,
        ),
        Err(DwrError::UnresolvedZeroIntegral {
            quantity: "primal QoI",
            cell: 1,
        })
    ));
}

#[test]
fn dwr_consumes_canonical_mms_class_and_problem_identities() {
    let normalized = MmsClass::new("dwr-identity", poly(vec![-0.0, 1.0, -1.0, -0.0]))
        .expect("canonical admitted class");
    let ordinary = MmsClass::new("dwr-identity", poly(vec![0.0, 1.0, -1.0]))
        .expect("same canonical admitted class");
    assert_eq!(normalized.identity(), ordinary.identity());
    assert_eq!(normalized.canonical_bytes(), ordinary.canonical_bytes());
    assert_eq!(normalized.forcing().coefficients(), &[2.0]);

    let coarse = MmsProblem::from_class(normalized.clone(), vec![-0.0, 0.5, 1.0])
        .expect("canonicalized coarse problem");
    let same =
        MmsProblem::from_class(ordinary, vec![0.0, 0.5, 1.0]).expect("same canonical problem");
    let refined = coarse
        .with_mesh(vec![0.0, 0.25, 0.5, 1.0])
        .expect("same class on a refined mesh");
    assert_eq!(coarse.identity(), same.identity());
    assert_eq!(coarse.class().identity(), normalized.identity());
    assert_eq!(refined.class().identity(), normalized.identity());
    assert_ne!(coarse.identity(), refined.identity());

    let renamed = MmsClass::new("dwr-identity-renamed", poly(vec![0.0, 1.0, -1.0]))
        .expect("renamed admitted class");
    let rescaled =
        MmsClass::new("dwr-identity", poly(vec![0.0, 2.0, -2.0])).expect("rescaled admitted class");
    assert_ne!(normalized.identity(), renamed.identity());
    assert_ne!(normalized.identity(), rescaled.identity());
}

#[test]
#[allow(clippy::cast_precision_loss)]
fn g4_hostile_maximum_dwr_shape_preflights_before_initial_cancellation() {
    let denominator = (MAX_DWR_MESH_NODES - 1) as f64;
    let mesh = (0..MAX_DWR_MESH_NODES)
        .map(|node| node as f64 / denominator)
        .collect::<Vec<_>>();
    let mut coefficients = vec![0.0; MAX_DWR_POLY_COEFFICIENTS];
    coefficients[1] = -1.0;
    *coefficients
        .last_mut()
        .expect("the admitted coefficient maximum is non-zero") = 1.0;
    let maximum_name = "m".repeat(MAX_FEM1D_CLASS_NAME_BYTES);
    let problem = MmsProblem::new(&maximum_name, poly(coefficients), mesh)
        .expect("the hostile maximum public DWR shape is admitted");
    let candidate = vec![0.0; MAX_DWR_MESH_NODES];

    let cancelled = with_cx(
        true,
        ExecMode::Deterministic,
        Budget::INFINITE,
        default_stream(),
        |cx| dwr_integral_qoi_with_cx(&problem, &candidate, 0.25, 0.75, cx, &VirtualClock::new()),
    );
    // 45_004_590 here vs the in-file 45_004_592 MAX-envelope preflight:
    // this fixture's canonical identity sits 2 bytes under the envelope
    // cap; both values are deterministic const arithmetic.
    assert!(matches!(
        cancelled,
        Err(DwrError::Cancelled {
            phase: "dwr.initial",
            completed_work_units: 0,
            planned_work_units,
        }) if planned_work_units == 45_004_590
            && planned_work_units
                <= u128::try_from(MAX_DWR_WORK_UNITS).expect("usize widens to u128")
    ));
}

#[test]
fn g4_dwr_and_accept_cancellation_are_bounded_and_retryable() {
    let problem = quartic_problem(16);
    let candidate = solve_p1(&problem).expect("G4 primal fixture must solve");
    let baseline =
        dwr_integral_qoi(&problem, &candidate, 0.25, 0.75).expect("healthy DWR baseline");

    let pre_cancelled = with_cx(
        true,
        ExecMode::Deterministic,
        Budget::INFINITE,
        default_stream(),
        |cx| dwr_integral_qoi_with_cx(&problem, &candidate, 0.25, 0.75, cx, &VirtualClock::new()),
    );
    assert!(matches!(
        pre_cancelled,
        Err(DwrError::Cancelled {
            phase: "dwr.initial",
            completed_work_units: 0,
            planned_work_units,
        }) if planned_work_units > 0
    ));

    #[allow(clippy::cast_precision_loss)]
    let stride_mesh: Vec<f64> = (0..=257).map(|node| node as f64 / 257.0).collect();
    let stride_problem = admitted_problem("g4-stride", vec![0.0], stride_mesh);
    let stride_candidate = vec![0.0; 258];
    let stride_cancelled = with_cx(
        false,
        ExecMode::Deterministic,
        Budget::INFINITE.with_poll_quota(3),
        default_stream(),
        |cx| {
            dwr_integral_qoi_with_cx(
                &stride_problem,
                &stride_candidate,
                0.25,
                0.75,
                cx,
                &VirtualClock::new(),
            )
        },
    );
    assert!(matches!(
        stride_cancelled,
        Err(DwrError::BudgetRefused {
            refusal: BudgetRefusal::PollsExhausted {
                phase: "dwr.validate-candidate",
                quota: 3,
            },
            completed_work_units: 257,
            planned_work_units,
        }) if planned_work_units > 257
    ));

    let mut dwr_phases = std::collections::BTreeSet::new();
    let mut reached_success = false;
    for quota in 0..512 {
        let attempt = with_cx(
            false,
            ExecMode::Deterministic,
            Budget::INFINITE.with_poll_quota(quota),
            default_stream(),
            |cx| {
                dwr_integral_qoi_with_cx(&problem, &candidate, 0.25, 0.75, cx, &VirtualClock::new())
            },
        );
        match attempt {
            Err(DwrError::BudgetRefused {
                refusal: BudgetRefusal::PollsExhausted { phase, .. },
                completed_work_units,
                planned_work_units,
            }) => {
                assert!(completed_work_units <= planned_work_units);
                dwr_phases.insert(phase);
                let retry = dwr_integral_qoi(&problem, &candidate, 0.25, 0.75)
                    .expect("budget-refused DWR leaves no retained state");
                assert_same_dwr_semantics(&baseline, &retry);
            }
            Ok(output) => {
                assert_same_dwr_semantics(&baseline, &output);
                reached_success = true;
                break;
            }
            Err(error) => panic!("quota sweep produced a non-cancellation refusal: {error}"),
        }
    }
    assert!(reached_success, "finite healthy quota must complete DWR");
    for phase in [
        "dwr.dual-assembly",
        "dwr.thomas-forward",
        "dwr.residual",
        "dwr.identity",
        "dwr.publish",
    ] {
        assert!(dwr_phases.contains(phase), "quota sweep missed {phase}");
    }

    let query = DwrQuery {
        qoi: "accept-finalization-".repeat(20),
        schema_version: DWR_QUERY_SCHEMA_VERSION,
        window: (0.25, 0.75),
        endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
        semantics: QoiSemantics::WindowedIntegral,
        dims: [0; 6],
        tolerance: 1.0,
    };
    let baseline_accept = accept(&query, &baseline, None);
    let pre_cancelled_accept = with_cx(
        true,
        ExecMode::Deterministic,
        Budget::INFINITE,
        default_stream(),
        |cx| accept_with_cx(&query, &baseline, None, cx, &VirtualClock::new()),
    );
    assert!(matches!(
        pre_cancelled_accept,
        Err(DwrError::Cancelled {
            phase: "dwr-accept.initial",
            completed_work_units: 0,
            planned_work_units,
        }) if planned_work_units
            == u128::try_from(query.qoi.len()).expect("bounded QoI length") + 2
    ));

    let mut accept_phases = std::collections::BTreeSet::new();
    let mut accept_reached_success = false;
    for quota in 0..32 {
        let attempt = with_cx(
            false,
            ExecMode::Deterministic,
            Budget::INFINITE.with_poll_quota(quota),
            default_stream(),
            |cx| accept_with_cx(&query, &baseline, None, cx, &VirtualClock::new()),
        );
        match attempt {
            Err(DwrError::BudgetRefused {
                refusal: BudgetRefusal::PollsExhausted { phase, .. },
                completed_work_units,
                planned_work_units,
            }) => {
                assert!(completed_work_units <= planned_work_units);
                accept_phases.insert(phase);
                let retry = accept(&query, &baseline, None);
                assert_same_accept_semantics(&baseline_accept, &retry);
            }
            Ok(outcome) => {
                assert_same_accept_semantics(&baseline_accept, &outcome);
                accept_reached_success = true;
                break;
            }
            Err(error) => panic!("accept quota sweep produced a non-cancellation refusal: {error}"),
        }
    }
    assert!(
        accept_reached_success,
        "finite healthy quota must complete accept"
    );
    assert!(accept_phases.contains("dwr-accept.identity"));
    assert!(accept_phases.contains("dwr-accept.publish"));
}

#[test]
#[allow(clippy::cast_precision_loss)]
fn g4_identity_hashing_resets_the_poll_stride_between_variable_streams() {
    let mesh: Vec<f64> = (0..=257).map(|node| node as f64 / 257.0).collect();
    let problem = admitted_problem("g4-identity-stride", vec![0.0], mesh);
    let candidate = vec![0.0; 258];
    assert!(problem.canonical_bytes().len() > DWR_POLL_STRIDE_ITEMS);
    assert!(candidate.len() > DWR_POLL_STRIDE_ITEMS);

    let mut identity_positions = Vec::new();
    let mut completed = false;
    for quota in 0..512 {
        let attempt = with_cx(
            false,
            ExecMode::Deterministic,
            Budget::INFINITE.with_poll_quota(quota),
            default_stream(),
            |cx| {
                dwr_integral_qoi_with_cx(&problem, &candidate, 0.25, 0.75, cx, &VirtualClock::new())
            },
        );
        match attempt {
            Err(DwrError::BudgetRefused {
                refusal:
                    BudgetRefusal::PollsExhausted {
                        phase: "dwr.identity",
                        ..
                    },
                completed_work_units,
                ..
            }) => identity_positions.push(completed_work_units),
            Err(DwrError::BudgetRefused {
                refusal: BudgetRefusal::PollsExhausted { .. },
                ..
            }) => {}
            Ok(_) => {
                completed = true;
                break;
            }
            Err(error) => panic!("identity-stride sweep produced an unexpected refusal: {error}"),
        }
    }
    assert!(
        completed,
        "a finite quota must complete the identity-stride fixture"
    );
    identity_positions.sort_unstable();
    identity_positions.dedup();
    assert!(
        identity_positions.len() >= 3,
        "the fixture must exercise multiple variable-length identity chunks"
    );
    for positions in identity_positions.windows(2) {
        assert!(
            positions[1] - positions[0] <= DWR_POLL_STRIDE_ITEMS as u128,
            "identity work advanced {} items between checkpoints",
            positions[1] - positions[0]
        );
    }
}

#[test]
fn g4_bracket_cancellation_drains_each_nested_verifier_phase_and_is_retryable() {
    let problem = quartic_problem(16);
    let candidate = solve_p1(&problem).expect("G4 bracket fixture must solve");
    let baseline = cauchy_schwarz(&problem, &candidate, &problem, &candidate)
        .expect("healthy bracket baseline");

    let pre_cancelled = with_cx(
        true,
        ExecMode::Deterministic,
        Budget::INFINITE,
        default_stream(),
        |cx| {
            Bracket::cauchy_schwarz(
                &problem,
                &candidate,
                &problem,
                &candidate,
                cx,
                &VirtualClock::new(),
            )
        },
    );
    assert!(matches!(
        pre_cancelled,
        Err(BracketError::Cancelled {
            phase: "dwr-bracket.initial",
            completed_work_units: 0,
            planned_work_units,
        }) if planned_work_units > 0
    ));

    let mut phases = std::collections::BTreeSet::new();
    let mut reached_success = false;
    for quota in 0..64 {
        let attempt = with_cx(
            false,
            ExecMode::Deterministic,
            Budget::INFINITE.with_poll_quota(quota),
            default_stream(),
            |cx| {
                Bracket::cauchy_schwarz(
                    &problem,
                    &candidate,
                    &problem,
                    &candidate,
                    cx,
                    &VirtualClock::new(),
                )
            },
        );
        match attempt {
            Err(BracketError::BudgetRefused {
                refusal: BudgetRefusal::PollsExhausted { phase, .. },
                completed_work_units,
                planned_work_units,
            }) => {
                assert!(completed_work_units <= planned_work_units);
                phases.insert(phase);
                let retry = cauchy_schwarz(&problem, &candidate, &problem, &candidate)
                    .expect("budget-refused nested verifier retains no partial bracket");
                assert_same_bracket_semantics(&baseline, &retry);
                assert_eq!(baseline.evidence_identity(), retry.evidence_identity());
            }
            Ok(bracket) => {
                assert_same_bracket_semantics(&baseline, &bracket);
                assert_ne!(baseline.evidence_identity(), bracket.evidence_identity());
                let replay = with_cx(
                    false,
                    ExecMode::Deterministic,
                    Budget::INFINITE.with_poll_quota(quota),
                    default_stream(),
                    |cx| {
                        Bracket::cauchy_schwarz(
                            &problem,
                            &candidate,
                            &problem,
                            &candidate,
                            cx,
                            &VirtualClock::new(),
                        )
                    },
                )
                .expect("the same finite healthy quota must replay");
                assert_same_bracket_semantics(&bracket, &replay);
                assert_eq!(bracket.evidence_identity(), replay.evidence_identity());
                reached_success = true;
                break;
            }
            Err(error) => panic!("bracket quota sweep produced a non-cancellation error: {error}"),
        }
    }
    assert!(
        reached_success,
        "finite healthy quota must complete bracket"
    );
    for phase in [
        "dwr-bracket.primal-verifier.validation",
        "dwr-bracket.primal-verifier.tightness",
        "dwr-bracket.primal-verifier.equilibrated",
        "dwr-bracket.primal-verifier.hash",
        "dwr-bracket.primal-verifier.finalization",
        "dwr-bracket.dual-verifier.validation",
        "dwr-bracket.dual-verifier.tightness",
        "dwr-bracket.dual-verifier.equilibrated",
        "dwr-bracket.dual-verifier.hash",
        "dwr-bracket.dual-verifier.finalization",
        "dwr-bracket.identity",
        "dwr-bracket.publish",
    ] {
        assert!(phases.contains(phase), "quota sweep missed {phase}");
    }
}

#[test]
fn g4_bracket_cancels_at_nested_verifier_global_work_boundary() {
    let problem = quartic_problem(64);
    let candidate = solve_p1(&problem).expect("G4 boundary fixture must solve");
    let verifier_plan = fs_verify::estimator::VerifierWorkPlan::for_inputs(&problem, &candidate)
        .expect("G4 boundary fixture has an admitted verifier plan");
    // fs-verify's 4aeaa7f (replay-admitted production receipts) grew the
    // verifier plan's tightness and total components; the plan identity is
    // fs-verify's authority, re-pinned at batch-verify.
    assert_eq!(verifier_plan.identity_fields(), [229, 64, 64, 144, 1, 502]);

    let baseline = cauchy_schwarz(&problem, &candidate, &problem, &candidate)
        .expect("healthy >256-work bracket baseline");
    // Six successful polls reach the primal tightness phase entry. The seventh
    // callback is the nested verifier's invocation-global 256-unit boundary.
    let cancelled = with_cx(
        false,
        ExecMode::Deterministic,
        Budget::INFINITE.with_poll_quota(6),
        default_stream(),
        |cx| {
            Bracket::cauchy_schwarz(
                &problem,
                &candidate,
                &problem,
                &candidate,
                cx,
                &VirtualClock::new(),
            )
        },
    );
    assert!(matches!(
        cancelled,
        Err(BracketError::BudgetRefused {
            refusal: BudgetRefusal::PollsExhausted {
                phase: "dwr-bracket.primal-verifier.tightness",
                quota: 6,
            },
            completed_work_units: 450,
            planned_work_units,
        }) if planned_work_units > 450
    ));

    let retry = cauchy_schwarz(&problem, &candidate, &problem, &candidate)
        .expect("nested boundary cancellation retains no partial bracket");
    assert_same_bracket_semantics(&baseline, &retry);
    assert_eq!(baseline.evidence_identity(), retry.evidence_identity());
}

#[test]
fn g5_bracket_identity_binds_execution_and_complete_nested_work_policy() {
    let problem = quartic_problem(8);
    let candidate = solve_p1(&problem).expect("G5 bracket fixture must solve");
    let run = |mode: ExecMode, budget: Budget, stream: StreamKey| {
        with_cx(false, mode, budget, stream, |cx| {
            Bracket::cauchy_schwarz(
                &problem,
                &candidate,
                &problem,
                &candidate,
                cx,
                &VirtualClock::new(),
            )
            .expect("G5 bracket execution")
        })
    };
    let stream = default_stream();
    let baseline = run(ExecMode::Deterministic, Budget::INFINITE, stream);
    let repeat = run(ExecMode::Deterministic, Budget::INFINITE, stream);
    assert_same_bracket_semantics(&baseline, &repeat);
    assert_eq!(baseline.evidence_identity(), repeat.evidence_identity());

    let variants = [
        ("mode", ExecMode::Fast, Budget::INFINITE, stream),
        (
            "deadline",
            ExecMode::Deterministic,
            Budget {
                deadline: Budget::with_deadline_at_ns(123).deadline,
                ..Budget::INFINITE
            },
            stream,
        ),
        (
            "poll quota",
            ExecMode::Deterministic,
            Budget::INFINITE.with_poll_quota(10_000),
            stream,
        ),
        (
            "cost quota",
            ExecMode::Deterministic,
            Budget::INFINITE.with_cost_quota(9_000_000),
            stream,
        ),
        (
            "priority",
            ExecMode::Deterministic,
            Budget::INFINITE.with_priority(7),
            stream,
        ),
        (
            "stream",
            ExecMode::Deterministic,
            Budget::INFINITE,
            StreamKey {
                iteration: stream.iteration + 1,
                ..stream
            },
        ),
    ];
    for (field, mode, budget, stream) in variants {
        let changed = run(mode, budget, stream);
        assert_same_bracket_semantics(&baseline, &changed);
        assert_ne!(
            baseline.evidence_identity(),
            changed.evidence_identity(),
            "bracket identity omitted {field}"
        );
    }

    let coarse = admitted_problem("g5-bracket-shape", vec![0.0], vec![0.0, 0.5, 1.0]);
    let refined = admitted_problem("g5-bracket-shape", vec![0.0], vec![0.0, 0.25, 0.5, 1.0]);
    let coarse_bracket =
        cauchy_schwarz(&coarse, &[0.0; 3], &coarse, &[0.0; 3]).expect("coarse zero bracket");
    let refined_bracket =
        cauchy_schwarz(&refined, &[0.0; 4], &refined, &[0.0; 4]).expect("refined zero bracket");
    for bracket in [&coarse_bracket, &refined_bracket] {
        assert!(bracket.bound().is_finite());
        assert!(bracket.bound() >= 0.0);
    }
    assert_ne!(
        coarse_bracket.evidence_identity(),
        refined_bracket.evidence_identity(),
        "nested verifier work shape and exact problem inputs are semantic"
    );
    assert_eq!(fs_verify::estimator::VERIFIER_WORK_PLAN_VERSION, 2);
    assert_eq!(fs_verify::estimator::VERIFIER_POLL_POLICY_VERSION, 1);
    assert_eq!(fs_verify::estimator::VERIFIER_POLL_STRIDE_WORK_UNITS, 256);
}

#[test]
fn g5_execution_identity_binds_mode_budget_stream_and_work_shape() {
    let problem = quartic_problem(12);
    let candidate = solve_p1(&problem).expect("G5 primal fixture must solve");
    let query = DwrQuery {
        qoi: "g5-integral".to_string(),
        schema_version: DWR_QUERY_SCHEMA_VERSION,
        window: (0.25, 0.75),
        endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
        semantics: QoiSemantics::WindowedIntegral,
        dims: [0; 6],
        tolerance: 1.0,
    };
    let run = |mode: ExecMode, budget: Budget, stream: StreamKey| {
        with_cx(false, mode, budget, stream, |cx| {
            let output = dwr_integral_qoi_with_cx(
                &problem,
                &candidate,
                0.25,
                0.75,
                cx,
                &VirtualClock::new(),
            )
            .expect("G5 execution must remain scientifically valid");
            let outcome = accept_with_cx(&query, &output, None, cx, &VirtualClock::new())
                .expect("G5 acceptance must remain valid");
            (output, outcome)
        })
    };

    let stream = default_stream();
    let (baseline_output, baseline_accept) = run(ExecMode::Deterministic, Budget::INFINITE, stream);
    let (repeat_output, repeat_accept) = run(ExecMode::Deterministic, Budget::INFINITE, stream);
    assert_same_dwr_semantics(&baseline_output, &repeat_output);
    assert_same_accept_semantics(&baseline_accept, &repeat_accept);
    assert_eq!(
        baseline_output.evidence_identity(),
        repeat_output.evidence_identity()
    );
    assert_eq!(
        baseline_accept.evidence_identity(),
        repeat_accept.evidence_identity()
    );

    let variants = [
        ("mode", ExecMode::Fast, Budget::INFINITE, stream),
        (
            "deadline",
            ExecMode::Deterministic,
            Budget {
                deadline: Budget::with_deadline_at_ns(123).deadline,
                ..Budget::INFINITE
            },
            stream,
        ),
        (
            "poll quota",
            ExecMode::Deterministic,
            Budget::INFINITE.with_poll_quota(10_000),
            stream,
        ),
        (
            "cost quota",
            ExecMode::Deterministic,
            Budget::INFINITE.with_cost_quota(9_000_000),
            stream,
        ),
        (
            "priority",
            ExecMode::Deterministic,
            Budget::INFINITE.with_priority(7),
            stream,
        ),
        (
            "seed",
            ExecMode::Deterministic,
            Budget::INFINITE,
            StreamKey {
                seed: stream.seed + 1,
                ..stream
            },
        ),
        (
            "kernel",
            ExecMode::Deterministic,
            Budget::INFINITE,
            StreamKey {
                kernel_id: stream.kernel_id + 1,
                ..stream
            },
        ),
        (
            "tile",
            ExecMode::Deterministic,
            Budget::INFINITE,
            StreamKey {
                tile: stream.tile + 1,
                ..stream
            },
        ),
        (
            "iteration",
            ExecMode::Deterministic,
            Budget::INFINITE,
            StreamKey {
                iteration: stream.iteration + 1,
                ..stream
            },
        ),
    ];
    for (field, mode, budget, stream) in variants {
        let (output, outcome) = run(mode, budget, stream);
        assert_same_dwr_semantics(&baseline_output, &output);
        assert_same_accept_semantics(&baseline_accept, &outcome);
        assert_ne!(
            baseline_output.evidence_identity(),
            output.evidence_identity(),
            "DWR identity omitted {field}"
        );
        assert_ne!(
            baseline_accept.evidence_identity(),
            outcome.evidence_identity(),
            "accept identity omitted {field}"
        );
    }

    let renamed_problem = admitted_problem(
        "g5-renamed-but-numerically-identical",
        vec![0.0, 1.0, -1.0, 2.0, -2.0],
        problem.mesh().to_vec(),
    );
    let renamed_output = dwr_integral_qoi(&renamed_problem, &candidate, 0.25, 0.75)
        .expect("renaming preserves the numerical problem");
    assert_same_dwr_semantics(&baseline_output, &renamed_output);
    assert_ne!(
        baseline_output.evidence_identity(),
        renamed_output.evidence_identity(),
        "exact canonical problem bytes are semantic"
    );

    let mut changed_candidate = candidate.clone();
    changed_candidate[candidate.len() / 2] += f64::EPSILON;
    let candidate_output = dwr_integral_qoi(&problem, &changed_candidate, 0.25, 0.75)
        .expect("finite interior candidate mutation remains admitted");
    assert_ne!(
        baseline_output.evidence_identity(),
        candidate_output.evidence_identity(),
        "candidate values are semantic"
    );
    let window_output = dwr_integral_qoi(&problem, &candidate, 0.2, 0.8)
        .expect("alternate finite QoI window remains admitted");
    assert_ne!(
        baseline_output.evidence_identity(),
        window_output.evidence_identity(),
        "QoI window is semantic"
    );

    let same_length_qoi = accept(
        &DwrQuery {
            qoi: "g5-integrax".to_string(),
            schema_version: DWR_QUERY_SCHEMA_VERSION,
            window: (0.25, 0.75),
            endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
            semantics: QoiSemantics::WindowedIntegral,
            dims: [0; 6],
            tolerance: query.tolerance,
        },
        &baseline_output,
        None,
    );
    assert_same_accept_semantics(&baseline_accept, &same_length_qoi);
    assert_ne!(
        baseline_accept.evidence_identity(),
        same_length_qoi.evidence_identity(),
        "QoI label content is semantic independently of work shape"
    );
    let changed_tolerance = accept(
        &DwrQuery {
            qoi: query.qoi.clone(),
            schema_version: DWR_QUERY_SCHEMA_VERSION,
            window: (0.25, 0.75),
            endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
            semantics: QoiSemantics::WindowedIntegral,
            dims: [0; 6],
            tolerance: 2.0,
        },
        &baseline_output,
        None,
    );
    assert_ne!(
        baseline_accept.evidence_identity(),
        changed_tolerance.evidence_identity(),
        "tolerance is semantic"
    );
    // A shifted scalar is unforgeable under the sealed accept (sj31i.1):
    // a genuinely different estimate (here, a different retained problem
    // and mesh) must move the acceptance identity instead.
    let changed_estimate = accept(&query, &solved_estimate(), None);
    assert_ne!(
        baseline_accept.evidence_identity(),
        changed_estimate.evidence_identity(),
        "DWR estimate is semantic"
    );

    let short = accept(
        &DwrQuery {
            qoi: "q".to_string(),
            schema_version: DWR_QUERY_SCHEMA_VERSION,
            window: (0.25, 0.75),
            endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
            semantics: QoiSemantics::WindowedIntegral,
            dims: [0; 6],
            tolerance: query.tolerance,
        },
        &baseline_output,
        None,
    );
    let long = accept(
        &DwrQuery {
            qoi: "same-science-longer-provenance-label".to_string(),
            schema_version: DWR_QUERY_SCHEMA_VERSION,
            window: (0.25, 0.75),
            endpoint_convention: EndpointConvention::HomogeneousDirichletZero,
            semantics: QoiSemantics::WindowedIntegral,
            dims: [0; 6],
            tolerance: query.tolerance,
        },
        &baseline_output,
        None,
    );
    assert_same_accept_semantics(&short, &long);
    assert_ne!(short.evidence_identity(), long.evidence_identity());
    assert_eq!(DWR_WORK_PLAN_VERSION, 2);
    assert_eq!(DWR_POLL_POLICY_VERSION, 3);
    assert_eq!(DWR_EVIDENCE_IDENTITY_VERSION, 7);
    assert_eq!(DWR_QUERY_SCHEMA_VERSION, 1);
    assert_eq!(DWR_POLL_STRIDE_ITEMS, 256);
}
