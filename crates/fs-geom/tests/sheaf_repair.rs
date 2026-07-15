//! Sheaf-repair conformance (the wqd.14 bead; runs under the
//! `sheaf-repair` feature). Acceptance: exact-component defects
//! produce a budget-eligible gauge step and a passing nominal residual;
//! coexact seeding retains a localized, explicitly non-causal circulation
//! diagnostic;
//! harmonic seeding retains a closed, non-exact witness and its full
//! interface support; predicted post-repair norms match actuals; the retained
//! first-fit fixture matches a dense reference; converged re-planning and
//! repair are budget-safe. These fixtures do not confer generic convergence or
//! orthogonality authority on the fixed-iteration routine.
#![cfg(feature = "sheaf-repair")]

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::router::{ConverterSpec, ErrorModel, MemoryCostOracle, RouteRequest, Router};
use fs_geom::sheaf::{Interface, SheafComplex};
use fs_geom::sheaf_repair::{
    AdmittedSheafSkeleton, COMPONENT_FLOOR, SheafRepairBudget, SheafRepairError, SheafSkeleton,
    SheafSkeletonError, apply_gauge, hodge_decompose, hodge_decompose_bounded, plan_repair,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-geom/sheaf-repair\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn with_gate_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5348_4541_4652_4550,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gate_cx(&gate, f)
}

/// A 3-patch triangle complex (one triple junction).
fn triangle() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (1, 2), (0, 2)],
        triangles: vec![(0, 1, 2)],
    }
}

/// A 4-patch ring (cycle, NO triangles): H¹ is nontrivial by design.
fn ring() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 4,
        edges: vec![(0, 1), (1, 2), (2, 3), (0, 3)],
        triangles: vec![],
    }
}

fn norm_inf(v: &[f64]) -> f64 {
    assert!(
        v.iter().all(|value| value.is_finite()),
        "test norm requires finite values: {v:?}"
    );
    v.iter().fold(0.0f64, |a, &b| a.max(b.abs()))
}

/// Dense least-squares oracle: minimize ‖m − A x‖² by normal equations
/// solved with Gaussian elimination (partial pivot) — an independent
/// code path from the module's Gauss–Seidel.
fn dense_projection(m: &[f64], columns: &[Vec<f64>]) -> Vec<f64> {
    let n = columns.len();
    let mut ata = vec![vec![0.0f64; n]; n];
    let mut atb = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            ata[i][j] = columns[i].iter().zip(&columns[j]).map(|(a, b)| a * b).sum();
        }
        atb[i] = columns[i].iter().zip(m).map(|(a, b)| a * b).sum();
    }
    // Ridge the (rank-deficient) gauge direction minimally.
    for (i, row) in ata.iter_mut().enumerate() {
        row[i] += 1e-12;
    }
    // Gaussian elimination.
    let mut aug: Vec<Vec<f64>> = ata
        .iter()
        .zip(&atb)
        .map(|(row, &b)| {
            let mut r = row.clone();
            r.push(b);
            r
        })
        .collect();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&a, &b| aug[a][col].abs().total_cmp(&aug[b][col].abs()))
            .expect("rows");
        aug.swap(col, pivot);
        let p = aug[col][col];
        if p.abs() < 1e-300 {
            continue;
        }
        for r in 0..n {
            if r != col {
                let f = aug[r][col] / p;
                let pivot_row = aug[col].clone();
                for (k, cell) in aug[r].iter_mut().enumerate().skip(col) {
                    *cell -= f * pivot_row[k];
                }
            }
        }
    }
    (0..n)
        .map(|i| {
            if aug[i][i].abs() < 1e-300 {
                0.0
            } else {
                aug[i][n] / aug[i][i]
            }
        })
        .collect()
}

#[test]
fn sr_001_sequential_fit_matches_dense_fixture() {
    let sk = triangle();
    // A mixed cochain: gauge part + circulation part.
    let gauge_part = sk.d0(&[0.0, 0.7, -0.3]);
    let circ_part = sk.d1t(&[0.4]);
    let m: Vec<f64> = gauge_part
        .iter()
        .zip(&circ_part)
        .map(|(a, b)| a + b)
        .collect();
    let split = hodge_decompose(&sk, &m);
    // Oracle: dense projection onto im δ⁰ (columns = δ⁰ of unit vertex
    // vectors, vertex 0 pinned) and im δ¹ᵀ.
    let d0_cols: Vec<Vec<f64>> = (1..sk.n_patches)
        .map(|i| {
            let mut e = vec![0.0; sk.n_patches];
            e[i] = 1.0;
            sk.d0(&e)
        })
        .collect();
    let c_oracle = dense_projection(&m, &d0_cols);
    let exact_oracle = {
        let mut full = vec![0.0; sk.n_patches];
        full[1..].copy_from_slice(&c_oracle);
        sk.d0(&full)
    };
    for (got, want) in split.exact.iter().zip(&exact_oracle) {
        assert!(
            (got - want).abs() < 1e-8,
            "exact component vs dense oracle: {got} vs {want}"
        );
    }
    // On a triangle (contractible), harmonic must vanish.
    assert!(
        norm_inf(&split.harmonic) < 1e-8,
        "contractible complex has no harmonic part: {:?}",
        split.harmonic
    );
    // Orthogonality residuals: δ⁰ᵀh ≈ 0 and δ¹h ≈ 0.
    assert!(norm_inf(&sk.d0t(&split.harmonic)) < 1e-8);
    assert!(norm_inf(&sk.d1(&split.harmonic)) < 1e-8);
    // This fixture also verifies near-orthogonality, so its diagnostic ratios
    // sum to approximately one. The API does not claim that law generically.
    let (fe, fc, fh) = split.fractions;
    assert!(
        (fe + fc + fh - 1.0).abs() < 1e-6,
        "this verified fixture approximately partitions energy: {fe} + {fc} + {fh}"
    );
    verdict(
        "sr-001",
        "the fitted exact component matches the dense-reference projection; the \
         fixture remainder vanishes and its checked ratios approximately sum to one",
    );
}

#[test]
fn sr_002_exact_defect_auto_repairs_within_budget() {
    let sk = triangle();
    // Seed a pure gauge defect: patch 2 drifted by +0.012.
    let mismatch = sk.d0(&[0.0, 0.0, 0.012]);
    let budgets = [0.02, 0.02, 0.02];
    let plan = plan_repair(&sk, &mismatch, &budgets, None);
    assert!(
        plan.gauge_step_eligible,
        "within budgets: gauge step eligible"
    );
    assert!(plan.split.fractions.0 > 0.999, "pure exact defect");
    assert!(plan.harmonic_support.is_empty(), "no harmonic remainder");
    // Predicted-vs-actual: apply the gauge, re-measure.
    let predicted = plan.proposals[0].expected_post_norm;
    let repaired = apply_gauge(&sk, &mismatch, &plan.gauge);
    let actual = norm_inf(&repaired);
    assert!(
        (predicted - actual).abs() < 1e-9,
        "prediction {predicted} vs actual {actual}"
    );
    assert!(actual < 1e-9, "nominal residual passes after repair");
    // Repair SAFETY: offsets stay within each chart's declared budget.
    for (off, b) in plan.gauge.iter().zip(&budgets) {
        assert!(off.abs() <= *b, "repair never exceeds a budget");
    }
    // Converged re-planning: planning from the repaired state yields a
    // near-zero follow-up gauge. Applying the original nonzero gauge twice is
    // deliberately not claimed to be idempotent.
    let plan2 = plan_repair(&sk, &repaired, &budgets, None);
    assert!(
        norm_inf(&plan2.gauge) < 1e-9,
        "no residual gauge on a passing model: {:?}",
        plan2.gauge
    );
    let repaired2 = apply_gauge(&sk, &repaired, &plan2.gauge);
    assert!(
        (norm_inf(&repaired2) - actual).abs() < 1e-12,
        "no-op repair"
    );
    // Over-budget variant: the SAME defect with a tight budget must NOT
    // auto-apply (needs explicit acceptance).
    let tight = [0.001, 0.001, 0.001];
    let gated = plan_repair(&sk, &mismatch, &tight, None);
    assert!(
        !gated.gauge_step_eligible,
        "budget gate blocks silent distortion"
    );
    assert!(
        gated.proposals[0].action.contains("EXCEEDS"),
        "the proposal says so: {}",
        gated.proposals[0].action
    );
    verdict(
        "sr-002",
        "gauge defect repaired to ~0 with exact prediction; converged re-planning is a \
         no-op; budget gate blocks over-budget auto-apply",
    );
}

#[test]
fn sr_002a_component_gauge_shift_finds_feasible_budget_representative() {
    let sk = SheafSkeleton {
        n_patches: 4,
        edges: vec![(0, 1), (2, 3)],
        triangles: Vec::new(),
    };
    // The pinned least-squares representative is [0,2,0,4], which appears to
    // violate the budgets. Independent component shifts produce [-1,1,-2,2]
    // without changing either coboundary correction.
    let mismatch = vec![2.0, 4.0];
    let budgets = [1.0, 1.0, 2.0, 2.0];
    let plan = plan_repair(&sk, &mismatch, &budgets, None);
    assert!(
        plan.gauge_step_eligible,
        "a feasible gauge representative exists"
    );
    assert_eq!(plan.gauge, vec![-1.0, 1.0, -2.0, 2.0]);
    assert_eq!(sk.d0(&plan.gauge), mismatch);
    assert_eq!(apply_gauge(&sk, &mismatch, &plan.gauge), vec![0.0; 2]);

    let impossible = plan_repair(&sk, &mismatch, &[0.9, 0.9, 1.9, 1.9], None);
    assert!(
        !impossible.gauge_step_eligible,
        "a component difference larger than the sum of its patch budgets must refuse auto-apply"
    );

    let one_edge = SheafSkeleton {
        n_patches: 2,
        edges: vec![(0, 1)],
        triangles: Vec::new(),
    };
    let slack = plan_repair(&one_edge, &[-2.0], &[f64::INFINITY, 1.0], None);
    assert_eq!(
        slack.gauge,
        vec![2.0, 0.0],
        "feasible shift interval [1,3] uses its deterministic maximum-slack midpoint"
    );
    assert_eq!(one_edge.d0(&slack.gauge), vec![-2.0]);

    let centered = plan_repair(&one_edge, &[2.0], &[100.0, 100.0], None);
    assert_eq!(
        centered.gauge,
        vec![-1.0, 1.0],
        "feasible interval [-100,98] uses its maximum-slack midpoint, not zero"
    );
    assert_eq!(one_edge.d0(&centered.gauge), vec![2.0]);
}

#[test]
fn sr_002aa_skeleton_extraction_refuses_unvalidated_public_complex() {
    let malformed = SheafComplex {
        n_patches: 2,
        interfaces: vec![Interface {
            patches: (1, 0),
            samples: Vec::new(),
        }],
        triples: Vec::new(),
        sampling_clip: None,
    };
    assert!(
        SheafSkeleton::of(&malformed).is_err(),
        "public malformed indices must not be copied into later incidence panics"
    );
}

#[test]
#[should_panic(expected = "one gauge budget per patch")]
fn sr_002b_short_budget_vector_hits_documented_precondition_guard() {
    // Regression: a budgets slice shorter than n_patches was silently truncated
    // by `potential.iter().zip(budgets)`, leaving trailing patches unchecked so
    // `gauge_step_eligible` could bless an over-budget distortion — the one thing
    // the planner promises never to do. Must fail closed.
    let sk = triangle(); // n_patches = 3
    let mismatch = vec![0.0; sk.edges.len()];
    let _ = plan_repair(&sk, &mismatch, &[1.0, 1.0], None); // only 2 budgets
}

#[test]
fn sr_003_coexact_seeding_retains_noncausal_diagnostic() {
    let sk = triangle();
    // Seed a pure circulation (the flipped-orientation signature): the
    // image of δ¹ᵀ.
    let mismatch = sk.d1t(&[0.05]);
    let plan = plan_repair(&sk, &mismatch, &[1.0; 3], None);
    assert!(
        plan.split.fractions.1 > 0.999,
        "pure coexact defect: {:?}",
        plan.split.fractions
    );
    let circulation_proposal = plan
        .proposals
        .iter()
        .find(|p| p.action.contains("coexact circulation candidate"))
        .expect("coexact circulation diagnostic present");
    assert!(
        circulation_proposal.action.contains("(0, 1, 2)"),
        "localized to the triple junction: {}",
        circulation_proposal.action
    );
    assert!(
        circulation_proposal
            .action
            .contains("algebra alone does not assign cause"),
        "diagnostic must not turn one hypothesis into a causal conclusion: {}",
        circulation_proposal.action
    );
    // Gauge repair CANNOT fix circulation: applying it leaves the norm.
    let repaired = apply_gauge(&sk, &mismatch, &plan.gauge);
    assert!(
        norm_inf(&repaired) > 0.9 * norm_inf(&mismatch),
        "circulation is not gauge-repairable"
    );
    verdict(
        "sr-003",
        "circulation seeding is >99.9% coexact, localized at the retained triangle \
         without assigning cause, and this fixture is not gauge-fixable",
    );
}

#[test]
fn sr_004_harmonic_seeding_retains_closed_nonexact_witness() {
    let sk = ring();
    // A circulation around the 4-cycle: with no 2-cells, nothing coexact
    // exists and no gauge kills a loop sum — genuinely harmonic.
    // Orientation: edges (0,1),(1,2),(2,3) run low→high along the loop;
    // (0,3) runs AGAINST it, so the loop cochain is (ε, ε, ε, −ε).
    let eps = 0.03;
    let mismatch = vec![eps, eps, eps, -eps];
    let plan = plan_repair(&sk, &mismatch, &[1.0; 4], None);
    assert!(
        norm_inf(&plan.split.harmonic) > 0.9 * eps,
        "the retained harmonic witness must be nonzero: {:?}",
        plan.split.harmonic
    );
    assert!(
        norm_inf(&sk.d1(&plan.split.harmonic)) < 1e-12,
        "the retained mismatch cochain must be closed"
    );
    assert!(
        norm_inf(&sk.d0t(&plan.split.harmonic)) < 1e-12,
        "the harmonic witness must be orthogonal to every exact cochain"
    );
    // Exact edge cochains telescope to zero around this oriented cycle.
    // A nonzero cycle pairing therefore witnesses that the retained closed
    // cochain is not in im(delta0), which is the extra evidence required
    // before calling an interface mismatch an H1 obstruction.
    let cycle_pairing = plan.split.harmonic[0] + plan.split.harmonic[1] + plan.split.harmonic[2]
        - plan.split.harmonic[3];
    assert!(
        (cycle_pairing - 4.0 * eps).abs() < 1e-12,
        "nonzero cycle pairing witnesses non-exactness: {cycle_pairing}"
    );
    assert!(
        plan.split.fractions.2 > 0.999,
        "pure harmonic: {:?}",
        plan.split.fractions
    );
    assert!(
        !plan.gauge_step_eligible,
        "a retained non-exact fixture must not be auto-repairable"
    );
    assert!(
        plan.split.fractions.0 < 1e-9,
        "the ring fixture must not acquire a material exact component: {:?}",
        plan.split.fractions
    );
    assert_eq!(
        plan.harmonic_support.len(),
        4,
        "the whole cycle is retained as harmonic support"
    );
    let remainder = plan
        .proposals
        .iter()
        .find(|p| p.action.contains("no generic exactness or topology claim"))
        .expect("honest candidate-remainder proposal");
    assert!(
        remainder.cost_s.is_infinite(),
        "no fabricated repair-cost claim"
    );
    // And indeed gauge repair achieves nothing.
    let repaired = apply_gauge(&sk, &mismatch, &plan.gauge);
    assert!(norm_inf(&repaired) > 0.9 * eps, "harmonic survives gauge");
    verdict(
        "sr-004",
        "retained cycle circulation is closed, non-exact by nonzero cycle pairing, \
         >99.9% harmonic, and outside the patch-gauge repair class with full support",
    );
}

#[test]
fn sr_004a_subfloor_remainder_is_retained_but_not_promoted_to_support() {
    let sk = ring();
    let exact = sk.d0(&[0.0, 1.0, 0.0, -1.0]);
    let eps = 1e-8;
    let cycle = [eps, eps, eps, -eps];
    let mismatch: Vec<f64> = exact.iter().zip(cycle).map(|(a, b)| a + b).collect();
    let plan = plan_repair(&sk, &mismatch, &[2.0; 4], None);
    assert!(
        plan.split.harmonic.iter().any(|value| *value != 0.0),
        "the raw diagnostic split retains the nonzero remainder"
    );
    assert!(
        plan.split.fractions.2 <= COMPONENT_FLOOR,
        "fixture must stay below component admission: {:?}",
        plan.split.fractions
    );
    assert!(
        plan.harmonic_support.is_empty(),
        "a component cannot bootstrap significance from its own maximum"
    );
    assert!(
        plan.proposals
            .iter()
            .all(|proposal| !proposal.action.contains("retained harmonic remainder")),
        "sub-floor residue must not create a scary +inf remainder proposal"
    );
}

#[test]
fn sr_005_router_reroute_proposal_ranks_by_expected_norm() {
    let sk = triangle();
    let mismatch = sk.d0(&[0.0, 0.0, 0.012]);
    // A router with one declared conversion available for the worst patch's
    // chart kind. The declaration is not authenticated certificate authority.
    let mut router = Router::new();
    router
        .register(ConverterSpec {
            name: "sdf->mesh/dc-interval".to_string(),
            from: "sdf".to_string(),
            to: "mesh".to_string(),
            base_cost_s: 2.0,
            error: ErrorModel::AdditiveAbs(5e-7),
            certified: true,
        })
        .expect("register");
    let oracle = MemoryCostOracle::new();
    let req = RouteRequest {
        from: "sdf".to_string(),
        to: "mesh".to_string(),
        scale: 1.0,
        max_abs_error: 1e-3,
        max_cost_s: 100.0,
    };
    let plan = plan_repair(&sk, &mismatch, &[1.0; 3], Some((&router, &oracle, &req)));
    let reroute = plan
        .proposals
        .iter()
        .find(|p| p.action.contains("reroute"))
        .expect("router proposal present");
    assert!(reroute.action.contains("dc-interval"), "{}", reroute.action);
    assert!((reroute.cost_s - 2.0).abs() < 1e-9, "router-modeled cost");
    // The route's composed representation error is not a post-repair seam
    // norm. It remains explicitly unavailable until the reroute is
    // constructively applied and re-evaluated.
    assert!(reroute.expected_post_norm.is_infinite());
    // Available constructive predictions sort before unavailable ones.
    for pair in plan.proposals.windows(2) {
        assert!(
            pair[0].expected_post_norm <= pair[1].expected_post_norm + 1e-12,
            "proposals ranked by expected norm"
        );
    }
    verdict(
        "sr-005",
        "router reroute proposal carries the planned route + modeled cost; ranking is \
         by expected post-repair norm",
    );
}

#[test]
fn sr_006_admitted_skeleton_rejects_the_first_structural_error() {
    assert_eq!(
        AdmittedSheafSkeleton::try_new(0, Vec::new(), Vec::new()),
        Err(SheafSkeletonError::EmptyComplex)
    );
    assert_eq!(
        AdmittedSheafSkeleton::try_new(3, vec![(0, 2), (0, 1)], Vec::new()),
        Err(SheafSkeletonError::InvalidEdge { index: 1 })
    );
    assert_eq!(
        AdmittedSheafSkeleton::try_new(3, vec![(0, 1), (0, 1)], Vec::new()),
        Err(SheafSkeletonError::InvalidEdge { index: 1 })
    );
    assert_eq!(
        AdmittedSheafSkeleton::try_new(3, vec![(1, 0)], Vec::new()),
        Err(SheafSkeletonError::InvalidEdge { index: 0 })
    );
    assert_eq!(
        AdmittedSheafSkeleton::try_new(3, vec![(0, 3)], Vec::new()),
        Err(SheafSkeletonError::InvalidEdge { index: 0 })
    );
    assert_eq!(
        AdmittedSheafSkeleton::try_new(3, vec![(0, 1), (0, 2)], vec![(0, 1, 2)],),
        Err(SheafSkeletonError::InvalidTriangle { index: 0 })
    );
    verdict(
        "sr-006",
        "opaque skeleton admission refuses empty, non-canonical, duplicate, reversed, out-of-range, and incomplete incidence deterministically",
    );
}

#[test]
fn sr_007_admitted_incidence_is_total_and_finite() {
    let skeleton = AdmittedSheafSkeleton::try_new(3, vec![(0, 1), (0, 2), (1, 2)], vec![(0, 1, 2)])
        .expect("canonical triangle admits");
    assert_eq!(skeleton.n_patches(), 3);
    assert_eq!(skeleton.edges(), &[(0, 1), (0, 2), (1, 2)]);
    assert_eq!(skeleton.triangles(), &[(0, 1, 2)]);

    assert_eq!(
        skeleton.d0(&[0.0, 1.0]),
        Err(SheafSkeletonError::CochainLength {
            role: "vertex",
            expected: 3,
            actual: 2,
        })
    );
    assert_eq!(
        skeleton.d0(&[0.0, f64::NAN, 1.0]),
        Err(SheafSkeletonError::NonFiniteCochain {
            role: "vertex",
            index: 1,
        })
    );
    assert_eq!(
        skeleton.d0(&[-f64::MAX, f64::MAX, 0.0]),
        Err(SheafSkeletonError::NumericalOverflow { stage: "d0" })
    );
    assert_eq!(
        skeleton.d1(&[1.0, 2.0]),
        Err(SheafSkeletonError::CochainLength {
            role: "edge",
            expected: 3,
            actual: 2,
        })
    );
    assert_eq!(
        skeleton.d1t(&[]),
        Err(SheafSkeletonError::CochainLength {
            role: "triangle",
            expected: 1,
            actual: 0,
        })
    );

    let vertex = [2.0, -3.0, 5.0];
    let edge = skeleton.d0(&vertex).expect("finite coboundary");
    assert_eq!(
        skeleton.d1(&edge).expect("finite boundary composition"),
        vec![0.0],
        "delta-one after delta-zero is exact in the admitted scalar complex"
    );
    assert_eq!(
        skeleton.d0t(&edge).expect("finite transpose"),
        vec![2.0, -13.0, 11.0]
    );
    assert_eq!(
        skeleton.d1t(&[4.0]).expect("finite triangle transpose"),
        vec![4.0, -4.0, 4.0]
    );
    verdict(
        "sr-007",
        "admitted incidence rejects shape, non-finite, and overflow inputs without panic and preserves delta-one composed with delta-zero",
    );
}

#[test]
fn sr_008_bounded_decomposition_accounts_before_work() {
    let skeleton = AdmittedSheafSkeleton::try_new(3, vec![(0, 1), (0, 2), (1, 2)], vec![(0, 1, 2)])
        .expect("canonical triangle admits");
    let mismatch = skeleton
        .d0(&[0.0, 0.5, -0.25])
        .expect("finite exact fixture");
    let budget = SheafRepairBudget {
        sweeps: 8,
        max_operator_evaluations: 56,
        max_scalar_slots: 42,
        poll_stride: 1,
    };
    with_cx(|cx| {
        let bounded = hodge_decompose_bounded(&skeleton, &mismatch, budget, cx)
            .expect("exact admitted schedule fits its envelope");
        assert_eq!(bounded.budget, budget);
        assert_eq!(bounded.usage.completed_sweeps, 16);
        assert_eq!(bounded.usage.operator_evaluations, 56);
        assert_eq!(bounded.usage.admitted_scalar_slots, 42);
        assert!(
            bounded
                .split
                .exact
                .iter()
                .chain(&bounded.split.potential)
                .chain(&bounded.split.coexact)
                .chain(&bounded.split.harmonic)
                .all(|value| value.is_finite())
        );

        assert_eq!(
            hodge_decompose_bounded(
                &skeleton,
                &mismatch,
                SheafRepairBudget {
                    max_operator_evaluations: 55,
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::WorkBudgetExceeded {
                required: 56,
                cap: 55,
            })
        );
        assert_eq!(
            hodge_decompose_bounded(
                &skeleton,
                &mismatch,
                SheafRepairBudget {
                    max_scalar_slots: 41,
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::MemoryBudgetExceeded {
                required: 42,
                cap: 41,
            })
        );
        assert_eq!(
            hodge_decompose_bounded(&skeleton, &[0.0, f64::NAN, 0.0], budget, cx,),
            Err(SheafRepairError::Skeleton(
                SheafSkeletonError::NonFiniteCochain {
                    role: "edge",
                    index: 1,
                }
            ))
        );
        assert_eq!(
            hodge_decompose_bounded(
                &skeleton,
                &mismatch,
                SheafRepairBudget {
                    sweeps: 0,
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::InvalidBudget { field: "sweeps" })
        );
        assert_eq!(
            hodge_decompose_bounded(
                &skeleton,
                &mismatch,
                SheafRepairBudget {
                    poll_stride: 0,
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::InvalidBudget {
                field: "poll_stride",
            })
        );
    });
    verdict(
        "sr-008",
        "bounded decomposition retains exact work and memory consumption and refuses under-budget or non-finite requests before publication",
    );
}

#[test]
fn sr_009_bounded_decomposition_observes_pre_cancellation() {
    let skeleton =
        AdmittedSheafSkeleton::try_new(2, vec![(0, 1)], Vec::new()).expect("single edge admits");
    let budget = SheafRepairBudget {
        sweeps: 2,
        max_operator_evaluations: 32,
        max_scalar_slots: 32,
        poll_stride: 1,
    };
    let gate = CancelGate::new();
    gate.request();
    with_gate_cx(&gate, |cx| {
        assert_eq!(
            hodge_decompose_bounded(&skeleton, &[1.0], budget, cx),
            Err(SheafRepairError::Cancelled {
                stage: "admission",
                completed_sweeps: 0,
                operator_evaluations: 0,
            })
        );
    });
    verdict(
        "sr-009",
        "pre-cancelled diagnostics refuse before allocation or operator work and retain zero consumption",
    );
}
