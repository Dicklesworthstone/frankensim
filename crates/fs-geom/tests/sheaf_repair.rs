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
use fs_exec::{BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::router::{ConverterSpec, ErrorModel, MemoryCostOracle, RouteRequest, Router};
use fs_geom::sheaf::{Interface, SheafComplex};
use fs_geom::sheaf_repair::{
    AdmittedSheafSkeleton, COMPONENT_FLOOR, SHEAF_NUMERICS_NORMALIZATION_V1,
    SHEAF_SPECTRUM_NORMALIZATION_V1, SheafNumericsOutcome, SheafNumericsStoppingReason,
    SheafRepairBudget, SheafRepairError, SheafRepairPlanBudget, SheafSkeleton, SheafSkeletonError,
    SheafSpectrumScope, apply_gauge, assess_hodge_decomposition_bounded, hodge_decompose,
    hodge_decompose_bounded, plan_repair, plan_repair_bounded, try_apply_gauge,
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

fn with_budget_cx<R>(budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5348_4541_4652_4550,
                kernel_id: 2,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// A 3-patch triangle complex (one triple junction).
fn triangle() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (1, 2), (0, 2)],
        triangles: vec![(0, 1, 2)],
    }
}

fn canonical_triangle() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (0, 2), (1, 2)],
        triangles: vec![(0, 1, 2)],
    }
}

fn admitted_triangle() -> AdmittedSheafSkeleton {
    AdmittedSheafSkeleton::try_new(3, vec![(0, 1), (0, 2), (1, 2)], vec![(0, 1, 2)])
        .expect("canonical triangle admits")
}

fn bounded_plan_budget(sweeps: usize) -> SheafRepairPlanBudget {
    let max_operator_evaluations = if sweeps == 8 { 56 } else { 2_408 };
    SheafRepairPlanBudget {
        repair: SheafRepairBudget {
            sweeps,
            max_operator_evaluations,
            max_work_items: if sweeps == 8 { 24_576 } else { 300_000 },
            max_scalar_slots: 42,
            poll_stride: 8,
        },
        max_plan_bytes: 8_192,
        max_action_bytes: 4_096,
        max_proposals: 4,
        max_harmonic_support: 3,
    }
}

fn numerics_budget(sweeps: usize) -> SheafRepairBudget {
    SheafRepairBudget {
        sweeps,
        max_operator_evaluations: 256,
        max_work_items: 32_768,
        max_scalar_slots: 128,
        poll_stride: 8,
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
    let gauge_part = sk.d0(&[0.0, 0.7, -0.3]).expect("valid gauge fixture");
    let circ_part = sk.d1t(&[0.4]).expect("valid circulation fixture");
    let m: Vec<f64> = gauge_part
        .iter()
        .zip(&circ_part)
        .map(|(a, b)| a + b)
        .collect();
    let split = hodge_decompose(&sk, &m).expect("valid mixed decomposition fixture");
    // Oracle: dense projection onto im δ⁰ (columns = δ⁰ of unit vertex
    // vectors, vertex 0 pinned) and im δ¹ᵀ.
    let d0_cols: Vec<Vec<f64>> = (1..sk.n_patches)
        .map(|i| {
            let mut e = vec![0.0; sk.n_patches];
            e[i] = 1.0;
            sk.d0(&e).expect("valid dense-reference basis")
        })
        .collect();
    let c_oracle = dense_projection(&m, &d0_cols);
    let exact_oracle = {
        let mut full = vec![0.0; sk.n_patches];
        full[1..].copy_from_slice(&c_oracle);
        sk.d0(&full).expect("valid dense-reference projection")
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
    assert!(norm_inf(&sk.d0t(&split.harmonic).expect("valid harmonic transpose")) < 1e-8);
    assert!(norm_inf(&sk.d1(&split.harmonic).expect("valid harmonic boundary")) < 1e-8);
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
    let mismatch = sk
        .d0(&[0.0, 0.0, 0.012])
        .expect("valid exact-defect fixture");
    let budgets = [0.02, 0.02, 0.02];
    let plan = plan_repair(&sk, &mismatch, &budgets, None).expect("valid repair plan");
    assert!(
        plan.gauge_step_eligible,
        "within budgets: gauge step eligible"
    );
    assert!(plan.split.fractions.0 > 0.999, "pure exact defect");
    assert!(plan.harmonic_support.is_empty(), "no harmonic remainder");
    // Predicted-vs-actual: apply the gauge, re-measure.
    let predicted = plan.proposals[0].expected_post_norm;
    let repaired = try_apply_gauge(&sk, &mismatch, &plan.gauge).expect("valid gauge application");
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
    let plan2 = plan_repair(&sk, &repaired, &budgets, None).expect("valid follow-up plan");
    assert!(
        norm_inf(&plan2.gauge) < 1e-9,
        "no residual gauge on a passing model: {:?}",
        plan2.gauge
    );
    let repaired2 = try_apply_gauge(&sk, &repaired, &plan2.gauge).expect("valid follow-up gauge");
    assert!(
        (norm_inf(&repaired2) - actual).abs() < 1e-12,
        "no-op repair"
    );
    // Over-budget variant: the SAME defect with a tight budget must NOT
    // auto-apply (needs explicit acceptance).
    let tight = [0.001, 0.001, 0.001];
    let gated = plan_repair(&sk, &mismatch, &tight, None).expect("valid gated plan");
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
    let plan = plan_repair(&sk, &mismatch, &budgets, None).expect("valid component plan");
    assert!(
        plan.gauge_step_eligible,
        "a feasible gauge representative exists"
    );
    assert_eq!(plan.gauge, vec![-1.0, 1.0, -2.0, 2.0]);
    assert_eq!(
        sk.d0(&plan.gauge).expect("valid component coboundary"),
        mismatch
    );
    assert_eq!(
        try_apply_gauge(&sk, &mismatch, &plan.gauge).expect("valid component gauge"),
        vec![0.0; 2]
    );

    let impossible =
        plan_repair(&sk, &mismatch, &[0.9, 0.9, 1.9, 1.9], None).expect("valid infeasible plan");
    assert!(
        !impossible.gauge_step_eligible,
        "a component difference larger than the sum of its patch budgets must refuse auto-apply"
    );

    let one_edge = SheafSkeleton {
        n_patches: 2,
        edges: vec![(0, 1)],
        triangles: Vec::new(),
    };
    let slack = plan_repair(&one_edge, &[-2.0], &[10.0, 1.0], None).expect("valid slack plan");
    assert_eq!(
        slack.gauge,
        vec![2.0, 0.0],
        "feasible shift interval [1,3] uses its deterministic maximum-slack midpoint"
    );
    assert_eq!(
        one_edge.d0(&slack.gauge).expect("valid slack coboundary"),
        vec![-2.0]
    );

    let centered =
        plan_repair(&one_edge, &[2.0], &[100.0, 100.0], None).expect("valid centered plan");
    assert_eq!(
        centered.gauge,
        vec![-1.0, 1.0],
        "feasible interval [-100,98] uses its maximum-slack midpoint, not zero"
    );
    assert_eq!(
        one_edge
            .d0(&centered.gauge)
            .expect("valid centered coboundary"),
        vec![2.0]
    );
}

#[test]
fn sr_002aa_skeleton_extraction_refuses_unvalidated_public_complex() {
    let empty = SheafComplex {
        n_patches: 0,
        interfaces: Vec::new(),
        triples: Vec::new(),
        sampling_clip: None,
    };
    assert_eq!(
        SheafSkeleton::of(&empty),
        Err(SheafSkeletonError::EmptyComplex)
    );

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
fn sr_002ab_raw_incidence_and_hodge_refuse_adversarial_inputs() {
    let skeleton = triangle();

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
        skeleton.d0t(&[0.0, 0.0]),
        Err(SheafSkeletonError::CochainLength {
            role: "edge",
            expected: 3,
            actual: 2,
        })
    );
    assert_eq!(
        skeleton.d0t(&[0.0, f64::INFINITY, 0.0]),
        Err(SheafSkeletonError::NonFiniteCochain {
            role: "edge",
            index: 1,
        })
    );
    assert_eq!(
        skeleton.d0t(&[f64::MAX, f64::MAX, f64::MAX]),
        Err(SheafSkeletonError::NumericalOverflow { stage: "d0t" })
    );

    assert_eq!(
        skeleton.d1(&[0.0, 0.0]),
        Err(SheafSkeletonError::CochainLength {
            role: "edge",
            expected: 3,
            actual: 2,
        })
    );
    assert_eq!(
        skeleton.d1(&[0.0, f64::NEG_INFINITY, 0.0]),
        Err(SheafSkeletonError::NonFiniteCochain {
            role: "edge",
            index: 1,
        })
    );
    assert_eq!(
        skeleton.d1(&[f64::MAX, f64::MAX, -f64::MAX]),
        Err(SheafSkeletonError::NumericalOverflow { stage: "d1" })
    );

    assert_eq!(
        skeleton.d1t(&[]),
        Err(SheafSkeletonError::CochainLength {
            role: "triangle",
            expected: 1,
            actual: 0,
        })
    );
    assert_eq!(
        skeleton.d1t(&[f64::NAN]),
        Err(SheafSkeletonError::NonFiniteCochain {
            role: "triangle",
            index: 0,
        })
    );
    let shared_edge_triangles = SheafSkeleton {
        n_patches: 4,
        edges: vec![(0, 1), (1, 2), (0, 2), (1, 3), (0, 3)],
        triangles: vec![(0, 1, 2), (0, 1, 3)],
    };
    assert_eq!(
        shared_edge_triangles.d1t(&[f64::MAX, f64::MAX]),
        Err(SheafSkeletonError::NumericalOverflow { stage: "d1t" })
    );

    let malformed_triangle = SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (1, 2)],
        triangles: vec![(0, 1, 2)],
    };
    assert_eq!(
        malformed_triangle.d1(&[0.0, 0.0]),
        Err(SheafSkeletonError::InvalidTriangle { index: 0 })
    );
    assert_eq!(
        hodge_decompose(&malformed_triangle, &[0.0, 0.0]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::InvalidTriangle { index: 0 }
        ))
    );
    assert_eq!(
        hodge_decompose(&skeleton, &[0.0, 0.0]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "mismatch",
                expected: 3,
                actual: 2,
            }
        ))
    );
    assert_eq!(
        hodge_decompose(&skeleton, &[0.0, f64::NAN, 0.0]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::NonFiniteCochain {
                role: "mismatch",
                index: 1,
            }
        ))
    );
    assert!(matches!(
        hodge_decompose(&skeleton, &[f64::MAX; 3]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::NumericalOverflow { .. }
        )) | Err(SheafRepairError::NumericalOverflow { .. })
    ));

    let vertex = [2.0, -3.0, 5.0];
    let edge = skeleton.d0(&vertex).expect("valid raw coboundary");
    assert_eq!(edge, vec![-5.0, 8.0, 3.0]);
    assert_eq!(
        skeleton.d1(&edge).expect("valid raw boundary composition"),
        vec![0.0],
        "raw incidence preserves caller edge-coordinate order"
    );
    assert_eq!(
        skeleton.d0t(&edge).expect("valid raw transpose"),
        vec![2.0, -13.0, 11.0]
    );
    assert_eq!(
        skeleton.d1t(&[4.0]).expect("valid raw triangle transpose"),
        vec![4.0, 4.0, -4.0]
    );
}

#[test]
fn sr_002b_malformed_repair_inputs_refuse_without_partial_plan() {
    let sk = triangle(); // n_patches = 3
    let mismatch = vec![0.0; sk.edges.len()];
    let baseline = plan_repair(&sk, &mismatch, &[1.0; 3], None).expect("valid baseline plan");

    let malformed = SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 3)],
        triangles: Vec::new(),
    };
    assert_eq!(
        plan_repair(&malformed, &[], &[], None),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::InvalidEdge { index: 0 }
        )),
        "topology refuses before cochain and budget validation"
    );
    assert_eq!(
        try_apply_gauge(&malformed, &[], &[]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::InvalidEdge { index: 0 }
        ))
    );

    assert_eq!(
        plan_repair(&sk, &mismatch, &[1.0, 1.0], None),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "gauge-budget",
                expected: 3,
                actual: 2,
            }
        ))
    );
    assert_eq!(
        plan_repair(&sk, &mismatch[..2], &[1.0; 3], None),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "mismatch",
                expected: 3,
                actual: 2,
            }
        ))
    );
    assert_eq!(
        plan_repair(&sk, &[0.0; 4], &[1.0; 3], None),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "mismatch",
                expected: 3,
                actual: 4,
            }
        ))
    );
    assert_eq!(
        plan_repair(&sk, &mismatch, &[1.0; 4], None),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "gauge-budget",
                expected: 3,
                actual: 4,
            }
        ))
    );
    assert_eq!(
        plan_repair(&sk, &[0.0, f64::NAN, 0.0], &[1.0; 3], None),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::NonFiniteCochain {
                role: "mismatch",
                index: 1,
            }
        ))
    );
    assert_eq!(
        plan_repair(&sk, &mismatch, &[1.0, f64::INFINITY, 1.0], None),
        Err(SheafRepairError::InvalidGaugeBudget { index: 1 })
    );
    assert_eq!(
        plan_repair(&sk, &mismatch, &[-1.0, 1.0, 1.0], None),
        Err(SheafRepairError::InvalidGaugeBudget { index: 0 })
    );

    assert_eq!(
        try_apply_gauge(&sk, &mismatch[..2], &[0.0; 3]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "mismatch",
                expected: 3,
                actual: 2,
            }
        ))
    );
    assert_eq!(
        try_apply_gauge(&sk, &mismatch, &[0.0; 2]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "gauge",
                expected: 3,
                actual: 2,
            }
        ))
    );
    assert_eq!(
        try_apply_gauge(&sk, &[0.0; 4], &[0.0; 3]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "mismatch",
                expected: 3,
                actual: 4,
            }
        ))
    );
    assert_eq!(
        try_apply_gauge(&sk, &mismatch, &[0.0; 4]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "gauge",
                expected: 3,
                actual: 4,
            }
        ))
    );
    assert_eq!(
        try_apply_gauge(&sk, &mismatch, &[0.0, f64::NAN, 0.0]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::NonFiniteCochain {
                role: "gauge",
                index: 1,
            }
        ))
    );
    assert_eq!(
        try_apply_gauge(&sk, &mismatch, &[f64::MAX, -f64::MAX, 0.0]),
        Err(SheafRepairError::NumericalOverflow {
            stage: "apply-gauge",
        })
    );
    assert_eq!(
        apply_gauge(&sk, &mismatch[..2], &[0.0; 3]),
        Err(SheafRepairError::Skeleton(
            SheafSkeletonError::CochainLength {
                role: "mismatch",
                expected: 3,
                actual: 2,
            }
        )),
        "the compatibility name preserves the typed refusal"
    );

    let replay = plan_repair(&sk, &mismatch, &[1.0; 3], None).expect("valid retry plan");
    assert_eq!(
        replay, baseline,
        "refusal cannot publish or retain a partial plan"
    );
}

#[test]
fn sr_003_coexact_seeding_retains_noncausal_diagnostic() {
    let sk = triangle();
    // Seed a pure circulation (the flipped-orientation signature): the
    // image of δ¹ᵀ.
    let mismatch = sk.d1t(&[0.05]).expect("valid coexact-seeding fixture");
    let plan = plan_repair(&sk, &mismatch, &[1.0; 3], None).expect("valid coexact plan");
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
    let repaired = try_apply_gauge(&sk, &mismatch, &plan.gauge).expect("valid coexact gauge");
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
    let plan = plan_repair(&sk, &mismatch, &[1.0; 4], None).expect("valid harmonic plan");
    assert!(
        norm_inf(&plan.split.harmonic) > 0.9 * eps,
        "the retained harmonic witness must be nonzero: {:?}",
        plan.split.harmonic
    );
    assert!(
        norm_inf(
            &sk.d1(&plan.split.harmonic)
                .expect("valid retained harmonic boundary")
        ) < 1e-12,
        "the retained mismatch cochain must be closed"
    );
    assert!(
        norm_inf(
            &sk.d0t(&plan.split.harmonic)
                .expect("valid retained harmonic transpose")
        ) < 1e-12,
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
    let repaired = try_apply_gauge(&sk, &mismatch, &plan.gauge).expect("valid harmonic gauge");
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
    let exact = sk
        .d0(&[0.0, 1.0, 0.0, -1.0])
        .expect("valid subfloor exact fixture");
    let eps = 1e-8;
    let cycle = [eps, eps, eps, -eps];
    let mismatch: Vec<f64> = exact.iter().zip(cycle).map(|(a, b)| a + b).collect();
    let plan = plan_repair(&sk, &mismatch, &[2.0; 4], None).expect("valid subfloor plan");
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
    let mismatch = sk
        .d0(&[0.0, 0.0, 0.012])
        .expect("valid reroute mismatch fixture");
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
    let plan = plan_repair(&sk, &mismatch, &[1.0; 3], Some((&router, &oracle, &req)))
        .expect("valid reroute plan");
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
        max_work_items: 8_192,
        max_scalar_slots: 42,
        poll_stride: 1,
    };
    with_cx(|cx| {
        let bounded = hodge_decompose_bounded(&skeleton, &mismatch, budget, cx)
            .expect("exact admitted schedule fits its envelope");
        assert_eq!(bounded.budget, budget);
        assert_eq!(bounded.usage.completed_sweeps, 16);
        assert_eq!(bounded.usage.operator_evaluations, 56);
        assert_eq!(bounded.usage.admitted_work_items, 5_760);
        assert!(bounded.usage.work_items <= bounded.usage.admitted_work_items);
        assert_eq!(bounded.usage.admitted_scalar_slots, 42);
        assert_eq!(bounded.usage.ambient_budget.refusal, None);
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
                    max_work_items: 5_759,
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::WorkItemBudgetExceeded {
                stage: "hodge-work-preflight",
                required: 5_760,
                cap: 5_759,
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
        "bounded decomposition retains exact operator work plus its conservative memory envelope and refuses under-budget or non-finite requests before publication",
    );
}

#[test]
fn sr_009_bounded_decomposition_observes_pre_cancellation() {
    let skeleton =
        AdmittedSheafSkeleton::try_new(2, vec![(0, 1)], Vec::new()).expect("single edge admits");
    let budget = SheafRepairBudget {
        sweeps: 2,
        max_operator_evaluations: 32,
        max_work_items: 1_024,
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
                work_items: 0,
            })
        );
    });
    verdict(
        "sr-009",
        "pre-cancelled diagnostics refuse before allocation or operator work and retain zero consumption",
    );
}

#[test]
fn sr_010_bounded_plan_retains_complete_budget_and_usage() {
    let skeleton = admitted_triangle();
    let mismatch = skeleton
        .d0(&[0.0, 0.5, -0.25])
        .expect("finite exact fixture");
    let budget = bounded_plan_budget(8);
    with_cx(|cx| {
        let bounded = plan_repair_bounded(&skeleton, &mismatch, &[1.0; 3], None, budget, cx)
            .expect("complete admitted plan");
        assert_eq!(bounded.budget, budget);
        assert_eq!(bounded.usage.repair.completed_sweeps, 16);
        assert_eq!(bounded.usage.repair.operator_evaluations, 56);
        assert_eq!(bounded.usage.repair.admitted_scalar_slots, 42);
        assert!(bounded.usage.repair.work_items > 0);
        assert!(
            bounded.usage.repair.work_items <= bounded.usage.repair.admitted_work_items,
            "measured scalar/graph/string work must stay inside preflight"
        );
        assert_eq!(
            bounded.usage.repair.ambient_budget.cost_charged,
            (bounded.usage.repair.work_items + bounded.usage.repair.operator_evaluations) as u64
        );
        assert_eq!(bounded.usage.repair.ambient_budget.refusal, None);
        assert!(bounded.usage.plan_memory_envelope <= budget.max_plan_bytes);
        assert!(bounded.usage.reserved_plan_bytes <= bounded.usage.plan_memory_envelope);
        assert!(bounded.usage.action_bytes <= budget.max_action_bytes);
        assert_eq!(bounded.usage.proposals, bounded.plan.proposals.len());
        assert_eq!(
            bounded.usage.harmonic_support,
            bounded.plan.harmonic_support.len()
        );
        assert!(bounded.plan.gauge_step_eligible);

        let required_work = bounded.usage.repair.admitted_work_items;
        assert_eq!(
            plan_repair_bounded(
                &skeleton,
                &mismatch,
                &[1.0; 3],
                None,
                SheafRepairPlanBudget {
                    repair: SheafRepairBudget {
                        max_work_items: required_work - 1,
                        ..budget.repair
                    },
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::WorkItemBudgetExceeded {
                stage: "plan-work-preflight",
                required: required_work as u128,
                cap: required_work - 1,
            })
        );

        let required_memory = bounded.usage.plan_memory_envelope;
        assert_eq!(
            plan_repair_bounded(
                &skeleton,
                &mismatch,
                &[1.0; 3],
                None,
                SheafRepairPlanBudget {
                    max_plan_bytes: required_memory - 1,
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::PlanMemoryBudgetExceeded {
                required: required_memory as u128,
                cap: required_memory - 1,
            })
        );

        assert!(matches!(
            plan_repair_bounded(
                &skeleton,
                &mismatch,
                &[1.0; 3],
                None,
                SheafRepairPlanBudget {
                    max_proposals: 0,
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::OutputBudgetExceeded {
                resource: "proposals",
                required: 1,
                cap: 0,
            })
        ));
        assert!(matches!(
            plan_repair_bounded(
                &skeleton,
                &mismatch,
                &[1.0; 3],
                None,
                SheafRepairPlanBudget {
                    max_action_bytes: 1,
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::OutputBudgetExceeded {
                resource: "action-bytes",
                cap: 1,
                ..
            })
        ));
        assert_eq!(
            plan_repair_bounded(
                &skeleton,
                &[0.0, f64::NAN, 0.0],
                &[1.0, -1.0, 1.0],
                None,
                budget,
                cx,
            ),
            Err(SheafRepairError::Skeleton(
                SheafSkeletonError::NonFiniteCochain {
                    role: "mismatch",
                    index: 1,
                }
            ))
        );
        assert_eq!(
            plan_repair_bounded(&skeleton, &mismatch, &[1.0, -1.0, 1.0], None, budget, cx,),
            Err(SheafRepairError::InvalidGaugeBudget { index: 1 })
        );
    });
    verdict(
        "sr-010",
        "bounded planning retains enforced work, memory, output, and ambient consumption and refuses every undersized envelope before publication",
    );
}

#[test]
fn sr_011_bounded_plan_matches_legacy_arithmetic_and_replays() {
    let raw = canonical_triangle();
    let admitted = admitted_triangle();
    let mismatch = raw
        .d0(&[0.0, -0.375, 0.625])
        .expect("valid legacy-parity fixture");
    let gauge_budgets = [1.0; 3];
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
        .expect("register bounded reroute fixture");
    let oracle = MemoryCostOracle::new();
    let request = RouteRequest {
        from: "sdf".to_string(),
        to: "mesh".to_string(),
        scale: 1.0,
        max_abs_error: 1e-3,
        max_cost_s: 100.0,
    };
    let route = router
        .plan(&request, &oracle)
        .expect("route is admitted before bounded repair planning");
    let legacy = plan_repair(
        &raw,
        &mismatch,
        &gauge_budgets,
        Some((&router, &oracle, &request)),
    )
    .expect("legacy plan over canonical admitted meaning");
    let budget = bounded_plan_budget(400);
    with_cx(|cx| {
        let bounded = plan_repair_bounded(
            &admitted,
            &mismatch,
            &gauge_budgets,
            Some(&route),
            budget,
            cx,
        )
        .expect("bounded compatibility plan");
        assert_eq!(bounded.plan, legacy);
        assert_eq!(bounded.usage.repair.completed_sweeps, 800);
        assert_eq!(bounded.usage.repair.operator_evaluations, 2_408);

        let replay = plan_repair_bounded(
            &admitted,
            &mismatch,
            &gauge_budgets,
            Some(&route),
            budget,
            cx,
        )
        .expect("deterministic replay");
        assert_eq!(replay, bounded);
    });
    verdict(
        "sr-011",
        "the 400-sweep admitted planner preserves legacy arithmetic, proposal text, ranking, and deterministic replay while adding bounded authority",
    );
}

#[test]
fn sr_012_bounded_plan_refuses_pre_cancel_and_final_publication() {
    let skeleton = admitted_triangle();
    let mismatch = skeleton
        .d0(&[0.0, 0.5, -0.25])
        .expect("finite exact fixture");
    let plan_budget = bounded_plan_budget(8);

    let gate = CancelGate::new();
    gate.request();
    with_gate_cx(&gate, |cx| {
        assert_eq!(
            plan_repair_bounded(&skeleton, &mismatch, &[1.0; 3], None, plan_budget, cx),
            Err(SheafRepairError::Cancelled {
                stage: "plan-admission",
                completed_sweeps: 0,
                operator_evaluations: 0,
                work_items: 0,
            })
        );
    });

    let generous = Budget {
        deadline: None,
        poll_quota: 100_000,
        cost_quota: None,
        priority: 0,
    };
    let baseline = with_budget_cx(generous, |cx| {
        plan_repair_bounded(&skeleton, &mismatch, &[1.0; 3], None, plan_budget, cx)
            .expect("generous poll budget")
    });
    let planned_cost = baseline.usage.repair.ambient_budget.planned_cost;
    let cost_refusal_budget = Budget {
        cost_quota: Some(planned_cost - 1),
        ..generous
    };
    assert!(matches!(
        with_budget_cx(cost_refusal_budget, |cx| {
            plan_repair_bounded(&skeleton, &mismatch, &[1.0; 3], None, plan_budget, cx)
        }),
        Err(SheafRepairError::AmbientBudgetRefused {
            refusal: BudgetRefusal::CostPlanExceedsQuota { planned, quota },
            completed_sweeps: 0,
            operator_evaluations: 0,
            work_items: 0,
        }) if planned == planned_cost && quota == planned_cost - 1
    ));
    let final_poll = baseline.usage.repair.ambient_budget.polls_used;
    assert!(final_poll > 1, "fixture must cross multiple boundaries");
    let final_refusal_budget = Budget {
        poll_quota: final_poll - 1,
        ..generous
    };
    let refusal = with_budget_cx(final_refusal_budget, |cx| {
        plan_repair_bounded(&skeleton, &mismatch, &[1.0; 3], None, plan_budget, cx)
    });
    assert!(matches!(
        refusal,
        Err(SheafRepairError::AmbientBudgetRefused {
            refusal: BudgetRefusal::PollsExhausted {
                phase: "plan-publication",
                quota,
            },
            completed_sweeps: 16,
            operator_evaluations: 56,
            ..
        }) if quota == final_poll - 1
    ));

    let retry = with_budget_cx(generous, |cx| {
        plan_repair_bounded(&skeleton, &mismatch, &[1.0; 3], None, plan_budget, cx)
            .expect("healthy retry after refused publication")
    });
    assert_eq!(retry, baseline);
    verdict(
        "sr-012",
        "pre-cancellation and final-publication exhaustion return no plan, while a healthy retry deterministically reproduces the complete baseline",
    );
}

#[test]
fn sr_013_bounded_plan_accounts_mixed_components_and_route() {
    let skeleton = AdmittedSheafSkeleton::try_new(
        7,
        vec![(0, 1), (0, 2), (1, 2), (3, 4), (3, 6), (4, 5), (5, 6)],
        vec![(0, 1, 2)],
    )
    .expect("filled triangle plus disconnected ring admits");
    // Triangle entries are exactly d1-transpose(1); ring entries are one
    // consistently oriented cycle. The fixture therefore exercises retained
    // coexact and harmonic components without relying on convergence noise.
    let mismatch = [1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0];

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
        .expect("register mixed-plan reroute fixture");
    let oracle = MemoryCostOracle::new();
    let request = RouteRequest {
        from: "sdf".to_string(),
        to: "mesh".to_string(),
        scale: 1.0,
        max_abs_error: 1e-3,
        max_cost_s: 100.0,
    };
    let route = router
        .plan(&request, &oracle)
        .expect("route admits before bounded planning");
    let budget = SheafRepairPlanBudget {
        repair: SheafRepairBudget {
            sweeps: 8,
            max_operator_evaluations: 124,
            max_work_items: 65_536,
            max_scalar_slots: 90,
            poll_stride: 8,
        },
        max_plan_bytes: 8_192,
        max_action_bytes: 4_096,
        max_proposals: 4,
        max_harmonic_support: 7,
    };

    with_cx(|cx| {
        let bounded =
            plan_repair_bounded(&skeleton, &mismatch, &[2.0; 7], Some(&route), budget, cx)
                .expect("mixed bounded plan");
        assert!(bounded.plan.split.fractions.0 <= COMPONENT_FLOOR);
        assert!(bounded.plan.split.fractions.1 > COMPONENT_FLOOR);
        assert!(bounded.plan.split.fractions.2 > COMPONENT_FLOOR);
        assert_eq!(
            bounded.plan.harmonic_support,
            vec![((3, 4), 1.0), ((3, 6), 1.0), ((4, 5), 1.0), ((5, 6), 1.0),]
        );
        assert_eq!(bounded.usage.repair.completed_sweeps, 16);
        assert_eq!(bounded.usage.repair.operator_evaluations, 124);
        assert_eq!(bounded.usage.repair.admitted_scalar_slots, 90);
        assert!(bounded.usage.repair.work_items <= bounded.usage.repair.admitted_work_items);
        assert_eq!(bounded.usage.proposals, 3);
        assert_eq!(bounded.usage.harmonic_support, 4);
        assert!(
            bounded.plan.proposals[0]
                .action
                .contains("coexact circulation")
        );
        assert!(bounded.plan.proposals[1].action.contains("reroute"));
        assert!(
            bounded.plan.proposals[2]
                .action
                .contains("retained harmonic remainder")
        );

        assert_eq!(
            plan_repair_bounded(
                &skeleton,
                &mismatch,
                &[2.0; 7],
                Some(&route),
                SheafRepairPlanBudget {
                    max_harmonic_support: 3,
                    ..budget
                },
                cx,
            ),
            Err(SheafRepairError::OutputBudgetExceeded {
                resource: "harmonic-support",
                required: 4,
                cap: 3,
            })
        );
    });
    verdict(
        "sr-013",
        "comparison-heavy incidence lookup, coexact localization, harmonic support, route formatting, output caps, and stable proposal ranking share one bounded accountant",
    );
}

#[test]
fn sr_014_numerics_converges_with_one_gauge_root_per_component() {
    let skeleton = AdmittedSheafSkeleton::try_new(5, vec![(0, 1), (2, 3)], Vec::new())
        .expect("two edges plus one isolated patch admit");
    let budget = numerics_budget(4);
    let outcome =
        with_cx(|cx| assess_hodge_decomposition_bounded(&skeleton, &[2.0, 4.0], 1e-12, budget, cx));
    let SheafNumericsOutcome::Converged(converged) = outcome else {
        panic!("independent exact components must converge: {outcome:?}");
    };
    assert_eq!(converged.exact(), &[2.0, 4.0]);
    assert_eq!(converged.potential(), &[0.0, 2.0, 0.0, 4.0, 0.0]);
    assert!(converged.coexact().iter().all(|value| *value == 0.0));
    assert!(converged.harmonic().iter().all(|value| *value == 0.0));
    let receipt = converged.receipt();
    assert_eq!(receipt.normalization_id, SHEAF_NUMERICS_NORMALIZATION_V1);
    assert_eq!(
        receipt.stopping_reason,
        SheafNumericsStoppingReason::ResidualBoundsSatisfied
    );
    assert!(receipt.primal_normal_equation.normalized.hi() <= 1e-12);
    assert!(receipt.dual_normal_equation.normalized.hi() <= 1e-12);
    assert!(receipt.reconstruction.normalized.hi() <= 1e-12);
    assert_eq!(receipt.source.n_patches(), 5);
    assert_eq!(receipt.source.edges(), &[(0, 1), (2, 3)]);
    assert_eq!(receipt.source.triangles(), &[]);
    assert_eq!(receipt.source.mismatch(), &[2.0, 4.0]);
    match &receipt.spectrum {
        SheafSpectrumScope::Unknown(report) => {
            assert_eq!(report.normalization_id, SHEAF_SPECTRUM_NORMALIZATION_V1);
            assert_eq!(report.normalization_scale, 2.0);
            assert_eq!((report.nullspace.lower, report.nullspace.upper), (3, 3));
            assert_eq!(report.nullspace.component_roots, &[0, 2, 4]);
            assert_eq!(report.structural_zero_cluster.normalized_hull.lo(), 0.0);
            assert_eq!(report.structural_zero_cluster.normalized_hull.hi(), 0.0);
            assert_eq!(report.structural_zero_cluster.multiplicity_lower, 3);
            assert_eq!(report.structural_zero_cluster.multiplicity_upper, 3);
            assert!(report.candidate_clusters.is_empty());
            assert_eq!(report.requested_range, None);
            assert_eq!(report.covered_range, None);
            assert_eq!(report.unresolved_modes, 5);
        }
    }

    let near_max = f64::MAX / 8.0;
    let one_edge =
        AdmittedSheafSkeleton::try_new(2, vec![(0, 1)], Vec::new()).expect("edge admits");
    let scaled =
        with_cx(|cx| assess_hodge_decomposition_bounded(&one_edge, &[near_max], 1e-12, budget, cx));
    let SheafNumericsOutcome::Converged(scaled) = scaled else {
        panic!("scale-safe exact edge must converge near f64::MAX: {scaled:?}");
    };
    assert_eq!(scaled.exact(), &[near_max]);
    assert!(
        scaled
            .receipt()
            .primal_normal_equation
            .absolute
            .hi()
            .is_finite()
    );
    assert!(scaled.receipt().primal_normal_equation.normalized.hi() <= 1e-12);

    let relabeled_source =
        with_cx(|cx| assess_hodge_decomposition_bounded(&skeleton, &[2.0, 5.0], 1e-12, budget, cx));
    let SheafNumericsOutcome::Converged(relabeled_source) = relabeled_source else {
        panic!("independent source-binding fixture must converge: {relabeled_source:?}");
    };
    assert_ne!(
        &receipt.source,
        &relabeled_source.receipt().source,
        "same-shape mismatch values cannot share a convergence source binding"
    );
    verdict(
        "sr-014",
        "the converged token retains isolated-component gauge roots, explicit spectral Unknown, replay residuals, and scale-safe near-MAX arithmetic",
    );
}

#[test]
fn sr_015_numerics_keeps_false_convergence_indeterminate_or_refused() {
    let path =
        AdmittedSheafSkeleton::try_new(3, vec![(0, 1), (1, 2)], Vec::new()).expect("path admits");
    let budget = numerics_budget(1);
    let outcome =
        with_cx(|cx| assess_hodge_decomposition_bounded(&path, &[1.0, 1.0], 1e-12, budget, cx));
    let SheafNumericsOutcome::Indeterminate(report) = outcome else {
        panic!("one sweep must not launder a partial fit: {outcome:?}");
    };
    assert_eq!(
        report.receipt().stopping_reason,
        SheafNumericsStoppingReason::SweepLimitReached
    );
    assert!(
        report.receipt().primal_normal_equation.normalized.hi() > 1e-12,
        "the outward residual upper bound exposes the unfinished normal equation"
    );
    assert_eq!(report.coboundary_candidate(), &[0.0, 1.0]);
    assert_eq!(report.remainder_candidate(), &[1.0, 0.0]);

    let invalid =
        with_cx(|cx| assess_hodge_decomposition_bounded(&path, &[1.0, 1.0], f64::NAN, budget, cx));
    assert_eq!(
        invalid,
        SheafNumericsOutcome::Refused(SheafRepairError::InvalidTolerance {
            field: "relative_tolerance",
        })
    );
    let underfunded = with_cx(|cx| {
        assess_hodge_decomposition_bounded(
            &path,
            &[1.0, 1.0],
            1e-12,
            SheafRepairBudget {
                max_operator_evaluations: 1,
                ..budget
            },
            cx,
        )
    });
    assert!(matches!(
        underfunded,
        SheafNumericsOutcome::Refused(SheafRepairError::WorkBudgetExceeded { .. })
    ));

    let gate = CancelGate::new();
    gate.request();
    let cancelled = with_gate_cx(&gate, |cx| {
        assess_hodge_decomposition_bounded(&path, &[1.0, 1.0], 1e-12, budget, cx)
    });
    assert!(matches!(
        cancelled,
        SheafNumericsOutcome::Refused(SheafRepairError::Cancelled {
            stage: "numerics-admission",
            completed_sweeps: 0,
            operator_evaluations: 0,
            work_items: 0,
        })
    ));
    verdict(
        "sr-015",
        "a fixed-sweep false convergence publishes only candidate-named Indeterminate evidence, while invalid tolerance, preflight work caps, and pre-cancellation refuse",
    );
}

#[test]
fn sr_016_numerics_tolerance_relabeling_and_scale_metamorphics() {
    let path =
        AdmittedSheafSkeleton::try_new(3, vec![(0, 1), (1, 2)], Vec::new()).expect("path admits");
    let one_sweep = numerics_budget(1);
    let baseline =
        with_cx(|cx| assess_hodge_decomposition_bounded(&path, &[1.0, 1.0], 0.0, one_sweep, cx));
    let SheafNumericsOutcome::Indeterminate(baseline) = baseline else {
        panic!("zero tolerance must expose the unfinished path fit: {baseline:?}");
    };
    let receipt = baseline.receipt();
    let boundary = [
        receipt.primal_normal_equation.normalized,
        receipt.dual_normal_equation.normalized,
        receipt.remainder_exact_orthogonality.normalized,
        receipt.coboundary_triangle_orthogonality.normalized,
        receipt.coboundary_remainder_orthogonality.normalized,
        receipt.triangle_remainder_orthogonality.normalized,
        receipt.reconstruction.normalized,
    ]
    .into_iter()
    .map(|bounds| bounds.hi())
    .fold(0.0, f64::max);
    assert!(boundary.is_finite() && boundary > 0.0);
    assert!(matches!(
        with_cx(|cx| assess_hodge_decomposition_bounded(
            &path,
            &[1.0, 1.0],
            boundary,
            one_sweep,
            cx,
        )),
        SheafNumericsOutcome::Converged(_)
    ));
    assert!(matches!(
        with_cx(|cx| assess_hodge_decomposition_bounded(
            &path,
            &[1.0, 1.0],
            boundary.next_down(),
            one_sweep,
            cx,
        )),
        SheafNumericsOutcome::Indeterminate(_)
    ));

    let unit_scale = 2.0f64.powi(400);
    let scaled = with_cx(|cx| {
        assess_hodge_decomposition_bounded(&path, &[unit_scale, unit_scale], 0.0, one_sweep, cx)
    });
    let SheafNumericsOutcome::Indeterminate(scaled) = scaled else {
        panic!("unit scaling cannot promote an unfinished fit: {scaled:?}");
    };
    let baseline_primal = baseline.receipt().primal_normal_equation.normalized;
    let scaled_primal = scaled.receipt().primal_normal_equation.normalized;
    assert!(
        baseline_primal.lo() <= scaled_primal.hi() && scaled_primal.lo() <= baseline_primal.hi(),
        "power-of-two unit scaling must preserve overlapping normalized residual evidence"
    );

    // Relabel old vertices [0, 1, 2, 3, 4] as [4, 2, 3, 0, 1]. Both
    // retained edge orientations reverse and canonical edge order swaps, so
    // the transported cochain is [-4, -2].
    let relabeled = AdmittedSheafSkeleton::try_new(5, vec![(0, 3), (2, 4)], Vec::new())
        .expect("relabel admits");
    let relabeled_outcome = with_cx(|cx| {
        assess_hodge_decomposition_bounded(&relabeled, &[-4.0, -2.0], 1e-12, numerics_budget(4), cx)
    });
    let SheafNumericsOutcome::Converged(relabeled_outcome) = relabeled_outcome else {
        panic!("transported exact components must remain converged: {relabeled_outcome:?}");
    };
    assert_eq!(relabeled_outcome.exact(), &[-4.0, -2.0]);
    assert_eq!(relabeled_outcome.potential(), &[0.0, 0.0, 0.0, -4.0, -2.0]);
    match &relabeled_outcome.receipt().spectrum {
        SheafSpectrumScope::Unknown(report) => {
            assert_eq!(report.nullspace.component_roots, &[0, 1, 2]);
        }
    }
    verdict(
        "sr-016",
        "inclusive tolerance gating, power-of-two unit scaling, and orientation-changing vertex relabeling preserve the bounded numerical authority boundary",
    );
}

#[test]
fn sr_017_edgeless_spectrum_has_one_set_valued_zero_cluster() {
    let edgeless =
        AdmittedSheafSkeleton::try_new(4, Vec::new(), Vec::new()).expect("isolates admit");
    let outcome = with_cx(|cx| {
        assess_hodge_decomposition_bounded(&edgeless, &[], 0.0, numerics_budget(1), cx)
    });
    let SheafNumericsOutcome::Converged(outcome) = outcome else {
        panic!("the empty cochain on isolated patches must converge: {outcome:?}");
    };
    let SheafSpectrumScope::Unknown(report) = &outcome.receipt().spectrum;
    assert_eq!(report.normalization_scale, 1.0);
    assert_eq!((report.nullspace.lower, report.nullspace.upper), (4, 4));
    assert_eq!(report.nullspace.component_roots, &[0, 1, 2, 3]);
    assert_eq!(report.structural_zero_cluster.normalized_hull.lo(), 0.0);
    assert_eq!(report.structural_zero_cluster.normalized_hull.hi(), 0.0);
    assert_eq!(report.structural_zero_cluster.multiplicity_lower, 4);
    assert_eq!(report.structural_zero_cluster.multiplicity_upper, 4);
    assert!(report.candidate_clusters.is_empty());
    assert_eq!(report.requested_range, None);
    assert_eq!(report.covered_range, None);
    assert_eq!(report.unresolved_modes, 4);

    let star = AdmittedSheafSkeleton::try_new(4, vec![(0, 1), (0, 2), (0, 3)], Vec::new())
        .expect("star admits");
    let star = with_cx(|cx| {
        assess_hodge_decomposition_bounded(&star, &[0.0; 3], 0.0, numerics_budget(1), cx)
    });
    let SheafNumericsOutcome::Converged(star) = star else {
        panic!("the zero cochain on a star must converge: {star:?}");
    };
    let SheafSpectrumScope::Unknown(star_report) = &star.receipt().spectrum;
    assert_eq!(star_report.normalization_scale, 6.0);
    assert_eq!(
        (star_report.nullspace.lower, star_report.nullspace.upper),
        (1, 1)
    );
    assert_eq!(star_report.nullspace.component_roots, &[0]);
    verdict(
        "sr-017",
        "edgeless and degree-three admitted graphs report exact structural nullity under the versioned Gershgorin scale while all numerical modes remain unresolved",
    );
}
