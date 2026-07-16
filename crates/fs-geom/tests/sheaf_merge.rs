//! Sheaf-merge conformance (the lmp4.12 crown jewel; runs under the
//! `sheaf-merge` feature). Acceptance: gauge-fit mismatches auto-reconcile
//! only after a nominal post-state residual check; dominant fixed-iteration
//! remainders remain candidate conflicts localized to the right cells with
//! both caller-supplied parent labels; type-level collisions are caught before any
//! decomposition; degraded-gap merges are flagged low-confidence; the
//! Sev-0 adversarial case (reconciliation cannot reach a passing state)
//! ESCALATES rather than reporting a resolved state falsely; trivial merges
//! take the fast paths; the candidate diagnostic measures retained conflicts.
#![cfg(feature = "sheaf-merge")]

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering};

use asupersync::time::TimeSource;
use asupersync::types::{Budget, Time};
use fs_exec::{BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::sheaf_merge::{
    BoundedCandidateRate, BoundedMergeOutcome, BranchState, CandidateRateBudget,
    CandidateRateError, Confidence, MergeOutcome, SheafMergeBudget, SheafMergeError,
    candidate_remainder_conflict_rate, candidate_remainder_conflict_rate_bounded, spectral_gap,
    three_way_merge, three_way_merge_bounded, type_conflicts,
};
use fs_geom::sheaf_repair::{
    AdmittedSheafSkeleton, SheafRepairBudget, SheafRepairError, SheafSkeleton, apply_gauge,
    hodge_decompose,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-geom/sheaf-merge\",\"case\":\"{case}\",\"verdict\":\"pass\",\
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
                seed: 0x5348_4541_464d_4552,
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
                seed: 0x5348_4541_464d_4552,
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

struct RequestingClock<'a> {
    gate: &'a CancelGate,
    reads: AtomicUsize,
    request_on: usize,
}

impl<'a> RequestingClock<'a> {
    fn new(gate: &'a CancelGate, request_on: usize) -> Self {
        Self {
            gate,
            reads: AtomicUsize::new(0),
            request_on,
        }
    }

    fn reads(&self) -> usize {
        self.reads.load(Ordering::SeqCst)
    }
}

impl TimeSource for RequestingClock<'_> {
    fn now(&self) -> Time {
        if self.reads.fetch_add(1, Ordering::SeqCst) + 1 == self.request_on {
            self.gate.request();
        }
        Time::ZERO
    }
}

fn with_clock_cx<R>(
    gate: &CancelGate,
    budget: Budget,
    clock: &dyn TimeSource,
    f: impl FnOnce(&Cx<'_>) -> R,
) -> R {
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x5348_4541_464d_4552,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            budget,
            ExecMode::Deterministic,
        )
        .with_time_source(clock);
        f(&cx)
    })
}

fn bounded_budget(sweeps: usize) -> SheafMergeBudget {
    SheafMergeBudget {
        repair: SheafRepairBudget {
            sweeps,
            max_operator_evaluations: 10_000,
            max_work_items: 10_000_000,
            max_scalar_slots: 10_000,
            poll_stride: 8,
        },
        spectral_sweeps: 64,
        max_scalar_slots: 10_000,
        max_output_bytes: 16_384,
        max_conflict_cells: 64,
        max_provenance_bytes: 512,
    }
}

fn candidate_budget(trials: usize, sweeps: usize) -> CandidateRateBudget {
    let merge = bounded_budget(sweeps);
    CandidateRateBudget {
        merge,
        max_trials: trials,
        max_operator_evaluations: 1_000_000,
        max_work_items: 1_000_000_000,
        max_scalar_slots: 100_000,
        max_total_output_bytes: merge
            .max_output_bytes
            .checked_mul(trials)
            .expect("small fixture output envelope"),
    }
}

fn canonical_triangle() -> (SheafSkeleton, AdmittedSheafSkeleton) {
    let raw = SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (0, 2), (1, 2)],
        triangles: vec![(0, 1, 2)],
    };
    let admitted = AdmittedSheafSkeleton::admit(raw.clone()).expect("canonical triangle admits");
    (raw, admitted)
}

fn canonical_ring() -> (SheafSkeleton, AdmittedSheafSkeleton) {
    let raw = SheafSkeleton {
        n_patches: 4,
        edges: vec![(0, 1), (0, 3), (1, 2), (2, 3)],
        triangles: Vec::new(),
    };
    let admitted = AdmittedSheafSkeleton::admit(raw.clone()).expect("canonical ring admits");
    (raw, admitted)
}

fn canonical_weighted_star() -> (SheafSkeleton, AdmittedSheafSkeleton) {
    let raw = SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (0, 2)],
        triangles: Vec::new(),
    };
    let admitted =
        AdmittedSheafSkeleton::admit(raw.clone()).expect("canonical weighted star admits");
    (raw, admitted)
}

/// 3-patch triangle (contractible: no harmonic space).
fn triangle() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (1, 2), (0, 2)],
        triangles: vec![(0, 1, 2)],
    }
}

/// 4-patch ring (one cycle: harmonic space is 1-dimensional).
fn ring() -> SheafSkeleton {
    SheafSkeleton {
        n_patches: 4,
        edges: vec![(0, 1), (1, 2), (2, 3), (0, 3)],
        triangles: vec![],
    }
}

fn branch(name: &str, mismatch: Vec<f64>) -> BranchState {
    BranchState {
        provenance: name.to_string(),
        mismatch,
        assignments: BTreeMap::new(),
    }
}

fn norm_inf(values: &[f64]) -> f64 {
    assert!(
        values.iter().all(|value| value.is_finite()),
        "test norm requires finite values: {values:?}"
    );
    values
        .iter()
        .fold(0.0f64, |largest, value| largest.max(value.abs()))
}

#[test]
fn sm_001_gauge_fit_auto_reconciles_with_checked_residual() {
    let sk = triangle();
    let base = vec![0.0; 3];
    // X re-gauges patch 1, Y re-gauges patch 2: pure coboundary edits.
    let x = branch(
        "agent-x@c1",
        sk.d0(&[0.0, 0.02, 0.0]).expect("valid agent-x coboundary"),
    );
    let y = branch(
        "agent-y@c2",
        sk.d0(&[0.0, 0.0, -0.015])
            .expect("valid agent-y coboundary"),
    );
    let out = three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6);
    match out {
        MergeOutcome::Resolved {
            merged,
            gauge,
            residual_receipt,
            confidence,
        } => {
            assert!(
                residual_receipt.post_norm <= residual_receipt.tol,
                "nominal residual was checked"
            );
            assert!(
                merged.iter().all(|v| v.abs() < 1e-9),
                "nominal residual-passing state: {merged:?}"
            );
            // This fixture's deterministic gauge fit recovered both offsets.
            assert!((gauge[1] - 0.02).abs() < 1e-9, "{gauge:?}");
            assert!((gauge[2] + 0.015).abs() < 1e-9, "{gauge:?}");
            assert!(matches!(confidence, Confidence::Normal { .. }));
        }
        other => panic!("coboundary edits must auto-resolve: {other:?}"),
    }
    verdict(
        "sm-001",
        "two gauge edits auto-reconcile; the nominal post-norm is checked under tol \
         and the recovered gauge matches both branches",
    );
}

#[test]
fn sm_002_retained_ring_witness_escalates_candidate_with_parent_labels() {
    let sk = ring();
    let base = vec![0.0; 4];
    // X and Y reinforce the same oriented circulation around the cycle.
    let x = branch("agent-x@c7", vec![0.03, 0.03, 0.03, -0.03]);
    let y = branch("agent-y@c9", vec![0.01, 0.01, 0.01, -0.01]);
    let union: Vec<f64> = x
        .mismatch
        .iter()
        .zip(&y.mismatch)
        .zip(&base)
        .map(|((x_value, y_value), base_value)| x_value + y_value - base_value)
        .collect();
    let split = hodge_decompose(&sk, &union).expect("valid ring decomposition");
    assert!(norm_inf(&split.harmonic) > 0.03, "nonzero retained witness");
    assert!(
        norm_inf(
            &sk.d1(&split.harmonic)
                .expect("valid ring harmonic boundary")
        ) < 1e-12,
        "ring witness is closed in the retained skeleton complex"
    );
    assert!(
        norm_inf(
            &sk.d0t(&split.harmonic)
                .expect("valid ring harmonic transpose")
        ) < 1e-12,
        "ring witness is orthogonal to patch-gauge cochains"
    );
    let cycle_pairing =
        split.harmonic[0] + split.harmonic[1] + split.harmonic[2] - split.harmonic[3];
    assert!(
        (cycle_pairing - 0.16).abs() < 1e-12,
        "nonzero cycle pairing witnesses non-exactness: {cycle_pairing}"
    );
    let out = three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6);
    match out {
        MergeOutcome::Conflicted {
            candidate_remainders,
            type_conflicts,
            ..
        } => {
            assert!(type_conflicts.is_empty());
            assert_eq!(candidate_remainders.len(), 1);
            let c = &candidate_remainders[0];
            assert_eq!(
                c.cells.len(),
                4,
                "all four cycle cells exceed this fixture tolerance"
            );
            assert_eq!(
                c.parents,
                ("agent-x@c7".to_string(), "agent-y@c9".to_string()),
                "both caller-supplied parent labels attached"
            );
        }
        other => panic!("the witnessed ring remainder must block auto-merge: {other:?}"),
    }
    verdict(
        "sm-002",
        "the retained ring cochain has a closed non-exact skeleton witness; runtime \
         output remains a candidate conflict with full support and both parents",
    );
}

#[test]
fn sm_003_gauge_edits_resolve_and_cycle_candidates_block_auto_merge() {
    // Property sweep: these seeded gauge-only edits converge and resolve; an
    // injected cycle component blocks automatic merge on the ring.
    let sk = ring();
    let base = vec![0.0; 4];
    let mut state = 0xd00d_u64;
    let mut lcg = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    for trial in 0..20 {
        let gx: Vec<f64> = (0..4).map(|_| 0.1 * lcg()).collect();
        let gy: Vec<f64> = (0..4).map(|_| 0.1 * lcg()).collect();
        let x = branch("x", sk.d0(&gx).expect("valid seeded x coboundary"));
        let y = branch("y", sk.d0(&gy).expect("valid seeded y coboundary"));
        let out = three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6);
        assert!(
            matches!(out, MergeOutcome::Resolved { .. }),
            "gauge-only trial {trial} must resolve"
        );
        // Now inject a cycle component into X.
        let mut mx = sk.d0(&gx).expect("valid tainted x coboundary");
        let eps = 0.02 + 0.05 * lcg().abs();
        for (k, v) in mx.iter_mut().enumerate() {
            *v += if k == 3 { -eps } else { eps };
        }
        let x2 = branch("x", mx);
        let y2 = branch("y", sk.d0(&gy).expect("valid repeated y coboundary"));
        let out2 = three_way_merge(&sk, &base, &x2, &y2, None, 1e-9, 1e-6);
        assert!(
            matches!(
                out2,
                MergeOutcome::Conflicted {
                    ref candidate_remainders,
                    ..
                } if !candidate_remainders.is_empty()
            ),
            "cycle-tainted trial {trial} must remain a candidate conflict"
        );
    }
    verdict(
        "sm-003",
        "20 seeded trials: converged gauge edits resolve; cycle-tainted remainders \
         stay candidate conflicts rather than topology claims",
    );
}

#[test]
fn sm_003b_exact_tree_remainder_is_candidate_not_topology() {
    let n_patches = 32;
    let sk = SheafSkeleton {
        n_patches,
        edges: (0..n_patches - 1).map(|patch| (patch, patch + 1)).collect(),
        triangles: Vec::new(),
    };
    assert_eq!(
        sk.edges.len() + 1,
        sk.n_patches,
        "connected tree has H1 = 0"
    );
    let mut x_potential = vec![0.0; n_patches];
    let mut y_potential = vec![0.0; n_patches];
    for patch in 1..n_patches {
        x_potential[patch] = x_potential[patch - 1] + 0.1 / 31.0;
        y_potential[patch] = y_potential[patch - 1] + 0.05 / 31.0;
    }
    let x = branch(
        "tree-x",
        sk.d0(&x_potential).expect("valid tree-x coboundary"),
    );
    let y = branch(
        "tree-y",
        sk.d0(&y_potential).expect("valid tree-y coboundary"),
    );
    let base = vec![0.0; sk.edges.len()];
    let union: Vec<f64> = x
        .mismatch
        .iter()
        .zip(&y.mismatch)
        .map(|(x_value, y_value)| x_value + y_value)
        .collect();
    let known_gauge: Vec<f64> = x_potential
        .iter()
        .zip(&y_potential)
        .map(|(x_value, y_value)| x_value + y_value)
        .collect();
    assert!(
        norm_inf(
            &apply_gauge(&sk, &union, &known_gauge).expect("valid explicit tree gauge repair")
        ) < 1e-15,
        "an explicit local patch-gauge repair exists"
    );

    match three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6) {
        MergeOutcome::Resolved {
            residual_receipt, ..
        } => assert!(
            residual_receipt.post_norm <= residual_receipt.tol,
            "a better converged solver may legitimately resolve this exact tree"
        ),
        MergeOutcome::Conflicted {
            candidate_remainders,
            type_conflicts,
            ..
        } => {
            assert!(type_conflicts.is_empty());
            assert!(
                !candidate_remainders.is_empty(),
                "bounded solve remainder must be retained only as a candidate"
            );
        }
        other => panic!(
            "the exact tree may resolve or retain an honest candidate, but must not take another path: {other:?}"
        ),
    }
}

#[test]
fn sm_004_sev0_escalates_instead_of_false_resolution() {
    // A COEXACT residue on the triangle: gauge reconciliation cannot
    // reach a passing state (circulation around the triple junction is
    // not a coboundary), and there is NO harmonic space here — the
    // naive path would report "resolved" despite a nominal residual above the
    // requested tolerance. The Sev-0 guard must escalate.
    let sk = triangle();
    let base = vec![0.0; 3];
    let circulation = sk.d1t(&[0.05]).expect("valid Sev-0 circulation fixture");
    let x = branch("x", circulation.clone());
    let y = branch("y", vec![0.0; 3]);
    // Wait: Y unchanged from base triggers the trivial path — perturb Y
    // slightly so the merge genuinely runs.
    let y = BranchState {
        mismatch: sk.d0(&[0.0, 1e-3, 0.0]).expect("valid Sev-0 y coboundary"),
        ..y
    };
    let out = three_way_merge(&sk, &base, &x, &y, None, 1e-6, 1e-6);
    match out {
        MergeOutcome::EscalatedUnresolved {
            post_norm,
            tol,
            fractions,
            confidence,
        } => {
            assert!(post_norm > tol, "the failure is real: {post_norm} > {tol}");
            assert!(
                fractions.1 > 0.5,
                "the residue is predominantly coexact: {fractions:?}"
            );
            assert!(matches!(confidence, Confidence::Normal { .. }));
        }
        other => panic!("Sev-0: must escalate, never report resolved falsely: {other:?}"),
    }
    // And the trivial fast paths themselves.
    let same = branch(
        "s",
        sk.d0(&[0.0, 0.01, 0.0])
            .expect("valid trivial-path coboundary"),
    );
    let t1 = three_way_merge(&sk, &base, &same, &same.clone(), None, 1e-9, 1e-6);
    assert!(matches!(t1, MergeOutcome::Trivial { reason, .. } if reason == "branches identical"));
    let unchanged = branch("u", base.clone());
    let t2 = three_way_merge(&sk, &base, &unchanged, &same, None, 1e-9, 1e-6);
    assert!(
        matches!(t2, MergeOutcome::Trivial { reason, .. } if reason == "X unchanged from base")
    );
    // Confidence must survive escalation too. Adding an isolated patch makes
    // the weighted graph disconnected (lambda_2 = 0) while leaving the
    // triangle circulation fixture intact.
    let low_gap = SheafSkeleton {
        n_patches: 4,
        edges: sk.edges.clone(),
        triangles: sk.triangles.clone(),
    };
    let x_low = branch(
        "x-low",
        low_gap.d1t(&[0.05]).expect("valid low-gap circulation"),
    );
    let y_low = branch(
        "y-low",
        low_gap
            .d0(&[0.0, 1e-3, 0.0, 0.0])
            .expect("valid low-gap coboundary"),
    );
    assert!(matches!(
        three_way_merge(
            &low_gap,
            &base,
            &x_low,
            &y_low,
            None,
            1e-6,
            1e-6,
        ),
        MergeOutcome::EscalatedUnresolved {
            confidence: Confidence::LowGap { gap, .. },
            ..
        } if gap.abs() < f64::EPSILON
    ));
    verdict(
        "sm-004",
        "a coexact residue escalates unresolved (never a false resolution); trivial \
         fast paths report edit equality without a residual receipt",
    );
}

#[test]
fn sm_005_assignment_authority_refusal_and_degraded_gap() {
    let sk = triangle();
    let base = vec![0.0; 3];
    // Pairwise-different values are a useful candidate diagnostic, but without
    // a base assignment map either branch might be unchanged. The merge must
    // refuse rather than manufacture three-way conflict authority.
    let mut x = branch(
        "x",
        sk.d0(&[0.0, 0.01, 0.0])
            .expect("valid assignment x coboundary"),
    );
    x.assignments
        .insert("loadcase/cruise".to_string(), "2.5g".to_string());
    let mut y = branch(
        "y",
        sk.d0(&[0.0, 0.0, 0.01])
            .expect("valid assignment y coboundary"),
    );
    y.assignments
        .insert("loadcase/cruise".to_string(), "3.0g".to_string());
    let candidates = type_conflicts(&x, &y);
    assert_eq!(candidates.len(), 1);
    assert_eq!(candidates[0].key, "loadcase/cruise");
    assert_eq!(candidates[0].x_value, "2.5g");
    assert_eq!(candidates[0].y_value, "3.0g");
    assert!(matches!(
        three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6),
        MergeOutcome::Refused { reason }
            if reason == "base-aware assignment merge is not represented"
    ));
    // Degraded gap: two clusters joined by ONE weak interface — the
    // weighted algebraic connectivity collapses and the merge is
    // flagged low-confidence (R5).
    let barbell = SheafSkeleton {
        n_patches: 6,
        edges: vec![(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5), (2, 3)],
        triangles: vec![(0, 1, 2), (3, 4, 5)],
    };
    let weights = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1e-4];
    let gap = spectral_gap(&barbell, Some(&weights));
    assert!(gap < 1e-3, "weak-link gap is tiny: {gap}");
    let xb = branch(
        "x",
        barbell
            .d0(&[0.0, 0.01, 0.0, 0.0, 0.0, 0.0])
            .expect("valid barbell x coboundary"),
    );
    let yb = branch(
        "y",
        barbell
            .d0(&[0.0, 0.0, 0.0, 0.0, -0.01, 0.0])
            .expect("valid barbell y coboundary"),
    );
    let base_b = vec![0.0; 7];
    let out_b = three_way_merge(&barbell, &base_b, &xb, &yb, Some(&weights), 1e-9, 1e-3);
    match out_b {
        MergeOutcome::Resolved { confidence, .. } => {
            assert!(
                matches!(confidence, Confidence::LowGap { gap, threshold }
                    if gap < threshold),
                "degraded-gap merge must be flagged: {confidence:?}"
            );
        }
        other => panic!("the barbell gauge merge itself resolves: {other:?}"),
    }
    // A healthy complex at the same threshold is Normal.
    let healthy_gap = spectral_gap(&sk, None);
    assert!(
        healthy_gap > 1e-3,
        "triangle is well-coupled: {healthy_gap}"
    );
    let disconnected = SheafSkeleton {
        n_patches: 4,
        edges: vec![(0, 1), (2, 3)],
        triangles: Vec::new(),
    };
    assert!(
        spectral_gap(&disconnected, None).abs() < f64::EPSILON,
        "lambda_2 remains zero when another positive eigenvalue exists"
    );
    let zero_weight_cut = SheafSkeleton {
        n_patches: 3,
        edges: vec![(0, 1), (1, 2)],
        triangles: Vec::new(),
    };
    assert!(
        spectral_gap(&zero_weight_cut, Some(&[1.0, 0.0])).abs() < f64::EPSILON,
        "a zero-weight cut disconnects the weighted graph"
    );
    verdict(
        "sm-005",
        "base-less assignment payloads refuse despite a pairwise-difference candidate; \
         the weak-link barbell merge resolves but is FLAGGED LowGap (R5)",
    );
}

#[test]
fn sm_005a_nonconflicting_assignments_refuse_instead_of_disappearing() {
    let sk = triangle();
    let base = vec![0.0; 3];
    let mismatch = sk
        .d0(&[0.0, 0.01, 0.0])
        .expect("valid nonconflicting assignment coboundary");
    let mut x = branch("x", mismatch.clone());
    x.assignments
        .insert("material/left".to_string(), "steel".to_string());
    let mut y = branch("y", mismatch);
    y.assignments
        .insert("loadcase/cruise".to_string(), "2.5g".to_string());
    assert!(matches!(
        three_way_merge(&sk, &base, &x, &y, None, 1e-9, 1e-6),
        MergeOutcome::Refused { reason }
            if reason == "base-aware assignment merge is not represented"
    ));
}

#[test]
fn sm_005b_nonfinite_or_malformed_inputs_refuse() {
    let sk = triangle();
    let base = vec![0.0; 3];
    let invalid = branch("invalid", vec![f64::NAN, 0.0, 0.0]);
    let valid = branch("valid", vec![0.0, 1e-3, 0.0]);
    assert!(spectral_gap(&sk, Some(&[f64::MAX; 3])).is_nan());
    assert!(matches!(
        three_way_merge(&sk, &base, &invalid, &valid, None, 1e-6, 1e-6),
        MergeOutcome::Refused { reason } if reason == "cochains must be finite"
    ));
    assert!(matches!(
        three_way_merge(
            &sk,
            &base,
            &valid,
            &valid,
            Some(&[1.0, 1.0]),
            1e-6,
            1e-6,
        ),
        MergeOutcome::Refused { reason }
            if reason == "weights must match edges and be finite non-negative values"
    ));
    assert!(matches!(
        three_way_merge(&sk, &base, &valid, &invalid, None, f64::INFINITY, 1e-6),
        MergeOutcome::Refused { .. }
    ));
    assert!(matches!(
        three_way_merge(&sk, &base, &valid, &invalid, None, 1e-6, 0.0),
        MergeOutcome::Refused { reason }
            if reason == "residual tolerance must be finite non-negative and gap threshold finite positive"
    ));
    let duplicate_triangle = SheafSkeleton {
        n_patches: 3,
        edges: sk.edges.clone(),
        triangles: vec![(0, 1, 2), (0, 1, 2)],
    };
    assert!(matches!(
        three_way_merge(
            &duplicate_triangle,
            &base,
            &valid,
            &valid,
            None,
            1e-6,
            1e-6,
        ),
        MergeOutcome::Refused { reason } if reason == "malformed sheaf skeleton"
    ));
    assert!(matches!(
        three_way_merge(&sk, &[0.0; 2], &valid, &valid, None, 1e-6, 1e-6),
        MergeOutcome::Refused { reason } if reason == "cochain length mismatch"
    ));

    let merged_overflow_x = branch("overflow-x", vec![f64::MAX, 0.0, 0.0]);
    let merged_overflow_y = branch("overflow-y", vec![f64::MAX, 1.0, 0.0]);
    assert!(matches!(
        three_way_merge(
            &sk,
            &base,
            &merged_overflow_x,
            &merged_overflow_y,
            None,
            1e-6,
            1e-6,
        ),
        MergeOutcome::Refused { reason }
            if reason == "merged cochain arithmetic is non-finite"
    ));

    let split_overflow_x = branch("split-overflow-x", vec![1e200, 0.0, 0.0]);
    let split_overflow_y = branch("split-overflow-y", vec![0.0, 1e200, 0.0]);
    assert!(matches!(
        three_way_merge(
            &sk,
            &base,
            &split_overflow_x,
            &split_overflow_y,
            None,
            1e-6,
            1e-6,
        ),
        MergeOutcome::Refused { reason }
            if reason == "decomposition arithmetic is non-finite"
    ));

    assert!(matches!(
        three_way_merge(
            &sk,
            &base,
            &valid,
            &branch("other", vec![1e-3, 0.0, 0.0]),
            Some(&[f64::MAX; 3]),
            1e-6,
            1e-6,
        ),
        MergeOutcome::Refused { reason }
            if reason == "spectral-gap arithmetic is non-finite"
    ));
}

#[test]
fn sm_006_candidate_diagnostic_harness() {
    // Candidate-remainder rate over seeded gauge-dominated edits. This is one
    // diagnostic input, not the full Proposal 10 unresolved-merge criterion.
    let ring4 = ring();
    let rate_ring = candidate_remainder_conflict_rate(&ring4, 60, 0.1, 0xfeed)
        .expect("well-formed seeded ring diagnostic");
    // Gauge-dominated edits on a cycle-bearing complex may leave a small
    // fixed-iteration candidate remainder. This fixture threshold catches
    // regressions in that narrow diagnostic; it is not the Proposal 10 kill
    // criterion.
    assert!(
        rate_ring < 0.25,
        "gauge-dominated merges must rarely retain candidate conflicts: {rate_ring}"
    );
    let tri = triangle();
    let rate_tri = candidate_remainder_conflict_rate(&tri, 60, 0.1, 0xfeed)
        .expect("well-formed seeded triangle diagnostic");
    assert!(
        rate_tri.abs() < f64::EPSILON,
        "no retained candidate conflicts in this triangle fixture: {rate_tri}"
    );
    println!(
        "{{\"metric\":\"candidate-remainder-conflict-rate\",\"ring\":{rate_ring:.3},\"triangle\":{rate_tri:.3},\
         \"fixture_threshold\":0.25}}"
    );
    verdict(
        "sm-006",
        "the candidate-conflict diagnostic runs: triangle rate 0 and ring rate below \
         this fixture's regression threshold",
    );

    assert_eq!(
        candidate_remainder_conflict_rate(&ring4, 0, 0.1, 0xfeed),
        Err(CandidateRateError::ZeroTrials)
    );
    assert_eq!(
        candidate_remainder_conflict_rate(&ring4, 1, f64::NAN, 0xfeed),
        Err(CandidateRateError::InvalidEditScale)
    );
    let malformed = SheafSkeleton {
        n_patches: 2,
        edges: vec![(1, 0)],
        triangles: Vec::new(),
    };
    assert_eq!(
        candidate_remainder_conflict_rate(&malformed, 1, 0.1, 0xfeed),
        Err(CandidateRateError::MalformedSkeleton)
    );
    assert!(matches!(
        candidate_remainder_conflict_rate(&ring4, 1, 1e200, 0xfeed),
        Err(CandidateRateError::TrialRefused { .. })
    ));
}

#[test]
#[allow(clippy::too_many_lines)] // One parity matrix across every successful merge verdict.
fn sm_007_bounded_merge_matches_legacy_outcomes_and_replays() {
    let budget = bounded_budget(400);

    let (triangle_raw, triangle_admitted) = canonical_triangle();
    let base_triangle = vec![0.0; 3];
    let resolved_x = branch(
        "agent-x@c1",
        triangle_admitted
            .d0(&[0.0, 0.02, 0.0])
            .expect("valid canonical x coboundary"),
    );
    let resolved_y = branch(
        "agent-y@c2",
        triangle_admitted
            .d0(&[0.0, 0.0, -0.015])
            .expect("valid canonical y coboundary"),
    );
    let legacy_resolved = three_way_merge(
        &triangle_raw,
        &base_triangle,
        &resolved_x,
        &resolved_y,
        None,
        1e-9,
        1e-6,
    );
    let bounded_resolved: BoundedMergeOutcome = with_cx(|cx| {
        three_way_merge_bounded(
            &triangle_admitted,
            &base_triangle,
            &resolved_x,
            &resolved_y,
            None,
            1e-9,
            1e-6,
            budget,
            cx,
        )
        .expect("bounded resolved fixture")
    });
    assert_eq!(bounded_resolved.outcome, legacy_resolved);
    assert!(matches!(
        bounded_resolved.outcome,
        MergeOutcome::Resolved { .. }
    ));
    assert!(bounded_resolved.usage.execution.work_items > 0);
    assert!(bounded_resolved.usage.spectral_sweeps_completed > 0);
    assert!(bounded_resolved.usage.reserved_output_bytes > 0);
    let replay = with_cx(|cx| {
        three_way_merge_bounded(
            &triangle_admitted,
            &base_triangle,
            &resolved_x,
            &resolved_y,
            None,
            1e-9,
            1e-6,
            budget,
            cx,
        )
        .expect("deterministic bounded replay")
    });
    assert_eq!(replay, bounded_resolved);

    let same = branch(
        "same",
        triangle_admitted
            .d0(&[0.0, 0.01, 0.0])
            .expect("valid identical coboundary"),
    );
    let legacy_trivial = three_way_merge(
        &triangle_raw,
        &base_triangle,
        &same,
        &same,
        None,
        1e-9,
        1e-6,
    );
    let bounded_trivial = with_cx(|cx| {
        three_way_merge_bounded(
            &triangle_admitted,
            &base_triangle,
            &same,
            &same,
            None,
            1e-9,
            1e-6,
            budget,
            cx,
        )
        .expect("bounded trivial fixture")
    });
    assert_eq!(bounded_trivial.outcome, legacy_trivial);
    assert_eq!(bounded_trivial.usage.spectral_sweeps_completed, 0);

    let unchanged = branch("base", base_triangle.clone());
    for (expected_reason, fast_x, fast_y) in [
        ("X unchanged from base", &unchanged, &same),
        ("Y unchanged from base", &same, &unchanged),
    ] {
        let legacy = three_way_merge(
            &triangle_raw,
            &base_triangle,
            fast_x,
            fast_y,
            None,
            1e-9,
            1e-6,
        );
        let bounded = with_cx(|cx| {
            three_way_merge_bounded(
                &triangle_admitted,
                &base_triangle,
                fast_x,
                fast_y,
                None,
                1e-9,
                1e-6,
                budget,
                cx,
            )
            .expect("bounded unchanged-branch fixture")
        });
        assert_eq!(bounded.outcome, legacy);
        assert!(matches!(
            bounded.outcome,
            MergeOutcome::Trivial { reason, .. } if reason == expected_reason
        ));
    }

    let mut unsupported = same.clone();
    unsupported
        .assignments
        .insert("loadcase/cruise".to_string(), "2.5g".to_string());
    with_cx(|cx| {
        assert_eq!(
            three_way_merge_bounded(
                &triangle_admitted,
                &base_triangle,
                &unsupported,
                &same,
                None,
                1e-9,
                1e-6,
                budget,
                cx,
            ),
            Err(SheafMergeError::AssignmentsUnsupported)
        );
    });

    let circulation = triangle_admitted
        .d1t(&[0.05])
        .expect("valid canonical circulation");
    let escalated_x = branch("x", circulation);
    let escalated_y = branch(
        "y",
        triangle_admitted
            .d0(&[0.0, 1e-3, 0.0])
            .expect("valid canonical perturbation"),
    );
    let legacy_escalated = three_way_merge(
        &triangle_raw,
        &base_triangle,
        &escalated_x,
        &escalated_y,
        None,
        1e-6,
        1e-6,
    );
    let bounded_escalated = with_cx(|cx| {
        three_way_merge_bounded(
            &triangle_admitted,
            &base_triangle,
            &escalated_x,
            &escalated_y,
            None,
            1e-6,
            1e-6,
            budget,
            cx,
        )
        .expect("bounded escalation fixture")
    });
    assert_eq!(bounded_escalated.outcome, legacy_escalated);

    let (ring_raw, ring_admitted) = canonical_ring();
    let base_ring = vec![0.0; 4];
    let conflict_x = branch("agent-x@c7", vec![0.03, -0.03, 0.03, 0.03]);
    let conflict_y = branch("agent-y@c9", vec![0.01, -0.01, 0.01, 0.01]);
    let legacy_conflict = three_way_merge(
        &ring_raw,
        &base_ring,
        &conflict_x,
        &conflict_y,
        None,
        1e-9,
        1e-6,
    );
    let bounded_conflict = with_cx(|cx| {
        three_way_merge_bounded(
            &ring_admitted,
            &base_ring,
            &conflict_x,
            &conflict_y,
            None,
            1e-9,
            1e-6,
            budget,
            cx,
        )
        .expect("bounded conflict fixture")
    });
    assert_eq!(bounded_conflict.outcome, legacy_conflict);
    assert_eq!(bounded_conflict.usage.conflict_cells, 4);
    assert_eq!(bounded_conflict.usage.provenance_bytes, 20);

    // This weighted star deliberately makes theta² overflow benignly inside
    // Jacobi while t=0, c=1, and s=0 stay finite. The bounded path must retain
    // the legacy finite-result behavior and LowGap classification.
    let (star_raw, star_admitted) = canonical_weighted_star();
    let base_star = vec![0.0; 2];
    let star_x = branch(
        "star-x",
        star_admitted
            .d0(&[0.0, 0.02, 0.0])
            .expect("valid weighted-star x coboundary"),
    );
    let star_y = branch(
        "star-y",
        star_admitted
            .d0(&[0.0, 0.0, -0.015])
            .expect("valid weighted-star y coboundary"),
    );
    let star_weights = [1e-200, 1.0];
    let legacy_star = three_way_merge(
        &star_raw,
        &base_star,
        &star_x,
        &star_y,
        Some(&star_weights),
        1e-9,
        1e-6,
    );
    let bounded_star = with_cx(|cx| {
        three_way_merge_bounded(
            &star_admitted,
            &base_star,
            &star_x,
            &star_y,
            Some(&star_weights),
            1e-9,
            1e-6,
            budget,
            cx,
        )
        .expect("bounded weighted-star fixture")
    });
    assert_eq!(bounded_star.outcome, legacy_star);
    assert!(matches!(
        bounded_star.outcome,
        MergeOutcome::Resolved {
            confidence: Confidence::LowGap { .. },
            ..
        }
    ));

    // A finite heavy edge makes the off-diagonal square +∞, but that sum is
    // only Jacobi's continue/break sentinel. Legacy and bounded paths must
    // still publish the same finite Normal-confidence result.
    let heavy_raw = SheafSkeleton {
        n_patches: 2,
        edges: vec![(0, 1)],
        triangles: Vec::new(),
    };
    let heavy_admitted =
        AdmittedSheafSkeleton::admit(heavy_raw.clone()).expect("canonical heavy edge admits");
    let base_heavy = [0.0];
    let heavy_x = branch(
        "heavy-x",
        heavy_admitted
            .d0(&[0.0, 0.02])
            .expect("valid heavy-edge x coboundary"),
    );
    let heavy_y = branch(
        "heavy-y",
        heavy_admitted
            .d0(&[0.0, -0.015])
            .expect("valid heavy-edge y coboundary"),
    );
    let heavy_weights = [1e200];
    let legacy_heavy = three_way_merge(
        &heavy_raw,
        &base_heavy,
        &heavy_x,
        &heavy_y,
        Some(&heavy_weights),
        1e-9,
        1e-6,
    );
    let bounded_heavy = with_cx(|cx| {
        three_way_merge_bounded(
            &heavy_admitted,
            &base_heavy,
            &heavy_x,
            &heavy_y,
            Some(&heavy_weights),
            1e-9,
            1e-6,
            budget,
            cx,
        )
        .expect("bounded heavy-edge fixture")
    });
    assert_eq!(bounded_heavy.outcome, legacy_heavy);
    assert!(matches!(
        bounded_heavy.outcome,
        MergeOutcome::Resolved {
            confidence: Confidence::Normal { gap },
            ..
        } if gap.is_finite()
    ));

    verdict(
        "sm-007",
        "one-accountant bounded merge exactly replays legacy trivial, resolved, conflicted, escalated, benign weighted-overflow, Normal, and LowGap numerical verdicts over identical canonical coordinates; invalid assignments remain typed errors",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Exact-boundary success and cap-minus-one refusal table.
fn sm_008_bounded_merge_refuses_each_undersized_envelope() {
    let (_, skeleton) = canonical_triangle();
    let base = vec![0.0; 3];
    let x = branch(
        "x",
        skeleton
            .d0(&[0.0, 0.02, 0.0])
            .expect("valid budget x coboundary"),
    );
    let y = branch(
        "y",
        skeleton
            .d0(&[0.0, 0.0, -0.015])
            .expect("valid budget y coboundary"),
    );
    let budget = bounded_budget(8);
    let baseline = with_cx(|cx| {
        three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx)
            .expect("baseline bounded envelope")
    });

    let required_operators = baseline.usage.execution.operator_evaluations;
    assert!(required_operators > 0);
    let required_work = baseline.usage.execution.admitted_work_items;
    let required_scalars = baseline.usage.execution.admitted_scalar_slots;
    let required_output = baseline.usage.admitted_output_bytes;
    let exact_budget = SheafMergeBudget {
        repair: SheafRepairBudget {
            max_operator_evaluations: required_operators,
            max_work_items: required_work,
            ..budget.repair
        },
        max_scalar_slots: required_scalars,
        max_output_bytes: required_output,
        ..budget
    };
    let exact = with_cx(|cx| {
        three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, exact_budget, cx)
            .expect("every exact admission boundary succeeds")
    });
    assert_eq!(exact.outcome, baseline.outcome);
    assert_eq!(exact.usage.execution.admitted_work_items, required_work);
    assert_eq!(
        exact.usage.execution.admitted_scalar_slots,
        required_scalars
    );
    assert_eq!(exact.usage.admitted_output_bytes, required_output);

    with_cx(|cx| {
        assert!(matches!(
            three_way_merge_bounded(
                &skeleton,
                &base,
                &x,
                &y,
                None,
                1e-9,
                1e-6,
                SheafMergeBudget {
                    repair: SheafRepairBudget {
                        max_operator_evaluations: required_operators - 1,
                        ..budget.repair
                    },
                    ..budget
                },
                cx,
            ),
            Err(SheafMergeError::Repair(
                SheafRepairError::WorkBudgetExceeded { required, cap }
            )) if required == required_operators as u128 && cap == required_operators - 1
        ));
    });

    with_cx(|cx| {
        assert!(matches!(
            three_way_merge_bounded(
                &skeleton,
                &base,
                &x,
                &y,
                None,
                1e-9,
                1e-6,
                SheafMergeBudget {
                    repair: SheafRepairBudget {
                        max_work_items: required_work - 1,
                        ..budget.repair
                    },
                    ..budget
                },
                cx,
            ),
            Err(SheafMergeError::Repair(
                SheafRepairError::WorkItemBudgetExceeded {
                    stage: "merge-work-preflight",
                    required,
                    cap,
                }
            )) if required == required_work as u128 && cap == required_work - 1
        ));
    });

    with_cx(|cx| {
        assert_eq!(
            three_way_merge_bounded(
                &skeleton,
                &base,
                &x,
                &y,
                None,
                1e-9,
                1e-6,
                SheafMergeBudget {
                    max_scalar_slots: required_scalars - 1,
                    ..budget
                },
                cx,
            ),
            Err(SheafMergeError::MemoryBudgetExceeded {
                required: required_scalars as u128,
                cap: required_scalars - 1,
            })
        );
    });

    with_cx(|cx| {
        assert_eq!(
            three_way_merge_bounded(
                &skeleton,
                &base,
                &x,
                &y,
                None,
                1e-9,
                1e-6,
                SheafMergeBudget {
                    max_output_bytes: required_output - 1,
                    ..budget
                },
                cx,
            ),
            Err(SheafMergeError::OutputBudgetExceeded {
                resource: "output-bytes",
                required: required_output as u128,
                cap: required_output - 1,
            })
        );
    });

    let huge = AdmittedSheafSkeleton::try_new(4_096, vec![(0, 1)], Vec::new())
        .expect("maximum admitted patch count with one edge");
    let huge_budget = SheafMergeBudget {
        repair: SheafRepairBudget {
            sweeps: 1,
            max_operator_evaluations: 100_000,
            max_work_items: 1_000_000_000,
            max_scalar_slots: 100_000,
            poll_stride: 8,
        },
        max_scalar_slots: 100_000,
        ..budget
    };
    with_cx(|cx| {
        assert_eq!(
            three_way_merge_bounded(
                &huge,
                &[0.0],
                &branch("huge-x", vec![1.0]),
                &branch("huge-y", vec![2.0]),
                None,
                1e-9,
                1e-6,
                huge_budget,
                cx,
            ),
            Err(SheafMergeError::MemoryBudgetExceeded {
                required: 4_096u128 * 4_096,
                cap: 100_000,
            })
        );
    });

    with_cx(|cx| {
        assert_eq!(
            three_way_merge_bounded(
                &skeleton,
                &base,
                &x,
                &y,
                None,
                1e-9,
                1e-6,
                SheafMergeBudget {
                    spectral_sweeps: 0,
                    ..budget
                },
                cx,
            ),
            Err(SheafMergeError::InvalidBudget {
                field: "spectral_sweeps"
            })
        );
        assert_eq!(
            three_way_merge_bounded(
                &skeleton,
                &base,
                &x,
                &y,
                None,
                1e-9,
                1e-6,
                SheafMergeBudget {
                    repair: SheafRepairBudget {
                        poll_stride: 0,
                        ..budget.repair
                    },
                    ..budget
                },
                cx,
            ),
            Err(SheafMergeError::Repair(SheafRepairError::InvalidBudget {
                field: "poll_stride"
            }))
        );
    });

    verdict(
        "sm-008",
        "exact operator, total-work, scalar, and logical-output caps succeed; cap-minus-one, invalid sweep/poll, and maximum-patch dense-spectrum envelopes refuse before allocation or publication",
    );
}

#[test]
fn sm_009_bounded_conflict_caps_are_exact_and_fail_closed() {
    let (_, skeleton) = canonical_ring();
    let base = vec![0.0; 4];
    let x = branch("agent-x@c7", vec![0.03, -0.03, 0.03, 0.03]);
    let y = branch("agent-y@c9", vec![0.01, -0.01, 0.01, 0.01]);
    let budget = SheafMergeBudget {
        max_conflict_cells: 4,
        max_provenance_bytes: 20,
        ..bounded_budget(400)
    };
    let baseline = with_cx(|cx| {
        three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx)
            .expect("exact conflict output caps")
    });
    assert_eq!(baseline.usage.conflict_cells, 4);
    assert_eq!(baseline.usage.provenance_bytes, 20);

    with_cx(|cx| {
        assert_eq!(
            three_way_merge_bounded(
                &skeleton,
                &base,
                &x,
                &y,
                None,
                1e-9,
                1e-6,
                SheafMergeBudget {
                    max_conflict_cells: 3,
                    ..budget
                },
                cx,
            ),
            Err(SheafMergeError::OutputBudgetExceeded {
                resource: "conflict-cells",
                required: 4,
                cap: 3,
            })
        );
        assert_eq!(
            three_way_merge_bounded(
                &skeleton,
                &base,
                &x,
                &y,
                None,
                1e-9,
                1e-6,
                SheafMergeBudget {
                    max_provenance_bytes: 19,
                    ..budget
                },
                cx,
            ),
            Err(SheafMergeError::OutputBudgetExceeded {
                resource: "provenance-bytes",
                required: 20,
                cap: 19,
            })
        );
    });

    verdict(
        "sm-009",
        "candidate conflict cells and UTF-8 parent bytes publish exactly at cap and refuse without truncation at cap minus one",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Pre/mid/final cancellation and deterministic retry transaction.
fn sm_010_bounded_merge_cancellation_and_ambient_refusal_publish_nothing() {
    let (_, skeleton) = canonical_triangle();
    let base = vec![0.0; 3];
    let x = branch(
        "x",
        skeleton
            .d0(&[0.0, 0.02, 0.0])
            .expect("valid cancellation x coboundary"),
    );
    let y = branch(
        "y",
        skeleton
            .d0(&[0.0, 0.0, -0.015])
            .expect("valid cancellation y coboundary"),
    );
    let budget = bounded_budget(8);

    let gate = CancelGate::new();
    gate.request();
    with_gate_cx(&gate, |cx| {
        assert_eq!(
            three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx,),
            Err(SheafMergeError::Repair(SheafRepairError::Cancelled {
                stage: "merge-admission",
                completed_sweeps: 0,
                operator_evaluations: 0,
                work_items: 0,
            }))
        );
    });

    let generous = Budget {
        deadline: None,
        poll_quota: 1_000_000,
        cost_quota: None,
        priority: 0,
    };
    let baseline = with_budget_cx(generous, |cx| {
        three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx)
            .expect("generous bounded merge")
    });
    let planned_cost = baseline.usage.execution.ambient_budget.planned_cost;
    let cost_refusal = Budget {
        cost_quota: Some(planned_cost - 1),
        ..generous
    };
    assert!(matches!(
        with_budget_cx(cost_refusal, |cx| {
            three_way_merge_bounded(
                &skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx,
            )
        }),
        Err(SheafMergeError::Repair(
            SheafRepairError::AmbientBudgetRefused {
                refusal: BudgetRefusal::CostPlanExceedsQuota { planned, quota },
                completed_sweeps: 0,
                operator_evaluations: 0,
                work_items: 0,
            }
        )) if planned == planned_cost && quota == planned_cost - 1
    ));

    let final_poll = baseline.usage.execution.ambient_budget.polls_used;
    assert!(final_poll > 1);
    let final_refusal = Budget {
        poll_quota: final_poll - 1,
        ..generous
    };
    assert!(matches!(
        with_budget_cx(final_refusal, |cx| {
            three_way_merge_bounded(
                &skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx,
            )
        }),
        Err(SheafMergeError::Repair(
            SheafRepairError::AmbientBudgetRefused {
                refusal: BudgetRefusal::PollsExhausted {
                    phase: "merge-publication",
                    quota,
                },
                ..
            }
        )) if quota == final_poll - 1
    ));
    let retry = with_budget_cx(generous, |cx| {
        three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx)
            .expect("healthy retry after publication refusal")
    });
    assert_eq!(retry, baseline);

    let timed_budget = Budget {
        deadline: Some(Time::MAX),
        ..generous
    };
    let timed_gate = CancelGate::new();
    let healthy_clock = RequestingClock::new(&timed_gate, usize::MAX);
    let timed_baseline = with_clock_cx(&timed_gate, timed_budget, &healthy_clock, |cx| {
        three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx)
            .expect("timed cancellation baseline")
    });
    let timed_polls = timed_baseline.usage.execution.ambient_budget.polls_used as usize;
    assert!(timed_polls > 2);
    assert_eq!(healthy_clock.reads(), timed_polls + 1);

    let mid_gate = CancelGate::new();
    let mid_request_on = (timed_polls / 2).max(2).min(timed_polls - 1);
    let mid_clock = RequestingClock::new(&mid_gate, mid_request_on);
    let mid_refusal = with_clock_cx(&mid_gate, timed_budget, &mid_clock, |cx| {
        three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx)
    });
    assert!(matches!(
        mid_refusal,
        Err(SheafMergeError::Repair(SheafRepairError::Cancelled {
            work_items,
            ..
        })) if work_items > 0 && work_items < timed_baseline.usage.execution.work_items
    ));

    // Admission reads the clock once; requesting on read P occurs during the
    // penultimate successful checkpoint, so the final publication boundary is
    // the first point that observes cancellation.
    let publication_gate = CancelGate::new();
    let publication_clock = RequestingClock::new(&publication_gate, timed_polls);
    let publication_refusal =
        with_clock_cx(&publication_gate, timed_budget, &publication_clock, |cx| {
            three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx)
        });
    assert!(matches!(
        publication_refusal,
        Err(SheafMergeError::Repair(SheafRepairError::Cancelled {
            stage: "merge-publication",
            ..
        }))
    ));
    assert_eq!(publication_clock.reads(), timed_polls);

    let timed_retry_gate = CancelGate::new();
    let retry_clock = RequestingClock::new(&timed_retry_gate, usize::MAX);
    let timed_retry = with_clock_cx(&timed_retry_gate, timed_budget, &retry_clock, |cx| {
        three_way_merge_bounded(&skeleton, &base, &x, &y, None, 1e-9, 1e-6, budget, cx)
            .expect("healthy retry after deterministic cancellation")
    });
    assert_eq!(timed_retry, timed_baseline);

    verdict(
        "sm-010",
        "pre-cancel, deterministic mid-run cancellation, true final-publication cancellation, ambient cost admission, and final poll exhaustion return no outcome; fresh retries reproduce both baselines",
    );
}

#[test]
fn sm_011_bounded_candidate_rate_matches_raw_seeded_parity_and_replays() {
    const TRIALS: usize = 3;
    const SCALE: f64 = 0.1;
    let (triangle_raw, triangle_admitted) = canonical_triangle();
    let (ring_raw, ring_admitted) = canonical_ring();

    for (raw, admitted, seed) in [
        (triangle_raw, triangle_admitted, 0xfeed),
        (ring_raw, ring_admitted, 0xdecafbad),
    ] {
        let legacy = candidate_remainder_conflict_rate(&raw, TRIALS, SCALE, seed)
            .expect("canonical raw seeded diagnostic");
        let budget = candidate_budget(TRIALS, 400);
        let bounded: BoundedCandidateRate = with_cx(|cx| {
            candidate_remainder_conflict_rate_bounded(&admitted, TRIALS, SCALE, seed, budget, cx)
                .expect("bounded seeded diagnostic")
        });
        assert_eq!(bounded.rate.to_bits(), legacy.to_bits());
        assert_eq!(bounded.seed, seed);
        assert_eq!(bounded.usage.trials_completed, TRIALS);
        assert!(bounded.usage.conflicts <= TRIALS);
        #[allow(clippy::cast_precision_loss)]
        let measured = bounded.usage.conflicts as f64 / TRIALS as f64;
        assert_eq!(bounded.rate.to_bits(), measured.to_bits());
        assert!(bounded.usage.requested_output_bytes <= bounded.usage.admitted_output_bytes);

        let replay = with_cx(|cx| {
            candidate_remainder_conflict_rate_bounded(&admitted, TRIALS, SCALE, seed, budget, cx)
                .expect("fresh deterministic seeded retry")
        });
        assert_eq!(replay, bounded);
    }

    verdict(
        "sm-011",
        "one-accountant bounded seeded diagnostics match canonical raw rate bits and replay the complete aggregate report",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Exact aggregate relations and cap-minus-one refusal table.
fn sm_012_bounded_candidate_rate_aggregate_caps_are_exact() {
    const TRIALS: usize = 3;
    const SCALE: f64 = 0.1;
    const SEED: u64 = 0x5eed_cafe;
    let (_, skeleton) = canonical_triangle();
    let single_budget = candidate_budget(1, 8);
    let single = with_cx(|cx| {
        candidate_remainder_conflict_rate_bounded(&skeleton, 1, SCALE, SEED, single_budget, cx)
            .expect("single-trial aggregate baseline")
    });
    let generous = candidate_budget(TRIALS, 8);
    with_cx(|cx| {
        assert_eq!(
            candidate_remainder_conflict_rate_bounded(&skeleton, 0, SCALE, SEED, generous, cx,),
            Err(CandidateRateError::ZeroTrials)
        );
        assert_eq!(
            candidate_remainder_conflict_rate_bounded(
                &skeleton,
                TRIALS,
                f64::NAN,
                SEED,
                generous,
                cx,
            ),
            Err(CandidateRateError::InvalidEditScale)
        );
    });
    let multi = with_cx(|cx| {
        candidate_remainder_conflict_rate_bounded(&skeleton, TRIALS, SCALE, SEED, generous, cx)
            .expect("multi-trial aggregate baseline")
    });

    assert_eq!(
        multi.usage.admitted_operator_evaluations,
        single.usage.admitted_operator_evaluations * TRIALS
    );
    assert_eq!(
        multi.usage.execution.operator_evaluations,
        single.usage.execution.operator_evaluations * TRIALS
    );
    let base_work = skeleton.edges().len();
    assert_eq!(
        multi.usage.execution.admitted_work_items,
        base_work + (single.usage.execution.admitted_work_items - base_work) * TRIALS
    );
    assert_eq!(
        multi.usage.execution.admitted_scalar_slots,
        single.usage.execution.admitted_scalar_slots
    );
    assert_eq!(
        multi.usage.admitted_output_bytes,
        single.usage.admitted_output_bytes * TRIALS
    );
    assert_eq!(
        multi.usage.spectral_sweeps_completed,
        single.usage.spectral_sweeps_completed * TRIALS
    );
    assert!(multi.usage.requested_output_bytes > 0);
    assert!(multi.usage.requested_output_bytes <= multi.usage.admitted_output_bytes);

    let required_operators = multi.usage.admitted_operator_evaluations;
    let required_work = multi.usage.execution.admitted_work_items;
    let required_scalars = multi.usage.execution.admitted_scalar_slots;
    let required_output = multi.usage.admitted_output_bytes;
    let exact = CandidateRateBudget {
        max_trials: TRIALS,
        max_operator_evaluations: required_operators,
        max_work_items: required_work,
        max_scalar_slots: required_scalars,
        max_total_output_bytes: required_output,
        ..generous
    };
    let exact_outcome = with_cx(|cx| {
        candidate_remainder_conflict_rate_bounded(&skeleton, TRIALS, SCALE, SEED, exact, cx)
            .expect("every exact aggregate boundary succeeds")
    });
    assert_eq!(exact_outcome.rate.to_bits(), multi.rate.to_bits());
    assert_eq!(exact_outcome.usage, multi.usage);

    with_cx(|cx| {
        assert_eq!(
            candidate_remainder_conflict_rate_bounded(
                &skeleton,
                TRIALS,
                SCALE,
                SEED,
                CandidateRateBudget {
                    max_trials: TRIALS - 1,
                    ..generous
                },
                cx,
            ),
            Err(CandidateRateError::TrialBudgetExceeded {
                requested: TRIALS,
                cap: TRIALS - 1,
            })
        );
    });
    with_cx(|cx| {
        assert!(matches!(
            candidate_remainder_conflict_rate_bounded(
                &skeleton,
                TRIALS,
                SCALE,
                SEED,
                CandidateRateBudget {
                    max_operator_evaluations: required_operators - 1,
                    ..generous
                },
                cx,
            ),
            Err(CandidateRateError::Bounded(SheafMergeError::Repair(
                SheafRepairError::WorkBudgetExceeded { required, cap }
            ))) if required == required_operators as u128 && cap == required_operators - 1
        ));
    });
    with_cx(|cx| {
        assert!(matches!(
            candidate_remainder_conflict_rate_bounded(
                &skeleton,
                TRIALS,
                SCALE,
                SEED,
                CandidateRateBudget {
                    max_work_items: required_work - 1,
                    ..generous
                },
                cx,
            ),
            Err(CandidateRateError::Bounded(SheafMergeError::Repair(
                SheafRepairError::WorkItemBudgetExceeded {
                    stage: "candidate-rate-work-preflight",
                    required,
                    cap,
                }
            ))) if required == required_work as u128 && cap == required_work - 1
        ));
    });
    with_cx(|cx| {
        assert_eq!(
            candidate_remainder_conflict_rate_bounded(
                &skeleton,
                TRIALS,
                SCALE,
                SEED,
                CandidateRateBudget {
                    max_scalar_slots: required_scalars - 1,
                    ..generous
                },
                cx,
            ),
            Err(CandidateRateError::Bounded(
                SheafMergeError::MemoryBudgetExceeded {
                    required: required_scalars as u128,
                    cap: required_scalars - 1,
                }
            ))
        );
    });
    with_cx(|cx| {
        assert_eq!(
            candidate_remainder_conflict_rate_bounded(
                &skeleton,
                TRIALS,
                SCALE,
                SEED,
                CandidateRateBudget {
                    max_total_output_bytes: required_output - 1,
                    ..generous
                },
                cx,
            ),
            Err(CandidateRateError::Bounded(
                SheafMergeError::OutputBudgetExceeded {
                    resource: "candidate-rate-output-bytes",
                    required: required_output as u128,
                    cap: required_output - 1,
                }
            ))
        );
    });

    let planned_cost = multi.usage.execution.ambient_budget.planned_cost;
    let cost_exact = Budget {
        deadline: None,
        poll_quota: 1_000_000,
        cost_quota: Some(planned_cost),
        priority: 0,
    };
    with_budget_cx(cost_exact, |cx| {
        candidate_remainder_conflict_rate_bounded(&skeleton, TRIALS, SCALE, SEED, exact, cx)
            .expect("exact ambient aggregate cost succeeds");
    });
    let cost_refusal = Budget {
        cost_quota: Some(planned_cost - 1),
        ..cost_exact
    };
    assert!(matches!(
        with_budget_cx(cost_refusal, |cx| {
            candidate_remainder_conflict_rate_bounded(
                &skeleton, TRIALS, SCALE, SEED, exact, cx,
            )
        }),
        Err(CandidateRateError::Bounded(SheafMergeError::Repair(
            SheafRepairError::AmbientBudgetRefused {
                refusal: BudgetRefusal::CostPlanExceedsQuota { planned, quota },
                completed_sweeps: 0,
                operator_evaluations: 0,
                work_items: 0,
            }
        ))) if planned == planned_cost && quota == planned_cost - 1
    ));

    verdict(
        "sm-012",
        "aggregate operator, work, scalar-peak, total-output, trial, and ambient-cost boundaries succeed exactly and refuse at cap minus one before partial publication",
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Pre/mid/final cancellation and fresh deterministic retry.
fn sm_013_bounded_candidate_rate_cancellation_is_atomic_and_retryable() {
    const TRIALS: usize = 3;
    const SCALE: f64 = 0.1;
    const SEED: u64 = 0xcafe_f00d;
    let (_, skeleton) = canonical_triangle();
    let budget = candidate_budget(TRIALS, 8);

    let gate = CancelGate::new();
    gate.request();
    with_gate_cx(&gate, |cx| {
        assert_eq!(
            candidate_remainder_conflict_rate_bounded(&skeleton, TRIALS, SCALE, SEED, budget, cx,),
            Err(CandidateRateError::Bounded(SheafMergeError::Repair(
                SheafRepairError::Cancelled {
                    stage: "candidate-rate-admission",
                    completed_sweeps: 0,
                    operator_evaluations: 0,
                    work_items: 0,
                }
            )))
        );
    });

    let timed_budget = Budget {
        deadline: Some(Time::MAX),
        poll_quota: 1_000_000,
        cost_quota: None,
        priority: 0,
    };
    let healthy_gate = CancelGate::new();
    let healthy_clock = RequestingClock::new(&healthy_gate, usize::MAX);
    let baseline = with_clock_cx(&healthy_gate, timed_budget, &healthy_clock, |cx| {
        candidate_remainder_conflict_rate_bounded(&skeleton, TRIALS, SCALE, SEED, budget, cx)
            .expect("timed seeded diagnostic baseline")
    });
    let polls = baseline.usage.execution.ambient_budget.polls_used as usize;
    assert!(polls > 2);
    assert_eq!(healthy_clock.reads(), polls + 1);

    let mid_gate = CancelGate::new();
    let mid_request_on = (polls / 2).max(2).min(polls - 1);
    let mid_clock = RequestingClock::new(&mid_gate, mid_request_on);
    let mid_refusal = with_clock_cx(&mid_gate, timed_budget, &mid_clock, |cx| {
        candidate_remainder_conflict_rate_bounded(&skeleton, TRIALS, SCALE, SEED, budget, cx)
    });
    assert!(matches!(
        mid_refusal,
        Err(CandidateRateError::Bounded(SheafMergeError::Repair(
            SheafRepairError::Cancelled { work_items, .. }
        ))) if work_items > 0 && work_items < baseline.usage.execution.work_items
    ));

    let publication_gate = CancelGate::new();
    let publication_clock = RequestingClock::new(&publication_gate, polls);
    let publication_refusal =
        with_clock_cx(&publication_gate, timed_budget, &publication_clock, |cx| {
            candidate_remainder_conflict_rate_bounded(&skeleton, TRIALS, SCALE, SEED, budget, cx)
        });
    assert!(matches!(
        publication_refusal,
        Err(CandidateRateError::Bounded(SheafMergeError::Repair(
            SheafRepairError::Cancelled {
                stage: "candidate-rate-publication",
                operator_evaluations,
                work_items,
                ..
            }
        ))) if operator_evaluations == baseline.usage.execution.operator_evaluations
            && work_items == baseline.usage.execution.work_items
    ));
    assert_eq!(publication_clock.reads(), polls);

    let retry_gate = CancelGate::new();
    let retry_clock = RequestingClock::new(&retry_gate, usize::MAX);
    let retry = with_clock_cx(&retry_gate, timed_budget, &retry_clock, |cx| {
        candidate_remainder_conflict_rate_bounded(&skeleton, TRIALS, SCALE, SEED, budget, cx)
            .expect("healthy retry after deterministic candidate cancellation")
    });
    assert_eq!(retry, baseline);

    verdict(
        "sm-013",
        "one-accountant seeded diagnostics publish no rate on pre, mid, or true final cancellation and a fresh retry reproduces the complete report",
    );
}
