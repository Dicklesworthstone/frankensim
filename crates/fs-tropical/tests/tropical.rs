//! Battery for tropical critical-path analytics (fs-tropical). Covers the
//! max-plus semiring, CPM critical path + slack + tune-next recommendations,
//! the tropical eigenvalue (Karp) cross-checked against brute-force cycle
//! enumeration, cycle detection, and the promotion audit.

use fs_tropical::{
    AuditEntry, MAX_BRUTE_FORCE_NODES, MAX_RECOMMENDATION_AUDIT_ENTRIES, MAX_TASK_DAG_NODES,
    MAX_TROPICAL_MATRIX_NODES, NEG_INF, RecommendationAudit, TaskDag, TropicalError,
    brute_force_max_cycle, critical_circuit, max_cycle_mean, oplus, otimes,
};

#[test]
fn the_max_plus_semiring_operates() {
    assert!((oplus(3.0, 5.0) - 5.0).abs() < 1e-12); // ⊕ = max
    assert!((otimes(3.0, 5.0) - 8.0).abs() < 1e-12); // ⊗ = +
    assert_eq!(otimes(NEG_INF, 5.0).to_bits(), NEG_INF.to_bits()); // −∞ absorbs
}

#[test]
fn the_critical_path_is_the_longest_path() {
    // 0→{1,2}→3; latencies [2,3,1,4]. Path 0-1-3 = 9 dominates 0-2-3 = 7.
    let dag = TaskDag::new(vec![2.0, 3.0, 1.0, 4.0])
        .with_edge(0, 1)
        .with_edge(0, 2)
        .with_edge(1, 3)
        .with_edge(2, 3);
    let cp = dag.critical_path().unwrap();
    assert!((cp.makespan - 9.0).abs() < 1e-12);
    assert!(cp.makespan_lo <= 9.0 && cp.makespan_hi >= 9.0);
    assert!(cp.path_is_unique);
    assert_eq!(cp.path, vec![0, 1, 3]);
    // task 2 is OFF the critical path -> positive slack; the rest are critical.
    assert!((cp.slack[2] - 2.0).abs() < 1e-12);
    for i in [0, 1, 3] {
        assert!(cp.slack[i].abs() < 1e-12);
    }
}

#[test]
fn tune_next_ranks_the_critical_bottleneck_first() {
    let dag = TaskDag::new(vec![2.0, 3.0, 1.0, 4.0])
        .with_edge(0, 1)
        .with_edge(0, 2)
        .with_edge(1, 3)
        .with_edge(2, 3);
    // critical tasks by latency: 3 (4) > 1 (3) > 0 (2); task 2 is off-path.
    assert_eq!(dag.tune_next().unwrap(), vec![3, 1, 0]);
    // the true bottleneck is the highest-latency critical task.
    assert_eq!(dag.bottleneck(), Ok(Some(3)));
}

#[test]
fn a_cyclic_dag_has_no_critical_path() {
    let cyclic = TaskDag::new(vec![1.0, 1.0]).with_edge(0, 1).with_edge(1, 0);
    assert_eq!(cyclic.critical_path(), Err(TropicalError::Cyclic));
}

#[test]
fn an_empty_dag_has_an_empty_critical_path() {
    // regression: must not panic (previously indexed best_pred[0] on an empty Vec).
    let cp = TaskDag::new(vec![]).critical_path().unwrap();
    assert!((cp.makespan - 0.0).abs() < 1e-12);
    assert_eq!((cp.makespan_lo, cp.makespan_hi), (0.0, 0.0));
    assert!(!cp.path_is_unique, "an empty graph has no unique path");
    assert!(cp.path.is_empty() && cp.slack.is_empty());
    assert!(TaskDag::new(vec![]).tune_next().unwrap().is_empty());
    assert_eq!(TaskDag::new(vec![]).bottleneck(), Ok(None));
}

#[test]
fn a_zero_latency_source_stays_on_the_path() {
    // regression: a predecessor with earliest-finish 0 was dropped from the
    // back-traced path by the old strict `ef[p] > es` (es started at 0.0).
    let dag = TaskDag::new(vec![0.0, 10.0, 5.0])
        .with_edge(0, 1)
        .with_edge(1, 2);
    let cp = dag.critical_path().unwrap();
    assert!((cp.makespan - 15.0).abs() < 1e-12);
    assert_eq!(cp.path, vec![0, 1, 2]); // task 0 (zero latency) is retained
}

#[test]
fn a_zero_latency_terminal_stays_on_the_source_to_sink_path() {
    // The global maximum finish is shared by tasks 1 and 0. The witness must
    // nevertheless end at terminal task 0 instead of truncating at task 1.
    let dag = TaskDag::new(vec![0.0, 1.0]).with_edge(1, 0);
    let cp = dag.critical_path().expect("finite zero-tail path");
    assert_eq!(cp.path, vec![1, 0]);
    assert!(cp.path_is_unique);
    assert_eq!(dag.tune_next(), Ok(vec![1]));
    assert_eq!(dag.bottleneck(), Ok(Some(1)));
}

#[test]
fn zero_work_and_tied_top_latencies_have_no_single_bottleneck() {
    let zero = TaskDag::new(vec![0.0]);
    let cp = zero.critical_path().expect("single finite task");
    assert!(cp.path_is_unique);
    assert_eq!(cp.path, vec![0]);
    assert!(zero.tune_next().unwrap().is_empty());
    assert_eq!(zero.bottleneck(), Ok(None));

    let tied = TaskDag::new(vec![2.0, 2.0, 1.0])
        .with_edge(0, 1)
        .with_edge(1, 2);
    assert!(tied.critical_path().unwrap().path_is_unique);
    assert_eq!(tied.tune_next(), Ok(vec![0, 1, 2]));
    assert_eq!(tied.bottleneck(), Ok(None));
}

#[test]
fn floating_makespan_is_outward_bounded_and_ties_have_no_single_target() {
    let fractional = TaskDag::new(vec![0.1, 0.2]).with_edge(0, 1);
    let cp = fractional.critical_path().expect("finite fractional path");
    assert!(cp.makespan_lo < cp.makespan);
    assert!(cp.makespan_hi > cp.makespan);
    assert!(cp.path_is_unique);

    let tied = TaskDag::new(vec![9.0, 1.0, 6.0, 4.0])
        .with_edge(0, 1)
        .with_edge(2, 3);
    let tied_cp = tied.critical_path().expect("two finite paths");
    assert!(!tied_cp.path_is_unique);
    assert!(tied.tune_next().expect("analysis succeeds").is_empty());
    assert_eq!(tied.bottleneck(), Ok(None));
}

#[test]
fn the_tropical_eigenvalue_matches_brute_force() {
    // 0→1 (3), 1→0 (5): the max cycle mean is (3+5)/2 = 4.
    let m = vec![vec![NEG_INF, 3.0], vec![5.0, NEG_INF]];
    let karp = max_cycle_mean(&m).unwrap().unwrap();
    assert!((karp - 4.0).abs() < 1e-9);
    let (cycle, mean) = brute_force_max_cycle(&m).unwrap().unwrap();
    assert!((karp - mean).abs() < 1e-9);
    assert_eq!(cycle, vec![0, 1]);
    assert_eq!(critical_circuit(&m), Ok(Some(vec![0, 1])));

    let close_loops = vec![vec![0.0, NEG_INF], vec![NEG_INF, 5e-13]];
    assert_eq!(
        brute_force_max_cycle(&close_loops)
            .expect("valid matrix")
            .expect("self loops")
            .0,
        vec![1]
    );
}

#[test]
fn a_faster_cycle_dominates_a_self_loop() {
    // 0→0 (2 self-loop), 0→1 (3), 1→0 (5): the 2-cycle (mean 4) beats the loop (2).
    let m = vec![vec![2.0, 3.0], vec![5.0, NEG_INF]];
    assert!((max_cycle_mean(&m).unwrap().unwrap() - 4.0).abs() < 1e-9);
    // a 3-cycle mean.
    let tri = vec![
        vec![NEG_INF, 1.0, NEG_INF],
        vec![NEG_INF, NEG_INF, 2.0],
        vec![3.0, NEG_INF, NEG_INF],
    ];
    assert!((max_cycle_mean(&tri).unwrap().unwrap() - 2.0).abs() < 1e-9); // (1+2+3)/3
    // an acyclic matrix has no cycle.
    let acyclic = vec![vec![NEG_INF, 1.0], vec![NEG_INF, NEG_INF]];
    assert_eq!(max_cycle_mean(&acyclic), Ok(None));
}

#[test]
fn malformed_or_oversized_inputs_fail_closed() {
    for latency in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0] {
        assert!(matches!(
            TaskDag::new(vec![latency]).critical_path(),
            Err(TropicalError::InvalidLatency { node: 0, .. })
        ));
    }
    assert!(matches!(
        TaskDag::new(vec![1.0; MAX_TASK_DAG_NODES + 1]).critical_path(),
        Err(TropicalError::ResourceLimit {
            resource: "task nodes",
            ..
        })
    ));
    assert_eq!(
        TaskDag::new(vec![]).with_edge(0, 0).critical_path(),
        Err(TropicalError::BadEdge { node: 0 })
    );
    assert!(matches!(
        max_cycle_mean(&[vec![0.0], vec![0.0, 0.0]]),
        Err(TropicalError::NonSquareMatrix { row: 0, .. })
    ));
    for invalid in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            max_cycle_mean(&[vec![invalid]]),
            Err(TropicalError::InvalidWeight { row: 0, col: 0 })
        );
    }
    assert!(matches!(
        max_cycle_mean(&vec![
            vec![NEG_INF; MAX_TROPICAL_MATRIX_NODES + 1];
            MAX_TROPICAL_MATRIX_NODES + 1
        ]),
        Err(TropicalError::ResourceLimit {
            resource: "matrix nodes",
            ..
        })
    ));
    assert!(matches!(
        brute_force_max_cycle(&vec![
            vec![NEG_INF; MAX_BRUTE_FORCE_NODES + 1];
            MAX_BRUTE_FORCE_NODES + 1
        ]),
        Err(TropicalError::ResourceLimit {
            resource: "brute-force cycle nodes",
            ..
        })
    ));
    assert!(matches!(
        max_cycle_mean(&[vec![f64::MAX, f64::MAX], vec![f64::MAX, f64::MAX]]),
        Err(TropicalError::NumericalOverflow { .. })
    ));
}

#[test]
fn the_promotion_audit_gates_on_realized_outcomes() {
    let mut audit = RecommendationAudit::new();
    audit
        .record(AuditEntry {
            recommended: 3,
            predicted_gain: 1.0,
            realized_gain: 0.9,
        })
        .unwrap();
    audit
        .record(AuditEntry {
            recommended: 1,
            predicted_gain: 0.5,
            realized_gain: 0.4,
        })
        .unwrap();
    audit
        .record(AuditEntry {
            recommended: 0,
            predicted_gain: 0.2,
            realized_gain: -0.1,
        })
        .unwrap(); // miss
    assert!((audit.hit_rate() - 2.0 / 3.0).abs() < 1e-9);
    // enough samples + a good hit rate -> promotable; a stricter bar is not met.
    assert_eq!(audit.promoted(3, 0.5), Ok(true));
    assert_eq!(audit.promoted(3, 0.8), Ok(false));
    assert_eq!(RecommendationAudit::new().promoted(3, 0.5), Ok(false));
    assert!(matches!(
        audit.promoted(0, 0.0),
        Err(TropicalError::InvalidAuditField {
            field: "min_samples",
            ..
        })
    ));
    assert!(matches!(
        audit.promoted(1, f64::NAN),
        Err(TropicalError::InvalidAuditField {
            field: "min_hit_rate",
            ..
        })
    ));

    let mut bounded = RecommendationAudit::new();
    for _ in 0..MAX_RECOMMENDATION_AUDIT_ENTRIES {
        bounded
            .record(AuditEntry {
                recommended: 0,
                predicted_gain: 0.0,
                realized_gain: 0.0,
            })
            .expect("exact audit boundary");
    }
    assert!(matches!(
        bounded.record(AuditEntry {
            recommended: 0,
            predicted_gain: 0.0,
            realized_gain: 0.0,
        }),
        Err(TropicalError::ResourceLimit {
            resource: "recommendation audit entries",
            ..
        })
    ));
    assert!(matches!(
        RecommendationAudit::new().record(AuditEntry {
            recommended: 0,
            predicted_gain: f64::NAN,
            realized_gain: 0.0,
        }),
        Err(TropicalError::InvalidAuditField {
            field: "predicted_gain",
            ..
        })
    ));
}
