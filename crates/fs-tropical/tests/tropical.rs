//! Battery for tropical critical-path analytics (fs-tropical). Covers the
//! max-plus semiring, CPM critical path + slack + tune-next recommendations,
//! the tropical eigenvalue (Karp) cross-checked against brute-force cycle
//! enumeration, cycle detection, and the promotion audit.

use fs_tropical::{
    AuditEntry, NEG_INF, RecommendationAudit, TaskDag, TropicalError, brute_force_max_cycle,
    critical_circuit, max_cycle_mean, oplus, otimes,
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
    let cp = dag.critical_path().unwrap();
    // critical tasks by latency: 3 (4) > 1 (3) > 0 (2); task 2 is off-path.
    assert_eq!(dag.tune_next(&cp), vec![3, 1, 0]);
    // the true bottleneck is the highest-latency critical task.
    assert_eq!(dag.bottleneck(&cp), Some(3));
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
    assert!(cp.path.is_empty() && cp.slack.is_empty());
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
fn the_tropical_eigenvalue_matches_brute_force() {
    // 0→1 (3), 1→0 (5): the max cycle mean is (3+5)/2 = 4.
    let m = vec![vec![NEG_INF, 3.0], vec![5.0, NEG_INF]];
    let karp = max_cycle_mean(&m).unwrap();
    assert!((karp - 4.0).abs() < 1e-9);
    let (cycle, mean) = brute_force_max_cycle(&m).unwrap();
    assert!((karp - mean).abs() < 1e-9);
    assert_eq!(cycle, vec![0, 1]);
    assert_eq!(critical_circuit(&m), Some(vec![0, 1]));
}

#[test]
fn a_faster_cycle_dominates_a_self_loop() {
    // 0→0 (2 self-loop), 0→1 (3), 1→0 (5): the 2-cycle (mean 4) beats the loop (2).
    let m = vec![vec![2.0, 3.0], vec![5.0, NEG_INF]];
    assert!((max_cycle_mean(&m).unwrap() - 4.0).abs() < 1e-9);
    // a 3-cycle mean.
    let tri = vec![
        vec![NEG_INF, 1.0, NEG_INF],
        vec![NEG_INF, NEG_INF, 2.0],
        vec![3.0, NEG_INF, NEG_INF],
    ];
    assert!((max_cycle_mean(&tri).unwrap() - 2.0).abs() < 1e-9); // (1+2+3)/3
    // an acyclic matrix has no cycle.
    let acyclic = vec![vec![NEG_INF, 1.0], vec![NEG_INF, NEG_INF]];
    assert_eq!(max_cycle_mean(&acyclic), None);
}

#[test]
fn the_promotion_audit_gates_on_realized_outcomes() {
    let mut audit = RecommendationAudit::new();
    audit.record(AuditEntry {
        recommended: 3,
        predicted_gain: 1.0,
        realized_gain: 0.9,
    });
    audit.record(AuditEntry {
        recommended: 1,
        predicted_gain: 0.5,
        realized_gain: 0.4,
    });
    audit.record(AuditEntry {
        recommended: 0,
        predicted_gain: 0.2,
        realized_gain: -0.1,
    }); // miss
    assert!((audit.hit_rate() - 2.0 / 3.0).abs() < 1e-9);
    // enough samples + a good hit rate -> promotable; a stricter bar is not met.
    assert!(audit.promoted(3, 0.5));
    assert!(!audit.promoted(3, 0.8));
    assert!(!RecommendationAudit::new().promoted(3, 0.5)); // no data
}
