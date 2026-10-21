use super::*;
use super::super::{DifferentialEvolutionState, DifferentialEvolutionStop};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use std::{cell::Cell, convert::Infallible};
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};

fn config(budget: usize) -> DifferentialEvolutionConfig {
    DifferentialEvolutionConfig { population_size: 8, differential_weight: 0.8,
        crossover_probability: 0.9, seed: 17, max_evaluations: budget,
        max_population_entries: 16, target: None }
}
fn state(budget: usize) -> ConstrainedDifferentialEvolutionState {
    ConstrainedDifferentialEvolutionState::new_with(&[[0.0, 1.0]; 2], Some(&[0.0; 2]),
        config(budget), 2, 0.0, &mut || true).unwrap()
}
fn evaluation(x: &[f64]) -> ConstrainedEvaluation {
    ConstrainedEvaluation { objective: x[0]*x[0] + x[1]*x[1],
        inequalities: vec![0.6-x[0], x[0]+x[1]-1.2] }
}
fn with_cx<T>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> T) -> T {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&gate, &Cx::new(
        &gate, arena, StreamKey { seed: 17, kernel_id: 55, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}

#[test]
fn g0_feasibility_dominates_every_objective_scale_without_a_penalty() {
    let feasible = score(f64::MAX, &[0.0, -f64::MAX], 2, 0.0).unwrap();
    let infeasible = score(-f64::MAX, &[f64::MIN_POSITIVE, -1.0], 2, 0.0).unwrap();
    assert!(feasible.rank(infeasible).is_lt());
    assert!(!infeasible.reaches(f64::MAX));
    let small = score(100.0, &[0.1, -1.0], 2, 0.0).unwrap();
    let large = score(-100.0, &[0.2, 0.2], 2, 0.0).unwrap();
    assert!(small.rank(large).is_lt());
    let tie = score(99.0, &[0.1, 0.05], 2, 0.0).unwrap();
    assert!(tie.rank(small).is_lt());
    assert!(tie.accepts(tie));
    let admitted = score(1.0, &[0.01], 1, 0.01).unwrap();
    assert!(admitted.is_feasible());
    assert_eq!(admitted.max_scaled_violation(), 0.01);
    assert!(admitted.reaches(1.0));
    let extreme = score(1.0, &[f64::MAX, f64::MAX], 2, 0.0).unwrap();
    assert_eq!(extreme.max_scaled_violation(), f64::MAX);
}

#[test]
fn g1_search_finds_a_feasible_boundary_but_never_calls_impossible_targets_recovered() {
    let mut search = state(6000);
    let report = search.advance_with(6000, &mut |x| Ok::<_, Infallible>(evaluation(x)), &mut || true).unwrap();
    let best = report.best.unwrap();
    assert!(best.value.is_feasible());
    assert!(best.coordinates[0] >= 0.6 && best.coordinates[0] < 0.6001);
    assert!(best.coordinates[1] < 0.001);
    assert!((best.value.objective()-0.36).abs() < 0.0002);
    let mut cfg = config(128); cfg.target = Some(1.0);
    let mut impossible = ConstrainedDifferentialEvolutionState::new_with(&[[0.0, 1.0]],
        Some(&[0.0]), cfg, 1, 0.0, &mut || true).unwrap();
    let report = impossible.advance_with(128, &mut |_| Ok::<_, Infallible>(ConstrainedEvaluation {
        objective: 0.0, inequalities: vec![1.0] }), &mut || true).unwrap();
    assert_eq!(report.stop, DifferentialEvolutionStop::EvaluationLimit);
    assert!(!report.best.unwrap().value.is_feasible());
}

#[test]
fn g3_all_feasible_search_matches_scalar_trace_and_resume_keeps_hard_caps() {
    let mut scalar = DifferentialEvolutionState::new_with(&[[0.0, 1.0]; 2], Some(&[0.0; 2]),
        config(128), &mut || true).unwrap();
    let mut expected = Vec::new();
    scalar.advance_with(128, &mut |x| { expected.push(x.to_vec());
        Ok::<_, Infallible>((x[0]-0.3).powi(2)+(x[1]-0.7).powi(2)) }, &mut || true).unwrap();
    for width in [1, 3, 8, 17] {
        let mut search = state(0);
        let mut actual = Vec::new();
        for total in 1..=128 {
            search.extend_evaluation_budget(total).unwrap();
            let report = search.advance_with(width, &mut |x| { actual.push(x.to_vec());
                Ok::<_, Infallible>(ConstrainedEvaluation { objective: (x[0]-0.3).powi(2)+(x[1]-0.7).powi(2),
                    inequalities: vec![-1.0, 0.0] }) }, &mut || true).unwrap();
            assert_eq!(report.evaluations, total);
            assert_eq!(actual, expected[..total]);
        }
        assert_eq!(search.best().unwrap().coordinates, scalar.best().unwrap().coordinates);
        assert_eq!(search.best().unwrap().value.objective(), scalar.best().unwrap().value);
    }
}

#[test]
fn g4_cancelled_success_is_retained_and_malformed_or_failed_callbacks_stay_errors() {
    let permit = Cell::new(true);
    let mut search = state(10);
    let result = search.advance_with(1, &mut |x| { permit.set(false);
        Ok::<_, Infallible>(evaluation(x)) }, &mut || permit.get());
    assert!(matches!(result, Err(DifferentialEvolutionRunError::Solver(Error::Cancelled))));
    assert_eq!(search.evaluations(), 1);
    assert!(search.best().is_none());
    let report = search.advance_with(0, &mut |_| -> Result<ConstrainedEvaluation, Infallible> {
        panic!("completed constrained physics was repeated") }, &mut || true).unwrap();
    assert_eq!(report.initialized_members, 1);
    assert!(!report.best.unwrap().value.is_feasible());
    let before = search.population().to_vec();
    let mut attempts = Vec::new();
    for rows in [vec![], vec![0.0, f64::NAN], vec![f64::INFINITY, 0.0]] {
        assert!(search.advance_with(1, &mut |x| { attempts.push(x.to_vec());
            Ok::<_, Infallible>(ConstrainedEvaluation { objective: 0.0, inequalities: rows.clone() })
        }, &mut || true).is_err());
        assert_eq!(search.population(), before);
    }
    let result = search.advance_with(1, &mut |x| { attempts.push(x.to_vec());
        Err::<ConstrainedEvaluation, _>("physical failure") }, &mut || true);
    assert!(matches!(result, Err(DifferentialEvolutionRunError::Evaluation("physical failure"))));
    search.advance_with(1, &mut |x| { attempts.push(x.to_vec());
        Ok::<_, Infallible>(evaluation(x)) }, &mut || true).unwrap();
    assert_eq!(search.evaluations(), 6);
    assert!(attempts.windows(2).all(|p| p[0] == p[1]));
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(ConstrainedDifferentialEvolutionState::new_with(&[[0.0, 1.0]], None,
            config(1), 1, tolerance, &mut || true).is_err());
    }
    assert!(ConstrainedDifferentialEvolutionState::new_with(&[[0.0, 1.0]], None,
        config(1), 0, 0.0, &mut || true).is_err());
}

#[test]
fn g5_constrained_parallel_selection_matches_serial_for_workers_waves_and_pauses() {
    with_cx(|_, cx| {
        let mut serial = state(257);
        let expected = serial.advance(cx, 257, &mut |x| Ok::<_, Infallible>(evaluation(x))).unwrap();
        let lease = OperationMemoryLease::bounded(64 << 20);
        for workers in [1, 4] {
            let pool = TilePool::for_host(workers, 17);
            for wave in [1, 3, 8] {
                let mut parallel = state(37);
                let first = parallel.advance_parallel(cx, &pool, &lease, 37, wave,
                    &|_, x| Ok::<_, Infallible>(evaluation(x))).unwrap();
                assert_eq!(first.evaluations, 37);
                parallel.extend_evaluation_budget(257).unwrap();
                parallel.advance(cx, 11, &mut |x| Ok::<_, Infallible>(evaluation(x))).unwrap();
                let actual = parallel.advance_parallel(cx, &pool, &lease, 1000, wave,
                    &|_, x| Ok::<_, Infallible>(evaluation(x))).unwrap();
                assert_eq!(actual, expected);
                assert_eq!(parallel.population(), serial.population());
            }
        }
    });
}

#[test]
fn g4_failed_parallel_member_preserves_feasible_and_infeasible_siblings() {
    with_cx(|_, cx| {
        let pool = TilePool::for_host(4, 17);
        let lease = OperationMemoryLease::bounded(64 << 20);
        let calls = AtomicUsize::new(0);
        let mut search = state(16);
        let result = search.advance_parallel(cx, &pool, &lease, 8, 8, &|_, x| {
            calls.fetch_add(1, AtomicOrdering::Relaxed);
            if x == [0.0, 0.0] { Err("material refusal") } else { Ok(evaluation(x)) }
        });
        assert!(matches!(result, Err(DifferentialEvolutionParallelError::Evaluation {
            member: 0, source: "material refusal" })));
        assert_eq!(search.evaluations(), 8);
        let report = search.advance(cx, 1, &mut |x| {
            calls.fetch_add(1, AtomicOrdering::Relaxed);
            assert_eq!(x, [0.0, 0.0]);
            Ok::<_, Infallible>(evaluation(x))
        }).unwrap();
        assert_eq!(report.initialized_members, 8);
        assert_eq!(calls.load(AtomicOrdering::Relaxed), 9);
        let mut reference = state(8);
        reference.advance(cx, 8, &mut |x| Ok::<_, Infallible>(evaluation(x))).unwrap();
        assert_eq!(search.population(), reference.population());
        assert_eq!(search.best(), reference.best());
    });
}
