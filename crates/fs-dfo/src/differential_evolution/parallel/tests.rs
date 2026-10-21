use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use std::convert::Infallible;
use std::sync::atomic::AtomicBool;

fn with_cx<T>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> T) -> T {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&gate, &Cx::new(
        &gate, arena, StreamKey { seed: 17, kernel_id: 55, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn state(cx: &Cx<'_>, limit: usize) -> DifferentialEvolutionState {
    DifferentialEvolutionState::new(cx, &[[-5.0, 5.0]; 3], Some(&[0.0; 3]),
        super::super::DifferentialEvolutionConfig { population_size: 8,
            differential_weight: 0.8, crossover_probability: 0.9, seed: 17,
            max_evaluations: limit, max_population_entries: 24, target: None }).unwrap()
}
fn objective(x: &[f64]) -> f64 { x.iter().map(|x| (x-0.25).powi(2)).sum() }

#[test]
fn g5_parallel_waves_match_serial_populations_at_every_budget_and_worker_count() {
    with_cx(|_, cx| {
        for threads in [1, 2, 4] {
            let pool = TilePool::for_host(threads, 17);
            let lease = OperationMemoryLease::bounded(64 << 20);
            for width in [1, 3, 8, 17] {
                let mut serial = state(cx, 37);
                let expected = serial.advance(cx, 37, &mut |x| Ok::<_, Infallible>(objective(x))).unwrap();
                let mut parallel = state(cx, 37);
                let actual = parallel.advance_parallel(cx, &pool, &lease, 37, width,
                    &|_, x| Ok::<_, Infallible>(objective(x))).unwrap();
                assert_eq!(actual, expected);
                assert_eq!(parallel.population(), serial.population());
                for total in [38, 48, 59, 64] {
                    parallel.extend_evaluation_budget(total).unwrap();
                    serial.extend_evaluation_budget(total).unwrap();
                    let a = parallel.advance_parallel(cx, &pool, &lease, usize::MAX, width,
                        &|_, x| Ok::<_, Infallible>(objective(x))).unwrap();
                    let b = serial.advance(cx, usize::MAX, &mut |x| Ok::<_, Infallible>(objective(x))).unwrap();
                    assert_eq!(a, b);
                    assert_eq!(parallel.population(), serial.population());
                }
            }
        }
    });
}

#[test]
fn g0_parallel_callback_limits_include_partial_initialization_and_generation() {
    with_cx(|_, cx| {
        let pool = TilePool::for_host(4, 17);
        let lease = OperationMemoryLease::bounded(64 << 20);
        for limit in 0..=33 {
            let calls = AtomicUsize::new(0);
            let mut search = state(cx, limit);
            let report = search.advance_parallel(cx, &pool, &lease, usize::MAX, 8, &|_, x| {
                calls.fetch_add(1, Ordering::Relaxed);
                Ok::<_, Infallible>(objective(x))
            }).unwrap();
            assert_eq!(calls.load(Ordering::Relaxed), limit);
            assert_eq!(report.evaluations, limit);
            assert_eq!(report.stop, DifferentialEvolutionStop::EvaluationLimit);
        }
        let mut search = state(cx, 64);
        for width in [1, 5, 2, 7, 11] {
            let before = search.evaluations();
            let report = search.advance_parallel(cx, &pool, &lease, width, 8,
                &|_, x| Ok::<_, Infallible>(objective(x))).unwrap();
            assert_eq!(report.evaluations-before, width);
            assert_eq!(report.stop, DifferentialEvolutionStop::Paused);
        }
    });
}

#[test]
fn g4_callback_failure_keeps_successful_siblings_and_serial_retry_uses_one_attempt() {
    with_cx(|_, cx| {
        let pool = TilePool::for_host(4, 17);
        let lease = OperationMemoryLease::bounded(64 << 20);
        let calls = AtomicUsize::new(0);
        let mut search = state(cx, 32);
        let result = search.advance_parallel(cx, &pool, &lease, 8, 8, &|_, x| {
            calls.fetch_add(1, Ordering::Relaxed);
            if x == [0.0; 3] { Err("physical refusal") } else { Ok(objective(x)) }
        });
        assert!(matches!(result, Err(DifferentialEvolutionParallelError::Evaluation {
            member: 0, source: "physical refusal" })));
        assert_eq!(search.evaluations(), 8);
        assert_eq!(calls.load(Ordering::Relaxed), 8);
        assert!(search.population().is_empty());
        let report = search.advance(cx, 1, &mut |x| {
            assert_eq!(x, [0.0; 3]);
            calls.fetch_add(1, Ordering::Relaxed);
            Ok::<_, Infallible>(objective(x))
        }).unwrap();
        assert_eq!(report.initialized_members, 8);
        assert_eq!(report.evaluations, 9);
        assert_eq!(calls.load(Ordering::Relaxed), 9);
        let mut reference = state(cx, 8);
        reference.advance(cx, 8, &mut |x| Ok::<_, Infallible>(objective(x))).unwrap();
        assert_eq!(search.population(), reference.population());
        assert_eq!(search.best(), reference.best());
    });
}

#[test]
fn g4_cancelled_wave_caches_completed_physics_and_resumes_on_a_fresh_context() {
    let calls = AtomicUsize::new(0);
    let mut search = with_cx(|gate, cx| {
        let pool = TilePool::for_host(1, 17);
        let lease = OperationMemoryLease::bounded(64 << 20);
        let mut search = state(cx, 32);
        let result = search.advance_parallel(cx, &pool, &lease, 8, 8, &|_, x| {
            calls.fetch_add(1, Ordering::Relaxed);
            gate.request();
            Ok::<_, Infallible>(objective(x))
        });
        assert!(matches!(result, Err(DifferentialEvolutionParallelError::Solver(
            DifferentialEvolutionError::Cancelled))));
        assert_eq!(search.evaluations(), 1);
        assert!(search.population().is_empty());
        search
    });
    with_cx(|_, cx| {
        let report = search.advance(cx, 0, &mut |_| -> Result<f64, Infallible> {
            panic!("completed callback was repeated")
        }).unwrap();
        assert_eq!(report.initialized_members, 1);
        assert_eq!(calls.load(Ordering::Relaxed), 1);
        let pool = TilePool::for_host(4, 17);
        let lease = OperationMemoryLease::bounded(64 << 20);
        search.advance_parallel(cx, &pool, &lease, 31, 8, &|_, x| {
            calls.fetch_add(1, Ordering::Relaxed);
            Ok::<_, Infallible>(objective(x))
        }).unwrap();
        let mut reference = state(cx, 32);
        reference.advance(cx, 32, &mut |x| Ok::<_, Infallible>(objective(x))).unwrap();
        assert_eq!(search.population(), reference.population());
        assert_eq!(search.best(), reference.best());
        assert_eq!(calls.load(Ordering::Relaxed), 32);
    });
}

#[test]
fn g0_nonfinite_outcome_and_executor_refusal_never_become_fitness() {
    with_cx(|_, cx| {
        let pool = TilePool::for_host(2, 17);
        let lease = OperationMemoryLease::bounded(64 << 20);
        let mut search = state(cx, 16);
        let result = search.advance_parallel(cx, &pool, &lease, 8, 8,
            &|_, x| Ok::<_, Infallible>(if x == [0.0; 3] { f64::NAN } else { objective(x) }));
        assert!(matches!(result, Err(DifferentialEvolutionParallelError::Solver(
            DifferentialEvolutionError::NonFinite(_)))));
        assert!(search.best().is_none());
        assert!(search.population().is_empty());
        assert_eq!(search.evaluations(), 8);
        let report = search.advance_parallel(cx, &pool, &lease, 1, 1,
            &|_, x| Ok::<_, Infallible>(objective(x))).unwrap();
        assert_eq!(report.initialized_members, 8);
        assert_eq!(report.evaluations, 9);
        let mut refused = state(cx, 16);
        let called = AtomicBool::new(false);
        assert!(refused.advance_parallel(cx, &pool, &OperationMemoryLease::bounded(0), 8, 8,
            &|_, _| { called.store(true, Ordering::Relaxed); Ok::<_, Infallible>(0.0) }).is_err());
        assert!(!called.load(Ordering::Relaxed));
        assert_eq!(refused.evaluations(), 0);
        assert!(refused.best().is_none());
    });
}

#[test]
fn g1_target_stops_between_waves_and_invalid_width_spends_nothing() {
    with_cx(|_, cx| {
        let pool = TilePool::for_host(4, 17);
        let lease = OperationMemoryLease::bounded(64 << 20);
        let mut search = state(cx, 32);
        search.config.target = Some(0.0);
        let report = search.advance_parallel(cx, &pool, &lease, 32, 8,
            &|_, x| Ok::<_, Infallible>(x.iter().map(|x| x*x).sum())).unwrap();
        assert_eq!(report.stop, DifferentialEvolutionStop::TargetReached);
        assert_eq!(report.evaluations, 8);
        assert_eq!(report.best.unwrap().value, 0.0);
        let before = search.evaluations();
        assert!(search.advance_parallel(cx, &pool, &lease, 32, 0,
            &|_, _| -> Result<f64, Infallible> { panic!("invalid wave invoked physics") }).is_err());
        assert_eq!(search.evaluations(), before);
    });
}
