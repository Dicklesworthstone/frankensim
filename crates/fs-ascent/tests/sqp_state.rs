//! Actual callback, checkpoint and constrained-solve regressions.
use fs_ascent::sqp::{SqpError, SqpSample, SqpState, SqpStop};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use std::cell::Cell;

type SampleResult = Result<Option<SqpSample>, &'static str>;

fn constrained(x: &[f64]) -> SampleResult {
    Ok(Some(SqpSample {
        f: (x[0] - 2.0).powi(2) + (x[1] - 1.0).powi(2),
        gradient: vec![2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)],
        ce: vec![x[0] + x[1] - 2.0], ci: vec![x[0] - 1.2],
        je: vec![1.0, 1.0], ji: vec![1.0, 0.0],
    }))
}

fn parabola(x: &[f64]) -> SampleResult {
    Ok(Some(SqpSample {
        f: x[0] * x[0], gradient: vec![2.0 * x[0]],
        ce: vec![], ci: vec![], je: vec![], ji: vec![],
    }))
}

fn rosenbrock(x: &[f64]) -> SampleResult {
    let r = x[1] - x[0] * x[0];
    Ok(Some(SqpSample {
        f: (1.0 - x[0]).powi(2) + 100.0 * r * r,
        gradient: vec![2.0 * (x[0] - 1.0) - 400.0 * x[0] * r, 200.0 * r],
        ce: vec![], ci: vec![], je: vec![], ji: vec![],
    }))
}

fn bits(v: &[f64]) -> Vec<u64> { v.iter().map(|v| v.to_bits()).collect() }

#[test]
fn constrained_optimum_uses_cached_derivatives_for_its_certificate() {
    let calls = Cell::new(0);
    let mut oracle = |x: &[f64]| { calls.set(calls.get() + 1); constrained(x) };
    let mut state = SqpState::try_new(&[0.0, 0.0], 16, &mut oracle, None).unwrap();
    let report = state.try_run(&mut oracle, 1e-8, 20, 30, None).unwrap();
    assert_eq!(report.stop, SqpStop::Converged);
    assert!(report.solution.kkt.within_tolerance(1e-8));
    assert!((state.point()[0] - 1.2).abs() < 1e-8);
    assert!((state.point()[1] - 0.8).abs() < 1e-8);
    assert!((report.solution.lambda[0] - 0.4).abs() < 1e-8);
    assert!((report.solution.nu[0] - 1.2).abs() < 1e-8);
    assert_eq!(calls.get(), state.evaluations());
    assert_eq!(state.evaluations(), 2);
}

#[test]
fn budget_in_backtracking_retains_the_accepted_point_and_never_overshoots() {
    let mut state = SqpState::try_new(&[2.0], 8, &mut parabola, None).unwrap();
    let report = state.try_run(&mut parabola, 1e-10, 20, 2, None).unwrap();
    assert_eq!(report.stop, SqpStop::EvaluationLimit);
    assert_eq!(state.point(), &[2.0]);
    assert_eq!(state.evaluations(), 2);
    assert_eq!(state.iterations(), 0);
    let report = state.try_run(&mut |_| panic!("exhausted budget evaluated"), 1e-10, 20, 2, None)
        .unwrap_or_else(|_: SqpError<&str>| panic!("unexpected error"));
    assert_eq!(report.stop, SqpStop::EvaluationLimit);
    let report = state.try_run(&mut parabola, 1e-10, 20, 4, None).unwrap();
    assert_eq!(report.stop, SqpStop::EvaluationLimit);
    assert!(report.solution.converged);
    assert_eq!(state.point(), &[0.0]);
    assert_eq!(state.evaluations(), 4);
}

#[test]
fn every_complete_step_split_replays_the_same_state_and_accounting() {
    let initial = SqpState::try_new(&[-1.2, 1.0], 8, &mut rosenbrock, None).unwrap();
    let mut straight = initial.clone();
    straight.try_run(&mut rosenbrock, 1e-12, 12, 1000, None).unwrap();
    assert!(straight.iterations() > 4);
    for split in 0..=12 {
        let mut resumed = initial.clone();
        resumed.try_run(&mut rosenbrock, 1e-12, split, 1000, None).unwrap();
        resumed.try_run(&mut rosenbrock, 1e-12, 12 - split, 1000, None).unwrap();
        assert_eq!(bits(straight.point()), bits(resumed.point()), "split {split}");
        assert_eq!(bits(straight.history()), bits(resumed.history()));
        assert_eq!(straight.sample(), resumed.sample());
        assert_eq!(straight.evaluations(), resumed.evaluations());
        assert_eq!(straight.iterations(), resumed.iterations());
        // Continuation also exercises retained curvature and penalty state.
        let mut a = straight.clone();
        a.try_run(&mut rosenbrock, 1e-12, 3, 1000, None).unwrap();
        resumed.try_run(&mut rosenbrock, 1e-12, 3, 1000, None).unwrap();
        assert_eq!(bits(a.point()), bits(resumed.point()));
        assert_eq!(a.evaluations(), resumed.evaluations());
    }
}

#[test]
fn original_callback_error_is_counted_and_recoverable() {
    let mut state = SqpState::try_new(&[2.0], 8, &mut parabola, None).unwrap();
    let sample = state.sample().clone();
    let error = state.try_run(&mut |_| Err("offline"), 1e-10, 20, 100, None).unwrap_err();
    assert_eq!(error, SqpError::Evaluation("offline"));
    assert_eq!(state.point(), &[2.0]);
    assert_eq!(state.sample(), &sample);
    assert_eq!(state.evaluations(), 2);
    assert_eq!(state.history(), &[4.0]);
    let report = state.try_run(&mut parabola, 1e-10, 20, 100, None).unwrap();
    assert_eq!(report.stop, SqpStop::Converged);
    assert_eq!(state.point(), &[0.0]);
    assert_eq!(state.evaluations(), 4);
}

#[test]
fn unavailable_domain_trial_shrinks_without_fabricating_a_gradient() {
    let mut state = SqpState::try_new(&[2.0], 8, &mut parabola, None).unwrap();
    let mut oracle = |x: &[f64]| if x[0] < -0.5 { Ok(None) } else { parabola(x) };
    let report = state.try_run(&mut oracle, 1e-10, 10, 10, None).unwrap();
    assert_eq!(report.stop, SqpStop::Converged);
    assert_eq!(state.point(), &[0.0]);
    assert_eq!(state.rejected_trials(), 1);
    assert_eq!(state.evaluations(), 3);
}

#[test]
fn malformed_trial_never_replaces_the_accepted_sample() {
    let mut state = SqpState::try_new(&[2.0], 8, &mut parabola, None).unwrap();
    let mut malformed = |x: &[f64]| {
        let mut sample = parabola(x)?.unwrap();
        sample.ci.push(-1.0);
        Ok::<_, &'static str>(Some(sample))
    };
    let error = state.try_run(&mut malformed, 1e-10, 10, 10, None).unwrap_err();
    assert!(matches!(error, SqpError::Shape { field: "inequalities", .. }));
    assert_eq!(state.point(), &[2.0]);
    assert_eq!(state.sample().ci.len(), 0);
    assert_eq!(state.evaluations(), 2);
}

#[test]
fn incompatible_linearization_is_a_stall_not_an_infeasibility_proof() {
    let mut impossible = |x: &[f64]| {
        let mut sample = parabola(x)?.unwrap();
        sample.ce = vec![x[0] - 1.0, x[0] - 2.0];
        sample.je = vec![1.0, 1.0];
        Ok::<_, &'static str>(Some(sample))
    };
    let mut state = SqpState::try_new(&[0.0], 8, &mut impossible, None).unwrap();
    let report = state.try_run(&mut impossible, 1e-10, 10, 10, None).unwrap();
    assert_eq!(report.stop, SqpStop::Stalled);
    assert!(!report.solution.converged);
    assert_eq!(state.evaluations(), 1);
}

#[test]
fn cancellation_after_a_trial_keeps_an_accepted_checkpoint() {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 0, kernel_id: 1, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let mut state = SqpState::try_new(&[0.0, 0.0], 16, &mut constrained, Some(&cx)).unwrap();
        let mut cancelled = |x: &[f64]| { gate.request(); constrained(x) };
        let error = state.try_run(&mut cancelled, 1e-8, 20, 100, Some(&cx)).unwrap_err();
        assert_eq!(error, SqpError::Cancelled);
        assert_eq!(state.point(), &[0.0, 0.0]);
        assert_eq!(state.iterations(), 0);
        assert_eq!(state.evaluations(), 2);
        let report = state.try_run(&mut constrained, 1e-8, 20, 100, None).unwrap();
        assert_eq!(report.stop, SqpStop::Converged);
        assert!(report.solution.kkt.within_tolerance(1e-8));
        assert_eq!(state.evaluations(), 3);
    });
}

#[test]
fn invalid_start_and_dense_caps_refuse_before_solving() {
    let error = SqpState::try_new(&[f64::NAN], 8, &mut |_| panic!("invalid start evaluated"), None)
        .unwrap_err();
    assert!(matches!(error, SqpError::<&str>::Invalid(_)));
    let error = SqpState::try_new(&[0.0, 0.0], 3, &mut constrained, None).unwrap_err();
    assert!(matches!(error, SqpError::Invalid(_)));
}
