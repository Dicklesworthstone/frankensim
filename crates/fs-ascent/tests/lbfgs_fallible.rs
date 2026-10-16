//! Direct behavior checks for fallible callbacks and hard evaluation ceilings.
use fs_ascent::lbfgs::{LbfgsError, LbfgsState};
use fs_ascent::{StopReason, StopRule};
use fs_ascent::wolfe::try_strong_wolfe_with_budget;
use std::cell::Cell;

fn rosenbrock(x: &[f64]) -> (f64, Vec<f64>) {
    let a = x[1] - x[0] * x[0];
    (100.0 * a * a + (1.0 - x[0]).powi(2),
     vec![-400.0 * x[0] * a - 2.0 * (1.0 - x[0]), 200.0 * a])
}

fn quadratic(x: &[f64]) -> Result<(f64, Vec<f64>), &'static str> {
    Ok((x[0] * x[0], vec![2.0 * x[0]]))
}

fn assert_same(a: &LbfgsState, b: &LbfgsState) {
    // Includes private curvature pairs, not just the final point.
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    assert_eq!(a.f.to_bits(), b.f.to_bits());
    for (x, y) in a.x.iter().zip(&b.x) { assert_eq!(x.to_bits(), y.to_bits()); }
}

#[test]
fn legacy_and_fallible_paths_have_identical_probes_and_checkpoints() {
    let rule = StopRule::GradNorm(0.0);
    let mut old_trace = Vec::new();
    let mut legacy_fg = |x: &[f64]| { old_trace.push(x.to_vec()); rosenbrock(x) };
    let mut old = LbfgsState::new(&[-1.2, 1.0], 7, &mut legacy_fg);
    old.run(&mut legacy_fg, &rule, 12);
    let mut new_trace = Vec::new();
    let mut fg = |x: &[f64]| -> Result<_, &'static str> {
        new_trace.push(x.to_vec()); Ok(rosenbrock(x))
    };
    let mut new = LbfgsState::try_new(&[-1.2, 1.0], 7, &mut fg).unwrap();
    new.try_run(&mut fg, &rule, 12, usize::MAX).unwrap();
    assert_eq!(old_trace, new_trace);
    assert_same(&old, &new);
}

#[test]
fn accepted_iteration_splits_retain_curvature_and_accounting() {
    let mut fg = |x: &[f64]| Ok::<_, &'static str>(rosenbrock(x));
    let initial = LbfgsState::try_new(&[-1.2, 1.0], 7, &mut fg).unwrap();
    let rule = StopRule::GradNorm(0.0);
    let mut straight = initial.clone();
    straight.try_run(&mut fg, &rule, 12, usize::MAX).unwrap();
    for split in 0..=12 {
        let mut resumed = initial.clone();
        resumed.try_run(&mut fg, &rule, split, usize::MAX).unwrap();
        resumed.try_run(&mut fg, &rule, 12 - split, usize::MAX).unwrap();
        assert_same(&straight, &resumed);
    }
}

#[test]
fn budget_stops_during_zoom_without_spending_or_accepting_an_extra_trial() {
    let calls = Cell::new(0);
    let mut fg = |x: &[f64]| { calls.set(calls.get() + 1); quadratic(x) };
    let mut state = LbfgsState::try_new(&[1.0], 4, &mut fg).unwrap();
    let report = state.try_run(&mut fg, &StopRule::GradNorm(0.0), 10, 2).unwrap();
    assert_eq!(report.reason, StopReason::Budget);
    assert_eq!(calls.get(), 2);
    assert_eq!(state.evals, 2);
    assert_eq!(state.x, vec![1.0]);
    assert_eq!(state.iters, 0);
    assert_eq!(state.history, vec![1.0]);
    let saved = state.clone();
    state.try_run(&mut fg, &StopRule::GradNorm(0.0), 10, 2).unwrap();
    assert_same(&saved, &state);
    assert_eq!(calls.get(), 2);
}

#[test]
fn final_budgeted_probe_can_land_a_valid_step() {
    let mut state = LbfgsState::try_new(&[1.0], 4, &mut quadratic).unwrap();
    let report = state.try_run(&mut quadratic, &StopRule::GradNorm(0.0), 10, 3).unwrap();
    assert_eq!(report.reason, StopReason::Budget);
    assert_eq!(report.evals, 3);
    assert_eq!(state.x, vec![0.0]);
    assert_eq!(state.history, vec![1.0, 0.0]);
    assert_eq!(state.iters, 1);
}

#[test]
fn all_combinator_cannot_override_a_hard_budget_leaf() {
    let mut state = LbfgsState::try_new(&[1.0], 4, &mut quadratic).unwrap();
    let rule = StopRule::All(vec![
        StopRule::ObjectiveBelow(-1.0),
        StopRule::Any(vec![StopRule::Budget(1), StopRule::Budget(5)]),
    ]);
    let report = state.try_run(
        &mut |_| -> Result<_, &'static str> { panic!("exhausted budget evaluated objective") },
        &rule, 10, 100,
    ).unwrap();
    assert_eq!(report.reason, StopReason::Budget);
    assert_eq!(report.evals, 1);
}

#[test]
fn failed_zoom_callback_keeps_checkpoint_and_counts_attempts() {
    let calls = Cell::new(0usize);
    let mut fg = |x: &[f64]| -> Result<_, &'static str> {
        calls.set(calls.get() + 1);
        if calls.get() == 3 { return Err("cancelled inside zoom"); }
        Ok((50.0 * x[0] * x[0], vec![100.0 * x[0]]))
    };
    let mut state = LbfgsState::try_new(&[1.0], 4, &mut fg).unwrap();
    let initial = state.clone();
    let error = state.try_run(&mut fg, &StopRule::GradNorm(1e-10), 20, 100).unwrap_err();
    assert_eq!(error, LbfgsError::Evaluation("cancelled inside zoom"));
    assert_eq!(state.evals, 3);
    assert_eq!(calls.get(), 3);
    let mut numerical = state.clone(); numerical.evals = initial.evals;
    assert_same(&numerical, &initial);
    let report = state.try_run(&mut fg, &StopRule::GradNorm(1e-10), 20, 100).unwrap();
    assert_eq!(report.reason, StopReason::GradNorm);
    assert!(state.x[0].abs() < 1e-10);
    assert_eq!(state.evals, calls.get());
}

#[test]
fn malformed_gradients_return_errors_not_truncated_zip_steps() {
    let mut state = LbfgsState::try_new(&[1.0], 4, &mut quadratic).unwrap();
    let error = state.try_run(
        &mut |_| Ok::<_, &'static str>((0.0, vec![])),
        &StopRule::GradNorm(0.0), 1, 10,
    ).unwrap_err();
    assert_eq!(error, LbfgsError::GradientLength { expected: 1, actual: 0 });
    assert_eq!(state.x, vec![1.0]);
    assert_eq!(state.evals, 2);
    assert!(matches!(LbfgsState::try_new(&[1.0], 0, &mut quadratic), Err(LbfgsError::InvalidInput(_))));
}

#[test]
fn fallible_wolfe_propagates_errors_in_zoom_and_respects_zero_budget() {
    let mut calls = 0;
    let mut phi = |_alpha| -> Result<_, &'static str> {
        calls += 1;
        if calls == 2 { Err("domain service unavailable") } else { Ok((5.0, -1.0)) }
    };
    assert_eq!(try_strong_wolfe_with_budget(&mut phi, 1.0, -1.0, 1.0, 1e-4, 0.9, 10).unwrap_err(), "domain service unavailable");
    assert_eq!(calls, 2);
    let outcome = try_strong_wolfe_with_budget(
        &mut |_| -> Result<_, &'static str> { panic!("zero budget called phi") },
        1.0, -1.0, 1.0, 1e-4, 0.9, 0,
    ).unwrap();
    assert!(!outcome.success);
    assert_eq!(outcome.evals, 0);
}
