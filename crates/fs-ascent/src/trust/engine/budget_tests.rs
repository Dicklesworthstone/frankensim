//! G0/G4/G5: real Krylov, final-model and objective callbacks under hard ceilings.
use super::*;
use std::cell::Cell;

fn fg(x: &[f64]) -> (f64, Vec<f64>) {
    (0.5*x[0]*x[0] + x[1]*x[1], vec![x[0], 2.0*x[1]])
}
fn hv(_: &[f64], v: &[f64]) -> Vec<f64> { vec![v[0], 2.0*v[1]] }
fn incumbent_equal(a: &TrustRegionState, b: &TrustRegionState) {
    assert_eq!(a.x, b.x);
    assert_eq!(a.f.to_bits(), b.f.to_bits());
    assert_eq!(a.g, b.g);
    assert_eq!(a.delta.to_bits(), b.delta.to_bits());
    assert_eq!(a.history, b.history);
    assert_eq!(a.iters, b.iters);
    assert_eq!(a.negative_curvature_hits, b.negative_curvature_hits);
    assert_eq!(a.stalled, b.stalled);
}

#[test]
fn zero_partial_and_exact_hessian_allowances_cover_every_product_role() {
    let initial = TrustRegionState::new(&[0.1, 0.1], &mut fg);
    let rule = StopRule::GradNorm(1e-12);
    let mut straight = initial.clone();
    straight.run(&mut fg, &mut hv, &rule, 1);
    assert_eq!(straight.hv_evals, 3); // two CG directions plus final model
    for cap in 0..=3 {
        let mut state = initial.clone();
        let calls = Cell::new(0);
        let result = state.run_bounded(&mut fg, &mut |x, v| {
            assert!(calls.get() < cap, "callback ran past its ceiling");
            calls.set(calls.get() + 1); hv(x, v)
        }, &rule, 1, cap);
        assert_eq!(calls.get(), cap);
        assert_eq!(result.solution.hv_evals, cap);
        if cap < 3 {
            assert_eq!(result.progress, TrustRegionProgress::HessianBudget);
            incumbent_equal(&initial, &state);
            assert_eq!(state.evals, 1);
        } else {
            incumbent_equal(&straight, &state);
            assert_eq!(state.evals, straight.evals);
            assert_ne!(result.progress, TrustRegionProgress::HessianBudget);
        }
    }
}

#[test]
fn exhausted_retries_spend_nothing_and_extension_restarts_without_refunds() {
    let initial = TrustRegionState::new(&[0.1, 0.1], &mut fg);
    let rule = StopRule::GradNorm(1e-12);
    let mut straight = initial.clone();
    straight.run(&mut fg, &mut hv, &rule, 10);
    for spent in 1..=2 {
        let mut state = initial.clone();
        let result = state.run_bounded(&mut fg, &mut hv, &rule, 10, spent);
        assert_eq!(result.progress, TrustRegionProgress::HessianBudget);
        for cap in [spent, 0] {
            let result = state.run_bounded(
                &mut |_| panic!("exhausted objective"), &mut |_, _| panic!("exhausted Hessian"),
                &rule, 10, cap,
            );
            assert_eq!(result.progress, TrustRegionProgress::HessianBudget);
            assert_eq!(state.hv_evals, spent);
            incumbent_equal(&initial, &state);
        }
        state.run_bounded(&mut fg, &mut hv, &rule, 10, spent + straight.hv_evals);
        incumbent_equal(&straight, &state);
        assert_eq!(state.evals, straight.evals);
        assert_eq!(state.hv_evals, straight.hv_evals + spent);
    }
}

#[test]
fn negative_curvature_requires_its_own_final_model_allowance() {
    let mut sample = |x: &[f64]| (-x[0]*x[0], vec![-2.0*x[0]]);
    let initial = TrustRegionState::new(&[1.0], &mut sample);
    let mut state = initial.clone();
    let rule = StopRule::GradNorm(0.0);
    let mut product = |_: &[f64], v: &[f64]| vec![-2.0*v[0]];
    assert_eq!(state.run_bounded(&mut sample, &mut product, &rule, 1, 1).progress,
        TrustRegionProgress::HessianBudget);
    incumbent_equal(&initial, &state);
    let result = state.run_bounded(&mut sample, &mut product, &rule, 1, 3);
    assert_eq!(result.solution.x, vec![2.0]);
    assert_eq!(result.solution.negative_curvature_hits, 1);
    assert_eq!(result.solution.hv_evals, 3);
    assert_eq!(state.radius(), 2.0);
}

#[test]
fn cancelled_last_permitted_product_is_charged_without_publication() {
    let initial = TrustRegionState::new(&[0.1, 0.1], &mut fg);
    for cap in 1..=3 {
        let mut state = initial.clone();
        let calls = Cell::new(0);
        let result = state.run_with_limits(&mut |_| panic!("cancelled trial"),
            &mut |x, v| { calls.set(calls.get()+1); hv(x, v) },
            &StopRule::GradNorm(1e-12), 1, cap, &mut || calls.get() == cap);
        assert_eq!(result.progress, TrustRegionProgress::Paused);
        assert_eq!(state.hv_evals, cap);
        incumbent_equal(&initial, &state);
    }
}

#[test]
fn exact_limit_can_publish_a_rejected_trial_but_not_its_point() {
    let mut state = TrustRegionState::new(&[0.1, 0.1], &mut fg);
    let result = state.run_bounded(&mut |_| (f64::INFINITY, vec![0.0; 2]), &mut hv,
        &StopRule::GradNorm(0.0), 1, 3);
    assert_eq!(result.progress, TrustRegionProgress::Stopped(StopReason::IterationCap));
    assert_eq!(result.solution.x, vec![0.1, 0.1]);
    assert_eq!(result.solution.evals, 2);
    assert_eq!(result.solution.hv_evals, 3);
    assert_eq!(state.radius(), 0.25);
    assert_eq!(state.history().len(), 2);
}

#[test]
fn objective_limits_and_zero_iteration_segments_do_not_spend_hessian_work() {
    for rule in [StopRule::Budget(1), StopRule::All(vec![StopRule::GradNorm(0.0), StopRule::Budget(1)])] {
        let mut state = TrustRegionState::new(&[0.1, 0.1], &mut fg);
        let result = state.run_bounded(&mut |_| panic!("objective ceiling"),
            &mut |_, _| panic!("objective ceiling must precede Hessian"), &rule, 10, 0);
        assert_eq!(result.progress, TrustRegionProgress::Stopped(StopReason::Budget));
        assert_eq!(result.solution.hv_evals, 0);
    }
    let mut state = TrustRegionState::new(&[0.1, 0.1], &mut fg);
    let result = state.run_bounded(&mut |_| panic!("zero iterations"),
        &mut |_, _| panic!("zero iterations"), &StopRule::GradNorm(0.0), 0, 0);
    assert_eq!(result.progress, TrustRegionProgress::Stopped(StopReason::IterationCap));
    assert_eq!(result.solution.hv_evals, 0);
}

#[test]
fn unwinding_and_malformed_callbacks_retain_incumbent_and_charge_every_role() {
    use std::panic::{catch_unwind, AssertUnwindSafe};
    let initial = TrustRegionState::new(&[0.1, 0.1], &mut fg);
    let rule = StopRule::GradNorm(1e-12);
    let mut direct = initial.clone();
    direct.run(&mut fg, &mut hv, &rule, 10);
    // Two Krylov calls, the final-model HVP, then the objective/gradient.
    for malformed in [false, true] {
        for fail_at in 1..=4 {
            let mut state = initial.clone();
            let calls = Cell::new(0);
            let result = catch_unwind(AssertUnwindSafe(|| {
                state.run_bounded(&mut |x| {
                    calls.set(calls.get() + 1);
                    if calls.get() == fail_at {
                        if malformed { return (0.0, Vec::new()); }
                        panic!("objective callback fault");
                    }
                    fg(x)
                }, &mut |x, v| {
                    calls.set(calls.get() + 1);
                    if calls.get() == fail_at {
                        if malformed { return Vec::new(); }
                        panic!("Hessian callback fault");
                    }
                    hv(x, v)
                }, &StopRule::Budget(2), 10, 3);
            }));
            assert!(result.is_err());
            assert_eq!(calls.get(), fail_at);
            incumbent_equal(&initial, &state);
            let spent_hv = fail_at.min(3);
            let spent_fg = usize::from(fail_at == 4);
            assert_eq!(state.hv_evals, spent_hv);
            assert_eq!(state.evals, 1 + spent_fg);
            let result = state.run_bounded(&mut fg, &mut hv, &rule, 10,
                spent_hv + direct.hv_evals);
            assert_eq!(result.progress, TrustRegionProgress::Stopped(StopReason::GradNorm));
            incumbent_equal(&direct, &state);
            assert_eq!(state.hv_evals, direct.hv_evals + spent_hv);
            assert_eq!(state.evals, direct.evals + spent_fg);
        }
    }
}

#[test]
fn exhausted_counter_ranges_refuse_before_callbacks_in_both_profiles() {
    let mut state = TrustRegionState::new(&[0.1, 0.1], &mut fg);
    state.evals = usize::MAX;
    let result = state.run(&mut |_| panic!("evaluation counter overflow"),
        &mut |_, _| panic!("evaluation overflow reached Hessian"),
        &StopRule::GradNorm(0.0), 10);
    assert_eq!(result.progress, TrustRegionProgress::Stopped(StopReason::Budget));
    state.evals = 1;
    state.hv_evals = usize::MAX;
    let result = state.run(&mut |_| panic!("Hessian limit reached objective"),
        &mut |_, _| panic!("Hessian counter overflow"),
        &StopRule::GradNorm(0.0), 10);
    assert_eq!(result.progress, TrustRegionProgress::HessianBudget);
    assert_eq!(result.solution.hv_evals, usize::MAX);
    assert_eq!(result.solution.iters, 0);
}
