use super::*;
use super::super::super::{CouplingMethod, tests::controls};
use std::cell::Cell;

fn settings() -> AdaptiveSettings {
    AdaptiveSettings { initial_step_s: 0.5, minimum_step_s: 1.0e-6,
        maximum_step_s: 0.5, method_order: 1 }
}
fn evolution() -> AdaptiveEvolution<f64> {
    AdaptiveEvolution::new(1.0, vec![1.0], 0.0, vec![1.0], controls(1, 8), settings()).unwrap()
}
fn decay(old: &f64, interval: StepInterval, _: &[f64]) -> Result<CouplingTrial<f64>, &'static str> {
    let next = old / (1.0 + interval.duration_s());
    Ok(CouplingTrial { state: next, image: vec![next], balance_residuals: vec![] })
}
fn distance(_: &f64, coarse: &f64, fine: &f64, _: StepInterval) -> Result<f64, &'static str> {
    Ok((coarse - fine).abs() / 1.0e-5)
}

#[test]
fn tighter_temporal_tolerance_reduces_error_against_independent_decay() {
    let mut loose = evolution();
    let a = loose.advance(20_000, &mut decay,
        &mut |_, coarse, fine, _| Ok((coarse - fine).abs() / 1e-3), &mut || false).unwrap();
    let mut tight = evolution();
    let b = tight.advance(20_000, &mut decay, &mut distance, &mut || false).unwrap();
    assert!(a.complete && b.complete);
    let exact = fs_math::det::exp(-1.0);
    assert!((tight.state() - exact).abs() < (loose.state() - exact).abs() / 4.0);
    assert!((tight.state() - exact).abs() < 0.002);
    assert!(b.accepted.len() > a.accepted.len());
    assert!(b.rejected > 0);
    assert!(b.accepted.iter().all(|step| step.error_ratio <= 1.0));
    assert_eq!(tight.time_s(), 1.0);
}

#[test]
fn rejection_publishes_neither_coarse_nor_intermediate_half_state() {
    let mut actual = evolution();
    let report = actual.advance(1, &mut decay, &mut distance, &mut || false).unwrap();
    assert_eq!(report.attempts, 1);
    assert_eq!(report.rejected, 1);
    assert!(report.accepted.is_empty() && !report.complete);
    assert!(report.evaluations > 0);
    assert_eq!(*actual.state(), 1.0);
    assert_eq!(actual.interface(), &[1.0]);
    assert_eq!(actual.time_s(), 0.0);
    assert_eq!(actual.next_step_s(), 0.25);
}

#[test]
fn publish_two_half_steps_not_the_unmeasured_richardson_extrapolation() {
    let mut actual = evolution();
    let report = actual.advance(1, &mut decay,
        &mut |_, _, _, _| Ok(0.0), &mut || false).unwrap();
    assert_eq!(report.accepted.len(), 1);
    assert_eq!(actual.state().to_bits(), ((1.0_f64 / 1.25) / 1.25).to_bits());
    assert_eq!(actual.time_s(), 0.5);
}

#[test]
fn nonconvergent_large_intervals_are_retried_without_corrupting_history() {
    let mut policy = controls(1, 4);
    policy.method = CouplingMethod::RelaxedPicard;
    policy.relaxation = 1.0;
    let mut actual = AdaptiveEvolution::new(1.0, vec![1.0], 0.0, vec![0.5], policy, settings()).unwrap();
    let calls = Cell::new(0);
    let report = actual.advance(1000, &mut |old, interval, x| {
        calls.set(calls.get() + 1);
        if interval.duration_s() > 0.125 {
            Ok(CouplingTrial { state: 1.0e6, image: vec![x[0] + 1.0], balance_residuals: vec![] })
        } else { decay(old, interval, x) }
    }, &mut |_, coarse, fine, _| Ok((coarse - fine).abs() / 0.1), &mut || false).unwrap();
    assert!(report.complete && report.rejected >= 2);
    assert_eq!(report.evaluations, calls.get());
    assert!(report.accepted.iter().all(|step| step.interval.duration_s() <= 0.125));
    assert!(*actual.state() < 1.0 && *actual.state() > 0.6);
}

#[test]
fn producer_refusal_is_not_disguised_as_a_retryable_time_error() {
    let mut actual = evolution();
    let before = actual.clone();
    let error = actual.advance(100, &mut |_, _, _| Err("invalid material"),
        &mut distance, &mut || false).unwrap_err();
    assert!(matches!(error.reason, AdaptiveFailure::Coupling(CouplingError {
        reason: CouplingFailure::Operator("invalid material"), ..
    })));
    assert_eq!(error.report.attempts, 1);
    assert_eq!(error.report.evaluations, 1);
    assert_eq!(actual, before);
}

#[test]
fn cancellation_during_second_half_discards_both_trial_advances() {
    let mut actual = evolution();
    let before = actual.clone();
    let stop = Cell::new(false);
    let error = actual.advance(100, &mut |old, interval, x| {
        if interval.start_s > 0.0 { stop.set(true); }
        decay(old, interval, x)
    }, &mut distance, &mut || stop.get()).unwrap_err();
    assert!(matches!(error.reason, AdaptiveFailure::Coupling(CouplingError {
        reason: CouplingFailure::Cancelled, ..
    })));
    assert_eq!(actual, before);
}

#[test]
fn estimator_error_nonfinite_negative_and_cancellation_all_prevent_publication() {
    for value in [f64::NAN, f64::INFINITY, -1.0] {
        let mut actual = evolution();
        let before = actual.clone();
        let error = actual.advance(1, &mut decay,
            &mut |_, _, _, _| Ok(value), &mut || false).unwrap_err();
        assert!(matches!(error.reason, AdaptiveFailure::InvalidDistance(_)));
        assert_eq!(actual, before);
    }
    let mut actual = evolution();
    let before = actual.clone();
    let error = actual.advance(1, &mut decay,
        &mut |_, _, _, _| Err("unavailable goal"), &mut || false).unwrap_err();
    assert!(matches!(error.reason, AdaptiveFailure::Estimator("unavailable goal")));
    assert_eq!(actual, before);
    let stop = Cell::new(false);
    let error = actual.advance(1, &mut decay,
        &mut |_, _, _, _| { stop.set(true); Ok(0.0) }, &mut || stop.get()).unwrap_err();
    assert!(matches!(error.reason, AdaptiveFailure::Cancelled));
    assert_eq!(actual, before);
}

#[test]
fn minimum_step_refuses_instead_of_silently_relaxing_accuracy() {
    let mut cfg = settings(); cfg.minimum_step_s = 0.25;
    let mut actual = AdaptiveEvolution::new(1.0, vec![1.0], 0.0, vec![1.0], controls(1, 8), cfg).unwrap();
    let error = actual.advance(100, &mut decay,
        &mut |_, _, _, _| Ok(2.0), &mut || false).unwrap_err();
    assert!(matches!(error.reason, AdaptiveFailure::AccuracyFloor { error_ratio: 2.0, .. }));
    assert_eq!(error.report.attempts, 2);
    assert_eq!(error.report.rejected, 1);
    assert_eq!(*actual.state(), 1.0);
    assert_eq!(actual.time_s(), 0.0);
}

#[test]
fn attempt_boundary_checkpoint_including_rejection_replays_bitwise() {
    let mut full = evolution();
    let report = full.advance(20_000, &mut decay, &mut distance, &mut || false).unwrap();
    let mut prefix = evolution();
    let first = prefix.advance(1, &mut decay, &mut distance, &mut || false).unwrap();
    assert_eq!(first.rejected, 1);
    let mut resumed = prefix.clone();
    let mut accepted = first.accepted;
    let mut evaluations = first.evaluations;
    let mut attempts = first.attempts;
    let mut rejected = first.rejected;
    for _ in 0..20_000 {
        let part = resumed.advance(3, &mut decay, &mut distance, &mut || false).unwrap();
        accepted.extend(part.accepted);
        evaluations += part.evaluations; attempts += part.attempts; rejected += part.rejected;
        if part.complete { break; }
    }
    assert_eq!(resumed, full);
    assert_eq!(accepted, report.accepted);
    assert_eq!((evaluations, attempts, rejected), (report.evaluations, report.attempts, report.rejected));
}

#[test]
fn forcing_breakpoints_are_hit_exactly_and_never_crossed_by_substeps() {
    let mut actual = AdaptiveEvolution::new(0.0, vec![0.0], 0.0, vec![0.17, 0.41, 0.45],
        controls(1, 8), settings()).unwrap();
    let report = actual.advance(100, &mut |old, interval, _| {
        assert!(![0.17, 0.41].iter().any(|&t| interval.start_s < t && t < interval.end_s));
        let power = if interval.start_s < 0.17 { 2.0 } else if interval.start_s < 0.41 { -1.0 } else { 0.0 };
        let next = old + power * interval.duration_s();
        Ok(CouplingTrial { state: next, image: vec![next], balance_residuals: vec![] })
    }, &mut distance, &mut || false).unwrap();
    assert!(report.complete);
    assert!(report.accepted.iter().any(|row| row.interval.end_s == 0.17));
    assert!(report.accepted.iter().any(|row| row.interval.end_s == 0.41));
    assert_eq!(actual.time_s(), 0.45);
    assert!((*actual.state() - 0.10).abs() < 1.0e-14);
}

#[test]
fn short_final_interval_is_allowed_but_an_unrepresentable_split_is_not() {
    let mut cfg = settings(); cfg.minimum_step_s = 0.25;
    let mut actual = AdaptiveEvolution::new(1.0, vec![1.0], 0.0, vec![0.1], controls(1, 8), cfg).unwrap();
    assert!(actual.advance(1, &mut decay, &mut |_, _, _, _| Ok(0.0), &mut || false).unwrap().complete);
    let start = 9_007_199_254_740_992.0_f64;
    let mut cfg = settings(); cfg.initial_step_s = 2.0; cfg.maximum_step_s = 2.0;
    let mut actual = AdaptiveEvolution::new(1.0, vec![1.0], start, vec![start + 2.0], controls(1, 8), cfg).unwrap();
    let error = actual.advance(1,
        &mut |_, _, _| -> Result<CouplingTrial<f64>, &'static str> { panic!("no representable split") },
        &mut distance, &mut || false).unwrap_err();
    assert!(matches!(error.reason, AdaptiveFailure::TimeResolution(_)));
    assert_eq!(error.report.evaluations, 0);
}

#[test]
fn zero_budget_completed_schedule_and_invalid_policy_do_no_work() {
    let mut actual = evolution();
    let before = actual.clone();
    let report = actual.advance(0,
        &mut |_, _, _| -> Result<CouplingTrial<f64>, &'static str> { panic!("no work") },
        &mut distance, &mut || false).unwrap();
    assert!(!report.complete && report.attempts == 0);
    assert_eq!(actual, before);
    actual.advance(20_000, &mut decay, &mut distance, &mut || false).unwrap();
    assert!(actual.advance(usize::MAX,
        &mut |_, _, _| -> Result<CouplingTrial<f64>, &'static str> { panic!("already done") },
        &mut distance, &mut || false).unwrap().complete);
    for i in 0..4 {
        let mut cfg = settings();
        match i { 0 => cfg.method_order = 0, 1 => cfg.minimum_step_s = 0.0,
            2 => cfg.initial_step_s = f64::NAN, _ => cfg.maximum_step_s = 0.01 }
        assert!(AdaptiveEvolution::new(1.0, vec![1.0], 0.0, vec![1.0], controls(1, 8), cfg).is_err());
    }
    let mut actual = evolution();
    let error = actual.advance(usize::MAX, &mut decay, &mut distance, &mut || false).unwrap_err();
    assert!(matches!(error.reason, AdaptiveFailure::WorkBudgetOverflow));
    assert_eq!(error.report.evaluations, 0);
}
