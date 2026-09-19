//! Product regressions for scale-aware constrained SQP, not golden receipts.
use fs_ascent::auglag::ConstrainedProblem;
use fs_ascent::{SqpReport, sqp};
use std::cell::Cell;

fn costly_feasibility(inequality: bool, scale: f64) -> SqpReport {
    let mut fg = |x: &[f64]| {
        (1e6 * x[0] + 0.5 * x[0] * x[0], vec![1e6 + x[0]])
    };
    let ce = |x: &[f64]| {
        if inequality { vec![] } else { vec![scale * (x[0] - 2.0)] }
    };
    let ce_jt = |_: &[f64], w: &[f64]| {
        vec![if inequality { 0.0 } else { scale * w[0] }]
    };
    let ci = |x: &[f64]| {
        if inequality { vec![scale * (2.0 - x[0])] } else { vec![] }
    };
    let ci_jt = |_: &[f64], w: &[f64]| {
        vec![if inequality { -scale * w[0] } else { 0.0 }]
    };
    let mut problem = ConstrainedProblem {
        fg: &mut fg, ce: &ce, ce_jt: &ce_jt, ci: &ci, ci_jt: &ci_jt,
    };
    sqp(&mut problem, &[1.0], 1e-7, 20)
}

#[test]
fn equality_can_restore_feasibility_against_a_large_objective_cost() {
    let r = costly_feasibility(false, 1.0);
    assert!(r.converged, "{r:?}");
    assert!((r.x[0] - 2.0).abs() < 1e-10);
    assert!((r.lambda[0] + 1_000_002.0).abs() < 1e-7);
    assert!(r.kkt.within_tolerance(1e-7));
}

#[test]
fn active_inequality_can_restore_feasibility_against_a_large_objective_cost() {
    let r = costly_feasibility(true, 1.0);
    assert!(r.converged, "{r:?}");
    assert!((r.x[0] - 2.0).abs() < 1e-10);
    assert!((r.nu[0] - 1_000_002.0).abs() < 1e-7);
    assert!(r.kkt.within_tolerance(1e-7));
}

#[test]
fn constraint_rescaling_does_not_disable_feasibility_restoration() {
    for scale in [1e-3, 1.0, 1e3] {
        for inequality in [false, true] {
            let r = costly_feasibility(inequality, scale);
            assert!(r.converged, "scale={scale}, inequality={inequality}: {r:?}");
            assert!((r.x[0] - 2.0).abs() < 1e-9);
            let physical_multiplier = if inequality { r.nu[0] } else { -r.lambda[0] } * scale;
            assert!((physical_multiplier - 1_000_002.0).abs() < 1e-7);
        }
    }
}

#[test]
fn small_objective_improvements_are_not_blocked_by_an_absolute_merit_floor() {
    let mut fg = |x: &[f64]| (0.5 * x[0] * x[0], vec![x[0]]);
    let empty = |_: &[f64]| vec![];
    let zero = |_: &[f64], _: &[f64]| vec![0.0];
    let mut problem = ConstrainedProblem {
        fg: &mut fg, ce: &empty, ce_jt: &zero, ci: &empty, ci_jt: &zero,
    };
    let r = sqp(&mut problem, &[1e-7], 1e-12, 10);
    assert!(r.converged, "{r:?}");
    assert!(r.x[0].abs() < 1e-12);
}

#[test]
fn objective_accounting_includes_validation_trials_and_terminal_certificates() {
    for max_iters in [0, 1, 10] {
        let calls = Cell::new(0usize);
        let mut fg = |x: &[f64]| {
            calls.set(calls.get() + 1);
            (0.5 * x[0] * x[0], vec![x[0]])
        };
        let empty = |_: &[f64]| vec![];
        let zero = |_: &[f64], _: &[f64]| vec![0.0];
        let mut problem = ConstrainedProblem {
            fg: &mut fg, ce: &empty, ce_jt: &zero, ci: &empty, ci_jt: &zero,
        };
        let r = sqp(&mut problem, &[1.0], 1e-12, max_iters);
        assert_eq!(r.evals, calls.get(), "max_iters={max_iters}");
    }
}
