//! Exact dyadic checks and an independent two-by-two elimination oracle for
//! the discrete goal certificate. No residual implementation is its oracle.

use fs_solver::goal::{
    GoalBoundStatus, GoalResidualError, GoalResidualLimits, GoalResidualReport,
    ScalarEnclosure, enclose_goal_error,
};
use fs_sparse::{Coo, Csr};

fn matrix(a: [[f64; 2]; 2]) -> Csr {
    let mut coo = Coo::new(2, 2);
    for (i, row) in a.into_iter().enumerate() {
        for (j, value) in row.into_iter().enumerate() { coo.push(i, j, value); }
    }
    coo.assemble()
}

const LIMITS: GoalResidualLimits = GoalResidualLimits { max_rows: 4096, max_nonzeros: 65536 };

fn contains(bound: ScalarEnclosure, value: f64) {
    assert!(bound.lower() <= value && value <= bound.upper(), "{value} outside {bound:?}");
}

fn evaluate(a: &Csr, b: &[f64], x: &[f64], g: &[f64], z: &[f64]) -> GoalResidualReport {
    enclose_goal_error(a, b, x, g, z, None, LIMITS, || true).expect("admitted system")
}

#[test]
fn nonsymmetric_transpose_sign_and_true_goal_error() {
    let a = matrix([[4.0, -1.0], [2.0, 3.0]]);
    // A^T [1/2, -1/2] = [1, -2]. The EXACT goal error is -1/4,
    // although neither coordinate of x* = [5/14, 3/7] is dyadic.
    let report = evaluate(&a, &[1.0, 2.0], &[0.25, 0.25], &[1.0, -2.0], &[0.5, -0.5]);
    contains(report.weighted_residual(), -0.25);
    assert!(report.weighted_residual().upper() < 0.0, "a residual sign flip must fail");
    assert!(report.dual_residual_one_upper() < 1e-13, "must apply A^T, not A");
    assert_eq!(report.status(), GoalBoundStatus::Enclosed);
    contains(report.goal_error().unwrap(), -0.25);
    assert!(report.dual_error_upper().unwrap() < 1e-13);
    assert!(report.evaluation_roundoff_upper() < 1e-13);
}

#[test]
fn an_inexact_dual_cannot_hide_goal_error() {
    let a = matrix([[4.0, -1.0], [2.0, 3.0]]);
    let report = evaluate(&a, &[1.0, 2.0], &[0.25, 0.25], &[1.0, -2.0], &[0.0, 0.0]);
    assert_eq!(report.weighted_residual().magnitude_upper(), 0.0);
    assert!(report.dual_error_upper().unwrap() >= 0.25);
    contains(report.goal_error().unwrap(), -0.25);
}

#[test]
fn checked_scaling_proves_a_bound_when_unscaled_dominance_cannot() {
    let a = matrix([[2.0, -2.0], [-1.0, 2.0]]);
    let unscaled = evaluate(&a, &[1.0, 1.0], &[0.0, 0.0], &[1.0, 0.0], &[0.0, 0.0]);
    assert_eq!(unscaled.status(), GoalBoundStatus::InverseBoundUnavailable);
    assert!(unscaled.goal_error().is_none());
    assert!(unscaled.inverse_infinity_upper().is_none());
    // A*[4,3] = [2,2]; the stored M-matrix has ||A^-1||_inf = 2.
    let scaled = enclose_goal_error(&a, &[1.0, 1.0], &[0.0, 0.0], &[1.0, 0.0], &[0.0, 0.0], Some(&[4.0, 3.0]), LIMITS, || true).unwrap();
    assert_eq!(scaled.status(), GoalBoundStatus::Enclosed);
    let inverse = scaled.inverse_infinity_upper().unwrap();
    assert!(inverse >= 2.0 && inverse < 2.0 + 1e-13);
    contains(scaled.goal_error().unwrap(), 2.0);
    // Scaling is independently checked, not an authority-bearing input bit.
    let wrong = enclose_goal_error(&a, &[1.0, 1.0], &[0.0, 0.0], &[1.0, 0.0], &[0.0, 0.0], Some(&[1.0, 10.0]), LIMITS, || true).unwrap();
    assert_eq!(wrong.status(), GoalBoundStatus::InverseBoundUnavailable);
    assert!(wrong.goal_error().is_none());
}

#[test]
fn singular_and_nondominant_systems_keep_only_the_residual_enclosure() {
    for a in [matrix([[1.0, -1.0], [-1.0, 1.0]]), matrix([[1.0, 2.0], [2.0, 1.0]])] {
        let report = evaluate(&a, &[1.0, 2.0], &[0.25, 0.5], &[1.0, 0.0], &[0.5, 0.25]);
        assert_eq!(report.status(), GoalBoundStatus::InverseBoundUnavailable);
        assert!(report.goal_error().is_none());
        assert!(report.dual_error_upper().is_none());
        assert!(report.weighted_residual().lower().is_finite());
    }
}

#[test]
fn exact_dyadic_family_and_independent_dense_goal_oracle() {
    let entries = [[4.0, -0.5], [1.0, 3.0]];
    let a = matrix(entries);
    for i in -16..=16 {
        for j in -8..=8 {
            let b = [1.0, 2.0];
            let x = [f64::from(i) / 8.0, f64::from(j) / 4.0];
            let g = [0.5, -1.0];
            let z = [f64::from(j) / 16.0, f64::from(i) / 16.0];
            let report = evaluate(&a, &b, &x, &g, &z);
            // All operations here are EXACT for this bounded dyadic family.
            let residual = [b[0] - (4.0*x[0] - 0.5*x[1]), b[1] - (x[0] + 3.0*x[1])];
            let weighted = z[0]*residual[0] + z[1]*residual[1];
            contains(report.weighted_residual(), weighted);
            assert!((weighted - report.nominal_correction()).abs() <= report.evaluation_roundoff_upper());
            // Independent 2x2 elimination, not the iterative solver or the
            // enclosure's residual/dominance code.
            let determinant = entries[0][0]*entries[1][1] - entries[0][1]*entries[1][0];
            let exact_x = [(entries[1][1]*b[0] - entries[0][1]*b[1])/determinant, (entries[0][0]*b[1] - entries[1][0]*b[0])/determinant];
            let error = g[0]*(exact_x[0] - x[0]) + g[1]*(exact_x[1] - x[1]);
            contains(report.goal_error().unwrap(), error);
        }
    }
}

#[test]
fn underflow_is_enclosed_and_not_reported_as_exact_zero() {
    let a = Csr::from_parts(1, 1, vec![0, 1], vec![0], vec![1.0]);
    let tiny = f64::from_bits(1);
    let report = evaluate(&a, &[tiny], &[0.0], &[0.5], &[0.5]);
    // The exact correction is 2^-1075, not representable in binary64.
    assert_eq!(report.nominal_correction(), 0.0);
    assert!(report.weighted_residual().upper() >= tiny);
    assert!(report.evaluation_roundoff_upper() >= tiny);
    assert!(report.goal_error().unwrap().upper() >= tiny);
}

#[test]
fn zero_goal_is_exact_and_replay_is_bitwise_stable() {
    let a = matrix([[4.0, -1.0], [2.0, 3.0]]);
    let b = [1.0, 2.0];
    let x = [0.1, 0.7];
    let report = evaluate(&a, &b, &x, &[0.0, 0.0], &[0.0, 0.0]);
    assert_eq!(report.goal_error().unwrap().magnitude_upper(), 0.0);
    assert_eq!(report.evaluation_roundoff_upper(), 0.0);
    let first = evaluate(&a, &b, &x, &[1.0, -2.0], &[0.5, -0.5]);
    let second = evaluate(&a, &b, &x, &[1.0, -2.0], &[0.5, -0.5]);
    assert_eq!(first, second);
    assert_eq!(first.weighted_residual().lower().to_bits(), second.weighted_residual().lower().to_bits());
    assert_eq!(first.weighted_residual().upper().to_bits(), second.weighted_residual().upper().to_bits());
}

#[test]
fn power_of_two_system_rescaling_preserves_the_mathematical_goal() {
    for power in [-500, -100, 0, 100, 500] {
        let factor = 2.0_f64.powi(power);
        let a = matrix([[4.0*factor, -factor], [2.0*factor, 3.0*factor]]);
        let report = evaluate(&a, &[factor, 2.0*factor], &[0.25, 0.25], &[1.0, -2.0], &[0.5/factor, -0.5/factor]);
        contains(report.weighted_residual(), -0.25);
        contains(report.goal_error().unwrap(), -0.25);
    }
}

#[test]
fn shape_limits_finiteness_scaling_and_arithmetic_fail_closed() {
    let a = matrix([[4.0, -1.0], [2.0, 3.0]]);
    let zero = [0.0, 0.0];
    for limits in [GoalResidualLimits { max_rows: 1, max_nonzeros: 4 }, GoalResidualLimits { max_rows: 2, max_nonzeros: 3 }] {
        assert!(matches!(enclose_goal_error(&a, &zero, &zero, &zero, &zero, None, limits, || true), Err(GoalResidualError::Limit { .. })));
    }
    let exact = GoalResidualLimits { max_rows: 2, max_nonzeros: 4 };
    assert!(enclose_goal_error(&a, &zero, &zero, &zero, &zero, None, exact, || true).is_ok());
    assert!(matches!(enclose_goal_error(&a, &[], &zero, &zero, &zero, None, LIMITS, || true), Err(GoalResidualError::Length { field: "rhs", .. })));
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(enclose_goal_error(&a, &[bad, 0.0], &zero, &zero, &zero, None, LIMITS, || true), Err(GoalResidualError::NonFinite { field: "rhs", .. })));
        let bad_a = matrix([[bad, 0.0], [0.0, 1.0]]);
        assert!(matches!(enclose_goal_error(&bad_a, &zero, &zero, &zero, &zero, None, LIMITS, || true), Err(GoalResidualError::NonFinite { field: "matrix row", .. })));
    }
    for scale in [[0.0, 1.0], [-1.0, 1.0]] {
        assert!(matches!(enclose_goal_error(&a, &zero, &zero, &zero, &zero, Some(&scale), LIMITS, || true), Err(GoalResidualError::NonPositiveScaling { .. })));
    }
    let overflow = Csr::from_parts(1, 1, vec![0, 1], vec![0], vec![f64::MAX]);
    assert_eq!(enclose_goal_error(&overflow, &[0.0], &[2.0], &[1.0], &[1.0], None, LIMITS, || true), Err(GoalResidualError::ArithmeticRange));
    let empty = Csr::from_parts(0, 0, vec![0], vec![], vec![]);
    assert!(matches!(enclose_goal_error(&empty, &[], &[], &[], &[], None, LIMITS, || true), Err(GoalResidualError::Shape { .. })));
}

#[test]
fn cancellation_refuses_without_modifying_any_input() {
    let a = matrix([[4.0, -1.0], [2.0, 3.0]]);
    let b = [1.0, 2.0];
    let x = [0.25, 0.25];
    assert_eq!(enclose_goal_error(&a, &b, &x, &b, &x, None, LIMITS, || false), Err(GoalResidualError::Cancelled));
    assert_eq!(x, [0.25, 0.25]);
    // A long row must not postpone checkpoints until its end.
    let n = 1024;
    let mut coo = Coo::new(n, n);
    for j in 0..n { coo.push(0, j, if j == 0 { 2048.0 } else { -1.0 }); }
    for i in 1..n { coo.push(i, i, 1.0); }
    let large = coo.assemble();
    let values = vec![1.0; n];
    let mut calls = 0;
    let result = enclose_goal_error(&large, &values, &values, &values, &values, None, LIMITS, || { calls += 1; calls < 16 });
    assert_eq!(result, Err(GoalResidualError::Cancelled));
    assert_eq!(calls, 16);
    assert!(values.iter().all(|&value| value == 1.0));
}
