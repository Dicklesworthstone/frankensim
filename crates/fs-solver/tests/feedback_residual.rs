//! Coupled stored-system error checked against independent small dense solves.
use fs_solver::goal::{GoalResidualError, GoalResidualLimits};
use fs_solver::goal::feedback::{
    FeedbackBoundStatus, FeedbackResidualLimits, enclose_affine_feedback_error,
};
use fs_sparse::{Coo, Csr};

fn matrix(rows: &[&[f64]]) -> Csr {
    let mut coo = Coo::new(rows.len(), rows[0].len());
    for (i, row) in rows.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() { if value != 0.0 { coo.push(i, j, value); } }
    }
    coo.assemble()
}
fn limits() -> FeedbackResidualLimits {
    FeedbackResidualLimits {
        solid: GoalResidualLimits { max_rows: 32, max_nonzeros: 1024 },
        max_ports: 8, max_transfer_nonzeros: 1024,
        max_response_entries: 256, max_verification_entries: 10_000,
    }
}

#[test]
fn freezing_a_reference_would_miss_the_whole_error() {
    let a = matrix(&[&[2.0]]);
    let b = matrix(&[&[1.0]]);
    let c = matrix(&[&[0.5]]);
    // x=0 solves A*x=0 exactly, but not (A-B*C)*x=B*d=1.
    let report = enclose_affine_feedback_error(&a, &[0.0], &[0.0], &b, &c, &[1.0],
        Some(&[vec![0.5]]), None, limits(), || true).unwrap();
    assert_eq!(report.status(), FeedbackBoundStatus::Enclosed);
    assert!(report.residual_infinity_upper() >= 1.0);
    assert!(report.gain_infinity_upper().unwrap() >= 0.25);
    let error = report.state_error_infinity_upper().unwrap();
    assert!(3.0 * error >= 2.0 && error < 0.667);
    assert_eq!(report.response_residual_infinity_upper().len(), 1);
}

#[test]
fn checked_response_columns_unlock_a_bound_the_norm_only_route_cannot() {
    let a = matrix(&[&[2.0, -1.0], &[-1.0, 2.0]]);
    let b = matrix(&[&[1.0], &[0.0]]);
    let c = matrix(&[&[0.0, 1.25]]);
    let run = |responses: Option<&[Vec<f64>]>| enclose_affine_feedback_error(
        &a, &[1.0, 0.0], &[0.0, 0.0], &b, &c, &[0.0], responses, None, limits(), || true,
    ).unwrap();
    let norm_only = run(None);
    assert_eq!(norm_only.status(), FeedbackBoundStatus::ContractionNotEstablished);
    assert!(norm_only.state_error_infinity_upper().is_none());
    let checked = run(Some(&[vec![2.0 / 3.0, 1.0 / 3.0]]));
    assert_eq!(checked.status(), FeedbackBoundStatus::Enclosed);
    assert!(checked.gain_infinity_upper().unwrap() < 0.834);
    assert!(7.0 * checked.state_error_infinity_upper().unwrap() >= 8.0);
    // A missing response is NOT an exact zero inverse column.
    let bad = run(Some(&[vec![0.0; 2]]));
    assert_eq!(bad.status(), FeedbackBoundStatus::ContractionNotEstablished);
    assert!(bad.response_residual_infinity_upper()[0] >= 1.0);
    assert!(bad.state_error_infinity_upper().is_none());
}

#[test]
fn no_contraction_is_unknown_not_a_singularity_claim() {
    let a = matrix(&[&[1.0]]);
    let b = matrix(&[&[1.0]]);
    for coefficient in [1.0, 2.0, -2.0] {
        let c = matrix(&[&[coefficient]]);
        let r = enclose_affine_feedback_error(&a, &[1.0], &[0.0], &b, &c, &[0.0],
            Some(&[vec![1.0]]), None, limits(), || true).unwrap();
        assert_eq!(r.status(), FeedbackBoundStatus::ContractionNotEstablished);
        assert!(r.coupled_inverse_infinity_upper().is_none());
    }
    let singular = matrix(&[&[0.0]]);
    let c = matrix(&[&[0.0]]);
    let r = enclose_affine_feedback_error(&singular, &[1.0], &[0.0], &b, &c, &[0.0],
        None, None, limits(), || true).unwrap();
    assert_eq!(r.status(), FeedbackBoundStatus::SolidInverseUnavailable);
    assert!(r.state_error_infinity_upper().is_none());
}

// Pivoted elimination is deliberately disjoint from production CG and its reports.
fn exact_small(a: &Csr, b: &Csr, c: &Csr, rhs: &[f64], d: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let mut rows = vec![vec![0.0; n + 1]; n];
    for i in 0..n {
        rows[i][n] = rhs[i];
        for k in 0..d.len() { rows[i][n] += b.get(i, k) * d[k]; }
        for j in 0..n {
            rows[i][j] = a.get(i, j);
            for k in 0..d.len() { rows[i][j] -= b.get(i, k) * c.get(k, j); }
        }
    }
    for k in 0..n {
        let pivot = (k..n).max_by(|&i, &j| rows[i][k].abs().total_cmp(&rows[j][k].abs())).unwrap();
        rows.swap(k, pivot);
        let diagonal = rows[k][k];
        assert!(diagonal.abs() > 1e-12);
        for j in k..=n { rows[k][j] /= diagonal; }
        for i in 0..n { if i != k {
            let value = rows[i][k];
            for j in k..=n { let pivot_value = rows[k][j]; rows[i][j] -= value * pivot_value; }
        } }
    }
    rows.iter().map(|row| row[n]).collect()
}

#[test]
fn signed_nonsymmetric_feedback_and_system_scaling_match_dense_oracles() {
    let b0 = [[0.5, -0.25], [0.0, 0.25]];
    let c = matrix(&[&[0.125, -0.25], &[0.25, 0.125]]);
    let d = [1.5, -0.75];
    let primal = [0.375, -0.125];
    for scale in [0.0625, 1.0, 16.0] {
        let a = matrix(&[&[3.0 * scale, scale], &[-scale, 2.0 * scale]]);
        let b = matrix(&[&b0[0].map(|v| v * scale), &b0[1].map(|v| v * scale)]);
        let rhs = [scale, -scale];
        let solved = exact_small(&a, &b, &c, &rhs, &d);
        let r = enclose_affine_feedback_error(&a, &rhs, &primal, &b, &c, &d,
            None, None, limits(), || true).unwrap();
        assert_eq!(r.status(), FeedbackBoundStatus::Enclosed);
        for i in 0..2 { assert!((solved[i] - primal[i]).abs() <= r.state_error_infinity_upper().unwrap()); }
    }
}

#[test]
fn zero_feedback_is_not_zero_residual_and_tiny_products_do_not_disappear() {
    let a = matrix(&[&[1.0]]);
    let b = matrix(&[&[1.0]]);
    let c = matrix(&[&[0.0]]);
    let tiny = f64::from_bits(1);
    let r = enclose_affine_feedback_error(&a, &[0.0], &[0.0], &b, &c, &[tiny],
        None, None, limits(), || true).unwrap();
    assert_eq!(r.gain_infinity_upper(), Some(0.0));
    assert!(r.residual_infinity_upper() >= tiny);
    assert!(r.state_error_infinity_upper().unwrap() >= tiny);
}

#[test]
fn shape_and_work_limits_refuse_before_repeated_verification() {
    let a = matrix(&[&[2.0]]);
    let b = matrix(&[&[1.0]]);
    let c = matrix(&[&[0.25]]);
    let columns = [vec![0.5]];
    let mut policies = Vec::new();
    let mut l = limits(); l.max_ports = 0; policies.push(l);
    let mut l = limits(); l.max_transfer_nonzeros = 1; policies.push(l);
    let mut l = limits(); l.max_response_entries = 0; policies.push(l);
    let mut l = limits(); l.max_verification_entries = 3; policies.push(l);
    for l in policies {
        let mut calls = 0;
        let r = enclose_affine_feedback_error(&a, &[0.0], &[0.0], &b, &c, &[0.0],
            Some(&columns), None, l, || { calls += 1; true });
        assert!(matches!(r, Err(GoalResidualError::Limit { .. })));
        assert_eq!(calls, 1);
    }
    let r = enclose_affine_feedback_error(&a, &[0.0], &[0.0], &b, &c, &[0.0],
        Some(&[vec![]]), None, limits(), || true);
    assert!(matches!(r, Err(GoalResidualError::Length { field: "response column", .. })));
    let bad = matrix(&[&[f64::NAN]]);
    assert!(matches!(enclose_affine_feedback_error(&a, &[0.0], &[0.0], &b, &bad, &[0.0],
        None, None, limits(), || true), Err(GoalResidualError::NonFinite { .. })));
}

#[test]
fn cancellation_including_the_publication_boundary_is_atomic_and_replayable() {
    let a = matrix(&[&[2.0]]);
    let b = matrix(&[&[1.0]]);
    let c = matrix(&[&[0.5]]);
    let columns = [vec![0.5]];
    let mut count = 0;
    let want = enclose_affine_feedback_error(&a, &[0.0], &[0.0], &b, &c, &[1.0],
        Some(&columns), None, limits(), || { count += 1; true }).unwrap();
    for stop in [1, count / 2, count] {
        let mut calls = 0;
        assert_eq!(enclose_affine_feedback_error(&a, &[0.0], &[0.0], &b, &c, &[1.0],
            Some(&columns), None, limits(), || { calls += 1; calls != stop }),
            Err(GoalResidualError::Cancelled));
    }
    let got = enclose_affine_feedback_error(&a, &[0.0], &[0.0], &b, &c, &[1.0],
        Some(&columns), None, limits(), || true).unwrap();
    assert_eq!(want, got);
}
