//! Independent small-system checks for the opt-in port-Schur inverse route.
use fs_solver::goal::{GoalResidualError, GoalResidualLimits};
use fs_solver::goal::feedback::{
    FeedbackBoundStatus, FeedbackInverseMethod, FeedbackResidualLimits,
    enclose_affine_feedback_error, enclose_affine_feedback_error_with_schur,
};
use fs_sparse::{Coo, Csr};

fn matrix(rows: &[&[f64]]) -> Csr {
    let mut coo = Coo::new(rows.len(), rows[0].len());
    for (i, row) in rows.iter().enumerate() {
        for (j, &value) in row.iter().enumerate() {
            if value != 0.0 { coo.push(i, j, value); }
        }
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
fn strong_feedback_can_be_nonsingular_without_whole_state_contraction() {
    let a = matrix(&[&[1.0]]);
    let b = matrix(&[&[1.0]]);
    for value in [-2.0, 2.0] {
        let c = matrix(&[&[value]]);
        let responses = [vec![1.0]];
        let old = enclose_affine_feedback_error(&a, &[1.0], &[0.0], &b, &c, &[0.0],
            Some(&responses), None, limits(), || true).unwrap();
        assert_eq!(old.status(), FeedbackBoundStatus::ContractionNotEstablished);
        let new = enclose_affine_feedback_error_with_schur(&a, &[1.0], &[0.0], &b, &c, &[0.0],
            Some(&responses), None, limits(), || true).unwrap();
        assert_eq!(new.status(), FeedbackBoundStatus::Enclosed);
        assert_eq!(new.inverse_method(), Some(FeedbackInverseMethod::PortSchurDominance));
        assert!(new.gain_infinity_upper().unwrap() >= 2.0);
        let exact = (1.0 - value).recip();
        assert!(new.schur_inverse_infinity_upper().unwrap() >= exact.abs());
        assert!(new.coupled_inverse_infinity_upper().unwrap() >= exact.abs());
        assert!(new.state_error_infinity_upper().unwrap() >= exact.abs());
        assert!(new.state_error_infinity_upper().unwrap() < 3.01);
    }
}

#[test]
fn singular_schur_and_inaccurate_response_columns_cannot_mint_authority() {
    let a = matrix(&[&[1.0]]);
    let b = matrix(&[&[1.0]]);
    let singular = matrix(&[&[1.0]]);
    // A dishonest approximate response could make I-C*V look invertible.
    // Its checked residual must widen the Schur interval back over zero.
    for response in [0.0, 0.5, 1.0, 2.0, -2.0] {
        let new = enclose_affine_feedback_error_with_schur(&a, &[1.0], &[0.0],
            &b, &singular, &[0.0], Some(&[vec![response]]), None, limits(), || true).unwrap();
        assert!(new.coupled_inverse_infinity_upper().is_none());
        assert!(new.schur_inverse_infinity_upper().is_none());
        assert!(new.state_error_infinity_upper().is_none());
        assert_eq!(new.inverse_method(), None);
    }
    let c = matrix(&[&[-2.0]]);
    let inaccurate = enclose_affine_feedback_error_with_schur(&a, &[1.0], &[0.0],
        &b, &c, &[0.0], Some(&[vec![0.0]]), None, limits(), || true).unwrap();
    assert!(inaccurate.state_error_infinity_upper().is_none());
    let missing = enclose_affine_feedback_error_with_schur(&a, &[1.0], &[0.0],
        &b, &c, &[0.0], None, None, limits(), || true).unwrap();
    assert_eq!(missing.inverse_method(), None);
}

#[test]
fn port_schur_uses_only_remaining_verification_work() {
    let a = matrix(&[&[1.0]]);
    let b = matrix(&[&[1.0]]);
    let c = matrix(&[&[-2.0]]);
    // Base: (1 row+1 nonzero)*2 passes = 4. Extra:
    // 2*C.nnz + n*p + p*p + 3*p = 7. Total exact proxy: 11.
    for (budget, admitted) in [(4, false), (10, false), (11, true)] {
        let mut l = limits();
        l.max_verification_entries = budget;
        let new = enclose_affine_feedback_error_with_schur(&a, &[1.0], &[0.0],
            &b, &c, &[0.0], Some(&[vec![1.0]]), None, l, || true).unwrap();
        assert_eq!(new.state_error_infinity_upper().is_some(), admitted);
    }
}

#[test]
fn already_enclosed_results_remain_exactly_the_original_results() {
    let a = matrix(&[&[2.0]]);
    let b = matrix(&[&[1.0]]);
    let c = matrix(&[&[0.5]]);
    let responses = [vec![0.5]];
    let mut l = limits(); l.max_verification_entries = 4;
    let old = enclose_affine_feedback_error(&a, &[0.0], &[0.0], &b, &c, &[1.0],
        Some(&responses), None, l, || true).unwrap();
    let new = enclose_affine_feedback_error_with_schur(&a, &[0.0], &[0.0], &b, &c, &[1.0],
        Some(&responses), None, l, || true).unwrap();
    assert_eq!(new, old);
    assert_eq!(new.inverse_method(), Some(FeedbackInverseMethod::StateContraction));
    assert!(new.schur_inverse_infinity_upper().is_none());
}

// Partial-pivot elimination is independent of the production response solver
// and the interval Schur checker. Also reconstruct the inverse from unit RHSs.
fn solve(mut a: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Vec<f64> {
    let n = rhs.len();
    for k in 0..n {
        let pivot = (k..n).max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs())).unwrap();
        a.swap(k, pivot); rhs.swap(k, pivot);
        let d = a[k][k]; assert!(d != 0.0);
        for j in k..n { a[k][j] /= d; }
        rhs[k] /= d;
        for i in 0..n { if i != k {
            let factor = a[i][k];
            for j in k..n { let v = a[k][j]; a[i][j] -= factor * v; }
            let pivot_rhs = rhs[k];
            rhs[i] -= factor * pivot_rhs;
        } }
    }
    rhs
}

#[test]
fn multiport_nonsymmetric_system_and_inverse_match_independent_elimination() {
    let c = matrix(&[&[-2.0, 0.125, 0.0], &[0.25, -3.0, 0.0]]);
    let responses = [vec![1.0, 0.0, 1.0], vec![0.0, 1.0, 1.0]];
    let primal = [0.5, -0.25, 0.125];
    for exponent in [-200, 0, 200] {
        let scale = 2.0_f64.powi(exponent);
        let a = matrix(&[&[2.0*scale, 0.0, 0.0], &[0.0, 4.0*scale, 0.0], &[0.0, 0.0, scale]]);
        let b = matrix(&[&[2.0*scale, 0.0], &[0.0, 4.0*scale], &[scale, scale]]);
        let rhs = [scale, -scale, 2.0*scale];
        let d = [0.5, -0.25];
        let new = enclose_affine_feedback_error_with_schur(&a, &rhs, &primal, &b, &c, &d,
            Some(&responses), None, limits(), || true).unwrap();
        assert_eq!(new.inverse_method(), Some(FeedbackInverseMethod::PortSchurDominance));
        let mut coupled = vec![vec![0.0; 3]; 3];
        let mut load = rhs.to_vec();
        for i in 0..3 {
            for k in 0..2 { load[i] += b.get(i, k)*d[k]; }
            for j in 0..3 {
                coupled[i][j] = a.get(i, j);
                for k in 0..2 { coupled[i][j] -= b.get(i, k)*c.get(k, j); }
            }
        }
        let exact = solve(coupled.clone(), load);
        for (&reference, &candidate) in exact.iter().zip(&primal) {
            assert!((reference-candidate).abs() <= new.state_error_infinity_upper().unwrap());
        }
        let mut row_norms = [0.0; 3];
        for column in 0..3 {
            let mut rhs = vec![0.0; 3]; rhs[column] = 1.0;
            for (sum, value) in row_norms.iter_mut().zip(solve(coupled.clone(), rhs)) {
                *sum += value.abs();
            }
        }
        assert!(row_norms.into_iter().fold(0.0_f64, f64::max)
            <= new.coupled_inverse_infinity_upper().unwrap());
    }
}

#[test]
fn cancellation_through_schur_publication_refuses_and_replay_is_identical() {
    let a = matrix(&[&[1.0]]);
    let b = matrix(&[&[1.0]]);
    let c = matrix(&[&[-2.0]]);
    let mut count = 0;
    let want = enclose_affine_feedback_error_with_schur(&a, &[1.0], &[0.0], &b, &c, &[0.0],
        Some(&[vec![1.0]]), None, limits(), || { count += 1; true }).unwrap();
    for stop in 1..=count {
        let mut calls = 0;
        assert_eq!(enclose_affine_feedback_error_with_schur(&a, &[1.0], &[0.0], &b, &c, &[0.0],
            Some(&[vec![1.0]]), None, limits(), || { calls += 1; calls != stop }),
            Err(GoalResidualError::Cancelled));
    }
    let got = enclose_affine_feedback_error_with_schur(&a, &[1.0], &[0.0], &b, &c, &[0.0],
        Some(&[vec![1.0]]), None, limits(), || true).unwrap();
    assert_eq!(want, got);
}
