//! Independent dyadic manufactured solution and rational inverse checks.

use fs_solver::goal::{
    GoalBoundStatus, GoalResidualError, GoalResidualLimits, enclose_goal_error,
    inverse::enclose_goal_error_with_inverse,
};
use fs_sparse::{Coo, Csr};

const LIMITS: GoalResidualLimits = GoalResidualLimits {
    max_rows: 4096,
    max_nonzeros: 65536,
};
const ZERO: [f64; 3] = [0.0; 3];

fn matrix<const N: usize>(values: [[f64; N]; N]) -> Csr {
    let mut coo = Coo::new(N, N);
    for (i, row) in values.into_iter().enumerate() {
        for (j, a) in row.into_iter().enumerate() {
            coo.push(i, j, a);
        }
    }
    coo.assemble()
}

fn nondominant() -> Csr {
    // SPD eigenvalues 1/4, 1/4, 5/2; comparison matrix is not an M-matrix.
    matrix([[1.0, 0.75, 0.75], [0.75, 1.0, 0.75], [0.75, 0.75, 1.0]])
}

fn inverse() -> Vec<Vec<f64>> {
    // A = (I + 3 J)/4, so A^-1 = 4 I - 6 J/5 exactly.
    vec![
        vec![2.8, -1.2, -1.2],
        vec![-1.2, 2.8, -1.2],
        vec![-1.2, -1.2, 2.8],
    ]
}

#[test]
fn approximate_inverse_encloses_true_goal_and_retains_dual_defect() {
    let a = nondominant();
    let (rhs, primal, goal, dual) = (
        [1.75, 2.0, 1.25],
        [0.25, 0.5, 0.0],
        [1.0, -2.0, 0.5],
        [0.5, -1.0, 0.25],
    );
    let base = enclose_goal_error(&a, &rhs, &primal, &goal, &dual, None, LIMITS, || true).unwrap();
    assert_eq!(base.status(), GoalBoundStatus::InverseBoundUnavailable);
    let report = enclose_goal_error_with_inverse(
        &a,
        &rhs,
        &primal,
        &goal,
        &dual,
        None,
        &inverse(),
        LIMITS,
        || true,
    )
    .unwrap();
    assert_eq!(report.status(), GoalBoundStatus::Enclosed);
    assert_eq!(report.weighted_residual(), base.weighted_residual());
    let bound = report.inverse_infinity_upper().unwrap();
    assert!(bound >= 5.2 && bound < 5.2 + 1e-12);
    assert!(report.dual_error_upper().unwrap() > 0.0);
    // Exact dyadic manufactured solution [1, 2, -1] gives error -11/4.
    let error = report.goal_error().unwrap();
    assert!(error.lower() <= -2.75 && error.upper() >= -2.75);
    assert!(report.weighted_residual().lower() > -2.75);
}

#[test]
fn nonsymmetric_inverse_columns_are_checked_in_the_declared_orientation() {
    let a = matrix([[1.0, 2.0], [3.0, 4.0]]);
    let columns = vec![vec![-2.0, 1.5], vec![1.0, -0.5]];
    let report = enclose_goal_error_with_inverse(
        &a,
        &[5.0, 11.0],
        &[0.0; 2],
        &[1.0, 0.0],
        &[0.0; 2],
        None,
        &columns,
        LIMITS,
        || true,
    )
    .unwrap();
    let norm = report.inverse_infinity_upper().unwrap();
    assert!(norm >= 3.0 && norm < 3.0 + 1e-12);
    let transposed = vec![vec![-2.0, 1.0], vec![1.5, -0.5]];
    let wrong = enclose_goal_error_with_inverse(
        &a,
        &[5.0, 11.0],
        &[0.0; 2],
        &[1.0, 0.0],
        &[0.0; 2],
        None,
        &transposed,
        LIMITS,
        || true,
    )
    .unwrap();
    assert_eq!(wrong.status(), GoalBoundStatus::InverseBoundUnavailable);
}

#[test]
fn bad_and_singular_proposals_never_supply_authority() {
    for (a, columns, expected) in [
        (
            nondominant(),
            vec![vec![0.0; 3]; 3],
            GoalBoundStatus::InverseBoundUnavailable,
        ),
        (
            nondominant(),
            vec![
                vec![1.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0],
                vec![0.0, 0.0, 1.0],
            ],
            GoalBoundStatus::InverseBoundUnavailable,
        ),
        (
            matrix([[1.0; 3]; 3]),
            inverse(),
            GoalBoundStatus::InverseBoundUnavailable,
        ),
        (
            nondominant(),
            vec![vec![f64::MAX; 3]; 3],
            GoalBoundStatus::BoundNotRepresentable,
        ),
    ] {
        let report = enclose_goal_error_with_inverse(
            &a,
            &ZERO,
            &ZERO,
            &ZERO,
            &ZERO,
            None,
            &columns,
            LIMITS,
            || true,
        )
        .unwrap();
        assert_eq!(report.status(), expected);
        assert!(report.inverse_infinity_upper().is_none());
        assert!(report.goal_error().is_none());
        assert!(report.dual_error_upper().is_none());
    }
}

#[test]
fn existing_bound_is_preserved_but_malformed_proposals_are_refused() {
    let a = matrix([[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]);
    let base = enclose_goal_error(&a, &ZERO, &ZERO, &ZERO, &ZERO, None, LIMITS, || true).unwrap();
    let zero_columns = vec![vec![0.0; 3]; 3];
    assert_eq!(
        enclose_goal_error_with_inverse(
            &a,
            &ZERO,
            &ZERO,
            &ZERO,
            &ZERO,
            None,
            &zero_columns,
            LIMITS,
            || true
        )
        .unwrap(),
        base
    );
    for columns in [vec![], vec![vec![]; 3]] {
        assert!(matches!(
            enclose_goal_error_with_inverse(
                &a,
                &ZERO,
                &ZERO,
                &ZERO,
                &ZERO,
                None,
                &columns,
                LIMITS,
                || true
            ),
            Err(GoalResidualError::Length { .. })
        ));
    }
    let mut invalid = inverse();
    invalid[1][2] = f64::NAN;
    assert!(matches!(
        enclose_goal_error_with_inverse(
            &a,
            &ZERO,
            &ZERO,
            &ZERO,
            &ZERO,
            None,
            &invalid,
            LIMITS,
            || true
        ),
        Err(GoalResidualError::NonFinite { .. })
    ));
    let limits = GoalResidualLimits {
        max_rows: 3,
        max_nonzeros: 35,
    };
    assert!(matches!(
        enclose_goal_error_with_inverse(
            &nondominant(),
            &ZERO,
            &ZERO,
            &ZERO,
            &ZERO,
            None,
            &inverse(),
            limits,
            || true
        ),
        Err(GoalResidualError::Limit {
            field: "inverse verification entries",
            ..
        })
    ));
}

#[test]
fn cancellation_is_checked_during_dense_verification_and_before_publication() {
    let mut entries = [[0.75; 16]; 16];
    for (i, row) in entries.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let a = matrix(entries);
    let zero = [0.0; 16];
    let mut base_calls = 0;
    enclose_goal_error(&a, &zero, &zero, &zero, &zero, None, LIMITS, || {
        base_calls += 1;
        true
    })
    .unwrap();
    let mut calls = 0;
    let result = enclose_goal_error_with_inverse(
        &a,
        &zero,
        &zero,
        &zero,
        &zero,
        None,
        &vec![vec![0.0; 16]; 16],
        LIMITS,
        || {
            calls += 1;
            calls < base_calls + 2
        },
    );
    assert_eq!(result, Err(GoalResidualError::Cancelled));
    assert_eq!(calls, base_calls + 2);
    let a = nondominant();
    base_calls = 0;
    enclose_goal_error(&a, &ZERO, &ZERO, &ZERO, &ZERO, None, LIMITS, || {
        base_calls += 1;
        true
    })
    .unwrap();
    calls = 0;
    assert_eq!(
        enclose_goal_error_with_inverse(
            &a,
            &ZERO,
            &ZERO,
            &ZERO,
            &ZERO,
            None,
            &inverse(),
            LIMITS,
            || {
                calls += 1;
                calls < base_calls + 2
            }
        ),
        Err(GoalResidualError::Cancelled)
    );
}
