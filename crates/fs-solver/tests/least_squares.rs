//! The sparse and regularized paths exercise the production LSQR driver.
use fs_solver::{
    RectLinearOp,
    least_squares::{CsrRectOp, LeastSquaresBuildError, RegularizedLeastSquares},
    lsqr::{LsqrConfig, LsqrStop, lsqr},
};
use fs_sparse::{Coo, Csr};

fn csr(rows: usize, cols: usize, values: &[f64]) -> Csr {
    assert_eq!(values.len(), rows * cols);
    let mut coo = Coo::new(rows, cols);
    for (i, &value) in values.iter().enumerate() {
        if value != 0.0 {
            coo.push(i / cols, i % cols, value);
        }
    }
    coo.assemble()
}

#[test]
fn sparse_rectangular_apply_and_true_transpose_drive_an_inconsistent_fit() {
    let matrix = csr(3, 2, &[1.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    let operator = CsrRectOp::new(&matrix);
    let mut transpose = [f64::NAN; 2];
    operator.apply_transpose(&[2.0, 3.0, 4.0], &mut transpose);
    assert!((transpose[0] - 6.0).abs() < 1e-14);
    assert!((transpose[1] - 7.0).abs() < 1e-14);
    let result = lsqr(&operator, &[1.0, 2.0, 2.0], LsqrConfig::default()).unwrap();
    assert!(result.report.converged());
    assert!((result.x[0] - 2.0 / 3.0).abs() < 1e-10);
    assert!((result.x[1] - 5.0 / 3.0).abs() < 1e-10);
}

#[test]
fn residual_weights_are_squared_by_the_objective_and_zero_excludes_an_outlier() {
    let matrix = csr(3, 1, &[1.0, 1.0, 1.0]);
    let operator = CsrRectOp::new(&matrix);
    let problem = RegularizedLeastSquares::new(
        &operator,
        &[0.0, 10.0, 1000.0],
        &[2.0, 1.0, 0.0],
        &[0.0],
        &[0.0],
    )
    .unwrap();
    let result = problem.solve(LsqrConfig::default()).unwrap();
    assert_eq!(result.report.stop, Some(LsqrStop::LeastSquares));
    // (4 * 0 + 1 * 10) / (4 + 1), NOT (2 * 0 + 1 * 10) / 3.
    assert!((result.x[0] - 2.0).abs() < 1e-10);
}

#[test]
fn damping_is_about_the_explicit_nonzero_prior() {
    let matrix = csr(1, 1, &[1.0]);
    let operator = CsrRectOp::new(&matrix);
    let problem = RegularizedLeastSquares::damped(&operator, &[10.0], 2.0, &[3.0]).unwrap();
    let result = problem.solve(LsqrConfig::default()).unwrap();
    // (10 + 2^2 * 3) / (1 + 2^2).
    assert!((result.x[0] - 4.4).abs() < 1e-10);
    assert_eq!(result.report.stop, Some(LsqrStop::LeastSquares));
}

#[test]
fn diagonal_regularization_resolves_a_rank_deficient_observation() {
    let matrix = csr(1, 2, &[1.0, 1.0]);
    let operator = CsrRectOp::new(&matrix);
    let problem = RegularizedLeastSquares::new(
        &operator,
        &[2.0],
        &[1.0],
        &[1.0, 2.0],
        &[0.0, 0.0],
    )
    .unwrap();
    let result = problem.solve(LsqrConfig::default()).unwrap();
    assert!(result.report.converged());
    assert!((result.x[0] - 8.0 / 9.0).abs() < 1e-10);
    assert!((result.x[1] - 2.0 / 9.0).abs() < 1e-10);
}

#[test]
fn augmented_transpose_obeys_the_dot_identity_and_resume_matches() {
    let matrix = csr(3, 2, &[1.0, 2.0, -3.0, 4.0, 0.5, -1.0]);
    let operator = CsrRectOp::new(&matrix);
    let problem = RegularizedLeastSquares::new(
        &operator,
        &[1.0, 2.0, 3.0],
        &[2.0, 0.0, 3.0],
        &[0.5, 2.0],
        &[1.0, -1.0],
    )
    .unwrap();
    let x = [0.2, -0.7];
    let y = [1.0, -2.0, 0.5, 4.0, -3.0];
    let mut ax = [0.0; 5];
    let mut aty = [0.0; 2];
    problem.apply(&x, &mut ax);
    problem.apply_transpose(&y, &mut aty);
    assert!((fs_solver::dot(&ax, &y) - fs_solver::dot(&x, &aty)).abs() < 1e-12);
    let mut state = problem.start(LsqrConfig::default()).unwrap();
    assert!(state.step(&problem).unwrap());
    let mut resumed = state.clone();
    state.run(&problem).unwrap();
    resumed.run(&problem).unwrap();
    assert_eq!(state, resumed);
    assert!(resumed.report().converged());
}

#[test]
fn invalid_weight_prior_dimension_and_overflow_are_refused() {
    let matrix = csr(1, 1, &[1.0]);
    let operator = CsrRectOp::new(&matrix);
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            RegularizedLeastSquares::damped(&operator, &[1.0], value, &[0.0]),
            Err(LeastSquaresBuildError::InvalidValue { .. })
        ));
        assert!(
            RegularizedLeastSquares::new(&operator, &[1.0], &[value], &[0.0], &[0.0]).is_err()
        );
    }
    assert!(matches!(
        RegularizedLeastSquares::damped(&operator, &[1.0], 1.0, &[]),
        Err(LeastSquaresBuildError::Dimension { .. })
    ));
    assert!(RegularizedLeastSquares::damped(&operator, &[1.0], 0.0, &[f64::NAN]).is_err());
    assert!(RegularizedLeastSquares::damped(&operator, &[1.0], 1e200, &[1e200]).is_err());
    assert!(
        RegularizedLeastSquares::new(&operator, &[1e200], &[1e200], &[0.0], &[0.0]).is_err()
    );
}
