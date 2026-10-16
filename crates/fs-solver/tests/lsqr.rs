//! Rectangular LSQR regressions with analytic answers and true residuals.
use fs_solver::{
    RectLinearOp,
    lsqr::{LsqrConfig, LsqrError, LsqrState, LsqrStop, lsqr},
};
use std::cell::Cell;

struct Dense {
    rows: usize,
    cols: usize,
    values: Vec<f64>,
}

impl Dense {
    fn new(rows: usize, cols: usize, values: &[f64]) -> Self {
        assert_eq!(values.len(), rows * cols);
        Self {
            rows,
            cols,
            values: values.to_vec(),
        }
    }
}

impl RectLinearOp for Dense {
    fn rows(&self) -> usize {
        self.rows
    }
    fn cols(&self) -> usize {
        self.cols
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.cols);
        assert_eq!(y.len(), self.rows);
        for (row, yi) in self.values.chunks_exact(self.cols).zip(y) {
            *yi = row.iter().zip(x).map(|(a, b)| a * b).sum();
        }
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.rows);
        assert_eq!(y.len(), self.cols);
        y.fill(0.0);
        for (row, xi) in self.values.chunks_exact(self.cols).zip(x) {
            for (a, yi) in row.iter().zip(y.iter_mut()) {
                *yi += a * xi;
            }
        }
    }
}

fn close(x: &[f64], expected: &[f64]) {
    assert_eq!(x.len(), expected.len());
    for (x, expected) in x.iter().zip(expected) {
        assert!((x - expected).abs() < 1e-9, "{x} != {expected}");
    }
}

#[test]
fn inconsistent_fit_stops_on_stationarity_not_a_fabricated_zero_residual() {
    let a = Dense::new(3, 2, &[1.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    let result = lsqr(&a, &[1.0, 2.0, 2.0], LsqrConfig::default()).unwrap();
    close(&result.x, &[2.0 / 3.0, 5.0 / 3.0]);
    assert_eq!(result.report.stop, Some(LsqrStop::LeastSquares));
    let r = result.report.history.last().unwrap();
    assert!((r.residual_norm - (1.0_f64 / 3.0).sqrt()).abs() < 1e-12);
    assert!(r.normal_residual_norm < 1e-10);
}

#[test]
fn underdetermined_solution_has_minimum_norm() {
    let a = Dense::new(2, 3, &[1.0, 0.0, 1.0, 0.0, 1.0, 1.0]);
    let result = lsqr(&a, &[1.0, 1.0], LsqrConfig::default()).unwrap();
    close(&result.x, &[1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0]);
    assert_eq!(result.report.stop, Some(LsqrStop::Compatible));
}

#[test]
fn rank_deficient_operator_does_not_require_invertible_normal_matrix() {
    let a = Dense::new(3, 2, &[1.0, 2.0, 2.0, 4.0, 3.0, 6.0]);
    let result = lsqr(&a, &[1.0, 2.0, 3.0], LsqrConfig::default()).unwrap();
    close(&result.x, &[0.2, 0.4]);
    assert!(result.report.converged());
}

#[test]
fn zero_rhs_and_zero_operator_are_distinct_successes() {
    let a = Dense::new(2, 1, &[0.0, 0.0]);
    let zero = lsqr(&a, &[0.0, 0.0], LsqrConfig::default()).unwrap();
    assert_eq!(zero.report.stop, Some(LsqrStop::Compatible));
    assert_eq!(zero.report.iterations, 0);
    let inconsistent = lsqr(&a, &[3.0, 4.0], LsqrConfig::default()).unwrap();
    assert_eq!(inconsistent.report.stop, Some(LsqrStop::LeastSquares));
    assert!((inconsistent.report.history[0].residual_norm - 5.0).abs() < 1e-14);
    close(&inconsistent.x, &[0.0]);
}

#[test]
fn finite_extreme_scales_do_not_overflow_or_underflow_the_norm() {
    for scale in [1e-200, 1e200] {
        let a = Dense::new(1, 1, &[scale]);
        let result = lsqr(&a, &[3.0], LsqrConfig::default()).unwrap();
        assert!(result.report.converged());
        assert!((scale * result.x[0] - 3.0).abs() < 1e-12);
        assert!(result.report.history[0].normal_residual_norm > 0.0);
    }
}

#[test]
fn checkpoints_resume_bitwise_and_budget_exhaustion_is_not_success() {
    let a = Dense::new(3, 3, &[1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0]);
    let b = [1.0, 1.0, 1.0];
    let mut straight = LsqrState::new(&a, &b, LsqrConfig::default()).unwrap();
    assert!(straight.step(&a).unwrap());
    assert_eq!(straight.report().stop, None);
    let mut resumed = straight.clone();
    straight.run(&a).unwrap();
    resumed.run(&a).unwrap();
    assert_eq!(straight, resumed);
    for (x, y) in straight.x().iter().zip(resumed.x()) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    assert!(straight.report().converged());
    close(straight.x(), &[1.0, 0.5, 1.0 / 3.0]);
    let stopped = straight.clone();
    assert!(!straight.step(&a).unwrap());
    assert_eq!(straight, stopped);
    for max_iters in [0, 1] {
        let config = LsqrConfig {
            max_iters,
            ..LsqrConfig::default()
        };
        let result = lsqr(&a, &b, config).unwrap();
        assert_eq!(result.report.stop, Some(LsqrStop::IterationLimit));
        assert!(!result.report.converged());
        assert_eq!(result.report.history.len(), max_iters + 1);
    }
}

struct Poisonable {
    a: Dense,
    poison: Cell<bool>,
}

impl RectLinearOp for Poisonable {
    fn rows(&self) -> usize {
        self.a.rows()
    }
    fn cols(&self) -> usize {
        self.a.cols()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.a.apply(x, y);
        if self.poison.get() {
            y[0] = f64::NAN;
        }
    }
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.a.apply_transpose(x, y);
    }
}

#[test]
fn failed_iteration_and_changed_shape_preserve_the_checkpoint() {
    let a = Poisonable {
        a: Dense::new(2, 2, &[1.0, 0.0, 0.0, 2.0]),
        poison: Cell::new(false),
    };
    let mut state = LsqrState::new(&a, &[1.0, 1.0], LsqrConfig::default()).unwrap();
    let checkpoint = state.clone();
    a.poison.set(true);
    assert!(matches!(state.step(&a), Err(LsqrError::NonFinite(_))));
    assert_eq!(state, checkpoint);
    let other = Dense::new(1, 1, &[1.0]);
    assert!(matches!(
        state.step(&other),
        Err(LsqrError::Dimension { .. })
    ));
    assert_eq!(state, checkpoint);
    a.poison.set(false);
    state.run(&a).unwrap();
    assert!(state.report().converged());
}

#[test]
fn invalid_inputs_are_refused() {
    let a = Dense::new(1, 1, &[1.0]);
    assert!(matches!(
        LsqrState::new(&a, &[], LsqrConfig::default()),
        Err(LsqrError::Dimension { .. })
    ));
    assert!(matches!(
        LsqrState::new(&a, &[f64::NAN], LsqrConfig::default()),
        Err(LsqrError::NonFinite(_))
    ));
    let empty = Dense::new(0, 1, &[]);
    assert!(matches!(
        LsqrState::new(&empty, &[], LsqrConfig::default()),
        Err(LsqrError::Empty)
    ));
    for tolerance in [f64::NAN, f64::INFINITY, -1.0, 1.0] {
        let config = LsqrConfig {
            relative_tolerance: tolerance,
            ..LsqrConfig::default()
        };
        assert!(matches!(
            LsqrState::new(&a, &[1.0], config),
            Err(LsqrError::InvalidTolerance(_))
        ));
    }
}
