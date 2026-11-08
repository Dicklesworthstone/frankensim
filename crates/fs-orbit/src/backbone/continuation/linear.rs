//! Numerical admission around fs-la's original rectangular QR, not a second
//! factorization. Unit-length columns keep its squared Householder sums finite.

use super::{BackboneContinuationError, OrbitError, Result, invalid};
use fs_la::factor::{Qr, qr};
use fs_math::det;

pub(super) struct Admitted {
    factor: Qr,
    maxima: Vec<f64>,
    norms: Vec<f64>,
    pub(super) condition: f64,
}
impl Admitted {
    pub(super) fn solve(&self, residual: &[f64]) -> Result<Vec<f64>> {
        let mut delta = self.factor.solve_ls(residual);
        for ((value, &norm), &maximum) in delta.iter_mut().zip(&self.norms).zip(&self.maxima) {
            // The product maximum*norm may overflow despite a finite correction.
            *value = (*value / norm) / maximum;
        }
        if delta.iter().any(|v| !v.is_finite()) {
            return Err(invalid("QR correction is not representable"));
        }
        Ok(delta)
    }
}

pub(super) fn admit(
    matrix: &mut [f64],
    rows: usize,
    cols: usize,
    maximum_condition: f64,
) -> Result<Admitted> {
    if matrix.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "original corrector Jacobian produced nonfinite values",
        ));
    }
    let mut maxima = vec![0.0_f64; cols];
    let mut norms = vec![0.0; cols];
    for col in 0..cols {
        for row in 0..rows {
            maxima[col] = maxima[col].max(matrix[row * cols + col].abs());
        }
        if maxima[col] == 0.0 {
            return Err(OrbitError::SingularJacobian.into());
        }
        let mut square = 0.0;
        for row in 0..rows {
            matrix[row * cols + col] /= maxima[col];
            square += matrix[row * cols + col] * matrix[row * cols + col];
        }
        norms[col] = det::sqrt(square);
        for row in 0..rows {
            matrix[row * cols + col] /= norms[col];
        }
    }
    let factor = qr(matrix, rows, cols);
    let floor = 64.0 * f64::EPSILON * rows as f64;
    for col in 0..cols {
        for row in 0..=col {
            if !factor.r(row, col).is_finite() {
                return Err(OrbitError::SingularJacobian.into());
            }
        }
        if factor.r(col, col).abs() <= floor {
            return Err(OrbitError::SingularJacobian.into());
        }
    }
    let condition = condition_1(&factor, cols);
    if !condition.is_finite() || condition > maximum_condition {
        return Err(BackboneContinuationError::IllConditioned {
            estimate: condition,
            maximum: maximum_condition,
        });
    }
    Ok(Admitted {
        factor,
        maxima,
        norms,
        condition,
    })
}

// The same bounded Hager/Higham iteration used by fs-la::Lu::condition_1,
// applied directly to this QR's triangular R. No extra matrix or normal
// equation is formed. Lowest-index tie breaking is deterministic.
fn condition_1(factor: &Qr, n: usize) -> f64 {
    let norm = (0..n)
        .map(|j| (0..=j).map(|i| factor.r(i, j).abs()).sum::<f64>())
        .fold(0.0_f64, f64::max);
    let mut x = vec![1.0 / n as f64; n];
    let mut z = vec![0.0; n];
    let mut estimate = 0.0;
    for _ in 0..5 {
        for i in (0..n).rev() {
            let mut value = x[i];
            for (j, &entry) in x.iter().enumerate().skip(i + 1) {
                value -= factor.r(i, j) * entry;
            }
            x[i] = value / factor.r(i, i);
        }
        let next: f64 = x.iter().map(|v| v.abs()).sum();
        if !next.is_finite() {
            return f64::INFINITY;
        }
        if next <= estimate {
            break;
        }
        estimate = next;
        for i in 0..n {
            let mut value = if x[i] >= 0.0 { 1.0 } else { -1.0 };
            for (j, &entry) in z.iter().enumerate().take(i) {
                value -= factor.r(j, i) * entry;
            }
            z[i] = value / factor.r(i, i);
        }
        let (mut largest, mut index) = (0.0, 0);
        for (j, &value) in z.iter().enumerate() {
            if !value.is_finite() {
                return f64::INFINITY;
            }
            if value.abs() > largest {
                (largest, index) = (value.abs(), j);
            }
        }
        x.fill(0.0);
        x[index] = 1.0;
    }
    estimate * norm
}
