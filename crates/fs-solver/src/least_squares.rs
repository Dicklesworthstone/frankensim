//! Weighted observations and diagonal Tikhonov priors for matrix-free LSQR.
//!
//! The problem is `min_x ||W(Ax-b)||^2 + ||D(x-prior)||^2`, where the
//! supplied diagonal entries of W and D multiply residuals, NOT squared
//! residuals. Use `W[i] = 1 / sigma[i]` for independent observation standard
//! deviations. The augmented operator `[W A; D]` is applied without building
//! an augmented sparse matrix or the normal equations. Residuals reported by
//! LSQR are for this AUGMENTED problem, not the unweighted physical data.

use crate::RectLinearOp;
use crate::lsqr::{LsqrConfig, LsqrError, LsqrSolution, LsqrState, lsqr};
use fs_sparse::Csr;
use std::fmt;

/// Borrowed rectangular CSR adapter; the transpose is applied without a copy.
#[derive(Debug, Clone, Copy)]
pub struct CsrRectOp<'a> {
    matrix: &'a Csr,
}

impl<'a> CsrRectOp<'a> {
    /// Borrow a canonical CSR matrix of any shape.
    #[must_use]
    pub const fn new(matrix: &'a Csr) -> Self {
        Self { matrix }
    }
}

impl RectLinearOp for CsrRectOp<'_> {
    fn rows(&self) -> usize {
        self.matrix.nrows()
    }

    fn cols(&self) -> usize {
        self.matrix.ncols()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.matrix.spmv(x, y);
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.rows(), "CSR transpose input length mismatch");
        assert_eq!(y.len(), self.cols(), "CSR transpose output length mismatch");
        y.fill(0.0);
        // For each output column, contributions arrive in ascending source-row
        // order. No hash iteration, atomics, or thread-dependent scatter order.
        for (row, xr) in x.iter().enumerate() {
            let (columns, values) = self.matrix.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                y[column] = value.mul_add(*xr, y[column]);
            }
        }
    }
}

/// Malformed weighted/regularized problem, refused before any operator apply.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LeastSquaresBuildError {
    /// The base operator is empty or its augmented row count would overflow.
    InvalidShape,
    /// A supplied vector has the wrong length.
    Dimension {
        /// Vector name.
        name: &'static str,
        /// Required number of entries.
        expected: usize,
        /// Supplied number of entries.
        actual: usize,
    },
    /// Non-finite data/prior/product, or a negative/non-finite diagonal entry.
    InvalidValue {
        /// Vector or product name.
        name: &'static str,
        /// Offending entry index.
        index: usize,
    },
}

impl fmt::Display for LeastSquaresBuildError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidShape => f.write_str("invalid augmented least-squares shape"),
            Self::Dimension {
                name,
                expected,
                actual,
            } => write!(
                f,
                "least-squares {name}: expected {expected}, observed {actual}"
            ),
            Self::InvalidValue { name, index } => {
                write!(f, "invalid least-squares {name} at index {index}")
            }
        }
    }
}

impl std::error::Error for LeastSquaresBuildError {}

/// Validated weighted observations and a diagonal prior penalty.
///
/// The underlying operator must remain immutable throughout a solve/resume.
/// Zero observation weights exclude those observations; zero penalties leave
/// those coordinates unregularized. These choices need not give a unique fit.
pub struct RegularizedLeastSquares<'a> {
    operator: &'a dyn RectLinearOp,
    weights: Vec<f64>,
    penalties: Vec<f64>,
    rhs: Vec<f64>,
}

impl<'a> RegularizedLeastSquares<'a> {
    /// Build `min ||diag(weights)(Ax-b)||^2 + ||diag(penalties)(x-prior)||^2`.
    /// All vectors are explicit: observation vectors have `rows` entries,
    /// and parameter vectors have `cols` entries. No values are silently clipped.
    pub fn new(
        operator: &'a dyn RectLinearOp,
        b: &[f64],
        weights: &[f64],
        penalties: &[f64],
        prior: &[f64],
    ) -> Result<Self, LeastSquaresBuildError> {
        let m = operator.rows();
        let n = operator.cols();
        let rows = m
            .checked_add(n)
            .ok_or(LeastSquaresBuildError::InvalidShape)?;
        if m == 0 || n == 0 {
            return Err(LeastSquaresBuildError::InvalidShape);
        }
        validate("observations", b, m, false)?;
        validate("weights", weights, m, true)?;
        validate("penalties", penalties, n, true)?;
        validate("prior", prior, n, false)?;
        let mut rhs = Vec::with_capacity(rows);
        rhs.extend(b.iter().zip(weights).map(|(value, weight)| value * weight));
        rhs.extend(
            prior
                .iter()
                .zip(penalties)
                .map(|(value, penalty)| value * penalty),
        );
        validate("augmented right-hand side", &rhs, rows, false)?;
        Ok(Self {
            operator,
            weights: weights.to_vec(),
            penalties: penalties.to_vec(),
            rhs,
        })
    }

    /// Scalar damping about an explicit prior: `||Ax-b||^2 + damping^2||x-prior||^2`.
    pub fn damped(
        operator: &'a dyn RectLinearOp,
        b: &[f64],
        damping: f64,
        prior: &[f64],
    ) -> Result<Self, LeastSquaresBuildError> {
        validate("damping", &[damping], 1, true)?;
        Self::new(
            operator,
            b,
            &vec![1.0; b.len()],
            &vec![damping; prior.len()],
            prior,
        )
    }

    /// The augmented right-hand side `[W b; D prior]` for a custom LSQR driver.
    #[must_use]
    pub fn right_hand_side(&self) -> &[f64] {
        &self.rhs
    }

    /// Start a checkpointable solve of the augmented system.
    pub fn start(&self, config: LsqrConfig) -> Result<LsqrState, LsqrError> {
        LsqrState::new(self, &self.rhs, config)
    }

    /// Solve with measured augmented residuals; inspect the report for convergence.
    pub fn solve(&self, config: LsqrConfig) -> Result<LsqrSolution, LsqrError> {
        lsqr(self, &self.rhs, config)
    }
}

impl RectLinearOp for RegularizedLeastSquares<'_> {
    fn rows(&self) -> usize {
        self.rhs.len()
    }

    fn cols(&self) -> usize {
        self.penalties.len()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.cols(), "augmented input length mismatch");
        assert_eq!(y.len(), self.rows(), "augmented output length mismatch");
        let (data, penalty) = y.split_at_mut(self.weights.len());
        self.operator.apply(x, data);
        for (value, weight) in data.iter_mut().zip(&self.weights) {
            *value *= weight;
        }
        for ((value, coefficient), xi) in penalty.iter_mut().zip(&self.penalties).zip(x) {
            *value = coefficient * xi;
        }
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.rows(), "augmented transpose input mismatch");
        assert_eq!(y.len(), self.cols(), "augmented transpose output mismatch");
        let (data, penalty) = x.split_at(self.weights.len());
        let weighted: Vec<f64> = data
            .iter()
            .zip(&self.weights)
            .map(|(x, w)| x * w)
            .collect();
        self.operator.apply_transpose(&weighted, y);
        for ((value, coefficient), xi) in y.iter_mut().zip(&self.penalties).zip(penalty) {
            *value += coefficient * xi;
        }
    }
}

fn validate(
    name: &'static str,
    values: &[f64],
    expected: usize,
    nonnegative: bool,
) -> Result<(), LeastSquaresBuildError> {
    if values.len() != expected {
        return Err(LeastSquaresBuildError::Dimension {
            name,
            expected,
            actual: values.len(),
        });
    }
    if let Some(index) = values
        .iter()
        .position(|v| !v.is_finite() || (nonnegative && *v < 0.0))
    {
        return Err(LeastSquaresBuildError::InvalidValue { name, index });
    }
    Ok(())
}
