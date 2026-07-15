//! Typed block composition and exact-order block preconditioner plumbing.
//!
//! The existing [`crate::LinearOp`] is deliberately square. Multiphysics and
//! saddle systems also need rectangular coupling blocks, so this module adds a
//! narrow rectangular apply trait and const-arity square block operators. It
//! does not materialize matrices or infer physical meaning from block names.

use crate::LinearOp;
use fs_sparse::precond::Precond;
use std::fmt;

/// A rectangular matrix-free linear map with an explicit transpose.
pub trait RectLinearOp {
    /// Output dimension.
    fn rows(&self) -> usize;
    /// Input dimension.
    fn cols(&self) -> usize;
    /// `y = A x`.
    fn apply(&self, x: &[f64], y: &mut [f64]);
    /// `y = A^T x`.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]);
}

/// Explicit rectangular view of a square [`LinearOp`].
///
/// The adapter is intentional: a blanket implementation would give every
/// square operator two visible `apply` methods whenever both traits are in
/// scope, making ordinary method calls ambiguous for downstream users.
#[derive(Clone, Copy)]
pub struct SquareBlock<'a> {
    operator: &'a dyn LinearOp,
}

impl<'a> SquareBlock<'a> {
    /// Borrow a square operator as a rectangular block.
    #[must_use]
    pub const fn new(operator: &'a dyn LinearOp) -> Self {
        Self { operator }
    }
}

impl RectLinearOp for SquareBlock<'_> {
    fn rows(&self) -> usize {
        self.operator.n()
    }

    fn cols(&self) -> usize {
        self.operator.n()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        LinearOp::apply(self.operator, x, y);
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        LinearOp::apply_transpose(self.operator, x, y);
    }
}

/// An explicit rectangular zero block.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ZeroBlock {
    rows: usize,
    cols: usize,
}

impl ZeroBlock {
    /// Construct a zero map of shape `rows x cols`.
    #[must_use]
    pub const fn new(rows: usize, cols: usize) -> Self {
        Self { rows, cols }
    }
}

impl RectLinearOp for ZeroBlock {
    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.cols, "zero-block input length mismatch");
        assert_eq!(y.len(), self.rows, "zero-block output length mismatch");
        y.fill(0.0);
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.rows, "zero-block transpose input mismatch");
        assert_eq!(y.len(), self.cols, "zero-block transpose output mismatch");
        y.fill(0.0);
    }
}

/// Shape or arithmetic refusal from block construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlockError {
    /// A zero-by-zero block operator was requested.
    Empty,
    /// Blocks in one block row disagreed on output dimension.
    RowMismatch {
        /// Block row.
        row: usize,
        /// Expected output dimension.
        expected: usize,
        /// Observed output dimension.
        actual: usize,
    },
    /// Blocks in one block column disagreed on input dimension.
    ColumnMismatch {
        /// Block column.
        column: usize,
        /// Expected input dimension.
        expected: usize,
        /// Observed input dimension.
        actual: usize,
    },
    /// Total block rows and columns differed.
    NotSquare {
        /// Total rows.
        rows: usize,
        /// Total columns.
        cols: usize,
    },
    /// Dimension addition overflowed.
    DimensionOverflow,
    /// A real/imaginary operator pair differed in dimension.
    ComplexPairMismatch {
        /// Real-operator dimension.
        real: usize,
        /// Imaginary-operator dimension.
        imaginary: usize,
    },
    /// A Schur coupling block did not match the declared split.
    SchurCouplingMismatch {
        /// Name of the coupling (`lower` or `upper`).
        coupling: &'static str,
        /// Expected rows.
        expected_rows: usize,
        /// Expected columns.
        expected_cols: usize,
        /// Actual rows.
        actual_rows: usize,
        /// Actual columns.
        actual_cols: usize,
    },
}

impl fmt::Display for BlockError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => f.write_str("block operator must contain at least one block"),
            Self::RowMismatch {
                row,
                expected,
                actual,
            } => write!(
                f,
                "block row {row} expected {expected} rows, observed {actual}"
            ),
            Self::ColumnMismatch {
                column,
                expected,
                actual,
            } => write!(
                f,
                "block column {column} expected {expected} columns, observed {actual}"
            ),
            Self::NotSquare { rows, cols } => {
                write!(f, "block operator is {rows}x{cols}, not square")
            }
            Self::DimensionOverflow => f.write_str("block dimension arithmetic overflowed"),
            Self::ComplexPairMismatch { real, imaginary } => write!(
                f,
                "real-equivalent pair dimensions differ: real {real}, imaginary {imaginary}"
            ),
            Self::SchurCouplingMismatch {
                coupling,
                expected_rows,
                expected_cols,
                actual_rows,
                actual_cols,
            } => write!(
                f,
                "Schur {coupling} coupling expected {expected_rows}x{expected_cols}, \
                 observed {actual_rows}x{actual_cols}"
            ),
        }
    }
}

impl core::error::Error for BlockError {}

/// A square, const-arity block operator.
///
/// `N` is the number of block rows and columns, while individual blocks may be
/// rectangular. Construction verifies every partition once; apply then uses a
/// fixed row-major accumulation order.
pub struct BlockOperator<'a, const N: usize> {
    blocks: [[&'a dyn RectLinearOp; N]; N],
    row_sizes: [usize; N],
    col_sizes: [usize; N],
    row_offsets: [usize; N],
    col_offsets: [usize; N],
    n: usize,
}

/// Typed 2-by-2 block operator.
pub type BlockOperator2<'a> = BlockOperator<'a, 2>;
/// Typed 3-by-3 block operator.
pub type BlockOperator3<'a> = BlockOperator<'a, 3>;

impl<'a, const N: usize> BlockOperator<'a, N> {
    /// Validate and retain a block table.
    #[must_use]
    pub fn new(blocks: [[&'a dyn RectLinearOp; N]; N]) -> Result<Self, BlockError> {
        if N == 0 {
            return Err(BlockError::Empty);
        }
        let mut row_sizes = [0usize; N];
        let mut col_sizes = [0usize; N];
        for row in 0..N {
            row_sizes[row] = blocks[row][0].rows();
            for column in 1..N {
                let actual = blocks[row][column].rows();
                if actual != row_sizes[row] {
                    return Err(BlockError::RowMismatch {
                        row,
                        expected: row_sizes[row],
                        actual,
                    });
                }
            }
        }
        for column in 0..N {
            col_sizes[column] = blocks[0][column].cols();
            for row in 1..N {
                let actual = blocks[row][column].cols();
                if actual != col_sizes[column] {
                    return Err(BlockError::ColumnMismatch {
                        column,
                        expected: col_sizes[column],
                        actual,
                    });
                }
            }
        }
        let rows = checked_total(&row_sizes)?;
        let cols = checked_total(&col_sizes)?;
        if rows != cols {
            return Err(BlockError::NotSquare { rows, cols });
        }
        Ok(Self {
            blocks,
            row_sizes,
            col_sizes,
            row_offsets: offsets(&row_sizes)?,
            col_offsets: offsets(&col_sizes)?,
            n: rows,
        })
    }

    /// Block-row dimensions.
    #[must_use]
    pub const fn row_sizes(&self) -> &[usize; N] {
        &self.row_sizes
    }

    /// Block-column dimensions.
    #[must_use]
    pub const fn col_sizes(&self) -> &[usize; N] {
        &self.col_sizes
    }
}

fn checked_total<const N: usize>(sizes: &[usize; N]) -> Result<usize, BlockError> {
    sizes.iter().try_fold(0usize, |total, size| {
        total
            .checked_add(*size)
            .ok_or(BlockError::DimensionOverflow)
    })
}

fn offsets<const N: usize>(sizes: &[usize; N]) -> Result<[usize; N], BlockError> {
    let mut result = [0usize; N];
    let mut next = 0usize;
    for (index, size) in sizes.iter().copied().enumerate() {
        result[index] = next;
        next = next
            .checked_add(size)
            .ok_or(BlockError::DimensionOverflow)?;
    }
    Ok(result)
}

impl<const N: usize> LinearOp for BlockOperator<'_, N> {
    fn n(&self) -> usize {
        self.n
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.n, "block input length mismatch");
        assert_eq!(y.len(), self.n, "block output length mismatch");
        y.fill(0.0);
        for row in 0..N {
            let y_start = self.row_offsets[row];
            let y_row = &mut y[y_start..y_start + self.row_sizes[row]];
            for column in 0..N {
                let x_start = self.col_offsets[column];
                let x_column = &x[x_start..x_start + self.col_sizes[column]];
                let mut contribution = vec![0.0; self.row_sizes[row]];
                self.blocks[row][column].apply(x_column, &mut contribution);
                for (value, addend) in y_row.iter_mut().zip(contribution) {
                    *value += addend;
                }
            }
        }
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.n, "block transpose input mismatch");
        assert_eq!(y.len(), self.n, "block transpose output mismatch");
        y.fill(0.0);
        for row in 0..N {
            let x_start = self.row_offsets[row];
            let x_row = &x[x_start..x_start + self.row_sizes[row]];
            for column in 0..N {
                let y_start = self.col_offsets[column];
                let y_column = &mut y[y_start..y_start + self.col_sizes[column]];
                let mut contribution = vec![0.0; self.col_sizes[column]];
                self.blocks[row][column].apply_transpose(x_row, &mut contribution);
                for (value, addend) in y_column.iter_mut().zip(contribution) {
                    *value += addend;
                }
            }
        }
    }
}

/// Real-f64 representation of the complex operator `A + iB`:
/// `[[A, -B], [B, A]]`.
pub struct RealEquivalentComplexOp<'a> {
    real: &'a dyn LinearOp,
    imaginary: &'a dyn LinearOp,
    n: usize,
}

impl<'a> RealEquivalentComplexOp<'a> {
    /// Validate a same-dimension real/imaginary pair.
    #[must_use]
    pub fn new(real: &'a dyn LinearOp, imaginary: &'a dyn LinearOp) -> Result<Self, BlockError> {
        if real.n() != imaginary.n() {
            return Err(BlockError::ComplexPairMismatch {
                real: real.n(),
                imaginary: imaginary.n(),
            });
        }
        let n = real
            .n()
            .checked_mul(2)
            .ok_or(BlockError::DimensionOverflow)?;
        Ok(Self { real, imaginary, n })
    }

    /// Dimension of one real or imaginary component.
    #[must_use]
    pub fn component_dimension(&self) -> usize {
        self.n / 2
    }
}

impl LinearOp for RealEquivalentComplexOp<'_> {
    fn n(&self) -> usize {
        self.n
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.n, "real-equivalent input mismatch");
        assert_eq!(y.len(), self.n, "real-equivalent output mismatch");
        let half = self.n / 2;
        let (xr, xi) = x.split_at(half);
        let (yr, yi) = y.split_at_mut(half);
        let mut ar = vec![0.0; half];
        let mut ai = vec![0.0; half];
        let mut br = vec![0.0; half];
        let mut bi = vec![0.0; half];
        LinearOp::apply(self.real, xr, &mut ar);
        LinearOp::apply(self.real, xi, &mut ai);
        LinearOp::apply(self.imaginary, xr, &mut br);
        LinearOp::apply(self.imaginary, xi, &mut bi);
        for index in 0..half {
            yr[index] = ar[index] - bi[index];
            yi[index] = br[index] + ai[index];
        }
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        assert_eq!(x.len(), self.n, "real-equivalent transpose input mismatch");
        assert_eq!(y.len(), self.n, "real-equivalent transpose output mismatch");
        let half = self.n / 2;
        let (xr, xi) = x.split_at(half);
        let (yr, yi) = y.split_at_mut(half);
        let mut atr = vec![0.0; half];
        let mut ati = vec![0.0; half];
        let mut btr = vec![0.0; half];
        let mut bti = vec![0.0; half];
        LinearOp::apply_transpose(self.real, xr, &mut atr);
        LinearOp::apply_transpose(self.real, xi, &mut ati);
        LinearOp::apply_transpose(self.imaginary, xr, &mut btr);
        LinearOp::apply_transpose(self.imaginary, xi, &mut bti);
        for index in 0..half {
            yr[index] = atr[index] + bti[index];
            yi[index] = ati[index] - btr[index];
        }
    }
}

/// Sign applied to the Schur solve in a block-LDU inverse.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SchurSolveSign {
    /// `p = S^-1 (r1 - B A^-1 r0)`.
    Positive,
    /// `p = -S^-1 (r1 - B A^-1 r0)`, as for `[A B^T; B -C]`
    /// with positive complement `S = C + B A^-1 B^T`.
    Negative,
}

/// Two-by-two block-LDU Schur preconditioner assembled from injected inverse
/// approximations and rectangular coupling actions.
pub struct BlockSchur2<'a> {
    first: usize,
    second: usize,
    first_inverse: &'a dyn Precond,
    schur_inverse: &'a dyn Precond,
    lower: &'a dyn RectLinearOp,
    upper: &'a dyn RectLinearOp,
    sign: SchurSolveSign,
}

impl<'a> BlockSchur2<'a> {
    /// Validate the split and retain injected inverse approximations.
    #[must_use]
    pub fn new(
        first: usize,
        second: usize,
        first_inverse: &'a dyn Precond,
        schur_inverse: &'a dyn Precond,
        lower: &'a dyn RectLinearOp,
        upper: &'a dyn RectLinearOp,
        sign: SchurSolveSign,
    ) -> Result<Self, BlockError> {
        check_coupling("lower", lower, second, first)?;
        check_coupling("upper", upper, first, second)?;
        first
            .checked_add(second)
            .ok_or(BlockError::DimensionOverflow)?;
        Ok(Self {
            first,
            second,
            first_inverse,
            schur_inverse,
            lower,
            upper,
            sign,
        })
    }
}

fn check_coupling(
    name: &'static str,
    coupling: &dyn RectLinearOp,
    rows: usize,
    cols: usize,
) -> Result<(), BlockError> {
    if coupling.rows() == rows && coupling.cols() == cols {
        Ok(())
    } else {
        Err(BlockError::SchurCouplingMismatch {
            coupling: name,
            expected_rows: rows,
            expected_cols: cols,
            actual_rows: coupling.rows(),
            actual_cols: coupling.cols(),
        })
    }
}

impl Precond for BlockSchur2<'_> {
    fn apply(&self, residual: &[f64], output: &mut [f64]) {
        let total = self.first + self.second;
        assert_eq!(residual.len(), total, "Schur residual length mismatch");
        assert_eq!(output.len(), total, "Schur output length mismatch");
        let (r0, r1) = residual.split_at(self.first);

        let mut z0 = vec![0.0; self.first];
        self.first_inverse.apply(r0, &mut z0);
        let mut lower_z0 = vec![0.0; self.second];
        self.lower.apply(&z0, &mut lower_z0);
        let schur_rhs: Vec<f64> = r1
            .iter()
            .zip(lower_z0)
            .map(|(right, coupled)| right - coupled)
            .collect();
        let mut z1 = vec![0.0; self.second];
        self.schur_inverse.apply(&schur_rhs, &mut z1);
        if self.sign == SchurSolveSign::Negative {
            for value in &mut z1 {
                *value = -*value;
            }
        }

        let mut upper_z1 = vec![0.0; self.first];
        self.upper.apply(&z1, &mut upper_z1);
        let mut correction = vec![0.0; self.first];
        self.first_inverse.apply(&upper_z1, &mut correction);
        for (value, remove) in z0.iter_mut().zip(correction) {
            *value -= remove;
        }
        output[..self.first].copy_from_slice(&z0);
        output[self.first..].copy_from_slice(&z1);
    }
}
