//! Algebraic goal error and residual-evaluation roundoff for a retained CSR
//! system. This is a DISCRETE, stored-coefficient statement, not an estimate
//! of geometry, assembly, discretization, material or model-form error.
//!
//! For `A x* = b`, `A^T z* = g`, `r = b - A x`, and `s = g - A^T z`,
//!
//! ```text
//! g^T (x* - x) = z^T r + s^T A^{-1} r.
//! ```
//!
//! The evaluator encloses `z^T r` with outward-rounded FMA, including
//! cancellation and subnormal products. It never mistakes an approximately
//! solved dual for an exact one. A full algebraic goal enclosure additionally
//! requires a proved inverse bound. We establish one by positive diagonal
//! scaling: if `delta = min_i (|A_ii| w_i - sum_{j != i} |A_ij| w_j) > 0`,
//! then `||A^{-1}||_infinity <= max(w)/delta`. To see this, choose an index
//! maximizing `|v_i|/w_i` and apply the triangle inequality to `(A v)_i`.
//! The same argument proves nonsingularity, without assuming symmetry or SPD.
//!
//! Every margin is rounded DOWN and every norm/error bound UP. Failure to
//! prove dominance leaves the residual enclosure usable but the full goal
//! enclosure absent. Callers must not fill that absence with zero. A supplied
//! positive scaling can establish dominance when unscaled row sums cannot;
//! it is evidence to CHECK, not a caller-declared stability constant.

use fs_sparse::Csr;

mod arithmetic;
use arithmetic::{add_up, down, finite, mul_up, up};

/// Work admission before allocation or numerical traversal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GoalResidualLimits {
    /// Maximum square-system dimension.
    pub max_rows: usize,
    /// Maximum number of stored coefficients, including explicit zeros.
    pub max_nonzeros: usize,
}

/// A refusal does not publish a partial residual or goal enclosure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GoalResidualError {
    /// A nonempty square system is required.
    Shape { rows: usize, columns: usize },
    /// An input vector has the wrong length.
    Length { field: &'static str, expected: usize, found: usize },
    /// A declared structural limit was exceeded.
    Limit { field: &'static str, required: usize, allowed: usize },
    /// An input scalar was nonfinite.
    NonFinite { field: &'static str, index: usize },
    /// The checked scaling must be strictly positive everywhere.
    NonPositiveScaling { index: usize },
    /// Finite outward endpoints could not be represented.
    ArithmeticRange,
    /// Scratch storage could not be admitted or allocated.
    Allocation,
    /// The caller's checkpoint requested cancellation.
    Cancelled,
}

impl std::fmt::Display for GoalResidualError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "algebraic goal enclosure refused: {self:?}")
    }
}

impl std::error::Error for GoalResidualError {}

/// An immutable finite enclosure produced by the evaluator.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScalarEnclosure {
    lower: f64,
    upper: f64,
}

impl ScalarEnclosure {
    /// Lower outward endpoint.
    #[must_use]
    pub const fn lower(self) -> f64 {
        self.lower
    }

    /// Upper outward endpoint.
    #[must_use]
    pub const fn upper(self) -> f64 {
        self.upper
    }

    /// Upper bound on the absolute value, with no additional arithmetic.
    #[must_use]
    pub fn magnitude_upper(self) -> f64 {
        self.lower.abs().max(self.upper.abs())
    }
}

/// Why a full algebraic goal enclosure is present or absent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GoalBoundStatus {
    /// The inverse bound and dual-error term were established.
    Enclosed,
    /// Positive scaled strict diagonal dominance was not established.
    InverseBoundUnavailable,
    /// The inverse or dual-remainder bound exceeds finite binary64 range.
    BoundNotRepresentable,
}

/// Immutable numerical results for exactly the supplied stored system and
/// vectors. This report carries no continuum or nonlinear-operator authority.
#[derive(Debug, Clone, PartialEq)]
pub struct GoalResidualReport {
    weighted_residual: ScalarEnclosure,
    nominal_correction: f64,
    evaluation_roundoff_upper: f64,
    primal_residual_infinity_upper: f64,
    dual_residual_one_upper: f64,
    inverse_infinity_upper: Option<f64>,
    dual_error_upper: Option<f64>,
    goal_error: Option<ScalarEnclosure>,
    status: GoalBoundStatus,
}

impl GoalResidualReport {
    /// Enclosure of `z^T (b - A x)` for the supplied approximate dual.
    #[must_use]
    pub const fn weighted_residual(&self) -> ScalarEnclosure {
        self.weighted_residual
    }
    /// Rounded CSR-order residual followed by a rounded FMA dot product.
    #[must_use]
    pub const fn nominal_correction(&self) -> f64 {
        self.nominal_correction
    }
    /// Error in evaluating that nominal correction, not total solver roundoff.
    #[must_use]
    pub const fn evaluation_roundoff_upper(&self) -> f64 {
        self.evaluation_roundoff_upper
    }
    /// Outward upper bound on `||b - A x||_infinity`.
    #[must_use]
    pub const fn primal_residual_infinity_upper(&self) -> f64 {
        self.primal_residual_infinity_upper
    }
    /// Outward upper bound on `||g - A^T z||_1` (the actual transpose).
    #[must_use]
    pub const fn dual_residual_one_upper(&self) -> f64 {
        self.dual_residual_one_upper
    }
    /// Verified stored-system inverse norm, never a caller-provided constant.
    #[must_use]
    pub const fn inverse_infinity_upper(&self) -> Option<f64> {
        self.inverse_infinity_upper
    }
    /// Upper bound on the omitted `s^T A^{-1} r` term, when established.
    #[must_use]
    pub const fn dual_error_upper(&self) -> Option<f64> {
        self.dual_error_upper
    }
    /// Enclosure of the TRUE discrete goal error, including dual-solve error.
    #[must_use]
    pub const fn goal_error(&self) -> Option<ScalarEnclosure> {
        self.goal_error
    }
    /// Absence of a full bound never invalidates or upgrades the residual bound.
    #[must_use]
    pub const fn status(&self) -> GoalBoundStatus {
        self.status
    }
}

struct Work<F> {
    checkpoint: F,
    left: usize,
}
impl<F: FnMut() -> bool> Work<F> {
    fn poll(&mut self) -> Result<(), GoalResidualError> {
        if (self.checkpoint)() {
            self.left = 512;
            Ok(())
        } else {
            Err(GoalResidualError::Cancelled)
        }
    }
    fn tick(&mut self) -> Result<(), GoalResidualError> {
        self.left -= 1;
        if self.left == 0 {
            self.poll()?;
        }
        Ok(())
    }
}

fn scratch<F: FnMut() -> bool>(
    n: usize,
    work: &mut Work<F>,
) -> Result<Vec<ScalarEnclosure>, GoalResidualError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| GoalResidualError::Allocation)?;
    for _ in 0..n {
        work.tick()?;
        values.push(ScalarEnclosure::point(0.0));
    }
    Ok(values)
}

/// Enclose the algebraic error of a linear goal on the supplied discrete
/// system. `scaling = None` tests unscaled strict diagonal dominance. A
/// positive supplied scaling is checked against THIS matrix on every call.
/// All inputs are immutable and all scratch is local, so cancellation cannot
/// publish partial authority or change a solver iterate.
///
/// The checkpoint is polled before work, after at most 512 scalar/entry
/// visits (including inside long sparse rows), and before publication. `false`
/// requests cancellation. Storage is one interval vector, O(rows); work is
/// O(rows + nonzeros), and no transpose matrix is materialized.
///
/// # Errors
/// Returns typed shape, budget, length, finite/scaling, allocation,
/// arithmetic-range or cancellation refusals. Inability to prove an inverse
/// bound is a reported no-bound state, not a fabricated certificate.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // One immutable primal/dual traversal keeps the rounding and stability checks together.
pub fn enclose_goal_error(
    matrix: &Csr,
    rhs: &[f64],
    primal: &[f64],
    goal: &[f64],
    dual: &[f64],
    scaling: Option<&[f64]>,
    limits: GoalResidualLimits,
    checkpoint: impl FnMut() -> bool,
) -> Result<GoalResidualReport, GoalResidualError> {
    let mut work = Work {
        checkpoint,
        left: 512,
    };
    work.poll()?;
    let n = matrix.nrows();
    if n == 0 || matrix.ncols() != n {
        return Err(GoalResidualError::Shape { rows: n, columns: matrix.ncols() });
    }
    for (field, required, allowed) in [
        ("rows", n, limits.max_rows),
        ("nonzeros", matrix.nnz(), limits.max_nonzeros),
    ] {
        if required > allowed { return Err(GoalResidualError::Limit { field, required, allowed }); }
    }
    for (field, values) in [("rhs", rhs), ("primal", primal), ("goal", goal), ("dual", dual)] {
        if values.len() != n {
            return Err(GoalResidualError::Length { field, expected: n, found: values.len() });
        }
        for (index, value) in values.iter().enumerate() {
            work.tick()?;
            if !value.is_finite() { return Err(GoalResidualError::NonFinite { field, index }); }
        }
    }
    let mut max_scale = 1.0_f64;
    if let Some(values) = scaling {
        if values.len() != n {
            return Err(GoalResidualError::Length { field: "scaling", expected: n, found: values.len() });
        }
        max_scale = 0.0;
        for (index, value) in values.iter().enumerate() {
            work.tick()?;
            if !value.is_finite() { return Err(GoalResidualError::NonFinite { field: "scaling", index }); }
            if *value <= 0.0 { return Err(GoalResidualError::NonPositiveScaling { index }); }
            max_scale = max_scale.max(*value);
        }
    }
    // Check representable scratch geometry before allocation.
    n.checked_mul(std::mem::size_of::<ScalarEnclosure>())
        .ok_or(GoalResidualError::Allocation)?;
    work.poll()?;
    let mut dual_residual = scratch(n, &mut work)?;
    for (slot, &value) in dual_residual.iter_mut().zip(goal) {
        work.tick()?;
        *slot = ScalarEnclosure::point(value);
    }
    let mut correction = ScalarEnclosure::point(0.0);
    let mut nominal = 0.0_f64;
    let mut primal_inf = 0.0_f64;
    let mut margin_lower = f64::INFINITY;
    let mut margin_representable = true;
    for row in 0..n {
        work.tick()?;
        let (columns, values) = matrix.row(row);
        let mut r = ScalarEnclosure::point(rhs[row]);
        let mut ax = 0.0_f64;
        let mut margin = 0.0_f64;
        for (&column, &a) in columns.iter().zip(values) {
            work.tick()?;
            if !a.is_finite() { return Err(GoalResidualError::NonFinite { field: "matrix row", index: row }); }
            r = r.add_scaled(ScalarEnclosure::point(primal[column]), -a)?;
            dual_residual[column] = dual_residual[column].add_scaled(ScalarEnclosure::point(dual[row]), -a)?;
            ax = finite(a.mul_add(primal[column], ax))?;
            if margin_representable && a != 0.0 {
                let w = scaling.map_or(1.0, |weights| weights[column]);
                let coefficient = if column == row { a.abs() } else { -a.abs() };
                match down(coefficient.mul_add(w, margin)) {
                    Ok(value) => margin = value,
                    Err(_) => margin_representable = false,
                }
            }
        }
        margin_lower = margin_lower.min(margin);
        primal_inf = primal_inf.max(r.magnitude_upper());
        correction = correction.add_scaled(r, dual[row])?;
        nominal = finite(dual[row].mul_add(finite(rhs[row] - ax)?, nominal))?;
    }
    let mut dual_one = 0.0;
    for value in dual_residual {
        work.tick()?;
        dual_one = add_up(dual_one, value.magnitude_upper())?;
    }
    let evaluation_roundoff_upper = correction.radius_about(nominal)?;
    let (mut inverse, mut dual_error, mut goal_error) = (None, None, None);
    let status = if !margin_representable {
        GoalBoundStatus::BoundNotRepresentable
    } else if margin_lower <= 0.0 {
        GoalBoundStatus::InverseBoundUnavailable
    } else {
        // The two residuals retain dual error. A small relative dual residual
        // by itself is not an inverse-stability or algebraic-error proof.
        let bound = (|| {
            let norm = up(max_scale / margin_lower)?;
            let remainder = mul_up(mul_up(dual_one, primal_inf)?, norm)?;
            let error = correction.widen(remainder)?;
            Ok::<_, GoalResidualError>((norm, remainder, error))
        })();
        match bound {
            Ok((norm, remainder, error)) => {
                inverse = Some(norm);
                dual_error = Some(remainder);
                goal_error = Some(error);
                GoalBoundStatus::Enclosed
            }
            Err(_) => GoalBoundStatus::BoundNotRepresentable,
        }
    };
    work.poll()?;
    Ok(GoalResidualReport {
        weighted_residual: correction,
        nominal_correction: nominal,
        evaluation_roundoff_upper,
        primal_residual_infinity_upper: primal_inf,
        dual_residual_one_upper: dual_one,
        inverse_infinity_upper: inverse,
        dual_error_upper: dual_error,
        goal_error,
        status,
    })
}
