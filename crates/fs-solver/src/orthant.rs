//! Bounded nonnegative quadratic minimization over a supplied Gram factor.
//!
//! Minimize `0.5*x^T*(A*A^T)*x + b^T*x`, subject to `x >= 0`.
//! The factor makes the mathematical Hessian positive semidefinite without
//! assuming independent rows. Deterministic projected coordinate sweeps are
//! a baseline, not a general conic solver or an infeasibility classifier.
//! Every returned point passes a freshly reconstructed KKT residual. All
//! iterates are scratch: cancellation or nonconvergence returns no solution.

/// Explicit allocation/work and gradient-residual limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GramOrthantConfig {
    /// Maximum rows/unknowns, in 1..=64.
    pub max_variables: usize,
    /// Maximum columns of the factor, in 1..=4096.
    pub max_columns: usize,
    /// Conservative setup-product ceiling `rows*rows*columns`.
    pub max_setup_products: usize,
    /// Maximum coordinate sweeps, in 1..=4096; each uses O(rows^2) work.
    pub max_sweeps: usize,
    /// Positive finite absolute KKT tolerance, in gradient units.
    pub gradient_tolerance: f64,
}

/// A checked nonnegative point, not a uniqueness or rank certificate.
#[derive(Clone, Debug, PartialEq)]
pub struct GramOrthantSolution {
    /// Nonnegative primal coordinates.
    pub point: Vec<f64>,
    /// Independently reconstructed `A*A^T*x + b` at the returned point.
    pub gradient: Vec<f64>,
    /// Maximum active stationarity/inactive dual-feasibility error.
    pub residual: f64,
    /// Completed coordinate sweeps (zero for an already optimal zero point).
    pub sweeps: usize,
}

/// A refusal never returns a partially solved point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GramOrthantError {
    /// Invalid caps, shapes, nonfinite data, or a zero/unrepresentable row norm.
    Invalid(&'static str),
    /// Caller cancellation at a bounded setup/coordinate/residual boundary.
    Cancelled,
    /// Finite inputs led to unrepresentable arithmetic.
    NonFinite,
    /// Budget exhausted; does NOT classify infeasibility or nonuniqueness.
    NotConverged { residual: f64, sweeps: usize },
}
impl core::fmt::Display for GramOrthantError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Gram orthant solve: {self:?}")
    }
}
impl std::error::Error for GramOrthantError {}

/// Solve a PSD nonnegative quadratic without inverting its possibly singular
/// Gram matrix. Rows are fixed in caller order. Cancellation is caller owned.
///
/// No diagonal regularization, hidden coefficient clipping, warm-start memory,
/// or rank-dependent uniqueness claim is introduced. Empty problems are legal.
pub fn solve_gram_orthant(
    factor: &[Vec<f64>],
    linear: &[f64],
    config: GramOrthantConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<GramOrthantSolution, GramOrthantError> {
    poll(&mut cancelled)?;
    if !(1..=64).contains(&config.max_variables)
        || !(1..=4096).contains(&config.max_columns)
        || !(1..=4096).contains(&config.max_sweeps)
        || !config.gradient_tolerance.is_finite() || config.gradient_tolerance <= 0.0 {
        return Err(GramOrthantError::Invalid("configuration"));
    }
    let n = factor.len();
    if n > config.max_variables || linear.len() != n || linear.iter().any(|v| !v.is_finite()) {
        return Err(GramOrthantError::Invalid("row count or linear term"));
    }
    let d = factor.first().map_or(0, Vec::len);
    if n != 0 && (d == 0 || d > config.max_columns) {
        return Err(GramOrthantError::Invalid("column count"));
    }
    let products = n.checked_mul(n).and_then(|v| v.checked_mul(d))
        .ok_or(GramOrthantError::Invalid("setup overflow"))?;
    if products > config.max_setup_products {
        return Err(GramOrthantError::Invalid("setup product ceiling"));
    }
    for row in factor {
        poll(&mut cancelled)?;
        if row.len() != d || row.iter().any(|v| !v.is_finite()) {
            return Err(GramOrthantError::Invalid("factor row"));
        }
    }
    let mut matrix = vec![0.0; n*n];
    for i in 0..n {
        for j in 0..=i {
            poll(&mut cancelled)?;
            let entry = finite(crate::dot(&factor[i], &factor[j]))?;
            matrix[i*n+j] = entry;
            matrix[j*n+i] = entry;
        }
        if matrix[i*n+i] <= 0.0 {
            return Err(GramOrthantError::Invalid("positive representable row norm required"));
        }
    }
    let mut point = vec![0.0; n];
    let mut gradient = linear.to_vec();
    let mut residual = kkt(&point, &gradient);
    if residual <= config.gradient_tolerance {
        poll(&mut cancelled)?;
        return Ok(GramOrthantSolution { point, gradient, residual, sweeps: 0 });
    }
    for sweep in 1..=config.max_sweeps {
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut off_diagonal = linear[i];
            for j in 0..n {
                if i != j { off_diagonal = finite(off_diagonal + matrix[i*n+j]*point[j])?; }
            }
            point[i] = finite(-off_diagonal / matrix[i*n+i])?.max(0.0);
        }
        // Reconstruct the full residual, never merely a last-coordinate delta.
        for i in 0..n {
            poll(&mut cancelled)?;
            gradient[i] = finite(linear[i] + crate::dot(&matrix[i*n..(i+1)*n], &point))?;
        }
        residual = kkt(&point, &gradient);
        if residual <= config.gradient_tolerance {
            poll(&mut cancelled)?;
            return Ok(GramOrthantSolution { point, gradient, residual, sweeps: sweep });
        }
    }
    Err(GramOrthantError::NotConverged { residual, sweeps: config.max_sweeps })
}

fn finite(v: f64) -> Result<f64, GramOrthantError> {
    if v.is_finite() { Ok(v) } else { Err(GramOrthantError::NonFinite) }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), GramOrthantError> {
    if cancelled() { Err(GramOrthantError::Cancelled) } else { Ok(()) }
}
fn kkt(point: &[f64], gradient: &[f64]) -> f64 {
    point.iter().zip(gradient).fold(0.0_f64, |r, (&x, &g)| {
        r.max(if x > 0.0 { g.abs() } else { (-g).max(0.0) })
    })
}
