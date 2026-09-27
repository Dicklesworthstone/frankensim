//! A sparse, whole-spectrum inverse bound, including obtuse-mesh operators.
//!
//! A numerical factor is only a proposal. For ANY real sparse columns L and
//! positive shift s, independently enclose E = A - s I - L L^T. If
//! eta >= max(||E||_1, ||E||_inf) and alpha = s - eta > 0, then
//!
//! ```text
//! v^T A v >= alpha ||v||_2^2,
//! ||A^-1||_inf <= sqrt(n) / alpha.
//! ```
//!
//! Indeed ||E||_2 <= sqrt(||E||_1 ||E||_inf) <= eta and L L^T is PSD.
//! This proves a lower bound for EVERY eigenvalue when A is symmetric, and
//! coercivity/nonsingularity even when stored assembly has small asymmetry.
//! An approximate Ritz pair alone cannot establish this proposition.
//!
//! The factor builder uses deterministic minimum-degree sparse elimination.
//! Neither successful numerical Cholesky nor its pivot sizes grant authority:
//! the separate outward Gram residual checks every stored coefficient AND
//! every fill entry. No dense inverse or dense n-by-n workspace is formed.

use std::collections::{BTreeMap, BTreeSet};

use fs_sparse::Csr;

use crate::goal::arithmetic::{add_up, down, mul_up, up};
use crate::goal::{
    GoalBoundStatus, GoalResidualError, GoalResidualLimits, GoalResidualReport,
    ScalarEnclosure, enclose_goal_error,
};

/// A sparse Gram column: distinct, strictly increasing row indices and values.
pub type GramColumn = Vec<(usize, f64)>;

/// Independent sparse storage/work admission, including unsuccessful shifts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SpectralInverseLimits {
    /// Ordinary input matrix envelope.
    pub system: GoalResidualLimits,
    /// Maximum simultaneously live logical sparse records plus reserved
    /// vector records. This is not a measured allocator/RSS bound; maps carry
    /// implementation overhead. The final owned CSR is included in admission.
    pub max_storage_entries: usize,
    /// Maximum scalar/edge visits across proposal, validation and verification.
    /// Ordered-map lookup overhead is not counted as additional scalar visits.
    pub max_work_entries: usize,
    /// Maximum shift attempts; each retry halves the previous proposal.
    pub max_shift_attempts: usize,
}

/// Why preparation stopped. No unsuccessful variant grants inverse authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectralStop {
    /// The full stored matrix passed the outward Gram-residual proof.
    Certified,
    /// No positive coercivity bound was established within the shift count.
    NotEstablished,
    /// The shared scalar/edge visit allowance ended before certification.
    WorkLimit,
    /// The live sparse-storage allowance ended before certification.
    StorageLimit,
    /// A finite candidate or outward endpoint could not be represented.
    ArithmeticRange,
}

/// Invalid input or explicit cancellation; neither publishes a certificate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpectralError {
    /// The ordinary matrix/vector/cancellation admission failed.
    Residual(GoalResidualError),
    /// A shift, attempt count or sparse proposal is malformed.
    InvalidProposal(&'static str),
}
impl From<GoalResidualError> for SpectralError {
    fn from(error: GoalResidualError) -> Self { Self::Residual(error) }
}
impl std::fmt::Display for SpectralError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "sparse spectral inverse refused: {self:?}")
    }
}
impl std::error::Error for SpectralError {}

/// Verified inverse bound owned together with the EXACT matrix it describes.
/// There is no public constructor from numbers, no mutable matrix accessor,
/// and no API applying this authority to a caller-selected replacement matrix.
/// Copying this object copies its matrix as well, not a detached stability
/// constant. The Gram workspace is discarded after successful verification.
#[derive(Debug, Clone)]
pub struct SpectralInverse {
    matrix: Csr,
    shift: f64,
    defect_upper: f64,
    coercivity_lower: f64,
    inverse_infinity_upper: f64,
}
impl SpectralInverse {
    /// Immutable, original stored coefficients, including explicit zeros.
    #[must_use]
    pub const fn matrix(&self) -> &Csr { &self.matrix }
    /// The numerical shift whose full residual was independently checked.
    #[must_use]
    pub const fn shift(&self) -> f64 { self.shift }
    /// Outward upper bound on max(||A-sI-LL^T||_1, ||A-sI-LL^T||_inf).
    #[must_use]
    pub const fn defect_upper(&self) -> f64 { self.defect_upper }
    /// Positive lower bound on coercivity of the exact stored operator.
    /// For a symmetric matrix this also bounds its smallest eigenvalue below.
    #[must_use]
    pub const fn coercivity_lower(&self) -> f64 { self.coercivity_lower }
    /// Verified infinity-norm inverse upper bound, for inspection only.
    #[must_use]
    pub const fn inverse_infinity_upper(&self) -> f64 { self.inverse_infinity_upper }

    /// Enclose a linear goal on THIS owned matrix without refactoring it.
    /// The ordinary outward primal/transpose-dual residuals are recomputed;
    /// the smaller proved inverse bound is used and dual error is retained.
    ///
    /// # Errors
    /// Ordinary residual input, budget, arithmetic and cancellation refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn enclose_goal_error(
        &self, rhs: &[f64], primal: &[f64], goal: &[f64], dual: &[f64],
        scaling: Option<&[f64]>, limits: GoalResidualLimits,
        mut checkpoint: impl FnMut() -> bool,
    ) -> Result<GoalResidualReport, GoalResidualError> {
        let mut report = enclose_goal_error(
            &self.matrix, rhs, primal, goal, dual, scaling, limits, &mut checkpoint,
        )?;
        if report.inverse_infinity_upper.is_none_or(|old| self.inverse_infinity_upper < old) {
            report.inverse_infinity_upper = Some(self.inverse_infinity_upper);
            let bounded = (|| {
                let remainder = mul_up(
                    mul_up(report.dual_residual_one_upper, report.primal_residual_infinity_upper)?,
                    self.inverse_infinity_upper,
                )?;
                Ok::<_, GoalResidualError>((remainder, report.weighted_residual.widen(remainder)?))
            })();
            match bounded {
                Ok((remainder, error)) => {
                    report.dual_error_upper = Some(remainder);
                    report.goal_error = Some(error);
                    report.status = GoalBoundStatus::Enclosed;
                }
                Err(GoalResidualError::ArithmeticRange) => {
                    report.dual_error_upper = None;
                    report.goal_error = None;
                    report.status = GoalBoundStatus::BoundNotRepresentable;
                }
                Err(error) => return Err(error),
            }
        }
        if !checkpoint() { return Err(GoalResidualError::Cancelled); }
        Ok(report)
    }
}

/// Preparation result, including paid work when a sufficient proof fails.
#[derive(Debug)]
pub struct SpectralPreparation {
    /// Present only after independent outward verification and final polling.
    pub certificate: Option<SpectralInverse>,
    /// Explicit successful or no-bound disposition.
    pub stop: SpectralStop,
    /// Actual charged visits; never resets on a smaller-shift retry.
    pub work_entries: usize,
    /// Peak admitted logical storage records, not measured bytes.
    pub peak_storage_entries: usize,
    /// Number of attempted shifts, including the final unsuccessful one.
    pub shift_attempts: usize,
}

#[derive(Debug)]
enum Halt { Work, Storage, Range, Fault(SpectralError) }
impl From<GoalResidualError> for Halt {
    fn from(error: GoalResidualError) -> Self {
        match error {
            GoalResidualError::ArithmeticRange => Self::Range,
            other => Self::Fault(other.into()),
        }
    }
}
struct Meter<F> {
    checkpoint: F,
    limits: SpectralInverseLimits,
    used: usize,
    live: usize,
    peak: usize,
    until_poll: usize,
}
impl<F: FnMut() -> bool> Meter<F> {
    fn poll(&mut self) -> Result<(), Halt> {
        self.until_poll = 256;
        if (self.checkpoint)() { Ok(()) }
        else { Err(Halt::Fault(GoalResidualError::Cancelled.into())) }
    }
    fn tick(&mut self) -> Result<(), Halt> {
        if self.used == self.limits.max_work_entries { return Err(Halt::Work); }
        self.used += 1;
        self.until_poll -= 1;
        if self.until_poll == 0 { self.poll()?; }
        Ok(())
    }
    fn grow(&mut self, count: usize) -> Result<(), Halt> {
        let next = self.live.checked_add(count).ok_or(Halt::Storage)?;
        if next > self.limits.max_storage_entries { return Err(Halt::Storage); }
        self.live = next;
        self.peak = self.peak.max(next);
        Ok(())
    }
    fn release(&mut self, count: usize) { self.live -= count; }
}
fn reserved<T>(n: usize) -> Result<Vec<T>, Halt> {
    let mut values = Vec::new();
    values.try_reserve_exact(n)
        .map_err(|_| Halt::Fault(GoalResidualError::Allocation.into()))?;
    Ok(values)
}
fn finite(value: f64) -> Result<f64, Halt> {
    if value.is_finite() { Ok(value) } else { Err(Halt::Range) }
}
fn admit<F: FnMut() -> bool>(a: &Csr, shift: f64, meter: &mut Meter<F>) -> Result<(), Halt> {
    meter.poll()?;
    let n = a.nrows();
    if n == 0 || a.ncols() != n {
        return Err(GoalResidualError::Shape { rows: n, columns: a.ncols() }.into());
    }
    for (field, required, allowed) in [
        ("rows", n, meter.limits.system.max_rows),
        ("nonzeros", a.nnz(), meter.limits.system.max_nonzeros),
    ] {
        if required > allowed {
            return Err(GoalResidualError::Limit { field, required, allowed }.into());
        }
    }
    if !shift.is_finite() || shift <= 0.0 || meter.limits.max_shift_attempts == 0 {
        return Err(Halt::Fault(SpectralError::InvalidProposal(
            "a positive finite shift and nonzero attempt count are required",
        )));
    }
    for row in 0..n {
        meter.tick()?;
        for &value in a.row(row).1 {
            meter.tick()?;
            if !value.is_finite() {
                return Err(GoalResidualError::NonFinite { field: "matrix row", index: row }.into());
            }
        }
    }
    Ok(())
}
fn finish<F: FnMut() -> bool>(
    mut meter: Meter<F>, attempts: usize, result: Result<Option<SpectralInverse>, Halt>,
) -> Result<SpectralPreparation, SpectralError> {
    // Cancellation takes precedence over an exhausted/inconclusive proposal.
    meter.poll().map_err(|halt| match halt {
        Halt::Fault(error) => error,
        _ => unreachable!("poll only observes cancellation"),
    })?;
    let (certificate, stop) = match result {
        Ok(Some(certificate)) => (Some(certificate), SpectralStop::Certified),
        Ok(None) => (None, SpectralStop::NotEstablished),
        Err(Halt::Work) => (None, SpectralStop::WorkLimit),
        Err(Halt::Storage) => (None, SpectralStop::StorageLimit),
        Err(Halt::Range) => (None, SpectralStop::ArithmeticRange),
        Err(Halt::Fault(error)) => return Err(error),
    };
    Ok(SpectralPreparation { certificate, stop, work_entries: meter.used,
        peak_storage_entries: meter.peak, shift_attempts: attempts })
}
fn meter<F: FnMut() -> bool>(limits: SpectralInverseLimits, checkpoint: F) -> Meter<F> {
    Meter { checkpoint, limits, used: 0, live: 0, peak: 0, until_poll: 256 }
}

/// Check arbitrary sparse Gram proposals without trusting their construction.
/// Every column must use canonical distinct row indices. Rank deficiency is
/// allowed: the positive shift and full residual, not rank claims, prove the
/// bound. Both halves of a nonsymmetric stored matrix are checked separately.
///
/// # Errors
/// Invalid matrix/shift/factor inputs and cancellation refuse. Work/storage
/// exhaustion, nonrepresentable arithmetic, or a large residual return an
/// explicit no-certificate result. No matrix or supplied factor is changed.
pub fn certify_shifted_gram(
    matrix: &Csr, shift: f64, columns: &[GramColumn], limits: SpectralInverseLimits,
    checkpoint: impl FnMut() -> bool,
) -> Result<SpectralPreparation, SpectralError> {
    let mut work = meter(limits, checkpoint);
    let result = (|| {
        admit(matrix, shift, &mut work)?;
        if columns.len() > matrix.nrows() {
            return Err(Halt::Fault(SpectralError::InvalidProposal("too many Gram columns")));
        }
        work.grow(matrix.nrows().checked_mul(8).ok_or(Halt::Storage)?)?;
        for column in columns {
            work.tick()?;
            work.grow(column.len())?;
            let mut previous = None;
            for &(row, value) in column {
                work.tick()?;
                if row >= matrix.nrows() || previous.is_some_and(|old| old >= row) || !value.is_finite() {
                    return Err(Halt::Fault(SpectralError::InvalidProposal(
                        "Gram columns require finite values at distinct increasing in-range rows",
                    )));
                }
                previous = Some(row);
            }
        }
        verify(matrix, shift, columns, &mut work)
    })();
    finish(work, 1, result)
}

/// Propose and independently verify a sparse shifted Cholesky factor.
/// The positive shift is a cost/conditioning proposal, NEVER a claimed
/// eigenvalue bound. Negative numerical pivots cause a smaller-shift retry;
/// scalar/edge work is shared across ALL attempts. Sparse fill and residual
/// expansion share the live-storage allowance. At most 256 charged visits
/// separate polls; ordered-map operations and the final bounded CSR clone
/// are opaque allocations/copies bracketed by polls.
///
/// # Errors
/// The same input/cancellation refusals as [`certify_shifted_gram`]. A failed
/// sufficient condition never asserts that the matrix is singular or indefinite.
pub fn prepare_spectral_inverse(
    matrix: &Csr, initial_shift: f64, limits: SpectralInverseLimits,
    checkpoint: impl FnMut() -> bool,
) -> Result<SpectralPreparation, SpectralError> {
    let mut work = meter(limits, checkpoint);
    let mut attempts = 0;
    let result = (|| {
        admit(matrix, initial_shift, &mut work)?;
        let mut shift = initial_shift;
        for _ in 0..limits.max_shift_attempts {
            attempts += 1;
            work.poll()?;
            // The preceding attempt's local maps/columns have been dropped.
            work.live = 0;
            work.grow(matrix.nrows().checked_mul(8).ok_or(Halt::Storage)?)?;
            if let Some(columns) = factor(matrix, shift, &mut work)? {
                if let Some(certificate) = verify(matrix, shift, &columns, &mut work)? {
                    return Ok(Some(certificate));
                }
            }
            shift *= 0.5;
            if shift == 0.0 { break; }
        }
        Ok(None)
    })();
    finish(work, attempts, result)
}

fn symmetric_add<F: FnMut() -> bool>(
    rows: &mut [BTreeMap<usize, f64>], i: usize, j: usize, value: f64,
    work: &mut Meter<F>,
) -> Result<(), Halt> {
    if !rows[i].contains_key(&j) { work.grow(2)?; }
    rows[i].insert(j, value);
    rows[j].insert(i, value);
    Ok(())
}

// Elimination uses the rounded symmetric part only as a factor PROPOSAL.
// Verification below always compares against the untouched ORIGINAL A.
fn factor<F: FnMut() -> bool>(
    a: &Csr, shift: f64, work: &mut Meter<F>,
) -> Result<Option<Vec<GramColumn>>, Halt> {
    let n = a.nrows();
    let mut rows = reserved(n)?;
    let mut diagonal = reserved(n)?;
    let mut columns = reserved(n)?;
    for _ in 0..n { work.tick()?; rows.push(BTreeMap::new()); diagonal.push(-shift); }
    for i in 0..n {
        for (&j, &value) in a.row(i).0.iter().zip(a.row(i).1) {
            work.tick()?;
            if i == j { diagonal[i] = finite(diagonal[i] + value)?; }
            else if value != 0.0 {
                let next = finite(value.mul_add(0.5, rows[i].get(&j).copied().unwrap_or(0.0)))?;
                symmetric_add(&mut rows, i, j, next, work)?;
            }
        }
    }
    let mut order = BTreeSet::new();
    for (i, row) in rows.iter().enumerate() { work.tick()?; order.insert((row.len(), i)); }
    while let Some((_, pivot)) = order.pop_first() {
        work.tick()?;
        let d = diagonal[pivot];
        if d <= 0.0 { return Ok(None); }
        let root = finite(fs_math::det::sqrt(d))?;
        if root <= 0.0 { return Ok(None); }
        let neighbors = std::mem::take(&mut rows[pivot]);
        // Detached neighbor records remain live until this pivot completes.
        work.grow(neighbors.len().checked_add(1).ok_or(Halt::Storage)?)?;
        let mut column = reserved(neighbors.len() + 1)?;
        let mut inserted_pivot = false;
        for (&i, &value) in &neighbors {
            work.tick()?;
            if !inserted_pivot && pivot < i { column.push((pivot, root)); inserted_pivot = true; }
            column.push((i, finite(value / root)?));
            order.remove(&(rows[i].len(), i));
            if rows[i].remove(&pivot).is_some() { work.release(1); }
        }
        if !inserted_pivot { column.push((pivot, root)); }
        for (slot, &(i, left)) in column.iter().enumerate() {
            if i == pivot { continue; }
            work.tick()?;
            diagonal[i] = finite((-left).mul_add(left, diagonal[i]))?;
            for &(j, right) in &column[slot + 1..] {
                if j == pivot { continue; }
                work.tick()?;
                let next = finite((-left).mul_add(right, rows[i].get(&j).copied().unwrap_or(0.0)))?;
                symmetric_add(&mut rows, i, j, next, work)?;
            }
        }
        for &i in neighbors.keys() { work.tick()?; order.insert((rows[i].len(), i)); }
        work.release(neighbors.len());
        columns.push(column);
    }
    work.poll()?;
    Ok(Some(columns))
}

fn residual_add<F: FnMut() -> bool>(
    rows: &mut [BTreeMap<usize, ScalarEnclosure>], i: usize, j: usize,
    value: f64, scale: f64, work: &mut Meter<F>,
) -> Result<(), Halt> {
    work.tick()?;
    if !rows[i].contains_key(&j) { work.grow(1)?; }
    let entry = rows[i].entry(j).or_insert_with(|| ScalarEnclosure::point(0.0));
    *entry = entry.add_scaled(ScalarEnclosure::point(value), scale)?;
    Ok(())
}
fn verify<F: FnMut() -> bool>(
    a: &Csr, shift: f64, columns: &[GramColumn], work: &mut Meter<F>,
) -> Result<Option<SpectralInverse>, Halt> {
    work.poll()?;
    let n = a.nrows();
    let mut residual = reserved(n)?;
    let mut column_sums = reserved(n)?;
    for _ in 0..n { work.tick()?; residual.push(BTreeMap::new()); column_sums.push(0.0); }
    for i in 0..n {
        for (&j, &value) in a.row(i).0.iter().zip(a.row(i).1) {
            work.tick()?;
            // Csr owns canonical unique coordinates; explicit zero entries
            // remain represented, and absent fill is independently introduced.
            work.grow(1)?;
            residual[i].insert(j, ScalarEnclosure::point(value));
        }
        residual_add(&mut residual, i, i, shift, -1.0, work)?;
    }
    for column in columns {
        work.tick()?;
        for (slot, &(i, left)) in column.iter().enumerate() {
            work.tick()?;
            for &(j, right) in &column[slot..] {
                residual_add(&mut residual, i, j, right, -left, work)?;
                if i != j { residual_add(&mut residual, j, i, right, -left, work)?; }
            }
        }
    }
    let mut row_norm = 0.0_f64;
    for row in &residual {
        work.tick()?;
        let mut sum = 0.0;
        for (&j, &entry) in row {
            work.tick()?;
            let absolute = entry.magnitude_upper();
            sum = add_up(sum, absolute)?;
            column_sums[j] = add_up(column_sums[j], absolute)?;
        }
        row_norm = row_norm.max(sum);
    }
    let mut defect = row_norm;
    for sum in column_sums { work.tick()?; defect = defect.max(sum); }
    let alpha = down(shift - defect)?;
    if alpha <= 0.0 { return Ok(None); }
    // Upward conversion handles n beyond binary64's exact integer range too.
    let n_upper = up(n as f64)?;
    let mut sqrt_n = up(fs_math::det::sqrt(n_upper))?;
    // Do not grant authority merely because a square-root backend returned
    // a plausible value: verify its squared lower endpoint against n.
    for _ in 0..4 {
        work.tick()?;
        if down(sqrt_n * sqrt_n)? >= n_upper { break; }
        sqrt_n = up(sqrt_n)?;
    }
    if down(sqrt_n * sqrt_n)? < n_upper { return Err(Halt::Range); }
    let inverse = up(sqrt_n / alpha)?;
    // Release residual records before retaining the exact CSR. Copy authority
    // and matrix together; later goal/feedback calls cannot substitute A.
    let mut released = 0;
    for row in &residual { work.tick()?; released += row.len(); }
    drop(residual);
    work.release(released);
    work.grow(a.nnz().checked_mul(2).and_then(|v| n.checked_add(1).and_then(|rows| v.checked_add(rows))).ok_or(Halt::Storage)?)?;
    work.poll()?;
    let matrix = a.clone();
    work.poll()?;
    Ok(Some(SpectralInverse { matrix, shift, defect_upper: defect,
        coercivity_lower: alpha, inverse_infinity_upper: inverse }))
}

#[cfg(test)]
mod tests;
