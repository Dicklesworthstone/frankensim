//! Bounded sparse interval LDL verification, without dense inverse columns.
//!
//! Reflect the stored lower triangle into H. Interval elimination encloses
//! H = L D L^T, and positive pivot lower endpoints establish its SPD property.
//! Two nonnegative triangular sweeps bound |H^-1| 1. The original matrix A
//! need not be exactly symmetric: its independently outward-bounded A-H
//! perturbation must also pass beta*||A-H||_inf < 1 before any A inverse bound
//! is published. This is stored-matrix authority only, not continuum error.

use std::collections::{BTreeMap, BTreeSet};
use fs_sparse::Csr;
use super::super::{
    GoalBoundStatus, GoalResidualError, GoalResidualLimits, GoalResidualReport,
    ScalarEnclosure, Work, enclose_goal_error,
};
use super::super::arithmetic::{add_up, down, mul_up, up};

/// Independent input, sparse-storage and elimination limits. Entry counts
/// are structural bounds, NOT allocator-byte or elapsed-time measurements.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SparseInverseLimits {
    /// Largest admitted square matrix dimension; no fixed 256-row ceiling.
    pub max_rows: usize,
    /// Largest original CSR storage, also retained for exact identity checks.
    pub max_input_nonzeros: usize,
    /// Active symmetric entries plus retained off-diagonal L entries.
    /// O(rows) diagonals, scratch and ordering nodes are additional.
    pub max_entries: usize,
    /// Total interval Schur-entry updates, including diagonal updates.
    pub max_updates: usize,
}

/// An inconclusive factor check is not a claim that the matrix is singular.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SparseInverseStatus {
    /// Positive interval pivots and the original-matrix perturbation passed.
    Enclosed,
    /// A pivot interval could not establish a strictly positive pivot.
    PivotNotPositive,
    /// The declared active/factor entry cap stopped further elimination.
    EntryLimit,
    /// The declared interval-update cap stopped further elimination.
    UpdateLimit,
    /// An outward arithmetic endpoint exceeded finite binary64 range.
    ArithmeticRange,
    /// The A-H perturbation was too large for the sufficient inverse check.
    PerturbationNotEstablished,
}

/// Work actually completed by one attempt, including unsuccessful attempts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SparseInverseWork {
    /// Completely eliminated columns.
    pub eliminated_rows: usize,
    /// Off-diagonal factor entries stored or reserved so far.
    pub factor_entries: usize,
    /// Completed interval Schur-entry updates.
    pub updates: usize,
    /// Peak active symmetric entries plus retained factor entries.
    pub peak_entries: usize,
}

/// Immutable inverse evidence tied to the EXACT stored CSR entries. Preparing
/// it may return an unavailable bound; inspect `status`/`inverse_upper`.
/// Reuse never repeats elimination, and a different matrix cannot borrow its
/// bound. The input snapshot costs O(nnz+rows); factor storage is released.
#[derive(Debug, Clone)]
pub struct VerifiedSparseInverse {
    matrix: Csr,
    inverse_upper: Option<f64>,
    perturbation_upper: Option<f64>,
    status: SparseInverseStatus,
    work: SparseInverseWork,
}

impl VerifiedSparseInverse {
    /// Factor with deterministic minimum-current-degree ordering, breaking
    /// ties by original row index. All interval updates and sparse-map visits
    /// are cancellable. Resource exhaustion during elimination retains only
    /// an explicit no-bound result and completed work, never a partial proof.
    ///
    /// # Errors
    /// Malformed/nonfinite input, input limits, failed scratch allocation or
    /// cancellation. Numerical/entry/update failures have typed statuses.
    pub fn prepare(
        matrix: &Csr, limits: SparseInverseLimits, checkpoint: impl FnMut() -> bool,
    ) -> Result<Self, GoalResidualError> {
        let mut work = Work { checkpoint, left: 512 };
        work.poll()?;
        let n = matrix.nrows();
        if n == 0 || matrix.ncols() != n {
            return Err(GoalResidualError::Shape { rows: n, columns: matrix.ncols() });
        }
        for (field, required, allowed) in [
            ("sparse inverse rows", n, limits.max_rows),
            ("sparse inverse input entries", matrix.nnz(), limits.max_input_nonzeros),
        ] {
            if required > allowed { return Err(GoalResidualError::Limit { field, required, allowed }); }
        }
        for row in 0..n {
            work.tick()?;
            for &value in matrix.row(row).1 {
                work.tick()?;
                if !value.is_finite() {
                    return Err(GoalResidualError::NonFinite { field: "sparse inverse matrix", index: row });
                }
            }
        }
        let mut stats = SparseInverseWork { eliminated_rows: 0, factor_entries: 0,
            updates: 0, peak_entries: 0 };
        let mut perturbation = None;
        let result = factor_bound(matrix, limits, &mut stats, &mut perturbation, &mut work);
        let (inverse_upper, status) = match result {
            Ok(Ok(bound)) => (Some(bound), SparseInverseStatus::Enclosed),
            Ok(Err(status)) => (None, status),
            Err(GoalResidualError::ArithmeticRange) => (None, SparseInverseStatus::ArithmeticRange),
            Err(error) => return Err(error),
        };
        // Copy with bounded visits/fallible vector reservation, rather than an
        // uninterruptible CSR clone after numerical work has been accepted.
        let snapshot = copy_matrix(matrix, &mut work)?;
        work.poll()?;
        Ok(Self { matrix: snapshot, inverse_upper, perturbation_upper: perturbation,
            status, work: stats })
    }

    /// Verified ||A^-1||_inf, not a supplied constant or a numerical estimate.
    #[must_use]
    pub const fn inverse_upper(&self) -> Option<f64> { self.inverse_upper }
    /// Why this attempt did or did not establish an inverse bound.
    #[must_use]
    pub const fn status(&self) -> SparseInverseStatus { self.status }
    /// Actual preparation work; subsequent enclosures perform no factorization.
    #[must_use]
    pub const fn work(&self) -> SparseInverseWork { self.work }
    /// Outward ||A-H||_inf, or None if preparation stopped before its check.
    #[must_use]
    pub const fn perturbation_upper(&self) -> Option<f64> { self.perturbation_upper }

    /// Enclose another goal/field on the same matrix, including the actual
    /// transposed dual defect. An existing tighter dominance bound is kept.
    /// The CSR identity comparison is exact (including signed zeros), with
    /// bounded checkpoints; no hash or caller assertion supplies authority.
    ///
    /// # Errors
    /// The ordinary residual refusals, a changed matrix, or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn enclose_goal(
        &self, matrix: &Csr, rhs: &[f64], primal: &[f64], goal: &[f64], dual: &[f64],
        scaling: Option<&[f64]>, limits: GoalResidualLimits,
        mut checkpoint: impl FnMut() -> bool,
    ) -> Result<GoalResidualReport, GoalResidualError> {
        let mut report = enclose_goal_error(matrix, rhs, primal, goal, dual, scaling,
            limits, &mut checkpoint)?;
        let mut work = Work { checkpoint, left: 512 };
        work.poll()?;
        if matrix.nrows() != self.matrix.nrows() || matrix.ncols() != self.matrix.ncols()
            || matrix.nnz() != self.matrix.nnz()
        { return Err(changed_matrix()); }
        for row in 0..matrix.nrows() {
            work.tick()?;
            let (indices, values) = matrix.row(row);
            let (expected_indices, expected_values) = self.matrix.row(row);
            if indices.len() != expected_indices.len() { return Err(changed_matrix()); }
            for ((&j, &a), (&k, &b)) in indices.iter().zip(values)
                .zip(expected_indices.iter().zip(expected_values))
            {
                work.tick()?;
                if j != k || a.to_bits() != b.to_bits() { return Err(changed_matrix()); }
            }
        }
        if let Some(bound) = self.inverse_upper
            && report.inverse_infinity_upper.is_none_or(|old| bound < old)
        {
            let enclosed = (|| {
                let radius = mul_up(mul_up(report.dual_residual_one_upper,
                    report.primal_residual_infinity_upper)?, bound)?;
                Ok::<_, GoalResidualError>((radius, report.weighted_residual.widen(radius)?))
            })();
            match enclosed
            {
                Ok((radius, error)) => {
                    report.inverse_infinity_upper = Some(bound);
                    report.dual_error_upper = Some(radius);
                    report.goal_error = Some(error);
                    report.status = GoalBoundStatus::Enclosed;
                }
                Err(GoalResidualError::ArithmeticRange) => {
                    // Preserve any independently established existing result.
                    if report.goal_error.is_none() { report.status = GoalBoundStatus::BoundNotRepresentable; }
                }
                Err(error) => return Err(error),
            }
        }
        work.poll()?;
        Ok(report)
    }
}

fn changed_matrix() -> GoalResidualError {
    GoalResidualError::MatrixMismatch
}
fn reserved<T>(count: usize) -> Result<Vec<T>, GoalResidualError> {
    let mut out = Vec::new();
    out.try_reserve_exact(count).map_err(|_| GoalResidualError::Allocation)?;
    Ok(out)
}
fn copy_matrix<F: FnMut() -> bool>(a: &Csr, work: &mut Work<F>) -> Result<Csr, GoalResidualError> {
    let mut ptr = reserved(a.nrows().checked_add(1).ok_or(GoalResidualError::Allocation)?)?;
    let mut indices = reserved(a.nnz())?;
    let mut values = reserved(a.nnz())?;
    ptr.push(0);
    for i in 0..a.nrows() {
        work.tick()?;
        let (columns, entries) = a.row(i);
        for (&j, &value) in columns.iter().zip(entries) {
            work.tick()?; indices.push(j); values.push(value);
        }
        ptr.push(indices.len());
    }
    Csr::try_from_parts_with_checkpoint(a.nrows(), a.ncols(), ptr, indices, values, || work.tick())?
        .ok_or(GoalResidualError::Allocation)
}
fn product(a: ScalarEnclosure, b: ScalarEnclosure) -> Result<ScalarEnclosure, GoalResidualError> {
    if (a.lower == 0.0 && a.upper == 0.0) || (b.lower == 0.0 && b.upper == 0.0) {
        return Ok(ScalarEnclosure::point(0.0));
    }
    let values = [a.lower*b.lower, a.lower*b.upper, a.upper*b.lower, a.upper*b.upper];
    let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Ok(ScalarEnclosure { lower: down(lo)?, upper: up(hi)? })
}
fn divide(a: ScalarEnclosure, d: ScalarEnclosure) -> Result<ScalarEnclosure, GoalResidualError> {
    // Called only after the positive-pivot gate; all endpoints are finite.
    if a.lower == 0.0 && a.upper == 0.0 { return Ok(a); }
    let values = [a.lower/d.lower, a.lower/d.upper, a.upper/d.lower, a.upper/d.upper];
    let lo = values.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    Ok(ScalarEnclosure { lower: down(lo)?, upper: up(hi)? })
}

type Attempt = Result<Result<f64, SparseInverseStatus>, GoalResidualError>;
#[allow(clippy::too_many_lines)]
fn factor_bound<F: FnMut() -> bool>(
    a: &Csr, limits: SparseInverseLimits, stats: &mut SparseInverseWork,
    perturbation: &mut Option<f64>, work: &mut Work<F>,
) -> Attempt {
    let n = a.nrows();
    let mut rows: Vec<BTreeMap<usize, ScalarEnclosure>> = reserved(n)?;
    let mut factors: Vec<Vec<(usize, f64)>> = reserved(n)?;
    let mut forward = reserved(n)?;
    let mut backward = reserved(n)?;
    let mut pivots = reserved(n)?;
    let mut difference = reserved(n)?;
    let mut order = reserved(n)?;
    for _ in 0..n {
        work.tick()?;
        rows.push(BTreeMap::new()); factors.push(Vec::new());
        forward.push(1.0); backward.push(0.0); pivots.push(0.0); difference.push(0.0);
    }
    let mut active = 0_usize;
    for i in 0..n {
        work.tick()?;
        if active == limits.max_entries { return Ok(Err(SparseInverseStatus::EntryLimit)); }
        rows[i].insert(i, ScalarEnclosure::point(a.get(i,i))); active += 1;
        let (columns, values) = a.row(i);
        for (&j, &value) in columns.iter().zip(values) {
            work.tick()?;
            if j < i {
                let upper = a.get(j,i);
                if upper != value { difference[j] = add_up(difference[j], up((upper-value).abs())?)?; }
                if value != 0.0 {
                    if limits.max_entries.saturating_sub(active) < 2 {
                        stats.peak_entries = active;
                        return Ok(Err(SparseInverseStatus::EntryLimit));
                    }
                    rows[i].insert(j, ScalarEnclosure::point(value));
                    rows[j].insert(i, ScalarEnclosure::point(value)); active += 2;
                }
            } else if j > i && a.row(j).0.binary_search(&i).is_err() {
                // An upper-only entry belongs to A-H even though H has no edge.
                difference[i] = add_up(difference[i], value.abs())?;
            }
        }
        stats.peak_entries = active;
    }
    let mut asymmetry = 0.0_f64;
    for value in difference { work.tick()?; asymmetry = asymmetry.max(value); }
    *perturbation = Some(asymmetry);
    let mut degrees = BTreeSet::new();
    for (i, row) in rows.iter().enumerate() { work.tick()?; degrees.insert((row.len()-1, i)); }
    while let Some((_, k)) = degrees.pop_first() {
        work.tick()?;
        let diagonal = rows[k][&k];
        if diagonal.lower <= 0.0 { return Ok(Err(SparseInverseStatus::PivotNotPositive)); }
        let m = rows[k].len()-1;
        let needed = m.checked_add(1).and_then(|v| m.checked_mul(v)).map(|v| v/2);
        if needed.is_none_or(|v| v > limits.max_updates.saturating_sub(stats.updates)) {
            return Ok(Err(SparseInverseStatus::UpdateLimit));
        }
        if m > limits.max_entries.saturating_sub(active.saturating_add(stats.factor_entries)) {
            return Ok(Err(SparseInverseStatus::EntryLimit));
        }
        let mut neighbors = reserved(m)?;
        let mut column = reserved(m)?;
        for (&i, &value) in &rows[k] {
            work.tick()?;
            if i == k { continue; }
            degrees.remove(&(rows[i].len()-1, i));
            let multiplier = divide(value, diagonal)?;
            neighbors.push((i, value, multiplier));
            column.push((i, multiplier.magnitude_upper()));
        }
        stats.factor_entries += m;
        stats.peak_entries = stats.peak_entries.max(active+stats.factor_entries);
        for (r, &(i, _, multiplier)) in neighbors.iter().enumerate() {
            work.tick()?;
            for &(j, value, _) in &neighbors[..=r] {
                work.tick()?;
                let old = rows[i].get(&j).copied();
                if old.is_none() && limits.max_entries
                    .saturating_sub(active+stats.factor_entries) < 2
                { return Ok(Err(SparseInverseStatus::EntryLimit)); }
                let updated = old.unwrap_or(ScalarEnclosure::point(0.0))
                    .add_scaled(product(multiplier, value)?, -1.0)?;
                rows[i].insert(j, updated);
                if i != j { rows[j].insert(i, updated); }
                if old.is_none() { active += 2; }
                stats.updates += 1;
                stats.peak_entries = stats.peak_entries.max(active+stats.factor_entries);
            }
        }
        for &(i, multiplier) in &column {
            work.tick()?;
            forward[i] = add_up(forward[i], mul_up(multiplier, forward[k])?)?;
            rows[i].remove(&k); active -= 1;
            degrees.insert((rows[i].len()-1, i));
        }
        active -= rows[k].len();
        while rows[k].pop_first().is_some() { work.tick()?; }
        factors[k] = column; pivots[k] = diagonal.lower; order.push(k);
        stats.eliminated_rows += 1;
    }
    let mut beta = 0.0_f64;
    for &k in order.iter().rev() {
        work.tick()?;
        let mut bound = up(forward[k]/pivots[k])?;
        for &(i, multiplier) in &factors[k] {
            work.tick()?;
            bound = add_up(bound, mul_up(multiplier, backward[i])?)?;
        }
        backward[k] = bound; beta = beta.max(bound);
    }
    if asymmetry == 0.0 { return Ok(Ok(beta)); }
    let gain = mul_up(beta, asymmetry)?;
    if gain >= 1.0 { return Ok(Err(SparseInverseStatus::PerturbationNotEstablished)); }
    let denominator = down(1.0-gain)?;
    if denominator <= 0.0 { return Ok(Err(SparseInverseStatus::PerturbationNotEstablished)); }
    Ok(Ok(up(beta/denominator)?))
}
