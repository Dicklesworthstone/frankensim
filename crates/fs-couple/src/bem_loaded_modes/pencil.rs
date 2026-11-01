//! Shared fixed balancing and original-coordinate audits for loaded boundaries.
//! Geometry-specific callers still assemble their own physical equations.
use super::*;

impl PhysicalPencil {
    /// Fixed two-sided scales from one declared reference operator. Frequency
    /// callbacks only APPLY these scales; they never rebalance a contour node.
    pub(super) fn fixed_scales(
        &self,
        n: usize,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<(Vec<f64>, Vec<f64>), BemModesError> {
        let mut row_scales = vec![1.0; n];
        let mut column_scales = vec![1.0; n];
        for i in 0..n {
            let row: f64 = self.envelope[i * n..(i + 1) * n].iter().sum();
            if !row.is_finite() || row <= 0.0 {
                return Err(bad("zero or unrepresentable reference row"));
            }
            row_scales[i] = row;
        }
        for j in 0..n {
            let col = (0..n)
                .map(|i| self.envelope[i * n + j] / row_scales[i])
                .fold(0.0, f64::max);
            let scale = 1.0 / col;
            if !scale.is_finite() || scale <= 0.0 {
                return Err(bad("zero or unrepresentable reference coordinate"));
            }
            column_scales[j] = scale;
        }
        for i in 0..n {
            let row: f64 = (0..n)
                .map(|j| self.envelope[i * n + j] * column_scales[j])
                .sum();
            if !row.is_finite() || row <= 0.0 {
                return Err(bad("unrepresentable balanced reference row"));
            }
            row_scales[i] = row;
        }
        poll(cancelled)?;
        Ok((row_scales, column_scales))
    }

    pub(super) fn into_balanced(
        mut self,
        row_scales: &[f64],
        column_scales: &[f64],
    ) -> Result<PoleEvaluation, BemModesError> {
        let n = row_scales.len();
        let mut scale = 0.0_f64;
        for i in 0..n {
            let mut row = 0.0;
            for j in 0..n {
                let k = i * n + j;
                let factor = column_scales[j] / row_scales[i];
                self.value.matrix[k] = self.value.matrix[k].scale(factor);
                self.value.derivative[k] = self.value.derivative[k].scale(factor);
                row += self.envelope[k] * factor;
            }
            scale = scale.max(row);
        }
        if !scale.is_finite()
            || scale <= 0.0
            || self
                .value
                .matrix
                .iter()
                .chain(&self.value.derivative)
                .any(|&z| !finite(z))
        {
            return Err(bad("balanced coupled operator is not representable"));
        }
        self.value.residual_scale = scale;
        Ok(self.value)
    }

    /// Shared compensated audit of original mechanical, boundary and memory
    /// equations. A source RHS acts on the first mechanical block only.
    ///
    /// Each row is judged by two clauses, the first that holds accepting it:
    ///   1. componentwise: `|sum_i| <= tolerance * (|rhs_i| + sum_j E_ij |v_j|)`, the row against its own terms;
    ///   2. rounding floor: `|sum_i| <= gamma_{n+3} * (|rhs_i| + sum_j E_ij * max_j |v_j|)`, the row against the
    ///      rounding of evaluating it at the vector's norm.
    ///
    /// The second clause exists for rows whose own terms are all below the norm-scale rounding: a mode that is
    /// dark to a coordinate (zero coupling) leaves that coordinate's entries as iteration noise, and their row
    /// divides noise by noise. `gamma_N = N u / (1 - N u)`, `u = 2^-53`, and `N = n + 3` counts: `n` for the
    /// summation of the `n + 1` terms of the row (the `n` products `A_ij v_j` and `-rhs_i`; a recursive sum of
    /// `n + 1` terms carries `gamma_n` of the sum of their magnitudes), `2` for the roundings inside each complex
    /// product, and `1` for the rounded stored entry of `A` or `v` that enters it. The contraction below is
    /// compensated, so the real error is smaller; the floor takes the plain bound. A row with signal passes
    /// clause 1 by itself, and a perturbation that matters moves `sum_i` by the row's own terms times its
    /// relative size, far above `gamma_{n+3}`.
    ///
    /// The reported value of a row accepted by clause 2 is `|sum_i| * tolerance / floor_i` (at most `tolerance`);
    /// a refused row reports its componentwise ratio.
    pub(super) fn audit_physical_rows(
        &self,
        modes: usize,
        boundary: usize,
        vector: &[C64],
        force: Option<&[C64]>,
        tolerance: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<Vec<RowAudit>, BemModesError> {
        let n = vector.len();
        if n.checked_mul(n) != Some(self.value.matrix.len())
            || self.envelope.len() != self.value.matrix.len()
            || modes.checked_add(boundary).is_none_or(|v| v > n)
            || force.is_some_and(|f| f.len() != modes)
        {
            return Err(bad(
                "physical residual needs complete state and force coordinates",
            ));
        }
        let units = (n + 3) as f64;
        let u = f64::EPSILON / 2.0;
        let gamma = units * u / (1.0 - units * u);
        let norm = vector.iter().map(|v| v.abs()).fold(0.0, f64::max);
        let mut rows = Vec::with_capacity(n);
        for i in 0..n {
            poll(cancelled)?;
            let rhs = force.and_then(|f| f.get(i)).copied().unwrap_or(C64::ZERO);
            let mut sum = -rhs;
            let mut correction = C64::ZERO;
            let mut denominator = rhs.abs();
            let mut row_envelope = 0.0_f64;
            for j in 0..n {
                let term = self.value.matrix[i * n + j] * vector[j] - correction;
                let next = sum + term;
                correction = (next - sum) - term;
                sum = next;
                denominator += self.envelope[i * n + j] * vector[j].abs();
                row_envelope += self.envelope[i * n + j];
            }
            if !finite(sum) || !denominator.is_finite() || !row_envelope.is_finite() {
                return Err(bad("physical residual contraction overflow"));
            }
            let relative = if denominator == 0.0 {
                if sum.abs() == 0.0 { 0.0 } else { f64::INFINITY }
            } else {
                sum.abs() / denominator
            };
            if !relative.is_finite() {
                return Err(bad("physical residual is not representable"));
            }
            let floor = gamma * (rhs.abs() + row_envelope * norm);
            let (clause, reported) = if relative <= tolerance {
                (RowClause::Componentwise, relative)
            } else if sum.abs() <= floor {
                (RowClause::RoundingFloor, sum.abs() * tolerance / floor)
            } else {
                (RowClause::Refused, relative)
            };
            let block = if i < modes {
                0
            } else if i < modes + boundary {
                1
            } else {
                2
            };
            rows.push(RowAudit {
                row: i,
                block,
                residual_abs: sum.abs(),
                floor,
                componentwise: relative,
                reported,
                clause,
            });
        }
        poll(cancelled)?;
        Ok(rows)
    }

    /// Per-block maxima of [`Self::audit_physical_rows`]; any refused row refuses the vector.
    pub(super) fn physical_residuals(
        &self,
        modes: usize,
        boundary: usize,
        vector: &[C64],
        force: Option<&[C64]>,
        tolerance: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<[f64; 3], BemModesError> {
        let rows = self.audit_physical_rows(modes, boundary, vector, force, tolerance, cancelled)?;
        let mut residuals = [0.0_f64; 3];
        for row in &rows {
            residuals[row.block] = residuals[row.block].max(row.reported);
        }
        if residuals.iter().any(|&r| r > tolerance) {
            return Err(BemModesError::PhysicalResidual {
                mechanical: residuals[0],
                boundary: residuals[1],
                memory: residuals[2],
            });
        }
        Ok(residuals)
    }
}

/// Which clause of the row audit accepted a row.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum RowClause {
    /// The row against its own terms, at the declared tolerance.
    Componentwise,
    /// The row against the rounding of evaluating it at the vector's norm.
    RoundingFloor,
    /// Neither: the vector is refused.
    Refused,
}

/// One audited row: `|sum_i|`, its rounding floor, the componentwise ratio, the value reported, and the clause.
#[derive(Clone, Copy, Debug)]
pub(super) struct RowAudit {
    pub(super) row: usize,
    pub(super) block: usize,
    pub(super) residual_abs: f64,
    pub(super) floor: f64,
    pub(super) componentwise: f64,
    pub(super) reported: f64,
    pub(super) clause: RowClause,
}

/// Original physical amplitudes, normalized by maximum structural displacement
/// or (for a zero-motion mode) maximum remaining original coordinate.
pub(super) fn physical_eigenstate(
    pole: &NonlinearPole,
    column_scales: &[f64],
    modes: usize,
) -> Result<Vec<C64>, BemModesError> {
    let mut vector: Vec<_> = pole
        .vector
        .iter()
        .zip(column_scales)
        .map(|(&v, &s)| v.scale(s))
        .collect();
    let mut amplitude = vector[..modes].iter().map(|v| v.abs()).fold(0.0, f64::max);
    if amplitude == 0.0 {
        amplitude = vector.iter().map(|v| v.abs()).fold(0.0, f64::max);
    }
    if !amplitude.is_finite() || amplitude <= 0.0 {
        return Err(bad("original eigenstate scale is not representable"));
    }
    for v in &mut vector {
        *v = divide(*v, amplitude);
    }
    if vector.iter().any(|&v| !finite(v)) {
        return Err(bad("original eigenstate overflow"));
    }
    Ok(vector)
}
