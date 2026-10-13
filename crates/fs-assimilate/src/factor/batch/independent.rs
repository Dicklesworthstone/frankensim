//! Structural operations on already admitted UD factors.

use crate::factor::{FactorBelief, packed_upper_index, packed_upper_len, try_f64_vec};
use crate::{
    AssimError, Belief, Cx, WorkPlan, WorkProgress, belief_validation_work_plan,
    canonicalize_belief_zeros, validate_belief_parts, validate_state_dimension,
};

impl FactorBelief {
    /// Form the joint belief of two DECLARED INDEPENDENT factor beliefs.
    /// Ordering is all of `self`, then all of `other`. Both diagonal blocks
    /// are copied exactly; cross-block entries are exactly zero. No enlarged
    /// dense PSD test or refactorization can erase an admitted null direction.
    /// Independence is a caller assumption, not inferred from marginal data.
    /// Do not use this to rejoin marginals of a previously correlated state.
    ///
    /// # Errors
    /// Refuses unsupported combined dimensions, allocation failure, exhausted
    /// budgets or cancellation. Inputs are borrowed and never changed.
    pub fn independent_product(&self, other: &Self, cx: &Cx<'_>) -> Result<Self, AssimError> {
        let left = self.dim();
        let n = left.checked_add(other.dim()).ok_or(AssimError::WorkPlanOverflow {
            phase: "independent factor shape",
        })?;
        validate_state_dimension(n)?;
        let nu = n as u128;
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(
            0, 0, 0, 0, 4 * nu * nu + 16 * nu + 16, 0, 0,
        )?)?;
        progress.checkpoint("independent factor admission")?;
        let mut mean = try_f64_vec(n, "independent factor mean", &mut progress)?;
        let mut diag = try_f64_vec(n, "independent factor diagonal", &mut progress)?;
        let mut upper = try_f64_vec(packed_upper_len(n)?, "independent factor upper", &mut progress)?;
        for (offset, block) in [(0, self), (left, other)] {
            for row in 0..block.dim() {
                mean[offset + row] = block.mean[row];
                diag[offset + row] = block.diag[row];
                progress.scalar("independent factor copy", 2)?;
                for column in row + 1..block.dim() {
                    upper[packed_upper_index(n, offset + row, offset + column)] =
                        block.upper[packed_upper_index(block.dim(), row, column)];
                    progress.scalar("independent factor copy", 1)?;
                }
            }
        }
        progress.checkpoint("independent factor publication")?;
        Ok(Self { mean, upper, diag })
    }

    /// Publish the checked dense marginal of the first `count` components.
    /// ALL factor columns contribute, including learned cross-block effects.
    /// Only the requested marginal is PSD-checked; a singular excluded nuisance
    /// block cannot prevent publishing an otherwise resolvable physical block.
    /// The retained joint factor remains unchanged. Cost is O(count² * dim)
    /// for reconstruction plus the existing marginal PSD admission.
    ///
    /// # Errors
    /// Refuses invalid dimensions, nonfinite/unresolved marginal covariance,
    /// allocation failure, cancellation or exhausted cumulative stage budgets.
    pub fn prefix_marginal(&self, count: usize, cx: &Cx<'_>) -> Result<Belief, AssimError> {
        validate_state_dimension(self.dim())?;
        validate_state_dimension(count)?;
        if count > self.dim() {
            return Err(AssimError::DimMismatch { state: self.dim(), operator: count });
        }
        let (nu, du) = (count as u128, self.dim() as u128);
        let validation = belief_validation_work_plan(count)?.total;
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(
            validation, 0, 0, 0, 4 * nu * nu * du + 8 * nu * nu + 16 * nu, 0, 0,
        )?)?;
        progress.checkpoint("factor marginal admission")?;
        let mut mean = try_f64_vec(count, "factor marginal mean", &mut progress)?;
        mean.copy_from_slice(&self.mean[..count]);
        progress.scalar("factor marginal mean", nu)?;
        let mut cov = Vec::new();
        cov.try_reserve_exact(count).map_err(|_| AssimError::AllocationRefused {
            stage: "factor marginal rows",
            requested_bytes: nu * core::mem::size_of::<Vec<f64>>() as u128,
        })?;
        for _ in 0..count {
            cov.push(try_f64_vec(count, "factor marginal row", &mut progress)?);
        }
        for row in 0..count {
            for column in row..count {
                let mut entry = 0.0;
                for k in column..self.dim() {
                    let left = if row == k { 1.0 } else { self.upper[packed_upper_index(self.dim(), row, k)] };
                    let right = if column == k { 1.0 } else { self.upper[packed_upper_index(self.dim(), column, k)] };
                    entry += left * self.diag[k] * right;
                    progress.scalar("factor marginal reconstruction", 3)?;
                }
                cov[row][column] = entry;
                cov[column][row] = entry;
            }
        }
        validate_belief_parts(&mean, &cov, "factor marginal validation", &mut progress)?;
        canonicalize_belief_zeros(&mut mean, &mut cov, &mut progress)?;
        progress.checkpoint("factor marginal publication")?;
        Ok(Belief { mean, cov })
    }
}
