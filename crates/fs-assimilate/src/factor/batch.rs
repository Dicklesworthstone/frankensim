//! Factor-retaining linear-Gaussian batches (sj31i.64).
//!
//! Whiten a correlated observation model once, then carry `U D U^T` through
//! every row. There is no dense Joseph update or dense PSD re-admission in
//! the row loop. Declaration order is deliberate: a whitening factor and
//! its rows form one ordered protocol, not an unordered set of sensors.

use super::{
    ContractionReceipt, FactorBelief, assimilate_scalar_planned, compensated_dot,
    factor_construction_work_plan, factor_update_work_plan, packed_upper_index, try_f64_vec,
};
use crate::{
    AssimError, Belief, Cx, MAX_DENSE_OBSERVATIONS, Observation, WorkPlan, WorkProgress,
    belief_validation_work_plan, canonicalize_belief_zeros, checked_work_add, checked_work_mul,
    validate_belief_parts, validate_state_dimension,
};

/// Diagnostics of one complete factor sequence, not a certificate of the
/// caller's statistical model or a general floating-point Loewner proof.
#[derive(Debug, Clone, PartialEq)]
pub struct FactorBatchDiagnostics {
    receipts: Vec<ContractionReceipt>,
    updated_row_indices: Vec<usize>,
    state_independent_row_indices: Vec<usize>,
    state_independent_misfit: f64,
    joint_misfit_before: f64,
    joint_misfit_after: f64,
}

impl FactorBatchDiagnostics {
    /// Scalar contraction evidence in sequence (whitened-row) order.
    #[must_use]
    pub fn receipts(&self) -> &[ContractionReceipt] {
        &self.receipts
    }

    /// Row indices aligned with `receipts()`. For a correlated batch these
    /// refer to the ordered whitened rows, not to individual physical sensors.
    #[must_use]
    pub fn updated_row_indices(&self) -> &[usize] {
        &self.updated_row_indices
    }

    /// Rows whose computed whitened operator is exactly zero. They add a
    /// constant likelihood residual but no state update. This numerical
    /// distinction is not an exact covariance-rank or information certificate.
    #[must_use]
    pub fn state_independent_row_indices(&self) -> &[usize] {
        &self.state_independent_row_indices
    }

    /// Constant squared whitened residual included in BOTH joint misfits.
    #[must_use]
    pub const fn state_independent_misfit(&self) -> f64 {
        self.state_independent_misfit
    }

    pub(crate) fn with_whitening_rows(
        mut self,
        updated_rows: Vec<usize>,
        independent_rows: Vec<usize>,
        constant_misfit: f64,
    ) -> Result<Self, AssimError> {
        debug_assert_eq!(updated_rows.len(), self.receipts.len());
        self.joint_misfit_before += constant_misfit;
        self.joint_misfit_after += constant_misfit;
        if !self.joint_misfit_before.is_finite() || !self.joint_misfit_after.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "joint whitened residual",
            });
        }
        self.updated_row_indices = updated_rows;
        self.state_independent_row_indices = independent_rows;
        self.state_independent_misfit = constant_misfit;
        Ok(self)
    }

    /// Full-batch weighted residual at the original prior mean.
    /// For whitened rows this is `(y - H mean)^T R^-1 (y - H mean)`.
    #[must_use]
    pub const fn joint_misfit_before(&self) -> f64 {
        self.joint_misfit_before
    }

    /// The SAME full-batch residual at the final posterior mean. This is not
    /// a sum of scalar receipts evaluated at different intermediate means.
    #[must_use]
    pub const fn joint_misfit_after(&self) -> f64 {
        self.joint_misfit_after
    }
}

/// A completed batch with its covariance still in factor form.
#[derive(Debug, Clone, PartialEq)]
pub struct FactorBatchAssimilation {
    belief: FactorBelief,
    diagnostics: FactorBatchDiagnostics,
}

impl FactorBatchAssimilation {
    /// Factor posterior, usable as the next batch's prior without refactoring.
    #[must_use]
    pub const fn belief(&self) -> &FactorBelief {
        &self.belief
    }

    /// Complete-batch diagnostics and per-row contraction evidence.
    #[must_use]
    pub const fn diagnostics(&self) -> &FactorBatchDiagnostics {
        &self.diagnostics
    }

    /// Move the posterior and diagnostics without copying the quadratic factor.
    #[must_use]
    pub fn into_parts(self) -> (FactorBelief, FactorBatchDiagnostics) {
        (self.belief, self.diagnostics)
    }
}

/// Fuse independent scalar rows in declaration order with one cumulative work
/// and poll budget. Correlated callers MUST whiten both values and operators
/// with their full observation covariance before calling this function.
///
/// Work is `O(m n^2)` with `O(n^2 + m)` retained storage (excluding the caller's
/// observations). No intermediate dense posterior is constructed. All input
/// dimensions are checked before updating; a failure publishes no partial
/// batch, and the borrowed prior is never mutated.
///
/// # Errors
/// Refuses empty/oversized batches, mismatched dimensions, non-finite
/// arithmetic, exhausted budgets, and cancellation. Per-row contraction
/// verdicts remain evidence, not statistical-model admission decisions.
pub fn assimilate_sequence(
    prior: &FactorBelief,
    observations: &[Observation],
    cx: &Cx<'_>,
) -> Result<FactorBatchAssimilation, AssimError> {
    if observations.is_empty() {
        return Err(AssimError::EmptyObservations);
    }
    if observations.len() > MAX_DENSE_OBSERVATIONS {
        return Err(AssimError::ObservationCountLimit {
            count: observations.len(),
            max: MAX_DENSE_OBSERVATIONS,
        });
    }
    let n = prior.dim();
    validate_state_dimension(n)?;
    let m = observations.len() as u128;
    let updates = checked_work_mul(factor_update_work_plan(n)?.total, m, "factor batch")?;
    let diagnostics = checked_work_mul(
        checked_work_add(16 * n as u128, 16, "factor batch misfit")?,
        m,
        "factor batch misfit",
    )?;
    let plan = WorkPlan::checked(0, m, 0, diagnostics, updates, 0, 0)?;
    let mut progress = WorkProgress::new(cx, plan)?;
    progress.checkpoint("factor batch admission")?;
    let mut updated_row_indices = Vec::new();
    updated_row_indices.try_reserve_exact(observations.len()).map_err(|_| {
        AssimError::AllocationRefused {
            stage: "factor batch row indices",
            requested_bytes: m * core::mem::size_of::<usize>() as u128,
        }
    })?;
    for (index, observation) in observations.iter().enumerate() {
        if observation.operator().len() != n {
            return Err(AssimError::DimMismatch {
                state: n,
                operator: observation.operator().len(),
            });
        }
        updated_row_indices.push(index);
        progress.records("factor batch dimensions", 1)?;
    }
    let mut receipts = Vec::new();
    receipts.try_reserve_exact(observations.len()).map_err(|_| {
        AssimError::AllocationRefused {
            stage: "factor batch receipts",
            requested_bytes: m * core::mem::size_of::<ContractionReceipt>() as u128,
        }
    })?;
    let joint_misfit_before = joint_misfit(prior.mean(), observations, &mut progress)?;
    let mut belief = prior.clone();
    for observation in observations {
        let update = assimilate_scalar_planned(&belief, observation, &mut progress)?;
        belief = update.belief;
        receipts.push(update.receipt);
    }
    let joint_misfit_after = joint_misfit(belief.mean(), observations, &mut progress)?;
    progress.checkpoint("factor batch publication")?;
    Ok(FactorBatchAssimilation {
        belief,
        diagnostics: FactorBatchDiagnostics {
            receipts,
            updated_row_indices,
            state_independent_row_indices: Vec::new(),
            state_independent_misfit: 0.0,
            joint_misfit_before,
            joint_misfit_after,
        },
    })
}

fn joint_misfit(
    mean: &[f64],
    observations: &[Observation],
    progress: &mut WorkProgress<'_, '_>,
) -> Result<f64, AssimError> {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for observation in observations {
        let (prediction, _) = compensated_dot(
            observation.operator(), mean, "batch prediction", "factor batch misfit", progress,
        )?;
        let residual = observation.value() - prediction;
        let scaled = residual / observation.noise_var().sqrt();
        let term = scaled * scaled;
        let next = sum + term;
        correction += if sum.abs() >= term.abs() {
            (sum - next) + term
        } else {
            (term - next) + sum
        };
        sum = next;
        progress.scalar("factor batch misfit", 8)?;
    }
    let total = sum + correction;
    if !total.is_finite() {
        return Err(AssimError::NonFiniteComputation { stage: "factor batch misfit" });
    }
    Ok(total)
}

impl FactorBelief {
    /// Reconstruct and validate ONE dense posterior at a publication boundary.
    /// Unlike `to_dense_covariance`, reconstruction polls cancellation at the
    /// fixed scalar stride. It and the ordinary PSD admission share one budget.
    /// This costs `O(n^3)` once, not once per observation.
    ///
    /// # Errors
    /// Refuses non-finite reconstruction, unresolved dense PSD admission,
    /// cancellation, allocation failure, or an exhausted budget. Keeping the
    /// factor result does not require this optional dense conversion.
    pub fn to_belief(&self, cx: &Cx<'_>) -> Result<Belief, AssimError> {
        let n = self.dim();
        validate_state_dimension(n)?;
        let construction = factor_construction_work_plan(n)?.total;
        let validation = belief_validation_work_plan(n)?.total;
        let plan = WorkPlan::checked(validation, 0, 0, 0, construction, 0, 0)?;
        let mut progress = WorkProgress::new(cx, plan)?;
        progress.checkpoint("factor dense publication")?;
        let mut cov = Vec::new();
        cov.try_reserve_exact(n).map_err(|_| AssimError::AllocationRefused {
            stage: "factor dense rows", requested_bytes: n as u128 * core::mem::size_of::<Vec<f64>>() as u128,
        })?;
        for _ in 0..n {
            cov.push(try_f64_vec(n, "factor dense row", &mut progress)?);
        }
        for i in 0..n {
            for j in i..n {
                let mut entry = 0.0;
                for k in j..n {
                    let u_ik = if k == i { 1.0 } else { self.upper[packed_upper_index(n, i, k)] };
                    let u_jk = if k == j { 1.0 } else { self.upper[packed_upper_index(n, j, k)] };
                    entry += u_ik * self.diag[k] * u_jk;
                    progress.scalar("factor dense reconstruction", 3)?;
                }
                cov[i][j] = entry;
                cov[j][i] = entry;
            }
        }
        let mut mean = self.mean.clone();
        validate_belief_parts(&mean, &cov, "factor dense validation", &mut progress)?;
        canonicalize_belief_zeros(&mut mean, &mut cov, &mut progress)?;
        progress.checkpoint("factor dense publication")?;
        Ok(Belief { mean, cov })
    }
}

mod independent;
