//! Mean-only controls in the physical coordinates of a Gaussian copula.
//!
//! An adjoint supplies derivatives in physical units, not latent-normal units.
//! Replay the exact existing copula transform before centering by the declared
//! uniform means. No physical solve, fit, new draw, or observation mutation is
//! performed by assessment. Numerical sampling/transform bias is not bounded.

use fs_blake3::{ContentHash, hash_domain};

use super::{GaussianCopulaExecution, physical_parameters};
use crate::product_checkpoint::execution_identity;
use crate::product_control_variate::{LinearControlEstimate, LinearControlVariate,
    UqControlError, assess_with_parameters, checkpoint};
use crate::product_execution::sample_parameters;
use crate::{CorrelationModel, UqCheckpointError, UqExecution, UqPlan, UqStatus};

/// Coefficients fixed before sampling and bound to the full physical copula.
///
/// Derivatives use the units/order returned by the execution's `marginals()`.
/// Expectations come from those marginals, never from the zero latent means.
/// Coefficients may increase variance; there is no favorable-result selection.
/// The caller must obtain them independently of this sample, for example from
/// one nominal adjoint. Freezing a fresh execution cannot prove that independence.
#[derive(Debug, Clone)]
pub struct CopulaLinearControlVariate {
    control: LinearControlVariate,
    sampling_identity: ContentHash,
}

impl GaussianCopulaExecution {
    fn control_identity(&self) -> ContentHash {
        let domain = hash_domain("org.frankensim.uq.copula-mean-control.v1", b"");
        execution_identity(self.execution.plan(), self.bound_model(domain))
    }

    /// Freeze physical-unit coefficients before any accepted or failed sample.
    /// The latent matrix, physical supports/units, seed and lifetime policy are
    /// all bound. No sampler draw or model evaluation is consumed.
    pub fn freeze_linear_control_variate(&self, gradient: &[f64])
        -> Result<CopulaLinearControlVariate, UqControlError>
    {
        if self.execution.status == UqStatus::Refused || self.execution.failure.is_some() {
            return Err(UqControlError::RefusedExecution);
        }
        if self.execution.attempted != 0 || !self.execution.values.is_empty() {
            return Err(UqControlError::AlreadySampled);
        }
        // Admission only: this independent plan supplies PHYSICAL expectations
        // and coefficient validation. It never draws an independent sample.
        let mut marginal_plan = self.execution.plan().clone();
        marginal_plan.parameters = self.marginals.clone();
        marginal_plan.correlation = CorrelationModel::Independent;
        let fresh = UqExecution::new(&marginal_plan).map_err(|_| UqControlError::PlanMismatch)?;
        Ok(CopulaLinearControlVariate {
            control: fresh.freeze_linear_control_variate(gradient)?,
            sampling_identity: self.control_identity(),
        })
    }

    /// Assess the fixed physical control on every accepted raw observation.
    /// Returns only a mean comparison, never adjusted quantiles or compliance.
    pub fn assess_linear_control_variate(&self, control: &CopulaLinearControlVariate)
        -> Result<Option<LinearControlEstimate>, UqControlError>
    {
        self.assess_linear_control_variate_interruptible(control, || false)
    }

    /// Bounded replay using the original latent draws AND physical transform.
    /// O(n*d) work and O(n+d) additional storage. Polls before each draw and at
    /// the existing statistics tiles. Refusal/cancellation returns no estimate
    /// and leaves the physical observations, reports and checkpoints untouched.
    /// Partial standard errors are descriptive, not optional-stopping bounds.
    pub fn assess_linear_control_variate_interruptible(
        &self, control: &CopulaLinearControlVariate, mut cancelled: impl FnMut() -> bool,
    ) -> Result<Option<LinearControlEstimate>, UqControlError> {
        if self.control_identity() != control.sampling_identity {
            return Err(UqControlError::PlanMismatch);
        }
        if self.execution.status == UqStatus::Refused || self.execution.failure.is_some() {
            return Err(UqControlError::RefusedExecution);
        }
        assess_with_parameters(&control.control, &self.execution.values, |ordinal| {
            let latent = sample_parameters(self.execution.plan(), self.execution.factor.as_deref(), ordinal)
                .map_err(|_| UqControlError::NumericalRange)?;
            physical_parameters(&self.marginals, &latent).map_err(|_| UqControlError::NumericalRange)
        }, &mut cancelled)
    }
}

impl CopulaLinearControlVariate {
    /// Analytic physical marginal means, in declaration order.
    #[must_use]
    pub fn parameter_means(&self) -> &[f64] { self.control.parameter_means() }

    /// Frozen physical-unit derivatives, not derivatives in latent coordinates.
    #[must_use]
    pub fn gradient(&self) -> &[f64] { self.control.gradient() }

    /// Save exact coefficient bits with the RAW copula observation prefix.
    /// The existing copula/raw checkpoint owners still enforce support, matrix,
    /// units, ordinal, status, checksum and model binding. Never checkpoint the
    /// read-only latent execution instead. Integrity is not authentication.
    pub fn checkpoint(&self, execution: &GaussianCopulaExecution, model: ContentHash)
        -> Result<Vec<u8>, UqCheckpointError>
    {
        if execution.control_identity() != self.sampling_identity {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let mut bytes = checkpoint::encode_header(self.gradient());
        let identity = checkpoint::bound_model(model, &bytes);
        bytes.extend_from_slice(&execution.checkpoint(identity)?);
        Ok(bytes)
    }

    /// Restore the original physical control without a nominal solve or refit.
    /// A different support, unit, latent matrix, model, policy or coefficient
    /// header refuses. Physical means are reconstructed, not read from the file.
    pub fn restore(plan: &UqPlan, matrix: &[Vec<f64>], model: ContentHash, bytes: &[u8])
        -> Result<(GaussianCopulaExecution, Self), UqCheckpointError>
    {
        let fresh = GaussianCopulaExecution::new(plan, matrix).map_err(UqCheckpointError::InvalidPlan)?;
        let (gradient, prefix) = checkpoint::decode_header(plan.parameters.len(), plan.budget_max_samples, bytes)?;
        let control = fresh.freeze_linear_control_variate(&gradient)
            .map_err(|_| UqCheckpointError::InvalidEncoding("invalid physical control coefficients"))?;
        let identity = checkpoint::bound_model(model, &bytes[..prefix]);
        let execution = GaussianCopulaExecution::restore(plan, matrix, identity, &bytes[prefix..])?;
        Ok((execution, control))
    }
}

#[cfg(test)]
mod tests;
