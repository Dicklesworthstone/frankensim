//! Frozen physical mean controls for replicated randomized quadrature.
//!
//! The numerical observable remains the original model output. A separate
//! estimator subtracts g^T(X-E[X]) using coefficients fixed before sampling.
//! Only complete equally sized scrambles enter either mean; their means, NOT
//! their dependent individual points, enter the descriptive standard error.
//! No variance-reduction guarantee, confidence interval or stopping rule is
//! implied. Finite-grid/normal-transform bias is not bounded by these errors.

use fs_blake3::{ContentHash, hash_domain};

use super::{QmcConfig, QmcEstimate, QmcExecution, scaled_mean_interruptible,
    summarize_replicates};
use crate::product_checkpoint::execution_identity;
use crate::product_control_variate::{LinearControlVariate, UqControlError,
    adjusted_samples, checkpoint, poll};
use crate::product_copula::{bound_model, physical_parameters};
use crate::{CorrelationModel, GaussianCopulaQmcExecution, ParameterUncertainty,
    PropagationMethod, UqCheckpointError, UqExecution, UqPlan, UqStatus};

/// Mean-only predictor in physical parameter units, frozen before sampling.
/// The complete sampler, dependence, physical marginals and net layout are
/// bound. Coefficients are not fitted to or selected using these observations.
/// The caller owns the adjoint/model binding and independence of the predictor;
/// constructing a fresh execution does not prove that independence.
#[derive(Debug, Clone)]
pub struct QmcLinearControlVariate {
    control: LinearControlVariate,
    sampling_identity: ContentHash,
}

/// Separate raw/adjusted mean estimates on exactly the same complete nets.
/// Adjusted responses are not a physical QoI distribution: no adjusted extrema,
/// quantiles, compliance probabilities or individual-point errors are exposed.
#[derive(Debug, Clone, PartialEq)]
pub struct QmcLinearControlEstimate {
    /// All accepted model responses, including an unfinished replicate.
    pub samples_accepted: usize,
    /// Responses actually used in this estimate; a whole number of nets.
    pub samples_in_estimate: usize,
    /// Independent completed scrambles, the statistical sample size.
    pub completed_replicates: usize,
    /// Original, immutable quadrature layout.
    pub config: QmcConfig,
    /// Raw replicate means in their original order.
    pub raw_replicate_means: Vec<f64>,
    /// Adjusted replicate means in that same order.
    pub controlled_replicate_means: Vec<f64>,
    /// Raw mean and between-replicate standard error. Absent before a full net.
    pub raw: Option<QmcEstimate>,
    /// Adjusted mean and between-replicate standard error. Absent before a full
    /// net; error remains absent until two complete nets exist.
    pub controlled: Option<QmcEstimate>,
    /// Adjusted/raw between-replicate variance, never clipped to imply benefit.
    /// Absent if error is unavailable, raw dispersion is zero or ratio overflows.
    pub variance_ratio: Option<f64>,
}

fn identity(source: &QmcExecution, physical: Option<&[ParameterUncertainty]>) -> ContentHash {
    let mut layout = Vec::with_capacity(17);
    layout.push(u8::from(physical.is_some()));
    layout.extend_from_slice(&(source.config.replicates as u64).to_le_bytes());
    layout.extend_from_slice(&(source.config.samples_per_replicate as u64).to_le_bytes());
    let model = hash_domain("org.frankensim.uq.qmc-linear-control.v1", &layout);
    let model = physical.map_or(model, |marginals| bound_model(marginals, model));
    execution_identity(&source.plan, model)
}

fn freeze(source: &QmcExecution, physical: Option<&[ParameterUncertainty]>, gradient: &[f64])
    -> Result<QmcLinearControlVariate, UqControlError>
{
    if source.status == UqStatus::Refused || source.failure.is_some() {
        return Err(UqControlError::RefusedExecution);
    }
    if source.attempted != 0 || !source.values.is_empty() {
        return Err(UqControlError::AlreadySampled);
    }
    // Admission/expectations ONLY. No Monte Carlo draws are used by this lane.
    let mut marginal_plan = source.plan.clone();
    marginal_plan.method = PropagationMethod::MonteCarlo;
    if let Some(marginals) = physical {
        marginal_plan.parameters = marginals.to_vec();
        marginal_plan.correlation = CorrelationModel::Independent;
    }
    let fresh = UqExecution::new(&marginal_plan).map_err(|_| UqControlError::PlanMismatch)?;
    Ok(QmcLinearControlVariate {
        control: fresh.freeze_linear_control_variate(gradient)?,
        sampling_identity: identity(source, physical),
    })
}

fn assess(
    source: &QmcExecution, physical: Option<&[ParameterUncertainty]>,
    control: &QmcLinearControlVariate, mut cancelled: impl FnMut() -> bool,
) -> Result<QmcLinearControlEstimate, UqControlError> {
    if identity(source, physical) != control.sampling_identity {
        return Err(UqControlError::PlanMismatch);
    }
    if source.status == UqStatus::Refused || source.failure.is_some() {
        return Err(UqControlError::RefusedExecution);
    }
    poll(&mut cancelled)?;
    let size = source.config.samples_per_replicate;
    let complete = source.values.len() / size;
    let count = complete * size;
    let values = &source.values[..count];
    let adjusted = adjusted_samples(&control.control, values, |ordinal| {
        let parameters = source.parameters(ordinal).map_err(|_| UqControlError::NumericalRange)?;
        match physical {
            None => Ok(parameters),
            Some(marginals) => physical_parameters(marginals, &parameters)
                .map_err(|_| UqControlError::NumericalRange),
        }
    }, &mut cancelled)?;
    let mut raw_means = Vec::with_capacity(complete);
    let mut controlled_means = Vec::with_capacity(complete);
    for (raw, controlled) in values.chunks_exact(size).zip(adjusted.chunks_exact(size)) {
        raw_means.push(scaled_mean_interruptible(raw, &mut cancelled)?);
        controlled_means.push(scaled_mean_interruptible(controlled, &mut cancelled)?);
    }
    poll(&mut cancelled)?;
    // At most 256 means: the existing owner supplies exactly the same raw
    // numerical result as QmcExecution::report, including finite-range handling.
    let raw = summarize_replicates(&raw_means);
    let controlled = summarize_replicates(&controlled_means);
    let variance_ratio = raw.as_ref().and_then(|r| r.standard_error)
        .zip(controlled.as_ref().and_then(|r| r.standard_error))
        .and_then(|(a, b)| {
            if a == 0.0 { return None; }
            let ratio = (b / a).powi(2);
            ratio.is_finite().then_some(ratio)
        });
    poll(&mut cancelled)?;
    Ok(QmcLinearControlEstimate {
        samples_accepted: source.values.len(), samples_in_estimate: count,
        completed_replicates: complete, config: source.config,
        raw_replicate_means: raw_means, controlled_replicate_means: controlled_means,
        raw, controlled, variance_ratio,
    })
}

impl QmcExecution {
    /// Freeze derivatives in the original parameter units before any sample.
    /// Marginal expectations are taken from the admitted probability plan.
    pub fn freeze_linear_control_variate(&self, gradient: &[f64])
        -> Result<QmcLinearControlVariate, UqControlError>
    { freeze(self, None, gradient) }

    /// Replay only the cheap sampler; preserve all physical observations.
    pub fn assess_linear_control_variate(&self, control: &QmcLinearControlVariate)
        -> Result<QmcLinearControlEstimate, UqControlError>
    { assess(self, None, control, || false) }

    /// Complete-net replay in O(n*d) work and O(n+d) additional storage. Poll
    /// before each input and every 512 reduction entries. Failed/cancelled
    /// assessments publish no partial estimate and never alter execution.
    pub fn assess_linear_control_variate_interruptible(
        &self, control: &QmcLinearControlVariate, cancelled: impl FnMut() -> bool,
    ) -> Result<QmcLinearControlEstimate, UqControlError>
    { assess(self, None, control, cancelled) }
}

impl GaussianCopulaQmcExecution {
    /// Freeze PHYSICAL derivatives, not latent-normal derivatives. Means are
    /// the declared uniform midpoints; matrix and physical transformation stay
    /// bound separately. Another sample stream's control cannot be attached.
    pub fn freeze_linear_control_variate(&self, gradient: &[f64])
        -> Result<QmcLinearControlVariate, UqControlError>
    { freeze(self.control_source(), Some(self.marginals()), gradient) }

    /// Assess physical-coordinate controls on complete independent scrambles.
    pub fn assess_linear_control_variate(&self, control: &QmcLinearControlVariate)
        -> Result<QmcLinearControlEstimate, UqControlError>
    { assess(self.control_source(), Some(self.marginals()), control, || false) }

    /// As the ordinary QMC assessment, replaying the exact original copula map.
    /// Cancellation never changes raw outcomes, compliance or checkpoints.
    pub fn assess_linear_control_variate_interruptible(
        &self, control: &QmcLinearControlVariate, cancelled: impl FnMut() -> bool,
    ) -> Result<QmcLinearControlEstimate, UqControlError>
    { assess(self.control_source(), Some(self.marginals()), control, cancelled) }
}

impl QmcLinearControlVariate {
    /// Analytic physical expectations in parameter declaration order.
    #[must_use]
    pub fn parameter_means(&self) -> &[f64] { self.control.parameter_means() }
    /// Original physical-unit coefficients, never a fitted sample regression.
    #[must_use]
    pub fn gradient(&self) -> &[f64] { self.control.gradient() }

    /// Save coefficients with RAW QMC observations, including a partial net.
    /// The original checkpoint owner binds layout/ordinal/status and model.
    pub fn checkpoint(&self, execution: &QmcExecution, model: ContentHash)
        -> Result<Vec<u8>, UqCheckpointError>
    {
        if identity(execution, None) != self.sampling_identity {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let mut bytes = checkpoint::encode_header(self.gradient());
        let bound = checkpoint::bound_model(model, &bytes);
        bytes.extend_from_slice(&execution.checkpoint(bound)?);
        Ok(bytes)
    }

    /// Save the physical-copula variant, retaining its physical supports/units
    /// as well as latent matrix and QMC layout. This is not an inner-latent save.
    pub fn checkpoint_copula(&self, execution: &GaussianCopulaQmcExecution, model: ContentHash)
        -> Result<Vec<u8>, UqCheckpointError>
    {
        if identity(execution.control_source(), Some(execution.marginals())) != self.sampling_identity {
            return Err(UqCheckpointError::IdentityMismatch);
        }
        let mut bytes = checkpoint::encode_header(self.gradient());
        let bound = checkpoint::bound_model(model, &bytes);
        bytes.extend_from_slice(&execution.checkpoint(bound)?);
        Ok(bytes)
    }

    /// Restore without refitting, re-solving a nominal point or replaying model
    /// evaluations. Coefficient substitution and any sampler/model change refuse.
    pub fn restore(plan: &UqPlan, config: QmcConfig, model: ContentHash, bytes: &[u8])
        -> Result<(QmcExecution, Self), UqCheckpointError>
    {
        let fresh = QmcExecution::new(plan, config).map_err(UqCheckpointError::InvalidPlan)?;
        let (gradient, prefix) = checkpoint::decode_header(plan.parameters.len(), plan.budget_max_samples, bytes)?;
        let control = fresh.freeze_linear_control_variate(&gradient)
            .map_err(|_| UqCheckpointError::InvalidEncoding("invalid QMC control coefficients"))?;
        let bound = checkpoint::bound_model(model, &bytes[..prefix]);
        let execution = QmcExecution::restore(plan, config, bound, &bytes[prefix..])?;
        Ok((execution, control))
    }

    /// Restore physical-copula controls and original raw observations together.
    /// Means are reconstructed from the admitted physical marginals. Checksums
    /// detect corruption but do not authenticate the data or its coefficient choice.
    pub fn restore_copula(plan: &UqPlan, matrix: &[Vec<f64>], config: QmcConfig,
        model: ContentHash, bytes: &[u8]) -> Result<(GaussianCopulaQmcExecution, Self), UqCheckpointError>
    {
        let fresh = GaussianCopulaQmcExecution::new(plan, matrix, config).map_err(UqCheckpointError::InvalidPlan)?;
        let (gradient, prefix) = checkpoint::decode_header(plan.parameters.len(), plan.budget_max_samples, bytes)?;
        let control = fresh.freeze_linear_control_variate(&gradient)
            .map_err(|_| UqCheckpointError::InvalidEncoding("invalid physical QMC control coefficients"))?;
        let bound = checkpoint::bound_model(model, &bytes[..prefix]);
        let execution = GaussianCopulaQmcExecution::restore(plan, matrix, config, bound, &bytes[prefix..])?;
        Ok((execution, control))
    }
}

#[cfg(test)]
mod tests;
