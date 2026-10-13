//! One-step, outcome-free sensor design under the retained joint posterior.
//!
//! Scores are reduction of a caller's linear physical target variance, not
//! EVPI, a calibrated utility, or a statistical-model certificate. Every probe
//! is evaluated against the same unchanged history; this is not joint selection
//! of a whole sensor set. Shared nuisance uncertainty and learned cross-terms
//! participate in exactly the same augmented observation used for acquisition.

use super::{SharedNoiseStream, copy_string, reserve, zeros};
use crate::factor::FactorBelief;
use crate::groups::RowId;
use crate::groups::jacobian::{JacobianRecord, SharedNoiseModel};
use crate::{AssimError, Cx, Observation, WorkPlan, WorkProgress, validate_observation_numeric_parts};
use fs_blake3::ContentHash;

/// A prospective reading with NO measured value or caller-supplied truth.
#[derive(Debug, Clone, PartialEq)]
pub struct SharedNoiseProbe {
    row: RowId,
    model_id: String,
    physical_dim: usize,
    operator: Vec<f64>,
    noise_variance: f64,
}

impl SharedNoiseProbe {
    /// Declare the sensor response and independent instrument variance.
    ///
    /// # Errors
    /// Refuses invalid operators/Jacobians, noise, dimensions or allocation.
    pub fn new(
        model: &SharedNoiseModel,
        row: RowId,
        physical_operator: Vec<f64>,
        jacobian: Vec<f64>,
        independent_noise_variance: f64,
    ) -> Result<Self, AssimError> {
        validate_observation_numeric_parts(&physical_operator, 0.0, independent_noise_variance)?;
        let physical_dim = physical_operator.len();
        // Reuse model-bound Jacobian admission, but store no pretend outcome.
        let checked = JacobianRecord::new(
            Observation::new(physical_operator, 0.0, independent_noise_variance, "design-probe")?,
            row, model, jacobian,
        )?;
        let d = physical_dim + checked.jacobian().len();
        crate::validate_state_dimension(d)?;
        let mut operator = reserve(d)?;
        operator.extend_from_slice(checked.observation().operator());
        operator.extend_from_slice(checked.jacobian());
        Ok(Self {
            row: checked.row().clone(), model_id: copy_string(model.content_id())?,
            physical_dim, operator, noise_variance: independent_noise_variance,
        })
    }

    /// Identity of the prospective NEW dataset row.
    #[must_use]
    pub const fn row(&self) -> &RowId { &self.row }

    /// Physical sensor operator.
    #[must_use]
    pub fn physical_operator(&self) -> &[f64] { &self.operator[..self.physical_dim] }

    /// Signed shared-error sensitivities in the admitted model's order.
    #[must_use]
    pub fn jacobian(&self) -> &[f64] { &self.operator[self.physical_dim..] }

    /// Independent instrument noise; shared uncertainty is retained separately.
    #[must_use]
    pub const fn noise_variance(&self) -> f64 { self.noise_variance }

    /// Attach a REAL acquired value to the exact response/noise/row declaration.
    /// Planning never calls this method or reads an oracle's truth.
    ///
    /// # Errors
    /// Refuses another model, invalid reading/identity, or allocation failure.
    pub fn record(&self, model: &SharedNoiseModel, value: f64, instrument: impl Into<String>) -> Result<JacobianRecord, AssimError> {
        if model.content_id() != self.model_id {
            return Err(AssimError::SharedSourceDeclaration {
                id: model.identity().to_owned(), reason: "probe differs from acquisition model",
            });
        }
        let mut operator = reserve(self.physical_dim)?;
        operator.extend_from_slice(self.physical_operator());
        let mut jacobian = reserve(self.jacobian().len())?;
        jacobian.extend_from_slice(self.jacobian());
        JacobianRecord::new(
            Observation::new(operator, value, self.noise_variance, instrument)?,
            self.row.clone(), model, jacobian,
        )
    }
}

/// Conditional Gaussian moments and a one-step variance-reduction score.
#[derive(Debug, Clone, PartialEq)]
pub struct ProbeScore {
    row: RowId,
    history: ContentHash,
    predicted_mean: f64,
    innovation_variance: f64,
    prior_target_variance: f64,
    posterior_target_variance: f64,
    variance_reduction: f64,
}

impl ProbeScore {
    /// Candidate row identity.
    #[must_use]
    pub const fn row(&self) -> &RowId { &self.row }
    /// Retained history at which the probe was evaluated.
    #[must_use]
    pub const fn history_id(&self) -> ContentHash { self.history }
    /// Mean of the prospective measurement, including learned nuisance mean.
    #[must_use]
    pub const fn predicted_mean(&self) -> f64 { self.predicted_mean }
    /// Measurement variance INCLUDING independent noise and joint uncertainty.
    #[must_use]
    pub const fn innovation_variance(&self) -> f64 { self.innovation_variance }
    /// Current variance of the caller's physical target.
    #[must_use]
    pub const fn prior_target_variance(&self) -> f64 { self.prior_target_variance }
    /// Predicted target variance after the reading, independent of its outcome.
    #[must_use]
    pub const fn posterior_target_variance(&self) -> f64 { self.posterior_target_variance }
    /// `Cov(target, reading)^2 / Var(reading)`, not a fixed-percent heuristic.
    #[must_use]
    pub const fn variance_reduction(&self) -> f64 { self.variance_reduction }
}

impl SharedNoiseStream {
    /// Rank prospective sensors by expected reduction of `target * x` variance.
    /// Ties use row identity. All probes refer to the SAME current posterior;
    /// re-rank after acquiring a reading. The stream is never mutated here.
    ///
    /// Computation uses `sqrt(D) U^T` actions, not dense covariance construction.
    /// The posterior target variance is evaluated as a Joseph sum of squares,
    /// avoiding catastrophic subtraction when the target is nearly determined.
    /// One work/poll accountant covers admission, all probes and ranking.
    ///
    /// # Errors
    /// Refuses invalid targets, duplicate/already-acquired rows, another model,
    /// an exhausted row cap, more than 256 probes, unresolved numerical range or
    /// contraction, allocation failure, cancellation or exhausted work/polls.
    pub fn rank_probes(&self, target: &[f64], probes: &[SharedNoiseProbe], cx: &Cx<'_>) -> Result<Vec<ProbeScore>, AssimError> {
        if target.len() != self.physical_dim {
            return Err(AssimError::DimMismatch { state: self.physical_dim, operator: target.len() });
        }
        validate_observation_numeric_parts(target, 0.0, 1.0)?;
        if probes.is_empty() { return Err(AssimError::EmptyObservations); }
        if probes.len() > 256 {
            return Err(AssimError::ObservationCountLimit { count: probes.len(), max: 256 });
        }
        if self.seen.len() == self.row_capacity {
            return Err(AssimError::ObservationCountLimit { count: self.seen.len() + 1, max: self.row_capacity });
        }
        let d = self.joint.dim();
        let (mu, du) = (probes.len() as u128, d as u128);
        let work = (mu + 1) * (8 * du * du + 256 * du + 16_384) + 256 * mu * mu + 1_024;
        let mut progress = WorkProgress::new(cx, WorkPlan::checked(0, 0, 0, 0, work, 0, 0)?)?;
        progress.checkpoint("shared probe admission")?;
        for (index, probe) in probes.iter().enumerate() {
            if probe.model_id != self.model_id {
                return Err(AssimError::SharedSourceDeclaration { id: probe.row.as_str().to_owned(), reason: "probe differs from retained shared model" });
            }
            if probe.physical_dim != self.physical_dim {
                return Err(AssimError::DimMismatch { state: self.physical_dim, operator: probe.physical_dim });
            }
            progress.comparison_bytes("shared probe row lookup", 13 * 256)?;
            if self.seen.binary_search_by(|seen| seen.row.cmp(&probe.row)).is_ok() {
                return Err(AssimError::DuplicateDatasetRow { id: probe.row.as_str().to_owned() });
            }
            for prior in &probes[..index] {
                progress.comparison_bytes("shared probe uniqueness", 256)?;
                if prior.row == probe.row {
                    return Err(AssimError::DuplicateDatasetRow { id: probe.row.as_str().to_owned() });
                }
            }
        }
        let mut augmented_target = zeros(d)?;
        augmented_target[..target.len()].copy_from_slice(target);
        let target_factor = project(&self.joint, &augmented_target, &mut progress)?;
        let target_root = norm(&target_factor, &mut progress)?;
        let prior_variance = square(target_root)?;
        let mut scores = reserve(probes.len())?;
        for probe in probes {
            let mut sensor_factor = project(&self.joint, &probe.operator, &mut progress)?;
            sensor_factor.push(probe.noise_variance.sqrt());
            let innovation_root = norm(&sensor_factor, &mut progress)?;
            let innovation_variance = square(innovation_root)?;
            let predicted_mean = dot(&probe.operator, self.joint.mean(), &mut progress)?;
            let (posterior, reduction) = if target_root == 0.0 { (0.0, 0.0) } else {
                let mut tn = reserve(d)?;
                let mut sn = reserve(d)?;
                for index in 0..d {
                    tn.push(target_factor[index] / target_root);
                    sn.push(sensor_factor[index] / innovation_root);
                    progress.scalar("shared probe normalization", 2)?;
                }
                let correlation = dot(&tn, &sn, &mut progress)?;
                if correlation.abs() > 1.0 { return Err(AssimError::CovarianceCertificationUnresolved); }
                let reduction = square(product(target_root, correlation)?)?;
                let posterior = if correlation == 0.0 { prior_variance } else {
                    let mut residual = reserve(d + 1)?;
                    for index in 0..d {
                        residual.push((-correlation).mul_add(sn[index], tn[index]));
                        progress.scalar("shared probe conditional variance", 2)?;
                    }
                    residual.push(correlation * (probe.noise_variance.sqrt() / innovation_root));
                    square(product(target_root, norm(&residual, &mut progress)?)?)?
                };
                if posterior > prior_variance || reduction > prior_variance {
                    return Err(AssimError::CovarianceCertificationUnresolved);
                }
                (posterior, reduction)
            };
            scores.push(ProbeScore { row: probe.row.clone(), history: self.history, predicted_mean,
                innovation_variance, prior_target_variance: prior_variance,
                posterior_target_variance: posterior, variance_reduction: reduction });
        }
        progress.comparison_bytes("shared probe ordering", mu * 32 * 256)?;
        scores.sort_unstable_by(|left, right| right.variance_reduction.total_cmp(&left.variance_reduction).then_with(|| left.row.cmp(&right.row)));
        progress.checkpoint("shared probe publication")?;
        Ok(scores)
    }
}

fn project(belief: &FactorBelief, operator: &[f64], progress: &mut WorkProgress<'_, '_>) -> Result<Vec<f64>, AssimError> {
    // Reserve the extra independent-noise coordinate used by probe norms.
    let mut result = reserve(belief.dim() + 1)?;
    for column in 0..belief.dim() {
        let mut value = operator[column];
        for (row, coefficient) in operator.iter().enumerate().take(column) {
            value = belief.upper_entry(row, column).ok_or(AssimError::CovarianceCertificationUnresolved)?.mul_add(*coefficient, value);
            progress.scalar("shared probe factor action", 2)?;
        }
        finite(value)?;
        let pivot = belief.diag(column).ok_or(AssimError::CovarianceCertificationUnresolved)?;
        let projected = value * pivot.sqrt();
        finite(projected)?;
        if projected == 0.0 && value != 0.0 && pivot != 0.0 { return Err(AssimError::CovarianceCertificationUnresolved); }
        result.push(projected);
        progress.scalar("shared probe factor action", 4)?;
    }
    Ok(result)
}

fn dot(left: &[f64], right: &[f64], progress: &mut WorkProgress<'_, '_>) -> Result<f64, AssimError> {
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for (a, b) in left.iter().zip(right) {
        let term = *a * *b;
        let next = sum + term;
        correction += if sum.abs() >= term.abs() { (sum - next) + term } else { (term - next) + sum };
        sum = next;
        progress.scalar("shared probe inner product", 8)?;
    }
    finite(sum + correction)
}

fn norm(values: &[f64], progress: &mut WorkProgress<'_, '_>) -> Result<f64, AssimError> {
    let mut scale = 0.0_f64;
    for value in values { scale = scale.max(value.abs()); progress.scalar("shared probe norm", 1)?; }
    if scale == 0.0 { return Ok(0.0); }
    let (mut sum, mut correction) = (0.0_f64, 0.0_f64);
    for value in values {
        let ratio = *value / scale;
        let term = ratio * ratio;
        let next = sum + term;
        correction += if sum >= term { (sum - next) + term } else { (term - next) + sum };
        sum = next;
        progress.scalar("shared probe norm", 8)?;
    }
    finite(scale * (sum + correction).sqrt())
}

fn product(left: f64, right: f64) -> Result<f64, AssimError> {
    let result = finite(left * right)?;
    if result == 0.0 && left != 0.0 && right != 0.0 {
        return Err(AssimError::CovarianceCertificationUnresolved);
    }
    Ok(result)
}

fn square(root: f64) -> Result<f64, AssimError> {
    let result = finite(root * root)?;
    if result == 0.0 && root != 0.0 { return Err(AssimError::CovarianceCertificationUnresolved); }
    Ok(result)
}

fn finite(value: f64) -> Result<f64, AssimError> {
    if value.is_finite() { Ok(value) }
    else { Err(AssimError::NonFiniteComputation { stage: "shared probe prediction" }) }
}
