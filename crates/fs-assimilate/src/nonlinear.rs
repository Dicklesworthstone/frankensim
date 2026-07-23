//! Receipted finite-difference linearization for nonlinear scalar observations.
//!
//! This module owns one bounded extended-Kalman seam: construct a local affine
//! observation at the prior mean, refuse numerically invalid probes, and
//! withhold publication when the standardized innovation exceeds a declared
//! validity gate. The actual covariance update remains owned by `assimilate`.

use fs_blake3::DomainHasher;
use fs_exec::Cx;

use super::{
    AssimError, Belief, Observation, assimilate, canonical_f64_bits, canonicalize_zero,
    validate_leaf_identity, validate_noise,
};

const LINEARIZATION_ID_DOMAIN: &str = "org.frankensim.fs-assimilate.nonlinear-fd-linearization.v1";
const LINEARIZATION_ID_PREFIX: &str = "nonlinear-fd-linearization:v1:";

/// Checked settings for scale-aware central finite differences and the
/// standardized-innovation validity gate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteDifferenceSettings {
    relative_step: f64,
    max_standardized_innovation: f64,
}

impl FiniteDifferenceSettings {
    /// Construct checked settings.
    ///
    /// The nominal perturbation for state component `x_i` is
    /// `relative_step * max(1, |x_i|)`. The actual representable plus/minus
    /// states and denominator are retained in the probe receipt.
    pub fn new(relative_step: f64, max_standardized_innovation: f64) -> Result<Self, AssimError> {
        if !relative_step.is_finite() || relative_step <= 0.0 {
            return Err(AssimError::InvalidLinearizationParameter {
                parameter: "relative finite-difference step",
            });
        }
        if !max_standardized_innovation.is_finite() || max_standardized_innovation <= 0.0 {
            return Err(AssimError::InvalidLinearizationParameter {
                parameter: "maximum standardized innovation",
            });
        }
        Ok(Self {
            relative_step: canonicalize_zero(relative_step),
            max_standardized_innovation: canonicalize_zero(max_standardized_innovation),
        })
    }

    /// Relative finite-difference scale.
    #[must_use]
    pub fn relative_step(self) -> f64 {
        self.relative_step
    }

    /// Largest standardized innovation admitted by the local linearization.
    #[must_use]
    pub fn max_standardized_innovation(self) -> f64 {
        self.max_standardized_innovation
    }
}

/// Checked scalar nonlinear observation declaration.
///
/// The model identity is semantic input to the receipt. The callback supplied
/// later must be a pure implementation of that identity; this crate cannot
/// authenticate arbitrary closure code.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearObservationSpec {
    reading: f64,
    noise_var: f64,
    instrument: String,
    model_id: String,
    settings: FiniteDifferenceSettings,
}

impl NonlinearObservationSpec {
    /// Construct a checked nonlinear observation declaration.
    pub fn new(
        reading: f64,
        noise_var: f64,
        instrument: impl Into<String>,
        model_id: impl Into<String>,
        settings: FiniteDifferenceSettings,
    ) -> Result<Self, AssimError> {
        if !reading.is_finite() {
            return Err(AssimError::NonFiniteObservationValue);
        }
        validate_noise(noise_var)?;
        let instrument = instrument.into();
        validate_leaf_identity("instrument", &instrument)?;
        let model_id = model_id.into();
        validate_leaf_identity("nonlinear_model", &model_id)?;
        Ok(Self {
            reading: canonicalize_zero(reading),
            noise_var: canonicalize_zero(noise_var),
            instrument,
            model_id,
            settings,
        })
    }

    /// Measured scalar value.
    #[must_use]
    pub fn reading(&self) -> f64 {
        self.reading
    }

    /// Declared measurement-noise variance.
    #[must_use]
    pub fn noise_var(&self) -> f64 {
        self.noise_var
    }

    /// Calibrated instrument identity.
    #[must_use]
    pub fn instrument(&self) -> &str {
        &self.instrument
    }

    /// Caller-declared nonlinear model identity.
    #[must_use]
    pub fn model_id(&self) -> &str {
        &self.model_id
    }

    /// Checked finite-difference and innovation-gate settings.
    #[must_use]
    pub fn settings(&self) -> FiniteDifferenceSettings {
        self.settings
    }
}

/// One retained central finite-difference experiment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteDifferenceProbe {
    component: usize,
    nominal_step: f64,
    plus_state: f64,
    minus_state: f64,
    plus_prediction: f64,
    minus_prediction: f64,
}

impl FiniteDifferenceProbe {
    /// Perturbed state component.
    #[must_use]
    pub fn component(self) -> usize {
        self.component
    }

    /// Requested scale-aware step before binary64 rounding.
    #[must_use]
    pub fn nominal_step(self) -> f64 {
        self.nominal_step
    }

    /// Actual representable plus state.
    #[must_use]
    pub fn plus_state(self) -> f64 {
        self.plus_state
    }

    /// Actual representable minus state.
    #[must_use]
    pub fn minus_state(self) -> f64 {
        self.minus_state
    }

    /// Model prediction at the plus state.
    #[must_use]
    pub fn plus_prediction(self) -> f64 {
        self.plus_prediction
    }

    /// Model prediction at the minus state.
    #[must_use]
    pub fn minus_prediction(self) -> f64 {
        self.minus_prediction
    }

    /// Actual finite-difference denominator.
    #[must_use]
    pub fn denominator(self) -> f64 {
        self.plus_state - self.minus_state
    }
}

/// Innovation-gate disposition for one nonlinear linearization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LinearizationDisposition {
    /// The standardized innovation was inside the declared gate and a checked
    /// affine observation was published.
    Admitted,
    /// The innovation was outside the declared local-validity region. No
    /// affine observation or posterior was published.
    DemotedLargeInnovation,
}

/// Immutable nonlinear finite-difference linearization report.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearLinearization {
    prediction: f64,
    innovation: f64,
    innovation_variance: f64,
    standardized_innovation: f64,
    jacobian: Vec<f64>,
    probes: Vec<FiniteDifferenceProbe>,
    observation: Option<Observation>,
    disposition: LinearizationDisposition,
    receipt_id: String,
}

impl NonlinearLinearization {
    /// Model prediction `h(mean)` at the linearization point.
    #[must_use]
    pub fn prediction(&self) -> f64 {
        self.prediction
    }

    /// Measurement innovation `reading - h(mean)`.
    #[must_use]
    pub fn innovation(&self) -> f64 {
        self.innovation
    }

    /// Local innovation variance `H P H^T + R`.
    #[must_use]
    pub fn innovation_variance(&self) -> f64 {
        self.innovation_variance
    }

    /// Absolute innovation divided by its local standard deviation.
    #[must_use]
    pub fn standardized_innovation(&self) -> f64 {
        self.standardized_innovation
    }

    /// Central finite-difference Jacobian row.
    #[must_use]
    pub fn jacobian(&self) -> &[f64] {
        &self.jacobian
    }

    /// Actual step and prediction receipts in state-component order.
    #[must_use]
    pub fn probes(&self) -> &[FiniteDifferenceProbe] {
        &self.probes
    }

    /// Checked local affine observation, present only when admitted.
    #[must_use]
    pub fn observation(&self) -> Option<&Observation> {
        self.observation.as_ref()
    }

    /// Validity-gate disposition.
    #[must_use]
    pub fn disposition(&self) -> LinearizationDisposition {
        self.disposition
    }

    /// Domain-separated receipt binding all numerical inputs and outputs.
    #[must_use]
    pub fn receipt_id(&self) -> &str {
        &self.receipt_id
    }
}

/// Nonlinear extended-Kalman seam result.
///
/// A demoted linearization contains no posterior. An admitted one is fused by
/// the established scalar Joseph updater.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearAssimilation {
    linearization: NonlinearLinearization,
    posterior: Option<Belief>,
}

impl NonlinearAssimilation {
    /// Complete linearization and gate receipt.
    #[must_use]
    pub fn linearization(&self) -> &NonlinearLinearization {
        &self.linearization
    }

    /// Updated belief, absent for a demoted linearization.
    #[must_use]
    pub fn posterior(&self) -> Option<&Belief> {
        self.posterior.as_ref()
    }
}

/// Linearize a caller-identified nonlinear scalar model at the prior mean.
///
/// A successful call evaluates the callback exactly `1 + 2 * prior.dim()`
/// times. It must be a deterministic, side-effect-free implementation of
/// `spec.model_id()`.
pub fn linearize_nonlinear_fd<F>(
    prior: &Belief,
    spec: &NonlinearObservationSpec,
    predict: F,
    cx: &Cx<'_>,
) -> Result<NonlinearLinearization, AssimError>
where
    F: Fn(&[f64]) -> f64,
{
    cx.checkpoint()
        .map_err(|_| nonlinear_cancelled("nonlinear preflight", 0, prior.dim() as u128 + 2))?;
    let prediction = predict(prior.mean());
    if !prediction.is_finite() {
        return Err(AssimError::NonFiniteModelPrediction {
            evaluation: "base",
            component: None,
        });
    }

    let mut state = prior.mean().to_vec();
    let mut jacobian = Vec::with_capacity(prior.dim());
    let mut probes = Vec::with_capacity(prior.dim());
    for component in 0..prior.dim() {
        cx.checkpoint().map_err(|_| {
            nonlinear_cancelled(
                "nonlinear finite differences",
                component as u128,
                prior.dim() as u128 + 2,
            )
        })?;
        let (derivative, probe) =
            central_difference_probe(&mut state, component, spec.settings.relative_step, &predict)?;
        jacobian.push(derivative);
        probes.push(probe);
    }

    let innovation = spec.reading - prediction;
    if !innovation.is_finite() {
        return Err(AssimError::NonFiniteComputation {
            stage: "nonlinear innovation",
        });
    }
    let innovation_variance = innovation_variance(prior, &jacobian, spec.noise_var)?;
    let standardized_innovation = innovation.abs() / innovation_variance.sqrt();
    if !standardized_innovation.is_finite() {
        return Err(AssimError::NonFiniteComputation {
            stage: "standardized nonlinear innovation",
        });
    }
    let disposition = if standardized_innovation <= spec.settings.max_standardized_innovation {
        LinearizationDisposition::Admitted
    } else {
        LinearizationDisposition::DemotedLargeInnovation
    };
    let observation = if disposition == LinearizationDisposition::Admitted {
        let affine_value = spec.reading - prediction + checked_dot(&jacobian, prior.mean())?;
        Some(Observation::new(
            jacobian.clone(),
            affine_value,
            spec.noise_var,
            spec.instrument.clone(),
        )?)
    } else {
        None
    };
    cx.checkpoint().map_err(|_| {
        nonlinear_cancelled(
            "nonlinear receipt",
            prior.dim() as u128 + 1,
            prior.dim() as u128 + 2,
        )
    })?;
    let receipt_id = linearization_receipt(
        prior,
        spec,
        prediction,
        innovation,
        innovation_variance,
        standardized_innovation,
        &jacobian,
        &probes,
        disposition,
    );
    Ok(NonlinearLinearization {
        prediction: canonicalize_zero(prediction),
        innovation: canonicalize_zero(innovation),
        innovation_variance: canonicalize_zero(innovation_variance),
        standardized_innovation: canonicalize_zero(standardized_innovation),
        jacobian,
        probes,
        observation,
        disposition,
        receipt_id,
    })
}

/// Linearize and, when the innovation gate admits, fuse the resulting affine
/// observation through the existing scalar Joseph updater.
pub fn assimilate_nonlinear_fd<F>(
    prior: &Belief,
    spec: &NonlinearObservationSpec,
    predict: F,
    cx: &Cx<'_>,
) -> Result<NonlinearAssimilation, AssimError>
where
    F: Fn(&[f64]) -> f64,
{
    let linearization = linearize_nonlinear_fd(prior, spec, predict, cx)?;
    let posterior = match linearization.observation() {
        Some(observation) => Some(assimilate(prior, observation, cx)?),
        None => None,
    };
    cx.checkpoint().map_err(|_| {
        nonlinear_cancelled(
            "nonlinear publication",
            prior.dim() as u128 + 1,
            prior.dim() as u128 + 2,
        )
    })?;
    Ok(NonlinearAssimilation {
        linearization,
        posterior,
    })
}

fn central_difference_probe<F>(
    state: &mut [f64],
    component: usize,
    relative_step: f64,
    predict: &F,
) -> Result<(f64, FiniteDifferenceProbe), AssimError>
where
    F: Fn(&[f64]) -> f64,
{
    let base = state[component];
    let nominal_step = relative_step * base.abs().max(1.0);
    if !nominal_step.is_finite() || nominal_step <= 0.0 {
        return Err(AssimError::FiniteDifferenceStepUnrepresentable { component });
    }
    let plus_state = base + nominal_step;
    let minus_state = base - nominal_step;
    if !plus_state.is_finite()
        || !minus_state.is_finite()
        || canonical_f64_bits(plus_state) == canonical_f64_bits(base)
        || canonical_f64_bits(minus_state) == canonical_f64_bits(base)
        || canonical_f64_bits(plus_state) == canonical_f64_bits(minus_state)
    {
        return Err(AssimError::FiniteDifferenceStepUnrepresentable { component });
    }
    let denominator = plus_state - minus_state;
    if !denominator.is_finite() || denominator <= 0.0 {
        return Err(AssimError::FiniteDifferenceStepUnrepresentable { component });
    }

    state[component] = plus_state;
    let plus_prediction = predict(state);
    state[component] = minus_state;
    let minus_prediction = predict(state);
    state[component] = base;
    if !plus_prediction.is_finite() {
        return Err(AssimError::NonFiniteModelPrediction {
            evaluation: "plus",
            component: Some(component),
        });
    }
    if !minus_prediction.is_finite() {
        return Err(AssimError::NonFiniteModelPrediction {
            evaluation: "minus",
            component: Some(component),
        });
    }
    let derivative = (plus_prediction - minus_prediction) / denominator;
    if !derivative.is_finite() {
        return Err(AssimError::NonFiniteObservationJacobian { component });
    }
    Ok((
        canonicalize_zero(derivative),
        FiniteDifferenceProbe {
            component,
            nominal_step,
            plus_state,
            minus_state,
            plus_prediction: canonicalize_zero(plus_prediction),
            minus_prediction: canonicalize_zero(minus_prediction),
        },
    ))
}

fn innovation_variance(
    prior: &Belief,
    jacobian: &[f64],
    noise_var: f64,
) -> Result<f64, AssimError> {
    let mut covariance_times_jacobian = Vec::with_capacity(prior.dim());
    for row in prior.covariance() {
        covariance_times_jacobian.push(checked_dot(row, jacobian)?);
    }
    let variance = checked_dot(jacobian, &covariance_times_jacobian)? + noise_var;
    if !variance.is_finite() {
        return Err(AssimError::NonFiniteComputation {
            stage: "nonlinear innovation variance",
        });
    }
    if variance <= 0.0 {
        return Err(AssimError::SingularInnovation);
    }
    Ok(variance)
}

fn checked_dot(left: &[f64], right: &[f64]) -> Result<f64, AssimError> {
    let mut sum = 0.0;
    for (&left, &right) in left.iter().zip(right) {
        sum = left.mul_add(right, sum);
        if !sum.is_finite() {
            return Err(AssimError::NonFiniteComputation {
                stage: "nonlinear linearization dot product",
            });
        }
    }
    Ok(sum)
}

#[allow(clippy::too_many_arguments)]
fn linearization_receipt(
    prior: &Belief,
    spec: &NonlinearObservationSpec,
    prediction: f64,
    innovation: f64,
    innovation_variance: f64,
    standardized_innovation: f64,
    jacobian: &[f64],
    probes: &[FiniteDifferenceProbe],
    disposition: LinearizationDisposition,
) -> String {
    let mut hasher = DomainHasher::new(LINEARIZATION_ID_DOMAIN);
    hash_belief(&mut hasher, prior);
    hash_f64(&mut hasher, spec.reading);
    hash_f64(&mut hasher, spec.noise_var);
    hash_bytes(&mut hasher, spec.instrument.as_bytes());
    hash_bytes(&mut hasher, spec.model_id.as_bytes());
    hash_f64(&mut hasher, spec.settings.relative_step);
    hash_f64(&mut hasher, spec.settings.max_standardized_innovation);
    hash_f64(&mut hasher, prediction);
    hash_f64(&mut hasher, innovation);
    hash_f64(&mut hasher, innovation_variance);
    hash_f64(&mut hasher, standardized_innovation);
    hash_u64(&mut hasher, jacobian.len() as u64);
    for derivative in jacobian {
        hash_f64(&mut hasher, *derivative);
    }
    hash_u64(&mut hasher, probes.len() as u64);
    for probe in probes {
        hash_u64(&mut hasher, probe.component as u64);
        hash_f64(&mut hasher, probe.nominal_step);
        hash_f64(&mut hasher, probe.plus_state);
        hash_f64(&mut hasher, probe.minus_state);
        hash_f64(&mut hasher, probe.plus_prediction);
        hash_f64(&mut hasher, probe.minus_prediction);
    }
    hasher.update(&[match disposition {
        LinearizationDisposition::Admitted => 0,
        LinearizationDisposition::DemotedLargeInnovation => 1,
    }]);
    format!("{LINEARIZATION_ID_PREFIX}{}", hasher.finalize())
}

fn hash_belief(hasher: &mut DomainHasher, belief: &Belief) {
    hash_u64(hasher, belief.dim() as u64);
    for value in belief.mean() {
        hash_f64(hasher, *value);
    }
    hash_u64(hasher, belief.covariance().len() as u64);
    for row in belief.covariance() {
        hash_u64(hasher, row.len() as u64);
        for value in row {
            hash_f64(hasher, *value);
        }
    }
}

fn hash_bytes(hasher: &mut DomainHasher, bytes: &[u8]) {
    hash_u64(hasher, bytes.len() as u64);
    hasher.update(bytes);
}

fn hash_f64(hasher: &mut DomainHasher, value: f64) {
    hasher.update(&canonical_f64_bits(value).to_le_bytes());
}

fn hash_u64(hasher: &mut DomainHasher, value: u64) {
    hasher.update(&value.to_le_bytes());
}

fn nonlinear_cancelled(phase: &'static str, completed: u128, planned: u128) -> AssimError {
    AssimError::Cancelled {
        phase,
        completed,
        planned,
    }
}
