//! Deterministic serial ensemble square-root analysis (Whitaker--Hamill, 2002).
//!
//! The state is an explicit, caller-supplied ensemble, not a dense covariance.
//! A scalar analysis costs O(state dimension * members) and allocates no
//! state-by-state matrix. The mean uses the ordinary sample Kalman gain; the
//! anomalies use its reduced square-root gain. Observations are not randomly
//! perturbed. The caller owns ensemble generation, units, model identity and
//! independent measurement-error assumptions. Nonlinear observation maps and
//! localization give approximate analyses, not exact Bayesian posteriors.
//!
//! IDs must increase within a timestamp; a successfully consumed reading cannot
//! be reapplied. Failed/cancelled calls leave every ensemble bit and the cursor
//! unchanged. Callback attempts remain charged in the caller-owned control.
//! A cancellation closure can bridge `fs_exec::Cx::checkpoint`; long model
//! callbacks must provide their own cancellation. No cross-ISA claim is made.

#[path = "ensemble/forecast.rs"]
pub mod forecast;

#[derive(Debug, Clone, PartialEq)]
pub enum EnsembleError {
    Invalid(&'static str),
    WorkspaceLimit { required: usize, limit: usize },
    ModelCallLimit,
    Allocation,
    NonFinite(&'static str),
    Model { member: usize, message: String },
    ObservationOrder,
    Cancelled,
}
impl std::fmt::Display for EnsembleError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ensemble estimation failed: {self:?}")
    }
}
impl std::error::Error for EnsembleError {}

/// Non-cloneable cumulative callback allowance. Numerical workspace excludes
/// the retained input ensemble and model-owned memory. Shape bounds also bound
/// all algebraic loops; callback-internal work needs its own model budget.
#[derive(Debug)]
pub struct EnsembleControl {
    max_calls: usize,
    calls: usize,
    workspace: usize,
}
impl EnsembleControl {
    pub fn new(max_model_calls: usize, max_workspace_components: usize) -> Self {
        Self { max_calls: max_model_calls, calls: 0, workspace: max_workspace_components }
    }
    pub fn model_calls(&self) -> usize { self.calls }
    /// Only raise limits; spent attempts are never reset or refunded.
    pub fn extend(&mut self, max_model_calls: usize, max_workspace_components: usize)
        -> Result<(), EnsembleError>
    {
        if max_model_calls < self.max_calls || max_workspace_components < self.workspace {
            return Err(EnsembleError::Invalid("limits may only increase"));
        }
        self.max_calls = max_model_calls;
        self.workspace = max_workspace_components;
        Ok(())
    }
    fn admit(&self, calls: usize, components: usize) -> Result<(), EnsembleError> {
        if components > self.workspace {
            return Err(EnsembleError::WorkspaceLimit { required: components, limit: self.workspace });
        }
        if calls > self.max_calls - self.calls { return Err(EnsembleError::ModelCallLimit); }
        Ok(())
    }
    fn charge(&mut self) -> Result<(), EnsembleError> {
        if self.calls == self.max_calls { return Err(EnsembleError::ModelCallLimit); }
        self.calls += 1;
        Ok(())
    }
}

/// One observation at the ensemble's current time. Sigma is a positive
/// standard deviation in the observation's units, not a variance or weight.
#[derive(Debug, Clone, Copy)]
pub struct EnsembleObservation {
    pub id: u64,
    pub time: f64,
    pub value: f64,
    pub sigma: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleAnalysis {
    pub observation_id: u64,
    pub predicted_mean: f64,
    pub predicted_std: f64,
    pub innovation_std: f64,
    pub standardized_innovation: f64,
    /// The anomaly-gain multiplier, in [1/2, 1].
    pub square_root_factor: f64,
    /// Zero means no supported sample cross-covariance (or localization zero),
    /// not a successful recovery of an unknown state.
    pub updated_components: usize,
}

/// Member-major finite state. Clone is a numerical checkpoint, including the
/// observation cursor; it does not clone/reset the external work allowance.
#[derive(Debug, Clone, PartialEq)]
pub struct Ensemble {
    time: f64,
    dimension: usize,
    count: usize,
    values: Vec<f64>,
    last_observation: Option<u64>,
}

fn poll<C: FnMut() -> bool>(cancelled: &mut C) -> Result<(), EnsembleError> {
    if cancelled() { Err(EnsembleError::Cancelled) } else { Ok(()) }
}
fn zeros(n: usize) -> Result<Vec<f64>, EnsembleError> {
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| EnsembleError::Allocation)?;
    values.resize(n, 0.0);
    Ok(values)
}
fn finite(x: f64, context: &'static str) -> Result<f64, EnsembleError> {
    if x.is_finite() { Ok(x) } else { Err(EnsembleError::NonFinite(context)) }
}
fn sum(values: impl Iterator<Item = f64>) -> Result<f64, EnsembleError> {
    let (mut total, mut correction) = (0.0_f64, 0.0_f64);
    for x in values {
        let next = finite(total + x, "reduction")?;
        correction += if total.abs() >= x.abs() { (total - next) + x } else { (x - next) + total };
        total = next;
    }
    finite(total + correction, "compensated reduction")
}
fn mean(values: impl Iterator<Item = f64> + Clone, count: usize) -> Result<f64, EnsembleError> {
    // Keep a constant ensemble exactly constant, and avoid summing its common
    // offset N times. Refuse arithmetic outside the representable envelope.
    let anchor = values.clone().next().ok_or(EnsembleError::Invalid("empty mean"))?;
    let denominator = count as f64;
    let correction = sum(values.map(|x| {
        let delta = x - anchor;
        if delta.is_finite() { delta / denominator } else { x / denominator - anchor / denominator }
    }))?;
    finite(anchor + correction, "mean")
}

impl Ensemble {
    /// Supply at least two complete members. No random initialization, implicit
    /// covariance, inflation, or rescaling is invented at this boundary.
    pub fn new(time: f64, dimension: usize, values: &[f64], max_components: usize)
        -> Result<Self, EnsembleError>
    {
        if !time.is_finite() || dimension == 0 || values.len() % dimension != 0
            || values.len() / dimension < 2 || values.len() > max_components
            || values.iter().any(|x| !x.is_finite())
        {
            return Err(EnsembleError::Invalid("finite time and at least two bounded complete members required"));
        }
        let mut owned = zeros(values.len())?;
        owned.copy_from_slice(values);
        Ok(Self { time, dimension, count: values.len() / dimension, values: owned, last_observation: None })
    }
    pub fn time(&self) -> f64 { self.time }
    pub fn dimension(&self) -> usize { self.dimension }
    pub fn member_count(&self) -> usize { self.count }
    pub fn member(&self, index: usize) -> Option<&[f64]> {
        (index < self.count).then(|| &self.values[index*self.dimension..(index+1)*self.dimension])
    }
    pub fn values(&self) -> &[f64] { &self.values }
    pub fn last_observation(&self) -> Option<u64> { self.last_observation }

    /// Analyze one independent-noise reading. A localization entry in [0,1]
    /// tapers that component's gain; `None` means global sample covariance.
    /// Zero localization leaves the component's original bits untouched.
    /// Unlocalized LINEAR observations reproduce the Kalman mean/covariance
    /// update for the supplied sample covariance, up to floating-point error.
    /// For nonlinear maps, predictions are reevaluated on the current ensemble
    /// for every subsequent reading; this is a serial approximation.
    ///
    /// Workspace is one candidate ensemble plus one member-sized vector.
    /// Stable member IDs are indices and must not be reordered between calls.
    pub fn assimilate_scalar<F, C>(
        &mut self, observation: EnsembleObservation, localization: Option<&[f64]>,
        predict: &mut F, control: &mut EnsembleControl, cancelled: &mut C,
    ) -> Result<EnsembleAnalysis, EnsembleError>
    where
        F: FnMut(usize, f64, &[f64]) -> Result<f64, String>,
        C: FnMut() -> bool,
    {
        poll(cancelled)?;
        if !observation.time.is_finite() || observation.time != self.time
            || !observation.value.is_finite() || !observation.sigma.is_finite() || observation.sigma <= 0.0
        { return Err(EnsembleError::Invalid("observation time, value or noise scale")); }
        if self.last_observation.is_some_and(|id| observation.id <= id) {
            return Err(EnsembleError::ObservationOrder);
        }
        if let Some(weights) = localization {
            if weights.len() != self.dimension || weights.iter().any(|w| !w.is_finite() || !(0.0..=1.0).contains(w)) {
                return Err(EnsembleError::Invalid("localization must contain one finite [0,1] weight per component"));
            }
        }
        let workspace = self.values.len().checked_add(self.count)
            .ok_or(EnsembleError::Invalid("workspace extent overflow"))?;
        control.admit(self.count, workspace)?;
        let mut anomalies = zeros(self.count)?;
        for (member, prediction) in anomalies.iter_mut().enumerate() {
            poll(cancelled)?;
            control.charge()?;
            let result = predict(member, self.time, self.member(member).expect("bounded member index"));
            poll(cancelled)?;
            *prediction = finite(result.map_err(|message| EnsembleError::Model { member, message })?, "prediction")?;
        }
        let observed_mean = mean(anomalies.iter().copied(), self.count)?;
        let mut scale = observation.sigma;
        for value in &mut anomalies {
            *value = finite(*value - observed_mean, "observation anomaly")?;
            scale = scale.max(value.abs());
        }
        for value in &mut anomalies { *value /= scale; }
        let divisor = (self.count - 1) as f64;
        let variance = sum(anomalies.iter().map(|x| (x / divisor) * x))?;
        let noise = observation.sigma / scale;
        let total_variance = finite(noise.mul_add(noise, variance), "innovation variance")?;
        if total_variance <= 0.0 { return Err(EnsembleError::NonFinite("vanishing innovation variance")); }
        let difference = observation.value - observed_mean;
        let innovation = finite(if difference.is_finite() { difference / scale }
            else { observation.value / scale - observed_mean / scale }, "scaled innovation")?;
        let alpha = 1.0 / (1.0 + noise / total_variance.sqrt());
        let mut report = EnsembleAnalysis { observation_id: observation.id,
            predicted_mean: observed_mean, predicted_std: finite(scale * variance.sqrt(), "forecast spread")?,
            innovation_std: finite(scale * total_variance.sqrt(), "innovation spread")?,
            standardized_innovation: finite(innovation / total_variance.sqrt(), "standardized innovation")?,
            square_root_factor: alpha, updated_components: 0 };
        let mut candidate = zeros(self.values.len())?;
        candidate.copy_from_slice(&self.values);
        for component in 0..self.dimension {
            poll(cancelled)?;
            let weight = localization.map_or(1.0, |weights| weights[component]);
            if weight == 0.0 || variance == 0.0 { continue; }
            let center = mean(self.values.chunks_exact(self.dimension).map(|x| x[component]), self.count)?;
            let cross = sum(self.values.chunks_exact(self.dimension).zip(&anomalies)
                .map(|(x, y)| ((x[component] - center) / divisor) * y))?;
            let gain = finite(weight * cross / total_variance, "scaled gain")?;
            if gain == 0.0 { continue; }
            let updated_mean = finite(gain.mul_add(innovation, center), "analysis mean")?;
            for (member, dy) in anomalies.iter().enumerate() {
                if member % 256 == 0 { poll(cancelled)?; }
                let offset = member * self.dimension + component;
                let anomaly = finite(self.values[offset] - center, "state anomaly")?;
                candidate[offset] = finite((-alpha * gain).mul_add(*dy, anomaly) + updated_mean, "analysis member")?;
            }
            report.updated_components += 1;
        }
        poll(cancelled)?;
        self.values = candidate;
        self.last_observation = Some(observation.id);
        Ok(report)
    }
}

#[cfg(test)]
#[path = "ensemble/tests.rs"]
mod tests;
