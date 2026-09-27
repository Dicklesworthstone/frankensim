//! Noise-scaled sensor objectives for the production transient calibration path.
//!
//! Each observed scalar contributes rho((prediction - reading) / sigma).
//! Quadratic loss is half weighted squared error. Huber loss is quadratic
//! inside an explicit dimensionless threshold and linear outside it; see
//! Huber (1964), and scipy.special.huber for the piecewise convention.
//! Sigma and the loss threshold are fixed data, not inferred parameters.
//! This diagonal weighting is not a correlated likelihood, a posterior, an
//! outlier classifier, or an identifiability/physical-validation certificate.

use super::{OdeVjp, SampleObjective, TransientFamily, TransientModel};
use std::sync::Arc;

/// Joint fitting with fixed cross-channel and cross-time error factors.
pub mod correlated;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SensorLoss {
    Quadratic,
    /// Strictly positive finite threshold in standardized-residual units.
    Huber { threshold: f64 },
}
impl SensorLoss {
    fn validate(self) -> Result<(), String> {
        if let Self::Huber { threshold } = self {
            if !threshold.is_finite() || threshold <= 0.0 {
                return Err("Huber threshold must be finite and positive".into());
            }
        }
        Ok(())
    }
    fn value_derivative(self, residual: f64) -> (f64, f64) {
        match self {
            Self::Huber { threshold } if residual.abs() > threshold =>
                (threshold * (residual.abs() - 0.5 * threshold), threshold.copysign(residual)),
            _ => ((0.5 * residual) * residual, residual),
        }
    }
}

/// One present scalar reading. Missing/censored readings must not be fabricated
/// as zeros: omit missing rows explicitly upstream; censored likelihoods are
/// outside this adapter. Repeated timestamps remain distinct objective terms.
#[derive(Debug, Clone, PartialEq)]
pub struct SensorReading {
    time: f64,
    channel: u64,
    value: f64,
    sigma: f64,
    loss: SensorLoss,
}
impl SensorReading {
    pub fn new(time: f64, channel: u64, value: f64, sigma: f64, loss: SensorLoss)
        -> Result<Self, String>
    {
        loss.validate()?;
        if !time.is_finite() || !value.is_finite() || !sigma.is_finite() || sigma <= 0.0 {
            return Err("reading time/value must be finite and sigma finite and positive".into());
        }
        Ok(Self { time, channel, value, sigma, loss })
    }
    pub fn time(&self) -> f64 { self.time }
    pub fn channel(&self) -> u64 { self.channel }
    pub fn value(&self) -> f64 { self.value }
    pub fn sigma(&self) -> f64 { self.sigma }
    pub fn loss(&self) -> SensorLoss { self.loss }

    /// Return (loss, dloss/dprediction), including BOTH factors of sigma in
    /// the quadratic derivative. No clipping of the prediction or observation.
    pub fn score(&self, prediction: f64) -> Result<(f64, f64), String> {
        if !prediction.is_finite() { return Err("non-finite sensor prediction".into()); }
        let difference = prediction - self.value;
        let residual = if difference.is_finite() { difference / self.sigma }
            else { prediction / self.sigma - self.value / self.sigma };
        if !residual.is_finite() { return Err("non-finite standardized sensor residual".into()); }
        let (value, derivative) = self.loss.value_derivative(residual);
        let derivative = derivative / self.sigma;
        if !value.is_finite() || !derivative.is_finite() {
            return Err("non-finite sensor objective or derivative".into());
        }
        Ok((value, derivative))
    }
}

/// Fixed data shared across parameter trials without copying every reading.
#[derive(Debug, Clone)]
pub struct SensorData {
    readings: Arc<[SensorReading]>,
    times: Arc<[f64]>,
}
impl SensorData {
    pub fn new(readings: &[SensorReading], max_readings: usize) -> Result<Self, String> {
        if readings.is_empty() || readings.len() > max_readings {
            return Err("nonempty sensor data must fit the reading cap".into());
        }
        if readings.windows(2).any(|pair| pair[0].time > pair[1].time) {
            return Err("sensor readings must be in nondecreasing time order".into());
        }
        Ok(Self { readings: Arc::from(readings),
            times: readings.iter().map(|r| r.time).collect::<Vec<_>>().into() })
    }
    pub fn readings(&self) -> &[SensorReading] { &self.readings }
    pub fn times(&self) -> &[f64] { &self.times }
}

/// Dynamics and a scalar sensor map at one immutable decision point. All
/// derivatives use the decision coordinates declared by SensorFamily::bounds.
/// The sensor VJP receives the already noise/loss-scaled seed; do not apply
/// the loss or noise scale a second time. All outputs must be overwritten.
pub trait SensorModel: OdeVjp {
    fn initial_values(&self) -> &[f64];
    fn initial_vjp(&self, initial_bar: &[f64], parameter_bar: &mut [f64]) -> Result<(), String>;
    fn predict(&self, channel: u64, time: f64, state: &[f64]) -> Result<f64, String>;
    #[allow(clippy::too_many_arguments)]
    fn prediction_vjp(
        &self, channel: u64, time: f64, state: &[f64], seed: f64,
        state_bar: &mut [f64], parameter_bar: &mut [f64],
    ) -> Result<(), String>;
}
pub trait SensorFamily {
    type Model: SensorModel;
    fn bounds(&self) -> &[[f64; 2]];
    fn instantiate(&self, point: &[f64]) -> Result<Self::Model, String>;
}

/// Plug noise-scaled measurements directly into TransientStudy (or a campaign).
/// Data and family are retained immutably; no parameter-dependent reweighting.
pub struct ObservedFamily<F> { family: F, data: SensorData }
impl<F> ObservedFamily<F> {
    pub fn new(family: F, data: SensorData) -> Self { Self { family, data } }
    pub fn data(&self) -> &SensorData { &self.data }
}
pub struct ObservedModel<M> { model: M, data: SensorData }
impl<F: SensorFamily> TransientFamily for ObservedFamily<F> {
    type Model = ObservedModel<F::Model>;
    fn bounds(&self) -> &[[f64; 2]] { self.family.bounds() }
    fn sample_times(&self) -> &[f64] { self.data.times() }
    fn instantiate(&self, point: &[f64]) -> Result<Self::Model, String> {
        Ok(ObservedModel { model: self.family.instantiate(point)?, data: self.data.clone() })
    }
}
impl<M: SensorModel> OdeVjp for ObservedModel<M> {
    fn dimension(&self) -> usize { self.model.dimension() }
    fn parameter_count(&self) -> usize { self.model.parameter_count() }
    fn rhs(&self, time: f64, state: &[f64], out: &mut [f64]) { self.model.rhs(time, state, out); }
    fn rhs_vjp(&self, time: f64, state: &[f64], seed: &[f64], x: &mut [f64], p: &mut [f64])
        -> Result<(), String>
    { self.model.rhs_vjp(time, state, seed, x, p) }
}
impl<M: SensorModel> SampleObjective for ObservedModel<M> {
    fn evaluate(&self, sample: usize, time: f64, state: &[f64], x: &mut [f64], p: &mut [f64])
        -> Result<f64, String>
    {
        let reading = self.data.readings.get(sample).ok_or("unknown sensor sample")?;
        if time != reading.time || state.len() != self.dimension()
            || x.len() != self.dimension() || p.len() != self.parameter_count()
        { return Err("sensor sample time or derivative dimensions do not match".into()); }
        let prediction = self.model.predict(reading.channel, time, state)?;
        let (value, seed) = reading.score(prediction)?;
        x.fill(f64::NAN); p.fill(f64::NAN);
        self.model.prediction_vjp(reading.channel, time, state, seed, x, p)?;
        if x.iter().chain(p.iter()).any(|v| !v.is_finite()) {
            return Err("non-finite or unwritten sensor derivative".into());
        }
        Ok(value)
    }
}
impl<M: SensorModel> TransientModel for ObservedModel<M> {
    fn initial_values(&self) -> &[f64] { self.model.initial_values() }
    fn initial_vjp(&self, initial_bar: &[f64], parameter_bar: &mut [f64]) -> Result<(), String> {
        self.model.initial_vjp(initial_bar, parameter_bar)
    }
}

#[cfg(test)]
mod tests;
