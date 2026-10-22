//! Finite-exposure measurements of the original thermal trajectory.
//!
//! A time mean is one measurement, not independent targets at every endpoint.
//! Trapezoidal weights integrate the piecewise-linear reconstruction of the
//! selected spatial readout and divide by the declared exposure duration. The
//! original physical time integrator and source sampling are unchanged. This is
//! a discrete quadrature, not a continuum-error bound or a sensor-response model.

use super::{ScanObservable, ScanSensitivityError, ScanStep, finite, invalid, poll, zeros};
use super::loss::ScanResidualLoss;
use super::spatial::{ScanObservationSpec, sealed};
use crate::ConductionMesh;
use fs_exec::Cx;

/// An instantaneous measurement or one time-averaged exposure. Both variants can
/// coexist in the same observation vector, with vertex or spatial readouts and
/// quadratic/Huber loss. The target and scale keep the inner observation's units.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanTemporalObservation<O> {
    observation: O,
    window: Option<Window>,
}

#[derive(Debug, Clone, PartialEq)]
struct Window {
    first: usize,
    times_s: Vec<f64>,
    weights: Vec<f64>,
}

impl<O: ScanObservationSpec> ScanTemporalObservation<O> {
    /// Preserve an existing measurement, including its loss and temporal layout.
    pub fn instantaneous(observation: O) -> Self { Self { observation, window: None } }

    /// Average from first_endpoint through observation.endpoint(), inclusive.
    /// The window starts at an accepted endpoint, NOT the run's initial state.
    /// At least two finite increasing times are required. All interior endpoints
    /// contribute; the exact selected grid is bound and rechecked before physics.
    /// Changing/refining that grid requires rebuilding the observation.
    ///
    /// max_samples bounds this row's retained clock/weights. Evaluation also
    /// charges every endpoint contribution against the shared max_observations
    /// allowance. Nested time averages refuse rather than averaging an average.
    pub fn trapezoidal_mean(cx: &Cx<'_>, observation: O, endpoints_s: &[f64],
        first_endpoint: usize, max_samples: usize) -> Result<Self, ScanSensitivityError>
    {
        poll(cx, 0)?;
        let last = observation.endpoint();
        let count = last.checked_sub(first_endpoint).and_then(|n| n.checked_add(1))
            .filter(|&n| n >= 2 && n <= max_samples)
            .ok_or_else(|| invalid("time mean needs at least two endpoints within its sample allowance"))?;
        if last >= endpoints_s.len() || observation.first_endpoint() != last {
            return Err(invalid("time mean requires an existing terminal endpoint and an instantaneous inner readout"));
        }
        observation.validate_times(cx, endpoints_s)?;
        let selected = &endpoints_s[first_endpoint..=last];
        let duration = selected[count-1] - selected[0];
        if !duration.is_finite() || duration <= 0.0 {
            return Err(invalid("time mean exposure must have a representable positive duration"));
        }
        let mut weights = zeros(count)?;
        for (i, pair) in selected.windows(2).enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            let dt = pair[1] - pair[0];
            if !pair[0].is_finite() || !pair[1].is_finite() || !dt.is_finite() || dt <= 0.0 {
                return Err(invalid("time mean samples must have finite strictly increasing seconds"));
            }
            let half = (dt / duration) * 0.5;
            if !half.is_finite() || half <= 0.0 {
                return Err(invalid("positive time-mean contribution is not representable"));
            }
            weights[i] = finite(weights[i] + half)?;
            weights[i+1] = finite(weights[i+1] + half)?;
        }
        // Only a common roundoff-scale renormalization: the relative weights
        // remain the declared nonuniform trapezoidal exposure, not sample counts.
        let mut sum = 0.0;
        for (i, &weight) in weights.iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            sum = finite(sum + weight)?;
        }
        for (i, weight) in weights.iter_mut().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            *weight = finite(*weight / sum)?;
            if *weight <= 0.0 { return Err(invalid("positive time-mean weight underflowed")); }
        }
        let mut times_s = zeros(count)?;
        for (i, (time, &original)) in times_s.iter_mut().zip(selected).enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            *time = original;
        }
        poll(cx, count)?;
        Ok(Self { observation, window: Some(Window { first: first_endpoint, times_s, weights }) })
    }

    pub fn observation(&self) -> &O { &self.observation }
    /// Exact bound samples and dimensionless normalized quadrature coefficients.
    /// None denotes the unchanged inner observation, not a zero-length average.
    pub fn time_weights(&self) -> Option<(&[f64], &[f64])> {
        self.window.as_ref().map(|w| (w.times_s.as_slice(), w.weights.as_slice()))
    }
}

impl<O: ScanObservationSpec> sealed::Sealed for ScanTemporalObservation<O> {}
impl<O: ScanObservationSpec> ScanObservationSpec for ScanTemporalObservation<O> {
    fn endpoint(&self) -> usize { self.observation.endpoint() }
    fn first_endpoint(&self) -> usize {
        self.window.as_ref().map_or_else(|| self.observation.first_endpoint(), |w| w.first)
    }
    fn time_weight(&self, endpoint: usize) -> f64 {
        match &self.window {
            None => self.observation.time_weight(endpoint),
            Some(w) => endpoint.checked_sub(w.first).and_then(|i| w.weights.get(i)).copied().unwrap_or(0.0),
        }
    }
    fn validate_times(&self, cx: &Cx<'_>, endpoints: &[f64]) -> Result<(), ScanSensitivityError> {
        self.observation.validate_times(cx, endpoints)?;
        if let Some(w) = &self.window {
            let actual = endpoints.get(w.first..=self.endpoint())
                .ok_or_else(|| invalid("time mean names an absent endpoint window"))?;
            if actual.len() != w.times_s.len() { return Err(invalid("time mean grid shape changed")); }
            for (i, (&a, &b)) in actual.iter().zip(&w.times_s).enumerate() {
                if i % 256 == 0 { poll(cx, i)?; }
                if a.to_bits() != b.to_bits() { return Err(invalid("time mean must be rebuilt for the selected physical grid")); }
            }
        }
        poll(cx, self.endpoint())
    }
    fn observable(&self) -> ScanObservable { self.observation.observable() }
    fn target(&self) -> f64 { self.observation.target() }
    fn scale(&self) -> f64 { self.observation.scale() }
    fn loss(&self) -> ScanResidualLoss { self.observation.loss() }
    fn validate_readout(&self, mesh: &ConductionMesh) -> Result<(), ScanSensitivityError> {
        self.observation.validate_readout(mesh)
    }
    fn prediction(&self, cx: &Cx<'_>, report: &ScanStep) -> Result<f64, ScanSensitivityError> {
        self.observation.prediction(cx, report)
    }
    fn add_seed(&self, cx: &Cx<'_>, seed: &mut [f64], value: f64) -> Result<(), ScanSensitivityError> {
        self.observation.add_seed(cx, seed, value)
    }
}

// One scoring implementation for full-storage, replayed and value-only runs.
// The residual/loss is formed once, after the entire measurement is accumulated.
pub(super) struct ObservationValues {
    pub(super) predictions: Vec<f64>,
    pub(super) residuals: Vec<f64>,
    pub(super) value: f64,
}
impl ObservationValues {
    pub(super) fn new(count: usize) -> Result<Self, ScanSensitivityError> {
        Ok(Self { predictions: zeros(count)?, residuals: zeros(count)?, value: 0.0 })
    }
    pub(super) fn record<O: ScanObservationSpec>(&mut self, observation: &O,
        row: usize, endpoint: usize, prediction: f64) -> Result<(), ScanSensitivityError>
    {
        let contribution = finite(finite(prediction)? * observation.time_weight(endpoint))?;
        self.predictions[row] = if endpoint == observation.first_endpoint() { contribution }
            else { finite(self.predictions[row] + contribution)? };
        if endpoint == observation.endpoint() {
            let residual = finite((self.predictions[row] - observation.target()) / observation.scale())?;
            self.value = finite(self.value + observation.loss().value_derivative(residual)?.0)?;
            self.residuals[row] = residual;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{ScanObservation, loss::ScanObservationWithLoss};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
        let gate = CancelGate::new();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            f(&gate, &Cx::new(&gate, arena,
                StreamKey { seed: 76, kernel_id: 51, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic));
        });
    }
    fn row() -> ScanObservation {
        ScanObservation { endpoint: 2, vertex: 0,
            observable: ScanObservable::SpecificEnthalpy, target: 7.0, scale: 2.0 }
    }

    #[test]
    fn g1_nonuniform_exposure_scores_one_average_and_retains_loss_wrapper_order() {
        with_cx(|_, cx| {
            let times = [1.0, 2.0, 4.0];
            let window = ScanTemporalObservation::trapezoidal_mean(cx, row(), &times, 0, 3).unwrap();
            let (_, weights) = window.time_weights().unwrap();
            for (&a, b) in weights.iter().zip([1.0/6.0, 0.5, 1.0/3.0]) {
                assert!((a-b).abs() < 1e-15);
            }
            let mut values = ObservationValues::new(1).unwrap();
            for (i, &time) in times.iter().enumerate() {
                values.record(&window, 0, i, 2.0*time + 1.0).unwrap();
                if i < 2 { assert_eq!(values.value, 0.0); }
            }
            assert!((values.predictions[0]-6.0).abs() < 1e-14);
            assert!((values.value-0.125).abs() < 1e-14);
            let loss = ScanResidualLoss::Huber { threshold: 0.25 };
            let inside = ScanTemporalObservation::trapezoidal_mean(cx,
                ScanObservationWithLoss::new(row(), loss).unwrap(), &times, 0, 3).unwrap();
            let outside = ScanObservationWithLoss::new(window, loss).unwrap();
            let (mut a, mut b) = (ObservationValues::new(1).unwrap(), ObservationValues::new(1).unwrap());
            for (i, &time) in times.iter().enumerate() {
                a.record(&inside, 0, i, 2.0*time+1.0).unwrap();
                b.record(&outside, 0, i, 2.0*time+1.0).unwrap();
            }
            assert_eq!(a.predictions, b.predictions);
            assert_eq!(a.residuals, b.residuals);
            assert_eq!(a.value, b.value);
            assert!((a.value-0.09375).abs() < 1e-14);
            assert_eq!(outside.first_endpoint(), 0);
            assert_eq!(outside.time_weight(3), 0.0);
        });
    }

    #[test]
    fn g0_windows_bind_actual_grid_and_refuse_underbudget_nested_or_cancelled_requests() {
        with_cx(|gate, cx| {
            for (times, first, cap) in [(vec![1.,2.,4.], 0, 2), (vec![1.,2.,4.], 2, 3),
                (vec![1.,2.,2.], 0, 3), (vec![1.,2.,f64::INFINITY], 0, 3)] {
                assert!(ScanTemporalObservation::trapezoidal_mean(cx, row(), &times, first, cap).is_err());
            }
            let window = ScanTemporalObservation::trapezoidal_mean(cx, row(), &[1.,2.,4.], 0, 3).unwrap();
            assert!(window.validate_times(cx, &[1.,2.5,4.]).is_err());
            assert!(window.validate_times(cx, &[1.,2.]).is_err());
            assert!(ScanTemporalObservation::trapezoidal_mean(cx, window, &[1.,2.,4.], 0, 3).is_err());
            gate.request();
            assert!(ScanTemporalObservation::trapezoidal_mean(cx, row(), &[1.,2.,4.], 0, 3).is_err());
        });
    }

    #[test]
    fn g3_instantaneous_wrapping_preserves_original_prediction_and_loss_bits() {
        let observation = row();
        let instant = ScanTemporalObservation::instantaneous(observation);
        for prediction in [-0.0_f64, 0.0, 1.25, 7.0, 1e100] {
            let mut values = ObservationValues::new(1).unwrap();
            values.record(&instant, 0, 2, prediction).unwrap();
            let residual = (prediction-observation.target)/observation.scale;
            assert_eq!(values.predictions[0].to_bits(), prediction.to_bits());
            assert_eq!(values.residuals[0].to_bits(), residual.to_bits());
            assert_eq!(values.value.to_bits(), ((0.5*residual)*residual).to_bits());
        }
    }
}
