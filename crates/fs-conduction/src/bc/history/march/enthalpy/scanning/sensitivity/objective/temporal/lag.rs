//! Analytic first-order sensor response to a piecewise-linear spatial readout.
use super::{ScanTemporalObservation, Window, ScanObservationSpec, ScanObservable,
    ScanSensitivityError, finite, invalid, poll, zeros};
use fs_exec::Cx;

impl<O: ScanObservationSpec> ScanTemporalObservation<O> {
    /// Sample tau*dy/dt + y = u at observation.endpoint(), with u the ORIGINAL
    /// spatial/vertex readout reconstructed linearly between accepted endpoints.
    /// initial_reading is y at first_endpoint, in the observable's units; it is
    /// fixed measurement metadata, NOT the simulated material's initial state.
    ///
    /// The time constant is positive seconds. Both it and the initial reading
    /// remain fixed during fitting. The start must be an accepted endpoint, not
    /// the run's initial state. At least two samples are required. This is a
    /// filtered terminal reading, not a time average, delay or coupled sensor.
    /// The physical solver, hard constraints and constitutive derivatives stay
    /// unchanged. Loss is applied once after this complete affine readout.
    ///
    /// Analytic segment integration is evaluated in binary64. Every positive
    /// coefficient and the surviving initial-state weight must be representable;
    /// overflow/underflow refuses rather than silently dropping sensor history.
    /// Nested temporal operators refuse. Grid/sample/cancellation admission is
    /// the same as for time means, including exact-grid binding on later use.
    #[allow(clippy::too_many_arguments)]
    pub fn first_order_lag(cx: &Cx<'_>, observation: O, endpoints_s: &[f64],
        first_endpoint: usize, time_constant_s: f64, initial_reading: f64,
        max_samples: usize) -> Result<Self, ScanSensitivityError>
    {
        poll(cx, 0)?;
        let last = observation.endpoint();
        let count = last.checked_sub(first_endpoint).and_then(|n| n.checked_add(1))
            .filter(|&n| n >= 2 && n <= max_samples)
            .ok_or_else(|| invalid("sensor lag needs at least two bounded endpoint samples"))?;
        if last >= endpoints_s.len() || observation.first_endpoint() != last
            || observation.prediction_offset() != 0.0
            || !time_constant_s.is_finite() || time_constant_s <= 0.0
            || !initial_reading.is_finite()
        { return Err(invalid("sensor lag needs an instantaneous readout, positive finite seconds and a finite initial reading")); }
        if observation.observable() == ScanObservable::LiquidMassFraction
            && !(0.0..=1.0).contains(&initial_reading)
        { return Err(invalid("initial sensor fraction must be dimensionless in [0,1], not percent")); }
        observation.validate_times(cx, endpoints_s)?;
        let selected = &endpoints_s[first_endpoint..=last];
        let mut weights = zeros(count)?;
        let mut times_s = zeros(count)?;
        let mut surviving = 1.0;
        for i in (0..count-1).rev() {
            poll(cx, i)?;
            let dt = selected[i+1] - selected[i];
            if !selected[i].is_finite() || !selected[i+1].is_finite()
                || !dt.is_finite() || dt <= 0.0
            { return Err(invalid("sensor lag needs finite strictly increasing sample times")); }
            let (decay, left, right) = coefficients(dt / time_constant_s)?;
            weights[i] = finite(weights[i] + positive(surviving * left)?)?;
            weights[i+1] = finite(weights[i+1] + positive(surviving * right)?)?;
            surviving = positive(surviving * decay)?;
        }
        let offset = finite(surviving * initial_reading)?;
        if initial_reading != 0.0 && offset == 0.0 {
            return Err(invalid("initial sensor contribution underflows"));
        }
        for (i, (time, &original)) in times_s.iter_mut().zip(selected).enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            *time = original;
        }
        poll(cx, count)?;
        Ok(Self { observation, window: Some(Window { first: first_endpoint,
            times_s, weights, offset }) })
    }
}

// y1 = a*y0 + b0*u0 + b1*u1, r=dt/tau, q=1-exp(-r).
// b1=1-q/r, b0=q-b1. For short intervals the first expression cancels;
// sum r/2-r^2/6+... instead. Ten terms at r<=1/8 have <2.5e-19 truncation.
// Deliberately no normalization: b0+b1=1-a, not one.
fn coefficients(r: f64) -> Result<(f64, f64, f64), ScanSensitivityError> {
    positive(r)?;
    let decay = positive((-r).exp())?;
    let q = -(-r).exp_m1();
    let right = if r <= 0.125 {
        let mut term = r * 0.5;
        let mut sum = term;
        for degree in 2..=10 {
            term *= -r / f64::from(degree + 1);
            sum += term;
        }
        sum
    } else { 1.0 - q / r };
    Ok((decay, positive(q - right)?, positive(right)?))
}

fn positive(value: f64) -> Result<f64, ScanSensitivityError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(invalid("sensor response coefficient overflows or underflows; no history is discarded"));
    }
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::ObservationValues;
    use super::super::super::ScanObservation;
    use super::super::super::loss::{ScanObservationWithLoss, ScanResidualLoss};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

    fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
        let gate = CancelGate::new();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| f(&gate,
            &Cx::new(&gate, arena, StreamKey {seed:76,kernel_id:51,tile:0,iteration:0},
                Budget::INFINITE, ExecMode::Deterministic)));
    }
    fn row() -> ScanObservation {
        ScanObservation {endpoint:2,vertex:0,observable:ScanObservable::SpecificEnthalpy,target:7.0,scale:2.0}
    }

    #[test]
    fn g1_lag_matches_constant_and_affine_input_with_nonuniform_clocks() {
        with_cx(|_, cx| {
            let times = [1.0, 1.125, 4.0];
            for tau in [0.1, 0.5, 1.0, 100.0] {
                for slope in [0.0, 2.0] {
                    let lag = ScanTemporalObservation::first_order_lag(cx, row(), &times, 0, tau, 6.0, 3).unwrap();
                    let mut values = ObservationValues::new(1).unwrap();
                    for (i, &t) in times.iter().enumerate() { values.record(&lag, 0, i, 1.0+slope*t).unwrap(); }
                    let expected = 1.0+slope*4.0-slope*tau
                        +(6.0-(1.0+slope)+slope*tau)*(-3.0/tau).exp();
                    assert!((values.predictions[0]-expected).abs() < 1e-12);
                    let sum: f64 = lag.time_weights().unwrap().1.iter().sum();
                    assert!((sum+(-3.0/tau).exp()-1.0).abs() < 1e-14);
                }
            }
            let (_, left, right) = coefficients(1e-14).unwrap();
            assert!((left/5e-15-1.0).abs() < 1e-13);
            assert!((right/5e-15-1.0).abs() < 1e-13);
        });
    }

    #[test]
    fn g3_robust_wrapper_order_keeps_initial_sensor_memory_and_single_loss() {
        with_cx(|_, cx| {
            let loss = ScanResidualLoss::Huber { threshold: 0.25 };
            let times = [1.0, 2.0, 4.0];
            let inner = ScanTemporalObservation::first_order_lag(cx,
                ScanObservationWithLoss::new(row(), loss).unwrap(), &times, 0, 2.0, 6.0, 3).unwrap();
            let outer = ScanObservationWithLoss::new(ScanTemporalObservation::first_order_lag(cx,
                row(), &times, 0, 2.0, 6.0, 3).unwrap(), loss).unwrap();
            let (mut a, mut b) = (ObservationValues::new(1).unwrap(), ObservationValues::new(1).unwrap());
            for (i, &t) in times.iter().enumerate() {
                a.record(&inner, 0, i, t).unwrap(); b.record(&outer, 0, i, t).unwrap();
            }
            assert_eq!(a.predictions, b.predictions);
            assert_eq!(a.residuals, b.residuals);
            assert_eq!(a.value.to_bits(), b.value.to_bits());
            assert_eq!(inner.prediction_offset().to_bits(), outer.prediction_offset().to_bits());
        });
    }

    #[test]
    fn g0_lag_refuses_invalid_policy_grid_nesting_underflow_and_cancellation() {
        with_cx(|gate, cx| {
            let times = [1.0, 2.0, 4.0];
            for tau in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e-300] {
                assert!(ScanTemporalObservation::first_order_lag(cx, row(), &times, 0, tau, 6.0, 3).is_err());
            }
            assert!(ScanTemporalObservation::first_order_lag(cx, row(), &times, 0, 1.0, f64::NAN, 3).is_err());
            assert!(ScanTemporalObservation::first_order_lag(cx, row(), &times, 0, 1.0, 6.0, 2).is_err());
            let lag = ScanTemporalObservation::first_order_lag(cx, row(), &times, 0, 1.0, 6.0, 3).unwrap();
            assert!(lag.validate_times(cx, &[1.0, 2.5, 4.0]).is_err());
            assert!(ScanTemporalObservation::first_order_lag(cx, lag, &times, 0, 1.0, 6.0, 3).is_err());
            let mut fraction = row(); fraction.observable = ScanObservable::LiquidMassFraction;
            assert!(ScanTemporalObservation::first_order_lag(cx, fraction, &times, 0, 1.0, 60.0, 3).is_err());
            gate.request();
            assert!(ScanTemporalObservation::first_order_lag(cx, row(), &times, 0, 1.0, 6.0, 3).is_err());
        });
    }
}
