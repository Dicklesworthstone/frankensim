//! Bulk absorption and dispersion using the original liquid or real-gas state.
//!
//! The existing observer owns energy-tail truncation, numerical precursor
//! padding, linear overlap-add convolution and streaming FIR arithmetic.
//! The classical constructor retains its minimum-phase approximation; explicit
//! thermoviscous selection preserves the original material's complex acoustic
//! wavenumber. Both use the same coherent state and numerical tolerance.
//!
//! The pressure history is filtered in source time. Geometric spreading and
//! the reported residual receiver delay remain the propagation consumer's
//! responsibility. `range/c0` is a phase reference: the parabolic
//! Navier–Stokes–Fourier model does not have a finite wavefront. This is one
//! homogeneous stationary acoustic branch, without boundary thermal modes,
//! molecular spectra or temperature evolution. Source applicability and the
//! explicit constant bulk viscosity remain on [`ClassicalBulkAbsorption`].

use super::{
    AcousticRealizeError, CLASSICAL_VALIDITY_OMEGA_TAU, FftC64, ObserverResponse,
    convolve_pressure_history, design_observer_response_with_budget, invalid, minimum_phase,
    process_fir,
};
use fs_material::fluid_acoustics::thermoviscous::absorption::ClassicalBulkAbsorption;
use fs_math::det;

/// Explicit propagation law for one original thermoviscous material state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FluidPropagationModel {
    /// Leading-order Stokes–Kirchhoff attenuation, with its finite-band
    /// minimum-phase approximation. This preserves the original constructor.
    ClassicalMinimumPhase,
    /// Exact spatial acoustic branch of the linear homogeneous
    /// Navier–Stokes–Fourier equations, retaining attenuation and dispersion.
    /// The observer currently admits the same `omega*tau <= 0.01` band as
    /// the classical path; this is a supported regime, not an error bound.
    ThermoviscousDispersion,
}

/// Numerical response and physical propagation bookkeeping for one frozen fluid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluidPropagationReceipt {
    /// Original coherent fluid coefficients and explicitly supplied bulk viscosity.
    pub absorption: ClassicalBulkAbsorption,
    /// Selected attenuation/dispersion law; never inferred from the fluid name.
    pub model: FluidPropagationModel,
    /// Physical path length [m].
    pub range_m: f64,
    /// Exact reciprocal of the supplied sample interval, without integer rounding.
    pub sample_rate_hz: f64,
    /// Retained causal FIR taps.
    pub taps: usize,
    /// Relative discarded impulse-response energy on the admitted design grid.
    /// This does not bound material/model error or maximum frequency error.
    pub truncated_energy: f64,
    /// Caller-selected numerical relative-energy ceiling; no ISO uncertainty.
    pub maximum_truncated_energy: f64,
    /// Worst `omega*tau` in the sampled band; at most the shared supported limit.
    /// Thermoviscous selection uses the full acoustic root within this band.
    pub maximum_weak_loss_parameter: f64,
    /// DC group delay of the actual retained FIR [s].
    pub low_frequency_delay_s: f64,
    /// Numerical Nyquist alignment and precursor padding [s].
    pub latency_s: f64,
    /// Remaining reference delay [s]: `range/c0 - low_frequency_delay_s`.
    /// May be negative; a causal receiver must add common numerical padding
    /// rather than silently omit this difference or count the delay twice.
    /// The parabolic thermoviscous model asserts no finite propagation front.
    pub receiver_delay_s: f64,
}

/// Prepared, allocation-free streaming observer on an exact sample clock.
#[derive(Clone, Debug)]
pub struct FluidObserverFilter {
    taps: Vec<f64>,
    history: Vec<f64>,
    head: usize,
    receipt: FluidPropagationReceipt,
}

impl FluidObserverFilter {
    /// Design classical bulk absorption over `range_m` on the exact clock
    /// `1/dt`. Requires the weak-loss condition through Nyquist, with
    /// `omega*tau <= CLASSICAL_VALIDITY_OMEGA_TAU`. All material evaluation is
    /// already complete on `absorption`; design performs no material queries.
    /// An explicit energy tolerance lies strictly between zero and one.
    pub fn new(
        absorption: ClassicalBulkAbsorption,
        range_m: f64,
        dt: f64,
        maximum_truncated_energy: f64,
    ) -> Result<Self, AcousticRealizeError> {
        Self::with_model(
            absorption,
            range_m,
            dt,
            maximum_truncated_energy,
            FluidPropagationModel::ClassicalMinimumPhase,
        )
    }

    /// Select the propagation law explicitly on the same exact sample clock.
    /// Thermoviscous propagation uses the original complex acoustic root;
    /// finite-band realization and its reported numerical latency remain
    /// approximations controlled by `maximum_truncated_energy`.
    /// Both models require `omega*tau <= CLASSICAL_VALIDITY_OMEGA_TAU` through
    /// Nyquist and propagate material/branch refusals without changing laws.
    pub fn with_model(
        absorption: ClassicalBulkAbsorption,
        range_m: f64,
        dt: f64,
        maximum_truncated_energy: f64,
        model: FluidPropagationModel,
    ) -> Result<Self, AcousticRealizeError> {
        let (response, receipt) =
            response(absorption, range_m, dt, maximum_truncated_energy, model)?;
        let history = vec![0.0; response.taps.len()];
        Ok(Self {
            taps: response.taps,
            history,
            head: 0,
            receipt,
        })
    }

    /// Immutable original material and numerical response report.
    pub const fn receipt(&self) -> &FluidPropagationReceipt {
        &self.receipt
    }

    /// Causal FIR in earliest-tap-first order, for original path composition.
    pub fn taps(&self) -> &[f64] {
        &self.taps
    }

    /// Filter one source-time pressure sample. No allocation or transcendental
    /// evaluation occurs here. A nonfinite sample or unrepresentable output
    /// refuses without advancing or poisoning the retained history.
    pub fn process(&mut self, pressure_pa: f64) -> Result<f64, AcousticRealizeError> {
        if !pressure_pa.is_finite() {
            return Err(invalid("fluid observer pressure must be finite"));
        }
        let previous_head = self.head;
        let previous_sample = self.history[self.head];
        let output = process_fir(&self.taps, &mut self.history, &mut self.head, pressure_pa);
        if !output.is_finite() {
            self.head = previous_head;
            self.history[self.head] = previous_sample;
            return Err(invalid("fluid observer output is unrepresentable"));
        }
        Ok(output)
    }
}

/// Apply the same causal FIR to an entire pressure history on its exact clock.
/// Zero prehistory is assumed, retained samples are unshifted source time,
/// and the tail beyond the supplied slice is discarded. Linear convolution
/// prevents late pressure impulses from wrapping into the record's start.
/// All admission/design/numeric refusals leave the supplied history unchanged.
pub fn propagate_pressure_history(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    absorption: ClassicalBulkAbsorption,
    maximum_truncated_energy: f64,
) -> Result<FluidPropagationReceipt, AcousticRealizeError> {
    propagate_pressure_history_with_model(
        pressure_pa,
        dt,
        range_m,
        absorption,
        maximum_truncated_energy,
        FluidPropagationModel::ClassicalMinimumPhase,
    )
}

/// Apply the explicitly selected material propagation to an entire history.
/// Shares the streaming response and the original linear convolution owner,
/// including zero prehistory, discarded output tail, exact clock and atomic
/// refusal. The receipt records physical-model selection and numerical delay.
pub fn propagate_pressure_history_with_model(
    pressure_pa: &mut [f64],
    dt: f64,
    range_m: f64,
    absorption: ClassicalBulkAbsorption,
    maximum_truncated_energy: f64,
    model: FluidPropagationModel,
) -> Result<FluidPropagationReceipt, AcousticRealizeError> {
    if pressure_pa.iter().any(|value| !value.is_finite()) {
        return Err(invalid("fluid observer pressure must be finite"));
    }
    let (response, receipt) =
        response(absorption, range_m, dt, maximum_truncated_energy, model)?;
    if range_m != 0.0 {
        convolve_pressure_history(pressure_pa, &response.taps)?;
    }
    Ok(receipt)
}

fn response(
    absorption: ClassicalBulkAbsorption,
    range_m: f64,
    dt: f64,
    maximum_truncated_energy: f64,
    model: FluidPropagationModel,
) -> Result<(ObserverResponse, FluidPropagationReceipt), AcousticRealizeError> {
    if !dt.is_finite() || dt <= 0.0 || !(1.0 / dt).is_finite() {
        return Err(invalid(
            "fluid observer needs a positive representable sample clock",
        ));
    }
    let rate = 1.0 / dt;
    let maximum_weak_loss_parameter = absorption
        .weak_loss_parameter(core::f64::consts::PI * rate)
        .map_err(|_| invalid("fluid observer weak-loss parameter is unrepresentable"))?;
    if maximum_weak_loss_parameter > CLASSICAL_VALIDITY_OMEGA_TAU {
        return Err(invalid(match model {
            FluidPropagationModel::ClassicalMinimumPhase => {
                "fluid observer Nyquist band exceeds classical weak-loss applicability"
            }
            FluidPropagationModel::ThermoviscousDispersion => {
                "fluid observer Nyquist band exceeds its supported thermoviscous regime"
            }
        }));
    }
    let response =
        design_observer_response_with_budget(range_m, rate, Some(maximum_truncated_energy), |n| {
            if model == FluidPropagationModel::ThermoviscousDispersion {
                return (0..=n / 2)
                    .map(|k| {
                        let omega = core::f64::consts::TAU * k as f64 * rate / n as f64;
                        thermoviscous_transfer(&absorption, omega, range_m)
                    })
                    .collect();
            }
            let log_magnitude: Result<Vec<f64>, AcousticRealizeError> = (0..=n / 2)
                .map(|k| {
                    let omega = core::f64::consts::TAU * k as f64 * rate / n as f64;
                    let alpha = absorption
                        .absorption_np_per_m(omega)
                        .map_err(|_| invalid("fluid observer absorption is unrepresentable"))?;
                    let value = -alpha * range_m;
                    if !value.is_finite() {
                        return Err(invalid(
                            "fluid observer attenuation over its path is unrepresentable",
                        ));
                    }
                    Ok(value)
                })
                .collect();
            let bins: Vec<FftC64> = minimum_phase(&log_magnitude?, n);
            Ok(bins)
        })?;
    let receiver_delay_s =
        range_m / absorption.state().acoustic().sound_speed_m_s() - response.low_frequency_delay_s;
    if !receiver_delay_s.is_finite() {
        return Err(invalid("fluid observer receiver delay is unrepresentable"));
    }
    let receipt = FluidPropagationReceipt {
        absorption,
        model,
        range_m,
        sample_rate_hz: rate,
        taps: response.taps.len(),
        truncated_energy: response.truncated_energy,
        maximum_truncated_energy,
        maximum_weak_loss_parameter,
        low_frequency_delay_s: response.low_frequency_delay_s,
        latency_s: response.latency_s,
        receiver_delay_s,
    };
    Ok((response, receipt))
}

// The material uses exp(i*k*r-i*omega*t); the FFT synthesizes exp(+i*omega*t).
// Conjugate the physical transfer and remove only the c0 reference delay.
// The common response owner supplies and reports any numerical precursor pad.
fn thermoviscous_transfer(
    absorption: &ClassicalBulkAbsorption,
    omega: f64,
    range_m: f64,
) -> Result<FftC64, AcousticRealizeError> {
    let wave = absorption
        .thermoviscous_wavenumber(omega)
        .map_err(|_| invalid("fluid thermoviscous acoustic branch could not be evaluated"))?;
    let log_magnitude = -wave.im * range_m;
    let phase = -(wave.re - omega / absorption.state().acoustic().sound_speed_m_s()) * range_m;
    if !log_magnitude.is_finite() || !phase.is_finite() {
        return Err(invalid("fluid thermoviscous path transfer is unrepresentable"));
    }
    let magnitude = det::exp(log_magnitude);
    let (sin, cos) = det::sin_cos(phase);
    Ok(FftC64::new(magnitude * cos, magnitude * sin))
}

#[cfg(test)]
mod tests;
