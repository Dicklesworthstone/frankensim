//! Causal classical absorption using the original liquid or real-gas state.
//!
//! The existing observer owns minimum-phase design, energy-tail truncation,
//! linear overlap-add convolution and streaming FIR arithmetic. This adapter
//! supplies the shared material's classical bulk law and an explicit numerical
//! tolerance. It never constructs a substitute ideal-gas state or assigns
//! terrestrial-air molecular relaxation/uncertainty to another fluid.
//!
//! The pressure history is filtered in source time. Geometric spreading and
//! the reported residual receiver delay remain the propagation consumer's
//! responsibility. The in-band minimum-phase realization is an approximation,
//! not a physical dispersion law, full Navier–Stokes solution, or source
//! uncertainty certificate. The fluid is uniform and stationary; source
//! applicability, constant bulk viscosity, and omitted relaxation remain
//! explicit on the supplied [`ClassicalBulkAbsorption`] owner.

use super::{
    AcousticRealizeError, CLASSICAL_VALIDITY_OMEGA_TAU, FftC64, ObserverResponse,
    convolve_pressure_history, design_observer_response_with_budget, invalid, minimum_phase,
    process_fir,
};
use fs_material::fluid_acoustics::thermoviscous::absorption::ClassicalBulkAbsorption;

/// Numerical response and physical propagation bookkeeping for one frozen fluid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FluidPropagationReceipt {
    /// Original coherent fluid coefficients and explicitly supplied bulk viscosity.
    pub absorption: ClassicalBulkAbsorption,
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
    /// Worst `omega*tau` in the sampled band; at most the shared weak-loss limit.
    pub maximum_weak_loss_parameter: f64,
    /// DC group delay of the actual retained FIR [s].
    pub low_frequency_delay_s: f64,
    /// Numerical Nyquist alignment and precursor padding [s].
    pub latency_s: f64,
    /// Remaining physical travel delay [s]: `range/c - low_frequency_delay_s`.
    /// May be negative; a causal receiver must add common numerical padding
    /// rather than silently omit this difference or count travel twice.
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
        let (response, receipt) = response(absorption, range_m, dt, maximum_truncated_energy)?;
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
    if pressure_pa.iter().any(|value| !value.is_finite()) {
        return Err(invalid("fluid observer pressure must be finite"));
    }
    let (response, receipt) = response(absorption, range_m, dt, maximum_truncated_energy)?;
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
        return Err(invalid(
            "fluid observer Nyquist band exceeds classical weak-loss applicability",
        ));
    }
    let response =
        design_observer_response_with_budget(range_m, rate, Some(maximum_truncated_energy), |n| {
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

#[cfg(test)]
mod tests;
