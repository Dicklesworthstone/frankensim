//! Joint amplitude estimation at supplied oscillatory decay poles.

use super::{
    FitDesign, HarmonicFit, HarmonicFitBudget, HarmonicFitError, HarmonicFitSummary, invalid,
};
use crate::C64;
use fs_math::det;

/// A supplied conjugate pole pair `-decay_rate_per_s ± i TAU frequency_hz`.
/// The decay rate acts on amplitude, not energy; zero retains a stationary tone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DampedTonePole {
    /// Positive oscillation frequency in hertz (the damped modal frequency).
    pub frequency_hz: f64,
    /// Nonnegative amplitude decay rate in reciprocal seconds.
    pub decay_rate_per_s: f64,
}

/// Reusable joint fit of DC and exponentially damped real tones:
/// `x(t) = dc + sum_k exp(-decay_k*t) * (a_k*cos(w_k*t) + b_k*sin(w_k*t))`.
///
/// Time zero is the first sample. Returned phasors are `a_k + i b_k` under
/// `exp(-i w_k t)`, in the caller's pole order and physical amplitude units.
/// They are initial amplitudes, not record-averaged amplitudes. The frequency
/// is the oscillation frequency after damping, not an undamped natural frequency.
/// Distinct decay rates may share a frequency, retaining physical double decay.
/// Every column, including DC, participates in the same Householder QR solve.
///
/// The caller supplies all poles. This plan does not identify decay rates or
/// frequencies, select mode count, filter aliasing, or certify a noise model.
/// An incorrect envelope or omitted mode can bias the fitted amplitudes and
/// remain in the sampled residual. Numerical conditioning does not certify
/// that close poles can be separated in noisy data.
///
/// Preparation costs O(N K^2 + K^3) work and O(N K) storage; reusing a prepared
/// plan costs O(N K + K^2) work. `fit_into` allocates nothing. All-zero decay
/// produces the same design and fit bits as the corresponding `ToneFitPlan`.
#[derive(Debug, Clone)]
pub struct DampedToneFitPlan {
    poles: Vec<DampedTonePole>,
    design: FitDesign,
}

impl DampedToneFitPlan {
    /// Prepare a uniformly sampled fit at the supplied distinct poles.
    /// Frequencies must be finite, positive and strictly below Nyquist; decay
    /// rates must be finite and nonnegative. Pole order is preserved, with no
    /// required sorting. Equal frequencies are allowed with distinct decays.
    ///
    /// At least `1 + 2 * poles.len()` samples are needed. The complete sampled
    /// basis must also pass the caller's conditioning floor: a pole that decays
    /// before its oscillation is resolved can be indistinguishable even when
    /// its frequency alone is below Nyquist. All storage is admitted before
    /// allocation, including retained pole declarations.
    pub fn new(
        sample_count: usize,
        timestep_s: f64,
        poles: &[DampedTonePole],
        budget: HarmonicFitBudget,
    ) -> Result<Self, HarmonicFitError> {
        if poles.is_empty() {
            return Err(invalid("at least one damped tone pole is required"));
        }
        for (k, pole) in poles.iter().enumerate() {
            if !pole.frequency_hz.is_finite()
                || pole.frequency_hz <= 0.0
                || !pole.decay_rate_per_s.is_finite()
                || pole.decay_rate_per_s < 0.0
            {
                return Err(invalid(
                    "damped poles need finite positive frequencies and nonnegative decay rates",
                ));
            }
            if poles[..k].contains(pole) {
                return Err(invalid("damped tone poles must be distinct"));
            }
        }
        // The shared admission already includes one f64 per retained tone;
        // damped declarations hold a frequency and a decay rate instead.
        let extra_workspace_bytes = poles
            .len()
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or_else(|| invalid("damped pole storage size overflows"))?;
        let design = FitDesign::with_envelopes(
            sample_count,
            timestep_s,
            poles.len(),
            |k| poles[k].frequency_hz * timestep_s,
            |j, k| {
                let decay = poles[k].decay_rate_per_s;
                if j == 0 || decay == 0.0 {
                    1.0
                } else {
                    // A positive exponent overflow means the envelope rounds
                    // to zero. Handle t=0 above to avoid an infinity * zero.
                    det::exp(-((decay * timestep_s) * j as f64))
                }
            },
            extra_workspace_bytes,
            budget,
        )?;
        Ok(Self {
            poles: poles.to_vec(),
            design,
        })
    }

    /// Supplied pole declarations in the same order as fitted phasors.
    #[must_use]
    pub fn poles(&self) -> &[DampedTonePole] {
        &self.poles
    }

    /// Required record and scratch length.
    #[must_use]
    pub const fn sample_count(&self) -> usize {
        self.design.samples
    }

    /// Uniform physical sampling interval in seconds.
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.design.timestep_s
    }

    /// Measured reciprocal 1-norm condition of the full triangular QR factor.
    #[must_use]
    pub const fn reciprocal_condition(&self) -> f64 {
        self.design.reciprocal_condition
    }

    /// Conservative numeric-array storage admitted at preparation.
    #[must_use]
    pub const fn workspace_bytes(&self) -> usize {
        self.design.workspace_bytes
    }

    /// Fit one record, allocating its scratch and returned first-sample phasors.
    pub fn fit(&self, samples: &[f64]) -> Result<HarmonicFit, HarmonicFitError> {
        self.design.fit(samples)
    }

    /// Fit without allocation, using length-N scratch and K phasor outputs.
    /// On error, phasors remain unchanged; scratch may have been modified.
    pub fn fit_into(
        &self,
        samples: &[f64],
        scratch: &mut [f64],
        phasors: &mut [C64],
    ) -> Result<HarmonicFitSummary, HarmonicFitError> {
        self.design.fit_into(samples, scratch, phasors)
    }
}
