//! Resolve an isolated finite resonance from the original complex response.
//!
//! The caller supplies a positive-frequency interval containing ONE magnitude
//! maximum, with monotonically rising/falling flanks and both half-power
//! crossings inside it. No initial frequency grid or modal fit is needed.
//! A closure can select a scalar duct impedance, one prepared modal response,
//! or an explicitly driven/observed network port, retaining its original
//! geometry, loss law, radiation, gas state and boundary conditions.
//!
//! This measures the peak of |response| and its -3.0103 dB amplitude crossings.
//! `Q_bandwidth = omega_peak / (omega_upper - omega_lower)` is an observable
//! response width, NOT a complex-pole quality factor, decay time, or natural
//! frequency. Background response and unresolved overlapping modes can change
//! it. Sampled violations of the one-peak assumption refuse; finite samples
//! cannot certify absence of an arbitrary hidden feature between queries.
//!
//! Golden-section maximum refinement and bracketed crossing refinement reuse
//! every frequency evaluation. Frequency and linewidth tolerances are separate:
//! the peak bracket must also be small relative to the measured linewidth.
//! Thus an unresolved lossless pole cannot acquire a fabricated finite Q from
//! the largest sampled amplitude. No damping or gain clipping is introduced.

use crate::DuctError;
use fs_math::{c64::C64, det};
use std::collections::BTreeMap;

/// Maximum retained original response evaluations for one refinement.
pub const MAX_RESONANCE_EVALUATIONS: usize = 65_536;

/// Numerical resolution of a caller-isolated response peak.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResonanceOptions {
    /// Final peak/crossing bracket widths relative to their angular frequency.
    pub relative_frequency_tolerance: f64,
    /// Crossing widths relative to their distance from the peak; the peak
    /// bracket is additionally at most this fraction of the full bandwidth.
    pub relative_bandwidth_tolerance: f64,
    /// Hard count of original callback evaluations, including endpoint probes.
    pub max_evaluations: usize,
}

impl Default for ResonanceOptions {
    fn default() -> Self {
        Self {
            relative_frequency_tolerance: 1.0e-10,
            relative_bandwidth_tolerance: 1.0e-6,
            max_evaluations: 256,
        }
    }
}

/// Measured finite response resonance and its actual half-power bandwidth.
#[derive(Debug, Clone, PartialEq)]
pub struct ResponseResonance {
    /// Angular frequency of the largest retained response in the final bracket.
    pub angular_frequency_rad_s: f64,
    /// Original complex response at that same frequency, including phase.
    pub response: C64,
    /// Final search interval around the isolated magnitude maximum.
    pub peak_bracket_rad_s: [f64; 2],
    /// Lower and upper frequencies where |response| = |peak| / sqrt(2).
    pub half_power_frequencies_rad_s: [f64; 2],
    /// Final sign-changing intervals of the two half-power searches.
    pub half_power_brackets_rad_s: [[f64; 2]; 2],
    /// Difference between the upper and lower half-power frequencies [rad/s].
    pub bandwidth_rad_s: f64,
    /// omega_peak / bandwidth. This is a response measure, not a pole/decay Q.
    pub bandwidth_quality_factor: f64,
    /// Number of original response evaluations; cached queries do not count.
    pub evaluations: usize,
}

/// A refusal never publishes an incomplete peak or a guessed linewidth.
#[derive(Debug, Clone, PartialEq)]
pub enum ResonanceError {
    /// Bounds, tolerances or work limits cannot admit a refinement.
    InvalidInput {
        /// Named invalid input.
        what: &'static str,
    },
    /// An original physical model refusal, unchanged and at its query frequency.
    Response {
        /// Original callback query [rad/s].
        omega: f64,
        /// Original physical-model refusal.
        source: DuctError,
    },
    /// The callback returned a nonfinite real or imaginary component.
    NonFiniteResponse {
        /// Original callback query [rad/s].
        omega: f64,
    },
    /// The interval does not enclose a positive finite interior magnitude peak.
    PeakNotBracketed,
    /// An observed reversal contradicts the caller's single-peak interval.
    NotIsolated,
    /// The named outer endpoint is still above the requested half-power level.
    HalfPowerNotBracketed {
        /// True for the lower-frequency flank, false for the upper flank.
        lower: bool,
    },
    /// Requested intervals cannot be resolved in represented arithmetic.
    ResolutionLimit,
    /// The callback budget was exhausted before all brackets were resolved.
    EvaluationLimit {
        /// Completed original response evaluations.
        evaluations: usize,
    },
    /// The caller cancelled before a complete result could be published.
    Cancelled {
        /// Completed original response evaluations.
        evaluations: usize,
    },
}

impl core::fmt::Display for ResonanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-DUCT-RESONANCE: {self:?}")
    }
}
impl std::error::Error for ResonanceError {}

#[derive(Clone, Copy)]
struct Sample {
    omega: f64,
    response: C64,
    log_amplitude: f64,
}

struct Samples<F, C> {
    evaluate: F,
    cancelled: C,
    limit: usize,
    values: BTreeMap<u64, Sample>,
}

impl<F, C> Samples<F, C>
where
    F: FnMut(f64) -> Result<C64, DuctError>,
    C: FnMut() -> bool,
{
    fn poll(&mut self) -> Result<(), ResonanceError> {
        if (self.cancelled)() {
            Err(ResonanceError::Cancelled {
                evaluations: self.values.len(),
            })
        } else {
            Ok(())
        }
    }

    fn at(&mut self, omega: f64) -> Result<Sample, ResonanceError> {
        self.poll()?;
        if let Some(value) = self.values.get(&omega.to_bits()) {
            return Ok(*value);
        }
        if self.values.len() == self.limit {
            return Err(ResonanceError::EvaluationLimit {
                evaluations: self.values.len(),
            });
        }
        let response =
            (self.evaluate)(omega).map_err(|source| ResonanceError::Response { omega, source })?;
        if !response.re.is_finite() || !response.im.is_finite() {
            return Err(ResonanceError::NonFiniteResponse { omega });
        }
        // Compare amplitudes without squaring, overflowing, or underflowing
        // the physical response. Constant unit/gain changes preserve the search.
        let scale = response.re.abs().max(response.im.abs());
        let log_amplitude = if scale == 0.0 {
            f64::NEG_INFINITY
        } else {
            let re = response.re / scale;
            let im = response.im / scale;
            det::ln(scale) + 0.5 * det::ln(re * re + im * im)
        };
        let value = Sample {
            omega,
            response,
            log_amplitude,
        };
        self.values.insert(omega.to_bits(), value);
        self.poll()?;
        Ok(value)
    }

    fn isolated(&self, peak: Sample) -> Result<(), ResonanceError> {
        let mut before: Option<Sample> = None;
        for &current in self.values.values() {
            if let Some(previous) = before {
                let noise = 128.0 * f64::EPSILON * peak.log_amplitude.abs().max(1.0);
                let wrong = if current.omega <= peak.omega {
                    previous.log_amplitude > current.log_amplitude + noise
                } else if previous.omega >= peak.omega {
                    current.log_amplitude > previous.log_amplitude + noise
                } else {
                    false
                };
                if wrong {
                    return Err(ResonanceError::NotIsolated);
                }
            }
            before = Some(current);
        }
        Ok(())
    }

    #[allow(clippy::float_cmp)] // Exact endpoint identity detects exhausted f64 resolution.
    fn crossing(
        &mut self,
        peak: Sample,
        lower: bool,
        options: &ResonanceOptions,
    ) -> Result<[f64; 2], ResonanceError> {
        let target = peak.log_amplitude - 0.5 * core::f64::consts::LN_2;
        let mut a = *self
            .values
            .first_key_value()
            .expect("endpoints evaluated")
            .1;
        let mut b = *self.values.last_key_value().expect("endpoints evaluated").1;
        if (if lower { a } else { b }).log_amplitude >= target {
            return Err(ResonanceError::HalfPowerNotBracketed { lower });
        }
        // Tightest already observed crossing bracket on the selected flank.
        // The one-peak check precedes this search and runs again after it.
        if lower {
            b = peak;
            for &sample in self.values.range(..=peak.omega.to_bits()).map(|(_, v)| v) {
                if sample.log_amplitude < target {
                    a = sample;
                } else {
                    b = sample;
                    break;
                }
            }
        } else {
            a = peak;
            for &sample in self.values.range(peak.omega.to_bits()..).map(|(_, v)| v) {
                if sample.log_amplitude > target {
                    a = sample;
                } else {
                    b = sample;
                    break;
                }
            }
        }
        loop {
            self.poll()?;
            let middle = f64::midpoint(a.omega, b.omega);
            let tolerance = (options.relative_frequency_tolerance * middle)
                .min(options.relative_bandwidth_tolerance * (middle - peak.omega).abs());
            if b.omega - a.omega <= tolerance {
                return Ok([a.omega, b.omega]);
            }
            if middle == a.omega || middle == b.omega {
                return Err(ResonanceError::ResolutionLimit);
            }
            let sample = self.at(middle)?;
            if (sample.log_amplitude < target) == lower {
                a = sample;
            } else {
                b = sample;
            }
        }
    }
}

struct Maximum {
    bracket: [f64; 2],
    left: Sample,
    right: Sample,
}

impl Maximum {
    fn refine<F, C>(
        &mut self,
        samples: &mut Samples<F, C>,
        relative_frequency_tolerance: f64,
        linewidth_tolerance: f64,
    ) -> Result<Sample, ResonanceError>
    where
        F: FnMut(f64) -> Result<C64, DuctError>,
        C: FnMut() -> bool,
    {
        // sqrt(5) is deterministic, shared with the original real math owner.
        let ratio = 0.5 * (det::sqrt(5.0) - 1.0);
        while self.bracket[1] - self.bracket[0]
            > (relative_frequency_tolerance * self.left.omega).min(linewidth_tolerance)
        {
            samples.poll()?;
            if self.left.log_amplitude > self.right.log_amplitude {
                self.bracket[1] = self.right.omega;
                self.right = self.left;
                let omega = self.bracket[1] - ratio * (self.bracket[1] - self.bracket[0]);
                if omega <= self.bracket[0] || omega >= self.right.omega {
                    return Err(ResonanceError::ResolutionLimit);
                }
                self.left = samples.at(omega)?;
            } else if self.left.log_amplitude < self.right.log_amplitude {
                self.bracket[0] = self.left.omega;
                self.left = self.right;
                let omega = self.bracket[0] + ratio * (self.bracket[1] - self.bracket[0]);
                if omega <= self.left.omega || omega >= self.bracket[1] {
                    return Err(ResonanceError::ResolutionLimit);
                }
                self.right = samples.at(omega)?;
            } else {
                // Equal samples retain the central interval of a unimodal peak.
                if self.left.omega <= self.bracket[0] && self.right.omega >= self.bracket[1] {
                    return Err(ResonanceError::ResolutionLimit);
                }
                self.bracket = [self.left.omega, self.right.omega];
                let width = self.bracket[1] - self.bracket[0];
                self.left = samples.at(self.bracket[1] - ratio * width)?;
                self.right = samples.at(self.bracket[0] + ratio * width)?;
            }
        }
        let middle = samples.at(f64::midpoint(self.bracket[0], self.bracket[1]))?;
        Ok([self.left, self.right, middle]
            .into_iter()
            .max_by(|a, b| a.log_amplitude.total_cmp(&b.log_amplitude))
            .unwrap())
    }
}

/// Refine one caller-isolated finite magnitude peak and both half-power flanks.
///
/// Frequencies are rad/s under the original callback's phasor convention. The
/// callback must describe the same model/drive/observer at every query and
/// return the SAME response on repeated queries. The supplied interval must
/// contain one unimodal peak and both half-power crossings. Wider multimode
/// sweeps must first isolate individual response peaks; this is no spectrum
/// completeness claim. No interpolation or continuation of a table is inferred.
///
/// Cancellation is checked before/after original evaluations and between
/// refinement steps. No callback is repeated at an identical frequency.
///
/// # Errors
/// Invalid bounds/tolerances, original solver refusal, nonfinite response,
/// missing/observably nonisolated peak, absent half-power flank, exhausted work,
/// cancellation, or insufficient floating-point resolution.
pub fn refine_response_resonance(
    interval_rad_s: [f64; 2],
    options: ResonanceOptions,
    evaluate: impl FnMut(f64) -> Result<C64, DuctError>,
    cancelled: impl FnMut() -> bool,
) -> Result<ResponseResonance, ResonanceError> {
    let [lo, hi] = interval_rad_s;
    if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || hi <= lo {
        return Err(ResonanceError::InvalidInput {
            what: "positive finite increasing frequency interval required",
        });
    }
    if [
        options.relative_frequency_tolerance,
        options.relative_bandwidth_tolerance,
    ]
    .iter()
    .any(|&v| !v.is_finite() || !(8.0 * f64::EPSILON..1.0).contains(&v))
        || options.max_evaluations == 0
        || options.max_evaluations > MAX_RESONANCE_EVALUATIONS
    {
        return Err(ResonanceError::InvalidInput {
            what: "representable relative tolerances and bounded positive evaluation count required",
        });
    }
    let mut samples = Samples {
        evaluate,
        cancelled,
        limit: options.max_evaluations,
        values: BTreeMap::new(),
    };
    let (low, high) = (samples.at(lo)?, samples.at(hi)?);
    let ratio = 0.5 * (det::sqrt(5.0) - 1.0);
    let mut maximum = Maximum {
        bracket: interval_rad_s,
        left: samples.at(hi - ratio * (hi - lo))?,
        right: samples.at(lo + ratio * (hi - lo))?,
    };
    let mut tolerance = f64::INFINITY;
    loop {
        let peak = maximum.refine(
            &mut samples,
            options.relative_frequency_tolerance,
            tolerance,
        )?;
        if !peak.log_amplitude.is_finite()
            || peak.log_amplitude <= low.log_amplitude
            || peak.log_amplitude <= high.log_amplitude
        {
            return Err(ResonanceError::PeakNotBracketed);
        }
        samples.isolated(peak)?;
        let lower = samples.crossing(peak, true, &options)?;
        let upper = samples.crossing(peak, false, &options)?;
        samples.isolated(peak)?;
        let half_power = [
            f64::midpoint(lower[0], lower[1]),
            f64::midpoint(upper[0], upper[1]),
        ];
        let bandwidth = half_power[1] - half_power[0];
        let q = peak.omega / bandwidth;
        if !bandwidth.is_finite() || bandwidth <= 0.0 || !q.is_finite() {
            return Err(ResonanceError::ResolutionLimit);
        }
        tolerance = (options.relative_frequency_tolerance * peak.omega)
            .min(options.relative_bandwidth_tolerance * bandwidth);
        if maximum.bracket[1] - maximum.bracket[0] <= tolerance {
            samples.poll()?;
            return Ok(ResponseResonance {
                angular_frequency_rad_s: peak.omega,
                response: peak.response,
                peak_bracket_rad_s: maximum.bracket,
                half_power_frequencies_rad_s: half_power,
                half_power_brackets_rad_s: [lower, upper],
                bandwidth_rad_s: bandwidth,
                bandwidth_quality_factor: q,
                evaluations: samples.values.len(),
            });
        }
    }
}

#[cfg(test)]
mod tests;
