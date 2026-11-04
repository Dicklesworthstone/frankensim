//! Local driven resonance and conventional half-power bandwidth of the original
//! complete system. A pole frequency alone does not locate an observed peak.
//!
//! The caller supplies an interval isolating the desired peak, a fixed physical
//! excitation and a frequency-independent real observation. Analytic response
//! derivatives bracket a change from increasing to decreasing amplitude. Every
//! trial factors the original descriptor pencil; no state, constraint, memory
//! branch or physical loss is removed, and no frequency grid is fitted.
//!
//! Brackets and residuals describe the numerical queries, not certified bounds
//! on a continuum or a guarantee that an interval contains only one peak.

use super::{HarmonicError, HarmonicModel, finite, poll};
use fs_math::{c64::C64, dd::Dd};

mod sensitivity;
pub use sensitivity::{DrivenResonanceParameterResponse, HalfPowerBandwidthDerivative};

/// A fixed real observation of original physical phasors, using the supplied
/// weights without normalization. Units and coordinate meaning belong to the
/// caller. Negative weights retain interference and antiresonances.
#[derive(Clone, Copy, Debug)]
pub enum HarmonicObservation<'a> {
    /// Weighted original differential states and algebraic multipliers.
    State(&'a [f64]),
    /// Weighted storage efforts and algebraic constraint reactions.
    Effort(&'a [f64]),
    /// Weighted original power-conjugate port outputs.
    Output(&'a [f64]),
}

/// Frequency resolution and total solve capacity for one local resonance query.
#[derive(Clone, Copy, Debug)]
pub struct ResonanceControls {
    /// Absolute numerical bracket width [rad/s], strictly positive.
    pub absolute_frequency_tolerance_rad_s: f64,
    /// Additional relative bracket width, in [0,1). The larger tolerance wins.
    pub relative_frequency_tolerance: f64,
    /// Maximum factor/response queries, including the peak and both crossings.
    /// Each derivative query additionally reuses its factor for one tangent solve.
    /// Must lie in 3..=4096; there is no unbounded search or implicit truncation.
    pub maximum_evaluations: usize,
    /// Also bracket the selected response's squared-amplitude half-height on
    /// both sides, inside the original caller interval. Missing crossings refuse.
    pub measure_half_power_bandwidth: bool,
}

/// A numerical frequency bracket resolved to the caller's requested width.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResonanceFrequency {
    /// Evaluated peak or midpoint estimate [rad/s].
    pub angular_frequency_rad_s: f64,
    /// Final sign/crossing bracket [rad/s], not an interval-arithmetic certificate.
    pub bracket_rad_s: (f64, f64),
}

/// Conventional -3 dB bandwidth of the SELECTED observation's amplitude squared.
/// This does not equate that observable with watts, stored energy, or material
/// heat. With overlapping peaks it need not equal a pole's decay quality factor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HalfPowerBandwidth {
    /// Crossing below the peak, at the selected squared-amplitude half-height.
    pub lower: ResonanceFrequency,
    /// Crossing above the peak, at the same selected half-height.
    pub upper: ResonanceFrequency,
    /// Upper minus lower crossing frequency [rad/s].
    pub bandwidth_rad_s: f64,
    /// Observed peak frequency divided by the measured bandwidth.
    pub response_quality_factor: f64,
}

/// One driven peak of the complete original model, with no state mutation.
#[derive(Clone, Debug, PartialEq)]
pub struct DrivenResonance {
    /// Local maximum resolved within the original caller interval.
    pub peak: ResonanceFrequency,
    /// Complex physical observation at the evaluated peak estimate.
    pub observation: C64,
    /// Derivative of that same observation with respect to angular frequency.
    pub observation_derivative: C64,
    /// Both resolved half-height crossings, only when requested.
    pub half_power: Option<HalfPowerBandwidth>,
    /// Actual number of factored physical response queries.
    pub evaluations: usize,
    /// Largest native response/tangent backward residual across all queries.
    /// Near resonance this is not a forward amplitude or frequency error bound.
    pub maximum_relative_residual: f64,
}

/// No partial peak/bandwidth is published on refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum DrivenResonanceError {
    /// Invalid controls, observation, excitation or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The caller's interval does not resolve the requested maximum/crossing.
    NotBracketed(&'static str),
    /// Caller-declared total evaluation capacity was exhausted.
    EvaluationLimit {
        /// Exhausted physical-query allowance.
        maximum: usize,
    },
    /// Floating-point frequency spacing cannot satisfy the requested tolerance.
    FrequencyResolution,
    /// Original system admission, solve, residual or cancellation refusal.
    Harmonic(HarmonicError),
}

impl core::fmt::Display for DrivenResonanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven resonance refused: {self:?}")
    }
}
impl std::error::Error for DrivenResonanceError {}
impl From<HarmonicError> for DrivenResonanceError {
    fn from(error: HarmonicError) -> Self {
        Self::Harmonic(error)
    }
}

#[derive(Clone, Copy)]
struct Sample {
    omega: f64,
    value: C64,
    derivative: C64,
    slope: f64,
}

impl HarmonicModel {
    /// Resolve a local maximum of a physical response under FIXED input phasors.
    /// The interval must have strictly positive amplitude slope at its lower
    /// endpoint and negative slope at its upper endpoint. If bandwidth is
    /// requested, both original endpoints must be below the selected peak's
    /// half-height. Caller knowledge must isolate the desired resonance: this
    /// search does not enumerate peaks or prove interval uniqueness.
    ///
    /// Frequencies are positive rad/s, observations use original units, and all
    /// native admission/cancellation/residual checks remain active. Exact poles,
    /// silent observations, unresolved stationary plateaus, missing brackets and
    /// exhausted evaluation capacity refuse. No linewidth or artificial damping
    /// is supplied. Existing harmonic preparation bounds every individual query.
    pub fn driven_resonance(
        &self,
        input: &[C64],
        observation: HarmonicObservation<'_>,
        interval_rad_s: (f64, f64),
        controls: ResonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenResonance, DrivenResonanceError> {
        poll(&mut cancelled)?;
        if input.len() != self.m || input.iter().any(|z| !finite(*z)) {
            return Err(DrivenResonanceError::Invalid(
                "one finite excitation per original port required",
            ));
        }
        let (weights, dimension) = match observation {
            HarmonicObservation::State(w) | HarmonicObservation::Effort(w) => (w, self.n),
            HarmonicObservation::Output(w) => (w, self.m),
        };
        if weights.len() != dimension
            || weights.iter().any(|w| !w.is_finite())
            || !weights.iter().any(|w| *w != 0.0)
        {
            return Err(DrivenResonanceError::Invalid(
                "finite nonzero observation on original coordinates required",
            ));
        }
        search_resonance(
            interval_rad_s,
            controls,
            cancelled,
            |omega, derivative, cancelled| {
                let factor = self.factor(omega, &mut *cancelled)?;
                let (response, slope, residual) = if derivative {
                    let result =
                        factor.response_with_frequency_derivative(input, None, &mut *cancelled)?;
                    let data = match observation {
                        HarmonicObservation::State(w) => (w, &result.state_derivative),
                        HarmonicObservation::Effort(w) => (w, &result.effort_derivative),
                        HarmonicObservation::Output(w) => (w, &result.output_derivative),
                    };
                    let slope = project(data.0, data.1.iter().copied(), cancelled)?;
                    (result.response, slope, result.relative_derivative_residual)
                } else {
                    (factor.response(input, &mut *cancelled)?, C64::ZERO, 0.0)
                };
                let data = match observation {
                    HarmonicObservation::State(w) => (w, &response.state),
                    HarmonicObservation::Effort(w) => (w, &response.effort),
                    HarmonicObservation::Output(w) => (w, &response.output),
                };
                let value = project(data.0, data.1.iter().copied(), cancelled)?;
                Ok((value, slope, residual.max(response.relative_residual)))
            },
        )
    }
}

// Both ordinary and wave-closed physical systems use one bracket/refinement
// owner. The evaluator returns the actual observed phasor, its optional exact
// tangent, and the largest independently checked physical equation residual.
// `derivative=false` avoids the second solve when refining half-height crossings.
pub(super) fn search_resonance<C, E>(
    interval_rad_s: (f64, f64),
    controls: ResonanceControls,
    mut cancelled: C,
    evaluate: E,
) -> Result<DrivenResonance, DrivenResonanceError>
where
    C: FnMut() -> bool,
    E: FnMut(f64, bool, &mut C) -> Result<(C64, C64, f64), DrivenResonanceError>,
{
    poll(&mut cancelled)?;
    let (low, high) = interval_rad_s;
    if !low.is_finite()
        || !high.is_finite()
        || low <= 0.0
        || high <= low
        || !controls.absolute_frequency_tolerance_rad_s.is_finite()
        || controls.absolute_frequency_tolerance_rad_s <= 0.0
        || !controls.relative_frequency_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.relative_frequency_tolerance)
        || !(3..=4096).contains(&controls.maximum_evaluations)
    {
        return Err(DrivenResonanceError::Invalid(
            "positive ordered interval and bounded resolution controls required",
        ));
    }
    let mut search = Search {
        evaluate,
        controls,
        cancelled,
        evaluations: 0,
        residual: 0.0,
    };
    let outer_low = search.sample(low, true)?;
    let outer_high = search.sample(high, true)?;
    let mut left = outer_low;
    let mut right = outer_high;
    if left.slope <= 0.0 || right.slope >= 0.0 {
        return Err(DrivenResonanceError::NotBracketed(
            "a peak requires increasing then decreasing nonzero response",
        ));
    }
    let mut stationary = None;
    while !search.resolved(left.omega, right.omega) {
        if let Some(middle) = stationary {
            // Reuse the exactly stationary evaluated response while its
            // two sign endpoints contract; do not refactor its frequency.
            let a = search.sample(interior_midpoint(left.omega, middle)?, true)?;
            let b = search.sample(interior_midpoint(middle, right.omega)?, true)?;
            if a.slope <= 0.0 || b.slope >= 0.0 {
                return Err(DrivenResonanceError::NotBracketed(
                    "stationary response does not retain a resolved maximum bracket",
                ));
            }
            left = a;
            right = b;
            continue;
        }
        let middle = search.sample(interior_midpoint(left.omega, right.omega)?, true)?;
        if middle.slope > 0.0 {
            left = middle;
        } else if middle.slope < 0.0 {
            right = middle;
        } else {
            // An exactly zero rounded slope is not a zero-width frequency
            // certificate. Retain an actual smaller sign bracket around it.
            stationary = Some(middle.omega);
        }
    }
    let peak_sample = search.sample(f64::midpoint(left.omega, right.omega), true)?;
    let peak = ResonanceFrequency {
        angular_frequency_rad_s: peak_sample.omega,
        bracket_rad_s: (left.omega, right.omega),
    };
    let half_power = if controls.measure_half_power_bandwidth {
        let target = peak_sample.value.abs() * core::f64::consts::FRAC_1_SQRT_2;
        if !target.is_finite() || target <= 0.0 {
            return Err(DrivenResonanceError::Invalid(
                "half-height amplitude is not representable",
            ));
        }
        if outer_low.value.abs() >= target || outer_high.value.abs() >= target {
            return Err(DrivenResonanceError::NotBracketed(
                "both half-height crossings must lie inside the supplied interval",
            ));
        }
        let lower = search.crossing(outer_low, peak_sample, target)?;
        let upper = search.crossing(outer_high, peak_sample, target)?;
        let bandwidth_rad_s = upper.angular_frequency_rad_s - lower.angular_frequency_rad_s;
        let response_quality_factor = peak_sample.omega / bandwidth_rad_s;
        if bandwidth_rad_s <= 0.0
            || !response_quality_factor.is_finite()
            || lower.bracket_rad_s.1 >= upper.bracket_rad_s.0
        {
            return Err(DrivenResonanceError::FrequencyResolution);
        }
        Some(HalfPowerBandwidth {
            lower,
            upper,
            bandwidth_rad_s,
            response_quality_factor,
        })
    } else {
        None
    };
    poll(&mut search.cancelled)?;
    Ok(DrivenResonance {
        peak,
        observation: peak_sample.value,
        observation_derivative: peak_sample.derivative,
        half_power,
        evaluations: search.evaluations,
        maximum_relative_residual: search.residual,
    })
}

struct Search<C, E> {
    evaluate: E,
    controls: ResonanceControls,
    cancelled: C,
    evaluations: usize,
    residual: f64,
}

impl<C, E> Search<C, E>
where
    C: FnMut() -> bool,
    E: FnMut(f64, bool, &mut C) -> Result<(C64, C64, f64), DrivenResonanceError>,
{
    fn resolved(&self, low: f64, high: f64) -> bool {
        high - low
            <= self
                .controls
                .absolute_frequency_tolerance_rad_s
                .max(self.controls.relative_frequency_tolerance * f64::midpoint(low, high))
    }

    fn sample(&mut self, omega: f64, derivative: bool) -> Result<Sample, DrivenResonanceError> {
        poll(&mut self.cancelled)?;
        if self.evaluations == self.controls.maximum_evaluations {
            return Err(DrivenResonanceError::EvaluationLimit {
                maximum: self.evaluations,
            });
        }
        self.evaluations += 1;
        let (value, slope, residual) = (self.evaluate)(omega, derivative, &mut self.cancelled)?;
        self.residual = self.residual.max(residual);
        if derivative && value == C64::ZERO {
            return Err(DrivenResonanceError::NotBracketed(
                "a silent observation or exact antiresonance is not a resonance peak",
            ));
        }
        let logarithmic_slope = if derivative { (slope / value).re } else { 0.0 };
        if !logarithmic_slope.is_finite() || !value.abs().is_finite() {
            return Err(DrivenResonanceError::Invalid(
                "observed amplitude or logarithmic slope is not representable",
            ));
        }
        Ok(Sample {
            omega,
            value,
            derivative: slope,
            slope: logarithmic_slope,
        })
    }

    fn crossing(
        &mut self,
        mut outside: Sample,
        mut inside: Sample,
        target: f64,
    ) -> Result<ResonanceFrequency, DrivenResonanceError> {
        loop {
            let low = outside.omega.min(inside.omega);
            let high = outside.omega.max(inside.omega);
            if self.resolved(low, high) {
                return Ok(ResonanceFrequency {
                    angular_frequency_rad_s: f64::midpoint(low, high),
                    bracket_rad_s: (low, high),
                });
            }
            let sample = self.sample(interior_midpoint(low, high)?, false)?;
            if sample.value.abs() < target {
                outside = sample;
            } else {
                inside = sample;
            }
        }
    }
}

fn interior_midpoint(low: f64, high: f64) -> Result<f64, DrivenResonanceError> {
    let middle = f64::midpoint(low, high);
    if middle <= low || middle >= high {
        Err(DrivenResonanceError::FrequencyResolution)
    } else {
        Ok(middle)
    }
}

pub(super) fn project(
    weights: &[f64],
    values: impl IntoIterator<Item = C64>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<C64, DrivenResonanceError> {
    let mut real = Dd::ZERO;
    let mut imaginary = Dd::ZERO;
    for (i, (&weight, value)) in weights.iter().zip(values).enumerate() {
        if i % 64 == 0 {
            poll(cancelled)?;
        }
        real = real + Dd::from_f64(weight) * Dd::from_f64(value.re);
        imaginary = imaginary + Dd::from_f64(weight) * Dd::from_f64(value.im);
    }
    let value = C64::new(real.to_f64(), imaginary.to_f64());
    if finite(value) {
        Ok(value)
    } else {
        Err(DrivenResonanceError::Invalid(
            "physical observation is not representable",
        ))
    }
}
