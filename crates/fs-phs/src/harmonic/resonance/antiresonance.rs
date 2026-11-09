//! Observed destructive-interference minima of the complete original response.
//! The same bounded refinement owner handles peaks, notches, and crossings.

use super::{
    DrivenResonanceError, HarmonicModel, HarmonicObservation, ResonanceControls,
    ResonanceFrequency, ResonanceSearchError, Search, StationaryKind, finite, poll, project,
    search_poll,
};
use fs_math::{c64::C64, det};

/// Frequency resolution and total physical-query capacity for one local notch.
#[derive(Clone, Copy, Debug)]
pub struct AntiresonanceControls {
    /// Absolute numerical bracket width [rad/s], strictly positive.
    pub absolute_frequency_tolerance_rad_s: f64,
    /// Additional relative bracket width in [0,1). The larger tolerance wins.
    pub relative_frequency_tolerance: f64,
    /// Total response queries, including tangents and requested recovery edges.
    /// Must lie in 3..=4096. Native tangents reuse their response factorization.
    pub maximum_evaluations: usize,
    /// Optional positive, finite, ABSOLUTE level of the selected |observation|^2
    /// at which to locate both recovery edges. Its units are the squared units
    /// of the caller's coherent observation; it is not necessarily watts.
    /// Both original endpoints must exceed this level and the minimum must be
    /// strictly below it. No baseline, relative notch depth, or Q is inferred.
    pub recovery_power: Option<f64>,
}

/// Two recovery crossings at a caller-declared absolute squared-amplitude level.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AntiresonanceRecovery {
    /// Caller-declared |observation|^2 recovery level, in original squared units.
    pub power_level: f64,
    /// Recovery crossing below the observed minimum.
    pub lower: ResonanceFrequency,
    /// Recovery crossing above the observed minimum.
    pub upper: ResonanceFrequency,
    /// Upper minus lower crossing frequency [rad/s]. This is a notch width,
    /// not a material loss, a free-pole linewidth, or a quality factor.
    pub bandwidth_rad_s: f64,
}

/// One local minimum of a coherent original physical response.
#[derive(Clone, Debug, PartialEq)]
pub struct DrivenAntiresonance {
    /// Minimum estimate and its strict decreasing-to-increasing sign bracket.
    pub minimum: ResonanceFrequency,
    /// Actual complex observation at the evaluated frequency. A finite bracket
    /// alone never rounds this to zero or certifies an exact cancellation.
    pub observation: C64,
    /// Analytic total derivative of that same observation per rad/s. This may
    /// be nonzero at an exact response zero; no log phase is formed there.
    pub observation_derivative: C64,
    /// Both recovery edges at the requested absolute level, when supplied.
    pub recovery: Option<AntiresonanceRecovery>,
    /// Actual physical response queries, including optional recovery edges.
    pub evaluations: usize,
    /// Largest native response/tangent backward residual across all queries.
    /// At cancellation this is not a relative forward error in the notch depth.
    pub maximum_relative_residual: f64,
}

/// Resolve one observed antiresonance without dividing by its response.
///
/// `evaluate(omega, need_derivative, cancelled)` has the same source contract as
/// [`super::search_resonance`]: the complete coherent complex observation, its
/// analytic TOTAL derivative per rad/s, and the largest independently checked
/// backward residual. Source, material, load, and observation frequency changes
/// belong in that derivative. Every query describes the same physical experiment.
/// The source keeps its original coordinates, constraints, memory, loss, physical
/// loads, budgets, cancellation, and residual admission. No operator is rebuilt
/// here, and no response or derivative is fitted from a frequency grid.
///
/// The caller must isolate one minimum in a positive interval whose endpoints
/// have strictly decreasing then increasing squared amplitude. The search uses
/// the sign of `Re(conj(H) H')`, scaled and accumulated in double-double precision
/// without squaring physical amplitudes or dividing by `H`. Exact zeros and
/// damped nonzero minima are both admitted. Silent responses, unresolved flat
/// minima, and missing brackets refuse. A bracket is numerical evidence, not a
/// continuum certificate, uniqueness proof, or bound on the minimum's depth.
///
/// With `recovery_power`, both crossings at that caller-declared absolute level
/// must lie inside the original interval. `need_derivative=false` permits the
/// source to skip its tangent for those crossings; that derivative slot is then
/// ignored. All work shares one total evaluation cap, polls cancellation before
/// and after physical queries, and returns no partial result on refusal.
pub fn search_antiresonance<E, C, F>(
    interval_rad_s: (f64, f64),
    controls: AntiresonanceControls,
    cancelled: C,
    evaluate: F,
) -> Result<DrivenAntiresonance, ResonanceSearchError<E>>
where
    C: FnMut() -> bool,
    F: FnMut(f64, bool, &mut C) -> Result<(C64, C64, f64), E>,
{
    let mut search = Search::new(
        interval_rad_s,
        ResonanceControls {
            absolute_frequency_tolerance_rad_s: controls.absolute_frequency_tolerance_rad_s,
            relative_frequency_tolerance: controls.relative_frequency_tolerance,
            maximum_evaluations: controls.maximum_evaluations,
            measure_half_power_bandwidth: false,
        },
        cancelled,
        evaluate,
    )?;
    if controls
        .recovery_power
        .is_some_and(|level| !level.is_finite() || level <= 0.0)
    {
        return Err(ResonanceSearchError::Invalid(
            "a positive finite absolute recovery power is required",
        ));
    }
    let outer_low = search.sample_kind(interval_rad_s.0, true, StationaryKind::Minimum)?;
    let outer_high = search.sample_kind(interval_rad_s.1, true, StationaryKind::Minimum)?;
    let (sample, bracket_rad_s) =
        search.stationary(outer_low, outer_high, StationaryKind::Minimum)?;
    let recovery = if let Some(power_level) = controls.recovery_power {
        // Compare amplitudes to sqrt(level), avoiding overflow/underflow of
        // |H|^2. The level remains in the caller's physical squared units.
        let target = det::sqrt(power_level);
        if sample.value.abs() >= target
            || outer_low.value.abs() <= target
            || outer_high.value.abs() <= target
        {
            return Err(ResonanceSearchError::NotBracketed(
                "both recovery crossings above the minimum must lie inside the supplied interval",
            ));
        }
        let lower = search.crossing_kind(outer_low, sample, target, StationaryKind::Minimum)?;
        let upper = search.crossing_kind(outer_high, sample, target, StationaryKind::Minimum)?;
        let bandwidth_rad_s = upper.angular_frequency_rad_s - lower.angular_frequency_rad_s;
        if bandwidth_rad_s <= 0.0 || lower.bracket_rad_s.1 >= upper.bracket_rad_s.0 {
            return Err(ResonanceSearchError::FrequencyResolution);
        }
        Some(AntiresonanceRecovery {
            power_level,
            lower,
            upper,
            bandwidth_rad_s,
        })
    } else {
        None
    };
    search_poll(&mut search.cancelled)?;
    Ok(DrivenAntiresonance {
        minimum: ResonanceFrequency {
            angular_frequency_rad_s: sample.omega,
            bracket_rad_s,
        },
        observation: sample.value,
        observation_derivative: sample.derivative,
        recovery,
        evaluations: search.evaluations,
        maximum_relative_residual: search.residual,
    })
}

impl HarmonicModel {
    /// Resolve an observed cancellation/minimum under fixed coherent port inputs
    /// in the complete original descriptor model. The interval must isolate one
    /// notch and bracket strictly decreasing then increasing squared amplitude.
    ///
    /// Original physical observation weights are used without normalization;
    /// signed weights preserve interference. Each trial reuses one original
    /// factor for response and exact tangent, retaining every constraint and
    /// material-memory state. A zero observed response does not imply zero
    /// internal state, energy, or loss. Optional recovery levels use the selected
    /// observation's squared units and imply no modal or material quality factor.
    pub fn driven_antiresonance(
        &self,
        input: &[C64],
        observation: HarmonicObservation<'_>,
        interval_rad_s: (f64, f64),
        controls: AntiresonanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenAntiresonance, DrivenResonanceError> {
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
        search_antiresonance(
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
        .map_err(DrivenResonanceError::from)
    }
}
