//! Discover observable driven peaks resolved by an explicit frequency grid.

use super::{
    DrivenResonance, ResonanceControls, ResonanceSearchError, finite, search_poll, search_resonance,
};
use fs_math::c64::C64;

mod physical;

/// Resolution and total source-query capacity for a driven-response survey.
#[derive(Clone, Copy, Debug)]
pub struct ResonanceSurveyControls {
    /// Strictly positive absolute peak-bracket width [rad/s].
    pub absolute_frequency_tolerance_rad_s: f64,
    /// Relative peak-bracket width in [0,1); the larger tolerance wins.
    pub relative_frequency_tolerance: f64,
    /// Total physical source calls, including the entire grid and all peak
    /// refinements. Must lie in 1..=1_048_576. Cache hits consume no calls.
    pub maximum_evaluations: usize,
    /// Maximum returned peaks, in 1..=4096. A further sign-bracketed candidate
    /// refuses the whole survey rather than silently truncating the list.
    pub maximum_peaks: usize,
}

/// Peaks observed under one physical excitation and observation over a grid.
#[derive(Clone, Debug, PartialEq)]
pub struct ResonanceSurvey {
    /// Ascending frequency; every `half_power` is None. Individual evaluation
    /// counts include new refinement calls only; shared grid calls belong to
    /// the survey total. Each peak still reports its own worst queried residual.
    pub peaks: Vec<DrivenResonance>,
    /// First and last supplied grid frequency [rad/s].
    pub interval_rad_s: (f64, f64),
    /// Number of supplied frequencies, all evaluated on a successful survey.
    pub grid_points: usize,
    /// Actual physical source calls, excluding reused grid samples.
    pub evaluations: usize,
    /// Worst original response/tangent backward residual, including grid points.
    /// This does not bound peak-frequency or forward-amplitude error.
    pub maximum_relative_residual: f64,
}

/// No partial list is published on source, refinement, cancellation or capacity
/// refusal. The native error remains inside `ResonanceSearchError::Source`.
#[derive(Clone, Debug, PartialEq)]
pub enum ResonanceSurveyError<E> {
    Search(ResonanceSearchError<E>),
    PeakLimit { maximum: usize },
}

impl<E> From<ResonanceSearchError<E>> for ResonanceSurveyError<E> {
    fn from(error: ResonanceSearchError<E>) -> Self {
        Self::Search(error)
    }
}

impl<E: core::fmt::Debug> core::fmt::Display for ResonanceSurveyError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven resonance survey refused: {self:?}")
    }
}

impl<E: std::error::Error + 'static> std::error::Error for ResonanceSurveyError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Search(error) => Some(error),
            Self::PeakLimit { .. } => None,
        }
    }
}

#[derive(Clone, Copy)]
struct GridSample {
    omega: f64,
    value: C64,
    derivative: C64,
    residual: f64,
    slope: f64,
}

/// Discover and refine the driven peaks resolved by a supplied frequency grid.
///
/// The grid must contain 2..=4096 strictly increasing finite positive rad/s
/// frequencies. Each point queries the same physical response and its complete
/// analytic frequency derivative, with the callback contract of
/// [`search_resonance`]. Every callback requests `need_derivative=true`.
/// No interpolation, fitted response or artificial loss is introduced.
///
/// Positive-to-negative amplitude-slope transitions identify candidates. A
/// stationary grid point between those signs belongs to one candidate, not two.
/// An exact zero separates candidates: silent or unobservable modes and sampled
/// antiresonances are not peaks. Band endpoints are not two-sided maxima.
/// Refinement uses the existing single-peak owner and reuses grid samples,
/// including candidate endpoints and a peak lying exactly on the grid.
///
/// This is explicitly a GRID-RESOLVED survey, not a free-mode census or proof of
/// spectral completeness. Unresolved peak/zero pairs can lie between samples;
/// callers must supply a grid that separates the features they need. A finite
/// grid does not prove uniqueness within a candidate interval. An empty result
/// means no maximum was resolved by these samples, not that the system has no
/// resonances. Flat candidates or native singularities still refuse.
///
/// Results contain peak frequencies, complex observations and analytic tangents.
/// Half-power crossings can span neighbouring peaks; use an independently
/// isolated [`search_resonance`] interval to request that observable's linewidth.
/// One global evaluation budget covers the grid and every refinement; the
/// single-peak owner's existing 4096-query bound also remains active per peak.
pub fn survey_resonances<E, C, F>(
    frequencies_rad_s: &[f64],
    controls: ResonanceSurveyControls,
    mut cancelled: C,
    evaluate: F,
) -> Result<ResonanceSurvey, ResonanceSurveyError<E>>
where
    C: FnMut() -> bool,
    F: FnMut(f64, bool, &mut C) -> Result<(C64, C64, f64), E>,
{
    search_poll(&mut cancelled)?;
    if !(2..=4096).contains(&frequencies_rad_s.len())
        || frequencies_rad_s
            .iter()
            .any(|omega| !omega.is_finite() || *omega <= 0.0)
        || frequencies_rad_s.windows(2).any(|pair| pair[0] >= pair[1])
        || !controls.absolute_frequency_tolerance_rad_s.is_finite()
        || controls.absolute_frequency_tolerance_rad_s <= 0.0
        || !controls.relative_frequency_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.relative_frequency_tolerance)
        || !(1..=1_048_576).contains(&controls.maximum_evaluations)
        || !(1..=4096).contains(&controls.maximum_peaks)
    {
        return Err(ResonanceSearchError::Invalid(
            "ordered positive grid, finite resolution and bounded survey capacity required",
        )
        .into());
    }
    let mut sampling = Sampling {
        evaluate,
        maximum_evaluations: controls.maximum_evaluations,
        evaluations: 0,
        maximum_relative_residual: 0.0,
    };
    let mut grid = Vec::with_capacity(frequencies_rad_s.len());
    for &omega in frequencies_rad_s {
        grid.push(sampling.query(omega, &mut cancelled)?);
    }
    let mut brackets = Vec::new();
    let mut rising = None;
    for (index, sample) in grid.iter().enumerate() {
        search_poll(&mut cancelled)?;
        if sample.value == C64::ZERO {
            rising = None;
        } else if sample.slope > 0.0 {
            rising = Some(index);
        } else if sample.slope < 0.0 {
            if let Some(left) = rising.take() {
                if brackets.len() == controls.maximum_peaks {
                    return Err(ResonanceSurveyError::PeakLimit {
                        maximum: controls.maximum_peaks,
                    });
                }
                brackets.push((left, index));
            }
        }
    }
    let refinement = ResonanceControls {
        absolute_frequency_tolerance_rad_s: controls.absolute_frequency_tolerance_rad_s,
        relative_frequency_tolerance: controls.relative_frequency_tolerance,
        maximum_evaluations: 4096,
        measure_half_power_bandwidth: false,
    };
    let mut peaks = Vec::with_capacity(brackets.len());
    for (left, right) in brackets {
        let before = sampling.evaluations;
        let mut peak = search_resonance(
            (grid[left].omega, grid[right].omega),
            refinement,
            &mut cancelled,
            |omega, _, cancelled| {
                let sample = match grid.binary_search_by(|sample| sample.omega.total_cmp(&omega)) {
                    Ok(index) => grid[index],
                    Err(_) => sampling.query(omega, *cancelled)?,
                };
                Ok((sample.value, sample.derivative, sample.residual))
            },
        )
        .map_err(flatten)?;
        peak.evaluations = sampling.evaluations - before;
        peaks.push(peak);
    }
    search_poll(&mut cancelled)?;
    Ok(ResonanceSurvey {
        peaks,
        interval_rad_s: (
            frequencies_rad_s[0],
            frequencies_rad_s[frequencies_rad_s.len() - 1],
        ),
        grid_points: grid.len(),
        evaluations: sampling.evaluations,
        maximum_relative_residual: sampling.maximum_relative_residual,
    })
}

struct Sampling<F> {
    evaluate: F,
    maximum_evaluations: usize,
    evaluations: usize,
    maximum_relative_residual: f64,
}

impl<F> Sampling<F> {
    fn query<E, C>(
        &mut self,
        omega: f64,
        cancelled: &mut C,
    ) -> Result<GridSample, ResonanceSearchError<E>>
    where
        C: FnMut() -> bool,
        F: FnMut(f64, bool, &mut C) -> Result<(C64, C64, f64), E>,
    {
        search_poll(cancelled)?;
        if self.evaluations == self.maximum_evaluations {
            return Err(ResonanceSearchError::EvaluationLimit {
                maximum: self.maximum_evaluations,
            });
        }
        self.evaluations += 1;
        let (value, derivative, residual) =
            (self.evaluate)(omega, true, cancelled).map_err(ResonanceSearchError::Source)?;
        search_poll(cancelled)?;
        if !finite(value) || !finite(derivative) || !residual.is_finite() || residual < 0.0 {
            return Err(ResonanceSearchError::Invalid(
                "finite physical observation/tangent and nonnegative backward residual required",
            ));
        }
        let slope = if value == C64::ZERO {
            0.0
        } else {
            (derivative / value).re
        };
        if !slope.is_finite() || !value.abs().is_finite() {
            return Err(ResonanceSearchError::Invalid(
                "observed amplitude or logarithmic slope is not representable",
            ));
        }
        self.maximum_relative_residual = self.maximum_relative_residual.max(residual);
        Ok(GridSample {
            omega,
            value,
            derivative,
            residual,
            slope,
        })
    }
}

// The callback also enforces the survey-wide budget. Unwrap exactly that layer
// so an original source error never becomes Source(Source(error)).
fn flatten<E>(error: ResonanceSearchError<ResonanceSearchError<E>>) -> ResonanceSearchError<E> {
    match error {
        ResonanceSearchError::Source(error) => error,
        ResonanceSearchError::Invalid(what) => ResonanceSearchError::Invalid(what),
        ResonanceSearchError::NotBracketed(what) => ResonanceSearchError::NotBracketed(what),
        ResonanceSearchError::EvaluationLimit { maximum } => {
            ResonanceSearchError::EvaluationLimit { maximum }
        }
        ResonanceSearchError::FrequencyResolution => ResonanceSearchError::FrequencyResolution,
        ResonanceSearchError::Cancelled => ResonanceSearchError::Cancelled,
    }
}
