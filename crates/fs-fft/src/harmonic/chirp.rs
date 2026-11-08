//! Joint amplitude and frequency-law estimation for linear chirps.

use super::refinement::{self, Controls, Domain};
use super::{C64, FitDesign, HarmonicFit, HarmonicFitBudget, HarmonicFitError, HarmonicFitSummary};
use fs_math::dd::Dd;

/// One linear instantaneous-frequency law across the sampled record.
///
/// With `T = (sample_count - 1) * timestep_s`, the phase in cycles is
/// `start_frequency_hz*t + (end_frequency_hz-start_frequency_hz)*t*t/(2*T)`.
/// Frequencies refer to the first and LAST samples, not to FFT bin edges or
/// the unsampled endpoint `sample_count * timestep_s`. Either drift sign is
/// allowed. Equal endpoints give a stationary tone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChirpedTone {
    /// Positive instantaneous frequency at the first sample [Hz].
    pub start_frequency_hz: f64,
    /// Positive instantaneous frequency at the last sample [Hz].
    pub end_frequency_hz: f64,
}

/// Joint DC/cosine/sine fit at supplied linear frequency trajectories.
///
/// All phases start at the first sample. The peak phasor `a + i*b` describes
/// `a*cos(phase) + b*sin(phase)` under `exp(-i phase)` in original input units.
/// Tone order is retained; trajectories may cross or share one endpoint.
/// Their complete sampled design must remain numerically distinguishable.
/// No frequency bins, resampling, or piecewise stationary windows are used.
///
/// This supplied-law plan does not estimate frequency, drift, damping, model
/// order or uncertainty. Nonlinear drift, envelopes and omitted tones remain
/// in the residual and may bias the fit. Preparation and reuse have the same
/// work/storage orders as [`super::ToneFitPlan`]; `fit_into` allocates nothing.
/// Zero drift preserves that plan's design and fit bits for the same ordering.
#[derive(Debug, Clone)]
pub struct ChirpedToneFitPlan {
    pub(super) tones: Vec<ChirpedTone>,
    pub(super) design: FitDesign,
}

impl ChirpedToneFitPlan {
    /// Prepare distinct linear chirps on one uniform physical clock.
    /// Both endpoint frequencies must be positive and strictly below Nyquist;
    /// linear interpolation then admits the entire frequency trajectory.
    /// At least `1 + 2 * tones.len()` samples and a finite duration are required.
    /// Retained endpoint pairs are included in the storage budget.
    pub fn new(
        sample_count: usize,
        timestep_s: f64,
        tones: &[ChirpedTone],
        budget: HarmonicFitBudget,
    ) -> Result<Self, HarmonicFitError> {
        if tones.is_empty() || sample_count < 2 {
            return Err(super::invalid(
                "chirp fitting needs tones and at least two samples",
            ));
        }
        let duration = (sample_count - 1) as f64 * timestep_s;
        if !(timestep_s.is_finite() && timestep_s > 0.0 && duration.is_finite() && duration > 0.0) {
            return Err(super::invalid(
                "chirp clock and duration must be positive and finite",
            ));
        }
        for (k, tone) in tones.iter().enumerate() {
            if ![tone.start_frequency_hz, tone.end_frequency_hz]
                .into_iter()
                .all(|f| f.is_finite() && f > 0.0 && f * timestep_s > 0.0 && f * timestep_s < 0.5)
            {
                return Err(super::invalid(
                    "chirp endpoint frequencies must be positive, finite and below Nyquist",
                ));
            }
            if tones[..k].contains(tone) {
                return Err(super::invalid(
                    "chirp frequency trajectories must be distinct",
                ));
            }
        }
        // The shared owner already budgets one retained scalar per tone.
        let extra_bytes = tones
            .len()
            .checked_mul(core::mem::size_of::<f64>())
            .ok_or_else(|| super::invalid("chirp endpoint storage size overflows"))?;
        let inverse_intervals = Dd::ONE / Dd::from_f64((sample_count - 1) as f64);
        let design = FitDesign::with_basis(
            sample_count,
            timestep_s,
            tones.len(),
            |k| tones[k].start_frequency_hz.max(tones[k].end_frequency_hz),
            |j, k| sample_sin_cos(tones[k], timestep_s, inverse_intervals, j),
            extra_bytes,
            budget,
        )?;
        Ok(Self {
            tones: tones.to_vec(),
            design,
        })
    }

    /// Supplied endpoint frequencies in the same order as fitted phasors.
    #[must_use]
    pub fn tones(&self) -> &[ChirpedTone] {
        &self.tones
    }

    /// Required record and scratch length.
    #[must_use]
    pub const fn sample_count(&self) -> usize {
        self.design.samples
    }

    /// Uniform physical sample interval [s].
    #[must_use]
    pub const fn timestep_s(&self) -> f64 {
        self.design.timestep_s
    }

    /// Reciprocal 1-norm condition of the complete linear amplitude design.
    #[must_use]
    pub const fn reciprocal_condition(&self) -> f64 {
        self.design.reciprocal_condition
    }

    /// Conservative peak numeric-array storage admitted at preparation.
    #[must_use]
    pub const fn workspace_bytes(&self) -> usize {
        self.design.workspace_bytes
    }

    /// Fit one record, allocating scratch and first-sample peak phasors.
    pub fn fit(&self, samples: &[f64]) -> Result<HarmonicFit, HarmonicFitError> {
        self.design.fit(samples)
    }

    /// Allocation-free fit; output phasors remain unchanged on refusal.
    pub fn fit_into(
        &self,
        samples: &[f64],
        scratch: &mut [f64],
        phasors: &mut [C64],
    ) -> Result<HarmonicFitSummary, HarmonicFitError> {
        self.design.fit_into(samples, scratch, phasors)
    }
}

// Form the integrated phase as a convex combination of endpoint frequencies,
// not a rounded drift rate times a large time squared. DD clock arithmetic and
// cycle reduction retain the original phase before the scalar trig call.
pub(super) fn sample_sin_cos(
    tone: ChirpedTone,
    dt: f64,
    inverse_intervals: Dd,
    index: usize,
) -> (f64, f64) {
    if tone.start_frequency_hz == tone.end_frequency_hz {
        return super::sample_sin_cos(tone.start_frequency_hz, dt, index);
    }
    let index = Dd::from_f64(index as f64);
    let half_fraction = index * inverse_intervals * Dd::from_f64(0.5);
    let average = Dd::from_f64(tone.start_frequency_hz) * (Dd::ONE - half_fraction)
        + Dd::from_f64(tone.end_frequency_hz) * half_fraction;
    super::sin_cos_cycles(average * Dd::from_f64(dt) * index)
}

/// Closed local rectangle for one chirp's first/last-sample frequencies [Hz].
/// Rectangles must be strictly disjoint in at least one endpoint coordinate.
/// This retains tone identity while allowing their trajectories to cross.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChirpedToneBounds {
    /// Positive, increasing search interval at the first sample [Hz].
    pub start_frequency_hz: [f64; 2],
    /// Positive, increasing search interval at the last sample [Hz].
    pub end_frequency_hz: [f64; 2],
}

/// Numerical resolution, work and storage admitted for local chirp inference.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChirpedToneRefinementBudget {
    /// Every untruncated endpoint-frequency correction must meet this limit.
    pub absolute_frequency_tolerance_hz: f64,
    /// All candidate trajectory sets, including rejected backtracking trials.
    pub max_evaluations: usize,
    /// Minimum reciprocal condition of both the full variable-projection
    /// Jacobian and its projected physical sensitivity. Coordinates are each
    /// endpoint frequency times duration, without independent normalization.
    /// Resolved sensitivity is required in the original common signal units;
    /// residual curvature cannot identify an absent tone. This is not uncertainty.
    pub min_chirp_reciprocal_condition: f64,
    /// Amplitude-design conditioning and total peak numeric-array storage,
    /// including refinement scratch and both candidates' channel fits.
    pub fit: HarmonicFitBudget,
}

/// Signed, untruncated endpoint corrections [Hz], not error bars.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChirpedToneCorrection {
    /// Signed first-sample frequency correction [Hz].
    pub start_frequency_hz: f64,
    /// Signed last-sample frequency correction [Hz].
    pub end_frequency_hz: f64,
}

/// Jointly inferred frequency trajectories and first-sample amplitudes.
#[derive(Debug, Clone, PartialEq)]
pub struct ChirpedToneFit {
    /// Inferred endpoint frequencies in the supplied seed order.
    pub tones: Vec<ChirpedTone>,
    /// DC, first-sample peak phasors and unexplained RMS in input units.
    pub fit: HarmonicFit,
    /// Evaluated candidate sets, including rejected trials.
    pub evaluations: usize,
    /// Remaining untruncated endpoint corrections in tone order.
    pub last_corrections: Vec<ChirpedToneCorrection>,
    /// Minimum condition of the full and projected frequency-law Jacobians.
    pub chirp_reciprocal_condition: f64,
}

/// Shared frequency trajectories with independent simultaneous channel phasors.
#[derive(Debug, Clone, PartialEq)]
pub struct ChirpedToneChannelsFit {
    /// Original seed order is retained, including when trajectories cross.
    pub tones: Vec<ChirpedTone>,
    /// Original channel order, input units, DC and first-sample phasors.
    pub fits: Vec<HarmonicFit>,
    /// Evaluated candidate sets; each shares one amplitude QR across channels.
    pub evaluations: usize,
    /// Remaining untruncated endpoint corrections in tone order.
    pub last_corrections: Vec<ChirpedToneCorrection>,
    /// Minimum condition of the full and projected frequency-law Jacobians.
    pub chirp_reciprocal_condition: f64,
}

/// Refusal to infer unresolved drift or publish unconverged chirp amplitudes.
#[derive(Debug, Clone, PartialEq)]
pub enum ChirpedToneRefinementError {
    /// The original amplitude owner refused data, conditioning or storage.
    Fit(HarmonicFitError),
    /// Invalid clock, domains, channel shapes or numerical controls.
    InvalidParameter {
        /// Admission condition that failed.
        what: &'static str,
    },
    /// Frequency-law derivatives lack resolved physical signal or rank.
    Unidentifiable {
        /// Refused full/projected Jacobian reciprocal 1-norm condition.
        reciprocal_condition: f64,
        /// Requested minimum reciprocal condition.
        minimum: f64,
    },
    /// Work exhaustion, active bounds or stagnation left excessive corrections.
    NotConverged {
        /// Best accepted trajectory set, without partial fitted amplitudes.
        tones: Vec<ChirpedTone>,
        /// Joint unexplained RMS in the channels' common physical units.
        residual_rms: f64,
        /// Evaluated candidate sets, including rejected trials.
        evaluations: usize,
        /// Remaining untruncated endpoint corrections at that candidate.
        last_corrections: Vec<ChirpedToneCorrection>,
    },
}

impl core::fmt::Display for ChirpedToneRefinementError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Fit(error) => error.fmt(f),
            Self::InvalidParameter { what } => write!(f, "chirped tones: {what}"),
            Self::Unidentifiable {
                reciprocal_condition,
                minimum,
            } => write!(
                f,
                "chirped tones are unidentifiable: reciprocal condition {reciprocal_condition:e} is below {minimum:e}"
            ),
            Self::NotConverged { evaluations, .. } => write!(
                f,
                "chirped tones did not converge after {evaluations} candidate sets"
            ),
        }
    }
}
impl std::error::Error for ChirpedToneRefinementError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Fit(error) => Some(error),
            _ => None,
        }
    }
}

fn tones(parameters: &[f64]) -> Vec<ChirpedTone> {
    parameters
        .chunks_exact(2)
        .map(|p| ChirpedTone {
            start_frequency_hz: p[0],
            end_frequency_hz: p[1],
        })
        .collect()
}
fn corrections(parameters: &[f64]) -> Vec<ChirpedToneCorrection> {
    parameters
        .chunks_exact(2)
        .map(|p| ChirpedToneCorrection {
            start_frequency_hz: p[0],
            end_frequency_hz: p[1],
        })
        .collect()
}
impl From<refinement::Error> for ChirpedToneRefinementError {
    fn from(error: refinement::Error) -> Self {
        match error {
            refinement::Error::Fit(error) => Self::Fit(error),
            refinement::Error::InvalidParameter { what } => Self::InvalidParameter { what },
            refinement::Error::Unidentifiable {
                reciprocal_condition,
                minimum,
            } => Self::Unidentifiable {
                reciprocal_condition,
                minimum,
            },
            refinement::Error::NotConverged {
                parameters,
                residual_rms,
                evaluations,
                corrections: delta,
            } => Self::NotConverged {
                tones: tones(&parameters),
                residual_rms,
                evaluations,
                last_corrections: corrections(&delta),
            },
        }
    }
}
fn invalid(what: &'static str) -> ChirpedToneRefinementError {
    ChirpedToneRefinementError::InvalidParameter { what }
}

/// Jointly refine locally seeded linear frequency trajectories and amplitudes.
///
/// Both instantaneous endpoint frequencies are inferred in hertz. Using
/// endpoints avoids poorly scaled frequency/Hz-per-second coordinates and
/// makes full-record Nyquist admission explicit. Frequency drift follows from
/// their difference divided by `(samples.len()-1)*timestep_s`.
///
/// The original shared QR optimizer retains the COMPLETE variable-projection
/// derivative `P_perp*D*c + (A^+)^T*D^T*r`, including the residual-dependent
/// term. Projected sensitivities are separately gated so residual curvature
/// cannot manufacture an observable chirp. Backtracking is monotone, every
/// attempted candidate is budgeted, and endpoint steps move the phase by at
/// most a quarter cycle across the record. Active bounds are not convergence.
///
/// This is local inference of a declared linear law and model order. It does
/// not discover tones, fit arbitrary frequency tracks or amplitude envelopes,
/// certify uncertainty, remove aliasing, or extrapolate beyond this record.
/// Nonlinear drift and unmodeled components remain in the residual and may
/// bias the result. These time-varying frequencies are not stationary modal poles.
pub fn refine_chirped_tones(
    samples: &[f64],
    timestep_s: f64,
    initial_tones: &[ChirpedTone],
    bounds: &[ChirpedToneBounds],
    budget: ChirpedToneRefinementBudget,
) -> Result<ChirpedToneFit, ChirpedToneRefinementError> {
    let mut result =
        refine_chirped_tones_channels(&[samples], timestep_s, initial_tones, bounds, budget)?;
    Ok(ChirpedToneFit {
        tones: result.tones,
        fit: result.fits.remove(0),
        evaluations: result.evaluations,
        last_corrections: result.last_corrections,
        chirp_reciprocal_condition: result.chirp_reciprocal_condition,
    })
}

/// Infer one set of trajectories from simultaneous common-clock channels.
/// Each channel keeps its own DC/phasors; original COMMON-unit amplitude weights
/// survive orthogonal frame changes. Constant channels retain their DC without
/// supplying frequency information. Convert unrelated units explicitly first.
///
/// One amplitude QR serves all channels of a candidate, and derivative blocks
/// fold through small QRs instead of normal equations or a channel-stack matrix.
/// Work is O(E C (N K² + K³)); numeric storage is O(N K + K² + C K).
/// This offline owner makes no mid-factorization cancellation claim.
pub fn refine_chirped_tones_channels(
    channels: &[&[f64]],
    timestep_s: f64,
    initial_tones: &[ChirpedTone],
    bounds: &[ChirpedToneBounds],
    budget: ChirpedToneRefinementBudget,
) -> Result<ChirpedToneChannelsFit, ChirpedToneRefinementError> {
    if initial_tones.is_empty() || initial_tones.len() != bounds.len() {
        return Err(invalid(
            "each chirped tone needs one endpoint search rectangle",
        ));
    }
    if !(timestep_s.is_finite() && timestep_s > 0.0)
        || !(budget.absolute_frequency_tolerance_hz.is_finite()
            && budget.absolute_frequency_tolerance_hz > 0.0)
        || budget.max_evaluations == 0
        || !(budget.min_chirp_reciprocal_condition.is_finite()
            && budget.min_chirp_reciprocal_condition > 0.0
            && budget.min_chirp_reciprocal_condition <= 1.0)
    {
        return Err(invalid(
            "invalid chirp clock, frequency tolerance, work or condition limit",
        ));
    }
    for (k, (seed, bound)) in initial_tones.iter().zip(bounds).enumerate() {
        for (value, [lo, hi]) in [
            (seed.start_frequency_hz, bound.start_frequency_hz),
            (seed.end_frequency_hz, bound.end_frequency_hz),
        ] {
            if ![value, lo, hi].into_iter().all(f64::is_finite)
                || !(lo > 0.0 && lo < hi && lo * timestep_s > 0.0 && hi * timestep_s < 0.5)
                || !(value >= lo && value <= hi)
            {
                return Err(invalid(
                    "chirp rectangles must contain their seeds and remain positive and sub-Nyquist at both endpoints",
                ));
            }
        }
        if bounds[..k].iter().any(|other| {
            !(bound.start_frequency_hz[1] < other.start_frequency_hz[0]
                || other.start_frequency_hz[1] < bound.start_frequency_hz[0]
                || bound.end_frequency_hz[1] < other.end_frequency_hz[0]
                || other.end_frequency_hz[1] < bound.end_frequency_hz[0])
        }) {
            return Err(invalid(
                "chirp rectangles must be strictly disjoint in at least one endpoint frequency",
            ));
        }
    }
    let result = refinement::refine_channels(
        channels,
        timestep_s,
        Domain::Chirped(initial_tones, bounds),
        Controls {
            tolerances: [budget.absolute_frequency_tolerance_hz; 2],
            max_evaluations: budget.max_evaluations,
            min_reciprocal_condition: budget.min_chirp_reciprocal_condition,
            fit: budget.fit,
        },
    )?;
    Ok(ChirpedToneChannelsFit {
        tones: tones(&result.parameters),
        fits: result.fits,
        evaluations: result.evaluations,
        last_corrections: corrections(&result.corrections),
        chirp_reciprocal_condition: result.reciprocal_condition,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_phase_retains_fractional_cycles_at_a_billion_sample_clock() {
        let last = 1_usize << 30;
        let p = ChirpedTone {
            start_frequency_hz: 1.0,
            end_frequency_hz: 3.0,
        };
        let (s, c) = sample_sin_cos(p, 0.125, Dd::ONE / Dd::from_f64(last as f64), last - 1);
        // Exact rational phase: i/8 + i²/(8L) at i=L-1 is
        // L/4 - 3/8 + 1/(8L). Remove its integer part independently.
        let fraction = -0.375 + 1.0 / (8.0 * last as f64);
        let (expected_s, expected_c) = (core::f64::consts::TAU * fraction).sin_cos();
        assert!((s - expected_s).abs() < 2e-15);
        assert!((c - expected_c).abs() < 2e-15);
    }
}
