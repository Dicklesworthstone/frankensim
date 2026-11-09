//! Runtime-independent refinement of one actual complex response extremum.
//!
//! This is the shared search owner for callback and executor-context entries.
//! It retains original evaluation refusals, strict numerical sign brackets,
//! bounded callback counts and the actual half-power line shape.

use fs_math::c64::C64;

/// One fixed-source, fixed-receiver transfer and its complete derivative.
/// Frequency-dependent sources, receivers and loads must include their slopes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicTransferSample {
    /// Complex observed amplitude, under `exp(-i omega t)`.
    pub amplitude: C64,
    /// Derivative of that same amplitude per rad/s.
    pub angular_derivative: C64,
}

impl HarmonicTransferSample {
    // Only its sign is used. Independent normalization avoids overflow from
    // squaring a high-Q response or multiplying its still larger derivative.
    fn slope(self) -> f64 {
        let a = component_scale(self.amplitude);
        let d = component_scale(self.angular_derivative);
        if a == 0.0 || d == 0.0 {
            return 0.0;
        }
        (self.amplitude.re / a).mul_add(
            self.angular_derivative.re / d,
            (self.amplitude.im / a) * (self.angular_derivative.im / d),
        )
    }

    fn finite_maximum_resolved(self, width: f64) -> bool {
        // A real-axis pole also changes the amplitude slope from + to -.
        // For H ~ 1/(omega - pole), |Re(H'/H)| * width is at least 2 at
        // the bracket midpoint; at a finite maximum it tends to zero.
        // The conservative 1/8 check is independent of response units and
        // phase rotation. It is not a forward-amplitude error certificate.
        self.amplitude != C64::ZERO
            && (self.angular_derivative / self.amplitude).re.abs() * width <= 0.125
    }
}

/// Frequency resolution and total callback-evaluation ceiling for one query.
#[derive(Clone, Copy, Debug)]
pub struct ResonanceControls {
    /// Absolute angular-frequency bracket width [rad/s], nonnegative.
    pub absolute_angular_tolerance_rad_s: f64,
    /// Relative bracket width, in [0, 1). At least one tolerance is positive.
    /// Resolution is the larger of this times the upper frequency and the
    /// absolute tolerance. Refinement also stops at irreducible floating-point
    /// spacing; the returned bracket makes that resolution visible. Neither
    /// tolerance certifies source/model accuracy. Maxima also require a resolved
    /// finite log-amplitude slope, so a narrow peak can require further refinement.
    pub relative_angular_tolerance: f64,
    /// Includes endpoints, interior samples and the final returned sample.
    /// A bandwidth query shares this ceiling across both flanks.
    pub maximum_evaluations: usize,
}

/// Classification of the local magnitude-squared extremum.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TransferExtremumKind {
    /// Positive amplitude slope on the left, negative on the right.
    Maximum,
    /// Negative amplitude slope on the left, positive on the right.
    Minimum,
}

/// A refined actual response peak or notch, not a free-decay pole.
#[derive(Clone, Copy, Debug)]
pub struct TransferExtremum {
    pub kind: TransferExtremumKind,
    pub angular_frequency_rad_s: f64,
    pub sample: HarmonicTransferSample,
    /// Final oriented numerical sign bracket [rad/s].
    pub bracket_rad_s: (f64, f64),
    /// Original evaluator calls. In a survey, this counts only new refinement
    /// calls; shared grid samples belong to the survey's total.
    pub evaluations: usize,
}

/// Observable peaks and notches resolved by a caller's frequency grid.
/// This is a driven-response survey, not a census of free vibration modes.
#[derive(Clone, Debug)]
pub struct TransferSurvey {
    /// Ascending frequency, retaining both maxima and minima of the same
    /// complex physical observation. No linewidth is inferred across neighbours.
    pub extrema: Vec<TransferExtremum>,
    /// First and last supplied angular frequency [rad/s].
    pub interval_rad_s: (f64, f64),
    /// Number of supplied frequencies, all evaluated on a successful query.
    pub grid_points: usize,
    /// Total original evaluator calls, including the grid and all refinements.
    /// Reusing a cached grid sample consumes neither a call nor this budget.
    pub evaluations: usize,
}

/// Width between two specified crossings of half the peak's `|H|^2`.
/// This is actual sensor-response bandwidth; it is not necessarily the modal
/// damping Q, a Lorentzian fit, or half the system's total dissipated power.
#[derive(Clone, Copy, Debug)]
pub struct HalfPowerBandwidth {
    pub lower_angular_frequency_rad_s: f64,
    pub upper_angular_frequency_rad_s: f64,
    pub angular_bandwidth_rad_s: f64,
    /// Peak frequency / measured bandwidth, not an inferred damping ratio.
    pub response_quality_factor: f64,
    pub lower_bracket_rad_s: (f64, f64),
    pub upper_bracket_rad_s: (f64, f64),
    pub evaluations: usize,
}

/// The original physical evaluator's error is retained without replacement.
#[derive(Debug, PartialEq)]
pub enum ResonanceError<E> {
    Invalid(&'static str),
    InvalidBracket(&'static str),
    /// An exactly stationary interior sample could not be isolated between
    /// strictly oriented slopes. A flat interval is not a strict resonance.
    UnresolvedStationaryInterval,
    Evaluation(E),
    EvaluationBudget {
        maximum_evaluations: usize,
    },
    /// The grid resolves more extrema than the caller admits. No truncated
    /// list is published, and no refinement is started before this check.
    ExtremumBudget {
        maximum_extrema: usize,
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for ResonanceError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-MODAL-RESONANCE: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ResonanceError<E> {}

/// Refine ONE isolated extremum of the supplied actual complex response.
/// Strict endpoint slopes are required: a zero-slope endpoint has only
/// one-sided evidence and refuses. A sampled lossless pole propagates its
/// evaluator's singularity error; an unresolved divergent slope refuses at
/// floating-point spacing. Neither becomes a finite response peak.
///
/// `cancelled` returning true refuses before publication. It is polled before
/// and after each evaluator call and during refinement. Long evaluator work
/// remains responsible for its own internal cancellation boundaries.
pub fn refine_transfer_extremum_with_cancel<E>(
    bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
) -> Result<TransferExtremum, ResonanceError<E>> {
    let mut query = Query::new(controls, cancelled, evaluate)?;
    valid_interval(bracket_rad_s)?;
    let left = query.sample(bracket_rad_s.0, HarmonicTransferSample::slope)?;
    let right = query.sample(bracket_rad_s.1, HarmonicTransferSample::slope)?;
    let kind = match (left.level > 0.0, right.level > 0.0) {
        (true, false) if right.level < 0.0 => TransferExtremumKind::Maximum,
        (false, true) if left.level < 0.0 => TransferExtremumKind::Minimum,
        _ => {
            return Err(ResonanceError::InvalidBracket(
                "strict opposite amplitude slopes required",
            ));
        }
    };
    let (point, bracket_rad_s) = query.refine(
        left,
        right,
        HarmonicTransferSample::slope,
        kind == TransferExtremumKind::Maximum,
    )?;
    Ok(TransferExtremum {
        kind,
        angular_frequency_rad_s: point.frequency,
        sample: point.sample,
        bracket_rad_s,
        evaluations: query.evaluations,
    })
}

/// Discover observable resonance peaks AND antiresonance notches over a band.
///
/// The supplied grid has at least two finite, strictly increasing, nonnegative
/// angular frequencies. All evaluations use the same original complex transfer
/// and its complete frequency tangent, including source, receiver and load
/// dispersion. Strict slope sign changes identify candidates; stationary grid
/// samples are retained between their two oriented neighbours. In particular,
/// an exact transfer zero can be a minimum and is never divided out.
///
/// Every candidate is refined by the same owner as
/// [`refine_transfer_extremum_with_cancel`]. The entire input grid is cached,
/// avoiding repeated physical solves for shared endpoints or stationary samples.
/// `controls.maximum_evaluations` covers the grid and ALL refinements together;
/// it also bounds cached storage before allocation. `maximum_extrema` must be
/// positive. Source refusal, singularity, cancellation, a flat unresolved
/// candidate or either exhausted budget publishes no partial result.
///
/// This finds features RESOLVED BY THE GRID, not every physical resonance.
/// A peak/notch pair can lie between samples with equal slope signs, and an
/// opposite-sign bracket need not contain just one extremum. The caller owns
/// that resolution and isolation obligation. Band endpoints are not two-sided
/// extrema; a silent transfer returns an empty survey. No pole is made finite
/// with invented damping. Each returned maximum can separately be passed to
/// [`refine_half_power_bandwidth_with_cancel`] with isolated flank brackets.
pub fn survey_transfer_extrema_with_cancel<E>(
    frequencies_rad_s: &[f64],
    controls: ResonanceControls,
    maximum_extrema: usize,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
) -> Result<TransferSurvey, ResonanceError<E>> {
    let mut query = Query::new(controls, cancelled, evaluate)?;
    if frequencies_rad_s.len() < 2
        || frequencies_rad_s.iter().any(|w| !w.is_finite() || *w < 0.0)
        || frequencies_rad_s.windows(2).any(|w| w[0] >= w[1])
        || maximum_extrema == 0
    {
        return Err(ResonanceError::Invalid(
            "ordered nonnegative frequency grid and positive extremum capacity required",
        ));
    }
    if frequencies_rad_s.len() > controls.maximum_evaluations {
        return Err(ResonanceError::EvaluationBudget {
            maximum_evaluations: controls.maximum_evaluations,
        });
    }
    query.grid.reserve_exact(frequencies_rad_s.len());
    for &frequency in frequencies_rad_s {
        let sample = query.sample(frequency, HarmonicTransferSample::slope)?;
        query.grid.push(sample);
    }

    // A zero slope alone is not evidence of an extremum. Skip it until two
    // strictly oriented neighbours establish either a peak or a notch. This
    // retains exact antiresonances and avoids double-counting a grid maximum.
    let mut previous: Option<usize> = None;
    let mut brackets = Vec::new();
    for index in 0..query.grid.len() {
        query.poll()?;
        let point = query.grid[index];
        if point.level == 0.0 {
            continue;
        }
        if let Some(left) = previous
            && query.grid[left].level.is_sign_positive() != point.level.is_sign_positive()
        {
            if brackets.len() == maximum_extrema {
                return Err(ResonanceError::ExtremumBudget { maximum_extrema });
            }
            brackets.push((left, index));
        }
        previous = Some(index);
    }

    let mut extrema = Vec::with_capacity(brackets.len());
    for (left, right) in brackets {
        let left = query.grid[left];
        let right = query.grid[right];
        let kind = if left.level > 0.0 {
            TransferExtremumKind::Maximum
        } else {
            TransferExtremumKind::Minimum
        };
        let before = query.evaluations;
        let (point, bracket_rad_s) = query.refine(
            left,
            right,
            HarmonicTransferSample::slope,
            kind == TransferExtremumKind::Maximum,
        )?;
        extrema.push(TransferExtremum {
            kind,
            angular_frequency_rad_s: point.frequency,
            sample: point.sample,
            bracket_rad_s,
            evaluations: query.evaluations - before,
        });
    }
    query.poll()?;
    Ok(TransferSurvey {
        extrema,
        interval_rad_s: (
            frequencies_rad_s[0],
            frequencies_rad_s[frequencies_rad_s.len() - 1],
        ),
        grid_points: frequencies_rad_s.len(),
        evaluations: query.evaluations,
    })
}

/// Refine a separately bracketed crossing on each side of this response peak.
/// Use the SAME immutable physical transfer/source/receiver as for `peak`.
/// The lower bracket must rise through half `|H_peak|^2`; the upper must fall
/// through it. The caller isolates each crossing from neighbouring resonances.
/// Missing crossings, a DC/zero peak, a notch and a singular evaluator refuse.
/// Cancellation has the same callback semantics as
/// [`refine_transfer_extremum_with_cancel`]; one evaluation cap covers both flanks.
pub fn refine_half_power_bandwidth_with_cancel<E>(
    peak: &TransferExtremum,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
) -> Result<HalfPowerBandwidth, ResonanceError<E>> {
    refine_half_power_bandwidth_samples_with_cancel(
        peak,
        lower_bracket_rad_s,
        upper_bracket_rad_s,
        controls,
        cancelled,
        evaluate,
    )
    .map(|(bandwidth, _)| bandwidth)
}

// Keep the final original samples for complete parameter-jet consistency checks.
// This shares the same search, budget and constant scratch as the public query.
pub(super) fn refine_half_power_bandwidth_samples_with_cancel<E>(
    peak: &TransferExtremum,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cancelled: impl FnMut() -> bool,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
) -> Result<(HalfPowerBandwidth, [HarmonicTransferSample; 2]), ResonanceError<E>> {
    let mut query = Query::new(controls, cancelled, evaluate)?;
    valid_interval(lower_bracket_rad_s)?;
    valid_interval(upper_bracket_rad_s)?;
    if peak.kind != TransferExtremumKind::Maximum
        || !peak.angular_frequency_rad_s.is_finite()
        || peak.angular_frequency_rad_s <= 0.0
        || !finite(peak.sample.amplitude)
        || component_scale(peak.sample.amplitude) == 0.0
        || lower_bracket_rad_s.1 > peak.angular_frequency_rad_s
        || upper_bracket_rad_s.0 < peak.angular_frequency_rad_s
    {
        return Err(ResonanceError::Invalid(
            "positive finite response peak and disjoint flanks required",
        ));
    }
    let level = |sample: HarmonicTransferSample| {
        // Common scaling retains the level sign at extreme response units,
        // including when directly forming |H_peak|^2 would overflow/underflow.
        let scale = component_scale(sample.amplitude).max(component_scale(peak.sample.amplitude));
        let h = C64::new(sample.amplitude.re / scale, sample.amplitude.im / scale);
        let p = C64::new(
            peak.sample.amplitude.re / scale,
            peak.sample.amplitude.im / scale,
        );
        h.norm_sq() - 0.5 * p.norm_sq()
    };
    let l0 = query.sample(lower_bracket_rad_s.0, level)?;
    let l1 = query.sample(lower_bracket_rad_s.1, level)?;
    let u0 = query.sample(upper_bracket_rad_s.0, level)?;
    let u1 = query.sample(upper_bracket_rad_s.1, level)?;
    if !(l0.level < 0.0 && l1.level > 0.0 && u0.level > 0.0 && u1.level < 0.0) {
        return Err(ResonanceError::InvalidBracket(
            "each flank must strictly bracket its half-power crossing",
        ));
    }
    let (lower, lower_bracket_rad_s) = query.refine(l0, l1, level, false)?;
    let (upper, upper_bracket_rad_s) = query.refine(u0, u1, level, false)?;
    if lower_bracket_rad_s.1 >= upper_bracket_rad_s.0 {
        // A coarse frequency tolerance can stop both flanks while they still
        // share the peak. Their midpoint difference is then a numerical step,
        // not a resolved positive linewidth of a narrow physical resonance.
        return Err(ResonanceError::InvalidBracket(
            "half-power crossing brackets do not resolve a positive linewidth",
        ));
    }
    let width = upper.frequency - lower.frequency;
    let quality = peak.angular_frequency_rad_s / width;
    if width <= 0.0 || !quality.is_finite() {
        return Err(ResonanceError::Invalid(
            "response linewidth or quality factor is not representable",
        ));
    }
    Ok((
        HalfPowerBandwidth {
            lower_angular_frequency_rad_s: lower.frequency,
            upper_angular_frequency_rad_s: upper.frequency,
            angular_bandwidth_rad_s: width,
            response_quality_factor: quality,
            lower_bracket_rad_s,
            upper_bracket_rad_s,
            evaluations: query.evaluations,
        },
        [lower.sample, upper.sample],
    ))
}

fn component_scale(z: C64) -> f64 {
    z.re.abs().max(z.im.abs())
}

pub(super) fn valid_interval<E>((lower, upper): (f64, f64)) -> Result<(), ResonanceError<E>> {
    if !lower.is_finite() || !upper.is_finite() || lower < 0.0 || upper <= lower {
        return Err(ResonanceError::Invalid(
            "finite increasing nonnegative frequency interval required",
        ));
    }
    Ok(())
}

// The refinement engine operates on original physical samples. Transfer
// magnitudes retain their finite-maximum guard; scalar observables can perform
// their own required curvature validation before publishing a result.
pub(super) trait QuerySample: Copy {
    const NONFINITE_MESSAGE: &'static str;
    fn finite(self) -> bool;
    fn finite_maximum_resolved(self, _width: f64) -> bool {
        true
    }
}

impl QuerySample for HarmonicTransferSample {
    const NONFINITE_MESSAGE: &'static str =
        "finite physical transfer and frequency derivative required";

    fn finite(self) -> bool {
        finite(self.amplitude) && finite(self.angular_derivative)
    }

    fn finite_maximum_resolved(self, width: f64) -> bool {
        HarmonicTransferSample::finite_maximum_resolved(self, width)
    }
}

#[derive(Clone, Copy)]
pub(super) struct Point<S> {
    pub(super) frequency: f64,
    pub(super) sample: S,
    pub(super) level: f64,
}

pub(super) struct Query<F, C, S> {
    controls: ResonanceControls,
    cancelled: C,
    evaluate: F,
    pub(super) evaluations: usize,
    // Only survey input samples are cached. Ordinary local refinement stays
    // allocation-free; survey memory is O(grid points + returned extrema),
    // independent of the number of source evaluations during refinement.
    grid: Vec<Point<S>>,
}

impl<F, C, S, E> Query<F, C, S>
where
    F: FnMut(f64) -> Result<S, E>,
    C: FnMut() -> bool,
    S: QuerySample,
{
    pub(super) fn new(
        controls: ResonanceControls,
        cancelled: C,
        evaluate: F,
    ) -> Result<Self, ResonanceError<E>> {
        if !controls.absolute_angular_tolerance_rad_s.is_finite()
            || controls.absolute_angular_tolerance_rad_s < 0.0
            || !controls.relative_angular_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_angular_tolerance)
            || controls.absolute_angular_tolerance_rad_s == 0.0
                && controls.relative_angular_tolerance == 0.0
            || controls.maximum_evaluations < 3
        {
            return Err(ResonanceError::Invalid(
                "finite positive frequency resolution and at least three evaluations required",
            ));
        }
        let mut query = Self {
            controls,
            cancelled,
            evaluate,
            evaluations: 0,
            grid: Vec::new(),
        };
        query.poll()?;
        Ok(query)
    }

    pub(super) fn poll(&mut self) -> Result<(), ResonanceError<E>> {
        if (self.cancelled)() {
            Err(ResonanceError::Cancelled)
        } else {
            Ok(())
        }
    }

    pub(super) fn sample(
        &mut self,
        frequency: f64,
        level: impl Fn(S) -> f64,
    ) -> Result<Point<S>, ResonanceError<E>> {
        self.poll()?;
        if let Ok(index) = self
            .grid
            .binary_search_by(|p| p.frequency.total_cmp(&frequency))
        {
            let mut point = self.grid[index];
            point.level = level(point.sample);
            return Ok(point);
        }
        if self.evaluations == self.controls.maximum_evaluations {
            return Err(ResonanceError::EvaluationBudget {
                maximum_evaluations: self.controls.maximum_evaluations,
            });
        }
        self.evaluations += 1;
        let sample = (self.evaluate)(frequency).map_err(ResonanceError::Evaluation)?;
        self.poll()?;
        if !sample.finite() {
            return Err(ResonanceError::Invalid(S::NONFINITE_MESSAGE));
        }
        Ok(Point {
            frequency,
            sample,
            level: level(sample),
        })
    }

    pub(super) fn refine(
        &mut self,
        mut left: Point<S>,
        mut right: Point<S>,
        level: impl Fn(S) -> f64 + Copy,
        finite_maximum: bool,
    ) -> Result<(Point<S>, (f64, f64)), ResonanceError<E>> {
        // An exactly stationary response stays between two shrinking sign
        // endpoints. Its original physical solve is still valid there; keep
        // just that sample instead of refactoring the same frequency on every
        // iteration. The cache is local to this level/transfer and O(1).
        let mut stationary: Option<Point<S>> = None;
        loop {
            self.poll()?;
            let width = right.frequency - left.frequency;
            let middle = left.frequency + 0.5 * width;
            let at_spacing = middle == left.frequency || middle == right.frequency;
            let point = if middle == left.frequency {
                left
            } else if middle == right.frequency {
                right
            } else if let Some(point) = stationary.filter(|p| p.frequency == middle) {
                point
            } else {
                self.sample(middle, level)?
            };
            let finite_response = !finite_maximum || point.sample.finite_maximum_resolved(width);
            let tolerance = self
                .controls
                .absolute_angular_tolerance_rad_s
                .max(self.controls.relative_angular_tolerance * right.frequency);
            if (width <= tolerance || at_spacing) && finite_response {
                return Ok((point, (left.frequency, right.frequency)));
            }
            if at_spacing {
                return Err(ResonanceError::InvalidBracket(
                    "finite response peak is unresolved at floating-point frequency spacing",
                ));
            }
            if point.level == 0.0 {
                stationary = Some(point);
                // Do not confuse a flat response or stationary inflection with
                // a peak. Keep two-sided strict orientation while shrinking
                // around an exactly stationary/crossing numerical sample.
                let a_frequency = left.frequency + 0.5 * (middle - left.frequency);
                let b_frequency = middle + 0.5 * (right.frequency - middle);
                if a_frequency == left.frequency
                    || a_frequency == middle
                    || b_frequency == middle
                    || b_frequency == right.frequency
                {
                    if !finite_response {
                        return Err(ResonanceError::InvalidBracket(
                            "finite response peak is unresolved at floating-point frequency spacing",
                        ));
                    }
                    return Ok((point, (left.frequency, right.frequency)));
                }
                let a = self.sample(a_frequency, level)?;
                let b = self.sample(b_frequency, level)?;
                if a.level == 0.0
                    || b.level == 0.0
                    || a.level.is_sign_positive() != left.level.is_sign_positive()
                    || b.level.is_sign_positive() != right.level.is_sign_positive()
                {
                    return Err(ResonanceError::UnresolvedStationaryInterval);
                }
                left = a;
                right = b;
            } else if point.level.is_sign_positive() == left.level.is_sign_positive() {
                left = point;
            } else {
                right = point;
            }
        }
    }
}

fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
