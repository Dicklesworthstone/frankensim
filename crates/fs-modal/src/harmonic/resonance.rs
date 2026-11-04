//! Resolve peaks and notches of an actual complex forced transfer.
//!
//! The stationary condition is `Re(conj(H) * dH/domega) = 0`. It retains
//! interference between modes, source phase, nonproportional damping and any
//! dispersive load already present in the evaluator. It is not a pole search:
//! displacement, velocity and pressure can peak at different frequencies.
//!
//! The caller supplies a continuous differentiable transfer and brackets ONE
//! isolated extremum or level crossing in each requested interval. Opposite
//! endpoint signs are checked; uniqueness and absence of unresolved poles are
//! not certified. Refinement is local, without a grid, spectral completeness
//! claim, damping floor or an assumed Lorentzian line shape. Each evaluation
//! keeps its original solve/domain/work limits, and cancellation publishes no
//! partial result. The returned bracket is numerical sign evidence, not an
//! interval enclosure accounting for the evaluator's forward error.

use super::{HarmonicFrequencyResponse, HarmonicResponseError, finite, invalid};
use fs_exec::Cx;
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
}

impl HarmonicFrequencyResponse {
    /// Observe a fixed real linear displacement receiver, including signed
    /// cancellation between modes. This uses the complete solved frequency
    /// tangent, not natural frequencies or modal amplitude magnitudes.
    pub fn displacement_transfer(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        project(
            receiver,
            &self.response.displacement,
            &self.displacement_derivative,
        )
    }

    /// Observe a fixed real linear velocity receiver. The derivative includes
    /// both `-i q` and `-i omega q_omega`, so its peak need not equal the
    /// displacement peak of the same structure and forcing.
    pub fn velocity_transfer(
        &self,
        receiver: &[f64],
    ) -> Result<HarmonicTransferSample, HarmonicResponseError> {
        project(receiver, &self.response.velocity, &self.velocity_derivative)
    }
}

fn project(
    receiver: &[f64],
    values: &[C64],
    derivatives: &[C64],
) -> Result<HarmonicTransferSample, HarmonicResponseError> {
    if receiver.len() != values.len()
        || derivatives.len() != values.len()
        || receiver.is_empty()
        || receiver.iter().any(|x| !x.is_finite())
    {
        return Err(invalid(
            "one finite fixed receiver weight per coordinate required",
        ));
    }
    let mut sums = [0.0; 4];
    let mut correction = [0.0; 4];
    for ((&weight, &value), &derivative) in receiver.iter().zip(values).zip(derivatives) {
        let terms = [
            weight * value.re,
            weight * value.im,
            weight * derivative.re,
            weight * derivative.im,
        ];
        for ((sum, error), term) in sums.iter_mut().zip(&mut correction).zip(terms) {
            let next = *sum + term;
            *error += if sum.abs() >= term.abs() {
                (*sum - next) + term
            } else {
                (term - next) + *sum
            };
            *sum = next;
        }
    }
    let amplitude = C64::new(sums[0] + correction[0], sums[1] + correction[1]);
    let angular_derivative = C64::new(sums[2] + correction[2], sums[3] + correction[3]);
    if !finite(amplitude) || !finite(angular_derivative) {
        return Err(invalid(
            "physical transfer or its derivative is not representable",
        ));
    }
    Ok(HarmonicTransferSample {
        amplitude,
        angular_derivative,
    })
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
    /// tolerance certifies source/model accuracy.
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
/// one-sided evidence and refuses. A lossless pole propagates its evaluator's
/// singularity error; it is never regularized into a finite response peak.
pub fn refine_transfer_extremum<E>(
    bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
) -> Result<TransferExtremum, ResonanceError<E>> {
    let mut query = Query::new(controls, cx, evaluate)?;
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
    let (point, bracket_rad_s) = query.refine(left, right, HarmonicTransferSample::slope)?;
    Ok(TransferExtremum {
        kind,
        angular_frequency_rad_s: point.frequency,
        sample: point.sample,
        bracket_rad_s,
        evaluations: query.evaluations,
    })
}

/// Refine a separately bracketed crossing on each side of this response peak.
/// Use the SAME immutable physical transfer/source/receiver as for `peak`.
/// The lower bracket must rise through half `|H_peak|^2`; the upper must fall
/// through it. The caller isolates each crossing from neighbouring resonances.
/// Missing crossings, a DC/zero peak, a notch and a singular evaluator refuse.
pub fn refine_half_power_bandwidth<E>(
    peak: &TransferExtremum,
    lower_bracket_rad_s: (f64, f64),
    upper_bracket_rad_s: (f64, f64),
    controls: ResonanceControls,
    cx: &Cx<'_>,
    evaluate: impl FnMut(f64) -> Result<HarmonicTransferSample, E>,
) -> Result<HalfPowerBandwidth, ResonanceError<E>> {
    let mut query = Query::new(controls, cx, evaluate)?;
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
    let (lower, lower_bracket_rad_s) = query.refine(l0, l1, level)?;
    let (upper, upper_bracket_rad_s) = query.refine(u0, u1, level)?;
    let width = upper.frequency - lower.frequency;
    let quality = peak.angular_frequency_rad_s / width;
    if width <= 0.0 || !quality.is_finite() {
        return Err(ResonanceError::Invalid(
            "response linewidth or quality factor is not representable",
        ));
    }
    Ok(HalfPowerBandwidth {
        lower_angular_frequency_rad_s: lower.frequency,
        upper_angular_frequency_rad_s: upper.frequency,
        angular_bandwidth_rad_s: width,
        response_quality_factor: quality,
        lower_bracket_rad_s,
        upper_bracket_rad_s,
        evaluations: query.evaluations,
    })
}

fn component_scale(z: C64) -> f64 {
    z.re.abs().max(z.im.abs())
}

fn valid_interval<E>((lower, upper): (f64, f64)) -> Result<(), ResonanceError<E>> {
    if !lower.is_finite() || !upper.is_finite() || lower < 0.0 || upper <= lower {
        return Err(ResonanceError::Invalid(
            "finite increasing nonnegative frequency interval required",
        ));
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct Point {
    frequency: f64,
    sample: HarmonicTransferSample,
    level: f64,
}

struct Query<'a, 'cx, F> {
    controls: ResonanceControls,
    cx: &'a Cx<'cx>,
    evaluate: F,
    evaluations: usize,
}

impl<'a, 'cx, F, E> Query<'a, 'cx, F>
where
    F: FnMut(f64) -> Result<HarmonicTransferSample, E>,
{
    fn new(
        controls: ResonanceControls,
        cx: &'a Cx<'cx>,
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
        let query = Self {
            controls,
            cx,
            evaluate,
            evaluations: 0,
        };
        query.poll()?;
        Ok(query)
    }

    fn poll(&self) -> Result<(), ResonanceError<E>> {
        self.cx.checkpoint().map_err(|_| ResonanceError::Cancelled)
    }

    fn sample(
        &mut self,
        frequency: f64,
        level: impl Fn(HarmonicTransferSample) -> f64,
    ) -> Result<Point, ResonanceError<E>> {
        self.poll()?;
        if self.evaluations == self.controls.maximum_evaluations {
            return Err(ResonanceError::EvaluationBudget {
                maximum_evaluations: self.controls.maximum_evaluations,
            });
        }
        self.evaluations += 1;
        let sample = (self.evaluate)(frequency).map_err(ResonanceError::Evaluation)?;
        self.poll()?;
        if !finite(sample.amplitude) || !finite(sample.angular_derivative) {
            return Err(ResonanceError::Invalid(
                "finite physical transfer and frequency derivative required",
            ));
        }
        Ok(Point {
            frequency,
            sample,
            level: level(sample),
        })
    }

    fn refine(
        &mut self,
        mut left: Point,
        mut right: Point,
        level: impl Fn(HarmonicTransferSample) -> f64 + Copy,
    ) -> Result<(Point, (f64, f64)), ResonanceError<E>> {
        loop {
            self.poll()?;
            let middle = left.frequency + 0.5 * (right.frequency - left.frequency);
            if middle == left.frequency || middle == right.frequency {
                return Ok((
                    if middle == left.frequency {
                        left
                    } else {
                        right
                    },
                    (left.frequency, right.frequency),
                ));
            }
            let point = self.sample(middle, level)?;
            let tolerance = self
                .controls
                .absolute_angular_tolerance_rad_s
                .max(self.controls.relative_angular_tolerance * right.frequency);
            if right.frequency - left.frequency <= tolerance {
                return Ok((point, (left.frequency, right.frequency)));
            }
            if point.level == 0.0 {
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
