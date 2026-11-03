//! Normal arrival from original mechanical endpoint responses.
//!
//! Every query starts from the same accepted body states. Only actual sampled
//! nonpenetrating endpoints can be returned; a secant estimate is never a body
//! state. A successful result is freshly replayed, retaining its original gap,
//! normal velocity and represented duration without projection.
//!
//! The bracket describes the supplied discrete response branch. It does not
//! certify the earliest exact-flow collision, a unique nonlinear root, or an
//! admissible impact/support handoff. The receiving owner must check every
//! contact and retain its original impulse, momentum and work gates.

use std::panic::{AssertUnwindSafe, catch_unwind};

/// Temporal, geometric and evaluation limits for one normal-arrival search.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArrivalControls {
    /// Smallest duration the original response may evaluate, in s.
    pub minimum_step_s: f64,
    /// Maximum bracket width for an approximate arrival, in s.
    pub time_tolerance_s: f64,
    /// Maximum actual nonnegative final gap, in m. Zero demands an exact
    /// numerical zero; it does not authorize rounding a position or gap.
    pub gap_tolerance_m: f64,
    /// New response evaluations, including the independent final replay.
    /// The supplied upper observation is already owned by the caller.
    pub maximum_trials: usize,
}

/// An actual replayed endpoint, still subject to all-contact physical checks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ArrivalReport {
    pub elapsed_s: f64,
    /// Duration coordinates from the same accepted predecessor. A point is
    /// reported only when the complete replay returns an exactly zero gap.
    pub bracket_s: [f64; 2],
    pub gap_m: f64,
    /// Original signed relative normal velocity; positive means separation.
    pub normal_velocity_m_s: f64,
    pub response_evaluations: usize,
}

#[derive(Debug)]
pub enum ArrivalError<E> {
    Invalid(&'static str),
    NonFinite(&'static str),
    Cancelled,
    Response(E),
    ResponsePanicked,
    NotBracketed {
        upper_gap_m: f64,
    },
    MinimumStep {
        bracket_s: [f64; 2],
        minimum_step_s: f64,
    },
    Unrepresentable {
        bracket_s: [f64; 2],
        requested_s: f64,
        actual_s: f64,
    },
    Budget {
        evaluations: usize,
        maximum: usize,
    },
    Replay {
        expected_elapsed_s: f64,
        actual_elapsed_s: f64,
        gap_m: f64,
    },
}

impl<E> ArrivalError<E> {
    /// Convert only the original-response error, preserving numerical refusal
    /// details. A receiving owner can first intercept its private typed earlier
    /// event marker and then map ordinary physical failures to its public error.
    pub fn map_response<F>(self, map: impl FnOnce(E) -> F) -> ArrivalError<F> {
        match self {
            Self::Invalid(reason) => ArrivalError::Invalid(reason),
            Self::NonFinite(operation) => ArrivalError::NonFinite(operation),
            Self::Cancelled => ArrivalError::Cancelled,
            Self::Response(error) => ArrivalError::Response(map(error)),
            Self::ResponsePanicked => ArrivalError::ResponsePanicked,
            Self::NotBracketed { upper_gap_m } => ArrivalError::NotBracketed { upper_gap_m },
            Self::MinimumStep {
                bracket_s,
                minimum_step_s,
            } => ArrivalError::MinimumStep {
                bracket_s,
                minimum_step_s,
            },
            Self::Unrepresentable {
                bracket_s,
                requested_s,
                actual_s,
            } => ArrivalError::Unrepresentable {
                bracket_s,
                requested_s,
                actual_s,
            },
            Self::Budget {
                evaluations,
                maximum,
            } => ArrivalError::Budget {
                evaluations,
                maximum,
            },
            Self::Replay {
                expected_elapsed_s,
                actual_elapsed_s,
                gap_m,
            } => ArrivalError::Replay {
                expected_elapsed_s,
                actual_elapsed_s,
                gap_m,
            },
        }
    }
}
impl<E: core::fmt::Display> core::fmt::Display for ArrivalError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Invalid(reason) => write!(f, "invalid normal-arrival search: {reason}"),
            Self::NonFinite(operation) => {
                write!(f, "non-finite normal-arrival operation: {operation}")
            }
            Self::Cancelled => write!(f, "normal-arrival search cancelled"),
            Self::Response(error) => write!(f, "original normal-arrival response: {error}"),
            Self::ResponsePanicked => write!(f, "original normal-arrival response panicked"),
            Self::NotBracketed { upper_gap_m } => write!(
                f,
                "normal arrival has no signed bracket: upper gap {upper_gap_m} m"
            ),
            Self::MinimumStep {
                bracket_s,
                minimum_step_s,
            } => write!(
                f,
                "normal-arrival bracket {bracket_s:?} s cannot resolve minimum step {minimum_step_s} s"
            ),
            Self::Unrepresentable {
                bracket_s,
                requested_s,
                actual_s,
            } => write!(
                f,
                "normal-arrival bracket {bracket_s:?} s cannot resolve request {requested_s} s (actual {actual_s} s)"
            ),
            Self::Budget {
                evaluations,
                maximum,
            } => write!(
                f,
                "normal-arrival response budget exhausted after {evaluations}/{maximum} calls"
            ),
            Self::Replay {
                expected_elapsed_s,
                actual_elapsed_s,
                gap_m,
            } => write!(
                f,
                "normal-arrival replay at {actual_elapsed_s} s (expected {expected_elapsed_s} s) has inadmissible gap {gap_m} m"
            ),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ArrivalError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        if let Self::Response(error) = self {
            Some(error)
        } else {
            None
        }
    }
}

#[derive(Clone, Copy)]
struct Sample {
    elapsed_s: f64,
    gap_m: f64,
    normal_velocity_m_s: f64,
}
impl Sample {
    fn new<E>(
        elapsed_s: f64,
        gap_m: f64,
        normal_velocity_m_s: f64,
    ) -> Result<Self, ArrivalError<E>> {
        if !positive(elapsed_s) || !gap_m.is_finite() || !normal_velocity_m_s.is_finite() {
            return Err(ArrivalError::NonFinite(
                "finite original endpoint at positive represented time required",
            ));
        }
        Ok(Self {
            elapsed_s,
            gap_m,
            normal_velocity_m_s,
        })
    }
}

/// Locate a normal arrival from a positive initial gap and a nonpositive upper.
///
/// `response` returns (actual represented duration, actual gap, actual signed
/// normal velocity). Its mechanical trajectory and source offsets must always
/// begin at the same unchanged accepted predecessor. The callback may retain
/// its complete native candidate; on success the independent replay is the last
/// callback. Its actual gap and velocity are returned verbatim.
///
/// An owner-specific error can signal a real earlier event, such as a tangential
/// stop. It is returned immediately as [`ArrivalError::Response`], without more
/// queries. Ordinary represented-clock rounding is an ordinary response value;
/// it must never be relabeled as an earlier physical event by this locator.
///
/// Secant trials alternate with bisection. Approximate success requires both a
/// narrow actual bracket and an original lower endpoint within the gap bound.
/// The known initial gap is never returned as a zero-duration step. A zero-gap
/// outward departure needs its own positive-gap bracket, not this entry point.
#[allow(clippy::too_many_arguments)]
pub(super) fn locate<E>(
    initial_gap_m: f64,
    upper_elapsed_s: f64,
    upper_gap_m: f64,
    upper_normal_velocity_m_s: f64,
    controls: ArrivalControls,
    response: &mut impl FnMut(f64, &mut dyn FnMut() -> bool) -> Result<(f64, f64, f64), E>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ArrivalReport, ArrivalError<E>> {
    poll(cancelled)?;
    if !positive(initial_gap_m)
        || !positive(controls.minimum_step_s)
        || !positive(controls.time_tolerance_s)
        || !controls.gap_tolerance_m.is_finite()
        || controls.gap_tolerance_m < 0.0
        || controls.maximum_trials == 0
        || !positive(upper_elapsed_s)
        || upper_elapsed_s < controls.minimum_step_s
    {
        return Err(ArrivalError::Invalid(
            "positive initial gap, finite controls, bounded work and positive original interval required",
        ));
    }
    let mut hi = Sample::new(upper_elapsed_s, upper_gap_m, upper_normal_velocity_m_s)?;
    let mut evaluations = 0;
    if hi.gap_m > 0.0 {
        return Err(ArrivalError::NotBracketed {
            upper_gap_m: hi.gap_m,
        });
    }
    if hi.gap_m == 0.0 {
        return replay(
            hi,
            [hi.elapsed_s; 2],
            true,
            controls,
            response,
            cancelled,
            &mut evaluations,
        );
    }
    // The initial gap is known, but no synthetic initial velocity or native
    // zero-duration candidate is manufactured for that bracket endpoint.
    let mut lo: Option<Sample> = None;
    let mut bisect = false;
    loop {
        poll(cancelled)?;
        let lo_s = lo.map_or(0.0, |sample| sample.elapsed_s);
        let lo_gap = lo.map_or(initial_gap_m, |sample| sample.gap_m);
        let bracket = [lo_s, hi.elapsed_s];
        let width = hi.elapsed_s - lo_s;
        if let Some(lower) = lo {
            if width <= controls.time_tolerance_s && lower.gap_m <= controls.gap_tolerance_m {
                return replay(
                    lower,
                    bracket,
                    false,
                    controls,
                    response,
                    cancelled,
                    &mut evaluations,
                );
            }
        }
        if hi.elapsed_s <= controls.minimum_step_s {
            return Err(ArrivalError::MinimumStep {
                bracket_s: bracket,
                minimum_step_s: controls.minimum_step_s,
            });
        }
        let mut request = f64::midpoint(lo_s, hi.elapsed_s);
        let mut used_bisection = true;
        if !bisect {
            // Normalize before adding opposite-sign gaps so a valid finite
            // bracket does not overflow just while forming its interpolation.
            let magnitude = lo_gap.max(-hi.gap_m);
            let positive = lo_gap / magnitude;
            let negative = -hi.gap_m / magnitude;
            let fraction = positive / (positive + negative);
            let secant = width.mul_add(fraction, lo_s);
            if secant > lo_s && secant < hi.elapsed_s {
                request = secant;
                used_bisection = false;
            }
        }
        request = request.max(controls.minimum_step_s);
        if !(request > lo_s && request < hi.elapsed_s) {
            return Err(ArrivalError::Unrepresentable {
                bracket_s: bracket,
                requested_s: request,
                actual_s: request,
            });
        }
        let sample = evaluate(request, controls, response, cancelled, &mut evaluations)?;
        if sample.elapsed_s < controls.minimum_step_s
            || sample.elapsed_s < lo_s
            || sample.elapsed_s > hi.elapsed_s
        {
            return Err(ArrivalError::Unrepresentable {
                bracket_s: bracket,
                requested_s: request,
                actual_s: sample.elapsed_s,
            });
        }
        if sample.gap_m == 0.0 {
            return replay(
                sample,
                [sample.elapsed_s; 2],
                true,
                controls,
                response,
                cancelled,
                &mut evaluations,
            );
        }
        if sample.elapsed_s == lo_s || sample.elapsed_s == hi.elapsed_s {
            // A close secant may round to an existing native clock. One
            // midpoint attempt distinguishes that from a bracket with no room.
            if used_bisection {
                return Err(ArrivalError::Unrepresentable {
                    bracket_s: bracket,
                    requested_s: request,
                    actual_s: sample.elapsed_s,
                });
            }
            bisect = true;
            continue;
        }
        if sample.gap_m > 0.0 {
            lo = Some(sample);
        } else {
            hi = sample;
        }
        bisect = !used_bisection;
    }
}

fn replay<E>(
    expected: Sample,
    bracket_s: [f64; 2],
    exact: bool,
    controls: ArrivalControls,
    response: &mut impl FnMut(f64, &mut dyn FnMut() -> bool) -> Result<(f64, f64, f64), E>,
    cancelled: &mut impl FnMut() -> bool,
    evaluations: &mut usize,
) -> Result<ArrivalReport, ArrivalError<E>> {
    let actual = evaluate(
        expected.elapsed_s,
        controls,
        response,
        cancelled,
        evaluations,
    )?;
    if actual.elapsed_s.to_bits() != expected.elapsed_s.to_bits()
        || actual.gap_m != expected.gap_m
        || actual.normal_velocity_m_s != expected.normal_velocity_m_s
        || actual.gap_m < 0.0
        || actual.gap_m > controls.gap_tolerance_m
        || (exact && actual.gap_m != 0.0)
    {
        return Err(ArrivalError::Replay {
            expected_elapsed_s: expected.elapsed_s,
            actual_elapsed_s: actual.elapsed_s,
            gap_m: actual.gap_m,
        });
    }
    poll(cancelled)?;
    Ok(ArrivalReport {
        elapsed_s: actual.elapsed_s,
        bracket_s,
        gap_m: actual.gap_m,
        normal_velocity_m_s: actual.normal_velocity_m_s,
        response_evaluations: *evaluations,
    })
}

fn evaluate<E>(
    duration_s: f64,
    controls: ArrivalControls,
    response: &mut impl FnMut(f64, &mut dyn FnMut() -> bool) -> Result<(f64, f64, f64), E>,
    cancelled: &mut impl FnMut() -> bool,
    evaluations: &mut usize,
) -> Result<Sample, ArrivalError<E>> {
    poll(cancelled)?;
    if *evaluations >= controls.maximum_trials {
        return Err(ArrivalError::Budget {
            evaluations: *evaluations,
            maximum: controls.maximum_trials,
        });
    }
    *evaluations += 1;
    let mut cancellation_seen = false;
    let returned = catch_unwind(AssertUnwindSafe(|| {
        response(duration_s, &mut || {
            cancellation_seen = cancellation_seen || cancelled();
            cancellation_seen
        })
    }));
    if cancellation_seen {
        return Err(ArrivalError::Cancelled);
    }
    let (actual_s, gap, velocity) = returned
        .map_err(|_| ArrivalError::ResponsePanicked)?
        .map_err(ArrivalError::Response)?;
    poll(cancelled)?;
    Sample::new(actual_s, gap, velocity)
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), ArrivalError<E>> {
    if cancelled() {
        Err(ArrivalError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "arrival/tests.rs"]
mod tests;
