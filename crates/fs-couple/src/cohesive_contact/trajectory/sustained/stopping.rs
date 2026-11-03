//! Locate a stopping endpoint on the original mechanical response branch.
//!
//! The scalar guard is the original tangential velocity projected onto its
//! initial direction. A zero guard alone is not a stop: the complete returned
//! velocity must also satisfy its physical tolerance. Every response starts
//! from the same accepted body states. This module never changes a velocity,
//! constructs a body, or admits a physical trajectory. The receiving owner
//! must still audit all contacts and all original work before adoption.
//!
//! A bracket describes the sampled discrete endpoint response, not an earliest
//! exact-flow event or a uniqueness theorem for nonlinear frictional dynamics.

use std::panic::{AssertUnwindSafe, catch_unwind};

/// Explicit temporal, physical and evaluation limits for one stopping search.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StoppingControls {
    /// Smallest duration the original mechanical owner may evaluate, in s.
    pub minimum_step_s: f64,
    /// Maximum bracket width for an approximate stop, in s.
    pub time_tolerance_s: f64,
    /// Maximum Euclidean norm of the complete relative tangent velocity.
    pub velocity_tolerance_m_s: f64,
    /// All new response calls, including the independent final replay.
    /// The supplied upper observation is not charged a second time here.
    pub maximum_trials: usize,
}

/// An actual replayed endpoint; physical adoption remains the caller's job.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StoppingReport {
    /// Actual represented duration returned by the original response.
    pub elapsed_s: f64,
    /// Signed-guard bracket, or a point only for an exactly zero full velocity.
    pub bracket_s: [f64; 2],
    /// Complete original relative velocity, without projection or zeroing.
    pub velocity_m_s: [f64; 2],
    pub response_evaluations: usize,
}

#[derive(Debug)]
pub enum StoppingError<E> {
    Invalid(&'static str),
    NonFinite(&'static str),
    Cancelled,
    Response(E),
    ResponsePanicked,
    NotSliding {
        initial_speed_m_s: f64,
    },
    NotBracketed {
        upper_projection_m_s: f64,
    },
    /// A scalar directional crossing with unresolved transverse motion.
    TransverseVelocity {
        elapsed_s: f64,
        velocity_m_s: [f64; 2],
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
        velocity_m_s: [f64; 2],
    },
}

impl<E: core::fmt::Display> core::fmt::Display for StoppingError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Response(error) => write!(f, "original stopping response: {error}"),
            Self::Invalid(reason) => write!(f, "invalid stopping search: {reason}"),
            Self::NonFinite(operation) => write!(f, "non-finite stopping operation: {operation}"),
            Self::Cancelled => write!(f, "stopping search cancelled"),
            Self::ResponsePanicked => write!(f, "original stopping response panicked"),
            Self::NotSliding { initial_speed_m_s } => write!(
                f,
                "stopping search starts at near-rest speed {initial_speed_m_s} m/s"
            ),
            Self::NotBracketed {
                upper_projection_m_s,
            } => write!(
                f,
                "stopping response has no signed bracket: upper guard {upper_projection_m_s} m/s"
            ),
            Self::TransverseVelocity {
                elapsed_s,
                velocity_m_s,
            } => write!(
                f,
                "directional crossing at {elapsed_s} s retains transverse velocity {velocity_m_s:?} m/s"
            ),
            Self::MinimumStep {
                bracket_s,
                minimum_step_s,
            } => write!(
                f,
                "stopping bracket {bracket_s:?} s lies below the original minimum step {minimum_step_s} s"
            ),
            Self::Unrepresentable {
                bracket_s,
                requested_s,
                actual_s,
            } => write!(
                f,
                "stopping bracket {bracket_s:?} s cannot resolve request {requested_s} s (actual {actual_s} s)"
            ),
            Self::Budget {
                evaluations,
                maximum,
            } => write!(
                f,
                "stopping response budget exhausted after {evaluations}/{maximum} calls"
            ),
            Self::Replay {
                expected_elapsed_s,
                actual_elapsed_s,
                velocity_m_s,
            } => write!(
                f,
                "stopping replay at {actual_elapsed_s} s (expected {expected_elapsed_s} s) retains velocity {velocity_m_s:?} m/s"
            ),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for StoppingError<E> {
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
    velocity_m_s: [f64; 2],
    projection_m_s: f64,
    speed_m_s: f64,
}
impl Sample {
    fn new<E>(
        elapsed_s: f64,
        velocity_m_s: [f64; 2],
        direction: [f64; 2],
    ) -> Result<Self, StoppingError<E>> {
        if !elapsed_s.is_finite() || velocity_m_s.iter().any(|v| !v.is_finite()) {
            return Err(StoppingError::NonFinite("original endpoint"));
        }
        let speed_m_s = velocity_m_s[0].hypot(velocity_m_s[1]);
        let projection_m_s = direction[0].mul_add(velocity_m_s[0], direction[1] * velocity_m_s[1]);
        if !speed_m_s.is_finite() || !projection_m_s.is_finite() {
            return Err(StoppingError::NonFinite("endpoint speed or guard"));
        }
        Ok(Self {
            elapsed_s,
            velocity_m_s,
            projection_m_s,
            speed_m_s,
        })
    }
}

/// Search a supplied positive-to-nonpositive directional bracket.
///
/// `response` returns the actual represented duration and full relative tangent
/// velocity from a fresh original body solve. Reflected high probes are private
/// branch observations, never accepted physical steps. Secant trials alternate
/// with bisection so one-sided interpolation cannot stagnate indefinitely.
///
/// Approximate stops require both time and full-vector velocity tolerances. A
/// strictly positive upper guard does not establish a bracket even if its
/// speed is tiny. An exact full-vector zero may identify a point bracket. Every
/// successful result comes from a fresh replay, which is the final callback.
/// No native response is requested at the known initial time zero.
#[allow(clippy::too_many_arguments)]
pub(super) fn locate<E>(
    start_velocity: [f64; 2],
    upper_elapsed_s: f64,
    upper_velocity: [f64; 2],
    controls: StoppingControls,
    response: &mut impl FnMut(f64, &mut dyn FnMut() -> bool) -> Result<(f64, [f64; 2]), E>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StoppingReport, StoppingError<E>> {
    poll(cancelled)?;
    if !positive(controls.minimum_step_s)
        || !positive(controls.time_tolerance_s)
        || !positive(controls.velocity_tolerance_m_s)
        || controls.maximum_trials == 0
        || !positive(upper_elapsed_s)
        || upper_elapsed_s < controls.minimum_step_s
        || start_velocity.iter().any(|v| !v.is_finite())
    {
        return Err(StoppingError::Invalid(
            "finite positive controls, duration and complete initial velocity required",
        ));
    }
    let speed = start_velocity[0].hypot(start_velocity[1]);
    if !speed.is_finite() {
        return Err(StoppingError::NonFinite("initial speed"));
    }
    if speed <= controls.velocity_tolerance_m_s {
        return Err(StoppingError::NotSliding {
            initial_speed_m_s: speed,
        });
    }
    let direction = [start_velocity[0] / speed, start_velocity[1] / speed];
    let mut lo = Sample::new(0.0, start_velocity, direction)?;
    let mut hi = Sample::new(upper_elapsed_s, upper_velocity, direction)?;
    let mut evaluations = 0;
    if hi.speed_m_s == 0.0 {
        return replay(
            hi,
            [hi.elapsed_s; 2],
            true,
            direction,
            controls,
            response,
            cancelled,
            &mut evaluations,
        );
    }
    if hi.projection_m_s > 0.0 {
        return Err(StoppingError::NotBracketed {
            upper_projection_m_s: hi.projection_m_s,
        });
    }
    if hi.projection_m_s == 0.0 && hi.speed_m_s > controls.velocity_tolerance_m_s {
        return Err(StoppingError::TransverseVelocity {
            elapsed_s: hi.elapsed_s,
            velocity_m_s: hi.velocity_m_s,
        });
    }
    let mut bisect = false;
    loop {
        poll(cancelled)?;
        let bracket = [lo.elapsed_s, hi.elapsed_s];
        let width = hi.elapsed_s - lo.elapsed_s;
        let best = if lo.elapsed_s >= controls.minimum_step_s && lo.speed_m_s <= hi.speed_m_s {
            lo
        } else {
            hi
        };
        if width <= controls.time_tolerance_s && best.speed_m_s <= controls.velocity_tolerance_m_s {
            return replay(
                best,
                bracket,
                false,
                direction,
                controls,
                response,
                cancelled,
                &mut evaluations,
            );
        }
        if hi.elapsed_s <= controls.minimum_step_s {
            return Err(StoppingError::MinimumStep {
                bracket_s: bracket,
                minimum_step_s: controls.minimum_step_s,
            });
        }
        let midpoint = f64::midpoint(lo.elapsed_s, hi.elapsed_s);
        let mut request = midpoint;
        let mut used_bisection = true;
        if !bisect {
            let magnitude = lo.projection_m_s.max(-hi.projection_m_s);
            let positive = lo.projection_m_s / magnitude;
            let negative = -hi.projection_m_s / magnitude;
            let fraction = positive / (positive + negative);
            let secant = width.mul_add(fraction, lo.elapsed_s);
            if secant > lo.elapsed_s && secant < hi.elapsed_s {
                request = secant;
                used_bisection = false;
            }
        }
        request = request.max(controls.minimum_step_s);
        if !(request > lo.elapsed_s && request < hi.elapsed_s) {
            return Err(StoppingError::Unrepresentable {
                bracket_s: bracket,
                requested_s: request,
                actual_s: request,
            });
        }
        let sample = evaluate(
            request,
            direction,
            controls,
            response,
            cancelled,
            &mut evaluations,
        )?;
        if sample.elapsed_s < controls.minimum_step_s
            || sample.elapsed_s < lo.elapsed_s
            || sample.elapsed_s > hi.elapsed_s
        {
            return Err(StoppingError::Unrepresentable {
                bracket_s: bracket,
                requested_s: request,
                actual_s: sample.elapsed_s,
            });
        }
        if sample.speed_m_s == 0.0 {
            return replay(
                sample,
                [sample.elapsed_s; 2],
                true,
                direction,
                controls,
                response,
                cancelled,
                &mut evaluations,
            );
        }
        if sample.elapsed_s == lo.elapsed_s || sample.elapsed_s == hi.elapsed_s {
            // A secant close to an endpoint may round to the same native
            // clock. Try a midpoint once before declaring no represented room.
            if used_bisection {
                return Err(StoppingError::Unrepresentable {
                    bracket_s: bracket,
                    requested_s: request,
                    actual_s: sample.elapsed_s,
                });
            }
            bisect = true;
            continue;
        }
        if sample.projection_m_s == 0.0 && sample.speed_m_s > controls.velocity_tolerance_m_s {
            return Err(StoppingError::TransverseVelocity {
                elapsed_s: sample.elapsed_s,
                velocity_m_s: sample.velocity_m_s,
            });
        }
        if sample.projection_m_s > 0.0 {
            lo = sample;
        } else {
            hi = sample;
        }
        bisect = !used_bisection;
    }
}

#[allow(clippy::too_many_arguments)]
fn replay<E>(
    expected: Sample,
    bracket_s: [f64; 2],
    exact: bool,
    direction: [f64; 2],
    controls: StoppingControls,
    response: &mut impl FnMut(f64, &mut dyn FnMut() -> bool) -> Result<(f64, [f64; 2]), E>,
    cancelled: &mut impl FnMut() -> bool,
    evaluations: &mut usize,
) -> Result<StoppingReport, StoppingError<E>> {
    let actual = evaluate(
        expected.elapsed_s,
        direction,
        controls,
        response,
        cancelled,
        evaluations,
    )?;
    if actual.elapsed_s.to_bits() != expected.elapsed_s.to_bits()
        || actual.speed_m_s > controls.velocity_tolerance_m_s
        || (exact && actual.speed_m_s != 0.0)
    {
        return Err(StoppingError::Replay {
            expected_elapsed_s: expected.elapsed_s,
            actual_elapsed_s: actual.elapsed_s,
            velocity_m_s: actual.velocity_m_s,
        });
    }
    poll(cancelled)?;
    Ok(StoppingReport {
        elapsed_s: actual.elapsed_s,
        bracket_s,
        velocity_m_s: actual.velocity_m_s,
        response_evaluations: *evaluations,
    })
}

fn evaluate<E>(
    duration_s: f64,
    direction: [f64; 2],
    controls: StoppingControls,
    response: &mut impl FnMut(f64, &mut dyn FnMut() -> bool) -> Result<(f64, [f64; 2]), E>,
    cancelled: &mut impl FnMut() -> bool,
    evaluations: &mut usize,
) -> Result<Sample, StoppingError<E>> {
    poll(cancelled)?;
    if *evaluations >= controls.maximum_trials {
        return Err(StoppingError::Budget {
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
        return Err(StoppingError::Cancelled);
    }
    let (actual_s, velocity) = returned
        .map_err(|_| StoppingError::ResponsePanicked)?
        .map_err(StoppingError::Response)?;
    poll(cancelled)?;
    if !positive(actual_s) {
        return Err(StoppingError::Invalid(
            "original response must advance positive finite time",
        ));
    }
    Sample::new(actual_s, velocity, direction)
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), StoppingError<E>> {
    if cancelled() {
        Err(StoppingError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
#[path = "stopping/tests.rs"]
mod tests;
