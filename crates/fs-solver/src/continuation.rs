//! Adaptive parameter continuation through an existing nonlinear corrector.
//!
//! This owner schedules bounded trials on a normalized parameter in [0,1].
//! The caller supplies every nonlinear solve and its independent acceptance.
//! Rejected trials never become seeds. No new residual, tangent, root algorithm,
//! physical time, stability selection or pseudo-arclength equation is introduced.

/// Explicit parameter resolution and number of nonlinear attempts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContinuationControls {
    pub initial_increment: f64,
    pub minimum_increment: f64,
    pub maximum_increment: f64,
    /// Includes refused correctors, excluding the caller's initial root at zero.
    pub maximum_attempts: usize,
}
impl ContinuationControls {
    pub fn validate(self) -> Result<(), &'static str> {
        if !self.minimum_increment.is_finite()
            || !self.initial_increment.is_finite()
            || !self.maximum_increment.is_finite()
            || self.minimum_increment <= 0.0
            || self.minimum_increment > self.initial_increment
            || self.initial_increment > self.maximum_increment
            || self.maximum_increment > 1.0
            || !(1..=256).contains(&self.maximum_attempts)
        {
            Err("finite 0 < minimum <= initial <= maximum <= 1 and 1..=256 attempts required")
        } else {
            Ok(())
        }
    }
}

/// Only the original owner decides whether a numerical refusal may be retried.
#[derive(Clone, Debug, PartialEq)]
pub enum CorrectorRefusal<E> {
    /// A smaller parameter increment may make the same corrector converge.
    Retry(E),
    /// Source, physical domain, control, cancellation or other terminal refusal.
    Stop(E),
}

/// Counts of actual parameter trials; these are not physical trajectory steps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContinuationReport {
    pub attempts: usize,
    pub accepted_steps: usize,
    pub rejected_steps: usize,
    pub smallest_accepted_increment: f64,
}

/// One accepted target root. No partially advanced state accompanies an error.
#[derive(Clone, Debug, PartialEq)]
pub struct Continued<T> {
    pub state: T,
    pub report: ContinuationReport,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ContinuationError<E> {
    Invalid(&'static str),
    Source {
        parameter: f64,
        source: E,
    },
    AttemptLimit {
        accepted_parameter: f64,
        attempts: usize,
    },
    IncrementLimit {
        accepted_parameter: f64,
        attempted_parameter: f64,
        source: E,
    },
    UnrepresentableParameter {
        accepted_parameter: f64,
    },
}
impl<E: core::fmt::Debug> core::fmt::Display for ContinuationError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "parameter continuation: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ContinuationError<E> {}

/// Follow a supplied, already admitted parameter-zero root to exactly one.
///
/// Each corrector receives the last accepted state, never a rejected candidate.
/// Success doubles the next increment up to the explicit maximum; a retryable
/// refusal halves the attempted increment. A final remaining interval may be
/// smaller than the minimum, but a refusal there cannot lower that minimum.
/// Every attempt and every acceptance polls the caller's checkpoint. Callbacks
/// must leave the borrowed accepted state and constitutive history unchanged;
/// their own side effects cannot be rolled back by this scheduling owner.
///
/// Every successful corrector must independently admit its complete target
/// residual. Continuation alone proves neither branch uniqueness nor smoothness,
/// and it does not cross a parameter fold by inventing an arclength constraint.
pub fn follow<T, E>(
    initial: T,
    controls: ContinuationControls,
    mut correct: impl FnMut(f64, &T) -> Result<T, CorrectorRefusal<E>>,
    mut checkpoint: impl FnMut() -> Result<(), E>,
) -> Result<Continued<T>, ContinuationError<E>> {
    checkpoint().map_err(|source| ContinuationError::Source {
        parameter: 0.0,
        source,
    })?;
    controls.validate().map_err(ContinuationError::Invalid)?;
    let mut state = initial;
    let mut parameter = 0.0;
    let mut increment = controls.initial_increment;
    let mut report = ContinuationReport {
        attempts: 0,
        accepted_steps: 0,
        rejected_steps: 0,
        smallest_accepted_increment: 1.0,
    };
    loop {
        checkpoint().map_err(|source| ContinuationError::Source { parameter, source })?;
        if report.attempts == controls.maximum_attempts {
            return Err(ContinuationError::AttemptLimit {
                accepted_parameter: parameter,
                attempts: report.attempts,
            });
        }
        let remaining = 1.0 - parameter;
        // Selecting one explicitly prevents a rounded endpoint from silently
        // replacing the exact caller's final problem by a nearly final problem.
        let target = if increment >= remaining {
            1.0
        } else {
            parameter + increment
        };
        if target <= parameter || target > 1.0 || !target.is_finite() {
            return Err(ContinuationError::UnrepresentableParameter {
                accepted_parameter: parameter,
            });
        }
        let attempted_increment = target - parameter;
        report.attempts += 1;
        match correct(target, &state) {
            Ok(candidate) => {
                checkpoint().map_err(|source| ContinuationError::Source {
                    parameter: target,
                    source,
                })?;
                state = candidate;
                parameter = target;
                report.accepted_steps += 1;
                report.smallest_accepted_increment =
                    report.smallest_accepted_increment.min(attempted_increment);
                if parameter == 1.0 {
                    return Ok(Continued { state, report });
                }
                increment = (2.0 * attempted_increment).min(controls.maximum_increment);
            }
            Err(CorrectorRefusal::Stop(source)) => {
                return Err(ContinuationError::Source {
                    parameter: target,
                    source,
                });
            }
            Err(CorrectorRefusal::Retry(source)) => {
                report.rejected_steps += 1;
                let smaller = 0.5 * attempted_increment;
                if smaller < controls.minimum_increment || smaller == 0.0 {
                    return Err(ContinuationError::IncrementLimit {
                        accepted_parameter: parameter,
                        attempted_parameter: target,
                        source,
                    });
                }
                increment = smaller;
            }
        }
    }
}

#[cfg(test)]
#[path = "continuation/tests.rs"]
mod tests;
