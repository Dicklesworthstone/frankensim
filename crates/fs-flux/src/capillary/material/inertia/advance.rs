//! Adaptive complete intervals over the same event-resolved inertial owner.
//!
//! Only actual conservative half-step candidates are accepted. Stop times are
//! temporal-error observables as well as control boundaries: an early stop
//! shared by coarse and fine trials must not hide unresolved pre-stop motion.

use super::{
    InertialCapillaryBoundary, InertialCapillaryIntervalControls, MaterialInertiaError,
    MaterialInertialCapillaryStep, SourceInertialCapillary, Sum, accumulate_interval, mechanics,
    poll, validate_boundary,
};
use fs_time::PiController;

/// Independent position, velocity and stopping-time error scales.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapillaryInertiaAdvanceControls {
    pub initial_step_s: f64,
    /// Smallest trial interval; the accepted half steps can be half this size.
    /// A final remainder may be shorter to reach the requested time exactly.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    pub maximum_attempts: usize,
    pub relative_tolerance: f64,
    pub length_absolute_tolerance_m: f64,
    pub velocity_absolute_tolerance_m_s: f64,
    pub event_time_absolute_tolerance_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CapillaryInertiaAdvanceError {
    Invalid(&'static str),
    Cancelled,
    Step {
        attempted_step_s: f64,
        error: MaterialInertiaError,
    },
    Attempts {
        attempts: usize,
        accepted_intervals: usize,
    },
    Tolerance {
        normalized_error: f64,
        step_s: f64,
    },
    EventResolution {
        event_time_s: f64,
        minimum_step_s: f64,
    },
    StepUnderflow,
}
impl core::fmt::Display for CapillaryInertiaAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive source inertial capillary: {self:?}")
    }
}
impl std::error::Error for CapillaryInertiaAdvanceError {}

/// Bounded-work counts include rejected trials; physical transfers do not.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CapillaryInertiaAdvanceReport {
    pub attempts: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub step_evaluations: usize,
    pub maximum_normalized_local_error: f64,
}

#[derive(Clone, Debug)]
pub struct CapillaryInertiaAdvanceCandidate<'a, 'd> {
    endpoint: SourceInertialCapillary<'a, 'd>,
    step: MaterialInertialCapillaryStep,
    report: CapillaryInertiaAdvanceReport,
}
impl<'a, 'd> CapillaryInertiaAdvanceCandidate<'a, 'd> {
    #[must_use]
    pub const fn endpoint(&self) -> &SourceInertialCapillary<'a, 'd> {
        &self.endpoint
    }
    #[must_use]
    pub const fn step(&self) -> &MaterialInertialCapillaryStep {
        &self.step
    }
    #[must_use]
    pub const fn report(&self) -> &CapillaryInertiaAdvanceReport {
        &self.report
    }
    #[must_use]
    pub fn into_parts(
        self,
    ) -> (
        SourceInertialCapillary<'a, 'd>,
        MaterialInertialCapillaryStep,
        CapillaryInertiaAdvanceReport,
    ) {
        (self.endpoint, self.step, self.report)
    }
}

impl<'a, 'd> SourceInertialCapillary<'a, 'd> {
    /// Advance the complete constant-boundary interval or retain the old owner.
    /// Position, velocity and resolved first-stop time are compared between one
    /// full step and two half steps. Smooth intervals use the midpoint factor
    /// of three; detected stop/branch transitions use the full discrepancy.
    /// Estimates control local error and are not a global error enclosure.
    pub fn advance_candidate(
        &self,
        input: InertialCapillaryBoundary,
        solve: InertialCapillaryIntervalControls,
        controls: CapillaryInertiaAdvanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CapillaryInertiaAdvanceCandidate<'a, 'd>, CapillaryInertiaAdvanceError> {
        use CapillaryInertiaAdvanceError as Error;
        let step_error = |h, error| match error {
            MaterialInertiaError::Cancelled => Error::Cancelled,
            other => Error::Step {
                attempted_step_s: h,
                error: other,
            },
        };
        poll(&mut cancelled).map_err(|e| step_error(0.0, e))?;
        validate_boundary(input, solve).map_err(|e| step_error(0.0, e))?;
        validate(controls)?;
        if input.duration_s == 0.0 {
            let candidate = self
                .step_candidate(input, solve, &mut cancelled)
                .map_err(|e| step_error(0.0, e))?;
            let (endpoint, step) = candidate.into_parts();
            return Ok(CapillaryInertiaAdvanceCandidate {
                endpoint,
                step,
                report: CapillaryInertiaAdvanceReport {
                    attempts: 0,
                    accepted_intervals: 0,
                    rejected_intervals: 0,
                    step_evaluations: 1,
                    maximum_normalized_local_error: 0.0,
                },
            });
        }
        let pi = PiController {
            k_p: 0.7 / 3.0,
            k_i: 0.4 / 3.0,
            ..PiController::default()
        };
        let mut endpoint = *self;
        let mut elapsed = 0.0;
        let mut proposal = controls.initial_step_s;
        let mut previous_error = 1.0_f64;
        let mut previous_motion = None;
        let mut sums = [Sum::default(); 10];
        let mut attempts = 0usize;
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        let mut evaluations = 0usize;
        let mut maximum_error = 0.0_f64;
        let mut stops = 0usize;
        let mut first_stop = None;
        let mut segments = 0usize;
        let mut iterations = 0usize;
        let mut max_pressure_residual = 0.0_f64;
        let mut initial_angle = None;
        let mut last_step = None;
        while elapsed < input.duration_s {
            poll(&mut cancelled).map_err(|e| step_error(0.0, e))?;
            if attempts == controls.maximum_attempts {
                return Err(Error::Attempts {
                    attempts,
                    accepted_intervals: accepted,
                });
            }
            attempts += 1;
            let remaining = input.duration_s - elapsed;
            let chosen = proposal.min(remaining);
            let end = if chosen == remaining {
                input.duration_s
            } else {
                elapsed + chosen
            };
            let h = end - elapsed;
            let half = 0.5 * h;
            if end <= elapsed || !h.is_finite() || half <= 0.0 {
                return Err(Error::StepUnderflow);
            }
            evaluations += 1;
            let coarse = match endpoint.step_candidate(
                InertialCapillaryBoundary {
                    duration_s: h,
                    ..input
                },
                solve,
                &mut cancelled,
            ) {
                Ok(value) => value,
                Err(error) if refinable(&error) && h > controls.minimum_step_s => {
                    rejected += 1;
                    proposal = shrink_to(h, elapsed, controls.minimum_step_s, 0.5 * h)?;
                    continue;
                }
                Err(error) => return Err(step_error(h, error)),
            };
            // Both a coarse trial and its first half can otherwise resolve the
            // identical unresolved stop, then stay pinned: their endpoint
            // discrepancy would be zero despite an inaccurate stopping time.
            if let Some(time) = coarse.step().first_stopping_time_s {
                if time < half {
                    rejected += 1;
                    // Keep a post-stop remainder. Ending exactly at the
                    // coarse midpoint stop can leave the two half steps with
                    // a tiny positive velocity, then repeat forever as the
                    // stopping horizon collapses cubically. This horizon puts
                    // the first half before the stop and leaves the second
                    // half room to publish actual rest or reversal.
                    let horizon = time + 0.5 * time;
                    if horizon < controls.minimum_step_s {
                        return Err(Error::EventResolution {
                            event_time_s: time,
                            minimum_step_s: controls.minimum_step_s,
                        });
                    }
                    proposal = shrink_to(h, elapsed, controls.minimum_step_s, horizon)?;
                    previous_error = 1.0;
                    continue;
                }
            }
            let fine = (|| {
                evaluations += 1;
                let first = endpoint.step_candidate(
                    InertialCapillaryBoundary {
                        duration_s: half,
                        ..input
                    },
                    solve,
                    &mut cancelled,
                )?;
                evaluations += 1;
                let second = first.endpoint().step_candidate(
                    InertialCapillaryBoundary {
                        duration_s: h - half,
                        ..input
                    },
                    solve,
                    &mut cancelled,
                )?;
                Ok::<_, MaterialInertiaError>((first, second))
            })();
            let (first, second) = match fine {
                Ok(value) => value,
                Err(error) if refinable(&error) && h > controls.minimum_step_s => {
                    rejected += 1;
                    proposal = shrink_to(h, elapsed, controls.minimum_step_s, 0.5 * h)?;
                    continue;
                }
                Err(error) => return Err(step_error(h, error)),
            };
            let changed = coarse.step().stopping_events > 0
                || first.step().stopping_events > 0
                || second.step().stopping_events > 0
                || coarse.step().final_motion != first.step().final_motion
                || first.step().final_motion != second.step().final_motion
                || previous_motion.is_some_and(|m| m != first.step().final_motion);
            let divisor = if changed { 1.0 } else { 3.0 };
            let mut error = state_error(
                endpoint.state(),
                coarse.endpoint().state(),
                second.endpoint().state(),
                controls,
                divisor,
            );
            let fine_stop = first
                .step()
                .first_stopping_time_s
                .or_else(|| second.step().first_stopping_time_s.map(|time| half + time));
            if let (Some(coarse_time), Some(fine_time)) =
                (coarse.step().first_stopping_time_s, fine_stop)
            {
                let scale = controls.event_time_absolute_tolerance_s
                    + controls.relative_tolerance * coarse_time.max(fine_time);
                let event_error = (fine_time - coarse_time).abs() / scale;
                error = error.max(event_error);
            }
            if error.is_nan() {
                return Err(Error::Invalid("nonfinite local error estimate"));
            }
            if changed {
                previous_error = 1.0;
            }
            if error > 1.0 {
                rejected += 1;
                if h <= controls.minimum_step_s {
                    return Err(Error::Tolerance {
                        normalized_error: error,
                        step_s: h,
                    });
                }
                let factor = (pi.safety * error.powf(-1.0 / 3.0)).clamp(pi.max_shrink, 0.9);
                proposal = shrink_to(h, elapsed, controls.minimum_step_s, h * factor)?;
                continue;
            }
            poll(&mut cancelled).map_err(|e| step_error(h, e))?;
            initial_angle.get_or_insert(first.step().initial_apparent_angle_radians);
            for (offset, step) in [(elapsed, first.step()), (elapsed + half, second.step())] {
                accumulate_interval(&mut sums, step);
                stops = stops
                    .checked_add(step.stopping_events)
                    .ok_or(Error::Invalid("stopping-event count overflow"))?;
                segments = segments
                    .checked_add(step.segments)
                    .ok_or(Error::Invalid("segment count overflow"))?;
                iterations = iterations
                    .checked_add(step.iterations)
                    .ok_or(Error::Invalid("iteration count overflow"))?;
                if let Some(time) = step.first_stopping_time_s {
                    first_stop.get_or_insert(offset + time);
                }
                max_pressure_residual =
                    max_pressure_residual.max(step.maximum_pressure_residual_pa);
            }
            let (next, physical) = second.into_parts();
            previous_motion = Some(physical.final_motion);
            last_step = Some(physical);
            endpoint = next;
            elapsed = end;
            accepted += 1;
            maximum_error = maximum_error.max(error);
            let safe_error = error.max(1e-300);
            let factor = pi.safety * safe_error.powf(-pi.k_p) * previous_error.powf(pi.k_i);
            if !factor.is_finite() {
                return Err(Error::Invalid("nonfinite PI step proposal"));
            }
            proposal = (h * factor.clamp(pi.max_shrink, if changed { 1.0 } else { pi.max_growth }))
                .clamp(controls.minimum_step_s, controls.maximum_step_s);
            previous_error = if changed { 1.0 } else { safe_error };
        }
        let last = last_step.expect("positive interval has an accepted step");
        let step = self
            .finish_step(
                endpoint,
                input,
                sums,
                solve.mechanics,
                initial_angle.expect("positive interval has an initial apparent angle"),
                last.final_apparent_angle_radians,
                last.final_motion,
                stops,
                first_stop,
                segments,
                iterations,
                max_pressure_residual,
            )
            .map_err(|e| step_error(input.duration_s, e))?;
        poll(&mut cancelled).map_err(|e| step_error(input.duration_s, e))?;
        Ok(CapillaryInertiaAdvanceCandidate {
            endpoint,
            step,
            report: CapillaryInertiaAdvanceReport {
                attempts,
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                step_evaluations: evaluations,
                maximum_normalized_local_error: maximum_error,
            },
        })
    }
}

fn validate(c: CapillaryInertiaAdvanceControls) -> Result<(), CapillaryInertiaAdvanceError> {
    if ![
        c.initial_step_s,
        c.minimum_step_s,
        c.maximum_step_s,
        c.length_absolute_tolerance_m,
        c.velocity_absolute_tolerance_m_s,
        c.event_time_absolute_tolerance_s,
    ]
    .into_iter()
    .all(|x| x.is_finite() && x > 0.0)
        || c.minimum_step_s > c.initial_step_s
        || c.initial_step_s > c.maximum_step_s
        || !c.relative_tolerance.is_finite()
        || !(0.0..=1.0).contains(&c.relative_tolerance)
        || c.maximum_attempts == 0
        || c.maximum_attempts > usize::MAX / 3
    {
        return Err(CapillaryInertiaAdvanceError::Invalid(
            "positive ordered timestep/error scales and a bounded attempt budget required",
        ));
    }
    Ok(())
}

fn state_error(
    old: super::InertialMeniscusState,
    coarse: super::InertialMeniscusState,
    fine: super::InertialMeniscusState,
    c: CapillaryInertiaAdvanceControls,
    divisor: f64,
) -> f64 {
    let mut error = 0.0_f64;
    for (a, b, d, atol) in [
        (
            old.wet_length_m,
            coarse.wet_length_m,
            fine.wet_length_m,
            c.length_absolute_tolerance_m,
        ),
        (
            old.velocity_m_s,
            coarse.velocity_m_s,
            fine.velocity_m_s,
            c.velocity_absolute_tolerance_m_s,
        ),
    ] {
        let scale = atol + c.relative_tolerance * a.abs().max(b.abs()).max(d.abs());
        if !scale.is_finite() || scale <= 0.0 {
            return f64::INFINITY;
        }
        let value = ((d - b).abs() / divisor) / scale;
        if !value.is_finite() {
            return f64::INFINITY;
        }
        error = error.max(value);
    }
    error
}

fn refinable(error: &MaterialInertiaError) -> bool {
    matches!(
        error,
        MaterialInertiaError::Mechanics(
            mechanics::InertialCapillaryError::StepTooLarge(_)
                | mechanics::InertialCapillaryError::NotConverged { .. }
        ) | MaterialInertiaError::StoppingBudget { .. }
    )
}

fn shrink_to(
    old: f64,
    elapsed: f64,
    minimum: f64,
    target: f64,
) -> Result<f64, CapillaryInertiaAdvanceError> {
    let proposed = target.max(minimum);
    let represented = (elapsed + proposed) - elapsed;
    if !proposed.is_finite()
        || proposed <= 0.0
        || proposed >= old
        || represented <= 0.0
        || represented >= old
    {
        return Err(CapillaryInertiaAdvanceError::StepUnderflow);
    }
    Ok(proposed)
}

#[cfg(test)]
mod tests {
    use super::{CapillaryInertiaAdvanceError, shrink_to};
    #[test]
    fn rounded_minimum_cannot_repeat_an_identical_rejected_interval() {
        let elapsed = 1.0;
        let minimum = 0.1;
        let h = (elapsed + minimum) - elapsed;
        assert!(h > minimum);
        assert_eq!(
            shrink_to(h, elapsed, minimum, 0.5 * h),
            Err(CapillaryInertiaAdvanceError::StepUnderflow)
        );
    }
}
