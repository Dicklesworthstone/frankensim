//! Local physical-time accuracy on the original simultaneous wall/film/gas solve.
//!
//! One coarse interval is compared with two conservative fine intervals. Only
//! the fine pair enters the shared window account and retained history. All
//! proposals spend one numerical work budget, including rejected comparisons.
//! Native wall enthalpy and phase fraction detect latent storage when T is flat.
//! This is a local error estimate, not a bound on the full trajectory's error.

use super::super::implicit::{WallFilmImplicitError, WallFilmImplicitReport, WorkMeter};
use super::*;
use fs_material::phase::EquilibriumPhaseState;
use fs_thermal::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

mod error;
pub use error::WallFilmTimeError;
use error::{Endpoint, compare};

/// Physical time accuracy plus the existing whole-window conservation/work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallFilmAdaptiveControls {
    pub integration: WallFilmImplicitWindowControls,
    /// Initial full comparison span. The window's min/max steps bound this span;
    /// the two accepted physical half intervals may be shorter than its minimum.
    pub initial_step_s: f64,
    pub maximum_comparisons: usize,
    /// Consecutive accuracy/admission reductions at one accepted physical time.
    pub maximum_reductions: usize,
    pub enthalpy_absolute_tolerance_j_kg: f64,
    /// Applies separately to the actual wall and gas temperatures.
    pub temperature_absolute_tolerance_k: f64,
    pub phase_fraction_absolute_tolerance: f64,
    /// Applies separately to liquid inventory and every original gas species.
    pub mass_absolute_tolerance_kg: f64,
    /// Scales CHANGE from the interval start, never absolute formation h or T.
    pub relative_tolerance: f64,
}
impl WallFilmAdaptiveControls {
    pub fn new(integration: WallFilmImplicitWindowControls) -> Self {
        Self {
            initial_step_s: integration.window.maximum_step_s,
            integration,
            maximum_comparisons: 65_536,
            maximum_reductions: 32,
            enthalpy_absolute_tolerance_j_kg: 1.0,
            temperature_absolute_tolerance_k: 1e-3,
            phase_fraction_absolute_tolerance: 1e-4,
            mass_absolute_tolerance_kg: 1e-10,
            relative_tolerance: 1e-3,
        }
    }
    fn validate(self) -> Result<(), WallFilmAdaptiveError> {
        if !positive(self.initial_step_s)
            || self.initial_step_s < self.integration.window.minimum_step_s
            || self.maximum_comparisons == 0
            || !positive(self.enthalpy_absolute_tolerance_j_kg)
            || !positive(self.temperature_absolute_tolerance_k)
            || !positive(self.phase_fraction_absolute_tolerance)
            || self.phase_fraction_absolute_tolerance > 1.0
            || !positive(self.mass_absolute_tolerance_kg)
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
        {
            return Err(WallFilmAdaptiveError::Invalid(
                "invalid wall-film time accuracy controls",
            ));
        }
        Ok(())
    }
}

/// Accepted fine trajectory and work spent on all coarse/fine proposals.
#[derive(Clone, Debug, PartialEq)]
pub struct WallFilmAdaptiveReport {
    /// Attempted/rejected physical intervals include discarded coarse solves
    /// and all unsuccessful fine proposals; samples contain accepted fine steps.
    pub window: WallFilmImplicitWindowReport,
    pub comparisons: usize,
    pub accepted_comparisons: usize,
    pub accuracy_rejections: usize,
    pub admission_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<WallFilmTimeError>,
    /// Explicit initial-step suggestion for a subsequent source window.
    pub suggested_step_s: f64,
}

#[derive(Debug)]
pub enum WallFilmAdaptiveError {
    Window(WallFilmWindowError),
    Invalid(&'static str),
    ComparisonBudget {
        accepted_time_s: f64,
        comparisons: usize,
    },
    TimeAccuracy {
        accepted_time_s: f64,
        span_s: f64,
        error: WallFilmTimeError,
    },
}
impl core::fmt::Display for WallFilmAdaptiveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive finite wall-film exchange: {self:?}")
    }
}
impl std::error::Error for WallFilmAdaptiveError {}
impl From<WallFilmWindowError> for WallFilmAdaptiveError {
    fn from(error: WallFilmWindowError) -> Self {
        Self::Window(error)
    }
}
impl From<WallFilmImplicitError> for WallFilmAdaptiveError {
    fn from(error: WallFilmImplicitError) -> Self {
        Self::Window(error.into())
    }
}

struct Proposal {
    candidate: UniformEnthalpyStep<WallFilmImplicitReport>,
    internal_heat_j: f64,
    gas_heat_j: f64,
}
impl Proposal {
    fn endpoint(&self) -> Endpoint<'_> {
        Endpoint {
            wall: self.candidate.state,
            gas: &self.candidate.report.next_gas,
            film: self.candidate.report.next_film,
        }
    }
}

impl WallFilmThermalEnvironment<'_> {
    /// Advance one uniform-source window with local h/T/phase/inventory control.
    ///
    /// Body and gas heat remain TOTAL window energies. Cumulative endpoints
    /// apportion them once, including the last partial comparison. Rejected
    /// work contributes no heat, mass, film event or elapsed physical time.
    /// Only the original explicit phase StepSizeLimit may trigger an admission
    /// retry; other source, nonlinear, conservation and cancellation errors
    /// propagate. No extrapolated phase state or artificial wetting seed is used.
    ///
    /// The first-order controller and unequal-half correction are the existing
    /// thermal owner's arithmetic. Its estimate does not certify a continuous
    /// dry-out event or accumulated global error. Known changes of source or
    /// environment require a new authored window at that physical boundary.
    pub fn advance_implicit_adaptive(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
        controls: WallFilmAdaptiveControls,
    ) -> Result<UniformEnthalpyStep<WallFilmAdaptiveReport>, WallFilmAdaptiveError> {
        poll(cx)?;
        controls.validate()?;
        let window = controls.integration.window;
        let mut account = WindowAccumulator::new(cx, self, input, window)?;
        let meter = WorkMeter::new(
            controls.integration.solve,
            controls.integration.maximum_gas_iterations,
        )?;
        let order = ConductionTimeOrder::First;
        let mut span = controls.initial_step_s.min(window.maximum_step_s);
        let mut attempts = 0usize;
        let mut comparisons = 0usize;
        let mut accepted_comparisons = 0usize;
        let mut accuracy_rejections = 0usize;
        let mut admission_rejections = 0usize;
        let mut reductions = 0usize;
        let mut maximum_accepted_normalized_error = 0.0_f64;
        let mut last_error = None;
        while account.time_s < input.duration_s {
            poll(cx)?;
            if comparisons >= controls.maximum_comparisons {
                return Err(WallFilmAdaptiveError::ComparisonBudget {
                    accepted_time_s: account.time_s,
                    comparisons,
                });
            }
            let start = account.time_s;
            let end = if span >= input.duration_s - start {
                input.duration_s
            } else {
                finite(start + span)?
            };
            let middle = f64::midpoint(start, end);
            if !(start < middle && middle < end && end <= input.duration_s) {
                return Err(WallFilmAdaptiveError::Invalid(
                    "adaptive window has no representable interior midpoint",
                ));
            }
            let dt = end - start;
            comparisons = count(comparisons, 1)?;
            let initial = Endpoint {
                wall: account.state,
                gas: &account.gas,
                film: account.film,
            };
            let proposal = (|| {
                let coarse = interval(
                    self,
                    cx,
                    input,
                    initial,
                    start,
                    end,
                    start,
                    controls.integration,
                    &meter,
                    &mut attempts,
                )?;
                let first = interval(
                    self,
                    cx,
                    input,
                    initial,
                    start,
                    middle,
                    start,
                    controls.integration,
                    &meter,
                    &mut attempts,
                )?;
                let fine = interval(
                    self,
                    cx,
                    input,
                    first.endpoint(),
                    middle,
                    end,
                    start,
                    controls.integration,
                    &meter,
                    &mut attempts,
                )?;
                Ok::<_, WallFilmWindowError>((coarse, first, fine))
            })();
            let (cause, wanted) = match proposal {
                Ok((coarse, first, fine)) => {
                    let correction = order
                        .fine_error_correction(middle - start, end - middle)
                        .ok_or(WallFilmAdaptiveError::Invalid(
                            "unrepresentable time error correction",
                        ))?;
                    let error = compare(
                        self,
                        cx,
                        initial,
                        coarse.endpoint(),
                        fine.endpoint(),
                        controls,
                        correction,
                    )?;
                    if error.normalized <= 1.0 {
                        account.record_step(
                            middle,
                            first.internal_heat_j,
                            first.gas_heat_j,
                            WindowStep::implicit(first.candidate),
                        )?;
                        account.record_step(
                            end,
                            fine.internal_heat_j,
                            fine.gas_heat_j,
                            WindowStep::implicit(fine.candidate),
                        )?;
                        accepted_comparisons = count(accepted_comparisons, 1)?;
                        maximum_accepted_normalized_error =
                            maximum_accepted_normalized_error.max(error.normalized);
                        last_error = Some(error);
                        let growth = order.accepted_growth_factor(error.normalized).ok_or(
                            WallFilmAdaptiveError::Invalid("invalid accepted time error"),
                        )?;
                        let suggested = if growth >= 1.0 {
                            dt.min(window.maximum_step_s / growth) * growth
                        } else {
                            dt * growth
                        };
                        span = suggested.clamp(window.minimum_step_s, window.maximum_step_s);
                        reductions = 0;
                        continue;
                    }
                    accuracy_rejections = count(accuracy_rejections, 1)?;
                    let factor = order.rejected_reduction_factor(error.normalized).ok_or(
                        WallFilmAdaptiveError::Invalid("invalid rejected time error"),
                    )?;
                    (
                        WallFilmAdaptiveError::TimeAccuracy {
                            accepted_time_s: start,
                            span_s: dt,
                            error,
                        },
                        dt * factor,
                    )
                }
                Err(
                    error @ WallFilmWindowError::StepSizeLimit {
                        maximum_duration_s, ..
                    },
                ) => {
                    admission_rejections = count(admission_rejections, 1)?;
                    let wanted = if positive(maximum_duration_s) {
                        maximum_duration_s.min(0.5 * dt)
                    } else {
                        0.0
                    };
                    (WallFilmAdaptiveError::Window(error), wanted)
                }
                Err(error) => return Err(error.into()),
            };
            if reductions >= controls.maximum_reductions
                || dt <= window.minimum_step_s
                || !positive(wanted)
            {
                return Err(cause);
            }
            // Accuracy can request the floor. An explicit source admission may
            // forbid even that interval, so never enlarge its stated maximum.
            let accuracy = matches!(cause, WallFilmAdaptiveError::TimeAccuracy { .. });
            let wanted = if accuracy {
                wanted.max(window.minimum_step_s)
            } else {
                if wanted < window.minimum_step_s {
                    return Err(cause);
                }
                wanted
            };
            let mut smaller = start + wanted;
            if !accuracy && smaller - start > wanted {
                smaller = smaller.next_down();
            }
            if accuracy && smaller - start < window.minimum_step_s {
                smaller = smaller.next_up();
            }
            if !(start < smaller && smaller < end)
                || smaller - start < window.minimum_step_s
                || (!accuracy && smaller - start > wanted)
            {
                return Err(cause);
            }
            span = smaller - start;
            reductions = count(reductions, 1)?;
        }
        let absolute_mass = account.sum_absolute_mass_residual_kg();
        let accepted_intervals =
            accepted_comparisons
                .checked_mul(2)
                .ok_or(WallFilmAdaptiveError::Invalid(
                    "adaptive accepted interval count overflows",
                ))?;
        let rejected =
            attempts
                .checked_sub(accepted_intervals)
                .ok_or(WallFilmAdaptiveError::Invalid(
                    "adaptive work count precedes accepted intervals",
                ))?;
        let physical = account.finish(cx, attempts, rejected)?;
        poll(cx)?;
        Ok(UniformEnthalpyStep {
            state: physical.state,
            external_heat_j: physical.external_heat_j,
            energy_residual_tolerance_j: physical.energy_residual_tolerance_j,
            report: WallFilmAdaptiveReport {
                window: WallFilmImplicitWindowReport {
                    physical: physical.report,
                    newton_iterations: meter.newton_iterations(),
                    residual_evaluations: meter.residual_evaluations(),
                    gas_iterations: meter.gas_iterations(),
                    charged_gas_iterations: meter.charged_gas_iterations(),
                    sum_absolute_mass_residual_kg: absolute_mass,
                },
                comparisons,
                accepted_comparisons,
                accuracy_rejections,
                admission_rejections,
                maximum_accepted_normalized_error,
                last_error,
                suggested_step_s: span,
            },
        })
    }
}

#[allow(clippy::too_many_arguments)]
fn interval(
    environment: &WallFilmThermalEnvironment<'_>,
    cx: &Cx<'_>,
    input: UniformEnthalpyStepInput<'_>,
    initial: Endpoint<'_>,
    start: f64,
    end: f64,
    accepted_time_s: f64,
    controls: WallFilmImplicitWindowControls,
    meter: &WorkMeter,
    attempts: &mut usize,
) -> Result<Proposal, WallFilmWindowError> {
    poll(cx)?;
    if *attempts >= controls.window.maximum_attempts {
        return Err(WallFilmWindowError::WorkLimit {
            attempted_intervals: *attempts,
            accepted_time_s,
        });
    }
    let internal_heat_j = increment(input.internal_heat_j, start, end, input.duration_s)?;
    let gas_heat_j = increment(environment.gas_heat_j, start, end, input.duration_s)?;
    let stage = WallFilmThermalEnvironment {
        gas: initial.gas,
        film: initial.film,
        gas_heat_j,
        ..environment.clone()
    };
    *attempts = count(*attempts, 1)?;
    let candidate = stage
        .advance_implicit_metered(
            cx,
            UniformEnthalpyStepInput {
                initial: initial.wall,
                internal_heat_j,
                duration_s: end - start,
                ..input
            },
            controls.solve,
            meter,
        )
        .map_err(|error| match error {
            WallFilmImplicitError::Physical(WallFilmThermalError::Phase(
                GasEvaporationError::StepSizeLimit {
                    maximum_duration_s,
                    reason,
                },
            )) => WallFilmWindowError::StepSizeLimit {
                accepted_time_s,
                requested_duration_s: end - start,
                maximum_duration_s,
                reason,
            },
            other => other.into(),
        })?;
    Ok(Proposal {
        candidate,
        internal_heat_j,
        gas_heat_j,
    })
}
