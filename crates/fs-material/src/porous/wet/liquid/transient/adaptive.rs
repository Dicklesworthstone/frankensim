//! Local temporal error control around the original conservative capillary solve.
//!
//! Backward Euler step doubling estimates the error of the fine path. Only its
//! two original physical endpoints and face transfers are accepted; Richardson
//! extrapolation never moves a saturation outside the original material domain.
//! This is a first-order evolution with a local error estimate, not a global
//! accuracy certificate. Spatial error and constitutive validity are separate.

use super::{
    LiquidSaturationBoundaryStepInput, LiquidSaturationCandidate, LiquidSaturationColumn,
    LiquidSaturationImplicitControls, LiquidSaturationImplicitError,
    LiquidSaturationImplicitReport, LiquidSaturationStepInput, compensated_sum, invalid, positive,
};
use crate::MaterialError;

/// Local temporal accuracy, whole-window water conservation and bounded work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidSaturationAdvanceControls {
    /// Coarse comparison interval bounds [s]. A shorter final remainder is
    /// allowed once; an unsuccessful remainder cannot be shortened below it.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// Total original implicit solves, including rejected coarse/fine probes.
    /// The admitted range is 3..=100_000.
    pub maximum_evaluations: usize,
    pub absolute_saturation_tolerance: f64,
    /// Applied to the larger coarse/fine CHANGE from the current state, not
    /// background saturation. Zero selects the absolute tolerance alone.
    pub relative_saturation_tolerance: f64,
    /// Fixed original-to-final liquid volume residual per transverse area [m].
    /// This budget is never multiplied by the number of accepted steps.
    pub volume_balance_tolerance_m: f64,
}

/// Accepted water transfers and actual computational work over a complete window.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidSaturationAdvanceReport {
    pub interval_s: [f64; 2],
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    pub linear_solves: usize,
    pub state_evaluations: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    /// Largest accepted local saturation-error ratio; at most one.
    pub maximum_error_ratio: f64,
    /// Signed integrated original Darcy flux at each face, positive left to
    /// right [m]. Rejected probes never contribute. There are cells+1 entries.
    pub face_volume_per_area_m: Vec<f64>,
    pub liquid_volume_per_area_before_m: f64,
    pub liquid_volume_per_area_after_m: f64,
    pub boundary_volume_per_area_m: f64,
    pub balance_residual_m: f64,
}

/// No partial state is returned by a refused or cancelled window.
#[derive(Debug)]
pub enum LiquidSaturationAdvanceError {
    Material(MaterialError),
    Step {
        elapsed_s: f64,
        step_s: f64,
        error: LiquidSaturationImplicitError,
    },
    Accuracy {
        elapsed_s: f64,
        step_s: f64,
        error_ratio: f64,
    },
    WorkBudget {
        evaluations: usize,
    },
    Conservation {
        residual_m: f64,
        maximum_m: f64,
    },
    Cancelled,
}
impl core::fmt::Display for LiquidSaturationAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive liquid saturation: {self:?}")
    }
}
impl std::error::Error for LiquidSaturationAdvanceError {}
impl From<MaterialError> for LiquidSaturationAdvanceError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct LiquidSaturationAdvanceCandidate {
    state: LiquidSaturationColumn,
    report: LiquidSaturationAdvanceReport,
}
impl LiquidSaturationAdvanceCandidate {
    pub const fn state(&self) -> &LiquidSaturationColumn {
        &self.state
    }
    pub const fn report(&self) -> &LiquidSaturationAdvanceReport {
        &self.report
    }
    pub fn into_parts(self) -> (LiquidSaturationColumn, LiquidSaturationAdvanceReport) {
        (self.state, self.report)
    }
}

impl LiquidSaturationAdvanceControls {
    fn validate(self) -> Result<(), MaterialError> {
        if !positive(self.minimum_step_s)
            || !positive(self.maximum_step_s)
            || self.minimum_step_s > self.maximum_step_s
            || !(3..=100_000).contains(&self.maximum_evaluations)
            || !positive(self.absolute_saturation_tolerance)
            || self.absolute_saturation_tolerance >= 1.0
            || !self.relative_saturation_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_saturation_tolerance)
            || !positive(self.volume_balance_tolerance_m)
        {
            return Err(invalid(
                "adaptive liquid diffusion needs ordered positive steps, 3..100000 evaluations, positive absolute saturation/volume tolerances and relative tolerance in [0,1)",
            ));
        }
        Ok(())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), LiquidSaturationAdvanceError> {
    if cancelled() {
        Err(LiquidSaturationAdvanceError::Cancelled)
    } else {
        Ok(())
    }
}

fn work(report: &mut LiquidSaturationAdvanceReport, native: &LiquidSaturationImplicitReport) {
    report.linear_solves += native.linear_solves;
    report.state_evaluations += native.state_evaluations;
}

fn evaluate(
    original: &LiquidSaturationColumn,
    input: LiquidSaturationBoundaryStepInput,
    nonlinear: LiquidSaturationImplicitControls,
    controls: LiquidSaturationAdvanceControls,
    report: &mut LiquidSaturationAdvanceReport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<LiquidSaturationCandidate, LiquidSaturationAdvanceError> {
    poll(cancelled)?;
    if report.evaluations == controls.maximum_evaluations {
        return Err(LiquidSaturationAdvanceError::WorkBudget {
            evaluations: report.evaluations,
        });
    }
    report.evaluations += 1;
    match original.step_candidate_implicit_with_boundaries(input, nonlinear, cancelled) {
        Ok((candidate, native)) => {
            work(report, &native);
            Ok(candidate)
        }
        Err(LiquidSaturationImplicitError::Cancelled) => {
            Err(LiquidSaturationAdvanceError::Cancelled)
        }
        Err(error) => {
            match &error {
                LiquidSaturationImplicitError::NotConverged { report: native, .. }
                | LiquidSaturationImplicitError::NoAdmissibleTrial { report: native, .. } => {
                    work(report, native);
                }
                _ => {}
            }
            Err(LiquidSaturationAdvanceError::Step {
                elapsed_s: original.elapsed_s,
                step_s: input.duration_s,
                error,
            })
        }
    }
}

impl LiquidSaturationColumn {
    /// Advance a complete pressure-driven interval with estimated temporal
    /// accuracy. Every probe uses the existing source-qualified implicit owner.
    pub fn advance_candidate_implicit(
        &self,
        input: LiquidSaturationStepInput,
        nonlinear: LiquidSaturationImplicitControls,
        controls: LiquidSaturationAdvanceControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<LiquidSaturationAdvanceCandidate, LiquidSaturationAdvanceError> {
        self.advance_candidate_implicit_with_boundaries(
            input.into(),
            nonlinear,
            controls,
            cancelled,
        )
    }

    /// Advance with the original sealed/flux/pressure/Robin boundaries held
    /// fixed. Supported nonlinear failures may shorten a trial; invalid original
    /// source coordinates remain fatal. Accepted half-step fluxes close one
    /// independently recomputed original-to-final water-volume budget.
    pub fn advance_candidate_implicit_with_boundaries(
        &self,
        input: LiquidSaturationBoundaryStepInput,
        nonlinear: LiquidSaturationImplicitControls,
        controls: LiquidSaturationAdvanceControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LiquidSaturationAdvanceCandidate, LiquidSaturationAdvanceError> {
        poll(&mut cancelled)?;
        input.validate()?;
        nonlinear.validate()?;
        controls.validate()?;
        let end = self.elapsed_s + input.duration_s;
        if !end.is_finite() || end <= self.elapsed_s {
            return Err(invalid("adaptive liquid time must advance representably").into());
        }
        let before = self.liquid_volume_per_area_m()?;
        let mut report = LiquidSaturationAdvanceReport {
            interval_s: [self.elapsed_s, end],
            accepted_intervals: 0,
            rejected_intervals: 0,
            evaluations: 0,
            linear_solves: 0,
            state_evaluations: 0,
            minimum_accepted_substep_s: f64::INFINITY,
            maximum_accepted_substep_s: 0.0,
            maximum_error_ratio: 0.0,
            face_volume_per_area_m: vec![0.0; self.cells.len() + 1],
            liquid_volume_per_area_before_m: before,
            liquid_volume_per_area_after_m: before,
            boundary_volume_per_area_m: 0.0,
            balance_residual_m: 0.0,
        };
        let mut corrections = vec![0.0; self.cells.len() + 1];
        let mut boundary_correction = 0.0;
        let mut state = self.clone();
        let mut step = controls.maximum_step_s.min(input.duration_s);
        while state.elapsed_s < end {
            poll(&mut cancelled)?;
            let remaining = end - state.elapsed_s;
            let next = if step >= remaining {
                end
            } else {
                state.elapsed_s + step
            };
            let duration = next - state.elapsed_s;
            let middle = state.elapsed_s + 0.5 * duration;
            let first = middle - state.elapsed_s;
            let second = next - middle;
            if !positive(first) || !positive(second) || next > end {
                return Err(
                    invalid("adaptive liquid half steps must advance representably").into(),
                );
            }
            let trial = (|| {
                let coarse = evaluate(
                    &state,
                    LiquidSaturationBoundaryStepInput {
                        duration_s: duration,
                        ..input
                    },
                    nonlinear,
                    controls,
                    &mut report,
                    &mut cancelled,
                )?;
                let half = evaluate(
                    &state,
                    LiquidSaturationBoundaryStepInput {
                        duration_s: first,
                        ..input
                    },
                    nonlinear,
                    controls,
                    &mut report,
                    &mut cancelled,
                )?;
                let fine = evaluate(
                    half.state(),
                    LiquidSaturationBoundaryStepInput {
                        duration_s: second,
                        ..input
                    },
                    nonlinear,
                    controls,
                    &mut report,
                    &mut cancelled,
                )?;
                // Actual representable halves can be unequal on a large clock.
                // The fine/coarse first-order error ratio is (a²+b²)/(2ab).
                let a = first / duration;
                let b = second / duration;
                let factor = (a * a + b * b) / (2.0 * a * b);
                let mut ratio = 0.0_f64;
                for (i, ((&base, &rough), &smooth)) in state
                    .saturations
                    .iter()
                    .zip(&coarse.state().saturations)
                    .zip(&fine.state().saturations)
                    .enumerate()
                {
                    if i % 128 == 0 {
                        poll(&mut cancelled)?;
                    }
                    let scale = controls.absolute_saturation_tolerance
                        + controls.relative_saturation_tolerance
                            * (rough - base).abs().max((smooth - base).abs());
                    ratio = ratio.max(factor * (smooth - rough).abs() / scale);
                }
                if !ratio.is_finite() || ratio > 1.0 {
                    return Err(LiquidSaturationAdvanceError::Accuracy {
                        elapsed_s: state.elapsed_s,
                        step_s: duration,
                        error_ratio: ratio,
                    });
                }
                Ok((half, fine, ratio))
            })();
            match trial {
                Ok((half, fine, ratio)) => {
                    for candidate in [&half, &fine] {
                        // Integrate the net original boundary transfer directly.
                        // Subtracting accumulated gross inlet/outlet throughflow
                        // would lose a small real uptake over a long window.
                        let y = candidate.step.boundary_volume_per_area_m - boundary_correction;
                        let sum = report.boundary_volume_per_area_m + y;
                        boundary_correction = (sum - report.boundary_volume_per_area_m) - y;
                        if !sum.is_finite() {
                            return Err(invalid(
                                "integrated liquid boundary transfer is unrepresentable",
                            )
                            .into());
                        }
                        report.boundary_volume_per_area_m = sum;
                        for (i, &flux) in candidate.step.face_flux_m_s.iter().enumerate() {
                            if i % 128 == 0 {
                                poll(&mut cancelled)?;
                            }
                            let term = candidate.step.duration_s * flux;
                            let y = term - corrections[i];
                            let sum = report.face_volume_per_area_m[i] + y;
                            corrections[i] = (sum - report.face_volume_per_area_m[i]) - y;
                            if !sum.is_finite() {
                                return Err(invalid(
                                    "integrated liquid face transfer is unrepresentable",
                                )
                                .into());
                            }
                            report.face_volume_per_area_m[i] = sum;
                        }
                    }
                    state = fine.state;
                    report.accepted_intervals += 1;
                    report.minimum_accepted_substep_s =
                        report.minimum_accepted_substep_s.min(first.min(second));
                    report.maximum_accepted_substep_s =
                        report.maximum_accepted_substep_s.max(first.max(second));
                    report.maximum_error_ratio = report.maximum_error_ratio.max(ratio);
                    // Dyadic growth is deterministic and avoids noisy growth
                    // around the acceptance threshold or a vanishing estimate.
                    step = if ratio < 0.25 {
                        (2.0 * duration).min(controls.maximum_step_s)
                    } else {
                        duration
                    };
                }
                Err(error) => {
                    let recoverable = matches!(
                        error,
                        LiquidSaturationAdvanceError::Accuracy { .. }
                            | LiquidSaturationAdvanceError::Step {
                                error: LiquidSaturationImplicitError::NotConverged { .. }
                                    | LiquidSaturationImplicitError::NoAdmissibleTrial { .. },
                                ..
                            }
                    );
                    let minimum = controls.minimum_step_s.min(remaining);
                    let reduced = (0.5 * duration).max(minimum);
                    if !recoverable || reduced >= duration {
                        return Err(error);
                    }
                    report.rejected_intervals += 1;
                    step = reduced;
                }
            }
        }
        report.liquid_volume_per_area_after_m = state.liquid_volume_per_area_m()?;
        let change =
            compensated_sum(self.cells.iter().enumerate().map(|(i, cell)| {
                cell.capacity_m() * (state.saturations[i] - self.saturations[i])
            }))?;
        report.balance_residual_m = change - report.boundary_volume_per_area_m;
        if !report.balance_residual_m.is_finite()
            || report.balance_residual_m.abs() > controls.volume_balance_tolerance_m
        {
            return Err(LiquidSaturationAdvanceError::Conservation {
                residual_m: report.balance_residual_m,
                maximum_m: controls.volume_balance_tolerance_m,
            });
        }
        poll(&mut cancelled)?;
        Ok(LiquidSaturationAdvanceCandidate { state, report })
    }
}
