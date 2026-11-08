//! Local time-error control for the original conservative SSPRK2 stages.
//!
//! Every comparison is exactly one full step versus two half steps. All three
//! use the same native CFL/domain admission, without internal subcycling. The
//! accepted state and every physical ledger come from the two fine steps.

use super::*;
use crate::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

/// Accuracy and bounded-work controls for explicit enthalpy transport.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveSsprk2Controls {
    /// Shared physical/CFL controls and stage-query budget for the whole call.
    pub step: Ssprk2Controls,
    pub initial_dt_s: f64,
    /// Floor for a complete coarse/fine comparison, not either fine half.
    /// An exact event/target remainder may be shorter; a rejected remainder
    /// cannot silently bypass this floor by further subdivision.
    pub minimum_dt_s: f64,
    /// Counts all attempted comparisons, including failed CFL/domain trials.
    pub maximum_trials: usize,
    /// Mandatory even on an isothermal latent plateau, in J/kg.
    pub enthalpy_absolute_tolerance_j_kg: f64,
    /// Optional additional norm using each cell's original T(h), in kelvin.
    pub temperature_absolute_tolerance_k: Option<f64>,
    /// Scales the change from the interval start, never an absolute h datum.
    pub relative_tolerance: f64,
}

/// Raw coarse/fine differences and the corrected fine local-error norm.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ssprk2TimeError {
    pub maximum_enthalpy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: Option<f64>,
    /// Maximum of corrected difference / (absolute + relative * change).
    /// The second-order correction is 1/3 for exactly equal fine spans.
    pub normalized: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct AdaptiveSsprk2Report {
    /// Accepted energy/flow work only; `energy.substeps` counts both fine steps.
    /// Stage queries and CFL/domain rejections include all discarded trials.
    pub work: Ssprk2Report,
    pub trials: usize,
    pub accepted_intervals: usize,
    pub accuracy_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<Ssprk2TimeError>,
    /// May be passed as the next call's initial dt; no hidden controller history.
    pub suggested_dt_s: f64,
}

/// Failure publishes neither a partial state nor speculative flow accounts.
#[derive(Clone, Debug, PartialEq)]
pub enum AdaptiveSsprk2Error {
    Deposition(DepositionError),
    Invalid(&'static str),
    MinimumStep {
        time_s: f64,
        dt_s: f64,
    },
    TimeAccuracy {
        time_s: f64,
        dt_s: f64,
        error: Ssprk2TimeError,
    },
}

impl From<DepositionError> for AdaptiveSsprk2Error {
    fn from(error: DepositionError) -> Self {
        Self::Deposition(error)
    }
}

impl Advected<()> {
    /// Error-controlled SSPRK2 with mandatory datum-invariant enthalpy error.
    /// Accepts the conservative fine trajectory; never extrapolates enthalpy.
    /// This is a smooth-solution local estimate, not a certified global bound.
    pub fn advance_ssprk2_adaptive_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: AdaptiveSsprk2Controls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, AdaptiveSsprk2Report), AdaptiveSsprk2Error> {
        self.advance_ssprk2_adaptive_with_coefficients(
            thermal, target_s, controls, unchanged, checkpoint,
        )
    }

    /// Repeated stage-local source/conductance queries have the same contract
    /// as [`Self::advance_ssprk2_with_coefficients`]. Native events and schedule
    /// knots cut comparisons exactly; speculative queries advance no history.
    pub fn advance_ssprk2_adaptive_with_coefficients<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: AdaptiveSsprk2Controls,
        coefficients: impl FnMut(
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, AdaptiveSsprk2Report), AdaptiveSsprk2Error> {
        advance_adaptive(self, thermal, target_s, controls, coefficients, checkpoint)
    }
}

impl ReconstructedAdvected<()> {
    /// Adaptive SSPRK2 for original limited fluxes and continuous pump/inlet
    /// schedules, with the same original-law h/T norms as donor-cell flow.
    pub fn advance_ssprk2_adaptive_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: AdaptiveSsprk2Controls,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, AdaptiveSsprk2Report), AdaptiveSsprk2Error> {
        self.advance_ssprk2_adaptive_with_coefficients(
            thermal, target_s, controls, unchanged, checkpoint,
        )
    }

    /// Repeatable stage-local coefficients; no Euler entropy or observer
    /// history is advanced at discarded coarse or fine stages.
    pub fn advance_ssprk2_adaptive_with_coefficients<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: AdaptiveSsprk2Controls,
        coefficients: impl FnMut(
            &DepositionSimulation<'m>,
        ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, Self, AdaptiveSsprk2Report), AdaptiveSsprk2Error> {
        advance_adaptive(self, thermal, target_s, controls, coefficients, checkpoint)
    }
}

fn time_error(
    initial: &DepositionSimulation<'_>,
    coarse: &DepositionSimulation<'_>,
    fine: &DepositionSimulation<'_>,
    correction: f64,
    controls: AdaptiveSsprk2Controls,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<Ssprk2TimeError, DepositionError> {
    let mut error = Ssprk2TimeError {
        maximum_temperature_difference_k: controls.temperature_absolute_tolerance_k.map(|_| 0.0),
        ..Ssprk2TimeError::default()
    };
    let component = |old: f64, a: f64, b: f64, absolute: f64| {
        let difference = finite(b - a)?.abs();
        let change = finite(a - old)?.abs().max(finite(b - old)?.abs());
        let allowance = finite(absolute + finite(controls.relative_tolerance * change)?)?;
        let ratio = (difference / allowance) * correction;
        // Positive infinity requests refinement; it is never fabricated into
        // a finite estimate. All physical coordinates remain finite.
        if ratio.is_nan() || ratio < 0.0 || allowance <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        Ok((difference, ratio))
    };
    for (index, (&id, old)) in initial.cells.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let (difference, ratio) = component(
            old.specific_h_j_kg,
            coarse.cells[&id].specific_h_j_kg,
            fine.cells[&id].specific_h_j_kg,
            controls.enthalpy_absolute_tolerance_j_kg,
        )?;
        error.maximum_enthalpy_difference_j_kg =
            error.maximum_enthalpy_difference_j_kg.max(difference);
        error.normalized = error.normalized.max(ratio);
        if let Some(absolute) = controls.temperature_absolute_tolerance_k {
            let (difference, ratio) = component(
                initial.temperature_k(id)?,
                coarse.temperature_k(id)?,
                fine.temperature_k(id)?,
                absolute,
            )?;
            error.maximum_temperature_difference_k = Some(
                error
                    .maximum_temperature_difference_k
                    .unwrap_or(0.0)
                    .max(difference),
            );
            error.normalized = error.normalized.max(ratio);
        }
    }
    Ok(error)
}

fn compare<'m, O: PureFlow>(
    initial: &DepositionSimulation<'m>,
    first: &mut Stage<O>,
    end: f64,
    unrepresented_heat_j: f64,
    controls: AdaptiveSsprk2Controls,
    work: &mut Ssprk2Report,
    coefficients: &mut impl FnMut(
        &DepositionSimulation<'m>,
    ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<Option<(Interval<'m, O>, Ssprk2TimeError)>, DepositionError> {
    let start = initial.time_s;
    let middle = start.midpoint(end);
    let (a, b) = (middle - start, end - middle);
    let correction = ConductionTimeOrder::Second
        .fine_error_correction(a, b)
        .ok_or(DepositionError::Arithmetic)?;
    let Some(coarse) = try_interval(
        initial,
        first,
        end,
        end - start,
        unrepresented_heat_j,
        controls.step,
        work,
        coefficients,
        checkpoint,
    )?
    else {
        return Ok(None);
    };
    let Some(mut half) = try_interval(
        initial,
        first,
        middle,
        a,
        unrepresented_heat_j,
        controls.step,
        work,
        coefficients,
        checkpoint,
    )?
    else {
        return Ok(None);
    };
    let mut midpoint = match evaluate(
        &mut half.thermal,
        &half.owner,
        controls.step,
        work,
        coefficients,
        checkpoint,
    ) {
        Ok(stage) => stage,
        Err(DepositionError::MaterialDomain) => return Ok(None),
        Err(error) => return Err(error),
    };
    if b > midpoint.rates.explicit_dt_s {
        return Ok(None);
    }
    let Some(mut fine) = try_interval(
        &half.thermal,
        &mut midpoint,
        end,
        b,
        half.unrepresented_heat_j,
        controls.step,
        work,
        coefficients,
        checkpoint,
    )?
    else {
        return Ok(None);
    };
    fine.external_heat.add_scaled(&half.external_heat, 1.0)?;
    let error = time_error(
        initial,
        &coarse.thermal,
        &fine.thermal,
        correction,
        controls,
        checkpoint,
    )?;
    Ok(Some((fine, error)))
}

fn advance_adaptive<'m, O: PureFlow>(
    owner: &O,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: AdaptiveSsprk2Controls,
    mut coefficients: impl FnMut(
        &DepositionSimulation<'m>,
    ) -> Result<BoundaryStepCoefficients<()>, DepositionError>,
    mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
) -> Result<(DepositionSimulation<'m>, O, AdaptiveSsprk2Report), AdaptiveSsprk2Error> {
    checkpoint()?;
    let positive = |x: f64| x.is_finite() && x > 0.0;
    if !target_s.is_finite()
        || target_s < thermal.time_s
        || !positive(controls.step.step.max_dt_s)
        || !positive(controls.initial_dt_s)
        || controls.initial_dt_s > controls.step.step.max_dt_s
        || !positive(controls.minimum_dt_s)
        || controls.minimum_dt_s > controls.initial_dt_s
        || !positive(controls.enthalpy_absolute_tolerance_j_kg)
        || controls
            .temperature_absolute_tolerance_k
            .is_some_and(|x| !positive(x))
        || !controls.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.relative_tolerance)
        || !controls.step.step.energy_tolerance_j.is_finite()
        || controls.step.step.energy_tolerance_j < 0.0
    {
        return Err(AdaptiveSsprk2Error::Invalid(
            "invalid adaptive SSPRK2 target, step or tolerance",
        ));
    }
    if owner.flow().entropy.is_some() {
        return Err(AdaptiveSsprk2Error::Invalid(
            "SSPRK2 cannot advance an Euler flow entropy history",
        ));
    }
    finite(thermal.energy_j())?;
    let mut next = thermal.clone();
    let mut owner = owner.clone();
    owner.begin()?;
    let mut report = AdaptiveSsprk2Report {
        work: Ssprk2Report {
            energy: empty_energy(),
            stage_evaluations: 0,
            rejected_intervals: 0,
        },
        trials: 0,
        accepted_intervals: 0,
        accuracy_rejections: 0,
        maximum_accepted_normalized_error: 0.0,
        last_error: None,
        suggested_dt_s: controls.initial_dt_s,
    };
    let mut external_heat = EnergySum::default();
    let mut unrepresented_heat = 0.0;
    loop {
        checkpoint()?;
        while next.cursor < next.schedule.len() && next.schedule[next.cursor].time_s == next.time_s
        {
            checkpoint()?;
            let event = next.schedule[next.cursor].event.clone();
            next.apply(event.clone(), &mut report.work.energy)?;
            next.cursor += 1;
            report.work.energy.events += 1;
            owner.event(&event, &next)?;
        }
        if next.time_s == target_s {
            break;
        }
        let mut stop = next
            .schedule
            .get(next.cursor)
            .map_or(target_s, |event| event.time_s.min(target_s));
        if let Some(transition) = owner.next_transition_time_s(&next, &mut checkpoint)? {
            if !transition.is_finite() || transition <= next.time_s {
                return Err(AdaptiveSsprk2Error::Invalid(
                    "invalid SSPRK2 source transition",
                ));
            }
            stop = stop.min(transition);
        }
        if next.cells.is_empty() {
            next.time_s = stop;
            continue;
        }
        if next
            .budget
            .max_substeps
            .saturating_sub(report.work.energy.substeps)
            < 2
        {
            return Err(DepositionError::Budget("SSPRK2 accepted intervals").into());
        }
        let mut first = evaluate(
            &mut next,
            &owner,
            controls.step,
            &mut report.work,
            &mut coefficients,
            &mut checkpoint,
        )?;
        let start = next.time_s;
        let bound = first.rates.explicit_dt_s;
        let mut proposed = report.suggested_dt_s.min(bound).min(stop - start);
        loop {
            checkpoint()?;
            if report.trials >= controls.maximum_trials {
                return Err(DepositionError::Budget("adaptive SSPRK2 trials").into());
            }
            let mut end = if proposed == stop - start {
                stop
            } else {
                finite(start + proposed)?.min(stop)
            };
            // Rounding the clock must not turn an admitted CFL span into a
            // longer one. All integrals and the estimator use actual spans.
            if end - start > bound {
                end = end.next_down();
            }
            let dt = end - start;
            if end <= start
                || start.midpoint(end) == start
                || start.midpoint(end) == end
                || (dt < controls.minimum_dt_s && end != stop)
            {
                return Err(AdaptiveSsprk2Error::MinimumStep {
                    time_s: start,
                    dt_s: dt,
                });
            }
            report.trials += 1;
            let trial = compare(
                &next,
                &mut first,
                end,
                unrepresented_heat,
                controls,
                &mut report.work,
                &mut coefficients,
                &mut checkpoint,
            )?;
            let factor = if let Some((fine, error)) = trial {
                report.last_error = Some(error);
                if error.normalized <= 1.0 {
                    external_heat.add_scaled(&fine.external_heat, 1.0)?;
                    next = fine.thermal;
                    owner = fine.owner;
                    unrepresented_heat = fine.unrepresented_heat_j;
                    report.work.energy.substeps += 2;
                    report.accepted_intervals += 1;
                    report.maximum_accepted_normalized_error = report
                        .maximum_accepted_normalized_error
                        .max(error.normalized);
                    let growth = ConductionTimeOrder::Second
                        .accepted_growth_factor(error.normalized)
                        .ok_or(DepositionError::Arithmetic)?;
                    report.suggested_dt_s =
                        (dt * growth).clamp(controls.minimum_dt_s, controls.step.step.max_dt_s);
                    break;
                }
                report.accuracy_rejections += 1;
                if dt <= controls.minimum_dt_s {
                    return Err(AdaptiveSsprk2Error::TimeAccuracy {
                        time_s: start,
                        dt_s: dt,
                        error,
                    });
                }
                ConductionTimeOrder::Second
                    .rejected_reduction_factor(error.normalized)
                    .ok_or(DepositionError::Arithmetic)?
            } else {
                let mut shorter = dt;
                halve(
                    &mut shorter,
                    &mut report.work,
                    controls.step.maximum_stage_evaluations,
                )?;
                if dt <= controls.minimum_dt_s {
                    return Err(AdaptiveSsprk2Error::MinimumStep {
                        time_s: start,
                        dt_s: dt,
                    });
                }
                0.5
            };
            proposed = (dt * factor).max(controls.minimum_dt_s);
            if proposed >= dt {
                return Err(AdaptiveSsprk2Error::MinimumStep {
                    time_s: start,
                    dt_s: dt,
                });
            }
        }
    }
    thermal.admit_advance_energy(
        &next,
        &external_heat,
        &mut report.work.energy,
        controls.step.step.energy_tolerance_j,
        &mut checkpoint,
    )?;
    Ok((next, owner, report))
}
