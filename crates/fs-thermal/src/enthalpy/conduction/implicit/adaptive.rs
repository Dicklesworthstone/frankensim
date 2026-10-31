//! Enthalpy-aware temporal accuracy around the original backward-Euler owner.
//!
//! One full interval is compared with two half intervals. Only the original
//! conservative half steps are accepted; comparison heat and entropy are never
//! published and no Richardson-extrapolated physical state is constructed.
//! Change-based error scales retain arbitrary caloric reference invariance and
//! resolve latent inventory even when every compared temperature is identical.
//! This follows the existing fs-conduction adaptive enthalpy semantics without
//! replacing the shared thermal residual or introducing a new material solver.
//! The inherited h-labelled coordinate measures the specific thermal energy
//! selected by the original law's storage basis. Enthalpy and fixed-volume
//! internal energy remain distinct declared potentials throughout accounting.

use super::*;

/// Local temporal estimates, separate from algebraic joule acceptance gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveImplicitConductionControls {
    /// Original endpoint controls and whole-call work ceilings. Minimum/maximum
    /// dt bound full comparison spans; accepted physical half steps are shorter.
    /// `maximum_halvings` bounds consecutive reductions at one accepted time.
    pub step: ImplicitConductionControls,
    /// At least the step floor; capped at the maximum and the next exact event.
    pub initial_dt_s: f64,
    /// Whole-call comparisons, including failed or rejected work. Each attempts
    /// at most three original endpoint solves; their residual budget is shared.
    pub maximum_trials: usize,
    /// Absolute error allowance on the original specific stored-energy
    /// coordinate, J/kg. Its inherited h label does not change the law's basis.
    pub enthalpy_absolute_tolerance_j_kg: f64,
    pub temperature_absolute_tolerance_k: f64,
    /// Some requires an original phase observable for every active cell.
    pub phase_fraction_absolute_tolerance: Option<f64>,
    /// Applied to each coordinate's CHANGE from the accepted state, never its
    /// absolute Kelvin or arbitrary enthalpy datum. Must be in [0,1).
    pub relative_tolerance: f64,
}
impl Default for AdaptiveImplicitConductionControls {
    fn default() -> Self {
        Self {
            step: ImplicitConductionControls::default(),
            initial_dt_s: 1.0,
            maximum_trials: 65_536,
            enthalpy_absolute_tolerance_j_kg: 1.0,
            temperature_absolute_tolerance_k: 1e-3,
            phase_fraction_absolute_tolerance: None,
            relative_tolerance: 1e-3,
        }
    }
}

/// Actual coarse/fine discrepancies. `normalized` includes the first-order
/// correction for unequal representable half spans; it is not a global bound.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ImplicitConductionTimeError {
    pub maximum_enthalpy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: f64,
    pub maximum_phase_fraction_difference: Option<f64>,
    pub normalized: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveImplicitConductionReport {
    /// Energy, entropy, events and accepted substeps describe only the fine
    /// trajectory. Nonlinear iterations and residual evaluations include ALL
    /// work. Rejected intervals include accuracy and nonlinear rejections.
    pub solve: ImplicitConductionReport,
    pub comparisons: usize,
    pub accepted_comparisons: usize,
    pub accuracy_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: Option<ImplicitConductionTimeError>,
    /// Explicit controller state for a later call; there is no hidden PI memory.
    pub suggested_dt_s: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AdaptiveImplicitConductionError {
    Implicit(ImplicitConductionError),
    Invalid(&'static str),
    TrialBudget,
    TimeAccuracy {
        time_s: f64,
        span_s: f64,
        error: ImplicitConductionTimeError,
    },
}
impl From<ImplicitConductionError> for AdaptiveImplicitConductionError {
    fn from(error: ImplicitConductionError) -> Self {
        Self::Implicit(error)
    }
}
impl From<DepositionError> for AdaptiveImplicitConductionError {
    fn from(error: DepositionError) -> Self {
        Self::Implicit(error.into())
    }
}
impl core::fmt::Display for AdaptiveImplicitConductionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive implicit conduction: {self:?}")
    }
}
impl std::error::Error for AdaptiveImplicitConductionError {}

/// Endpoint adapter for the shared conservative step-doubling controller.
///
/// An implementation retains the original fixed cell identities, masses and
/// caloric laws. It supplies complete backward-Euler endpoints, their original
/// energy/entropy accounts and source-owned observables, never an alternative
/// temperature evolution. Every residual evaluation and Newton iteration,
/// including failed trials, must be charged to `work` by `solve_interval`.
/// Events are applied by the original timeline through the supplied observer;
/// `admit` must refuse changes outside the adapter's fixed physical binding.
pub trait AdaptiveConductionStepper: Sized {
    type Workspace;
    type Observer: DepositionObserver<Error = DepositionError>;

    fn cell_ids(&self) -> &[CellId];
    fn workspace_cells(workspace: &Self::Workspace) -> usize;
    fn observer(&self) -> Self::Observer;
    fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError>;

    /// Solve exactly this interval, without internal physical subdivisions.
    /// The residual ceiling in `controls` is shared through cumulative `work`.
    fn solve_interval<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut Self::Workspace,
        work: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, EnergyAudit, EntropyAudit, f64), ImplicitConductionError>;

    /// Sensible owners need only the original caloric temperature. A phase
    /// owner must explicitly supply its own phase observable when requested.
    fn observable(
        &self,
        thermal: &DepositionSimulation<'_>,
        id: CellId,
    ) -> Result<ImplicitConductionObservables, DepositionError> {
        Ok(ImplicitConductionObservables {
            temperature_k: thermal.temperature_k(id)?,
            phase_mass_fraction: None,
        })
    }

    /// Advance with local h/T/phase accuracy estimates and original heat/entropy.
    ///
    /// Exact scheduled boundaries cap each comparison before either trajectory
    /// is computed. All candidate work shares one residual-evaluation ceiling.
    /// Only Newton nonconvergence, a singular Newton system or excessive time
    /// error causes a retry. Material, energy, cancellation and work refusals
    /// propagate unchanged. No endpoint is silently internally subdivided.
    ///
    /// The fine pair consumes two accepted-step slots and publishes atomically.
    /// A late failure returns no partial prefix. A final/event-clipped span below
    /// the floor can be tried once; it is never accepted merely for being small.
    /// The local estimator needs independent temporal and spatial validation,
    /// particularly at constitutive corners. It is not a certified global bound.
    fn advance_adaptive_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: AdaptiveImplicitConductionControls,
        workspace: &mut Self::Workspace,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (DepositionSimulation<'m>, AdaptiveImplicitConductionReport),
        AdaptiveImplicitConductionError,
    > {
        checkpoint()?;
        validate_advance(
            thermal.time_s,
            target_s,
            self.cell_ids().len(),
            controls.step,
            Self::workspace_cells(workspace),
        )?;
        if !controls.initial_dt_s.is_finite()
            || controls.initial_dt_s < controls.step.minimum_dt_s
            || controls.maximum_trials == 0
            || !positive(controls.enthalpy_absolute_tolerance_j_kg)
            || !positive(controls.temperature_absolute_tolerance_k)
            || !controls.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_tolerance)
            || controls
                .phase_fraction_absolute_tolerance
                .is_some_and(|x| !positive(x) || x > 1.0)
        {
            return Err(AdaptiveImplicitConductionError::Invalid(
                "invalid adaptive accuracy, initial step or comparison budget",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let mut next = thermal.clone();
        let mut observer = self.observer();
        let mut old_observables = self.time_observables(&next, controls, &mut checkpoint)?;
        let mut span = controls.initial_dt_s.min(controls.step.maximum_dt_s);
        let mut report = AdaptiveImplicitConductionReport {
            solve: empty_report(),
            comparisons: 0,
            accepted_comparisons: 0,
            accuracy_rejections: 0,
            maximum_accepted_normalized_error: 0.0,
            last_error: None,
            suggested_dt_s: span,
        };
        let mut supplied = Sum::default();
        let mut reductions = 0;
        loop {
            checkpoint()?;
            let (with_events, observed, events) = next.advance_observed(
                &observer,
                next.time_s,
                StepConfig {
                    max_dt_s: controls.step.maximum_dt_s,
                    energy_tolerance_j: controls.step.total_tolerance_j,
                },
                &mut checkpoint,
            )?;
            next = with_events;
            observer = observed;
            report.solve.energy.events += events.events;
            self.admit(&next, &mut checkpoint)?;
            if next.time_s == target_s {
                break;
            }
            let maximum = controls.step.maximum_steps.min(next.budget.max_substeps);
            if maximum.saturating_sub(report.solve.energy.substeps) < 2 {
                return Err(DepositionError::Budget("adaptive implicit half steps").into());
            }
            if report.comparisons >= controls.maximum_trials {
                return Err(AdaptiveImplicitConductionError::TrialBudget);
            }
            let start = next.time_s;
            let stop = next
                .schedule
                .get(next.cursor)
                .map_or(target_s, |event| event.time_s.min(target_s));
            let end = if span >= stop - start {
                stop
            } else {
                finite(start + span)?
            };
            let middle = f64::midpoint(start, end);
            if !(start < middle && middle < end && end <= stop) {
                return Err(AdaptiveImplicitConductionError::Invalid(
                    "adaptive interval has no representable interior midpoint",
                ));
            }
            let dt = end - start;
            report.comparisons += 1;
            let candidate = (|| {
                let (coarse, _, _, _) = self.solve_interval(
                    &next,
                    end,
                    controls.step,
                    workspace,
                    &mut report.solve,
                    &mut checkpoint,
                )?;
                let (half, e1, s1, r1) = self.solve_interval(
                    &next,
                    middle,
                    controls.step,
                    workspace,
                    &mut report.solve,
                    &mut checkpoint,
                )?;
                let (fine, e2, s2, r2) = self.solve_interval(
                    &half,
                    end,
                    controls.step,
                    workspace,
                    &mut report.solve,
                    &mut checkpoint,
                )?;
                Ok::<_, ImplicitConductionError>((coarse, fine, e1, e2, s1, s2, r1.max(r2)))
            })();
            let (cause, factor) = match candidate {
                Ok((coarse, fine, e1, e2, s1, s2, residual)) => {
                    let coarse_observables =
                        self.time_observables(&coarse, controls, &mut checkpoint)?;
                    let fine_observables =
                        self.time_observables(&fine, controls, &mut checkpoint)?;
                    // First-order fine error / (coarse-fine) =
                    // (h1²+h2²)/(2*h1*h2), exactly one for equal halves.
                    let (h1, h2) = (middle - start, end - middle);
                    let correction = finite(0.5 * (h1 / h2 + h2 / h1))?;
                    let error = self.time_error(
                        &next,
                        &coarse,
                        &fine,
                        [&old_observables, &coarse_observables, &fine_observables],
                        controls,
                        correction,
                        &mut checkpoint,
                    )?;
                    if error.normalized <= 1.0 {
                        supplied.add(e1.external_heat_j)?;
                        supplied.add(e2.external_heat_j)?;
                        accumulate_entropy(&mut report.solve.entropy, s1)?;
                        accumulate_entropy(&mut report.solve.entropy, s2)?;
                        report.solve.maximum_cell_residual_j =
                            report.solve.maximum_cell_residual_j.max(residual);
                        report.solve.energy.substeps += 2;
                        report.accepted_comparisons += 1;
                        report.maximum_accepted_normalized_error = report
                            .maximum_accepted_normalized_error
                            .max(error.normalized);
                        report.last_error = Some(error);
                        let growth = if error.normalized == 0.0 {
                            2.0
                        } else {
                            (0.9 / error.normalized.sqrt()).clamp(0.5, 2.0)
                        };
                        let suggested = if growth >= 1.0 {
                            dt.min(controls.step.maximum_dt_s / growth) * growth
                        } else {
                            dt * growth
                        };
                        span =
                            suggested.clamp(controls.step.minimum_dt_s, controls.step.maximum_dt_s);
                        report.suggested_dt_s = span;
                        old_observables = fine_observables;
                        next = fine;
                        reductions = 0;
                        continue;
                    }
                    report.accuracy_rejections += 1;
                    (
                        AdaptiveImplicitConductionError::TimeAccuracy {
                            time_s: start,
                            span_s: dt,
                            error,
                        },
                        (0.9 / error.normalized.sqrt()).clamp(0.1, 0.5),
                    )
                }
                Err(
                    error @ (ImplicitConductionError::NotConverged { .. }
                    | ImplicitConductionError::KrylovNotConverged { .. }
                    | ImplicitConductionError::Linear(TridiagonalError::Singular { .. })),
                ) => (error.into(), 0.5),
                Err(error) => return Err(error.into()),
            };
            report.solve.rejected_intervals += 1;
            if reductions >= controls.step.maximum_halvings || dt <= controls.step.minimum_dt_s {
                return Err(cause);
            }
            let wanted = (dt * factor).max(controls.step.minimum_dt_s);
            let mut smaller = start + wanted;
            if smaller - start < controls.step.minimum_dt_s {
                smaller = smaller.next_up();
            }
            if !(start < smaller && smaller < end) {
                return Err(cause);
            }
            span = smaller - start;
            reductions += 1;
        }
        let mut storage = Sum::default();
        for &id in self.cell_ids() {
            checkpoint()?;
            storage.add(
                thermal.cells[&id].mass_kg
                    * (next.cells[&id].specific_h_j_kg - thermal.cells[&id].specific_h_j_kg),
            )?;
        }
        report.solve.energy.external_heat_j = supplied.total()?;
        report.solve.energy.stored_energy_change_j = storage.total()?;
        report.solve.energy.residual_j = finite(
            report.solve.energy.stored_energy_change_j - report.solve.energy.external_heat_j,
        )?;
        if report.solve.energy.residual_j.abs() > controls.step.total_tolerance_j {
            return Err(DepositionError::EnergyBalance(report.solve.energy.residual_j).into());
        }
        checkpoint()?;
        Ok((next, report))
    }

    fn time_observables(
        &self,
        thermal: &DepositionSimulation<'_>,
        controls: AdaptiveImplicitConductionControls,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<ImplicitConductionObservables>, AdaptiveImplicitConductionError> {
        let mut values = Vec::with_capacity(self.cell_ids().len());
        for (i, &id) in self.cell_ids().iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let value = self.observable(thermal, id)?;
            if !positive(value.temperature_k)
                || value
                    .phase_mass_fraction
                    .is_some_and(|x| !x.is_finite() || !(0.0..=1.0).contains(&x))
            {
                return Err(DepositionError::MaterialDomain.into());
            }
            if controls.phase_fraction_absolute_tolerance.is_some()
                && value.phase_mass_fraction.is_none()
            {
                return Err(AdaptiveImplicitConductionError::Invalid(
                    "phase accuracy requested without an original material phase observable",
                ));
            }
            values.push(value);
        }
        checkpoint()?;
        Ok(values)
    }

    fn time_error(
        &self,
        old: &DepositionSimulation<'_>,
        coarse: &DepositionSimulation<'_>,
        fine: &DepositionSimulation<'_>,
        values: [&[ImplicitConductionObservables]; 3],
        controls: AdaptiveImplicitConductionControls,
        correction: f64,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitConductionTimeError, AdaptiveImplicitConductionError> {
        let mut error = ImplicitConductionTimeError {
            maximum_phase_fraction_difference: controls
                .phase_fraction_absolute_tolerance
                .map(|_| 0.0),
            ..ImplicitConductionTimeError::default()
        };
        for (i, &id) in self.cell_ids().iter().enumerate() {
            if i % 256 == 0 {
                checkpoint()?;
            }
            let (difference, ratio) = component_error(
                old.cells[&id].specific_h_j_kg,
                coarse.cells[&id].specific_h_j_kg,
                fine.cells[&id].specific_h_j_kg,
                controls.enthalpy_absolute_tolerance_j_kg,
                controls.relative_tolerance,
                correction,
            )?;
            error.maximum_enthalpy_difference_j_kg =
                error.maximum_enthalpy_difference_j_kg.max(difference);
            error.normalized = error.normalized.max(ratio);
            let (difference, ratio) = component_error(
                values[0][i].temperature_k,
                values[1][i].temperature_k,
                values[2][i].temperature_k,
                controls.temperature_absolute_tolerance_k,
                controls.relative_tolerance,
                correction,
            )?;
            error.maximum_temperature_difference_k =
                error.maximum_temperature_difference_k.max(difference);
            error.normalized = error.normalized.max(ratio);
            if let Some(tolerance) = controls.phase_fraction_absolute_tolerance {
                let phase = |which: usize| {
                    values[which][i]
                        .phase_mass_fraction
                        .expect("phase admission precedes error comparison")
                };
                let (difference, ratio) = component_error(
                    phase(0),
                    phase(1),
                    phase(2),
                    tolerance,
                    controls.relative_tolerance,
                    correction,
                )?;
                error.maximum_phase_fraction_difference = Some(
                    error
                        .maximum_phase_fraction_difference
                        .unwrap_or(0.0)
                        .max(difference),
                );
                error.normalized = error.normalized.max(ratio);
            }
        }
        checkpoint()?;
        Ok(error)
    }
}

impl<'a> AdaptiveConductionStepper for ImplicitConductionChain<'a> {
    type Workspace = ImplicitConductionWorkspace;
    type Observer = ConductionObserver<'a>;

    fn cell_ids(&self) -> &[CellId] {
        &self.ids
    }
    fn workspace_cells(workspace: &Self::Workspace) -> usize {
        workspace.cells()
    }
    fn observer(&self) -> Self::Observer {
        self.conduction.clone()
    }

    fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        ImplicitConductionChain::admit(self, thermal, checkpoint)
    }

    fn observable(
        &self,
        thermal: &DepositionSimulation<'_>,
        id: CellId,
    ) -> Result<ImplicitConductionObservables, DepositionError> {
        let cell = &thermal.cells[&id];
        self.laws[&cell.material_id].time_observables(cell.specific_h_j_kg)
    }

    fn solve_interval<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: ImplicitConductionControls,
        workspace: &mut Self::Workspace,
        work: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, EnergyAudit, EntropyAudit, f64), ImplicitConductionError>
    {
        self.solve_step(thermal, target_s, controls, workspace, work, checkpoint)
    }
}

impl ImplicitConductionChain<'_> {
    /// Use the shared original-law step-doubling controller. See
    /// [`AdaptiveConductionStepper::advance_adaptive_to`] for error/event semantics.
    pub fn advance_adaptive_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: AdaptiveImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (DepositionSimulation<'m>, AdaptiveImplicitConductionReport),
        AdaptiveImplicitConductionError,
    > {
        AdaptiveConductionStepper::advance_adaptive_to(
            self, thermal, target_s, controls, workspace, checkpoint,
        )
    }
}

fn component_error(
    old: f64,
    coarse: f64,
    fine: f64,
    absolute: f64,
    relative: f64,
    correction: f64,
) -> Result<(f64, f64), DepositionError> {
    let difference = finite(fine - coarse)?.abs();
    let change = finite(coarse - old)?.abs().max(finite(fine - old)?.abs());
    let allowance = finite(absolute + finite(relative * change)?)?;
    let ratio = (difference / allowance) * correction;
    // An infinite ratio causes rejection; it must not become a fabricated
    // finite error. Every physical coordinate and difference remains finite.
    if ratio.is_nan() || ratio < 0.0 || !positive(allowance) {
        return Err(DepositionError::Arithmetic);
    }
    Ok((difference, ratio))
}

#[cfg(test)]
#[path = "adaptive/tests.rs"]
mod tests;
