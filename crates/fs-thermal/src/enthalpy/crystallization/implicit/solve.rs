//! Existing Newton--Krylov on the original kinetic thermal-network residual.
//!
//! Every residual trial starts from one frozen physical history. The existing
//! fallible finite-difference owner supplies numerical Jacobian actions; the
//! original thermal and phase owners still decide physical acceptance.

use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

use fs_math::dd::Dd;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};

use super::super::{
    CrystallizationAdaptiveAudit, CrystallizationAudit, CrystallizationControls,
    CrystallizationSimulation, CrystallizationTimelineError, DepositionError, empty_audit, finite,
    positive,
};
use super::CrystallizationBackwardEulerStep;
use crate::enthalpy::backward_euler::{BoundaryStepCoefficients, StepCoefficients};

/// Solver controls are separate from the original caloric and energy gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CrystallizationImplicitControls {
    pub phase: CrystallizationControls,
    pub newton: NewtonKrylovConfig,
    /// Maximum-component perturbation in CELL joules for numerical Jacobian
    /// actions. It must resolve the original stored coordinate and caloric solve.
    pub difference_step_j: f64,
    /// Whole-call limits, including every physical endpoint reconstruction.
    pub maximum_iterations: usize,
    pub maximum_evaluations: usize,
}
impl Default for CrystallizationImplicitControls {
    fn default() -> Self {
        Self {
            phase: CrystallizationControls::default(),
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-9,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            difference_step_j: 1e-2,
            maximum_iterations: 10_000,
            maximum_evaluations: 100_000,
        }
    }
}
impl CrystallizationImplicitControls {
    pub fn validate(self) -> Result<(), CrystallizationImplicitError> {
        self.newton.validate()?;
        self.phase
            .caloric
            .validate()
            .map_err(CrystallizationTimelineError::from)?;
        if !positive(self.phase.maximum_dt_s)
            || !positive(self.difference_step_j)
            || !self.phase.cell_energy_tolerance_j.is_finite()
            || self.phase.cell_energy_tolerance_j < 0.0
            || !self.phase.total_energy_tolerance_j.is_finite()
            || self.phase.total_energy_tolerance_j < 0.0
            || self.maximum_iterations == 0
            || self.maximum_evaluations == 0
        {
            return Err(CrystallizationImplicitError::Invalid(
                "invalid implicit crystallization controls",
            ));
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum CrystallizationImplicitError {
    Physical(CrystallizationTimelineError),
    Numerical(NewtonError),
    Invalid(&'static str),
    NonConverged(NewtonReport),
    Interrupted(SolverRunProgress),
}
impl From<CrystallizationTimelineError> for CrystallizationImplicitError {
    fn from(value: CrystallizationTimelineError) -> Self {
        Self::Physical(value)
    }
}
impl From<DepositionError> for CrystallizationImplicitError {
    fn from(value: DepositionError) -> Self {
        Self::Physical(value.into())
    }
}
impl From<NewtonError> for CrystallizationImplicitError {
    fn from(value: NewtonError) -> Self {
        Self::Numerical(value)
    }
}
impl core::fmt::Display for CrystallizationImplicitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit crystallization: {self:?}")
    }
}
impl std::error::Error for CrystallizationImplicitError {}

/// Physical accounts include accepted intervals only. Numerical work includes
/// the initial guess, every uncached trial/stencil and final physical recheck.
#[derive(Clone, Debug, PartialEq)]
pub struct CrystallizationImplicitReport {
    pub phase: CrystallizationAudit,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    pub maximum_accepted_residual_norm_j: f64,
}

impl<'a> CrystallizationSimulation<'a> {
    /// Solve one complete interval over all original active cells and links.
    /// Ordinary finite bodies retain their own caloric law and storage basis;
    /// kinetic cells require their explicitly selected adaptive history controls.
    ///
    /// The supplied coefficient callback sees actual proposed kinetic states,
    /// not all-amorphous inverse temperatures. It may supply original endpoint
    /// sources, link conductances and occupied boundary-channel replacements.
    /// Its physical meaning must remain immutable during speculative calls.
    /// No contact topology, coverage, mass or source history may change here.
    ///
    /// Start events must already be processed. Events at target remain pending
    /// until the caller accepts all coupled owners; an interior event refuses.
    /// Newton convergence never replaces per-cell caloric/energy acceptance.
    /// Source-domain and stencil errors propagate rather than being converted
    /// into zero residuals or silently retried on a different material model.
    #[allow(clippy::type_complexity)]
    pub fn solve_implicit_step<T>(
        &self,
        target_s: f64,
        controls: CrystallizationImplicitControls,
        mut resolve: impl FnMut(
            &CrystallizationSimulation<'a>,
        )
            -> Result<BoundaryStepCoefficients<T>, CrystallizationTimelineError>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, T, CrystallizationImplicitReport), CrystallizationImplicitError> {
        checkpoint()?;
        controls.validate()?;
        let step =
            CrystallizationBackwardEulerStep::new(self, target_s, controls.phase, &mut checkpoint)?;
        let poll = RefCell::new(&mut checkpoint);
        let source = RefCell::new(&mut resolve);
        let evaluations = Cell::new(0usize);
        // Numerical Jv repeatedly requests the same base point. Two bounded
        // entries retain it alongside the last stencil and avoid reconstructing
        // every material history twice per direction. Only residuals are cached;
        // the accepted physical state and source payload are recomputed below.
        let cache: RefCell<Vec<(Vec<f64>, Vec<f64>)>> = RefCell::new(Vec::with_capacity(2));
        let problem = FallibleFiniteDifference::new(
            step.cell_ids().len(),
            controls.difference_step_j,
            |x: &[f64]| {
                (poll.borrow_mut())()?;
                let coordinates = step.resolvable_increments_j(x);
                {
                    let mut entries = cache.borrow_mut();
                    if let Some(index) = entries
                        .iter()
                        .position(|(point, _)| same_point(point, &coordinates))
                    {
                        let entry = entries.remove(index);
                        let residual = entry.1.clone();
                        entries.push(entry);
                        return Ok(residual);
                    }
                }
                spend_evaluation(&evaluations, controls.maximum_evaluations)?;
                let trial = step.evaluate(
                    &coordinates,
                    |state| (source.borrow_mut())(state),
                    || (poll.borrow_mut())(),
                )?;
                let residual = trial.residual_j().to_vec();
                let mut entries = cache.borrow_mut();
                if entries.len() == 2 {
                    entries.remove(0);
                }
                entries.push((coordinates, residual.clone()));
                Ok::<_, CrystallizationTimelineError>(residual)
            },
        )
        .map_err(difference_error)?;
        let initialization =
            NewtonKrylovState::new(&problem, vec![0.0; step.cell_ids().len()], controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        let mut solver = initialization?;
        let mut checkpoint_error = None;
        let run = solver.run_controlled(&problem, controls.maximum_iterations, &mut || {
            if problem.has_failed() {
                return true;
            }
            match (poll.borrow_mut())() {
                Ok(()) => false,
                Err(error) => {
                    checkpoint_error = Some(error);
                    true
                }
            }
        });
        if let Some(error) = problem.take_failure() {
            return Err(difference_error(error));
        }
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        if run.progress != SolverRunProgress::Complete {
            return Err(CrystallizationImplicitError::Interrupted(run.progress));
        }
        if !run.report.converged {
            return Err(CrystallizationImplicitError::NonConverged(run.report));
        }
        spend_evaluation(&evaluations, controls.maximum_evaluations)?;
        let trial = step.evaluate(
            &step.resolvable_increments_j(&solver.x),
            |state| (source.borrow_mut())(state),
            || (poll.borrow_mut())(),
        )?;
        let (next, payload, phase) = trial.accept(|| (poll.borrow_mut())())?;
        (poll.borrow_mut())()?;
        Ok((
            next,
            payload,
            CrystallizationImplicitReport {
                phase,
                newton_iterations: run.report.iterations,
                residual_evaluations: evaluations.get(),
                maximum_accepted_residual_norm_j: run.report.residual_norm,
            },
        ))
    }

    /// Advance the existing deposition timeline with implicit kinetic heat
    /// transport and its authored sources/conductances. No explicit thermal CFL
    /// is imposed; maximum_dt_s remains an independent physical time-resolution
    /// requirement, not a nonlinear accuracy certificate.
    ///
    /// Birth/contact/separation/boundary events are consumed by the original
    /// timeline in authored order, including events at the requested endpoint.
    /// Each interval freezes its own accepted histories. Work and noncancelling
    /// energy budgets are shared across the entire call. A late failure returns
    /// no partial simulation, including events or newly deposited inventory.
    pub fn advance_implicit_to(
        &self,
        target_s: f64,
        controls: CrystallizationImplicitControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, CrystallizationImplicitReport), CrystallizationImplicitError> {
        checkpoint()?;
        controls.validate()?;
        if !target_s.is_finite() || target_s < self.time_s() {
            return Err(CrystallizationImplicitError::Invalid(
                "implicit target must be finite and not earlier",
            ));
        }
        let mut next = self.clone();
        let mut report = CrystallizationImplicitReport {
            phase: CrystallizationAudit {
                thermal: empty_audit(),
                crystallized_mass_kg: 0.0,
                maximum_cell_energy_residual_j: 0.0,
                absolute_energy_residual_j: 0.0,
                adaptive_phase: CrystallizationAdaptiveAudit::default(),
            },
            newton_iterations: 0,
            residual_evaluations: 0,
            maximum_accepted_residual_norm_j: 0.0,
        };
        let mut sums = [Dd::ZERO; 7];
        loop {
            checkpoint()?;
            if next
                .pending_events()
                .first()
                .is_some_and(|event| event.time_s == next.time_s())
            {
                let (after_events, events) =
                    next.advance_to(next.time_s(), controls.phase, &mut checkpoint)?;
                add_audit(&mut report.phase, events, &mut sums)?;
                next = after_events;
            }
            if next.time_s() == target_s {
                break;
            }
            let stop = next
                .pending_events()
                .first()
                .map_or(target_s, |event| event.time_s.min(target_s));
            if next.cells().is_empty() {
                // The original timeline owns empty-time traversal and births.
                let (after_empty, events) =
                    next.advance_to(stop, controls.phase, &mut checkpoint)?;
                add_audit(&mut report.phase, events, &mut sums)?;
                next = after_empty;
                continue;
            }
            if report.phase.thermal.substeps >= self.thermal.budget.max_substeps {
                return Err(
                    DepositionError::Budget("implicit crystallization thermal intervals").into(),
                );
            }
            let remaining = finite(stop - next.time_s())?;
            let endpoint = if controls.phase.maximum_dt_s >= remaining {
                stop
            } else {
                finite(next.time_s() + controls.phase.maximum_dt_s)?
            };
            if endpoint <= next.time_s() {
                return Err(DepositionError::Arithmetic.into());
            }
            let mut local = controls;
            local.maximum_iterations = controls
                .maximum_iterations
                .saturating_sub(report.newton_iterations);
            local.maximum_evaluations = controls
                .maximum_evaluations
                .saturating_sub(report.residual_evaluations);
            if local.maximum_iterations == 0 || local.maximum_evaluations == 0 {
                return Err(DepositionError::Budget(
                    "implicit crystallization whole-call solver work",
                )
                .into());
            }
            let (candidate, (), accepted) = next.solve_implicit_step(
                endpoint,
                local,
                |_| {
                    Ok(StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: BTreeMap::new(),
                        payload: (),
                    }
                    .into())
                },
                &mut checkpoint,
            )?;
            report.newton_iterations =
                add_count(report.newton_iterations, accepted.newton_iterations)?;
            report.residual_evaluations =
                add_count(report.residual_evaluations, accepted.residual_evaluations)?;
            report.maximum_accepted_residual_norm_j = report
                .maximum_accepted_residual_norm_j
                .max(accepted.maximum_accepted_residual_norm_j);
            add_audit(&mut report.phase, accepted.phase, &mut sums)?;
            if report.phase.absolute_energy_residual_j > controls.phase.total_energy_tolerance_j {
                return Err(DepositionError::EnergyBalance(
                    report.phase.absolute_energy_residual_j,
                )
                .into());
            }
            next = candidate;
        }
        if report.phase.absolute_energy_residual_j > controls.phase.total_energy_tolerance_j {
            return Err(
                DepositionError::EnergyBalance(report.phase.absolute_energy_residual_j).into(),
            );
        }
        checkpoint()?;
        Ok((next, report))
    }
}

fn spend_evaluation(count: &Cell<usize>, limit: usize) -> Result<(), CrystallizationTimelineError> {
    if count.get() >= limit {
        return Err(
            DepositionError::Budget("implicit crystallization residual evaluations").into(),
        );
    }
    count.set(count.get() + 1);
    Ok(())
}
fn same_point(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits())
}
fn difference_error(
    error: FiniteDifferenceError<CrystallizationTimelineError>,
) -> CrystallizationImplicitError {
    match error {
        FiniteDifferenceError::Source(error) => error.into(),
        FiniteDifferenceError::Invalid(message) => CrystallizationImplicitError::Invalid(message),
    }
}
fn add_count(a: usize, b: usize) -> Result<usize, DepositionError> {
    a.checked_add(b).ok_or(DepositionError::Budget(
        "implicit crystallization work counter",
    ))
}
fn add_audit(
    total: &mut CrystallizationAudit,
    part: CrystallizationAudit,
    sums: &mut [Dd; 7],
) -> Result<(), DepositionError> {
    for (sum, value) in sums.iter_mut().zip([
        part.thermal.incoming_mass_kg,
        part.thermal.incoming_energy_j,
        part.thermal.external_heat_j,
        part.thermal.stored_energy_change_j,
        part.thermal.residual_j,
        part.crystallized_mass_kg,
        part.absolute_energy_residual_j,
    ]) {
        *sum = *sum + Dd::from_f64(value);
        finite(sum.hi)?;
        finite(sum.lo)?;
    }
    total.thermal.incoming_mass_kg = sums[0].to_f64();
    total.thermal.incoming_energy_j = sums[1].to_f64();
    total.thermal.external_heat_j = sums[2].to_f64();
    total.thermal.stored_energy_change_j = sums[3].to_f64();
    total.thermal.residual_j = sums[4].to_f64();
    total.crystallized_mass_kg = sums[5].to_f64();
    total.absolute_energy_residual_j = sums[6].to_f64();
    total.thermal.substeps = add_count(total.thermal.substeps, part.thermal.substeps)?;
    total.thermal.events = add_count(total.thermal.events, part.thermal.events)?;
    total.maximum_cell_energy_residual_j = total
        .maximum_cell_energy_residual_j
        .max(part.maximum_cell_energy_residual_j);
    let (a, b) = (&mut total.adaptive_phase, part.adaptive_phase);
    a.accepted_substeps = add_count(a.accepted_substeps, b.accepted_substeps)?;
    a.rejected_trials = add_count(a.rejected_trials, b.rejected_trials)?;
    a.attempted_trials = add_count(a.attempted_trials, b.attempted_trials)?;
    a.maximum_temperature_discrepancy_k = a
        .maximum_temperature_discrepancy_k
        .max(b.maximum_temperature_discrepancy_k);
    a.maximum_crystalline_mass_fraction_discrepancy = a
        .maximum_crystalline_mass_fraction_discrepancy
        .max(b.maximum_crystalline_mass_fraction_discrepancy);
    a.maximum_exposure_discrepancy = a
        .maximum_exposure_discrepancy
        .max(b.maximum_exposure_discrepancy);
    Ok(())
}
