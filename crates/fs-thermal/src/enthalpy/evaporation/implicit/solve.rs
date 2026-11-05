//! Automatic Newton--Krylov on the original coupled heat/evaporation residual.
//!
//! Numerical differences and globalization belong to `fs_solver`. The original
//! thermal and liquid owners retain source domains, inventory and independent
//! physical acceptance. Every speculative evaluation starts from the same
//! initial state; no heat-only intermediate state or substitute phase law is
//! constructed. Nonlinear convergence is separate from timestep accuracy.

use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};

use super::{
    EvaporatingBackwardEulerAudit, EvaporatingCellBinding, EvaporatingNetworkBackwardEulerAudit,
    EvaporatingNetworkBackwardEulerStep,
};
use crate::enthalpy::backward_euler::BoundaryExchangeStepCoefficients;
use crate::enthalpy::evaporation::ThermalEvaporationError;
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, nonnegative, positive};
use crate::evaporation::{EvaporationLaw, GasFilm, LiquidSurfaceState};

/// Common nonlinear work and physical acceptance controls for the heat network.
/// Every liquid domain belongs to its own `EvaporatingCellBinding`; there is no
/// global range that could silently override heterogeneous source support.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvaporationNetworkImplicitControls {
    pub cell_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    pub mass_tolerance_kg: f64,
    pub newton: NewtonKrylovConfig,
    pub difference_step_j: f64,
    pub maximum_iterations: usize,
    pub maximum_evaluations: usize,
}
impl Default for EvaporationNetworkImplicitControls {
    fn default() -> Self {
        Self {
            cell_energy_tolerance_j: 1e-8,
            total_energy_tolerance_j: 1e-6,
            mass_tolerance_kg: 1e-12,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-9,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            difference_step_j: 1e-2,
            maximum_iterations: 64,
            maximum_evaluations: 10_000,
        }
    }
}
impl EvaporationNetworkImplicitControls {
    pub(super) fn validate(self) -> Result<(), EvaporationImplicitError> {
        self.newton.validate()?;
        if !nonnegative(self.cell_energy_tolerance_j)
            || !nonnegative(self.total_energy_tolerance_j)
            || !nonnegative(self.mass_tolerance_kg)
            || !positive(self.difference_step_j)
            || self.maximum_iterations == 0
            || self.maximum_evaluations == 0
        {
            return Err(EvaporationImplicitError::Invalid(
                "invalid implicit evaporation controls",
            ));
        }
        Ok(())
    }
}

/// Numerical work controls and independent physical acceptance budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvaporationImplicitControls {
    /// Explicitly admitted pure-liquid branch on the original caloric datum.
    pub pure_liquid_enthalpy_range_j_kg: [f64; 2],
    /// Largest accepted local energy defect, including the independent surface account [J].
    pub cell_energy_tolerance_j: f64,
    /// Largest accepted native/global open-system energy defect [J].
    pub total_energy_tolerance_j: f64,
    /// Largest accepted whole-network and represented-rate mass defect [kg].
    pub mass_tolerance_kg: f64,
    /// Original solver's linear-work, nonlinear-convergence and globalization policy.
    pub newton: NewtonKrylovConfig,
    /// Maximum-component numerical Jacobian perturbation in CELL joules.
    /// It must resolve the stored enthalpy coordinate and remain in source support.
    pub difference_step_j: f64,
    /// Maximum nonlinear attempts in this call, including rejected attempts.
    pub maximum_iterations: usize,
    /// Maximum uncached physical residual evaluations, including the initial
    /// guess, every stencil/line-search trial and the final physical recheck.
    pub maximum_evaluations: usize,
}
impl EvaporationImplicitControls {
    /// Numerical defaults with an explicitly supplied physical liquid domain.
    /// No default material interval is inferred. Time refinement and a suitable
    /// Jacobian difference scale remain caller choices.
    pub fn new(pure_liquid_enthalpy_range_j_kg: [f64; 2]) -> Self {
        let common = EvaporationNetworkImplicitControls::default();
        Self {
            pure_liquid_enthalpy_range_j_kg,
            cell_energy_tolerance_j: common.cell_energy_tolerance_j,
            total_energy_tolerance_j: common.total_energy_tolerance_j,
            mass_tolerance_kg: common.mass_tolerance_kg,
            newton: common.newton,
            difference_step_j: common.difference_step_j,
            maximum_iterations: common.maximum_iterations,
            maximum_evaluations: common.maximum_evaluations,
        }
    }

    /// Explicitly copy numerical work and acceptance gates for a network solve.
    /// The old single-cell domain is separate: the compatibility entry point
    /// always retains it in that cell's binding. Network callers supply every
    /// domain through their binding slice.
    pub const fn network_controls(self) -> EvaporationNetworkImplicitControls {
        EvaporationNetworkImplicitControls {
            cell_energy_tolerance_j: self.cell_energy_tolerance_j,
            total_energy_tolerance_j: self.total_energy_tolerance_j,
            mass_tolerance_kg: self.mass_tolerance_kg,
            newton: self.newton,
            difference_step_j: self.difference_step_j,
            maximum_iterations: self.maximum_iterations,
            maximum_evaluations: self.maximum_evaluations,
        }
    }

    pub(super) fn validate(self) -> Result<(), EvaporationImplicitError> {
        self.network_controls().validate()?;
        let [lo, hi] = self.pure_liquid_enthalpy_range_j_kg;
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return Err(EvaporationImplicitError::Invalid(
                "invalid implicit evaporation controls",
            ));
        }
        Ok(())
    }
}

/// Original physical refusal or a bounded numerical solve that did not finish.
#[derive(Debug)]
pub enum EvaporationImplicitError {
    Physical(ThermalEvaporationError),
    Numerical(NewtonError),
    Invalid(&'static str),
    NonConverged(NewtonReport),
    Interrupted(SolverRunProgress),
}
impl From<ThermalEvaporationError> for EvaporationImplicitError {
    fn from(error: ThermalEvaporationError) -> Self {
        Self::Physical(error)
    }
}
impl From<DepositionError> for EvaporationImplicitError {
    fn from(error: DepositionError) -> Self {
        Self::Physical(error.into())
    }
}
impl From<NewtonError> for EvaporationImplicitError {
    fn from(error: NewtonError) -> Self {
        Self::Numerical(error)
    }
}
impl core::fmt::Display for EvaporationImplicitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit evaporation solve: {self:?}")
    }
}
impl std::error::Error for EvaporationImplicitError {}

/// Accepted physical interval and the numerical work spent obtaining it.
#[derive(Clone, Debug, PartialEq)]
pub struct EvaporationImplicitReport {
    pub physical: EvaporatingBackwardEulerAudit,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    /// Euclidean residual norm [J] reported by the original numerical solver.
    /// Independent local, total and mass gates also passed before publication.
    pub maximum_accepted_residual_norm_j: f64,
}

/// One accepted heat-network interval with every selected surface receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct EvaporationNetworkImplicitReport {
    pub physical: EvaporatingNetworkBackwardEulerAudit,
    pub newton_iterations: usize,
    pub residual_evaluations: usize,
    pub maximum_accepted_residual_norm_j: f64,
}

impl<'m> DepositionSimulation<'m> {
    /// The original single-cell entry point, delegated to the shared network
    /// solve with exactly one binding. The supplied liquid range, typed source
    /// callback, final receipt and `L: ?Sized` support retain their meaning.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn solve_implicit_evaporating_step<L: EvaporationLaw + ?Sized, T>(
        &self,
        id: CellId,
        law: &L,
        film: GasFilm,
        target_s: f64,
        controls: EvaporationImplicitControls,
        mut resolve: impl FnMut(
            &DepositionSimulation<'m>,
            LiquidSurfaceState,
        )
            -> Result<BoundaryExchangeStepCoefficients<T>, ThermalEvaporationError>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, T, EvaporationImplicitReport), EvaporationImplicitError> {
        checkpoint()?;
        controls.validate()?;
        let binding = EvaporatingCellBinding {
            id,
            law,
            film,
            pure_liquid_enthalpy_range_j_kg: controls.pure_liquid_enthalpy_range_j_kg,
        };
        let (next, payload, report) = self.solve_implicit_evaporating_network_step(
            &[binding],
            target_s,
            controls.network_controls(),
            |thermal, liquids| resolve(thermal, liquids[&id]),
            &mut checkpoint,
        )?;
        Ok((
            next,
            payload,
            EvaporationImplicitReport {
                physical: EvaporatingBackwardEulerAudit::from_network(&report.physical, id),
                newton_iterations: report.newton_iterations,
                residual_evaluations: report.residual_evaluations,
                maximum_accepted_residual_norm_j: report.maximum_accepted_residual_norm_j,
            },
        ))
    }

    /// Automatically solve one original heat graph with all selected liquid
    /// surfaces simultaneously. Unselected cells exchange heat without phase
    /// mass transfer; empty bindings explicitly select a heat-only solve.
    ///
    /// Process start-time events first. Every selected liquid cell must be active;
    /// an interior event refuses and events exactly at target remain pending.
    /// The supplied callback has the original residual adapter's meanings: the
    /// thermal state retains initial mass for energy coordinates, while the
    /// canonical liquid map carries all proposed new masses. Return independent
    /// sources, contacts, boundary exchanges and conservative transfers only.
    ///
    /// Callbacks must keep physical history immutable during speculative calls.
    /// A two-entry residual cache avoids duplicate endpoint reconstructions in
    /// numerical Jacobian actions; the final physical endpoint and its payload
    /// are always evaluated afresh. Source/domain/stencil failures propagate
    /// unchanged and do not become zero residuals or different physical models.
    /// Finite differences supply numerical tangents, not analytic adjoints.
    ///
    /// No explicit thermal CFL is imposed, but timestep accuracy and a supported
    /// nonlinear path remain separate obligations. Failure or cancellation
    /// returns no partial timeline or phase state. All original acceptance gates
    /// remain authoritative even if the numerical solver reports convergence.
    #[allow(clippy::type_complexity)]
    pub fn solve_implicit_evaporating_network_step<L: EvaporationLaw + ?Sized, T>(
        &self,
        bindings: &[EvaporatingCellBinding<'_, L>],
        target_s: f64,
        controls: EvaporationNetworkImplicitControls,
        mut resolve: impl FnMut(
            &DepositionSimulation<'m>,
            &BTreeMap<CellId, LiquidSurfaceState>,
        )
            -> Result<BoundaryExchangeStepCoefficients<T>, ThermalEvaporationError>,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Self, T, EvaporationNetworkImplicitReport), EvaporationImplicitError> {
        checkpoint()?;
        controls.validate()?;
        let step =
            EvaporatingNetworkBackwardEulerStep::new(self, bindings, target_s, &mut checkpoint)?;
        let poll = RefCell::new(&mut checkpoint);
        let source = RefCell::new(&mut resolve);
        let evaluations = Cell::new(0usize);
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
                    |state, liquid| (source.borrow_mut())(state, liquid),
                    || (poll.borrow_mut())(),
                )?;
                let residual = trial.residual_j().to_vec();
                let mut entries = cache.borrow_mut();
                if entries.len() == 2 {
                    entries.remove(0);
                }
                entries.push((coordinates, residual.clone()));
                Ok::<_, ThermalEvaporationError>(residual)
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
            return Err(EvaporationImplicitError::Interrupted(run.progress));
        }
        if !run.report.converged {
            return Err(EvaporationImplicitError::NonConverged(run.report));
        }
        spend_evaluation(&evaluations, controls.maximum_evaluations)?;
        let trial = step.evaluate(
            &step.resolvable_increments_j(&solver.x),
            |state, liquid| (source.borrow_mut())(state, liquid),
            || (poll.borrow_mut())(),
        )?;
        let (next, payload, physical) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.total_energy_tolerance_j,
            controls.mass_tolerance_kg,
            || (poll.borrow_mut())(),
        )?;
        (poll.borrow_mut())()?;
        Ok((
            next,
            payload,
            EvaporationNetworkImplicitReport {
                physical,
                newton_iterations: run.report.iterations,
                residual_evaluations: evaluations.get(),
                maximum_accepted_residual_norm_j: run.report.residual_norm,
            },
        ))
    }
}

fn spend_evaluation(count: &Cell<usize>, limit: usize) -> Result<(), ThermalEvaporationError> {
    if count.get() >= limit {
        return Err(DepositionError::Budget("implicit evaporation residual evaluations").into());
    }
    count.set(count.get() + 1);
    Ok(())
}
fn same_point(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits())
}
fn difference_error(
    error: FiniteDifferenceError<ThermalEvaporationError>,
) -> EvaporationImplicitError {
    match error {
        FiniteDifferenceError::Source(error) => error.into(),
        FiniteDifferenceError::Invalid(message) => EvaporationImplicitError::Invalid(message),
    }
}
