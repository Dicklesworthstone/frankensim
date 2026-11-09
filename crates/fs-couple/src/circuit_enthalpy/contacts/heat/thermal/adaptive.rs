//! Complete nonlinear material heat through the original step-doubling owner.

use super::super::{HeatAssembly, as_circuit_error};
use super::{HeatAssemblySolve, step_with_entropy};
use crate::circuit_enthalpy::{CircuitEnthalpyError, entropy::CircuitEntropyAudit};
use fs_exec::Cx;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_solver::NewtonReport;
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionSimulation, EnergyAudit, ThermalLaw,
    conduction::implicit::{
        ImplicitConductionControls, ImplicitConductionError, ImplicitConductionObservables,
        ImplicitConductionReport,
        adaptive::{
            AdaptiveConductionStepper, AdaptiveImplicitConductionControls,
            AdaptiveImplicitConductionError, AdaptiveImplicitConductionReport,
        },
    },
    entropy::EntropyAudit,
};
use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

/// Physical-time accuracy and the existing thermal nonlinear policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveHeatAssemblyControls {
    /// The evaluation ceiling applies to the WHOLE advance: coarse, rejected,
    /// accepted fine and final-source queries all share it.
    pub step: HeatAssemblySolve,
    pub initial_dt_s: f64,
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    pub maximum_steps: usize,
    pub maximum_halvings: usize,
    pub maximum_trials: usize,
    pub enthalpy_absolute_tolerance_j_kg: f64,
    pub temperature_absolute_tolerance_k: f64,
    pub phase_fraction_absolute_tolerance: Option<f64>,
    /// Scales changes from the last accepted state, never the energy datum.
    pub relative_tolerance: f64,
}
impl AdaptiveHeatAssemblyControls {
    fn native(self) -> AdaptiveImplicitConductionControls {
        AdaptiveImplicitConductionControls {
            step: ImplicitConductionControls {
                maximum_dt_s: self.maximum_dt_s,
                minimum_dt_s: self.minimum_dt_s,
                maximum_steps: self.maximum_steps,
                maximum_newton_iterations: self.step.maximum_newton_iterations,
                // The endpoint's original Newton globalization owns this work.
                maximum_backtracks: 0,
                maximum_halvings: self.maximum_halvings,
                maximum_residual_evaluations: self.step.maximum_evaluations,
                cell_tolerance_j: self.step.cell_energy_tolerance_j,
                total_tolerance_j: self.step.closed_energy_tolerance_j,
            },
            initial_dt_s: self.initial_dt_s,
            maximum_trials: self.maximum_trials,
            enthalpy_absolute_tolerance_j_kg: self.enthalpy_absolute_tolerance_j_kg,
            temperature_absolute_tolerance_k: self.temperature_absolute_tolerance_k,
            phase_fraction_absolute_tolerance: self.phase_fraction_absolute_tolerance,
            relative_tolerance: self.relative_tolerance,
        }
    }
}

/// Only accepted fine intervals contribute physical energy and entropy.
pub struct AdaptiveHeatAssemblyReport<S> {
    /// Original fine-trajectory account and local error estimate. Its residual
    /// evaluation count includes the final complete-source query as source work.
    pub time: AdaptiveImplicitConductionReport,
    /// Fresh complete source at the returned endpoint, AFTER any target event.
    pub source: S,
}

#[derive(Debug)]
pub enum AdaptiveHeatAssemblyError<E> {
    Time(AdaptiveImplicitConductionError),
    /// Original material, source, nonlinear or work failure, without erasure.
    Solve(CircuitEnthalpyError<E>),
}
impl<E: core::fmt::Debug> core::fmt::Display for AdaptiveHeatAssemblyError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive complete material heat: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for AdaptiveHeatAssemblyError<E> {}

struct Stepper<'a, 'scope, 'm, A: HeatAssembly + ?Sized> {
    assembly: &'a A,
    cx: &'a Cx<'scope>,
    initial: &'a DepositionSimulation<'m>,
    controls: HeatAssemblySolve,
    ids: Vec<CellId>,
    curves: BTreeMap<CellId, &'a EquilibriumEnthalpyPhaseCurve>,
    transitions: &'a [f64],
    evaluations: Cell<usize>,
    failure: RefCell<Option<CircuitEnthalpyError<A::Error>>>,
    last_nonconvergence: RefCell<Option<CircuitEnthalpyError<A::Error>>>,
}
impl<A: HeatAssembly + ?Sized> Stepper<'_, '_, '_, A> {
    fn charge(&self) -> Result<(), CircuitEnthalpyError<A::Error>> {
        self.cx
            .checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let count = self.evaluations.get();
        if count >= self.controls.maximum_evaluations {
            return Err(CircuitEnthalpyError::CircuitBudget);
        }
        self.evaluations.set(count + 1);
        Ok(())
    }

    fn fatal(&self, error: CircuitEnthalpyError<A::Error>) -> ImplicitConductionError {
        *self.failure.borrow_mut() = Some(error);
        // The outer boundary restores the original typed failure immediately.
        // This nonretryable signal never becomes a numerical time-step retry.
        ImplicitConductionError::Invalid("complete heat assembly refused its original trial")
    }
}

fn record_work(
    work: &mut ImplicitConductionReport,
    report: &NewtonReport,
) -> Result<(), ImplicitConductionError> {
    work.nonlinear_iterations = work
        .nonlinear_iterations
        .checked_add(report.iterations)
        .ok_or(DepositionError::Arithmetic)?;
    for iteration in &report.history {
        work.linear_iterations = work
            .linear_iterations
            .checked_add(iteration.linear_iterations)
            .ok_or(DepositionError::Arithmetic)?;
    }
    Ok(())
}

// The private thermal-only driver has no electrical loss, reversible transfer
// or advected store. Preserve the original entropy terms and their explicit
// unclassified-source flag; ordinary scalar heat does not acquire bath entropy.
fn thermal_entropy(source: CircuitEntropyAudit) -> Result<EntropyAudit, DepositionError> {
    let temporal_defect_j_k = source.thermal_storage_change_j_k
        - source.reservoir_entropy_input_j_k
        - source.conduction_generation_j_k
        - source.unclassified_source_at_cell_j_k;
    if !temporal_defect_j_k.is_finite() {
        return Err(DepositionError::Arithmetic);
    }
    Ok(EntropyAudit {
        thermal_storage_change_j_k: source.thermal_storage_change_j_k,
        reservoir_entropy_input_j_k: source.reservoir_entropy_input_j_k,
        conduction_generation_j_k: source.conduction_generation_j_k,
        advection: None,
        unclassified_source_at_cell_j_k: source.unclassified_source_at_cell_j_k,
        has_unclassified_sources: source.has_unclassified_sources,
        temporal_defect_j_k,
    })
}

impl<A: HeatAssembly + ?Sized> AdaptiveConductionStepper for Stepper<'_, '_, '_, A> {
    type Workspace = Vec<f64>;
    type Observer = ();

    fn cell_ids(&self) -> &[CellId] {
        &self.ids
    }
    fn workspace_cells(workspace: &Vec<f64>) -> usize {
        workspace.len()
    }
    fn observer(&self) {}
    fn next_transition_time_s(&self, after_s: f64) -> Option<f64> {
        self.transitions
            .get(self.transitions.partition_point(|&time| time <= after_s))
            .copied()
    }

    fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        checkpoint()?;
        if thermal.cells().len() != self.ids.len() || thermal.links() != self.initial.links() {
            return Err(ImplicitConductionError::Invalid(
                "adaptive heat requires its original fixed topology",
            ));
        }
        for (&id, original) in self.initial.cells() {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(ImplicitConductionError::Invalid(
                    "adaptive heat cell identity changed",
                ))?;
            let curve = self
                .curves
                .get(&id)
                .ok_or(ImplicitConductionError::Invalid(
                    "adaptive heat needs every original material curve",
                ))?;
            if cell.material_id != original.material_id
                || cell.mass_kg != original.mass_kg
                || cell.geometry != original.geometry
                || cell.body_id != original.body_id
                || cell.segment_id != original.segment_id
                || !thermal.uses_material_law(cell.material_id, *curve as &dyn ThermalLaw)
            {
                return Err(ImplicitConductionError::Invalid(
                    "adaptive heat must retain original finite stores",
                ));
            }
        }
        for &link in self.initial.links().keys() {
            checkpoint()?;
            if thermal.contact_patch(link) != self.initial.contact_patch(link) {
                return Err(ImplicitConductionError::Invalid(
                    "adaptive heat contact coverage changed",
                ));
            }
        }
        Ok(())
    }

    fn observable(
        &self,
        thermal: &DepositionSimulation<'_>,
        id: CellId,
    ) -> Result<ImplicitConductionObservables, DepositionError> {
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(DepositionError::MaterialDomain)?;
        let state = self
            .curves
            .get(&id)
            .ok_or(DepositionError::MaterialDomain)?
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .map_err(|_| DepositionError::MaterialDomain)?;
        Ok(ImplicitConductionObservables {
            temperature_k: thermal.temperature_k(id)?,
            phase_mass_fraction: Some(state.liquid_mass_fraction()),
        })
    }

    fn solve_interval<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        _: ImplicitConductionControls,
        workspace: &mut Vec<f64>,
        work: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, EnergyAudit, EntropyAudit, f64), ImplicitConductionError>
    {
        checkpoint()?;
        let mut controls = self.controls;
        controls.maximum_evaluations = controls
            .maximum_evaluations
            .saturating_sub(self.evaluations.get());
        if controls.maximum_evaluations == 0 {
            return Err(self.fatal(CircuitEnthalpyError::CircuitBudget));
        }
        let result = step_with_entropy(
            self.assembly,
            self.cx,
            thermal,
            target_s,
            controls,
            workspace,
            true,
            || self.charge(),
        );
        work.residual_evaluations = self.evaluations.get();
        match result {
            Ok((next, report, entropy)) => {
                record_work(work, &report.nonlinear)?;
                self.last_nonconvergence.borrow_mut().take();
                let entropy = entropy.ok_or(ImplicitConductionError::Invalid(
                    "missing original thermal entropy account",
                ))?;
                let mut energy = report.thermal;
                energy.external_heat_j = report.external_heat_j;
                energy.residual_j = report.energy_residual_j;
                checkpoint()?;
                Ok((
                    next,
                    energy,
                    thermal_entropy(entropy)?,
                    report.maximum_cell_residual_j,
                ))
            }
            Err(CircuitEnthalpyError::NotConverged(report)) => {
                record_work(work, &report)?;
                // Only a retry signal: the original full numerical report is
                // restored on exhaustion. Its norm bounds each local residual.
                let maximum_residual_j = report.residual_norm;
                *self.last_nonconvergence.borrow_mut() =
                    Some(CircuitEnthalpyError::NotConverged(report));
                Err(ImplicitConductionError::NotConverged { maximum_residual_j })
            }
            Err(error) => Err(self.fatal(error)),
        }
    }
}

/// Advance a complete nonlinear heat assembly with local time-error control.
///
/// The ORIGINAL native controller compares one full interval with two half
/// intervals and accepts only the conservative fine trajectory. Every source,
/// wall/interface temperature and finite-fluid convection coefficient is
/// re-evaluated inside each native nonlinear solve. Coarse/rejected heat never
/// enters storage or the accepted energy/entropy report; no extrapolated state
/// or additional physical store is constructed.
///
/// Native source/bath boundary events are applied at their exact times. The
/// optional strictly increasing finite source-transition roster also clips
/// comparisons for continuous callback-owned slope changes. Discontinuous
/// source/bath changes use native events; callbacks must retain a repeatable
/// declared time convention and may not advance history.
/// Birth, separation, contact geometry or finite-inventory changes require a
/// new prepared owner and refuse within this fixed-topology advance. All active
/// original caloric curves provide h/T and optional liquid-fraction error scales.
///
/// The adaptive path computes original caloric entropy only at each interval's
/// final admission. Independent callback sources, including prescribed-fluid
/// powers, remain explicitly unclassified at-cell entropy; no bath temperature
/// or irreversible production is inferred. The ordinary one-step energy-only
/// API is unchanged. The source evaluation cap includes ALL trial/final queries.
///
/// Only native numerical nonconvergence and temporal error retry. Original
/// source/phase/energy failures, cancellation and exhausted work remain typed
/// refusals. Even a late failure returns no accepted prefix; source callbacks
/// themselves must be pure. Returned source evidence is fresh after any event
/// at the requested endpoint. Local error estimates are not global, spatial or
/// empirical-model accuracy guarantees; small-source and algebraic tolerances
/// must be chosen consistently with requested temporal accuracy.
#[allow(clippy::too_many_arguments)]
pub fn advance_heat_assembly_adaptive_candidate<'m, A: HeatAssembly + ?Sized>(
    assembly: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'m>,
    target_s: f64,
    controls: AdaptiveHeatAssemblyControls,
    source_transitions_s: &[f64],
) -> Result<
    (
        DepositionSimulation<'m>,
        AdaptiveHeatAssemblyReport<A::Sample>,
    ),
    AdaptiveHeatAssemblyError<A::Error>,
> {
    use AdaptiveHeatAssemblyError::{Solve, Time};
    cx.checkpoint()
        .map_err(|_| Solve(CircuitEnthalpyError::Cancelled))?;
    if controls.step.maximum_cells == 0
        || thermal.cells().len() > controls.step.maximum_cells
        || !controls.step.difference_step_j.is_finite()
        || controls.step.difference_step_j <= 0.0
        || !controls.step.thermal_energy_tolerance_j.is_finite()
        || controls.step.thermal_energy_tolerance_j < 0.0
        || source_transitions_s.len() > controls.maximum_trials
        || source_transitions_s.iter().any(|time| !time.is_finite())
        || source_transitions_s
            .windows(2)
            .any(|pair| pair[0] >= pair[1])
    {
        return Err(Solve(CircuitEnthalpyError::Invalid(
            "bounded original cells and a finite ordered source-transition roster are required",
        )));
    }
    controls
        .step
        .newton
        .validate()
        .map_err(CircuitEnthalpyError::from)
        .map_err(Solve)?;
    let material = assembly
        .prepared_heat_network()
        .thermal_material_network()
        .map_err(|error| Solve(as_circuit_error(error).map_source(A::Error::from)))?;
    let stepper = Stepper {
        assembly,
        cx,
        initial: thermal,
        controls: controls.step,
        ids: thermal.cells().keys().copied().collect(),
        curves: material
            .cells
            .iter()
            .map(|cell| (cell.cell, cell.curve))
            .collect(),
        transitions: source_transitions_s,
        evaluations: Cell::new(0),
        failure: RefCell::new(None),
        last_nonconvergence: RefCell::new(None),
    };
    let mut workspace = Vec::new();
    workspace
        .try_reserve_exact(stepper.ids.len())
        .map_err(|_| Time(DepositionError::Budget("adaptive heat guess storage").into()))?;
    workspace.resize(stepper.ids.len(), 0.0);
    let result =
        stepper.advance_adaptive_to(thermal, target_s, controls.native(), &mut workspace, || {
            cx.checkpoint().map_err(|_| DepositionError::Cancelled)
        });
    if let Some(error) = stepper.failure.borrow_mut().take() {
        return Err(Solve(error));
    }
    let (next, mut time) = match result {
        Err(
            error @ AdaptiveImplicitConductionError::Implicit(
                ImplicitConductionError::NotConverged { .. },
            ),
        ) => {
            return Err(match stepper.last_nonconvergence.borrow_mut().take() {
                Some(original) => Solve(original),
                None => Time(error),
            });
        }
        other => other.map_err(Time)?,
    };
    stepper.charge().map_err(Solve)?;
    let endpoint = assembly.resolve_heat(cx, &next).map_err(Solve)?;
    if !endpoint.bulk.resistors.is_empty() || !endpoint.contacts.is_empty() {
        return Err(Solve(CircuitEnthalpyError::Invalid(
            "adaptive thermal-only assembly cannot discard electrical coefficients",
        )));
    }
    time.solve.residual_evaluations = stepper.evaluations.get();
    cx.checkpoint()
        .map_err(|_| Solve(CircuitEnthalpyError::Cancelled))?;
    Ok((
        next,
        AdaptiveHeatAssemblyReport {
            time,
            source: endpoint.bulk.source,
        },
    ))
}
