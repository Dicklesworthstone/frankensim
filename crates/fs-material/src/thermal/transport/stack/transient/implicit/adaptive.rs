//! Material endpoints plugged into the existing thermal time-error controller.

use super::{
    BoundLayeredThermalBody, LayerConductionConfig, LayeredThermalAdaptiveControls,
    LayeredThermalAdaptiveError, LayeredThermalAdaptiveReport, LayeredThermalSolveControls,
    LayeredThermalSolveError, LayeredThermalStep, LayeredThermalWorkspace, solve::SolveWork,
};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    conduction::implicit::{
        ImplicitConductionControls, ImplicitConductionError, ImplicitConductionReport,
        adaptive::AdaptiveConductionStepper,
    },
    entropy::EntropyAudit,
};
use std::collections::BTreeMap;

impl BoundLayeredThermalBody<'_> {
    /// Advance with the original thermal owner's adaptive backward-Euler pair.
    ///
    /// Each comparison uses this body's exact integrated conductivity/contact
    /// fluxes and sparse Newton solve. Only the two conservative half steps
    /// contribute heat and entropy; no extrapolated physical state is invented.
    /// Error scales use changes in source-specific storage and temperature,
    /// independently of the original local and whole-advance joule gates.
    /// Sensible entropy must be enabled in the bound storage laws. Phase error
    /// control refuses because this sensible body has no phase observable.
    ///
    /// All cells/contacts must already be active. Scheduled Boundary events are
    /// applied at their exact times, including the target. Births, geometry,
    /// contact changes belong to an enclosing event owner and
    /// refuse here. In contrast, a single `LayeredThermalStep` leaves target
    /// events pending. Failure/cancellation returns no partial thermal prefix.
    /// The reported local estimate is not a certified global error or adjoint.
    pub fn advance_adaptive_to<'m>(
        &self,
        initial: &DepositionSimulation<'m>,
        target_s: f64,
        config: LayerConductionConfig,
        controls: LayeredThermalAdaptiveControls,
        workspace: &mut LayeredThermalWorkspace,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'m>, LayeredThermalAdaptiveReport), LayeredThermalAdaptiveError>
    {
        LayeredStepper {
            binding: self,
            ids: self.cell_ids(),
            config,
        }
        .advance_adaptive_to(initial, target_s, controls, workspace, checkpoint)
    }
}

struct LayeredStepper<'s, 'a> {
    binding: &'s BoundLayeredThermalBody<'a>,
    ids: Vec<CellId>,
    config: LayerConductionConfig,
}

/// A fixed-mass layered body permits boundary schedules, not topology changes.
#[derive(Clone)]
struct BoundaryEvents;
impl DepositionObserver for BoundaryEvents {
    type Error = DepositionError;
    fn event(&mut self, event: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        if matches!(event, Event::Boundary(..)) {
            Ok(())
        } else {
            Err(DepositionError::InvalidInput(
                "adaptive layered body only admits boundary events",
            ))
        }
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Err(DepositionError::InvalidInput(
            "layered event observer cannot advance thermal time",
        ))
    }
}

impl AdaptiveConductionStepper for LayeredStepper<'_, '_> {
    type Workspace = LayeredThermalWorkspace;
    type Observer = BoundaryEvents;

    fn cell_ids(&self) -> &[CellId] {
        &self.ids
    }
    fn workspace_cells(workspace: &Self::Workspace) -> usize {
        workspace.cells()
    }
    fn observer(&self) -> Self::Observer {
        BoundaryEvents
    }

    fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        if thermal.cells().len() != self.ids.len() || thermal.links().len() != self.ids.len() - 1 {
            return Err(ImplicitConductionError::Invalid(
                "adaptive layered solve requires exactly its bound cells and links",
            ));
        }
        self.binding
            .step_coefficients(thermal, self.config, checkpoint)?;
        Ok(())
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
        let remaining = controls
            .maximum_residual_evaluations
            .saturating_sub(work.residual_evaluations);
        if remaining == 0 {
            return Err(DepositionError::Budget("adaptive layered residual evaluations").into());
        }
        let step = LayeredThermalStep::new(
            self.binding,
            thermal,
            target_s,
            self.config,
            &mut *checkpoint,
        )?;
        let mut spent = SolveWork::default();
        let result = step.solve_inner(
            LayeredThermalSolveControls {
                maximum_newton_iterations: controls.maximum_newton_iterations,
                maximum_backtracks: controls.maximum_backtracks,
                maximum_residual_evaluations: remaining,
                cell_tolerance_j: controls.cell_tolerance_j,
                total_tolerance_j: controls.total_tolerance_j,
            },
            workspace,
            true,
            &mut spent,
            checkpoint,
        );
        work.residual_evaluations += spent.residual_evaluations;
        work.nonlinear_iterations = work
            .nonlinear_iterations
            .checked_add(spent.newton_iterations)
            .ok_or(DepositionError::Budget(
                "adaptive layered iteration counter",
            ))?;
        let (next, _, report) = result.map_err(|error| match error {
            LayeredThermalSolveError::Thermal(error) => ImplicitConductionError::Thermal(error),
            LayeredThermalSolveError::Linear(error) => ImplicitConductionError::Linear(error),
            LayeredThermalSolveError::NotConverged {
                maximum_cell_residual_j,
                ..
            } => ImplicitConductionError::NotConverged {
                maximum_residual_j: maximum_cell_residual_j,
            },
        })?;
        let entropy = report.entropy.ok_or(ImplicitConductionError::Invalid(
            "missing original layered entropy audit",
        ))?;
        Ok((next, report.energy, entropy, report.maximum_cell_residual_j))
    }
}
