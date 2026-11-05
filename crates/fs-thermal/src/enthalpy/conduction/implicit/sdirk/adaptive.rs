//! SDIRK2 endpoints in the existing h/T/phase step-doubling controller.
use super::super::adaptive::{
    AdaptiveConductionStepper, AdaptiveImplicitConductionControls, AdaptiveImplicitConductionError,
    AdaptiveImplicitConductionReport, ConductionTimeOrder,
};
use super::*;

// Private so callers cannot bypass the whole-window event preflight.
struct Sdirk2Stepper<'s, 'law> {
    solver: &'s ImplicitConductionChain<'law>,
}

impl<'s, 'law> AdaptiveConductionStepper for Sdirk2Stepper<'s, 'law> {
    type Workspace = ImplicitConductionWorkspace;
    type Observer = ConductionObserver<'law>;

    fn cell_ids(&self) -> &[CellId] {
        self.solver.cell_ids()
    }

    fn workspace_cells(workspace: &Self::Workspace) -> usize {
        workspace.cells()
    }

    fn observer(&self) -> Self::Observer {
        self.solver.conduction.clone()
    }

    fn time_order(&self) -> ConductionTimeOrder {
        ConductionTimeOrder::Second
    }

    fn admit(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ImplicitConductionError> {
        ImplicitConductionChain::admit(self.solver, thermal, checkpoint)
    }

    fn observable(
        &self,
        thermal: &DepositionSimulation<'_>,
        id: CellId,
    ) -> Result<ImplicitConductionObservables, DepositionError> {
        AdaptiveConductionStepper::observable(self.solver, thermal, id)
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
        self.solver
            .solve_sdirk2_step(thermal, target_s, controls, workspace, work, checkpoint)
    }
}

impl ImplicitConductionChain<'_> {
    /// Advance native SDIRK2 with the existing h/T/phase local-error controller.
    ///
    /// The entire requested window must satisfy [`Self::advance_sdirk2_to`]'s
    /// fixed-mass, constant-source/bath scope before the generic controller
    /// can process any events. Attached implicit advection and pending events
    /// through the target refuse even if the first comparison would end sooner.
    ///
    /// Each comparison computes one full and two half SDIRK2 intervals; all six
    /// stages, rejected comparisons and entropy evaluations consume the same
    /// work budgets. Only the two physical fine endpoints contribute heat,
    /// entropy and accepted-step counts. No extrapolated state is published.
    ///
    /// Actual half spans a,b give fine-error correction
    /// `(a³+b³)/(3ab(a+b))`, equal to 1/3 for equal halves. The controller uses
    /// the corresponding cubic local-error scaling, with bounded deterministic
    /// arithmetic. Caloric reference-invariant h changes and original T/phase
    /// observables retain their existing tolerances. This is a local estimate,
    /// not a global error bound or a positivity guarantee; constitutive corners
    /// still require independent refinement. Failure publishes no partial run.
    pub fn advance_sdirk2_adaptive_to<'m>(
        &self,
        thermal: &DepositionSimulation<'m>,
        target_s: f64,
        controls: AdaptiveImplicitConductionControls,
        workspace: &mut ImplicitConductionWorkspace,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (DepositionSimulation<'m>, AdaptiveImplicitConductionReport),
        AdaptiveImplicitConductionError,
    > {
        checkpoint()?;
        self.admit_sdirk2_window(thermal, target_s)?;
        AdaptiveConductionStepper::advance_adaptive_to(
            &Sdirk2Stepper { solver: self },
            thermal,
            target_s,
            controls,
            workspace,
            checkpoint,
        )
    }
}
