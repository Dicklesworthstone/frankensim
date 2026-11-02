//! Source-bound prescribed-pressure windows using the shared diffusion driver.

use super::slab::adaptive::{
    self, MembraneDiffusionAdvanceCandidate, MembraneDiffusionAdvanceControls,
    MembraneDiffusionAdvanceError, MembraneDiffusionAdvanceWorkspace, Native,
};
use super::{
    MembraneDiffusionSolve, MembraneDiffusionState, MembraneDiffusionStep,
    MembraneDiffusionWorkspace, ResolvedMembraneDiffusion, ResolvedMembraneDiffusionError,
};

impl<'a> Native for ResolvedMembraneDiffusion<'a> {
    type Workspace = MembraneDiffusionWorkspace;
    type Error = ResolvedMembraneDiffusionError;
    fn membrane(&self) -> &MembraneDiffusionState {
        self.state()
    }
    fn is_cancelled(error: &Self::Error) -> bool {
        *error == ResolvedMembraneDiffusionError::Cancelled
    }
    fn native_step(
        &self,
        reservoirs: [f64; 2],
        duration: f64,
        solve: MembraneDiffusionSolve,
        work: &mut Self::Workspace,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, MembraneDiffusionStep), Self::Error> {
        Ok(self
            .step_candidate(reservoirs, duration, solve, work, cancelled)?
            .into_parts())
    }
}

impl<'a> ResolvedMembraneDiffusion<'a> {
    pub fn adaptive_workspace(
        &self,
    ) -> Result<MembraneDiffusionAdvanceWorkspace, ResolvedMembraneDiffusionError> {
        Ok(self.state().adaptive_workspace()?)
    }
    /// Complete one constant prescribed-pressure window with the original
    /// selected D/S/interface laws checked on EVERY coarse and fine probe.
    /// Only accepted native halves advance state and boundary inventory. Source
    /// refusals are hard errors, never reasons to silently reduce the time step.
    /// The final state retains its original binding and final-step support;
    /// temperature and coefficients remain fixed. Split calls at boundary events.
    pub fn advance_candidate(
        &self,
        reservoirs_pa: [f64; 2],
        duration_s: f64,
        controls: MembraneDiffusionAdvanceControls,
        work: &mut MembraneDiffusionAdvanceWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        MembraneDiffusionAdvanceCandidate<Self>,
        MembraneDiffusionAdvanceError<ResolvedMembraneDiffusionError>,
    > {
        adaptive::advance(self, reservoirs_pa, duration_s, controls, work, cancelled)
    }
}
