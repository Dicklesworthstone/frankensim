//! Local open-grid SGS in the accepted spatial Guo force frame.

use super::*;
use crate::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;

impl<L: EquilibriumLattice3<Q> + SgsLattice3<Q>, const Q: usize>
    ForcedOpenBoundaryGrid3<L, Q>
{
    /// Numerical payload for the original open-grid SGS workspace.
    pub fn required_sgs_workspace_bytes(&self) -> Result<usize, SgsStepError<L::Error>> {
        self.inner.required_sgs_workspace_bytes()
    }

    /// Reserve one reusable workspace; no force or population state is copied.
    pub fn sgs_workspace(
        &self,
        maximum_workspace_bytes: usize,
    ) -> Result<OpenSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        self.inner.sgs_workspace(maximum_workspace_bytes)
    }

    /// Prepare local stress from the SAME accepted populations and spatial
    /// force. Physical velocity includes each cell's own Guo half-source;
    /// reservoir and wall derivatives retain their original native geometry.
    /// A future force field must never be used to interpret this opening state.
    ///
    /// The optional collision field is fresh molecular transport. Original
    /// bulk/higher-mode rates and all shared closure/ratio limits are retained.
    pub fn prepare_sgs(
        &self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut OpenSgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        self.inner.prepare_open_sgs(
            config,
            molecular,
            workspace,
            Some(&self.forces),
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Advance one native interval with the held spatial force and incoming SGS.
    /// Preparation changes scratch only. The original native transaction owns
    /// population/clock/wall publication and every late refusal.
    pub fn step_serial_with_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut OpenSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
    ) -> Result<OpenSgsStep3, SgsStepError<L::Error>> {
        let sgs = self.prepare_sgs(cx, config, molecular, workspace)?;
        let momentum = self.step_serial(cx, Some(workspace.collisions()), origin)?;
        Ok(OpenSgsStep3 { momentum, sgs })
    }

    /// Advance with a prescribed force at the NEXT accepted time. SGS samples
    /// only the retained opening force; native collision uses that same epoch,
    /// while streamed donors and reservoirs use the supplied endpoint force.
    /// The original trapezoidal body impulse and reservoir/wall accounts are
    /// unchanged. No force work or modeled transfer is deposited as heat.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        forces_at_end: Vec<[f64; 3]>,
        workspace: &mut OpenSgsWorkspace3<L::Collision>,
        origin: [f64; 3],
    ) -> Result<OpenSgsStep3, SgsStepError<L::Error>> {
        let sgs = self.prepare_sgs(cx, config, molecular, workspace)?;
        let momentum = self.advance_serial(
            cx,
            Some(workspace.collisions()),
            forces_at_end,
            origin,
        )?;
        Ok(OpenSgsStep3 { momentum, sgs })
    }
}
