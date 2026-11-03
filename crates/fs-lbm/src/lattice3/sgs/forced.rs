//! Exploratory SGS on the original spatial/time-dependent force owner.
//!
//! Every closure sees `(momentum + F_n/2)/rho` from the accepted population
//! and force epoch. A supplied `F_(n+1)` belongs only to the original force
//! transaction: collision uses `F_n`, endpoint admission/readback use
//! `F_(n+1)`, and populations, force history and clock publish together.
//! Reinterpreting the incoming field with the endpoint force would add an
//! artificial strain and alter both local and dynamic transport.
//!
//! These adapters reuse the existing gradients, fits, molecular conversion,
//! reports and workspaces. They allocate no additional preparation arrays and
//! do not expose the underlying grid without its force history. Molecular
//! inputs must remain original; never resubmit previously augmented rates.

use super::*;
use crate::lattice3::sources::{ForcedGrid3, ForcedLattice3};
use crate::sgs3::dynamic::DynamicSgsConfig;
use crate::sgs3::dynamic::channel::ChannelDynamicSgsConfig;
use crate::sgs3::{ChannelDynamicSgsWorkspace3, DynamicSgsStepReport, DynamicSgsWorkspace3};

fn gate_checkpoint<E>(gate: &CancelGate) -> Result<(), GridError3<E>> {
    if gate.is_requested() {
        Err(GridError3::Cancelled)
    } else {
        Ok(())
    }
}

impl<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Original local closure workspace; no force snapshot or second grid.
    pub fn sgs_workspace(&self) -> Result<SgsWorkspace3<L::Collision>, GridError3<L::Error>> {
        self.grid.sgs_workspace()
    }

    /// Original periodic dynamic workspace and complete numerical byte gate.
    pub fn dynamic_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<DynamicSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        self.grid.dynamic_sgs_workspace(maximum_bytes)
    }

    /// Original planar-channel workspace, wall topology and numerical byte gate.
    pub fn channel_dynamic_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<ChannelDynamicSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        self.grid.channel_dynamic_sgs_workspace(maximum_bytes)
    }

    /// Prepare Smagorinsky/WALE/Vreman from the CURRENT spatial force epoch.
    /// Existing solid and stationary halfway-wall gradients remain authoritative.
    pub fn prepare_sgs(
        &self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        self.grid
            .prepare_sgs_with_forces(config, molecular, workspace, Some(self.forces()), || {
                cx.checkpoint().map_err(|_| GridError3::Cancelled)
            })
    }

    /// Hold the current force while applying one local-SGS collision field.
    pub fn step_serial_with_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        let report = self.prepare_sgs(cx, config, molecular, workspace)?;
        self.step_serial(cx, workspace.collisions())?;
        Ok(report)
    }

    /// Fit at the current epoch, then atomically apply the NEXT force endpoint.
    pub fn advance_serial_with_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        let report = self.prepare_sgs(cx, config, molecular, workspace)?;
        self.advance_serial(cx, workspace.collisions(), forces_at_end)?;
        Ok(report)
    }

    /// Pooled original transaction with a held spatial force and local SGS.
    pub fn step_pooled_with_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        let report = self.grid.prepare_sgs_with_forces(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            || gate_checkpoint(gate),
        )?;
        self.step_pooled(pool, gate, workspace.collisions())?;
        Ok(report)
    }

    /// Pooled local-SGS advance; failure retains both populations and force epoch.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        let report = self.grid.prepare_sgs_with_forces(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            || gate_checkpoint(gate),
        )?;
        self.advance_pooled(pool, gate, workspace.collisions(), forces_at_end)?;
        Ok(report)
    }

    /// Periodic dynamic fit from current physical velocities and density.
    pub fn prepare_dynamic_sgs(
        &self,
        cx: &Cx<'_>,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        self.grid.prepare_dynamic_sgs_with_forces(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Hold the current force for one original step with periodic dynamic SGS.
    pub fn step_serial_with_dynamic_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.prepare_dynamic_sgs(cx, config, molecular, workspace)?;
        self.step_serial(cx, workspace.collisions())?;
        Ok(report)
    }

    /// Incoming-time dynamic fit with the original trapezoidal force advance.
    pub fn advance_serial_with_dynamic_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.prepare_dynamic_sgs(cx, config, molecular, workspace)?;
        self.advance_serial(cx, workspace.collisions(), forces_at_end)?;
        Ok(report)
    }

    /// Pooled original step with current-epoch periodic dynamic transport.
    pub fn step_pooled_with_dynamic_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.grid.prepare_dynamic_sgs_with_forces(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            || gate_checkpoint(gate),
        )?;
        self.step_pooled(pool, gate, workspace.collisions())?;
        Ok(report)
    }

    /// Pooled dynamic advance; the next force never enters the incoming fit.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_dynamic_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: DynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut DynamicSgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.grid.prepare_dynamic_sgs_with_forces(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            || gate_checkpoint(gate),
        )?;
        self.advance_pooled(pool, gate, workspace.collisions(), forces_at_end)?;
        Ok(report)
    }

    /// Channel fit with current spatial forces and original stationary walls.
    pub fn prepare_channel_dynamic_sgs(
        &self,
        cx: &Cx<'_>,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        self.grid.prepare_channel_dynamic_sgs_with_forces(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Hold the current force for one original channel-dynamic step.
    pub fn step_serial_with_channel_dynamic_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.prepare_channel_dynamic_sgs(cx, config, molecular, workspace)?;
        self.step_serial(cx, workspace.collisions())?;
        Ok(report)
    }

    /// Current-force channel transport followed by the original endpoint advance.
    pub fn advance_serial_with_channel_dynamic_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.prepare_channel_dynamic_sgs(cx, config, molecular, workspace)?;
        self.advance_serial(cx, workspace.collisions(), forces_at_end)?;
        Ok(report)
    }

    /// Pooled channel step with unchanged boundary and source ownership.
    pub fn step_pooled_with_channel_dynamic_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.grid.prepare_channel_dynamic_sgs_with_forces(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            || gate_checkpoint(gate),
        )?;
        self.step_pooled(pool, gate, workspace.collisions())?;
        Ok(report)
    }

    /// Pooled channel advance; both source endpoints retain their original roles.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_channel_dynamic_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.grid.prepare_channel_dynamic_sgs_with_forces(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            || gate_checkpoint(gate),
        )?;
        self.advance_pooled(pool, gate, workspace.collisions(), forces_at_end)?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
