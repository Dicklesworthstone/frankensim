//! Moving-wall SGS with the original spatial body-force history.
//!
//! Preparation uses the accepted F_n/2 physical velocity. The same complete
//! uniform wall pair feeds the local/dynamic closure and the actual native
//! stream; F_(n+1) is reserved for endpoint source admission. All wrappers use
//! the original workspaces and molecular baselines. No force snapshot, extra
//! population array, second collision or persistent SGS law is introduced.

use super::*;
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallMotion3};
use crate::lattice3::wall_loads::WallMomentumProbe3;

impl<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// Prepare local SGS from current spatial forces and both actual wall
    /// velocities. The complete uniform-plane admission is shared with Grid3.
    pub fn prepare_channel_sgs_with_moving_walls(
        &self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        self.grid.prepare_channel_sgs_with_motion(
            config, molecular, workspace, Some(self.forces()),
            WallMotion3 { probe, velocities },
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Hold the spatial force and advance local SGS with actual wall loads.
    pub fn step_serial_with_channel_sgs_and_moving_walls(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.prepare_channel_sgs_with_moving_walls(
            cx, config, molecular, probe, velocities, workspace,
        )?;
        let walls = self.step_serial_with_moving_wall_momentum(
            cx, probe, velocities, workspace.collisions(),
        )?;
        Ok((report, walls))
    }

    /// Fit at F_n, then publish the corrected stream, wall receipt and next
    /// force endpoint together. Failure retains the complete accepted epoch.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_channel_sgs_and_moving_walls(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.prepare_channel_sgs_with_moving_walls(
            cx, config, molecular, probe, velocities, workspace,
        )?;
        let walls = self.advance_serial_with_moving_wall_momentum(
            cx, probe, velocities, workspace.collisions(), forces_at_end,
        )?;
        Ok((report, walls))
    }

    /// Pooled local-SGS wall step with unchanged current body-force history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_channel_sgs_and_moving_walls<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.grid.prepare_channel_sgs_with_motion(
            config, molecular, workspace, Some(self.forces()),
            WallMotion3 { probe, velocities }, || gate_checkpoint(gate),
        )?;
        let walls = self.step_pooled_with_moving_wall_momentum(
            pool, gate, probe, velocities, workspace.collisions(),
        )?;
        Ok((report, walls))
    }

    /// Pooled current-epoch local fit and original trapezoidal source advance.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_channel_sgs_and_moving_walls<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.grid.prepare_channel_sgs_with_motion(
            config, molecular, workspace, Some(self.forces()),
            WallMotion3 { probe, velocities }, || gate_checkpoint(gate),
        )?;
        let walls = self.advance_pooled_with_moving_wall_momentum(
            pool, gate, probe, velocities, workspace.collisions(), forces_at_end,
        )?;
        Ok((report, walls))
    }

    /// Fit dynamic channel SGS in the current spatial-force frame, using the
    /// complete native wall pair for both filtering and shear dissipation.
    pub fn prepare_channel_dynamic_sgs_with_moving_walls(
        &self,
        cx: &Cx<'_>,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        self.grid.prepare_channel_dynamic_sgs_with_motion(
            config, molecular, workspace, Some(self.forces()),
            Some(WallMotion3 { probe, velocities }),
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Hold spatial body forces while advancing dynamic SGS and native walls.
    pub fn step_serial_with_channel_dynamic_sgs_and_moving_walls(
        &mut self,
        cx: &Cx<'_>,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<(DynamicSgsStepReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.prepare_channel_dynamic_sgs_with_moving_walls(
            cx, config, molecular, probe, velocities, workspace,
        )?;
        let walls = self.step_serial_with_moving_wall_momentum(
            cx, probe, velocities, workspace.collisions(),
        )?;
        Ok((report, walls))
    }

    /// Prepare the dynamic fit at F_n before advancing the force endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_channel_dynamic_sgs_and_moving_walls(
        &mut self,
        cx: &Cx<'_>,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<(DynamicSgsStepReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.prepare_channel_dynamic_sgs_with_moving_walls(
            cx, config, molecular, probe, velocities, workspace,
        )?;
        let walls = self.advance_serial_with_moving_wall_momentum(
            cx, probe, velocities, workspace.collisions(), forces_at_end,
        )?;
        Ok((report, walls))
    }

    /// Pooled dynamic channel SGS with the same held force and wall drive.
    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_channel_dynamic_sgs_and_moving_walls<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<(DynamicSgsStepReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.grid.prepare_channel_dynamic_sgs_with_motion(
            config, molecular, workspace, Some(self.forces()),
            Some(WallMotion3 { probe, velocities }), || gate_checkpoint(gate),
        )?;
        let walls = self.step_pooled_with_moving_wall_momentum(
            pool, gate, probe, velocities, workspace.collisions(),
        )?;
        Ok((report, walls))
    }

    /// Pooled dynamic wall step with atomic source-history publication.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_channel_dynamic_sgs_and_moving_walls<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<(DynamicSgsStepReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.grid.prepare_channel_dynamic_sgs_with_motion(
            config, molecular, workspace, Some(self.forces()),
            Some(WallMotion3 { probe, velocities }), || gate_checkpoint(gate),
        )?;
        let walls = self.advance_pooled_with_moving_wall_momentum(
            pool, gate, probe, velocities, workspace.collisions(), forces_at_end,
        )?;
        Ok((report, walls))
    }
}
