//! Material rheology, channel SGS and prescribed wall motion in one native step.
//!
//! Molecular viscosity is resolved from the incoming force-centered velocity and
//! the exact prescribed halfway-wall velocities. The selected channel SGS then
//! augments that molecular field using the same force frame and wall motion.
//! Finally one original native collision/stream applies the resulting rates and
//! returns the actual moving-wall momentum/work receipt. No second collision,
//! duplicate population buffer or persistent eddy viscosity is introduced.

use super::super::moving_walls::channel_walls;
use super::*;
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallMotion3};
use crate::lattice3::sources::{ForcedGrid3, ForcedLattice3};
use crate::lattice3::wall_loads::WallMomentumProbe3;
use crate::sgs3::dynamic::channel::ChannelDynamicSgsConfig;
use crate::sgs3::{ChannelDynamicSgsWorkspace3, SgsConfig, SgsWorkspace3};
use fs_exec::{CancelGate, Cx, KernelRunner};

/// Channel-only SGS choice for source-resolved moving-wall rheology.
pub enum ChannelRheologySgs3<'a, C> {
    Local {
        config: SgsConfig,
        workspace: &'a mut SgsWorkspace3<C>,
    },
    Dynamic {
        config: ChannelDynamicSgsConfig,
        workspace: &'a mut ChannelDynamicSgsWorkspace3<C>,
    },
}

fn prepare_channel_sgs<'a, L: SgsLattice3<Q>, const Q: usize>(
    grid: &Grid3<L, Q>,
    forces: Option<&[[f64; 3]]>,
    collisions: &mut [L::Collision],
    sgs: ChannelRheologySgs3<'a, L::Collision>,
    motion: WallMotion3<'_, L, Q>,
    mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
) -> Result<(RheologySgsReport3, &'a [L::Collision]), SgsStepError<L::Error>> {
    let (report, applied) = match sgs {
        ChannelRheologySgs3::Local { config, workspace } => {
            let report = grid.prepare_channel_sgs_with_motion(
                config,
                Some(&*collisions),
                workspace,
                forces,
                motion,
                &mut checkpoint,
            )?;
            (RheologySgsReport3::Local(report), workspace.collisions())
        }
        ChannelRheologySgs3::Dynamic { config, workspace } => {
            let report = grid.prepare_channel_dynamic_sgs_with_motion(
                config,
                Some(&*collisions),
                workspace,
                forces,
                Some(motion),
                &mut checkpoint,
            )?;
            (RheologySgsReport3::Dynamic(report), workspace.collisions())
        }
    };
    for (trial, field) in collisions.chunks_mut(64).zip(applied.chunks(64)) {
        checkpoint()?;
        trial.copy_from_slice(field);
    }
    checkpoint()?;
    Ok((report, applied))
}

impl<L: SgsLattice3<Q>, const Q: usize> Grid3<L, Q> {
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_channel_material_rheology_sgs_and_moving_walls<
        'w,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: ChannelRheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        (RheologySgsStep3<'w, L::Collision, F>, MovingWallMomentum3),
        RheologySgsError3<L::Error>,
    > {
        let motion = WallMotion3 { probe, velocities };
        let mut checkpoint = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let walls = channel_walls(self, &motion, &mut checkpoint)?;
        prepare_workspace_with_walls(
            self, reference, None, Some(walls), workspace, resolve, &mut checkpoint,
        )?;
        let (report, applied) = prepare_channel_sgs(
            self, None, &mut workspace.trial.collisions, sgs, motion, &mut checkpoint,
        )?;
        let wall =
            self.step_serial_with_moving_wall_momentum(cx, probe, velocities, Some(applied))?;
        Ok((publish(workspace, report), wall))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_channel_material_rheology_sgs_and_moving_walls<
        'w,
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: ChannelRheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        (RheologySgsStep3<'w, L::Collision, F>, MovingWallMomentum3),
        RheologySgsError3<L::Error>,
    > {
        let motion = WallMotion3 { probe, velocities };
        let mut checkpoint = || gate_poll(gate);
        let walls = channel_walls(self, &motion, &mut checkpoint)?;
        prepare_workspace_with_walls(
            self, reference, None, Some(walls), workspace, resolve, &mut checkpoint,
        )?;
        let (report, applied) = prepare_channel_sgs(
            self, None, &mut workspace.trial.collisions, sgs, motion, &mut checkpoint,
        )?;
        let wall = self.step_pooled_with_moving_wall_momentum(
            pool, gate, probe, velocities, Some(applied),
        )?;
        Ok((publish(workspace, report), wall))
    }
}

impl<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    #[allow(clippy::too_many_arguments)]
    pub fn step_serial_with_channel_material_rheology_sgs_and_moving_walls<
        'w,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: ChannelRheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        (RheologySgsStep3<'w, L::Collision, F>, MovingWallMomentum3),
        RheologySgsError3<L::Error>,
    > {
        let motion = WallMotion3 { probe, velocities };
        let mut checkpoint = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let walls = channel_walls(&self.grid, &motion, &mut checkpoint)?;
        prepare_workspace_with_walls(
            &self.grid, reference, Some(self.forces()), Some(walls), workspace, resolve,
            &mut checkpoint,
        )?;
        let (report, applied) = prepare_channel_sgs(
            &self.grid, Some(self.forces()), &mut workspace.trial.collisions, sgs, motion,
            &mut checkpoint,
        )?;
        let wall = self.step_serial_with_moving_wall_momentum(cx, probe, velocities, applied)?;
        Ok((publish(workspace, report), wall))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_channel_material_rheology_sgs_and_moving_walls<
        'w,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        cx: &Cx<'_>,
        reference: &MaterialScaling<'_, R>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        forces_at_end: Vec<[f64; 3]>,
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: ChannelRheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        (RheologySgsStep3<'w, L::Collision, F>, MovingWallMomentum3),
        RheologySgsError3<L::Error>,
    > {
        let motion = WallMotion3 { probe, velocities };
        let mut checkpoint = || cx.checkpoint().map_err(|_| GridError3::Cancelled);
        let walls = channel_walls(&self.grid, &motion, &mut checkpoint)?;
        prepare_workspace_with_walls(
            &self.grid, reference, Some(self.forces()), Some(walls), workspace, resolve,
            &mut checkpoint,
        )?;
        let (report, applied) = prepare_channel_sgs(
            &self.grid, Some(self.forces()), &mut workspace.trial.collisions, sgs, motion,
            &mut checkpoint,
        )?;
        let wall = self.advance_serial_with_moving_wall_momentum(
            cx, probe, velocities, applied, forces_at_end,
        )?;
        Ok((publish(workspace, report), wall))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_channel_material_rheology_sgs_and_moving_walls<
        'w,
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: ChannelRheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        (RheologySgsStep3<'w, L::Collision, F>, MovingWallMomentum3),
        RheologySgsError3<L::Error>,
    > {
        let motion = WallMotion3 { probe, velocities };
        let mut checkpoint = || gate_poll(gate);
        let walls = channel_walls(&self.grid, &motion, &mut checkpoint)?;
        prepare_workspace_with_walls(
            &self.grid, reference, Some(self.forces()), Some(walls), workspace, resolve,
            &mut checkpoint,
        )?;
        let (report, applied) = prepare_channel_sgs(
            &self.grid, Some(self.forces()), &mut workspace.trial.collisions, sgs, motion,
            &mut checkpoint,
        )?;
        let wall =
            self.step_pooled_with_moving_wall_momentum(pool, gate, probe, velocities, applied)?;
        Ok((publish(workspace, report), wall))
    }

    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_channel_material_rheology_sgs_and_moving_walls<
        'w,
        P: KernelRunner,
        R: FluidTransport + ?Sized,
        F: FluidTransport,
    >(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        reference: &MaterialScaling<'_, R>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        forces_at_end: Vec<[f64; 3]>,
        workspace: &'w mut RheologyWorkspace3<L::Collision, F>,
        sgs: ChannelRheologySgs3<'_, L::Collision>,
        resolve: impl FnMut(RheologyPoint3, &FluidTransportState) -> Result<F, MaterialError>,
    ) -> Result<
        (RheologySgsStep3<'w, L::Collision, F>, MovingWallMomentum3),
        RheologySgsError3<L::Error>,
    > {
        let motion = WallMotion3 { probe, velocities };
        let mut checkpoint = || gate_poll(gate);
        let walls = channel_walls(&self.grid, &motion, &mut checkpoint)?;
        prepare_workspace_with_walls(
            &self.grid, reference, Some(self.forces()), Some(walls), workspace, resolve,
            &mut checkpoint,
        )?;
        let (report, applied) = prepare_channel_sgs(
            &self.grid, Some(self.forces()), &mut workspace.trial.collisions, sgs, motion,
            &mut checkpoint,
        )?;
        let wall = self.advance_pooled_with_moving_wall_momentum(
            pool, gate, probe, velocities, applied, forces_at_end,
        )?;
        Ok((publish(workspace, report), wall))
    }
}
