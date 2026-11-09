//! Local SGS at the original fixed, arbitrarily selected moving wall links.
//!
//! Unlike a channel test filter, a local closure needs no homogeneous plane.
//! Each missing Cartesian neighbor uses its actual axial wall-link velocity;
//! diagonal-only contacts do not invent an axial wall. Kinematics and SGS share
//! this mapping, and native moving bounce-back retains the complete link field.

use super::*;
use crate::lattice3::moving_walls::MovingWallMomentum3;
use crate::lattice3::sources::{ForcedGrid3, ForcedLattice3};

impl<L: Lattice3<Q>, const Q: usize> WallMotion3<'_, L, Q> {
    /// Sample admitted motion in increasing cell order, using fixed storage.
    /// The native probe is sorted and unique, so the complete sweep visits each
    /// selected link at most once. Unselected axial links stay stationary.
    pub(super) fn cartesian_velocities(
        &self,
    ) -> Result<impl FnMut(usize, usize, bool) -> [f64; 3] + '_, GridError3<L::Error>> {
        let mut axial_face = [6usize; Q];
        for face in 0..6 {
            let mut incoming = [0; 3];
            incoming[face / 2] = if face % 2 == 0 { 1 } else { -1 };
            let q = L::E
                .iter()
                .position(|&e| e == incoming)
                .ok_or(GridError3::Input(
                    "moving-wall gradients require Cartesian lattice links",
                ))?;
            axial_face[q] = face;
        }
        let mut cursor = 0;
        let mut cached_cell = usize::MAX;
        let mut walls = [[0.0; 3]; 6];
        Ok(move |cell, axis, high| {
            if cached_cell != cell {
                cached_cell = cell;
                walls = [[0.0; 3]; 6];
                while cursor < self.probe.links.len() && self.probe.links[cursor].cell <= cell {
                    let link = &self.probe.links[cursor];
                    let face = axial_face[link.incoming];
                    if link.cell == cell && face < 6 {
                        walls[face] = self.velocities[cursor];
                    }
                    cursor += 1;
                }
            }
            walls[2 * axis + usize::from(high)]
        })
    }
}

impl<L: SgsLattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Prepare local Smagorinsky, WALE or Vreman from the accepted physical
    /// velocity and the SAME wall-link drives used by native moving streaming.
    /// Interior solids, multiple wall pairs, spatially varying tangential
    /// motion and selected surface patches keep their original geometry.
    /// Unselected reflected links stay stationary; no plane average is taken.
    ///
    /// Native Mach, tangency, topology and local impermeability admission runs
    /// first. The shared full-gradient stencil uses each axial link's velocity
    /// at its halfway location. Uniform co-motion adds no spurious wall shear.
    /// The original molecular field, independent collision rates and ratio
    /// refusal remain authoritative. Use fresh molecular values, never a prior
    /// SGS-augmented collision field.
    ///
    /// Reuses the ordinary SGS workspace; there is no step allocation or
    /// per-cell probe search. This is an exploratory resolved-wall model, not
    /// a wall function, moving geometry/remap or turbulent-flow qualification.
    pub fn prepare_sgs_with_moving_walls(
        &self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        self.prepare_sgs_with_motion(
            config,
            molecular,
            workspace,
            None,
            WallMotion3 { probe, velocities },
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Advance the local closure and original moving-wall stream together.
    /// Wall impulse/work comes from the actual SGS-augmented populations.
    /// The SGS report describes incoming resolved/modeled shear transfer; it
    /// is not deposited heat or an integrated total-energy balance.
    /// Preparation, collision, stream and receipt failures publish no state.
    pub fn step_serial_with_sgs_and_moving_walls(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self
            .prepare_sgs_with_moving_walls(cx, config, molecular, probe, velocities, workspace)?;
        let walls = self.step_serial_with_moving_wall_momentum(
            cx,
            probe,
            velocities,
            Some(workspace.collisions()),
        )?;
        Ok((report, walls))
    }

    /// The same serial coefficient preparation and original pooled collision,
    /// moving stream and canonical wall receipt. Changing worker count does
    /// not change the closure arithmetic, accepted populations or wall loads.
    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_sgs_and_moving_walls<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.prepare_sgs_with_motion(
            config,
            molecular,
            workspace,
            None,
            WallMotion3 { probe, velocities },
            || {
                if gate.is_requested() {
                    Err(GridError3::Cancelled)
                } else {
                    Ok(())
                }
            },
        )?;
        let walls = self.step_pooled_with_moving_wall_momentum(
            pool,
            gate,
            probe,
            velocities,
            Some(workspace.collisions()),
        )?;
        Ok((report, walls))
    }

    pub(in crate::lattice3) fn prepare_sgs_with_motion(
        &self,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
        forces: Option<&[[f64; 3]]>,
        motion: WallMotion3<'_, L, Q>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        checkpoint()?;
        config.validate()?;
        motion.admit(self, &mut checkpoint)?;
        self.prepare_sgs_with_forces_and_walls(
            config,
            molecular,
            workspace,
            forces,
            motion.cartesian_velocities()?,
            checkpoint,
        )
    }
}

impl<L: SgsLattice3<Q> + ForcedLattice3<Q>, const Q: usize> ForcedGrid3<L, Q> {
    /// The same local wall-aware closure in the CURRENT spatial force frame.
    /// Accepted F_n/2 enters the full velocity gradient before fitting. A
    /// later endpoint force never changes this incoming wall/shear evaluation.
    pub fn prepare_sgs_with_moving_walls(
        &self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        self.grid.prepare_sgs_with_motion(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            WallMotion3 { probe, velocities },
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Hold the original spatial forces while advancing arbitrary fixed moving
    /// walls and local SGS in the native population/load transaction.
    pub fn step_serial_with_sgs_and_moving_walls(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self
            .prepare_sgs_with_moving_walls(cx, config, molecular, probe, velocities, workspace)?;
        let walls = self.step_serial_with_moving_wall_momentum(
            cx,
            probe,
            velocities,
            workspace.collisions(),
        )?;
        Ok((report, walls))
    }

    /// Prepare at F_n, then atomically publish the native moving-wall stream
    /// and F_(n+1). A refused endpoint preserves populations, forces and clock.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_serial_with_sgs_and_moving_walls(
        &mut self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
        forces_at_end: Vec<[f64; 3]>,
    ) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self
            .prepare_sgs_with_moving_walls(cx, config, molecular, probe, velocities, workspace)?;
        let walls = self.advance_serial_with_moving_wall_momentum(
            cx,
            probe,
            velocities,
            workspace.collisions(),
            forces_at_end,
        )?;
        Ok((report, walls))
    }

    /// Original pooled fixed-force step with the same wall-aware preparation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_sgs_and_moving_walls<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<(SgsReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let report = self.grid.prepare_sgs_with_motion(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            WallMotion3 { probe, velocities },
            || {
                if gate.is_requested() {
                    Err(GridError3::Cancelled)
                } else {
                    Ok(())
                }
            },
        )?;
        let walls = self.step_pooled_with_moving_wall_momentum(
            pool,
            gate,
            probe,
            velocities,
            workspace.collisions(),
        )?;
        Ok((report, walls))
    }

    /// Original pooled trapezoidal force advance and actual wall receipt.
    /// Preparation uses only the accepted force history; endpoint admission
    /// and joint publication belong to the unchanged native source owner.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_pooled_with_sgs_and_moving_walls<P: KernelRunner>(
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
        let report = self.grid.prepare_sgs_with_motion(
            config,
            molecular,
            workspace,
            Some(self.forces()),
            WallMotion3 { probe, velocities },
            || {
                if gate.is_requested() {
                    Err(GridError3::Cancelled)
                } else {
                    Ok(())
                }
            },
        )?;
        let walls = self.advance_pooled_with_moving_wall_momentum(
            pool,
            gate,
            probe,
            velocities,
            workspace.collisions(),
            forces_at_end,
        )?;
        Ok((report, walls))
    }
}

#[cfg(test)]
mod tests;
