//! Accepted Lagrangian turbulence memory on the original lattice transaction.

use super::*;
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallMotion3};
use crate::lattice3::wall_loads::WallMomentumProbe3;
use crate::sgs3::dynamic::DynamicSgsError;
use crate::sgs3::dynamic::lagrangian::{
    LagrangianDynamicSgs, LagrangianSgsConfig, LagrangianSgsTrial,
};
use crate::sgs3::{DynamicSgsStepReport, LagrangianSgsWorkspace3};

impl<L: SgsLattice3<Q>, const Q: usize> Grid3<L, Q> {
    /// Allocate a trajectory's local collision fields and two turbulence-memory
    /// banks. All three axes must be periodic, with no solids. The declared
    /// byte bound covers actual numerical capacities; stepping adds no arrays.
    pub fn lagrangian_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<LagrangianSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        self.admit_dynamic_topology(|| Ok(()))?;
        self.new_lagrangian_sgs_workspace(None, maximum_bytes)
    }

    /// Allocate trajectory memory for a resolved planar channel. Exactly one
    /// axis has fixed halfway walls; both tangents are periodic and all cells
    /// are fluid. The workspace binds that normal axis and cannot be reused on
    /// a periodic grid or a differently oriented channel. Stepping reuses its
    /// two history banks and local collision fields without adding arrays.
    pub fn channel_lagrangian_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<LagrangianSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        let normal = self.admit_channel_sgs_topology(|| Ok(()))?;
        self.new_lagrangian_sgs_workspace(Some(normal), maximum_bytes)
    }

    fn new_lagrangian_sgs_workspace(
        &self,
        normal_axis: Option<usize>,
        maximum_bytes: usize,
    ) -> Result<LagrangianSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        let cells = self.f.len();
        let local_bytes = cells
            .checked_mul(4 * core::mem::size_of::<f64>() + core::mem::size_of::<L::Collision>())
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian native workspace size overflow",
            ))?;
        let required = local_bytes
            .checked_add(LagrangianDynamicSgs::required_bytes(cells)?)
            .ok_or(DynamicSgsError::Invalid(
                "Lagrangian native workspace size overflow",
            ))?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            }
            .into());
        }
        let result = LagrangianSgsWorkspace3 {
            local: self.sgs_workspace()?,
            filter: match normal_axis {
                Some(normal) => LagrangianDynamicSgs::new_channel(
                    self.dims,
                    [1.0; 3],
                    normal,
                    maximum_bytes - local_bytes,
                )?,
                None => {
                    LagrangianDynamicSgs::new(self.dims, [1.0; 3], maximum_bytes - local_bytes)?
                }
            },
        };
        let required = result.workspace_bytes();
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            }
            .into());
        }
        Ok(result)
    }

    /// Advance the native D3Q19/D3Q27 flow and its Lagrangian dynamic SGS history
    /// in one transaction. Memory time is in lattice timesteps. Fits use the
    /// accepted Guo-centered physical velocity; the constructor's molecular
    /// collision or the explicitly supplied molecular field remains unchanged.
    /// D3Q27 retains its independent trace and higher collision rates.
    ///
    /// A closure, collision, Mach or cancellation refusal publishes neither
    /// populations, step count nor history. The fit clock is the incoming
    /// `completed_steps`, so failed trials can be retried with identical memory.
    /// The workspace admits either periodic flow or a channel with stationary
    /// halfway walls. Channel filtering and pathline interpolation never wrap
    /// through a wall. This is exploratory resolved-wall LES, not a DNS claim.
    pub fn step_serial_with_lagrangian_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let (trial, dissipation) = self.prepare_lagrangian_sgs_with_forces(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            None,
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        self.step_serial_with_collisions(cx, workspace.local.collisions())?;
        let fit = trial.report();
        trial.commit();
        Ok(DynamicSgsStepReport { fit, dissipation })
    }

    /// Identical deterministic fit followed by the original pooled step.
    /// Accepted history is swapped only after that native owner publishes.
    pub fn step_pooled_with_lagrangian_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let (trial, dissipation) = self.prepare_lagrangian_sgs_with_forces(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            None,
            || {
                if gate.is_requested() {
                    Err(GridError3::Cancelled)
                } else {
                    Ok(())
                }
            },
        )?;
        self.step_pooled_with_collisions(pool, gate, workspace.local.collisions())?;
        let fit = trial.report();
        trial.commit();
        Ok(DynamicSgsStepReport { fit, dissipation })
    }

    /// Advance channel turbulence memory with the actual native moving-wall
    /// impulse and work transaction. The probe must include every reflected
    /// link on both walls; each plane has one uniform tangential velocity.
    /// Those same velocities define the fit, pathline boundary and dissipation.
    /// Any preparation, collision, stream or receipt refusal leaves accepted
    /// populations, clock and history unchanged.
    pub fn step_serial_with_channel_lagrangian_sgs_and_moving_walls(
        &mut self,
        cx: &Cx<'_>,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
    ) -> Result<(DynamicSgsStepReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let (trial, dissipation) = self.prepare_lagrangian_sgs_with_motion(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            None,
            Some(WallMotion3 { probe, velocities }),
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )?;
        let walls = self.step_serial_with_moving_wall_momentum(
            cx,
            probe,
            velocities,
            Some(workspace.local.collisions()),
        )?;
        let fit = trial.report();
        trial.commit();
        Ok((DynamicSgsStepReport { fit, dissipation }, walls))
    }

    /// Deterministic channel fit followed by the original pooled moving-wall
    /// transaction. Turbulence memory commits only after populations, the
    /// lattice clock and the canonical-order wall receipt have all succeeded.
    #[allow(clippy::too_many_arguments)]
    pub fn step_pooled_with_channel_lagrangian_sgs_and_moving_walls<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut LagrangianSgsWorkspace3<L::Collision>,
    ) -> Result<(DynamicSgsStepReport, MovingWallMomentum3), SgsStepError<L::Error>> {
        let (trial, dissipation) = self.prepare_lagrangian_sgs_with_motion(
            config,
            molecular,
            &mut workspace.local,
            &mut workspace.filter,
            None,
            Some(WallMotion3 { probe, velocities }),
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
            Some(workspace.local.collisions()),
        )?;
        let fit = trial.report();
        trial.commit();
        Ok((DynamicSgsStepReport { fit, dissipation }, walls))
    }

    // Keep the fit, departure velocity and dissipation in the same accepted
    // force frame. The endpoint force belongs to the native step alone.
    pub(super) fn prepare_lagrangian_sgs_with_forces<'a>(
        &self,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        local: &mut SgsWorkspace3<L::Collision>,
        filter: &'a mut LagrangianDynamicSgs,
        forces: Option<&[[f64; 3]]>,
        checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(LagrangianSgsTrial<'a>, SgsReport), SgsStepError<L::Error>> {
        self.prepare_lagrangian_sgs_with_motion(
            config, molecular, local, filter, forces, None, checkpoint,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::lattice3) fn prepare_lagrangian_sgs_with_motion<'a>(
        &self,
        config: LagrangianSgsConfig,
        molecular: Option<&[L::Collision]>,
        local: &mut SgsWorkspace3<L::Collision>,
        filter: &'a mut LagrangianDynamicSgs,
        forces: Option<&[[f64; 3]]>,
        motion: Option<WallMotion3<'_, L, Q>>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<(LagrangianSgsTrial<'a>, SgsReport), SgsStepError<L::Error>> {
        checkpoint()?;
        let walls = match filter.normal_axis() {
            Some(expected) => {
                let normal = self.admit_channel_sgs_topology(&mut checkpoint)?;
                if normal != expected {
                    return Err(GridError3::Input(
                        "Lagrangian SGS workspace wall-normal axis mismatch",
                    )
                    .into());
                }
                Some(match motion {
                    Some(motion) => {
                        self.channel_sgs_wall_velocities(normal, &motion, &mut checkpoint)?
                    }
                    None => [[0.0; 3]; 2],
                })
            }
            None => {
                self.admit_dynamic_topology(&mut checkpoint)?;
                if motion.is_some() {
                    return Err(DynamicSgsError::Invalid(
                        "moving-wall Lagrangian SGS requires a channel workspace",
                    )
                    .into());
                }
                None
            }
        };
        config.validate()?;
        let cells = self.f.len();
        if local.velocity.len() != cells
            || local.density.len() != cells
            || local.collisions.len() != cells
            || filter.shape() != self.dims
            || molecular.is_some_and(|field| field.len() != cells)
            || forces.is_some_and(|field| field.len() != cells)
        {
            return Err(GridError3::Input(
                "Lagrangian SGS workspace, molecular or force field mismatch",
            )
            .into());
        }
        // The public clock is integral. A rounded floating clock must never
        // make two distinct accepted steps appear to be the same fit time.
        if self.completed_steps() as u128 > (1u128 << 53) {
            return Err(DynamicSgsError::Invalid(
                "Lagrangian lattice clock exceeds exact f64 range",
            )
            .into());
        }
        for cell in 0..cells {
            if cell % 64 == 0 {
                checkpoint()?;
            }
            let force = forces.map_or(self.force, |field| field[cell]);
            let (density, velocity) =
                macro_cell::<L, Q>(&self.f[cell], force, cell, self.budget.maximum_mach)?;
            local.density[cell] = density;
            local.velocity[cell] = velocity;
        }
        let time = self.completed_steps() as f64;
        let trial = match walls {
            Some(walls) => filter.prepare_with_walls(
                config,
                time,
                &local.density,
                &local.velocity,
                walls,
                || checkpoint().is_err(),
            )?,
            None => filter.prepare(config, time, &local.density, &local.velocity, || {
                checkpoint().is_err()
            })?,
        };
        let dissipation = match walls {
            Some(walls) => self.apply_dynamic_sgs_viscosities_with_walls(
                config.maximum_eddy_to_molecular,
                molecular,
                trial.viscosities(),
                local,
                |_, _, positive| walls[usize::from(positive)],
                checkpoint,
            )?,
            None => self.apply_dynamic_sgs_viscosities(
                config.maximum_eddy_to_molecular,
                molecular,
                trial.viscosities(),
                local,
                checkpoint,
            )?,
        };
        Ok((trial, dissipation))
    }
}

#[cfg(test)]
mod channel_tests;
