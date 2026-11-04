//! Local and dynamic plane-channel viscosity on the original native lattice step owner.
//!
//! The unique nonperiodic axis defines two fixed halfway walls; both
//! tangents remain periodic and homogeneous. No solid or wall cell enters the
//! test filter. The shared model supplies only shear viscosity; D3Q27 retains
//! its independent trace and all higher central-moment/cumulant rates.
//! The explicit moving-wall entries bind the same uniform tangential wall
//! velocities to the filter, strain diagnostics and original native stream.

use super::*;
use crate::lattice3::moving_walls::{MovingWallMomentum3, WallMotion3};
use crate::lattice3::wall_loads::WallMomentumProbe3;
use crate::sgs3::dynamic::DynamicSgsError;
use crate::sgs3::dynamic::channel::{ChannelDynamicSgs, ChannelDynamicSgsConfig};
use crate::sgs3::{ChannelDynamicSgsWorkspace3, DynamicSgsStepReport};

impl<L: SgsLattice3<Q>, const Q: usize> Grid3<L, Q> {
    fn admit_channel_sgs_topology(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<usize, SgsStepError<L::Error>> {
        let mut normal = None;
        for (axis, &periodic) in self.periodic.iter().enumerate() {
            if !periodic && normal.replace(axis).is_some() {
                return Err(DynamicSgsError::Invalid(
                    "channel SGS requires exactly one fixed wall pair and two periodic axes",
                )
                .into());
            }
        }
        let normal = normal.ok_or(DynamicSgsError::Invalid(
            "channel SGS requires exactly one fixed wall pair and two periodic axes",
        ))?;
        for tile in self.solid.chunks(64) {
            checkpoint()?;
            if tile.iter().any(|&solid| solid) {
                return Err(DynamicSgsError::Invalid(
                    "channel SGS requires an all-fluid rectangular region between its external walls",
                ).into());
            }
        }
        Ok(normal)
    }

    // A native probe is an immutable, duplicate-free subset of reflected
    // links. On this admitted all-fluid rectangle, its exact full link count
    // therefore proves that both physical walls are represented completely.
    // No per-cell lookup or additional workspace is needed for uniform planes.
    fn channel_sgs_wall_velocities(
        &self,
        normal: usize,
        motion: &WallMotion3<'_, L, Q>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<[[f64; 3]; 2], SgsStepError<L::Error>> {
        motion.admit(self, &mut checkpoint)?;
        let layers = L::E.iter().try_fold(0usize, |sum, direction| {
            sum.checked_add((direction[normal].unsigned_abs() as usize).min(self.dims[normal]))
                .ok_or(DynamicSgsError::Invalid("channel wall link count overflow"))
        })?;
        let required = (self.f.len() / self.dims[normal])
            .checked_mul(layers)
            .ok_or(DynamicSgsError::Invalid("channel wall link count overflow"))?;
        if motion.probe.measured_links() != required {
            return Err(DynamicSgsError::Invalid(
                "moving channel SGS requires the complete native links of both walls",
            )
            .into());
        }
        let mut walls = [None; 2];
        for (index, (link, &velocity)) in
            motion.probe.links.iter().zip(motion.velocities).enumerate()
        {
            if index % 64 == 0 {
                checkpoint()?;
            }
            let side = usize::from(L::E[link.incoming][normal] < 0);
            if let Some(expected) = walls[side] {
                if velocity != expected {
                    return Err(DynamicSgsError::Invalid(
                        "channel SGS requires a uniform velocity on each wall plane",
                    )
                    .into());
                }
            } else {
                walls[side] = Some(velocity);
            }
        }
        checkpoint()?;
        let [Some(lower), Some(upper)] = walls else {
            return Err(DynamicSgsError::Invalid(
                "moving channel SGS requires both native wall planes",
            )
            .into());
        };
        Ok([lower, upper])
    }

    /// Prepare local Smagorinsky, WALE or Vreman transport with the SAME wall
    /// motion used by native streaming. The domain has exactly one nonperiodic
    /// axis, two periodic tangents and no interior solids. The probe must contain
    /// every reflected link on both walls, including a stationary partner.
    /// Each wall has one uniform tangential velocity in native probe order.
    ///
    /// Current force-centered velocity and the physical halfway-wall values
    /// supply the full gradient, local closure and shear-dissipation account.
    /// The original molecular field and independent collision rates remain
    /// authoritative. No test filter, extra workspace or step allocation is
    /// added; the ordinary SGS workspace is reused. Zero drives preserve the
    /// stationary arithmetic. This is an exploratory resolved-wall closure,
    /// without a wall function, unresolved-energy model or acceptance claim.
    pub fn prepare_channel_sgs_with_moving_walls(
        &self,
        cx: &Cx<'_>,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut SgsWorkspace3<L::Collision>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        self.prepare_channel_sgs_with_motion(
            config,
            molecular,
            workspace,
            WallMotion3 { probe, velocities },
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Advance the local closure and prescribed walls in one original native
    /// transaction. The returned wall impulse/work measures the actual stream;
    /// the SGS report is its pre-step shear diagnostic, not integrated heat.
    /// Preparation, collision, endpoint and receipt refusals leave accepted
    /// populations, force and clock unchanged.
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
            cx,
            probe,
            velocities,
            Some(workspace.collisions()),
        )?;
        Ok((report, walls))
    }

    /// Identical local preparation followed by the original pooled moving-wall
    /// owner. Canonical probe and cell order retain worker-count determinism;
    /// the collision augmentation belongs only to this step.
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
        let report = self.prepare_channel_sgs_with_motion(
            config,
            molecular,
            workspace,
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

    fn prepare_channel_sgs_with_motion(
        &self,
        config: SgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut SgsWorkspace3<L::Collision>,
        motion: WallMotion3<'_, L, Q>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<SgsReport, SgsStepError<L::Error>> {
        checkpoint()?;
        config.validate()?;
        let normal = self.admit_channel_sgs_topology(&mut checkpoint)?;
        let walls = self.channel_sgs_wall_velocities(normal, &motion, &mut checkpoint)?;
        self.prepare_sgs_with_forces_and_walls(
            config,
            molecular,
            workspace,
            None,
            |_, _, positive| walls[usize::from(positive)],
            checkpoint,
        )
    }

    /// Allocate the complete reusable macroscopic, collision and channel-filter
    /// payload. The wall normal is the grid's unique nonperiodic axis; at least
    /// two normal cells and three per periodic tangent are required. No solid
    /// cells are admitted. The byte gate covers predicted and actual numerical
    /// capacities, excluding fixed metadata and the grid's original buffers.
    pub fn channel_dynamic_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<ChannelDynamicSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        let normal = self.admit_channel_sgs_topology(|| Ok(()))?;
        let cells = self.f.len();
        let local = cells
            .checked_mul(4 * core::mem::size_of::<f64>() + core::mem::size_of::<L::Collision>())
            .ok_or(DynamicSgsError::Invalid("channel workspace size overflow"))?;
        let required = local
            .checked_add(ChannelDynamicSgs::required_bytes(cells)?)
            .ok_or(DynamicSgsError::Invalid("channel workspace size overflow"))?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            }
            .into());
        }
        let filter = ChannelDynamicSgs::new(self.dims, 1.0, normal, maximum_bytes - local)?;
        let result = ChannelDynamicSgsWorkspace3 {
            local: self.sgs_workspace()?,
            filter,
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

    /// Fit each homogeneous wall-parallel plane from current Guo-centered
    /// velocities and the actual stationary halfway walls. Molecular rates
    /// come from the supplied original field or constructor collision law;
    /// previously augmented rates must not be resubmitted as molecular input.
    ///
    /// Preparation allocates nothing and polls every 64 cells, including shared
    /// filtering. It changes no populations, force, accepted clock or molecular
    /// law. A zero fit preserves every original collision bit exactly. This is
    /// exploratory, resolved-wall low-Mach LES, not high-Re/DNS qualification.
    pub fn prepare_channel_dynamic_sgs(
        &self,
        cx: &Cx<'_>,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        self.prepare_channel_dynamic_sgs_impl(config, molecular, workspace, || {
            cx.checkpoint().map_err(|_| GridError3::Cancelled)
        })
    }

    /// Fit dynamic SGS using the same prescribed wall motion as native
    /// moving bounce-back. The probe must select EVERY reflected link on both
    /// channel walls, including the stationary wall when only one is driven.
    /// `velocities` is in probe order, with a spatially uniform tangential
    /// vector on each wall. The two wall values are derived from that exact
    /// native drive; partial probes and nonuniform plane velocities refuse.
    ///
    /// The existing channel workspace and probe are reused without allocation.
    /// Current Guo-centered flow, both wall velocities and the original local
    /// molecular law determine the fit and dissipation diagnostics together.
    /// This prepares scratch only; publication still belongs to native stepping.
    pub fn prepare_channel_dynamic_sgs_with_moving_walls(
        &self,
        cx: &Cx<'_>,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        probe: &WallMomentumProbe3<L, Q>,
        velocities: &[[f64; 3]],
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        self.prepare_channel_dynamic_sgs_with_motion(
            config,
            molecular,
            workspace,
            None,
            Some(WallMotion3 { probe, velocities }),
            || cx.checkpoint().map_err(|_| GridError3::Cancelled),
        )
    }

    /// Advance dynamic channel SGS and return the actual moving-wall impulse
    /// and wall work from the original corrected stream. The fit, wall drive
    /// and collision field share this step; none persists into a later
    /// ordinary step. Any preparation, stream or receipt refusal publishes
    /// neither populations nor the accepted clock.
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
            cx,
            probe,
            velocities,
            Some(workspace.collisions()),
        )?;
        Ok((report, walls))
    }

    /// Deterministic preparation followed by the original pooled moving-wall
    /// transaction. Wall loads use canonical probe order, independently of
    /// worker count. The complete probe and uniform-plane rules are identical
    /// to [`Self::prepare_channel_dynamic_sgs_with_moving_walls`].
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
        let report = self.prepare_channel_dynamic_sgs_with_motion(
            config,
            molecular,
            workspace,
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
            Some(workspace.collisions()),
        )?;
        Ok((report, walls))
    }

    /// The original serial collision/stream transaction applies the prepared
    /// field for one step. Every endpoint Mach/collision check and cancellation
    /// gate precedes the original population/clock publication.
    pub fn step_serial_with_channel_dynamic_sgs(
        &mut self,
        cx: &Cx<'_>,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.prepare_channel_dynamic_sgs(cx, config, molecular, workspace)?;
        self.step_serial_with_collisions(cx, workspace.collisions())?;
        Ok(report)
    }

    /// Identical deterministic preparation followed by the original pooled
    /// collision/stream owner. Failed scratch is rebuilt on retry, and the
    /// constructor collision remains unchanged for later ordinary steps.
    pub fn step_pooled_with_channel_dynamic_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        let report = self.prepare_channel_dynamic_sgs_impl(config, molecular, workspace, || {
            if gate.is_requested() {
                Err(GridError3::Cancelled)
            } else {
                Ok(())
            }
        })?;
        self.step_pooled_with_collisions(pool, gate, workspace.collisions())?;
        Ok(report)
    }

    fn prepare_channel_dynamic_sgs_impl(
        &self,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
        checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        self.prepare_channel_dynamic_sgs_with_forces(config, molecular, workspace, None, checkpoint)
    }

    pub(super) fn prepare_channel_dynamic_sgs_with_forces(
        &self,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
        forces: Option<&[[f64; 3]]>,
        checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        self.prepare_channel_dynamic_sgs_with_motion(
            config, molecular, workspace, forces, None, checkpoint,
        )
    }

    fn prepare_channel_dynamic_sgs_with_motion(
        &self,
        config: ChannelDynamicSgsConfig,
        molecular: Option<&[L::Collision]>,
        workspace: &mut ChannelDynamicSgsWorkspace3<L::Collision>,
        forces: Option<&[[f64; 3]]>,
        motion: Option<WallMotion3<'_, L, Q>>,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        checkpoint()?;
        let normal = self.admit_channel_sgs_topology(&mut checkpoint)?;
        config.validate()?;
        let walls = match motion {
            Some(motion) => self.channel_sgs_wall_velocities(normal, &motion, &mut checkpoint)?,
            None => [[0.0; 3]; 2],
        };
        let cells = self.f.len();
        let local = &mut workspace.local;
        if local.velocity.len() != cells
            || local.density.len() != cells
            || local.collisions.len() != cells
            || workspace.filter.shape() != self.dims
            || workspace.filter.normal_axis() != normal
            || molecular.is_some_and(|field| field.len() != cells)
            || forces.is_some_and(|field| field.len() != cells)
        {
            return Err(
                GridError3::Input("channel SGS workspace or molecular field mismatch").into(),
            );
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
        let fit =
            workspace
                .filter
                .prepare(config, &local.density, &local.velocity, walls, || {
                    checkpoint().is_err()
                })?;
        let dissipation = self.apply_dynamic_sgs_viscosities_with_walls(
            config.maximum_eddy_to_molecular,
            molecular,
            workspace.filter.viscosities(),
            local,
            |_, _, positive| walls[usize::from(positive)],
            checkpoint,
        )?;
        Ok(DynamicSgsStepReport { fit, dissipation })
    }
}

#[cfg(test)]
mod local_tests;
#[cfg(test)]
mod tests;
