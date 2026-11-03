//! Plane-channel dynamic viscosity on the original native lattice step owner.
//!
//! The unique nonperiodic axis defines two stationary halfway walls; both
//! tangents remain periodic and homogeneous. No solid or wall cell enters the
//! test filter. The shared model supplies only shear viscosity; D3Q27 retains
//! its independent trace and all higher central-moment/cumulant rates.

use super::*;
use crate::sgs3::dynamic::DynamicSgsError;
use crate::sgs3::dynamic::channel::{ChannelDynamicSgs, ChannelDynamicSgsConfig};
use crate::sgs3::{ChannelDynamicSgsWorkspace3, DynamicSgsStepReport};

impl<L: SgsLattice3<Q>, const Q: usize> Grid3<L, Q> {
    fn admit_channel_dynamic_topology(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<usize, SgsStepError<L::Error>> {
        let mut normal = None;
        for (axis, &periodic) in self.periodic.iter().enumerate() {
            if !periodic && normal.replace(axis).is_some() {
                return Err(DynamicSgsError::Invalid(
                    "channel SGS requires exactly one stationary wall pair and two periodic axes",
                )
                .into());
            }
        }
        let normal = normal.ok_or(DynamicSgsError::Invalid(
            "channel SGS requires exactly one stationary wall pair and two periodic axes",
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

    /// Allocate the complete reusable macroscopic, collision and channel-filter
    /// payload. The wall normal is the grid's unique nonperiodic axis; at least
    /// two normal cells and three per periodic tangent are required. No solid
    /// cells are admitted. The byte gate covers predicted and actual numerical
    /// capacities, excluding fixed metadata and the grid's original buffers.
    pub fn channel_dynamic_sgs_workspace(
        &self,
        maximum_bytes: usize,
    ) -> Result<ChannelDynamicSgsWorkspace3<L::Collision>, SgsStepError<L::Error>> {
        let normal = self.admit_channel_dynamic_topology(|| Ok(()))?;
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
        mut checkpoint: impl FnMut() -> Result<(), GridError3<L::Error>>,
    ) -> Result<DynamicSgsStepReport, SgsStepError<L::Error>> {
        checkpoint()?;
        let normal = self.admit_channel_dynamic_topology(&mut checkpoint)?;
        config.validate()?;
        let cells = self.f.len();
        let local = &mut workspace.local;
        if local.velocity.len() != cells
            || local.density.len() != cells
            || local.collisions.len() != cells
            || workspace.filter.shape() != self.dims
            || workspace.filter.normal_axis() != normal
            || molecular.is_some_and(|field| field.len() != cells)
        {
            return Err(
                GridError3::Input("channel SGS workspace or molecular field mismatch").into(),
            );
        }
        for cell in 0..cells {
            if cell % 64 == 0 {
                checkpoint()?;
            }
            let (density, velocity) =
                macro_cell::<L, Q>(&self.f[cell], self.force, cell, self.budget.maximum_mach)?;
            local.density[cell] = density;
            local.velocity[cell] = velocity;
        }
        let fit = workspace.filter.prepare(
            config,
            &local.density,
            &local.velocity,
            [[0.0; 3]; 2],
            || checkpoint().is_err(),
        )?;
        let dissipation = self.apply_dynamic_sgs_viscosities(
            config.maximum_eddy_to_molecular,
            molecular,
            workspace.filter.viscosities(),
            local,
            checkpoint,
        )?;
        Ok(DynamicSgsStepReport { fit, dissipation })
    }
}

#[cfg(test)]
mod tests;
