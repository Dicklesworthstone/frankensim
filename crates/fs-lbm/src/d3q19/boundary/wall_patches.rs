//! Spatial partitions of the original post-collision wall momentum exchange.

use super::{
    BoundaryGrid3, BoundaryLink3, BoundarySgsError3, BoundaryStepError3, CollisionModel3,
    WallMomentumExchange3,
};
use crate::sgs3::{SgsConfig, SgsReport, SgsWorkspace3};
use core::mem::size_of;
use fs_exec::CancelGate;

/// A fixed partition of a canonical wall-link selection, with reusable loads.
///
/// Every selected link belongs to exactly one nonempty patch. The same original
/// wall reflection supplies force, moment, mass exchange and wall work; this
/// does not infer patch area, pressure, an acoustic source or unresolved stress.
/// Moments retain the native grid-coordinate convention, including its
/// periodic-winding limitation. Origins need not lie on their patches.
///
/// Accepted loads change only after the corresponding fluid step succeeds.
/// Scratch is rebuilt on retry. Measurement visits the selected links once,
/// caches post-collision density per fluid cell and allocates nothing per step.
#[derive(Debug)]
pub struct WallPatchProbe3 {
    dimensions: [usize; 3],
    links: Vec<BoundaryLink3>,
    patch_ids: Vec<usize>,
    origins: Vec<[f64; 3]>,
    accepted: Vec<WallMomentumExchange3>,
    pending: Vec<WallMomentumExchange3>,
    has_sample: bool,
    workspace_bytes: usize,
}

impl WallPatchProbe3 {
    /// Loads of the last accepted step, in the authored patch order.
    #[must_use]
    pub fn exchanges(&self) -> Option<&[WallMomentumExchange3]> {
        self.has_sample.then_some(self.accepted.as_slice())
    }

    #[must_use]
    pub fn moment_origins(&self) -> &[[f64; 3]] {
        &self.origins
    }

    /// Requested heap payload, excluding allocator overhead and the probe value.
    #[must_use]
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    pub(super) fn validate(&self, grid: &BoundaryGrid3) -> Result<(), BoundaryStepError3> {
        if grid.dimensions() != self.dimensions {
            return Err(refused(
                "wall patch probe belongs to different grid dimensions",
            ));
        }
        grid.validate_wall_measurement(&self.links, [0.0; 3])
    }

    pub(super) fn stage(
        &mut self,
        grid: &BoundaryGrid3,
        gate: &CancelGate,
    ) -> Result<(), BoundaryStepError3> {
        self.pending.fill(WallMomentumExchange3::default());
        let pending = &mut self.pending;
        let patch_ids = &self.patch_ids;
        let origins = &self.origins;
        grid.visit_wall_momentum(
            &self.links,
            Some(gate),
            |i, link, outgoing, rho, velocity| {
                let patch = patch_ids[i];
                pending[patch].add_link(link, outgoing, rho, velocity, origins[patch]);
            },
        )?;
        for exchange in &self.pending {
            let _ = exchange.checked()?;
        }
        Ok(())
    }

    pub(super) fn publish(&mut self) {
        core::mem::swap(&mut self.accepted, &mut self.pending);
        self.has_sample = true;
    }
}

impl BoundaryGrid3 {
    /// Prepare spatial loads from an authored partition of wall links.
    ///
    /// `links` retains the unique canonical order of [`Self::boundary_links`];
    /// `patch_ids[i]` selects one origin for `links[i]`. Every patch must contain
    /// a link. The budget admits all copied links, IDs, origins and both load
    /// buffers before allocation. Topology is rechecked when stepping, so a
    /// probe cannot silently measure a different wall after geometry changes.
    ///
    /// # Errors
    /// Original link/origin validation, empty or invalid partitions, budget or
    /// size overflow and allocation refusal return the native wall-load error.
    pub fn wall_patch_probe(
        &self,
        links: &[BoundaryLink3],
        patch_ids: &[usize],
        moment_origins: &[[f64; 3]],
        maximum_workspace_bytes: usize,
    ) -> Result<WallPatchProbe3, BoundaryStepError3> {
        if moment_origins.is_empty()
            || links.len() != patch_ids.len()
            || moment_origins.len() > links.len()
            || moment_origins.iter().flatten().any(|x| !x.is_finite())
            || patch_ids.iter().any(|&id| id >= moment_origins.len())
        {
            return Err(refused(
                "wall patches require finite origins and one valid ID per link",
            ));
        }
        self.validate_wall_measurement(links, [0.0; 3])?;
        let workspace_bytes = links
            .len()
            .checked_mul(size_of::<BoundaryLink3>() + size_of::<usize>())
            .and_then(|bytes| {
                moment_origins
                    .len()
                    .checked_mul(size_of::<[f64; 3]>() + 2 * size_of::<WallMomentumExchange3>())
                    .and_then(|patch_bytes| bytes.checked_add(patch_bytes))
            })
            .filter(|&bytes| bytes <= maximum_workspace_bytes)
            .ok_or_else(|| refused("wall patch workspace exceeds byte budget"))?;
        let mut pending: Vec<WallMomentumExchange3> = Vec::new();
        pending
            .try_reserve_exact(moment_origins.len())
            .map_err(|_| refused("wall patch workspace allocation refused"))?;
        pending.resize(moment_origins.len(), WallMomentumExchange3::default());
        for &id in patch_ids {
            pending[id].measured_links += 1;
        }
        if pending.iter().any(|exchange| exchange.measured_links == 0) {
            return Err(refused("every wall patch must contain a selected link"));
        }
        pending.fill(WallMomentumExchange3::default());
        Ok(WallPatchProbe3 {
            dimensions: self.dimensions(),
            links: copied(links)?,
            patch_ids: copied(patch_ids)?,
            origins: copied(moment_origins)?,
            accepted: copied(&pending)?,
            pending,
            has_sample: false,
            workspace_bytes,
        })
    }

    /// Advance once and retain every selected spatial wall load from that
    /// actual molecular collision. The existing boundary reconstruction and
    /// transactional streaming owner publish the fluid before the probe.
    /// No full-domain force scan or extra collision is performed.
    ///
    /// # Errors
    /// Invalid topology, cancellation, collision, nonfinite exchange or open
    /// boundary refusal leaves accepted populations and probe loads unchanged.
    pub fn step_with_wall_patches(
        &mut self,
        gate: &CancelGate,
        probe: &mut WallPatchProbe3,
    ) -> Result<(), BoundaryStepError3> {
        probe.validate(self)?;
        self.topology_locked = true;
        self.try_collide_with_gate(Some(gate))?;
        probe.stage(self, gate)?;
        self.stream_collided_serial(gate)?;
        probe.publish();
        Ok(())
    }

    /// Advance with an explicitly authored SGS closure and retain each patch's
    /// loads from the same SGS-augmented collision that advances the fluid.
    /// The closure, physical velocity gradients and refusal limits retain their
    /// native owners. No unresolved random force or acoustic term is added.
    ///
    /// # Errors
    /// The original SGS and wall-load refusals preserve accepted populations
    /// and all previously accepted patch loads.
    pub fn step_with_sgs_and_wall_patches(
        &mut self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SgsWorkspace3<CollisionModel3>,
        probe: &mut WallPatchProbe3,
    ) -> Result<SgsReport, BoundarySgsError3> {
        probe.validate(self)?;
        let report = self.prepare_sgs(gate, config, workspace)?;
        self.step_sgs_prepared_measured(gate, workspace, |grid| probe.stage(grid, gate))?;
        probe.publish();
        Ok(report)
    }
}

fn refused(reason: &'static str) -> BoundaryStepError3 {
    BoundaryStepError3::WallMomentum { reason }
}

fn copied<T: Copy>(values: &[T]) -> Result<Vec<T>, BoundaryStepError3> {
    let mut output = Vec::new();
    output
        .try_reserve_exact(values.len())
        .map_err(|_| refused("wall patch workspace allocation refused"))?;
    output.extend_from_slice(values);
    Ok(output)
}

#[cfg(test)]
mod tests;
