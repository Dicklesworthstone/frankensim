//! Incoming resolved-shear material rheology with native sparse SGS.
//!
//! Evaluating a nonlinear viscosity on the resolved strain is an explicit LES
//! approximation: filtering the nonlinear molecular stress need not equal
//! evaluating that stress at the filtered strain. The SGS closure models eddy
//! transport only; it supplies no closure for that unresolved rheological
//! correlation. No constitutive-memory or implicit flow/source solve is added.

use super::{SparseMaterialSgsError3, checkpoint};
use crate::d3q19::sparse::sgs::{SparseSgsError3, SparseSgsWorkspace3};
use crate::d3q19::sparse::{SparseError3, SparseGrid3, TILE_CELLS, demorton3};
use crate::d3q19::{CollisionModel3, TILE};
use crate::lattice3::{D3q19, ViscousLattice3};
use crate::material::immersed::{MaterialImmersedPoint3, MaterialImmersedSample3, strain_rate};
use crate::material::{MaterialLbmError, MaterialScaling, plan_momentum_scaling};
use crate::sgs3::{SgsConfig, SgsReport};
use fs_exec::{CancelGate, KernelRunner};
use fs_material::{MaterialError, fluid_transport::FluidTransport};

/// Reusable sparse gradient/collision scratch and two original-source epochs.
/// The applied epoch changes only after a successful native flow publication.
/// Construct again after changing topology, physical scaling or collision law.
/// Neither a previous SGS field nor an earlier source viscosity is fed back
/// as the next step's molecular constitutive response.
pub struct SparseRheologySgsWorkspace3<M> {
    flow: SparseSgsWorkspace3,
    scaling: [u64; 3],
    model: CollisionModel3,
    applied: Vec<MaterialImmersedSample3<M>>,
    pending: Vec<MaterialImmersedSample3<M>>,
    report: Option<SgsReport>,
}

impl<M> SparseRheologySgsWorkspace3<M> {
    /// Original responses and exact incoming shear/transport snapshots from
    /// the last successful step, in Morton-tile then x-fastest lane order.
    /// Empty until the first success; a refused step leaves this epoch intact.
    pub fn applied(&self) -> &[MaterialImmersedSample3<M>] {
        &self.applied
    }

    /// Diagnostics from the same accepted step as `applied()`.
    pub fn last_report(&self) -> Option<SgsReport> {
        self.report
    }

    /// Retained buffer capacities, excluding allocator/struct overhead and
    /// any heap storage owned by the original material responses themselves.
    pub fn workspace_bytes(&self) -> usize {
        self.flow.workspace_bytes()
            + (self.applied.capacity() + self.pending.capacity())
                * core::mem::size_of::<MaterialImmersedSample3<M>>()
    }

    fn accept(&mut self, report: SgsReport) {
        core::mem::swap(&mut self.applied, &mut self.pending);
        self.report = Some(report);
    }
}

impl<R: FluidTransport + ?Sized> MaterialScaling<'_, R> {
    /// Allocate active-cell scratch once, before sampling any source. The
    /// explicit budget bounds both retained source epochs and the flow fields.
    /// Borrowed source responses and `&dyn FluidTransport` remain supported.
    pub fn sparse_rheology_sgs_workspace<M>(
        &self,
        grid: &SparseGrid3,
        max_active_cells: usize,
        gate: &CancelGate,
    ) -> Result<SparseRheologySgsWorkspace3<M>, SparseMaterialSgsError3> {
        checkpoint(gate)?;
        let count = grid
            .keys
            .len()
            .checked_mul(TILE_CELLS)
            .filter(|&n| n <= max_active_cells)
            .ok_or(SparseMaterialSgsError3::Invalid(
                "active rheology cell budget exceeded",
            ))?;
        let flow = grid.sgs_workspace(gate)?;
        let mut applied = Vec::new();
        let mut pending = Vec::new();
        applied
            .try_reserve_exact(count)
            .map_err(|_| SparseSgsError3::Allocation)?;
        pending
            .try_reserve_exact(count)
            .map_err(|_| SparseSgsError3::Allocation)?;
        checkpoint(gate)?;
        Ok(SparseRheologySgsWorkspace3 {
            flow,
            scaling: self.sparse_rheology_scaling(),
            model: grid.model,
            applied,
            pending,
            report: None,
        })
    }

    /// Advance the sparse native owner using a fresh original material
    /// response at each cell's current resolved shear, plus the selected SGS.
    ///
    /// The force-centered velocity and stationary halfway-wall stencil are
    /// shared with sparse SGS. One full gradient supplies BOTH source shear
    /// `sqrt(2 D:D) / dt` [1/s] and eddy viscosity; rigid rotation contributes
    /// no source shear. The physical velocity and distance scales cancel dx.
    /// Each active cell is queried once, and its returned `FluidTransport` is
    /// read once. The existing material scaling converts its mu/rho0 into the
    /// molecular tau; the local SGS/molecular ratio uses that same response.
    /// Zero SGS still applies the full heterogeneous molecular collision field.
    ///
    /// This uses incoming-time resolved-shear generalized-Newtonian rheology.
    /// Nonlinear filtered molecular stress requires additional correlations
    /// that this approximation does not close. Source temperature/composition
    /// stay caller-owned, reference density must remain identical, and Cp/k
    /// impose no thermal lattice bound. Maintain the scaling's scenario speed
    /// bound. Native family and non-shear collision controls remain unchanged.
    ///
    /// Refusal leaves populations, step count and `workspace.applied()` intact.
    /// Callbacks must not publish external constitutive history: their side
    /// effects cannot be rolled back. Adapter buffers allocate only at workspace
    /// construction; original response-owned storage remains source-owned.
    pub fn step_sparse_with_material_rheology_and_sgs<M: FluidTransport>(
        &self,
        grid: &mut SparseGrid3,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseRheologySgsWorkspace3<M>,
        resolve: impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<SgsReport, SparseMaterialSgsError3> {
        let report = self.prepare_sparse_rheology_sgs(grid, gate, config, workspace, resolve)?;
        grid.step_serial_inner(Some(gate), Some(workspace.flow.collisions()))
            .map_err(SparseSgsError3::from)?;
        workspace.accept(report);
        Ok(report)
    }

    /// The identical source/SGS preparation through the original pooled
    /// collision and streaming kernels. Sampling stays in canonical order;
    /// worker count cannot change source queries or collision arithmetic.
    /// Pooled scheduling retains its existing chunk-descriptor allocations.
    /// Source snapshots publish only after the native flow's final acceptance.
    #[allow(clippy::too_many_arguments)]
    pub fn step_sparse_pooled_with_material_rheology_and_sgs<M: FluidTransport, P: KernelRunner>(
        &self,
        grid: &mut SparseGrid3,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseRheologySgsWorkspace3<M>,
        resolve: impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<SgsReport, SparseMaterialSgsError3> {
        let report = self.prepare_sparse_rheology_sgs(grid, gate, config, workspace, resolve)?;
        grid.step_pooled_inner(pool, gate, false, Some(workspace.flow.collisions()))
            .map_err(SparseSgsError3::from)?;
        workspace.accept(report);
        Ok(report)
    }

    fn sparse_rheology_scaling(&self) -> [u64; 3] {
        [
            self.cell_size_m().to_bits(),
            self.timestep_s().to_bits(),
            self.transport().density_kg_m3().to_bits(),
        ]
    }

    fn prepare_sparse_rheology_sgs<M: FluidTransport>(
        &self,
        grid: &SparseGrid3,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseRheologySgsWorkspace3<M>,
        mut resolve: impl FnMut(MaterialImmersedPoint3) -> Result<M, MaterialError>,
    ) -> Result<SgsReport, SparseMaterialSgsError3> {
        checkpoint(gate)?;
        if workspace.scaling != self.sparse_rheology_scaling() || workspace.model != grid.model {
            return Err(SparseMaterialSgsError3::Invalid(
                "sparse rheology physical scaling or native collision law changed",
            ));
        }
        workspace.pending.clear();
        let pending = &mut workspace.pending;
        let mut origin = [0; 3];
        grid.prepare_sgs_sampled(gate, config, &mut workspace.flow, |slot, lane, gradient| {
            checkpoint(gate)?;
            let key = grid.keys[slot];
            // Canonical lane order decodes each tile only once, without a
            // second coordinate buffer or a sparse-map lookup per source.
            if lane == 0 {
                let (tx, ty, tz) = demorton3(key);
                origin = [tx, ty, tz].map(|p| p as usize * TILE);
            }
            let cell = [
                origin[0] + lane % TILE,
                origin[1] + (lane / TILE) % TILE,
                origin[2] + lane / (TILE * TILE),
            ];
            let error = |source| SparseMaterialSgsError3::Source { cell, source };
            let lattice_rate = strain_rate(gradient);
            let rate = lattice_rate / self.timestep_s();
            if !rate.is_finite() || (lattice_rate != 0.0 && rate == 0.0) {
                return Err(error(MaterialLbmError::Invalid(
                    "sparse physical shear rate is unrepresentable",
                )));
            }
            let point = MaterialImmersedPoint3 {
                cell,
                shear_rate_per_s: rate,
            };
            let material = resolve(point);
            checkpoint(gate)?;
            let material = material.map_err(|source| error(source.into()))?;
            let transport = material.transport_state();
            checkpoint(gate)?;
            let transport = transport.map_err(|source| error(source.into()))?;
            if transport.density_kg_m3().to_bits() != self.transport().density_kg_m3().to_bits() {
                return Err(error(MaterialLbmError::Invalid(
                    "sparse rheology requires unchanged reference density",
                )));
            }
            let tau = plan_momentum_scaling(
                transport.kinematic_viscosity_m2_s(),
                self.cell_size_m(),
                self.timestep_s(),
                0.0,
            )
            .map_err(error)?
            .0
            .tau;
            let collision = D3q19::with_shear_relaxation(grid.model, tau).map_err(|source| {
                SparseSgsError3::Flow(SparseError3::Collision {
                    tile_key: key,
                    lane,
                    source,
                })
            })?;
            pending.push(MaterialImmersedSample3 {
                point,
                material,
                transport,
                shear_relaxation_time: tau,
            });
            Ok(collision)
        })
    }
}

#[cfg(test)]
mod tests;
