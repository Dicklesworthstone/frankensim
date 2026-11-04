//! Original material responses under the native sparse SGS trajectory.
//!
//! One frozen source epoch supplies molecular viscosity. Resolved temperature,
//! composition and source-law evaluation stay with their existing owners; this
//! adapter does not invent a second material model or an energy equation.
//! [`rheology`] resamples those original sources at incoming resolved shear.

use super::{SparseGrid3, SparseSgsError3, SparseSgsWorkspace3, checkpoint};
use crate::d3q19::sparse::{SparseError3, demorton3};
use crate::d3q19::{CollisionModel3, TILE};
use crate::lattice3::{D3q19, ViscousLattice3};
use crate::material::{MaterialLbmError, MaterialScaling, plan_momentum_scaling};
use crate::sgs3::{SgsConfig, SgsReport};
use fs_exec::{CancelGate, KernelRunner};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};

/// Live source viscosity evaluated at the current resolved shear rate.
pub mod rheology;

/// Complete original response and the molecular relaxation applied at one
/// active cell. Borrowed concrete sources and `&dyn FluidTransport` work too.
pub struct SparseMaterialSample3<M> {
    pub cell: [usize; 3],
    pub material: M,
    pub transport: FluidTransportState,
    pub molecular_relaxation_time: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SparseMaterialSgsError3 {
    Sgs(SparseSgsError3),
    Invalid(&'static str),
    Source {
        cell: [usize; 3],
        source: MaterialLbmError,
    },
}

impl core::fmt::Display for SparseMaterialSgsError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-SPARSE-MATERIAL-SGS3: {self:?}")
    }
}
impl std::error::Error for SparseMaterialSgsError3 {}
impl From<SparseSgsError3> for SparseMaterialSgsError3 {
    fn from(error: SparseSgsError3) -> Self {
        Self::Sgs(error)
    }
}

/// Immutable active-cell material epoch. Original source responses and their
/// admitted transport snapshots remain together; stepping never resamples
/// them or reads a previously SGS-augmented viscosity as molecular input.
/// Rebind after the caller accepts new thermal/composition/source state.
pub struct SparseMaterialSgs3<M> {
    tile_dims: (usize, usize, usize),
    keys: Vec<u64>,
    model: CollisionModel3,
    molecular: Vec<CollisionModel3>,
    samples: Vec<SparseMaterialSample3<M>>,
}

impl<M> SparseMaterialSgs3<M> {
    /// Canonical Morton-tile/lane order, containing no inactive-cell padding.
    pub fn samples(&self) -> &[SparseMaterialSample3<M>] {
        &self.samples
    }

    /// Original molecular collisions in the same order as `samples()`.
    pub fn molecular_collisions(&self) -> &[CollisionModel3] {
        &self.molecular
    }

    fn admit(&self, grid: &SparseGrid3, gate: &CancelGate) -> Result<(), SparseMaterialSgsError3> {
        checkpoint(gate)?;
        if self.tile_dims != grid.tile_dims() || self.keys != grid.keys {
            return Err(SparseSgsError3::TopologyChanged.into());
        }
        if self.model != grid.model {
            return Err(SparseMaterialSgsError3::Invalid(
                "native collision law changed after material binding",
            ));
        }
        Ok(())
    }

    /// Prepare total viscosity and molecular/SGS shear diagnostics from this
    /// source epoch and the grid's CURRENT accepted velocity. The SGS/molecular
    /// ratio is checked against each cell's own sourced molecular viscosity.
    pub fn prepare_sgs(
        &self,
        grid: &SparseGrid3,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseSgsWorkspace3,
    ) -> Result<SgsReport, SparseMaterialSgsError3> {
        self.admit(grid, gate)?;
        Ok(grid.prepare_sgs_molecular(gate, config, workspace, Some(&self.molecular))?)
    }

    /// Advance the native sparse owner with this local molecular field plus
    /// the explicit SGS closure. Zero SGS retains the heterogeneous molecular
    /// field. No source sampling, allocation or fallible work follows the
    /// native publication boundary; refusals preserve populations and steps.
    pub fn step_serial(
        &self,
        grid: &mut SparseGrid3,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseSgsWorkspace3,
    ) -> Result<SgsReport, SparseMaterialSgsError3> {
        self.admit(grid, gate)?;
        Ok(grid.step_serial_with_sgs_molecular(gate, config, workspace, Some(&self.molecular))?)
    }

    /// The same retained source field and SGS arithmetic through the original
    /// pooled collision/streaming owner and its failure-atomic publication.
    pub fn step_pooled<P: KernelRunner>(
        &self,
        grid: &mut SparseGrid3,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseSgsWorkspace3,
    ) -> Result<SgsReport, SparseMaterialSgsError3> {
        self.admit(grid, gate)?;
        Ok(grid.step_pooled_with_sgs_molecular(
            pool,
            gate,
            config,
            workspace,
            Some(&self.molecular),
        )?)
    }
}

impl<R: FluidTransport + ?Sized> MaterialScaling<'_, R> {
    /// Resolve one original material response per active sparse cell and bind
    /// its molecular viscosity to this scaling's rho0, dx and dt. Only active
    /// cells are sampled/allocated, in Morton-tile then x-fastest lane order.
    /// Borrow an already resolved source to retain a held material epoch.
    ///
    /// Every response must retain the SAME reference density: this isothermal
    /// incompressible owner cannot remap mass, inertia or the pressure datum.
    /// The existing momentum-scaling routine converts mu/rho0 to relaxation;
    /// heat capacity and conductivity are retained but do not impose a thermal
    /// lattice bound or create heat evolution. The grid's original collision
    /// family and non-shear controls remain authoritative.
    ///
    /// Keep this reference scaling for physical interpretation and maintain
    /// its scenario speed bound. Temperature, species transport, nonlinear
    /// shear-law closure and source updates remain caller-owned. Callback
    /// side effects cannot be rolled back: do not publish external state from
    /// `resolve`. Admission/source/cancellation failures mutate no grid state
    /// and return no partial material binding.
    pub fn bind_sparse_sgs_materials<M: FluidTransport>(
        &self,
        grid: &SparseGrid3,
        max_active_cells: usize,
        gate: &CancelGate,
        mut resolve: impl FnMut([usize; 3]) -> Result<M, MaterialError>,
    ) -> Result<SparseMaterialSgs3<M>, SparseMaterialSgsError3> {
        checkpoint(gate)?;
        if grid.keys.is_empty() {
            return Err(SparseSgsError3::EmptyGrid.into());
        }
        let count = grid
            .keys
            .len()
            .checked_mul(TILE * TILE * TILE)
            .filter(|&n| n <= max_active_cells)
            .ok_or(SparseMaterialSgsError3::Invalid(
                "active material cell budget exceeded",
            ))?;
        grid.model.validate().map_err(|source| {
            SparseSgsError3::Flow(SparseError3::Collision {
                tile_key: grid.keys[0],
                lane: 0,
                source,
            })
        })?;
        let mut keys = super::allocated(grid.keys.len(), 0)?;
        let mut molecular = Vec::new();
        let mut samples = Vec::new();
        molecular
            .try_reserve_exact(count)
            .map_err(|_| SparseSgsError3::Allocation)?;
        samples
            .try_reserve_exact(count)
            .map_err(|_| SparseSgsError3::Allocation)?;
        for (slot, &key) in grid.keys.iter().enumerate() {
            checkpoint(gate)?;
            keys[slot] = key;
            let (tx, ty, tz) = demorton3(key);
            for lane in 0..TILE * TILE * TILE {
                let cell = [
                    tx as usize * TILE + lane % TILE,
                    ty as usize * TILE + (lane / TILE) % TILE,
                    tz as usize * TILE + lane / (TILE * TILE),
                ];
                let error = |source| SparseMaterialSgsError3::Source { cell, source };
                let material = resolve(cell).map_err(|source| error(source.into()))?;
                let transport = material
                    .transport_state()
                    .map_err(|source| error(source.into()))?;
                if transport.density_kg_m3().to_bits() != self.transport().density_kg_m3().to_bits()
                {
                    return Err(error(MaterialLbmError::Invalid(
                        "sparse material viscosity requires unchanged reference density",
                    )));
                }
                let (scaled, _) = plan_momentum_scaling(
                    transport.kinematic_viscosity_m2_s(),
                    self.cell_size_m(),
                    self.timestep_s(),
                    0.0,
                )
                .map_err(error)?;
                let collision =
                    D3q19::with_shear_relaxation(grid.model, scaled.tau).map_err(|source| {
                        SparseSgsError3::Flow(SparseError3::Collision {
                            tile_key: key,
                            lane,
                            source,
                        })
                    })?;
                molecular.push(collision);
                samples.push(SparseMaterialSample3 {
                    cell,
                    material,
                    transport,
                    molecular_relaxation_time: scaled.tau,
                });
            }
        }
        checkpoint(gate)?;
        Ok(SparseMaterialSgs3 {
            tile_dims: grid.tile_dims(),
            keys,
            model: grid.model,
            molecular,
            samples,
        })
    }
}

#[cfg(test)]
mod tests;
