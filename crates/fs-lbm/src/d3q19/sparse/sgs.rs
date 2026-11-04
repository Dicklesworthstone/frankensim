//! Explicit SGS viscosity on the original sparse collision/streaming owner.
//!
//! Only active cells are retained. Six neighbor slots per tile make gradient
//! reconstruction independent of sparse-map lookups in the per-cell sweep.
//! Inactive tiles and external links keep their stationary halfway-wall
//! meaning; they are not empty fluid, periodic wraparound, or free surfaces.

use super::{CollisionModel3, SparseError3, SparseGrid3, TILE, TILE_CELLS, demorton3, morton3};
use crate::lattice3::{D3q19, Lattice3, ViscousLattice3};
use crate::sgs3::{SgsConfig, SgsError, SgsReport, SgsWorkspace3, eddy_viscosity};
use fs_exec::{CancelGate, KernelRunner};

/// Retained source responses supplying spatial molecular viscosity under SGS.
pub mod material;

/// Refused preparation or native sparse step. Published populations and step
/// count stay intact; a partially prepared workspace is not a valid output.
#[derive(Clone, Debug, PartialEq)]
pub enum SparseSgsError3 {
    Flow(SparseError3),
    Closure(SgsError),
    Allocation,
    EmptyGrid,
    /// Active Morton keys or domain dimensions differ from the prepared plan.
    TopologyChanged,
    MolecularFieldLength {
        expected: usize,
        actual: usize,
    },
    InvalidState {
        tile_key: u64,
        lane: usize,
    },
}

impl core::fmt::Display for SparseSgsError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-SPARSE-SGS3: {self:?}")
    }
}
impl std::error::Error for SparseSgsError3 {}
impl From<SparseError3> for SparseSgsError3 {
    fn from(error: SparseError3) -> Self {
        Self::Flow(error)
    }
}
impl From<SgsError> for SparseSgsError3 {
    fn from(error: SgsError) -> Self {
        Self::Closure(error)
    }
}

/// Reusable fields and axial adjacency for one active set. Memory depends on
/// the active tiles, never the enclosing domain volume. Construct a new
/// workspace after activation/deactivation changes the set; equal cell counts
/// alone do not identify a topology.
pub struct SparseSgsWorkspace3 {
    tile_dims: (usize, usize, usize),
    keys: Vec<u64>,
    neighbors: Vec<[Option<usize>; 6]>,
    local: SgsWorkspace3<CollisionModel3>,
}

impl SparseSgsWorkspace3 {
    /// Morton keys corresponding to contiguous 64-cell collision blocks.
    pub fn tile_keys(&self) -> &[u64] {
        &self.keys
    }

    /// Valid after successful preparation, in `slot * 64 + lane` order.
    /// Lanes are x-fastest within each 4x4x4 tile.
    pub fn collisions(&self) -> &[CollisionModel3] {
        self.local.collisions()
    }

    /// Retained numerical capacities, excluding allocator and struct overhead.
    pub fn workspace_bytes(&self) -> usize {
        self.keys.capacity() * core::mem::size_of::<u64>()
            + self.neighbors.capacity() * core::mem::size_of::<[Option<usize>; 6]>()
            + self.local.velocity.capacity() * core::mem::size_of::<[f64; 3]>()
            + self.local.density.capacity() * core::mem::size_of::<f64>()
            + self.local.collisions.capacity() * core::mem::size_of::<CollisionModel3>()
    }

    fn neighbor(&self, slot: usize, lane: usize, axis: usize, positive: bool) -> Option<[f64; 3]> {
        let stride = [1, TILE, TILE * TILE][axis];
        let at = (lane / stride) % TILE;
        let (other_slot, other_lane) = if positive {
            if at + 1 < TILE {
                (slot, lane + stride)
            } else {
                (
                    self.neighbors[slot][2 * axis + 1]?,
                    lane - (TILE - 1) * stride,
                )
            }
        } else if at > 0 {
            (slot, lane - stride)
        } else {
            (self.neighbors[slot][2 * axis]?, lane + (TILE - 1) * stride)
        };
        Some(self.local.velocity[other_slot * TILE_CELLS + other_lane])
    }

    fn gradient(&self, slot: usize, lane: usize) -> [[f64; 3]; 3] {
        let mut gradient = [[0.0; 3]; 3];
        let center = self.local.velocity[slot * TILE_CELLS + lane];
        for axis in 0..3 {
            let low = self.neighbor(slot, lane, axis, false);
            let high = self.neighbor(slot, lane, axis, true);
            for component in 0..3 {
                // Same quadratic halfway-wall stencil as Grid3. The wall,
                // cell center and next fluid center lie at -1/2, 0 and 1.
                gradient[component][axis] = match (low, high) {
                    (Some(low), Some(high)) => 0.5 * (high[component] - low[component]),
                    (None, Some(high)) => center[component] + high[component] / 3.0,
                    (Some(low), None) => -center[component] - low[component] / 3.0,
                    (None, None) => 0.0,
                };
            }
        }
        gradient
    }
}

fn checkpoint(gate: &CancelGate) -> Result<(), SparseSgsError3> {
    if gate.is_requested() {
        Err(SparseError3::Cancelled.into())
    } else {
        Ok(())
    }
}

fn allocated<T: Clone>(count: usize, value: T) -> Result<Vec<T>, SparseSgsError3> {
    let mut out = Vec::new();
    out.try_reserve_exact(count)
        .map_err(|_| SparseSgsError3::Allocation)?;
    out.resize(count, value);
    Ok(out)
}

impl SparseGrid3 {
    /// Prepare reusable active-cell fields and tile-neighbor indices. No SGS
    /// model is enabled by constructing a workspace. Empty active sets refuse.
    pub fn sgs_workspace(&self, gate: &CancelGate) -> Result<SparseSgsWorkspace3, SparseSgsError3> {
        checkpoint(gate)?;
        if self.keys.is_empty() {
            return Err(SparseSgsError3::EmptyGrid);
        }
        let cells = self
            .keys
            .len()
            .checked_mul(TILE_CELLS)
            .ok_or(SparseSgsError3::Allocation)?;
        let mut keys = allocated(self.keys.len(), 0)?;
        let mut neighbors = allocated(self.keys.len(), [None; 6])?;
        let dims = [self.ntx, self.nty, self.ntz];
        for (slot, &key) in self.keys.iter().enumerate() {
            checkpoint(gate)?;
            keys[slot] = key;
            let (tx, ty, tz) = demorton3(key);
            neighbors[slot] = core::array::from_fn(|side| {
                let axis = side / 2;
                let mut p = [tx, ty, tz];
                let n = i64::from(p[axis]) + if side % 2 == 0 { -1 } else { 1 };
                if n < 0 || n as usize >= dims[axis] {
                    return None;
                }
                p[axis] = n as u32;
                self.index.get(&morton3(p[0], p[1], p[2])).copied()
            });
        }
        let local = SgsWorkspace3 {
            velocity: allocated(cells, [0.0; 3])?,
            density: allocated(cells, 0.0)?,
            collisions: allocated(cells, self.model)?,
        };
        checkpoint(gate)?;
        Ok(SparseSgsWorkspace3 {
            tile_dims: self.tile_dims(),
            keys,
            neighbors,
            local,
        })
    }

    /// Build total-viscosity collisions from the accepted, force-centered
    /// velocity. Uses the shared full-gradient Smagorinsky/WALE/Vreman laws,
    /// explicit filter width and viscosity-ratio refusal; no clipping occurs.
    /// Molecular viscosity always comes from this grid's native law, never a
    /// prior augmented field. Preparation allocates no memory and performs no
    /// sparse-map lookup per cell. Diagnostics follow canonical Morton order.
    pub fn prepare_sgs(
        &self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseSgsWorkspace3,
    ) -> Result<SgsReport, SparseSgsError3> {
        self.prepare_sgs_molecular(gate, config, workspace, None)
    }

    fn prepare_sgs_molecular(
        &self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseSgsWorkspace3,
        molecular: Option<&[CollisionModel3]>,
    ) -> Result<SgsReport, SparseSgsError3> {
        checkpoint(gate)?;
        config.validate()?;
        if workspace.tile_dims != self.tile_dims() || workspace.keys != self.keys {
            return Err(SparseSgsError3::TopologyChanged);
        }
        if let Some(field) = molecular {
            if field.len() != workspace.local.collisions.len() {
                return Err(SparseSgsError3::MolecularFieldLength {
                    expected: workspace.local.collisions.len(),
                    actual: field.len(),
                });
            }
        }
        self.model
            .validate()
            .map_err(|source| SparseError3::Collision {
                tile_key: self.keys[0],
                lane: 0,
                source,
            })?;
        let uniform_nu = self.model.kinematic_viscosity();
        if !uniform_nu.is_finite() || uniform_nu <= 0.0 {
            return Err(SgsError::InvalidMolecularViscosity { cell: 0 }.into());
        }
        for (slot, &key) in self.keys.iter().enumerate() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let f = self.populations(slot, lane);
                let rho: f64 = f.iter().sum();
                let u: [f64; 3] = core::array::from_fn(|a| {
                    (f.iter()
                        .enumerate()
                        .map(|(q, p)| f64::from(D3q19::E[q][a]) * p)
                        .sum::<f64>()
                        + 0.5 * self.force[a])
                        / rho
                });
                if !rho.is_finite() || rho <= 0.0 || f.iter().chain(&u).any(|v| !v.is_finite()) {
                    return Err(SparseSgsError3::InvalidState {
                        tile_key: key,
                        lane,
                    });
                }
                let cell = slot * TILE_CELLS + lane;
                workspace.local.velocity[cell] = u;
                workspace.local.density[cell] = rho;
            }
        }
        let mut report = SgsReport::default();
        for (slot, &key) in self.keys.iter().enumerate() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let cell = slot * TILE_CELLS + lane;
                let original = molecular.map_or(self.model, |field| field[cell]);
                let nu = if molecular.is_some() {
                    original
                        .validate()
                        .map_err(|source| SparseError3::Collision {
                            tile_key: key,
                            lane,
                            source,
                        })?;
                    original.kinematic_viscosity()
                } else {
                    uniform_nu
                };
                if !nu.is_finite() || nu <= 0.0 {
                    return Err(SgsError::InvalidMolecularViscosity { cell }.into());
                }
                let gradient = workspace.gradient(slot, lane);
                let eddy = eddy_viscosity(config.model, config.filter_width_lu, gradient)?;
                let ratio = eddy / nu;
                if !ratio.is_finite() || ratio > config.maximum_eddy_to_molecular {
                    return Err(SgsError::ViscosityRatioExceeded {
                        cell,
                        ratio,
                        maximum: config.maximum_eddy_to_molecular,
                    }
                    .into());
                }
                workspace.local.collisions[cell] = if eddy == 0.0 {
                    original
                } else {
                    D3q19::with_shear_relaxation(original, 0.5 + (nu + eddy) / D3q19::CS2).map_err(
                        |source| SparseError3::Collision {
                            tile_key: key,
                            lane,
                            source,
                        },
                    )?
                };
                let mut strain2 = 0.0;
                let dilation_third = (gradient[0][0] + gradient[1][1] + gradient[2][2]) / 3.0;
                for i in 0..3 {
                    for j in 0..3 {
                        let s = 0.5 * (gradient[i][j] + gradient[j][i])
                            - if i == j { dilation_third } else { 0.0 };
                        strain2 += s * s;
                    }
                }
                let factor = 2.0 * workspace.local.density[cell] * strain2;
                report.fluid_cells += 1;
                report.resolved_dissipation_lu += factor * nu;
                report.sgs_dissipation_lu += factor * eddy;
                report.maximum_eddy_viscosity_lu = report.maximum_eddy_viscosity_lu.max(eddy);
                report.maximum_eddy_to_molecular = report.maximum_eddy_to_molecular.max(ratio);
            }
        }
        checkpoint(gate)?;
        if !report.resolved_dissipation_lu.is_finite() || !report.sgs_dissipation_lu.is_finite() {
            return Err(SgsError::NonFiniteViscosity.into());
        }
        Ok(report)
    }

    /// Advance the original sparse collide/stream owner with prepared SGS
    /// viscosity. Pre-step dissipation is diagnostic, not added heat. Neither
    /// model nor viscosity persists into a subsequent ordinary step. Zero
    /// coefficient retains the native scalar-collision path and its bits.
    pub fn step_serial_with_sgs(
        &mut self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseSgsWorkspace3,
    ) -> Result<SgsReport, SparseSgsError3> {
        self.step_serial_with_sgs_molecular(gate, config, workspace, None)
    }

    fn step_serial_with_sgs_molecular(
        &mut self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseSgsWorkspace3,
        molecular: Option<&[CollisionModel3]>,
    ) -> Result<SgsReport, SparseSgsError3> {
        let report = self.prepare_sgs_molecular(gate, config, workspace, molecular)?;
        // Disabling SGS must retain a spatial molecular field. Only the
        // original uniform owner can return to its scalar-collision path.
        let field = (molecular.is_some() || config.model.coefficient() != 0.0)
            .then_some(workspace.collisions());
        self.step_serial_inner(Some(gate), field)?;
        Ok(report)
    }

    /// The same SGS preparation and original pooled collide/stream kernels.
    /// Worker count changes neither closure arithmetic nor populations. SGS
    /// introduces no step allocation; the existing pool's chunk descriptors
    /// remain its own scheduling cost. Cancellation publishes no candidate.
    pub fn step_pooled_with_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseSgsWorkspace3,
    ) -> Result<SgsReport, SparseSgsError3> {
        self.step_pooled_with_sgs_molecular(pool, gate, config, workspace, None)
    }

    fn step_pooled_with_sgs_molecular<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SparseSgsWorkspace3,
        molecular: Option<&[CollisionModel3]>,
    ) -> Result<SgsReport, SparseSgsError3> {
        let report = self.prepare_sgs_molecular(gate, config, workspace, molecular)?;
        let field = (molecular.is_some() || config.model.coefficient() != 0.0)
            .then_some(workspace.collisions());
        self.step_pooled_inner(pool, gate, false, field)?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
