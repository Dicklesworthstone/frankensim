//! Explicit SGS through the existing moving/open boundary and wall-load owner.
//!
//! Closure preparation uses the actual Guo-frame velocity, including a local
//! force field. It does not change any boundary, force, or collision default.
//! The regularized open-face scheme retains its existing accuracy envelope;
//! this integration does not qualify a turbulent inlet, jet or outlet model.

use super::{
    BoundaryGrid3, BoundaryLink3, BoundaryStepError3, CollisionModel3, E3, Face3, FaceBoundary3,
    Q3, StreamView, TILE_CELLS, WallMomentumExchange3, cell_force, collision_kernel3,
    effective_wall_velocity,
};
use crate::lattice3::{D3q19, ViscousLattice3};
use crate::sgs3::{SgsConfig, SgsError, SgsReport, SgsWorkspace3, eddy_viscosity};
use fs_exec::{CancelGate, KernelRunner};

/// Refusals leave the published populations unchanged. Scratch may contain
/// an incomplete preparation or collision pass and is rebuilt on retry.
#[derive(Clone, Debug, PartialEq)]
pub enum BoundarySgsError3 {
    /// Existing collision, boundary, load-measurement or cancellation error.
    Boundary(BoundaryStepError3),
    /// Invalid closure controls, acceptance use or viscosity-ratio refusal.
    Closure(SgsError),
    /// Workspace was allocated for a different cell count.
    Workspace,
    /// A reusable workspace buffer could not be allocated.
    Allocation,
    /// The resolved cell has nonfinite populations/velocity or invalid density.
    InvalidState { cell: [usize; 3] },
}

impl core::fmt::Display for BoundarySgsError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-BOUNDARY-SGS3: {self:?}")
    }
}
impl std::error::Error for BoundarySgsError3 {}
impl From<BoundaryStepError3> for BoundarySgsError3 {
    fn from(error: BoundaryStepError3) -> Self {
        Self::Boundary(error)
    }
}
impl From<SgsError> for BoundarySgsError3 {
    fn from(error: SgsError) -> Self {
        Self::Closure(error)
    }
}

/// Pre-step SGS diagnostics and optional loads from the actual same step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundarySgsReport3 {
    /// Instantaneous shear dissipation and viscosity extrema in lattice units.
    pub sgs: SgsReport,
    /// Wall impulses and work, when a link selection was requested.
    pub wall_momentum: Option<WallMomentumExchange3>,
}

#[derive(Clone, Copy)]
enum Neighbour {
    Fluid([f64; 3]),
    Wall([f64; 3]),
    Velocity,
    Pressure,
}

fn checkpoint(gate: &CancelGate) -> Result<(), BoundaryStepError3> {
    if gate.is_requested() {
        Err(BoundaryStepError3::Cancelled)
    } else {
        Ok(())
    }
}

impl BoundaryGrid3 {
    /// Allocate the velocity/density/collision cache once. Prepared collisions
    /// use this owner's native `tile * TILE_CELLS + lane` order.
    pub fn sgs_workspace(&self) -> Result<SgsWorkspace3<CollisionModel3>, BoundarySgsError3> {
        let count = self.f[0].len() * TILE_CELLS;
        let mut velocity = Vec::new();
        let mut density = Vec::new();
        let mut collisions = Vec::new();
        velocity
            .try_reserve_exact(count)
            .map_err(|_| BoundarySgsError3::Allocation)?;
        density
            .try_reserve_exact(count)
            .map_err(|_| BoundarySgsError3::Allocation)?;
        collisions
            .try_reserve_exact(count)
            .map_err(|_| BoundarySgsError3::Allocation)?;
        velocity.resize(count, [0.0; 3]);
        density.resize(count, 0.0);
        collisions.resize(count, self.collision_model);
        Ok(SgsWorkspace3 {
            velocity,
            density,
            collisions,
        })
    }

    /// Prepare Smagorinsky/WALE from the current full velocity gradient.
    /// Molecular viscosity always comes from the constructor collision law;
    /// a previous SGS field is never reused as molecular transport.
    ///
    /// Periodic/interior derivatives are centred. Halfway walls use quadratic
    /// reconstruction through their actual moving or stationary value. The
    /// wall/open and wall/wall corner convention is the streaming owner's.
    /// On-site velocity boundaries use a quadratic one-sided derivative when
    /// two interior fluid cells exist, otherwise a linear derivative. Pressure
    /// outlets use their declared zero-normal-velocity-gradient condition.
    /// These stencils add no allocations and poll at least every 64 cells.
    pub fn prepare_sgs(
        &self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SgsWorkspace3<CollisionModel3>,
    ) -> Result<SgsReport, BoundarySgsError3> {
        checkpoint(gate)?;
        config.validate()?;
        let count = self.f[0].len() * TILE_CELLS;
        if workspace.velocity.len() != count
            || workspace.density.len() != count
            || workspace.collisions.len() != count
        {
            return Err(BoundarySgsError3::Workspace);
        }
        let nu = self.collision_model.kinematic_viscosity();
        if !(nu.is_finite() && nu > 0.0) {
            return Err(SgsError::InvalidMolecularViscosity { cell: 0 }.into());
        }
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let cell = tile * TILE_CELLS + lane;
                if self.solid[tile] & (1u64 << lane) != 0 {
                    workspace.velocity[cell] = [0.0; 3];
                    workspace.density[cell] = 0.0;
                    continue;
                }
                let mut rho = 0.0;
                let mut momentum = [0.0; 3];
                for (q, e) in E3.iter().enumerate() {
                    let population = self.f[q][tile].0[lane];
                    if !population.is_finite() {
                        return Err(BoundarySgsError3::InvalidState {
                            cell: self.coords(tile, lane),
                        });
                    }
                    rho += population;
                    momentum[0] += f64::from(e.0) * population;
                    momentum[1] += f64::from(e.1) * population;
                    momentum[2] += f64::from(e.2) * population;
                }
                let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                let velocity = core::array::from_fn(|a| (momentum[a] + 0.5 * force[a]) / rho);
                if !(rho.is_finite() && rho > 0.0) || !velocity.iter().all(|v| v.is_finite()) {
                    return Err(BoundarySgsError3::InvalidState {
                        cell: self.coords(tile, lane),
                    });
                }
                workspace.velocity[cell] = velocity;
                workspace.density[cell] = rho;
            }
        }
        let mut report = SgsReport::default();
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let cell = tile * TILE_CELLS + lane;
                workspace.collisions[cell] = self.collision_model;
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let gradient = self.sgs_gradient(self.coords(tile, lane), &workspace.velocity);
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
                if eddy != 0.0 {
                    workspace.collisions[cell] =
                        D3q19::with_shear_relaxation(self.collision_model, 0.5 + 3.0 * (nu + eddy))
                            .map_err(|source| BoundaryStepError3::Collision {
                                cell: self.coords(tile, lane),
                                source,
                            })?;
                }
                let dilation = (gradient[0][0] + gradient[1][1] + gradient[2][2]) / 3.0;
                let mut strain2 = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        let s = 0.5 * (gradient[i][j] + gradient[j][i])
                            - if i == j { dilation } else { 0.0 };
                        strain2 += s * s;
                    }
                }
                let factor = 2.0 * workspace.density[cell] * strain2;
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

    /// Serial SGS step using the original collision, streaming, boundary and
    /// wall-load kernels. Pass a canonical wall-link selection to measure
    /// loads on the same SGS-augmented post-collision state. `None` skips it.
    ///
    /// Failure publishes neither populations nor loads; the existing stream
    /// scratch is reused. The zero coefficient preserves ordinary step bits.
    pub fn step_with_sgs(
        &mut self,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SgsWorkspace3<CollisionModel3>,
        measurement: Option<(&[BoundaryLink3], [f64; 3])>,
    ) -> Result<BoundarySgsReport3, BoundarySgsError3> {
        if let Some((links, origin)) = measurement {
            self.validate_wall_measurement(links, origin)?;
        }
        let sgs = self.prepare_sgs(gate, config, workspace)?;
        self.topology_locked = true;
        let collide = collision_kernel3();
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let f = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                let post = collide(f, workspace.collisions[tile * TILE_CELLS + lane], force)
                    .map_err(|source| BoundaryStepError3::Collision {
                        cell: self.coords(tile, lane),
                        source,
                    })?;
                for (field, value) in self.post.iter_mut().zip(post) {
                    field[tile].0[lane] = value;
                }
            }
        }
        let wall_momentum = match measurement {
            Some((links, origin)) => Some(self.wall_momentum(links, origin, Some(gate))?),
            None => None,
        };
        let mut next = self.next.take().unwrap_or_else(|| self.f.clone());
        let outcome: Result<(), BoundaryStepError3> = (|| {
            let view = StreamView {
                post: &self.post,
                pull: &self.pull,
                link_masks: &self.link_masks,
                open_link_masks: &self.open_link_masks,
                stationary_link_masks: &self.stationary_link_masks,
                boundaries: self.boundaries,
                dims: [self.nx, self.ny, self.nz],
            };
            for tile in 0..self.f[0].len() {
                checkpoint(gate)?;
                for lane in 0..TILE_CELLS {
                    if self.solid[tile] & (1u64 << lane) != 0 {
                        continue;
                    }
                    let values = view.pull_cell(tile, lane);
                    for (field, value) in next.iter_mut().zip(values) {
                        field[tile].0[lane] = value;
                    }
                }
            }
            Ok(())
        })();
        if let Err(error) = outcome {
            self.next = Some(next);
            return Err(error.into());
        }
        self.publish_next(next, gate, true)?;
        Ok(BoundarySgsReport3 { sgs, wall_momentum })
    }

    /// Pooled SGS collision/streaming and optional wall loads through the
    /// existing transactional step owner. Preparation and load reduction use
    /// the same canonical order as the serial path for every worker count.
    pub fn step_pooled_with_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: SgsConfig,
        workspace: &mut SgsWorkspace3<CollisionModel3>,
        measurement: Option<(&[BoundaryLink3], [f64; 3])>,
    ) -> Result<BoundarySgsReport3, BoundarySgsError3> {
        if let Some((links, origin)) = measurement {
            self.validate_wall_measurement(links, origin)?;
        }
        let sgs = self.prepare_sgs(gate, config, workspace)?;
        let exchange =
            self.step_pooled_inner(pool, gate, measurement, Some(&workspace.collisions))?;
        Ok(BoundarySgsReport3 {
            sgs,
            wall_momentum: measurement.map(|_| exchange),
        })
    }

    fn sgs_neighbour(
        &self,
        p: [usize; 3],
        axis: usize,
        positive: bool,
        velocity: &[[f64; 3]],
    ) -> Neighbour {
        let dims = [self.nx, self.ny, self.nz];
        let outside = if positive {
            p[axis] + 1 == dims[axis]
        } else {
            p[axis] == 0
        };
        let mut other = p;
        if outside {
            let face = Face3::ALL[2 * axis + usize::from(positive)];
            match self.boundaries.face(face) {
                FaceBoundary3::Periodic => other[axis] = if positive { 0 } else { dims[axis] - 1 },
                FaceBoundary3::Wall { .. } => {
                    let (tile, lane) = self.addr(p[0], p[1], p[2]);
                    let q = 2 * axis + if positive { 2 } else { 1 };
                    let wall = if self.stationary_link_masks[tile].0[lane] & (1u32 << q) != 0 {
                        [0.0; 3]
                    } else {
                        effective_wall_velocity(self.boundaries, dims, p)
                    };
                    return Neighbour::Wall(wall);
                }
                FaceBoundary3::Velocity { .. } => return Neighbour::Velocity,
                FaceBoundary3::Pressure { .. } => return Neighbour::Pressure,
            }
        } else {
            other[axis] = if positive { p[axis] + 1 } else { p[axis] - 1 };
        }
        if self.is_solid(other[0], other[1], other[2]) {
            Neighbour::Wall([0.0; 3])
        } else {
            let (tile, lane) = self.addr(other[0], other[1], other[2]);
            Neighbour::Fluid(velocity[tile * TILE_CELLS + lane])
        }
    }

    fn sgs_gradient(&self, p: [usize; 3], velocity: &[[f64; 3]]) -> [[f64; 3]; 3] {
        let (tile, lane) = self.addr(p[0], p[1], p[2]);
        let center = velocity[tile * TILE_CELLS + lane];
        let mut gradient = [[0.0; 3]; 3];
        for axis in 0..3 {
            let low = self.sgs_neighbour(p, axis, false, velocity);
            let high = self.sgs_neighbour(p, axis, true, velocity);
            for component in 0..3 {
                gradient[component][axis] = match (low, high) {
                    (Neighbour::Pressure, _) | (_, Neighbour::Pressure) => 0.0,
                    (Neighbour::Fluid(low), Neighbour::Fluid(high)) => {
                        0.5 * (high[component] - low[component])
                    }
                    (Neighbour::Wall(wall), Neighbour::Fluid(high)) => {
                        (center[component] - wall[component])
                            + (high[component] - wall[component]) / 3.0
                    }
                    (Neighbour::Fluid(low), Neighbour::Wall(wall)) => {
                        -(center[component] - wall[component])
                            - (low[component] - wall[component]) / 3.0
                    }
                    (Neighbour::Wall(low), Neighbour::Wall(high)) => {
                        high[component] - low[component]
                    }
                    (Neighbour::Velocity, Neighbour::Fluid(first)) => {
                        self.sgs_open_derivative(p, axis, false, component, center, first, velocity)
                    }
                    (Neighbour::Fluid(first), Neighbour::Velocity) => {
                        self.sgs_open_derivative(p, axis, true, component, center, first, velocity)
                    }
                    (Neighbour::Velocity, Neighbour::Wall(wall)) => {
                        2.0 * (wall[component] - center[component])
                    }
                    (Neighbour::Wall(wall), Neighbour::Velocity) => {
                        2.0 * (center[component] - wall[component])
                    }
                    (Neighbour::Velocity, Neighbour::Velocity) => 0.0,
                };
            }
        }
        gradient
    }

    fn sgs_open_derivative(
        &self,
        p: [usize; 3],
        axis: usize,
        upper: bool,
        component: usize,
        center: [f64; 3],
        first: [f64; 3],
        velocity: &[[f64; 3]],
    ) -> f64 {
        let mut second = p;
        second[axis] = if upper { p[axis] - 2 } else { p[axis] + 2 };
        let sign = if upper { -1.0 } else { 1.0 };
        let difference = first[component] - center[component];
        if self.is_solid(second[0], second[1], second[2]) {
            return sign * difference;
        }
        let (tile, lane) = self.addr(second[0], second[1], second[2]);
        sign * (2.0 * difference
            - 0.5 * (velocity[tile * TILE_CELLS + lane][component] - center[component]))
    }
}
