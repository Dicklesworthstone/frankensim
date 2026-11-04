//! Nonlocal dynamic closure through the native boundary collision/load owner.

use super::{
    BoundaryGrid3, BoundaryLink3, BoundarySgsError3, BoundaryStepError3, CollisionModel3, D3q19,
    Face3, FaceBoundary3, SgsError, SgsReport, SgsWorkspace3, TILE_CELLS, ViscousLattice3,
    WallMomentumExchange3, checkpoint,
};
use crate::sgs3::DynamicSgsStepReport;
use crate::sgs3::dynamic::channel::{ChannelDynamicSgs, ChannelDynamicSgsConfig};
use crate::sgs3::dynamic::{DynamicSgsConfig, DynamicSgsError, PeriodicDynamicSgs};
use fs_exec::{CancelGate, KernelRunner};

/// Explicit dynamic fit and viscosity limits for the actual boundary topology.
/// No closure is enabled by default, and acceptance rendering is refused.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BoundaryDynamicSgsConfig3 {
    /// Fully periodic box, with domain or local-box averaging as requested.
    Periodic(DynamicSgsConfig),
    /// One impermeable wall pair and two periodic homogeneous directions.
    /// Contractions are averaged within each wall-parallel plane separately.
    Channel(ChannelDynamicSgsConfig),
}

/// A failed preparation or step never publishes changed populations or loads.
/// Workspace outputs are usable only after a successful preparation or step.
#[derive(Clone, Debug, PartialEq)]
pub enum BoundaryDynamicSgsError3 {
    Boundary(BoundarySgsError3),
    Dynamic(DynamicSgsError),
}

impl core::fmt::Display for BoundaryDynamicSgsError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LBM-BOUNDARY-DYNAMIC-SGS3: {self:?}")
    }
}
impl std::error::Error for BoundaryDynamicSgsError3 {}
impl From<BoundarySgsError3> for BoundaryDynamicSgsError3 {
    fn from(error: BoundarySgsError3) -> Self {
        Self::Boundary(error)
    }
}
impl From<BoundaryStepError3> for BoundaryDynamicSgsError3 {
    fn from(error: BoundaryStepError3) -> Self {
        Self::Boundary(error.into())
    }
}
impl From<DynamicSgsError> for BoundaryDynamicSgsError3 {
    fn from(error: DynamicSgsError) -> Self {
        Self::Dynamic(error)
    }
}

enum DynamicFilter {
    Periodic(PeriodicDynamicSgs),
    Channel(ChannelDynamicSgs),
}

impl DynamicFilter {
    fn shape(&self) -> [usize; 3] {
        match self {
            Self::Periodic(filter) => filter.shape(),
            Self::Channel(filter) => filter.shape(),
        }
    }

    fn normal_axis(&self) -> Option<usize> {
        match self {
            Self::Periodic(_) => None,
            Self::Channel(filter) => Some(filter.normal_axis()),
        }
    }

    fn viscosities(&self) -> &[f64] {
        match self {
            Self::Periodic(filter) => filter.viscosities(),
            Self::Channel(filter) => filter.viscosities(),
        }
    }

    fn workspace_bytes(&self) -> usize {
        match self {
            Self::Periodic(filter) => filter.workspace_bytes(),
            Self::Channel(filter) => filter.workspace_bytes(),
        }
    }
}

/// Reused native macro/collision fields plus x-fastest filter inputs and fit.
/// Preparation and stepping allocate no numerical buffers. The domain shape
/// and periodic/channel topology are fixed; current wall velocities are read
/// from the grid for every preparation.
pub struct BoundaryDynamicSgsWorkspace3 {
    local: SgsWorkspace3<CollisionModel3>,
    density: Vec<f64>,
    velocity: Vec<[f64; 3]>,
    filter: DynamicFilter,
}

impl BoundaryDynamicSgsWorkspace3 {
    /// Minimum numerical bytes, excluding the grid, allocator overhead and
    /// small metadata. Admission also checks actual capacities after allocation.
    pub fn required_bytes(cells: usize) -> Result<usize, DynamicSgsError> {
        cells
            .checked_mul(25 * core::mem::size_of::<f64>() + core::mem::size_of::<CollisionModel3>())
            .ok_or(DynamicSgsError::Invalid(
                "boundary dynamic workspace size overflow",
            ))
    }

    /// Actual numerical capacities, in bytes, excluding allocator metadata.
    pub fn workspace_bytes(&self) -> usize {
        (3 * (self.local.velocity.capacity() + self.velocity.capacity())
            + self.local.density.capacity()
            + self.density.capacity())
            * core::mem::size_of::<f64>()
            + self.local.collisions.capacity() * core::mem::size_of::<CollisionModel3>()
            + self.filter.workspace_bytes()
    }

    /// Total-viscosity collisions in the boundary grid's native tile order.
    pub fn collisions(&self) -> &[CollisionModel3] {
        self.local.collisions()
    }

    /// Eddy viscosities in x-fastest order, independently of native tiling.
    pub fn eddy_viscosities(&self) -> &[f64] {
        self.filter.viscosities()
    }
}

/// Dynamic fit and pre-step shear dissipation, with actual same-step wall loads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryDynamicSgsReport3 {
    pub dynamic: DynamicSgsStepReport,
    pub wall_momentum: Option<WallMomentumExchange3>,
}

impl BoundaryGrid3 {
    // The shared filter only admits complete periodic boxes or planar channels.
    // Recheck each preparation because the workspace can outlive topology edits.
    fn dynamic_sgs_geometry(
        &self,
        gate: &CancelGate,
    ) -> Result<(Option<usize>, [[f64; 3]; 2]), BoundaryDynamicSgsError3> {
        checkpoint(gate)?;
        let mut normal = None;
        let mut walls = [[0.0; 3]; 2];
        for axis in 0..3 {
            match (
                self.boundaries.face(Face3::ALL[2 * axis]),
                self.boundaries.face(Face3::ALL[2 * axis + 1]),
            ) {
                (FaceBoundary3::Periodic, FaceBoundary3::Periodic) => {}
                (FaceBoundary3::Wall { velocity: low }, FaceBoundary3::Wall { velocity: high })
                    if normal.is_none() && low[axis] == 0.0 && high[axis] == 0.0 =>
                {
                    normal = Some(axis);
                    walls = [low, high];
                }
                _ => {
                    return Err(DynamicSgsError::Invalid(
                        "dynamic boundary SGS requires a periodic box or one impermeable planar wall pair",
                    )
                    .into());
                }
            }
        }
        for solid in &self.solid {
            checkpoint(gate)?;
            if *solid != 0 {
                return Err(DynamicSgsError::Invalid(
                    "dynamic boundary SGS does not admit interior solids",
                )
                .into());
            }
        }
        Ok((normal, walls))
    }

    /// Allocate a bounded dynamic workspace for this grid's actual topology.
    /// Open faces, interior obstacles and multiple wall pairs are refused.
    /// The filter uses unit lattice spacing; no external length is invented.
    pub fn dynamic_sgs_workspace(
        &self,
        gate: &CancelGate,
        maximum_bytes: usize,
    ) -> Result<BoundaryDynamicSgsWorkspace3, BoundaryDynamicSgsError3> {
        let (normal, _) = self.dynamic_sgs_geometry(gate)?;
        let count = self.f[0].len() * TILE_CELLS;
        let required = BoundaryDynamicSgsWorkspace3::required_bytes(count)?;
        if required > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required,
                maximum: maximum_bytes,
            }
            .into());
        }
        let filter = match normal {
            Some(axis) => DynamicFilter::Channel(ChannelDynamicSgs::new(
                self.dimensions(),
                1.0,
                axis,
                maximum_bytes,
            )?),
            None => DynamicFilter::Periodic(PeriodicDynamicSgs::new(
                self.dimensions(),
                1.0,
                maximum_bytes,
            )?),
        };
        let mut density = Vec::new();
        density
            .try_reserve_exact(count)
            .map_err(|_| DynamicSgsError::Allocation)?;
        density.resize(count, 0.0);
        let mut velocity = Vec::new();
        velocity
            .try_reserve_exact(count)
            .map_err(|_| DynamicSgsError::Allocation)?;
        velocity.resize(count, [0.0; 3]);
        let workspace = BoundaryDynamicSgsWorkspace3 {
            local: self.sgs_workspace()?,
            density,
            velocity,
            filter,
        };
        let actual = workspace.workspace_bytes();
        if actual > maximum_bytes {
            return Err(DynamicSgsError::WorkspaceBytes {
                required: actual,
                maximum: maximum_bytes,
            }
            .into());
        }
        checkpoint(gate)?;
        Ok(workspace)
    }

    /// Fit a Favre Germano–Lilly closure from the current half-force velocity
    /// and, for a channel, both actual moving-wall velocities. The shared fit
    /// reports backscatter suppression and refuses excessive coefficients;
    /// this owner additionally refuses excessive eddy/molecular viscosity.
    ///
    /// The original uniform molecular collision is the baseline at every call,
    /// never a previously augmented field. D3Q19 retains its existing tied
    /// second-order shear/bulk rate and independent higher-moment rates. This
    /// is an exploratory LES path, not a wall model or DNS calibration.
    pub fn prepare_dynamic_sgs(
        &self,
        gate: &CancelGate,
        config: BoundaryDynamicSgsConfig3,
        workspace: &mut BoundaryDynamicSgsWorkspace3,
    ) -> Result<DynamicSgsStepReport, BoundaryDynamicSgsError3> {
        let (normal, walls) = self.dynamic_sgs_geometry(gate)?;
        if workspace.filter.shape() != self.dimensions() || workspace.filter.normal_axis() != normal
        {
            return Err(BoundarySgsError3::Workspace.into());
        }
        let maximum_ratio = match (config, &workspace.filter) {
            (BoundaryDynamicSgsConfig3::Periodic(config), DynamicFilter::Periodic(_)) => {
                config.validate()?;
                config.maximum_eddy_to_molecular
            }
            (BoundaryDynamicSgsConfig3::Channel(config), DynamicFilter::Channel(_)) => {
                config.validate()?;
                config.maximum_eddy_to_molecular
            }
            _ => {
                return Err(DynamicSgsError::Invalid(
                    "dynamic controls do not match boundary topology",
                )
                .into());
            }
        };
        self.prepare_sgs_macros(gate, &mut workspace.local)?;
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let [x, y, z] = self.coords(tile, lane);
                let flat = (z * self.ny + y) * self.nx + x;
                let native = tile * TILE_CELLS + lane;
                workspace.density[flat] = workspace.local.density[native];
                workspace.velocity[flat] = workspace.local.velocity[native];
            }
        }
        let fit = match (config, &mut workspace.filter) {
            (BoundaryDynamicSgsConfig3::Periodic(config), DynamicFilter::Periodic(filter)) => {
                filter.prepare(config, &workspace.density, &workspace.velocity, || {
                    gate.is_requested()
                })?
            }
            (BoundaryDynamicSgsConfig3::Channel(config), DynamicFilter::Channel(filter)) => filter
                .prepare(
                    config,
                    &workspace.density,
                    &workspace.velocity,
                    walls,
                    || gate.is_requested(),
                )?,
            _ => unreachable!("configuration/topology admitted above"),
        };
        let nu = self.collision_model.kinematic_viscosity();
        let viscosities = workspace.filter.viscosities();
        let mut dissipation = SgsReport::default();
        for tile in 0..self.f[0].len() {
            checkpoint(gate)?;
            for lane in 0..TILE_CELLS {
                let cell = tile * TILE_CELLS + lane;
                let coordinate = self.coords(tile, lane);
                let [x, y, z] = coordinate;
                let eddy = viscosities[(z * self.ny + y) * self.nx + x];
                let ratio = eddy / nu;
                if !ratio.is_finite() || ratio > maximum_ratio {
                    return Err(BoundarySgsError3::from(SgsError::ViscosityRatioExceeded {
                        cell,
                        ratio,
                        maximum: maximum_ratio,
                    })
                    .into());
                }
                workspace.local.collisions[cell] = if eddy == 0.0 {
                    self.collision_model
                } else {
                    D3q19::with_shear_relaxation(self.collision_model, 0.5 + 3.0 * (nu + eddy))
                        .map_err(|source| BoundaryStepError3::Collision {
                            cell: coordinate,
                            source,
                        })?
                };
                let gradient = self.sgs_gradient(coordinate, &workspace.local.velocity);
                let dilation = (gradient[0][0] + gradient[1][1] + gradient[2][2]) / 3.0;
                let mut strain2 = 0.0;
                for i in 0..3 {
                    for j in 0..3 {
                        let s = 0.5 * (gradient[i][j] + gradient[j][i])
                            - if i == j { dilation } else { 0.0 };
                        strain2 += s * s;
                    }
                }
                let factor = 2.0 * workspace.local.density[cell] * strain2;
                dissipation.fluid_cells += 1;
                dissipation.resolved_dissipation_lu += factor * nu;
                dissipation.sgs_dissipation_lu += factor * eddy;
                dissipation.maximum_eddy_viscosity_lu =
                    dissipation.maximum_eddy_viscosity_lu.max(eddy);
                dissipation.maximum_eddy_to_molecular =
                    dissipation.maximum_eddy_to_molecular.max(ratio);
            }
        }
        checkpoint(gate)?;
        if !dissipation.resolved_dissipation_lu.is_finite()
            || !dissipation.sgs_dissipation_lu.is_finite()
        {
            return Err(BoundarySgsError3::from(SgsError::NonFiniteViscosity).into());
        }
        Ok(DynamicSgsStepReport { fit, dissipation })
    }

    /// Native serial collision/streaming with the dynamically fitted closure.
    /// Optional wall loads use the actual augmented post-collision populations
    /// before the sole stream/publication boundary, exactly as local SGS does.
    pub fn step_with_dynamic_sgs(
        &mut self,
        gate: &CancelGate,
        config: BoundaryDynamicSgsConfig3,
        workspace: &mut BoundaryDynamicSgsWorkspace3,
        measurement: Option<(&[BoundaryLink3], [f64; 3])>,
    ) -> Result<BoundaryDynamicSgsReport3, BoundaryDynamicSgsError3> {
        if let Some((links, origin)) = measurement {
            self.validate_wall_measurement(links, origin)?;
        }
        let dynamic = self.prepare_dynamic_sgs(gate, config, workspace)?;
        let wall_momentum =
            self.step_sgs_prepared_measured(gate, &workspace.local, |grid| match measurement {
                Some((links, origin)) => grid.wall_momentum(links, origin, Some(gate)).map(Some),
                None => Ok(None),
            })?;
        Ok(BoundaryDynamicSgsReport3 {
            dynamic,
            wall_momentum,
        })
    }

    /// Pooled native collision/streaming with the same canonical fit and load
    /// reduction as the serial path, independent of worker count.
    pub fn step_pooled_with_dynamic_sgs<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        config: BoundaryDynamicSgsConfig3,
        workspace: &mut BoundaryDynamicSgsWorkspace3,
        measurement: Option<(&[BoundaryLink3], [f64; 3])>,
    ) -> Result<BoundaryDynamicSgsReport3, BoundaryDynamicSgsError3> {
        if let Some((links, origin)) = measurement {
            self.validate_wall_measurement(links, origin)?;
        }
        let dynamic = self.prepare_dynamic_sgs(gate, config, workspace)?;
        let exchange =
            self.step_pooled_inner(pool, gate, measurement, Some(&workspace.local.collisions))?;
        Ok(BoundaryDynamicSgsReport3 {
            dynamic,
            wall_momentum: measurement.map(|_| exchange),
        })
    }
}
