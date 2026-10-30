//! Link-wise boundary conditions for the D3Q19 tile layout.
//!
//! This module deliberately keeps [`super::Duct`]'s frozen bit surface as the
//! body-force/periodic fixture, while
//! [`BoundaryGrid3`] owns the more general boundary semantics.  Solid links
//! use halfway bounce-back, moving planar walls add the standard momentum
//! correction, and open planar faces use a regularized non-equilibrium stress
//! reconstruction.  The latter is face-generic and preserves the prescribed
//! density/velocity moments without maintaining six hand-specialized D3Q19
//! formula tables.

use core::ops::ControlFlow;
use std::sync::Mutex;

use fs_exec::{CancelGate, Cancelled, Cx, KernelRunner, Reduce, TileKernel, TilePlan};

use super::{
    CollisionError3, CollisionModel3, E3, OPP3, Q3, TILE, TILE_CELLS, Tile, W3, collide_cell3,
    equilibrium3,
};
use crate::CS2;

mod sgs;
pub use sgs::{BoundarySgsError3, BoundarySgsReport3};
mod profile;
pub use profile::{BoundaryProfileError3, InletFlux3, VelocityProfileReport3};
mod partial_saturation;
pub use partial_saturation::{
    LubricatedPlanesExchange3, LubricationPlane3, PartialSaturationBodiesExchange3,
    PartialSaturationBody3, PartialSaturationBodyExchange3, PartialSaturationBodyResponse3,
    PartialSaturationError3, PartialSaturationExchange3, PartialSaturationForcing3,
    PlaneLubricationExchange3, SphereLubrication3, SphereLubricationExchange3,
};
mod body_motion;
pub use body_motion::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BodyFrameShape3, MovingLubricatedPlanesExchange3,
    MovingPartialSaturationBodiesExchange3, MovingPartialSaturationBody3,
    MovingPartialSaturationExchange3,
};

/// Bit-semantics version for the D3Q19 boundary surface.
///
/// This covers face ordering, link-mask construction, solid-cell
/// voxelization, moving-wall correction, regularized open-face
/// reconstruction, and deterministic traversal. Bump it whenever any of
/// those rules can move result bits.
pub const D3Q19_BOUNDARY_BIT_SEMANTICS_VERSION: u32 = 1;

/// Axis-aligned domain face in deterministic registry order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum Face3 {
    /// `x = 0` face.
    XMin = 0,
    /// `x = nx - 1` face.
    XMax = 1,
    /// `y = 0` face.
    YMin = 2,
    /// `y = ny - 1` face.
    YMax = 3,
    /// `z = 0` face.
    ZMin = 4,
    /// `z = nz - 1` face.
    ZMax = 5,
}

impl Face3 {
    /// All faces in the order used by boundary reconstruction.
    pub const ALL: [Face3; 6] = [
        Face3::XMin,
        Face3::XMax,
        Face3::YMin,
        Face3::YMax,
        Face3::ZMin,
        Face3::ZMax,
    ];

    #[inline]
    const fn index(self) -> usize {
        self as usize
    }

    #[inline]
    const fn axis(self) -> usize {
        self.index() / 2
    }

    #[inline]
    const fn is_min(self) -> bool {
        self.index().is_multiple_of(2)
    }

    #[inline]
    const fn opposite(self) -> Face3 {
        match self {
            Face3::XMin => Face3::XMax,
            Face3::XMax => Face3::XMin,
            Face3::YMin => Face3::YMax,
            Face3::YMax => Face3::YMin,
            Face3::ZMin => Face3::ZMax,
            Face3::ZMax => Face3::ZMin,
        }
    }

    /// Outward unit normal.
    #[must_use]
    pub const fn normal(self) -> [i32; 3] {
        let sign = if self.is_min() { -1 } else { 1 };
        match self.axis() {
            0 => [sign, 0, 0],
            1 => [0, sign, 0],
            _ => [0, 0, sign],
        }
    }
}

/// Boundary rule attached to one axis-aligned face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FaceBoundary3 {
    /// Wrap pulls to the opposite face. Periodicity must be paired on an axis.
    Periodic,
    /// Halfway bounce-back. `velocity` must be finite and tangential.
    Wall {
        /// Wall velocity in lattice units.
        velocity: [f64; 3],
    },
    /// Prescribed on-site velocity with density extrapolated from the first
    /// interior cell and regularized non-equilibrium stress reconstruction.
    Velocity {
        /// Target velocity in lattice units.
        velocity: [f64; 3],
    },
    /// Prescribed on-site density (isothermal pressure `p = c_s² rho`) with
    /// velocity extrapolated from the first interior cell.
    Pressure {
        /// Positive finite target density.
        density: f64,
    },
}

impl FaceBoundary3 {
    /// A stationary halfway bounce-back wall.
    #[must_use]
    pub const fn stationary_wall() -> FaceBoundary3 {
        FaceBoundary3::Wall { velocity: [0.0; 3] }
    }

    #[inline]
    const fn is_periodic(self) -> bool {
        matches!(self, FaceBoundary3::Periodic)
    }

    #[inline]
    const fn is_open(self) -> bool {
        matches!(
            self,
            FaceBoundary3::Velocity { .. } | FaceBoundary3::Pressure { .. }
        )
    }
}

/// Complete boundary specification for the six domain faces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BoundarySpec3 {
    faces: [FaceBoundary3; 6],
}

impl BoundarySpec3 {
    /// Construct from face rules in [`Face3::ALL`] order.
    #[must_use]
    pub const fn new(faces: [FaceBoundary3; 6]) -> BoundarySpec3 {
        BoundarySpec3 { faces }
    }

    /// Fully periodic domain.
    #[must_use]
    pub const fn periodic() -> BoundarySpec3 {
        BoundarySpec3::new([FaceBoundary3::Periodic; 6])
    }

    /// Square duct: stationary x/y walls and periodic z.
    #[must_use]
    pub const fn periodic_duct_z() -> BoundarySpec3 {
        let wall = FaceBoundary3::stationary_wall();
        BoundarySpec3::new([
            wall,
            wall,
            wall,
            wall,
            FaceBoundary3::Periodic,
            FaceBoundary3::Periodic,
        ])
    }

    /// Square duct with regularized velocity inlet at z-min and pressure
    /// outlet at z-max.
    #[must_use]
    pub const fn velocity_pressure_duct_z(
        inlet_velocity: [f64; 3],
        outlet_density: f64,
    ) -> BoundarySpec3 {
        let wall = FaceBoundary3::stationary_wall();
        BoundarySpec3::new([
            wall,
            wall,
            wall,
            wall,
            FaceBoundary3::Velocity {
                velocity: inlet_velocity,
            },
            FaceBoundary3::Pressure {
                density: outlet_density,
            },
        ])
    }

    /// Pressure-driven square duct with independently prescribed densities on
    /// the z faces.
    #[must_use]
    pub const fn pressure_duct_z(inlet_density: f64, outlet_density: f64) -> BoundarySpec3 {
        let wall = FaceBoundary3::stationary_wall();
        BoundarySpec3::new([
            wall,
            wall,
            wall,
            wall,
            FaceBoundary3::Pressure {
                density: inlet_density,
            },
            FaceBoundary3::Pressure {
                density: outlet_density,
            },
        ])
    }

    /// Closed cavity with a moving y-max lid. The velocity is validated as
    /// tangential when the grid is constructed.
    #[must_use]
    pub const fn lid_cavity(lid_velocity: [f64; 3]) -> BoundarySpec3 {
        let wall = FaceBoundary3::stationary_wall();
        BoundarySpec3::new([
            wall,
            wall,
            wall,
            FaceBoundary3::Wall {
                velocity: lid_velocity,
            },
            wall,
            wall,
        ])
    }

    /// Rule for one face.
    #[must_use]
    pub const fn face(self, face: Face3) -> FaceBoundary3 {
        self.faces[face.index()]
    }

    fn validate(self) {
        for min_face in [Face3::XMin, Face3::YMin, Face3::ZMin] {
            let max_face = min_face.opposite();
            assert_eq!(
                self.face(min_face).is_periodic(),
                self.face(max_face).is_periodic(),
                "periodic boundaries must be paired on the {min_face:?}/{max_face:?} axis"
            );
        }

        let mut open_axes = [false; 3];
        for face in Face3::ALL {
            match self.face(face) {
                FaceBoundary3::Periodic => {}
                FaceBoundary3::Wall { velocity } => {
                    assert!(
                        velocity.iter().all(|component| component.is_finite()),
                        "wall velocity on {face:?} must be finite"
                    );
                    let normal = face.normal();
                    let normal_velocity = velocity
                        .iter()
                        .zip(normal)
                        .map(|(u, n)| *u * f64::from(n))
                        .sum::<f64>();
                    assert!(
                        normal_velocity.abs() <= 16.0 * f64::EPSILON,
                        "wall velocity on {face:?} must be tangential"
                    );
                    let speed_sq = velocity.iter().map(|u| u * u).sum::<f64>();
                    assert!(
                        speed_sq < 0.03,
                        "wall velocity on {face:?} exceeds the low-Mach admission envelope"
                    );
                }
                FaceBoundary3::Velocity { velocity } => {
                    assert!(
                        velocity.iter().all(|component| component.is_finite()),
                        "velocity boundary on {face:?} must be finite"
                    );
                    let speed_sq = velocity.iter().map(|u| u * u).sum::<f64>();
                    assert!(
                        speed_sq < 0.03,
                        "velocity boundary on {face:?} exceeds the low-Mach admission envelope"
                    );
                    open_axes[face.axis()] = true;
                }
                FaceBoundary3::Pressure { density } => {
                    assert!(
                        density.is_finite() && density > 0.0,
                        "pressure-boundary density on {face:?} must be positive and finite"
                    );
                    open_axes[face.axis()] = true;
                }
            }
        }
        assert!(
            open_axes.into_iter().filter(|is_open| *is_open).count() <= 1,
            "open velocity/pressure faces may occupy only one axis per grid"
        );
    }
}

/// One aligned tile of per-cell D3Q19 wall-link masks.
///
/// Bit `q` is set when pull direction `q` crosses either a wall face or a
/// voxelized solid cell. Bit zero is never set because the rest population
/// does not cross a link.
#[derive(Clone, PartialEq, Eq)]
#[repr(align(128))]
pub struct LinkMaskTile3([u32; TILE_CELLS]);

impl LinkMaskTile3 {
    fn empty() -> LinkMaskTile3 {
        LinkMaskTile3([0; TILE_CELLS])
    }

    /// Masks in canonical local-lane order (x-fastest, then y, then z).
    #[must_use]
    pub const fn as_array(&self) -> &[u32; TILE_CELLS] {
        &self.0
    }
}

impl core::fmt::Debug for LinkMaskTile3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("LinkMaskTile3")
            .field(&self.0.as_slice())
            .finish()
    }
}

/// Canonically ordered wall link exposed for sparse-tile bookkeeping.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BoundaryLink3 {
    /// Tile index in x-fastest tile-major order.
    pub tile: usize,
    /// Lane within the 4x4x4 tile.
    pub lane: usize,
    /// Global cell coordinate.
    pub cell: [usize; 3],
    /// D3Q19 direction index.
    pub direction: usize,
}

/// Force, torque, and work transferred to selected D3Q19 wall links in one
/// lattice step. All quantities use lattice units (`dt = dx = 1`).
///
/// The boundary-relative impulse is `(c_out-u_wall) f_out -
/// (c_in-u_wall) f_in` (Wen et al., 2014,
/// <https://arxiv.org/abs/1303.0625>). The separate population and moving-mass
/// terms close `wall_impulse + fluid_population_impulse =
/// wall_velocity_mass_impulse`; the same identity holds for angular impulse.
/// Thus moving mass is not mistaken for an equal-and-opposite contact force.
/// Torque uses the halfway intersection measured from the requested origin;
/// cell centers have coordinates `[x+0.5, y+0.5, z+0.5]`, as in voxelization.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[must_use]
pub struct WallMomentumExchange3 {
    /// Hydrodynamic impulse delivered to the selected walls.
    pub wall_impulse: [f64; 3],
    /// Change in fluid-population momentum across the selected wall links.
    pub fluid_population_impulse: [f64; 3],
    /// Incoming-minus-outgoing population mass across the selected links.
    pub fluid_mass_change: f64,
    /// Sum of `u_wall * (f_in - f_out)` across the selected links.
    pub wall_velocity_mass_impulse: [f64; 3],
    /// Wall angular impulse about the requested origin.
    pub wall_angular_impulse: [f64; 3],
    /// Angular impulse of the fluid-population momentum change.
    pub fluid_population_angular_impulse: [f64; 3],
    /// Angular impulse of the moving-mass term.
    pub wall_velocity_mass_angular_impulse: [f64; 3],
    /// Work delivered to walls: sum of link impulse dotted with wall velocity.
    pub wall_work: f64,
    /// Number of measured links; each selected link is counted once.
    pub measured_links: usize,
}

impl WallMomentumExchange3 {
    fn add_link(
        &mut self,
        link: BoundaryLink3,
        outgoing: f64,
        rho: f64,
        velocity: [f64; 3],
        origin: [f64; 3],
    ) {
        let q = link.direction;
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let incoming = wall_population(outgoing, rho, q, velocity);
        let mass_change = incoming - outgoing;
        let fluid = e.map(|component| component * (outgoing + incoming));
        let moving_mass = velocity.map(|component| component * mass_change);
        let wall = core::array::from_fn::<_, 3, _>(|axis| {
            (-e[axis] - velocity[axis]) * outgoing - (e[axis] - velocity[axis]) * incoming
        });
        let offset = core::array::from_fn::<_, 3, _>(|axis| {
            link.cell[axis] as f64 + 0.5 - 0.5 * e[axis] - origin[axis]
        });
        for axis in 0..3 {
            let b = (axis + 1) % 3;
            let c = (axis + 2) % 3;
            self.wall_impulse[axis] += wall[axis];
            self.fluid_population_impulse[axis] += fluid[axis];
            self.wall_velocity_mass_impulse[axis] += moving_mass[axis];
            self.wall_angular_impulse[axis] += offset[b] * wall[c] - offset[c] * wall[b];
            self.fluid_population_angular_impulse[axis] +=
                offset[b] * fluid[c] - offset[c] * fluid[b];
            self.wall_velocity_mass_angular_impulse[axis] +=
                offset[b] * moving_mass[c] - offset[c] * moving_mass[b];
            self.wall_work += wall[axis] * velocity[axis];
        }
        self.fluid_mass_change += mass_change;
        self.measured_links += 1;
    }

    fn checked(self) -> Result<Self, BoundaryStepError3> {
        if ![
            self.wall_impulse,
            self.fluid_population_impulse,
            self.wall_velocity_mass_impulse,
            self.wall_angular_impulse,
            self.fluid_population_angular_impulse,
            self.wall_velocity_mass_angular_impulse,
        ]
        .into_iter()
        .flatten()
        .chain([self.fluid_mass_change, self.wall_work])
        .all(f64::is_finite)
        {
            return Err(BoundaryStepError3::WallMomentum {
                reason: "wall momentum exchange is non-finite",
            });
        }
        Ok(self)
    }
}

/// Inertia and applied impulse of one translating planar wall for one lattice
/// step. Supply the minimum face first, then the opposite maximum face.
///
/// The wall velocity is the corresponding [`FaceBoundary3::Wall`] velocity
/// stored in the grid. Mass is the total mass of the complete periodic face,
/// in the same lattice mass units as the sum of fluid densities (`dx = dt = 1`).
/// The wall has two tangential translation degrees of freedom; its normal
/// position and orientation are constrained.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiniteMassWall3 {
    /// Positive finite total wall mass, not an areal density.
    pub mass: f64,
    /// Applied impulse at the geometric face center during this step. A normal
    /// component is carried by the support reaction and does no work.
    pub external_impulse: [f64; 3],
}

/// Accepted two-way fluid/wall exchange for one finite-mass planar wall.
/// All quantities use lattice units. Angular impulses use the moment origin
/// supplied to the step; external impulse acts at the geometric face center.
///
/// Tangential wall momentum changes by the hydrodynamic plus external
/// impulse. Normal and rotational support reactions are reported explicitly.
/// The hydrodynamic and external work sum to the wall kinetic-energy change
/// to roundoff because the boundary uses the wall's midpoint velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
#[must_use]
pub struct FiniteMassWallReceipt3 {
    /// Face carrying this wall.
    pub face: Face3,
    /// Total wall mass used in the solve.
    pub mass: f64,
    /// Wall velocity before the step.
    pub velocity_before: [f64; 3],
    /// Wall velocity used by bounce-back during the step.
    pub velocity_midpoint: [f64; 3],
    /// Wall velocity published in the boundary specification after the step.
    pub velocity_after: [f64; 3],
    /// Actual link-wise fluid impulse, torque, mass exchange, and work.
    pub hydrodynamic: WallMomentumExchange3,
    /// Applied impulse at the geometric face center.
    pub external_impulse: [f64; 3],
    /// Support impulse on the wall, canceling the total normal impulse.
    pub constraint_impulse: [f64; 3],
    /// Support angular impulse on the wall: the moment of the normal support
    /// force plus the couple that prevents rotation about the face center.
    pub constraint_angular_impulse: [f64; 3],
    /// Change in the wall's tangential translational kinetic energy.
    pub kinetic_energy_change: f64,
    /// Applied impulse dotted with midpoint velocity.
    pub external_work: f64,
}

struct FiniteMassWallProposal3 {
    midpoint_boundaries: BoundarySpec3,
    final_boundaries: BoundarySpec3,
    receipts: [FiniteMassWallReceipt3; 2],
}

enum PullSource {
    Fluid { tile: usize, lane: usize },
    Wall { voxel: bool },
    Open,
}

/// General D3Q19 grid with deterministic link-wise boundary bookkeeping.
pub struct BoundaryGrid3 {
    nx: usize,
    ny: usize,
    nz: usize,
    collision_model: CollisionModel3,
    force: [f64; 3],
    boundaries: BoundarySpec3,
    f: [Vec<Tile>; Q3],
    post: [Vec<Tile>; Q3],
    /// One 64-bit solid occupancy word per 4x4x4 tile.
    solid: Vec<u64>,
    /// Wall links. At a wall/open seam, directions that cross a wall remain
    /// wall-owned while pure open-face directions remain reconstructable.
    link_masks: Vec<LinkMaskTile3>,
    /// Open links on face-interior cells reconstructed after streaming.
    open_link_masks: Vec<LinkMaskTile3>,
    /// Wall bits whose owner is a stationary voxel obstacle or a stationary
    /// wall/open seam, even if the same cell touches a moving exterior wall.
    stationary_link_masks: Vec<LinkMaskTile3>,
    /// Solid topology is immutable after it is defined or state evolution
    /// begins; this prevents silent mass/topology transitions.
    topology_locked: bool,
    /// Precompiled fluid pull sources: entry `(tile * TILE_CELLS + lane) *
    /// Q3 + q` is the flat source cell `tile * TILE_CELLS + lane` of
    /// direction `q`, or [`NO_FLUID_SOURCE`] for wall/open links and solid
    /// cells. A pure function of topology, compiled with the link masks, so
    /// streaming never re-derives neighbours (bit-identical by construction:
    /// it selects the same source the classifier would).
    pull: Vec<u32>,
    /// Stream destination for [`BoundaryGrid3::step_pooled`], allocated on
    /// first use as a copy of `f` (so never-written solid entries agree) and
    /// swapped with `f` only after a complete pass: a cancelled or refused
    /// pooled step publishes nothing.
    next: Option<[Vec<Tile>; Q3]>,
    /// Optional per-cell Guo body force added to the uniform `force`
    /// (flat index `tile * TILE_CELLS + lane`). `None` leaves every result
    /// bit-identical to the uniform-force grid.
    force_field: Option<Vec<[f64; 3]>>,
    /// Open-face cells in canonical reconstruction order (face, then z, y,
    /// x) with the interior cell each reads. Cached on the first step; the
    /// topology is immutable from then on, so the per-step pass visits only
    /// these cells instead of scanning the whole grid.
    open_cells: Option<Vec<([usize; 3], [usize; 3], Face3)>>,
    /// Explicit spatial velocity traces, sampled outside the lattice step.
    /// Each face owns two reusable plane buffers for atomic profile updates.
    velocity_profiles: [Option<profile::ProfileState3>; 6],
}

/// Lattice tiles per pooled kernel tile (512 cells).
const POOLED_GROUP_TILES: usize = 8;

/// Refusal from a pooled or finite-mass-wall [`BoundaryGrid3`] step. The
/// published state is the pre-step state in every case.
#[derive(Debug, Clone, PartialEq)]
pub enum BoundaryStepError3 {
    /// The cancel gate tripped during the collide or stream pass.
    Cancelled,
    /// The first (lowest canonical index) cell whose collision refused.
    Collision {
        /// Cell coordinates.
        cell: [usize; 3],
        /// The collision diagnostic.
        source: CollisionError3,
    },
    /// Pool-level refusal (worker fault or admission), rendered.
    Pool(String),
    /// A cell used by open-face reconstruction or finite-wall coupling has
    /// non-positive or non-finite density: the lattice left its physical regime.
    Unphysical {
        /// Cell coordinates.
        cell: [usize; 3],
    },
    /// Invalid wall-link selection/origin, or non-finite momentum exchange.
    WallMomentum {
        /// Stable refusal detail.
        reason: &'static str,
    },
    /// Inadmissible finite wall mass, impulse, geometry, or solved velocity.
    FiniteMassWall {
        /// Stable refusal detail.
        reason: &'static str,
    },
}

impl core::fmt::Display for BoundaryStepError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => write!(f, "pooled D3Q19 step cancelled"),
            Self::Collision { cell, source } => {
                write!(f, "D3Q19 collision refused at cell {cell:?}: {source:?}")
            }
            Self::Pool(detail) => write!(f, "pooled D3Q19 step refused by the pool: {detail}"),
            Self::Unphysical { cell } => {
                write!(
                    f,
                    "D3Q19 boundary source cell {cell:?} lost positive finite density"
                )
            }
            Self::WallMomentum { reason } => write!(f, "D3Q19 wall momentum refused: {reason}"),
            Self::FiniteMassWall { reason } => {
                write!(f, "D3Q19 finite-mass wall coupling refused: {reason}")
            }
        }
    }
}

impl std::error::Error for BoundaryStepError3 {}

/// Pull-table sentinel for links whose source is not a fluid cell.
const NO_FLUID_SOURCE: u32 = u32::MAX;

impl BoundaryGrid3 {
    /// Unit-density fluid at rest under the supplied face rules.
    ///
    /// Every dimension must be a positive multiple of [`TILE`]. `tau` must
    /// be finite and greater than 0.5, and `force` must be finite.
    ///
    /// # Panics
    /// Panics when dimensions, allocation geometry, relaxation time, force,
    /// or boundary parameters are inadmissible.
    #[must_use]
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        tau: f64,
        force: [f64; 3],
        boundaries: BoundarySpec3,
    ) -> BoundaryGrid3 {
        Self::with_collision_model(nx, ny, nz, CollisionModel3::Bgk { tau }, force, boundaries)
    }

    /// Unit-density fluid at rest under an explicit collision model.
    ///
    /// Every dimension must be a positive multiple of [`TILE`]. The collision
    /// model and force must be finite and admissible. BGK and moment-space
    /// models support body forces on periodic/walled domains; regularized
    /// open velocity/pressure faces still require a force-free domain.
    ///
    /// # Panics
    /// Panics when dimensions, allocation geometry, collision parameters,
    /// force, or boundary parameters are inadmissible.
    #[must_use]
    pub fn with_collision_model(
        nx: usize,
        ny: usize,
        nz: usize,
        collision_model: CollisionModel3,
        force: [f64; 3],
        boundaries: BoundarySpec3,
    ) -> BoundaryGrid3 {
        assert!(
            nx > 0
                && ny > 0
                && nz > 0
                && nx.is_multiple_of(TILE)
                && ny.is_multiple_of(TILE)
                && nz.is_multiple_of(TILE),
            "grid dimensions must be positive multiples of {TILE} (got {nx}x{ny}x{nz})"
        );
        collision_model
            .validate()
            .expect("D3Q19 collision model must be physically admissible");
        assert!(
            force.iter().all(|component| component.is_finite()),
            "body force must be finite"
        );
        assert!(
            collision_model.supports_body_force()
                || force
                    .iter()
                    .all(|value| matches!(value.classify(), core::num::FpCategory::Zero)),
            "the selected D3Q19 collision model does not admit body forcing"
        );
        boundaries.validate();
        let has_open_face = Face3::ALL
            .into_iter()
            .any(|face| boundaries.face(face).is_open());
        assert!(
            !has_open_face
                || force
                    .iter()
                    .all(|value| matches!(value.classify(), core::num::FpCategory::Zero)),
            "regularized velocity/pressure boundaries currently require zero body force"
        );
        let ntx = nx / TILE;
        let nty = ny / TILE;
        let ntz = nz / TILE;
        let tiles = ntx
            .checked_mul(nty)
            .and_then(|xy| xy.checked_mul(ntz))
            .expect("grid tile count overflow");
        let f0 = equilibrium3(1.0, [0.0; 3]);
        let f = core::array::from_fn(|i| vec![Tile::filled(f0[i]); tiles]);
        let post = core::array::from_fn(|i| vec![Tile::filled(f0[i]); tiles]);
        let mut grid = BoundaryGrid3 {
            nx,
            ny,
            nz,
            collision_model,
            force,
            boundaries,
            f,
            post,
            solid: vec![0; tiles],
            link_masks: vec![LinkMaskTile3::empty(); tiles],
            open_link_masks: vec![LinkMaskTile3::empty(); tiles],
            stationary_link_masks: vec![LinkMaskTile3::empty(); tiles],
            topology_locked: false,
            pull: Vec::new(),
            next: None,
            force_field: None,
            open_cells: None,
            velocity_profiles: core::array::from_fn(|_| None),
        };
        let (wall_masks, open_masks, stationary_masks) = grid.compile_link_masks(&grid.solid);
        grid.link_masks = wall_masks;
        grid.open_link_masks = open_masks;
        grid.stationary_link_masks = stationary_masks;
        grid.pull = grid.compile_pull_table(&grid.solid);
        grid
    }

    /// Grid dimensions in cells.
    #[must_use]
    pub const fn dimensions(&self) -> [usize; 3] {
        [self.nx, self.ny, self.nz]
    }

    /// Kinematic viscosity implied by the selected collision model.
    #[must_use]
    pub fn viscosity(&self) -> f64 {
        self.collision_model.kinematic_viscosity()
    }

    /// Collision model used for every fluid-cell collision.
    #[must_use]
    pub const fn collision_model(&self) -> CollisionModel3 {
        self.collision_model
    }

    /// Boundary specification used by the grid.
    #[must_use]
    pub const fn boundaries(&self) -> BoundarySpec3 {
        self.boundaries
    }

    /// Per-tile wall-link masks in canonical tile order.
    #[must_use]
    pub fn tile_link_masks(&self) -> &[LinkMaskTile3] {
        &self.link_masks
    }

    /// Per-tile open-link masks in canonical tile order. These are disjoint
    /// from [`BoundaryGrid3::tile_link_masks`]. At a wall/open rim, a
    /// direction that crosses both faces is wall-owned, while a direction
    /// that crosses only the open face remains in this mask.
    #[must_use]
    pub fn tile_open_link_masks(&self) -> &[LinkMaskTile3] {
        &self.open_link_masks
    }

    /// The link mask at one cell.
    #[must_use]
    pub fn link_mask(&self, x: usize, y: usize, z: usize) -> u32 {
        let (tile, lane) = self.addr_checked(x, y, z);
        self.link_masks[tile].0[lane]
    }

    /// Open velocity/pressure links at one face-interior cell.
    #[must_use]
    pub fn open_link_mask(&self, x: usize, y: usize, z: usize) -> u32 {
        let (tile, lane) = self.addr_checked(x, y, z);
        self.open_link_masks[tile].0[lane]
    }

    /// All wall links in deterministic tile, lane, direction order.
    #[must_use]
    pub fn boundary_links(&self) -> Vec<BoundaryLink3> {
        let mut links = Vec::new();
        for (tile, masks) in self.link_masks.iter().enumerate() {
            for lane in 0..TILE_CELLS {
                let mask = masks.0[lane];
                if mask == 0 {
                    continue;
                }
                let cell = self.coords(tile, lane);
                for direction in 1..Q3 {
                    if mask & (1u32 << direction) != 0 {
                        links.push(BoundaryLink3 {
                            tile,
                            lane,
                            cell,
                            direction,
                        });
                    }
                }
            }
        }
        links
    }

    /// Whether one lattice cell is voxelized solid.
    #[must_use]
    pub fn is_solid(&self, x: usize, y: usize, z: usize) -> bool {
        let (tile, lane) = self.addr_checked(x, y, z);
        self.solid[tile] & (1u64 << lane) != 0
    }

    /// Define immutable solid occupancy by sampling an SDF-like classifier at
    /// cell centers. Non-positive values are solid; all samples must be
    /// finite. Sampling, validation, and mask compilation happen against
    /// temporary storage, then commit atomically, so a caught panic cannot
    /// leave occupancy and runtime stencils inconsistent.
    ///
    /// This operation is initialization-only and may be called once, before a
    /// perturbation or time step. Dynamic topology requires a future explicit
    /// mass/topology transition receipt rather than mutating this grid.
    pub fn voxelize_sdf(&mut self, mut signed_distance: impl FnMut([f64; 3]) -> f64) {
        assert!(
            !self.topology_locked,
            "solid topology is immutable after initialization or state evolution"
        );
        let mut proposed_solid = vec![0u64; self.solid.len()];
        let mut fluid_cells = 0usize;
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    let sample = signed_distance([x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5]);
                    assert!(sample.is_finite(), "voxelized SDF samples must be finite");
                    let (tile, lane) = self.addr(x, y, z);
                    if sample <= 0.0 {
                        proposed_solid[tile] |= 1u64 << lane;
                    } else {
                        fluid_cells += 1;
                    }
                }
            }
        }
        assert!(
            fluid_cells > 0,
            "voxelized domain must retain at least one fluid cell"
        );
        self.validate_open_neighbors(&proposed_solid);
        let (wall_masks, open_masks, stationary_masks) = self.compile_link_masks(&proposed_solid);
        let pull = self.compile_pull_table(&proposed_solid);
        self.solid = proposed_solid;
        self.link_masks = wall_masks;
        self.open_link_masks = open_masks;
        self.stationary_link_masks = stationary_masks;
        self.pull = pull;
        self.topology_locked = true;
    }

    /// Density of a fluid cell.
    ///
    /// # Panics
    /// Panics for an out-of-range or solid cell.
    #[must_use]
    pub fn density(&self, x: usize, y: usize, z: usize) -> f64 {
        assert!(!self.is_solid(x, y, z), "solid cells have no fluid density");
        let (tile, lane) = self.addr(x, y, z);
        (0..Q3).map(|q| self.f[q][tile].0[lane]).sum()
    }

    /// Force-corrected velocity of a fluid cell.
    ///
    /// # Panics
    /// Panics for an out-of-range or solid cell.
    #[must_use]
    pub fn velocity(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        assert!(
            !self.is_solid(x, y, z),
            "solid cells have no fluid velocity"
        );
        let (rho, momentum) = self.raw_moments(x, y, z);
        let (tile, lane) = self.addr(x, y, z);
        let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
        core::array::from_fn(|axis| (momentum[axis] + 0.5 * force[axis]) / rho)
    }

    /// Total held body-force density (uniform plus per-cell, in lattice
    /// units) at a fluid cell. Physical momentum density is the first
    /// population moment plus half this force.
    ///
    /// # Panics
    /// Panics for an out-of-range or solid cell.
    #[must_use]
    pub fn force_density(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        assert!(!self.is_solid(x, y, z), "solid cells have no fluid force");
        let (tile, lane) = self.addr(x, y, z);
        cell_force(self.force, self.force_field.as_deref(), tile, lane)
    }

    /// Select the initial per-cell body-force frame (lattice units), ADDED
    /// to the uniform force. Samples every fluid cell, skips solids, and
    /// replaces any previous field without changing stored populations.
    ///
    /// This initializer deliberately changes the interpreted physical
    /// momentum by half the change in force. When replacing a held force
    /// on an evolved state, use [`Self::rebind_force_field`] to preserve
    /// physical momentum and stress instead.
    ///
    /// # Panics
    /// If the grid has an open velocity/pressure face (whose regularized
    /// reconstruction assumes a force-free fluid), the collision model does
    /// not admit body forcing, or a sample or its sum with the uniform force
    /// is non-finite. Refusal leaves the prior field and populations intact.
    pub fn set_force_field(&mut self, force: impl FnMut([usize; 3]) -> [f64; 3]) {
        self.force_field = Some(self.sample_force_field(force));
    }

    /// Replace a held per-cell body force without advancing time or adding
    /// an impulse. The uniform force remains part of each cell's total force.
    /// Use a zero field to switch off the per-cell contribution.
    ///
    /// Guo populations store a half-source shift: physical momentum is
    /// `sum_q(f_q c_q) + F/2`. Rebinding applies `S(u, (F_old-F_new)/2)` at
    /// the OLD force-corrected velocity. This preserves mass, physical
    /// momentum, and the full force-corrected nonequilibrium stress, rather
    /// than resetting the cell to equilibrium. It contributes no work.
    /// The next step then applies the new force for one lattice time step.
    ///
    /// All samples and proposed populations are checked before publication.
    /// The existing collision scratch holds the proposal, avoiding another
    /// lattice-sized population allocation. Unchanged forces preserve the
    /// corresponding population bits exactly; solids are never sampled.
    ///
    /// # Panics
    /// The same field-admission failures as [`Self::set_force_field`], or
    /// a changed cell with nonpositive/nonfinite density, nonfinite physical
    /// velocity, or a nonfinite population-frame transform. Even a callback
    /// panic or late-cell refusal leaves accepted populations and the held
    /// force field intact.
    pub fn rebind_force_field(&mut self, force: impl FnMut([usize; 3]) -> [f64; 3]) {
        let field = self.sample_force_field(force);
        for tile in 0..self.f[0].len() {
            for lane in 0..TILE_CELLS {
                let mut populations = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                if self.solid[tile] & (1u64 << lane) == 0 {
                    let previous = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                    let next = cell_force(self.force, Some(&field), tile, lane);
                    if previous != next {
                        let [x, y, z] = self.coords(tile, lane);
                        let (rho, momentum) = self.raw_moments(x, y, z);
                        let velocity = core::array::from_fn(|axis| {
                            (momentum[axis] + 0.5 * previous[axis]) / rho
                        });
                        assert!(
                            velocity.iter().all(|value| value.is_finite()),
                            "force rebinding requires finite physical velocity"
                        );
                        crate::thermal::forcing::reframe3(
                            &mut populations,
                            velocity,
                            previous,
                            next,
                        );
                        let density = populations.iter().sum::<f64>();
                        assert!(
                            populations.iter().all(|value| value.is_finite())
                                && density.is_finite()
                                && density > 0.0,
                            "force rebinding produced an inadmissible population frame"
                        );
                    }
                }
                for (proposed, value) in self.post.iter_mut().zip(populations) {
                    proposed[tile].0[lane] = value;
                }
            }
        }
        core::mem::swap(&mut self.f, &mut self.post);
        self.force_field = Some(field);
    }

    fn sample_force_field(&self, mut force: impl FnMut([usize; 3]) -> [f64; 3]) -> Vec<[f64; 3]> {
        assert!(
            !Face3::ALL
                .into_iter()
                .any(|face| self.boundaries.face(face).is_open()),
            "per-cell forcing requires a grid without open velocity/pressure faces"
        );
        assert!(
            self.collision_model.supports_body_force(),
            "the selected collision model does not admit body forcing"
        );
        let mut field = vec![[0.0; 3]; self.f[0].len() * TILE_CELLS];
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    if self.is_solid(x, y, z) {
                        continue;
                    }
                    let value = force([x, y, z]);
                    assert!(
                        value
                            .iter()
                            .zip(self.force)
                            .all(|(local, uniform)| local.is_finite()
                                && (local + uniform).is_finite()),
                        "per-cell force and total force must be finite"
                    );
                    let (tile, lane) = self.addr(x, y, z);
                    field[tile * TILE_CELLS + lane] = value;
                }
            }
        }
        field
    }

    /// Population vector of a fluid cell, copied in D3Q19 direction order.
    /// This is a read-only inspection surface for independent moment and
    /// stress conformance checks.
    ///
    /// # Panics
    /// Panics for an out-of-range or solid cell.
    #[must_use]
    pub fn populations(&self, x: usize, y: usize, z: usize) -> [f64; Q3] {
        assert!(
            !self.is_solid(x, y, z),
            "solid cells have no fluid populations"
        );
        let (tile, lane) = self.addr(x, y, z);
        core::array::from_fn(|q| self.f[q][tile].0[lane])
    }
    /// Overwrite one fluid cell's populations (fixture initialization;
    /// the write-side counterpart of [`BoundaryGrid3::populations`]).
    ///
    /// # Panics
    /// Panics for an out-of-range or solid cell, or any non-finite value.
    pub fn set_populations(&mut self, x: usize, y: usize, z: usize, values: &[f64; Q3]) {
        assert!(
            !self.is_solid(x, y, z),
            "solid cells have no fluid populations"
        );
        assert!(
            values.iter().all(|v| v.is_finite()),
            "set_populations refuses non-finite values"
        );
        let (tile, lane) = self.addr(x, y, z);
        for (field, value) in self.f.iter_mut().zip(values) {
            field[tile].0[lane] = *value;
        }
    }

    /// Total mass over fluid cells in canonical tile/lane order.
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        let mut total = 0.0;
        for tile in 0..self.f[0].len() {
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) == 0 {
                    for q in 0..Q3 {
                        total += self.f[q][tile].0[lane];
                    }
                }
            }
        }
        total
    }

    /// Deterministic density perturbation for replay/golden fixtures.
    pub fn perturb(&mut self, seed: u64, amplitude: f64) {
        assert!(
            amplitude.is_finite() && amplitude.abs() < 1.0,
            "perturbation amplitude magnitude must be finite and less than one"
        );
        self.topology_locked = true;
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    if self.is_solid(x, y, z) {
                        continue;
                    }
                    let mut h = seed
                        ^ (x as u64)
                            .wrapping_mul(0x9e37_79b9_7f4a_7c15)
                            .wrapping_add((y as u64).wrapping_mul(0xbf58_476d_1ce4_e5b9))
                            .wrapping_add((z as u64).wrapping_mul(0x94d0_49bb_1331_11eb));
                    h ^= h >> 30;
                    h = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
                    h ^= h >> 27;
                    let unit = (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0;
                    let rho = 1.0 + amplitude * unit;
                    assert!(rho > 0.0, "perturbation must retain positive density");
                    self.write_equilibrium(x, y, z, rho, [0.0; 3]);
                }
            }
        }
    }

    /// One selected collision, link-wise pull-stream, and open-face
    /// reconstruction step.
    pub fn step(&mut self) {
        self.topology_locked = true;
        self.collide();
        self.stream();
        self.apply_open_boundaries();
    }

    /// Advance the existing serial step and return the hydrodynamic loads on
    /// a selected wall or obstacle. Obtain links once from [`Self::boundary_links`]
    /// and retain any subset in its original order; immutable topology makes
    /// the selection reusable. Empty selections are valid.
    ///
    /// This reads only the selected links of the actual post-collision state,
    /// with no full-grid momentum scan and no change to stepping arithmetic.
    /// Open-boundary momentum flux and body forcing are separate exchanges.
    ///
    /// # Panics
    /// Invalid, duplicate, or out-of-order links, a non-finite origin, collision
    /// failure, or non-finite exchange. Invalid measurements are refused before
    /// populations advance, including a non-finite computed exchange.
    pub fn step_with_wall_momentum(
        &mut self,
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
    ) -> WallMomentumExchange3 {
        self.validate_wall_measurement(links, moment_origin)
            .expect("wall measurement requires valid ordered links and a finite origin");
        self.topology_locked = true;
        self.collide();
        let exchange = self
            .wall_momentum(links, moment_origin, None)
            .expect("wall momentum exchange must remain finite");
        self.stream();
        self.apply_open_boundaries();
        exchange
    }

    /// Advance the fluid and two finite-mass planar walls together. The grid
    /// must have exactly two opposite complete wall faces, periodic tangent
    /// directions, and no voxel solids. `walls` is in minimum/maximum face
    /// order. Initial velocities come from the boundary specification;
    /// accepted final velocities replace them for the next step.
    ///
    /// With the actual post-collision populations held fixed, each complete
    /// face has tangential impulse `J(u) = J(0) - D u`, where
    /// `D = sum_face(rho_post)/3`. The direct midpoint solve is
    /// `u_mid = (M u_old + (J(0) + J_external)/2)/(M + D/2)`.
    /// This couples even light walls without explicit impulse/velocity
    /// feedback, repeated collisions, or a fluid-state clone. The actual
    /// bounce-back populations are then evaluated at `u_mid`.
    ///
    /// Normal translation and rotation are constrained. This is tangential
    /// motion of fixed planar geometry, not moving-voxel topology. Body
    /// forces on the fluid remain separate external momentum input; periodic
    /// angular-momentum balances also require the periodic transport flux.
    ///
    /// # Errors
    /// Invalid geometry, nonpositive/nonfinite mass, nonfinite impulse or
    /// origin, collision failure, nonfinite exchange, or midpoint/final wall
    /// speed outside the existing low-Mach envelope. Every refusal leaves
    /// fluid populations and wall velocities unchanged.
    pub fn step_with_finite_mass_walls(
        &mut self,
        walls: [FiniteMassWall3; 2],
        moment_origin: [f64; 3],
    ) -> Result<[FiniteMassWallReceipt3; 2], BoundaryStepError3> {
        let faces = self.validate_finite_mass_walls(walls, moment_origin)?;
        self.try_collide()?;
        let proposal = self.finite_mass_wall_proposal(faces, walls, moment_origin, None)?;
        self.topology_locked = true;
        self.stream_with_boundaries(proposal.midpoint_boundaries);
        self.boundaries = proposal.final_boundaries;
        Ok(proposal.receipts)
    }

    /// Run a fixed number of deterministic steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// One step with the collide and pull-stream passes executed as tile
    /// kernels on `pool`; the open-face reconstruction (face cells only)
    /// stays serial. Every kernel tile writes a disjoint group of
    /// [`POOLED_GROUP_TILES`] lattice tiles from read-only inputs through the
    /// SAME per-cell functions as [`BoundaryGrid3::step`], so the published
    /// state is bit-identical to a serial step for any worker count.
    ///
    /// # Errors
    /// [`BoundaryStepError3`]; the pre-step state is intact in every case,
    /// and re-issuing the step is deterministic.
    pub fn step_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
    ) -> Result<(), BoundaryStepError3> {
        self.step_pooled_inner(pool, gate, None, None).map(|_| ())
    }

    /// Pooled counterpart of [`Self::step_with_wall_momentum`]. Loads reduce
    /// in the same selected-link order for every worker count. Measurement
    /// checkpoints occur at most 256 links apart and before publication.
    ///
    /// # Errors
    /// [`BoundaryStepError3`]. Invalid measurement, collision, cancellation,
    /// and non-finite exchange publish neither a new fluid state nor loads.
    pub fn step_pooled_with_wall_momentum<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        links: &[BoundaryLink3],
        moment_origin: [f64; 3],
    ) -> Result<WallMomentumExchange3, BoundaryStepError3> {
        self.validate_wall_measurement(links, moment_origin)?;
        self.step_pooled_inner(pool, gate, Some((links, moment_origin)), None)
    }

    /// Pooled counterpart of [`Self::step_with_finite_mass_walls`]. Both
    /// face reductions have a fixed order independent of worker count. A
    /// canceled or refused step publishes neither fluid nor wall motion.
    ///
    /// # Errors
    /// [`BoundaryStepError3`], including cancellation before publication.
    pub fn step_pooled_with_finite_mass_walls<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        walls: [FiniteMassWall3; 2],
        moment_origin: [f64; 3],
    ) -> Result<[FiniteMassWallReceipt3; 2], BoundaryStepError3> {
        let faces = self.validate_finite_mass_walls(walls, moment_origin)?;
        self.step_pooled_coupling_inner(pool, gate, None, None, Some((faces, walls, moment_origin)))
            .map(|(_, receipts)| receipts.expect("finite wall step returns its receipts"))
    }

    fn step_pooled_inner<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        measurement: Option<(&[BoundaryLink3], [f64; 3])>,
        models: Option<&[CollisionModel3]>,
    ) -> Result<WallMomentumExchange3, BoundaryStepError3> {
        self.step_pooled_coupling_inner(pool, gate, measurement, models, None)
            .map(|(exchange, _)| exchange)
    }

    fn step_pooled_coupling_inner<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        measurement: Option<(&[BoundaryLink3], [f64; 3])>,
        models: Option<&[CollisionModel3]>,
        finite_walls: Option<([Face3; 2], [FiniteMassWall3; 2], [f64; 3])>,
    ) -> Result<(WallMomentumExchange3, Option<[FiniteMassWallReceipt3; 2]>), BoundaryStepError3>
    {
        if finite_walls.is_none() {
            self.topology_locked = true;
        }
        {
            let kernel = CollidePass {
                f: &self.f,
                solid: &self.solid,
                post: chunk_groups(&mut self.post),
                model: self.collision_model,
                models,
                force: self.force,
                force_field: self.force_field.as_deref(),
            };
            match run_pass(pool, gate, &kernel)? {
                FirstCollision(None) => {}
                FirstCollision(Some((index, source))) => {
                    return Err(BoundaryStepError3::Collision {
                        cell: self.coords(index / TILE_CELLS, index % TILE_CELLS),
                        source,
                    });
                }
            }
        }
        let exchange = match measurement {
            Some((links, origin)) => self.wall_momentum(links, origin, Some(gate))?,
            None => WallMomentumExchange3::default(),
        };
        let proposal = match finite_walls {
            Some((faces, walls, origin)) => {
                Some(self.finite_mass_wall_proposal(faces, walls, origin, Some(gate))?)
            }
            None => None,
        };
        let mut next = self.next.take().unwrap_or_else(|| self.f.clone());
        let outcome = {
            let kernel = StreamPass {
                view: StreamView {
                    post: &self.post,
                    pull: &self.pull,
                    link_masks: &self.link_masks,
                    open_link_masks: &self.open_link_masks,
                    stationary_link_masks: &self.stationary_link_masks,
                    boundaries: proposal
                        .as_ref()
                        .map_or(self.boundaries, |proposal| proposal.midpoint_boundaries),
                    dims: [self.nx, self.ny, self.nz],
                },
                solid: &self.solid,
                next: chunk_groups(&mut next),
            };
            run_pass(pool, gate, &kernel)
        };
        if let Err(error) = outcome {
            self.next = Some(next);
            return Err(error);
        }
        self.publish_next(
            next,
            gate,
            measurement.is_some() || models.is_some() || proposal.is_some(),
        )?;
        let receipts = proposal.map(|proposal| {
            self.boundaries = proposal.final_boundaries;
            self.topology_locked = true;
            proposal.receipts
        });
        Ok((exchange, receipts))
    }

    /// Shared publication for the pooled owner and explicit serial SGS.
    /// All boundary reconstruction preconditions are checked before the swap.
    fn publish_next(
        &mut self,
        mut next: [Vec<Tile>; Q3],
        gate: &CancelGate,
        check_cancellation: bool,
    ) -> Result<(), BoundaryStepError3> {
        // The open-face reconstruction asserts a positive density at its
        // source cells; check them on the unpublished state so a lattice
        // that left its physical regime refuses instead of panicking.
        if self.open_cells.is_none() {
            self.open_cells = Some(self.compile_open_cells());
        }
        if let Some(cells) = &self.open_cells {
            for &(_, interior, _) in cells {
                let (tile, lane) = self.addr(interior[0], interior[1], interior[2]);
                let rho = (0..Q3).map(|q| next[q][tile].0[lane]).sum::<f64>();
                if !(rho.is_finite() && rho > 0.0) {
                    self.next = Some(next);
                    return Err(BoundaryStepError3::Unphysical { cell: interior });
                }
            }
        }
        if check_cancellation && gate.is_requested() {
            self.next = Some(next);
            return Err(BoundaryStepError3::Cancelled);
        }
        core::mem::swap(&mut self.f, &mut next);
        self.next = Some(next);
        self.apply_open_boundaries();
        Ok(())
    }

    /// Velocity component over one z section, row-major in `(y, x)`.
    #[must_use]
    pub fn velocity_section_z(&self, z: usize, component: usize) -> Vec<f64> {
        assert!(z < self.nz, "z section out of range");
        assert!(component < 3, "velocity component out of range");
        let mut section = Vec::with_capacity(self.nx * self.ny);
        for y in 0..self.ny {
            for x in 0..self.nx {
                section.push(if self.is_solid(x, y, z) {
                    0.0
                } else {
                    self.velocity(x, y, z)[component]
                });
            }
        }
        section
    }

    fn validate_finite_mass_walls(
        &self,
        walls: [FiniteMassWall3; 2],
        origin: [f64; 3],
    ) -> Result<[Face3; 2], BoundaryStepError3> {
        let refuse = |reason| BoundaryStepError3::FiniteMassWall { reason };
        if !origin.into_iter().all(f64::is_finite) {
            return Err(refuse("moment origin must be finite"));
        }
        for wall in walls {
            if !(wall.mass.is_finite() && wall.mass > 0.0) {
                return Err(refuse("wall mass must be positive and finite"));
            }
            if !wall.external_impulse.into_iter().all(f64::is_finite) {
                return Err(refuse("external wall impulse must be finite"));
            }
        }
        if self.solid.iter().any(|&word| word != 0) {
            return Err(refuse("finite planar walls do not admit voxel solids"));
        }
        let mut faces = None;
        for min in [Face3::XMin, Face3::YMin, Face3::ZMin] {
            let max = min.opposite();
            match (self.boundaries.face(min), self.boundaries.face(max)) {
                (FaceBoundary3::Periodic, FaceBoundary3::Periodic) => {}
                (FaceBoundary3::Wall { velocity: lo }, FaceBoundary3::Wall { velocity: hi })
                    if faces.is_none() =>
                {
                    if lo[min.axis()] != 0.0 || hi[min.axis()] != 0.0 {
                        return Err(refuse("finite planar walls require zero normal velocity"));
                    }
                    faces = Some([min, max]);
                }
                _ => {
                    return Err(refuse(
                        "finite walls require one opposite wall pair and periodic tangent axes",
                    ));
                }
            }
        }
        faces.ok_or_else(|| refuse("finite walls require one opposite wall pair"))
    }

    fn finite_mass_wall_proposal(
        &self,
        faces: [Face3; 2],
        walls: [FiniteMassWall3; 2],
        origin: [f64; 3],
        gate: Option<&CancelGate>,
    ) -> Result<FiniteMassWallProposal3, BoundaryStepError3> {
        let receipts = [
            self.finite_mass_wall_receipt(faces[0], walls[0], origin, gate)?,
            self.finite_mass_wall_receipt(faces[1], walls[1], origin, gate)?,
        ];
        let mut midpoint_boundaries = self.boundaries;
        let mut final_boundaries = self.boundaries;
        for receipt in receipts {
            midpoint_boundaries.faces[receipt.face.index()] = FaceBoundary3::Wall {
                velocity: receipt.velocity_midpoint,
            };
            final_boundaries.faces[receipt.face.index()] = FaceBoundary3::Wall {
                velocity: receipt.velocity_after,
            };
        }
        Ok(FiniteMassWallProposal3 {
            midpoint_boundaries,
            final_boundaries,
            receipts,
        })
    }

    fn finite_mass_wall_receipt(
        &self,
        face: Face3,
        wall: FiniteMassWall3,
        origin: [f64; 3],
        gate: Option<&CancelGate>,
    ) -> Result<FiniteMassWallReceipt3, BoundaryStepError3> {
        let dims = [self.nx, self.ny, self.nz];
        let normal = face.axis();
        let mut face_mass = 0.0;
        let mut stationary_impulse = [0.0; 3];
        for (index, cell) in planar_face_cells(dims, face).enumerate() {
            if index.is_multiple_of(32) && gate.is_some_and(CancelGate::is_requested) {
                return Err(BoundaryStepError3::Cancelled);
            }
            let (tile, lane) = self.addr(cell[0], cell[1], cell[2]);
            let rho = (0..Q3).map(|q| self.post[q][tile].0[lane]).sum::<f64>();
            if !(rho.is_finite() && rho > 0.0) {
                return Err(BoundaryStepError3::Unphysical { cell });
            }
            face_mass += rho;
            let mask = self.link_masks[tile].0[lane];
            for q in 1..Q3 {
                if mask & (1u32 << q) != 0 {
                    let outgoing = self.post[OPP3[q]][tile].0[lane];
                    let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
                    for axis in 0..3 {
                        stationary_impulse[axis] -= 2.0 * e[axis] * outgoing;
                    }
                }
            }
        }
        // Complete tangential link pairs cancel the moving-mass term, and
        // sum(2*w*e_t*e_t/c_s^2) = 1/3 on every planar boundary cell.
        let damping = face_mass / 3.0;
        let denominator = wall.mass + 0.5 * damping;
        if !(damping.is_finite() && denominator.is_finite() && denominator > 0.0) {
            return Err(BoundaryStepError3::FiniteMassWall {
                reason: "finite-wall inertia and fluid response must remain finite",
            });
        }
        let FaceBoundary3::Wall { velocity: before } = self.boundaries.face(face) else {
            unreachable!("finite-wall geometry was validated before collision")
        };
        let mut midpoint = [0.0; 3];
        let mut after = [0.0; 3];
        for axis in 0..3 {
            if axis != normal {
                // Solve for the increment: this avoids M*u overflow for a
                // heavy wall and never divides a load residual by a tiny M.
                let increment = (stationary_impulse[axis] + wall.external_impulse[axis]
                    - damping * before[axis])
                    / denominator;
                midpoint[axis] = before[axis] + 0.5 * increment;
                after[axis] = before[axis] + increment;
            }
        }
        for velocity in [midpoint, after] {
            if !velocity.into_iter().all(f64::is_finite)
                || velocity.into_iter().map(|u| u * u).sum::<f64>() >= 0.03
            {
                return Err(BoundaryStepError3::FiniteMassWall {
                    reason: "solved wall velocity exceeds the finite low-Mach envelope",
                });
            }
        }
        let mut hydrodynamic = WallMomentumExchange3::default();
        for (index, cell) in planar_face_cells(dims, face).enumerate() {
            if index.is_multiple_of(32) && gate.is_some_and(CancelGate::is_requested) {
                return Err(BoundaryStepError3::Cancelled);
            }
            let (tile, lane) = self.addr(cell[0], cell[1], cell[2]);
            let rho = (0..Q3).map(|q| self.post[q][tile].0[lane]).sum();
            let mask = self.link_masks[tile].0[lane];
            for (direction, opposite) in OPP3.into_iter().enumerate().skip(1) {
                if mask & (1u32 << direction) != 0 {
                    hydrodynamic.add_link(
                        BoundaryLink3 {
                            tile,
                            lane,
                            cell,
                            direction,
                        },
                        self.post[opposite][tile].0[lane],
                        rho,
                        midpoint,
                        origin,
                    );
                }
            }
        }
        let hydrodynamic = hydrodynamic.checked()?;
        let mut constraint_impulse = [0.0; 3];
        constraint_impulse[normal] =
            -(hydrodynamic.wall_impulse[normal] + wall.external_impulse[normal]);
        let arm = core::array::from_fn::<_, 3, _>(|axis| {
            let center = if axis == normal {
                if face.is_min() {
                    0.0
                } else {
                    dims[axis] as f64
                }
            } else {
                0.5 * dims[axis] as f64
            };
            center - origin[axis]
        });
        let constraint_angular_impulse = core::array::from_fn(|axis| {
            let b = (axis + 1) % 3;
            let c = (axis + 2) % 3;
            arm[b] * (hydrodynamic.wall_impulse[c] + constraint_impulse[c])
                - arm[c] * (hydrodynamic.wall_impulse[b] + constraint_impulse[b])
                - hydrodynamic.wall_angular_impulse[axis]
        });
        let kinetic_energy_change = (0..3)
            .map(|axis| {
                0.5 * wall.mass * (after[axis] - before[axis]) * (after[axis] + before[axis])
            })
            .sum();
        let external_work = (0..3)
            .map(|axis| wall.external_impulse[axis] * midpoint[axis])
            .sum();
        if !constraint_impulse
            .into_iter()
            .chain(constraint_angular_impulse)
            .chain([kinetic_energy_change, external_work])
            .all(f64::is_finite)
        {
            return Err(BoundaryStepError3::FiniteMassWall {
                reason: "wall constraint reaction and work must remain finite",
            });
        }
        if gate.is_some_and(CancelGate::is_requested) {
            return Err(BoundaryStepError3::Cancelled);
        }
        Ok(FiniteMassWallReceipt3 {
            face,
            mass: wall.mass,
            velocity_before: before,
            velocity_midpoint: midpoint,
            velocity_after: after,
            hydrodynamic,
            external_impulse: wall.external_impulse,
            constraint_impulse,
            constraint_angular_impulse,
            kinetic_energy_change,
            external_work,
        })
    }

    fn validate_wall_measurement(
        &self,
        links: &[BoundaryLink3],
        origin: [f64; 3],
    ) -> Result<(), BoundaryStepError3> {
        if !origin.into_iter().all(f64::is_finite) {
            return Err(BoundaryStepError3::WallMomentum {
                reason: "moment origin must be finite",
            });
        }
        let mut previous = None;
        for link in links {
            let key = (link.tile, link.lane, link.direction);
            if previous.is_some_and(|previous| previous >= key) {
                return Err(BoundaryStepError3::WallMomentum {
                    reason: "wall links must be unique and in canonical order",
                });
            }
            if link.tile >= self.link_masks.len()
                || link.lane >= TILE_CELLS
                || link.direction == 0
                || link.direction >= Q3
                || self.coords(link.tile, link.lane) != link.cell
                || self.link_masks[link.tile].0[link.lane] & (1u32 << link.direction) == 0
            {
                return Err(BoundaryStepError3::WallMomentum {
                    reason: "selection contains a link outside the grid's wall topology",
                });
            }
            previous = Some(key);
        }
        Ok(())
    }

    fn wall_momentum(
        &self,
        links: &[BoundaryLink3],
        origin: [f64; 3],
        gate: Option<&CancelGate>,
    ) -> Result<WallMomentumExchange3, BoundaryStepError3> {
        let mut exchange = WallMomentumExchange3::default();
        let mut previous_cell = None;
        let mut rho = 0.0;
        let mut exterior_velocity = [0.0; 3];
        for (index, link) in links.iter().enumerate() {
            if index.is_multiple_of(256) && gate.is_some_and(CancelGate::is_requested) {
                return Err(BoundaryStepError3::Cancelled);
            }
            let cell = (link.tile, link.lane);
            if previous_cell != Some(cell) {
                rho = (0..Q3).map(|q| self.post[q][link.tile].0[link.lane]).sum();
                exterior_velocity = effective_wall_velocity(
                    self.boundaries,
                    [self.nx, self.ny, self.nz],
                    link.cell,
                );
                previous_cell = Some(cell);
            }
            let velocity = if self.stationary_link_masks[link.tile].0[link.lane]
                & (1u32 << link.direction)
                != 0
            {
                [0.0; 3]
            } else {
                exterior_velocity
            };
            exchange.add_link(
                *link,
                self.post[OPP3[link.direction]][link.tile].0[link.lane],
                rho,
                velocity,
                origin,
            );
        }
        if gate.is_some_and(CancelGate::is_requested) {
            return Err(BoundaryStepError3::Cancelled);
        }
        exchange.checked()
    }

    fn addr(&self, x: usize, y: usize, z: usize) -> (usize, usize) {
        let (ntx, nty) = (self.nx / TILE, self.ny / TILE);
        let tile = (z / TILE * nty + y / TILE) * ntx + x / TILE;
        let lane = (z % TILE * TILE + y % TILE) * TILE + x % TILE;
        (tile, lane)
    }

    fn addr_checked(&self, x: usize, y: usize, z: usize) -> (usize, usize) {
        assert!(
            x < self.nx && y < self.ny && z < self.nz,
            "cell coordinate out of range"
        );
        self.addr(x, y, z)
    }

    fn coords(&self, tile: usize, lane: usize) -> [usize; 3] {
        tile_lane_coords([self.nx, self.ny, self.nz], tile, lane)
    }

    fn raw_moments(&self, x: usize, y: usize, z: usize) -> (f64, [f64; 3]) {
        let (tile, lane) = self.addr(x, y, z);
        let mut rho = 0.0;
        let mut momentum = [0.0; 3];
        for (q, velocity) in E3.iter().enumerate() {
            let fq = self.f[q][tile].0[lane];
            rho += fq;
            momentum[0] += f64::from(velocity.0) * fq;
            momentum[1] += f64::from(velocity.1) * fq;
            momentum[2] += f64::from(velocity.2) * fq;
        }
        assert!(
            rho.is_finite() && rho > 0.0,
            "fluid density must remain positive and finite"
        );
        (rho, momentum)
    }

    fn write_equilibrium(&mut self, x: usize, y: usize, z: usize, rho: f64, velocity: [f64; 3]) {
        let equilibrium = equilibrium3(rho, velocity);
        let (tile, lane) = self.addr(x, y, z);
        for (field, value) in self.f.iter_mut().zip(equilibrium) {
            field[tile].0[lane] = value;
        }
    }

    fn collide(&mut self) {
        self.try_collide()
            .expect("BoundaryGrid3 constructor and prior state admit selected collision");
    }

    fn try_collide(&mut self) -> Result<(), BoundaryStepError3> {
        for tile in 0..self.f[0].len() {
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let populations = core::array::from_fn(|direction| self.f[direction][tile].0[lane]);
                let force = cell_force(self.force, self.force_field.as_deref(), tile, lane);
                let post =
                    collide_cell3(populations, self.collision_model, force).map_err(|source| {
                        BoundaryStepError3::Collision {
                            cell: self.coords(tile, lane),
                            source,
                        }
                    })?;
                for (field, value) in self.post.iter_mut().zip(post) {
                    field[tile].0[lane] = value;
                }
            }
        }
        Ok(())
    }

    fn stream(&mut self) {
        self.stream_with_boundaries(self.boundaries);
    }

    fn stream_with_boundaries(&mut self, boundaries: BoundarySpec3) {
        // Canonical tile/lane traversal. Every destination population is a
        // pure read of the immutable post-collision state, so the visiting
        // order cannot move result bits.
        let view = StreamView {
            post: &self.post,
            pull: &self.pull,
            link_masks: &self.link_masks,
            open_link_masks: &self.open_link_masks,
            stationary_link_masks: &self.stationary_link_masks,
            boundaries,
            dims: [self.nx, self.ny, self.nz],
        };
        for tile in 0..self.f[0].len() {
            let solid = self.solid[tile];
            for lane in 0..TILE_CELLS {
                if solid & (1u64 << lane) != 0 {
                    continue;
                }
                let values = view.pull_cell(tile, lane);
                for (field, value) in self.f.iter_mut().zip(values) {
                    field[tile].0[lane] = value;
                }
            }
        }
    }

    fn compile_pull_table(&self, solid: &[u64]) -> Vec<u32> {
        let cells = self.f[0].len() * TILE_CELLS;
        assert!(
            u32::try_from(cells).is_ok_and(|n| n < NO_FLUID_SOURCE),
            "BoundaryGrid3 pull table supports fewer than 2^32 - 1 cells"
        );
        let mut pull = vec![NO_FLUID_SOURCE; cells * Q3];
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    if self.is_solid_in(solid, x, y, z) {
                        continue;
                    }
                    let (tile, lane) = self.addr(x, y, z);
                    let base = (tile * TILE_CELLS + lane) * Q3;
                    for q in 0..Q3 {
                        if let PullSource::Fluid { tile, lane } =
                            self.classify_pull(solid, x, y, z, q)
                        {
                            pull[base + q] = (tile * TILE_CELLS + lane) as u32;
                        }
                    }
                }
            }
        }
        pull
    }

    fn compile_open_cells(&self) -> Vec<([usize; 3], [usize; 3], Face3)> {
        let mut cells = Vec::new();
        for face in Face3::ALL {
            if !self.boundaries.face(face).is_open() {
                continue;
            }
            for z in 0..self.nz {
                for y in 0..self.ny {
                    for x in 0..self.nx {
                        if !self.on_face(x, y, z, face)
                            || self.is_solid(x, y, z)
                            || self.open_link_mask(x, y, z) == 0
                        {
                            continue;
                        }
                        let interior = Self::interior_neighbor(x, y, z, face);
                        assert!(
                            !self.is_solid(interior[0], interior[1], interior[2]),
                            "open boundary {face:?} requires a fluid first-interior neighbor"
                        );
                        cells.push(([x, y, z], interior, face));
                    }
                }
            }
        }
        cells
    }

    fn apply_open_boundaries(&mut self) {
        if self.open_cells.is_none() {
            self.open_cells = Some(self.compile_open_cells());
        }
        let cells = self.open_cells.take().unwrap_or_default();
        for &([x, y, z], interior, face) in &cells {
            let (neighbor_rho, neighbor_velocity, stress) =
                self.regularized_source(interior[0], interior[1], interior[2]);
            let (rho, velocity) = match self.boundaries.face(face) {
                FaceBoundary3::Velocity { velocity } => (
                    neighbor_rho,
                    self.profile_velocity(face, [x, y, z]).unwrap_or(velocity),
                ),
                FaceBoundary3::Pressure { density } => (density, neighbor_velocity),
                FaceBoundary3::Periodic | FaceBoundary3::Wall { .. } => unreachable!(),
            };
            let reconstructed = regularized_populations(rho, velocity, stress);
            let (tile, lane) = self.addr(x, y, z);
            let wall_mask = self.link_masks[tile].0[lane];
            let open_mask = self.open_link_masks[tile].0[lane];
            if wall_mask == 0 {
                for (field, value) in self.f.iter_mut().zip(reconstructed) {
                    field[tile].0[lane] = value;
                }
            } else {
                // Mixed wall/open cells retain their bounced wall
                // populations and streamed tangential populations.
                // Only directions whose source crosses the open
                // face alone receive the regularized closure.
                for (q, value) in reconstructed.into_iter().enumerate().skip(1) {
                    if open_mask & (1u32 << q) != 0 {
                        self.f[q][tile].0[lane] = value;
                    }
                }
            }
        }
        self.open_cells = Some(cells);
    }

    fn regularized_source(&self, x: usize, y: usize, z: usize) -> (f64, [f64; 3], [[f64; 3]; 3]) {
        let (rho, momentum) = self.raw_moments(x, y, z);
        let velocity = core::array::from_fn(|axis| (momentum[axis] + 0.5 * self.force[axis]) / rho);
        let equilibrium = equilibrium3(rho, velocity);
        let (tile, lane) = self.addr(x, y, z);
        let mut stress = [[0.0; 3]; 3];
        for q in 0..Q3 {
            let nonequilibrium = self.f[q][tile].0[lane] - equilibrium[q];
            let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
            for row in 0..3 {
                for column in 0..3 {
                    stress[row][column] += e[row] * e[column] * nonequilibrium;
                }
            }
        }
        (rho, velocity, stress)
    }

    fn compile_link_masks(
        &self,
        solid: &[u64],
    ) -> (Vec<LinkMaskTile3>, Vec<LinkMaskTile3>, Vec<LinkMaskTile3>) {
        let mut wall_masks = vec![LinkMaskTile3::empty(); self.link_masks.len()];
        let mut open_masks = vec![LinkMaskTile3::empty(); self.link_masks.len()];
        let mut stationary_masks = vec![LinkMaskTile3::empty(); self.link_masks.len()];
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    if self.is_solid_in(solid, x, y, z) {
                        continue;
                    }
                    let (tile, lane) = self.addr(x, y, z);
                    let mut wall_mask = 0u32;
                    let mut open_mask = 0u32;
                    let mut stationary_mask = 0u32;
                    for q in 1..Q3 {
                        match self.classify_pull(solid, x, y, z, q) {
                            PullSource::Wall { voxel } => {
                                wall_mask |= 1u32 << q;
                                if voxel {
                                    stationary_mask |= 1u32 << q;
                                }
                            }
                            PullSource::Open => open_mask |= 1u32 << q,
                            PullSource::Fluid { .. } => {}
                        }
                    }
                    // A mixed wall/open seam keeps disjoint per-link owners:
                    // a direction crossing both faces is wall-owned, while a
                    // pure open-normal direction is reconstructed. Wall-owned
                    // links on the mixed cell are stationary so a moving wall
                    // cannot inject an unpaired corner correction.
                    if wall_mask != 0 && open_mask != 0 {
                        stationary_mask |= wall_mask;
                    }
                    wall_masks[tile].0[lane] = wall_mask;
                    open_masks[tile].0[lane] = open_mask;
                    stationary_masks[tile].0[lane] = stationary_mask;
                }
            }
        }
        (wall_masks, open_masks, stationary_masks)
    }

    fn validate_open_neighbors(&self, solid: &[u64]) {
        for face in Face3::ALL {
            if !self.boundaries.face(face).is_open() {
                continue;
            }
            for z in 0..self.nz {
                for y in 0..self.ny {
                    for x in 0..self.nx {
                        if self.on_face(x, y, z, face) && !self.is_solid_in(solid, x, y, z) {
                            let neighbor = Self::interior_neighbor(x, y, z, face);
                            assert!(
                                !self.is_solid_in(solid, neighbor[0], neighbor[1], neighbor[2]),
                                "open boundary {face:?} requires a fluid first-interior neighbor"
                            );
                        }
                    }
                }
            }
        }
    }

    fn is_solid_in(&self, solid: &[u64], x: usize, y: usize, z: usize) -> bool {
        let (tile, lane) = self.addr(x, y, z);
        solid[tile] & (1u64 << lane) != 0
    }

    fn classify_pull(&self, solid: &[u64], x: usize, y: usize, z: usize, q: usize) -> PullSource {
        let mut coordinate = [x, y, z];
        let velocity = E3[q];
        let deltas = [-velocity.0, -velocity.1, -velocity.2];
        let dimensions = [self.nx, self.ny, self.nz];
        let mut crossed_wall = false;
        let mut open = false;
        for axis in 0..3 {
            let crossed = match deltas[axis] {
                -1 if coordinate[axis] == 0 => Some(Face3::ALL[2 * axis]),
                1 if coordinate[axis] + 1 == dimensions[axis] => Some(Face3::ALL[2 * axis + 1]),
                -1 => {
                    coordinate[axis] -= 1;
                    None
                }
                1 => {
                    coordinate[axis] += 1;
                    None
                }
                _ => None,
            };
            let Some(face) = crossed else {
                continue;
            };
            match self.boundaries.face(face) {
                FaceBoundary3::Periodic => {
                    coordinate[axis] = if face.is_min() {
                        dimensions[axis] - 1
                    } else {
                        0
                    };
                }
                FaceBoundary3::Wall { .. } => {
                    crossed_wall = true;
                }
                FaceBoundary3::Velocity { .. } | FaceBoundary3::Pressure { .. } => {
                    open = true;
                }
            }
        }
        if crossed_wall {
            return PullSource::Wall { voxel: false };
        }
        if open {
            return PullSource::Open;
        }
        let (tile, lane) = self.addr(coordinate[0], coordinate[1], coordinate[2]);
        if solid[tile] & (1u64 << lane) != 0 {
            PullSource::Wall { voxel: true }
        } else {
            PullSource::Fluid { tile, lane }
        }
    }

    fn on_face(&self, x: usize, y: usize, z: usize, face: Face3) -> bool {
        match face {
            Face3::XMin => x == 0,
            Face3::XMax => x + 1 == self.nx,
            Face3::YMin => y == 0,
            Face3::YMax => y + 1 == self.ny,
            Face3::ZMin => z == 0,
            Face3::ZMax => z + 1 == self.nz,
        }
    }

    fn interior_neighbor(x: usize, y: usize, z: usize, face: Face3) -> [usize; 3] {
        match face {
            Face3::XMin => [x + 1, y, z],
            Face3::XMax => [x - 1, y, z],
            Face3::YMin => [x, y + 1, z],
            Face3::YMax => [x, y - 1, z],
            Face3::ZMin => [x, y, z + 1],
            Face3::ZMax => [x, y, z - 1],
        }
    }
}

impl core::fmt::Debug for BoundaryGrid3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("BoundaryGrid3")
            .field("nx", &self.nx)
            .field("ny", &self.ny)
            .field("nz", &self.nz)
            .field("collision_model", &self.collision_model)
            .field("force", &self.force)
            .field("boundaries", &self.boundaries)
            .field(
                "solid_cells",
                &self.solid.iter().map(|bits| bits.count_ones()).sum::<u32>(),
            )
            .finish_non_exhaustive()
    }
}

/// The uniform force, plus the cell's field entry when a field is installed
/// (without a field the uniform value is returned untouched, bit for bit).
#[inline]
fn cell_force(uniform: [f64; 3], field: Option<&[[f64; 3]]>, tile: usize, lane: usize) -> [f64; 3] {
    match field {
        None => uniform,
        Some(field) => {
            let local = field[tile * TILE_CELLS + lane];
            [
                uniform[0] + local[0],
                uniform[1] + local[1],
                uniform[2] + local[2],
            ]
        }
    }
}

fn tile_lane_coords(dims: [usize; 3], tile: usize, lane: usize) -> [usize; 3] {
    let (ntx, nty) = (dims[0] / TILE, dims[1] / TILE);
    let tx = tile % ntx;
    let rem = tile / ntx;
    let ty = rem % nty;
    let tz = rem / nty;
    let lx = lane % TILE;
    let ly = (lane / TILE) % TILE;
    let lz = lane / (TILE * TILE);
    [tx * TILE + lx, ty * TILE + ly, tz * TILE + lz]
}

/// Face-only traversal in global coordinate order; no grid-volume scan or
/// temporary link list is needed by the finite-wall response calculation.
fn planar_face_cells(dims: [usize; 3], face: Face3) -> impl Iterator<Item = [usize; 3]> {
    let normal = face.axis();
    let tangents = match normal {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    };
    (0..dims[tangents[0]] * dims[tangents[1]]).map(move |index| {
        let mut cell = [0; 3];
        cell[normal] = if face.is_min() { 0 } else { dims[normal] - 1 };
        cell[tangents[0]] = index % dims[tangents[0]];
        cell[tangents[1]] = index / dims[tangents[0]];
        cell
    })
}

/// One velocity owns every exterior wall link of a boundary cell. At a
/// seam where incident wall faces disagree (for example, moving lid meets
/// stationary side wall), the entire cell is stationary. This keeps the
/// moving-wall corrections pairwise balanced instead of injecting mass
/// through an arbitrarily owned diagonal link.
fn effective_wall_velocity(
    boundaries: BoundarySpec3,
    dims: [usize; 3],
    cell: [usize; 3],
) -> [f64; 3] {
    let mut velocity: Option<[f64; 3]> = None;
    for face in Face3::ALL {
        let axis = face.axis();
        let on_face = if face.is_min() {
            cell[axis] == 0
        } else {
            cell[axis] + 1 == dims[axis]
        };
        if !on_face {
            continue;
        }
        let FaceBoundary3::Wall {
            velocity: face_velocity,
        } = boundaries.face(face)
        else {
            continue;
        };
        velocity = Some(match velocity {
            Some(existing) if !same_velocity(existing, face_velocity) => [0.0; 3],
            Some(existing) => existing,
            None => face_velocity,
        });
    }
    velocity.unwrap_or([0.0; 3])
}

/// Read-only inputs of the pull-stream, shared by the serial and pooled
/// passes so both evaluate the identical per-cell expression.
struct StreamView<'a> {
    post: &'a [Vec<Tile>; Q3],
    pull: &'a [u32],
    link_masks: &'a [LinkMaskTile3],
    open_link_masks: &'a [LinkMaskTile3],
    stationary_link_masks: &'a [LinkMaskTile3],
    boundaries: BoundarySpec3,
    dims: [usize; 3],
}

// One operation tree for the actual incoming population and its measured
// momentum exchange. Keep the existing streaming arithmetic bit-for-bit.
fn wall_population(outgoing: f64, rho: f64, q: usize, velocity: [f64; 3]) -> f64 {
    let e = E3[q];
    let eu_wall = f64::from(e.0).mul_add(
        velocity[0],
        f64::from(e.1).mul_add(velocity[1], f64::from(e.2) * velocity[2]),
    );
    outgoing + 2.0 * W3[q] * rho * eu_wall / CS2
}

impl StreamView<'_> {
    /// Post-stream populations of one fluid cell before open-face
    /// reconstruction.
    fn pull_cell(&self, tile: usize, lane: usize) -> [f64; Q3] {
        let post = self.post;
        let base = (tile * TILE_CELLS + lane) * Q3;
        let wall_mask = self.link_masks[tile].0[lane];
        let open_mask = self.open_link_masks[tile].0[lane];
        if wall_mask | open_mask == 0 {
            return core::array::from_fn(|q| {
                let source = self.pull[base + q] as usize;
                post[q][source / TILE_CELLS].0[source % TILE_CELLS]
            });
        }
        let rho_post = (0..Q3).map(|q| post[q][tile].0[lane]).sum::<f64>();
        let stationary_mask = self.stationary_link_masks[tile].0[lane];
        core::array::from_fn(|q| {
            let bit = 1u32 << q;
            if wall_mask & bit != 0 {
                let velocity = if stationary_mask & bit != 0 {
                    [0.0; 3]
                } else {
                    effective_wall_velocity(
                        self.boundaries,
                        self.dims,
                        tile_lane_coords(self.dims, tile, lane),
                    )
                };
                wall_population(post[OPP3[q]][tile].0[lane], rho_post, q, velocity)
            } else if open_mask & bit != 0 {
                // Every population on an open face-interior cell is replaced
                // by the regularized pass that follows streaming.
                post[OPP3[q]][tile].0[lane]
            } else {
                let source = self.pull[base + q];
                debug_assert_ne!(
                    source, NO_FLUID_SOURCE,
                    "compiled boundary masks must classify every non-fluid link"
                );
                let source = source as usize;
                post[q][source / TILE_CELLS].0[source % TILE_CELLS]
            }
        })
    }
}

/// Split the 19 SoA fields into aligned groups of [`POOLED_GROUP_TILES`]
/// tiles, one exclusive (uncontended) lock per kernel tile.
fn chunk_groups(fields: &mut [Vec<Tile>; Q3]) -> Vec<Mutex<[&mut [Tile]; Q3]>> {
    let groups = fields[0].len().div_ceil(POOLED_GROUP_TILES);
    let mut iters: Vec<_> = fields
        .iter_mut()
        .map(|field| field.chunks_mut(POOLED_GROUP_TILES))
        .collect();
    (0..groups)
        .map(|_| {
            Mutex::new(core::array::from_fn(|q| {
                iters[q]
                    .next()
                    .expect("every SoA field has the same tile count")
            }))
        })
        .collect()
}

/// First collision refusal in canonical (tile, lane) order.
struct FirstCollision(Option<(usize, CollisionError3)>);

impl Reduce for FirstCollision {
    fn identity() -> Self {
        FirstCollision(None)
    }

    fn merge(self, other: Self) -> Self {
        if self.0.is_some() { self } else { other }
    }
}

struct CollidePass<'a> {
    f: &'a [Vec<Tile>; Q3],
    solid: &'a [u64],
    post: Vec<Mutex<[&'a mut [Tile]; Q3]>>,
    model: CollisionModel3,
    models: Option<&'a [CollisionModel3]>,
    force: [f64; 3],
    force_field: Option<&'a [[f64; 3]]>,
}

impl TileKernel for CollidePass<'_> {
    type Out = FirstCollision;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("fs-lbm.boundary3.collide", self.post.len() as u64)
    }

    fn run(&self, group: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, FirstCollision> {
        if let Err(cancelled) = cx.checkpoint() {
            return ControlFlow::Break(cancelled);
        }
        let group = group as usize;
        let mut chunk = self.post[group]
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        for local in 0..chunk[0].len() {
            let tile = group * POOLED_GROUP_TILES + local;
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let populations = core::array::from_fn(|q| self.f[q][tile].0[lane]);
                let force = cell_force(self.force, self.force_field, tile, lane);
                let model = self
                    .models
                    .map_or(self.model, |field| field[tile * TILE_CELLS + lane]);
                match collide_cell3(populations, model, force) {
                    Ok(post) => {
                        for (field, value) in chunk.iter_mut().zip(post) {
                            field[local].0[lane] = value;
                        }
                    }
                    Err(source) => {
                        return ControlFlow::Continue(FirstCollision(Some((
                            tile * TILE_CELLS + lane,
                            source,
                        ))));
                    }
                }
            }
        }
        ControlFlow::Continue(FirstCollision(None))
    }
}

struct StreamPass<'a> {
    view: StreamView<'a>,
    solid: &'a [u64],
    next: Vec<Mutex<[&'a mut [Tile]; Q3]>>,
}

impl TileKernel for StreamPass<'_> {
    type Out = ();

    fn tiles(&self) -> TilePlan {
        TilePlan::new("fs-lbm.boundary3.stream", self.next.len() as u64)
    }

    fn run(&self, group: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, ()> {
        if let Err(cancelled) = cx.checkpoint() {
            return ControlFlow::Break(cancelled);
        }
        let group = group as usize;
        let mut chunk = self.next[group]
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        for local in 0..chunk[0].len() {
            let tile = group * POOLED_GROUP_TILES + local;
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let values = self.view.pull_cell(tile, lane);
                for (field, value) in chunk.iter_mut().zip(values) {
                    field[local].0[lane] = value;
                }
            }
        }
        ControlFlow::Continue(())
    }
}

fn run_pass<P: KernelRunner, K: TileKernel>(
    pool: &P,
    gate: &CancelGate,
    kernel: &K,
) -> Result<K::Out, BoundaryStepError3> {
    let (outcome, _report) = pool.run_with_gate(kernel, gate);
    outcome.map_err(|error| {
        if gate.is_requested() {
            BoundaryStepError3::Cancelled
        } else {
            BoundaryStepError3::Pool(format!("{error:?}"))
        }
    })
}

fn same_velocity(left: [f64; 3], right: [f64; 3]) -> bool {
    left.into_iter().zip(right).all(|(left, right)| {
        left.to_bits() == right.to_bits()
            || (matches!(left.classify(), core::num::FpCategory::Zero)
                && matches!(right.classify(), core::num::FpCategory::Zero))
    })
}

fn regularized_populations(rho: f64, velocity: [f64; 3], stress: [[f64; 3]; 3]) -> [f64; Q3] {
    let mut populations = equilibrium3(rho, velocity);
    let coefficient = 1.0 / (2.0 * CS2 * CS2);
    for q in 0..Q3 {
        let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
        let mut contraction = 0.0;
        for row in 0..3 {
            for column in 0..3 {
                let isotropic = if row == column { CS2 } else { 0.0 };
                contraction += (e[row] * e[column] - isotropic) * stress[row][column];
            }
        }
        populations[q] += coefficient * W3[q] * contraction;
    }
    populations
}
