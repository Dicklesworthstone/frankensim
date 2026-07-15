//! Link-wise boundary conditions for the D3Q19 tile layout.
//!
//! This module deliberately leaves [`super::Duct`] unchanged: its frozen
//! bit surface remains the body-force/periodic fixture, while
//! [`BoundaryGrid3`] owns the more general boundary semantics.  Solid links
//! use halfway bounce-back, moving planar walls add the standard momentum
//! correction, and open planar faces use a regularized non-equilibrium stress
//! reconstruction.  The latter is face-generic and preserves the prescribed
//! density/velocity moments without maintaining six hand-specialized D3Q19
//! formula tables.

use super::{
    CollisionModel3, E3, OPP3, Q3, TILE, TILE_CELLS, Tile, W3, collide_cell3, equilibrium3,
};
use crate::CS2;

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
    tau: f64,
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
}

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
        assert!(
            nx > 0
                && ny > 0
                && nz > 0
                && nx.is_multiple_of(TILE)
                && ny.is_multiple_of(TILE)
                && nz.is_multiple_of(TILE),
            "grid dimensions must be positive multiples of {TILE} (got {nx}x{ny}x{nz})"
        );
        assert!(
            tau.is_finite() && tau > 0.5,
            "relaxation time tau must be finite and greater than 0.5"
        );
        assert!(
            force.iter().all(|component| component.is_finite()),
            "body force must be finite"
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
            tau,
            force,
            boundaries,
            f,
            post,
            solid: vec![0; tiles],
            link_masks: vec![LinkMaskTile3::empty(); tiles],
            open_link_masks: vec![LinkMaskTile3::empty(); tiles],
            stationary_link_masks: vec![LinkMaskTile3::empty(); tiles],
            topology_locked: false,
        };
        let (wall_masks, open_masks, stationary_masks) = grid.compile_link_masks(&grid.solid);
        grid.link_masks = wall_masks;
        grid.open_link_masks = open_masks;
        grid.stationary_link_masks = stationary_masks;
        grid
    }

    /// Grid dimensions in cells.
    #[must_use]
    pub const fn dimensions(&self) -> [usize; 3] {
        [self.nx, self.ny, self.nz]
    }

    /// Kinematic viscosity `nu = (tau - 1/2)/3`.
    #[must_use]
    pub fn viscosity(&self) -> f64 {
        (self.tau - 0.5) / 3.0
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
        self.solid = proposed_solid;
        self.link_masks = wall_masks;
        self.open_link_masks = open_masks;
        self.stationary_link_masks = stationary_masks;
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
        core::array::from_fn(|axis| (momentum[axis] + 0.5 * self.force[axis]) / rho)
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

    /// One BGK/Guo collide, link-wise pull-stream, and open-face
    /// reconstruction step.
    pub fn step(&mut self) {
        self.topology_locked = true;
        self.collide();
        self.stream();
        self.apply_open_boundaries();
    }

    /// Run a fixed number of deterministic steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
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
        let (ntx, nty) = (self.nx / TILE, self.ny / TILE);
        let tx = tile % ntx;
        let rem = tile / ntx;
        let ty = rem % nty;
        let tz = rem / nty;
        let lx = lane % TILE;
        let ly = (lane / TILE) % TILE;
        let lz = lane / (TILE * TILE);
        [tx * TILE + lx, ty * TILE + ly, tz * TILE + lz]
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
        for tile in 0..self.f[0].len() {
            for lane in 0..TILE_CELLS {
                if self.solid[tile] & (1u64 << lane) != 0 {
                    continue;
                }
                let populations = core::array::from_fn(|direction| self.f[direction][tile].0[lane]);
                let post = collide_cell3(
                    populations,
                    CollisionModel3::Bgk { tau: self.tau },
                    self.force,
                )
                .expect("BoundaryGrid3 constructor and prior state admit BGK/Guo collision");
                for (field, value) in self.post.iter_mut().zip(post) {
                    field[tile].0[lane] = value;
                }
            }
        }
    }

    fn stream(&mut self) {
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    if self.is_solid(x, y, z) {
                        continue;
                    }
                    let (destination_tile, destination_lane) = self.addr(x, y, z);
                    let rho_post = (0..Q3)
                        .map(|q| self.post[q][destination_tile].0[destination_lane])
                        .sum::<f64>();
                    let wall_mask = self.link_masks[destination_tile].0[destination_lane];
                    let open_mask = self.open_link_masks[destination_tile].0[destination_lane];
                    let stationary_mask =
                        self.stationary_link_masks[destination_tile].0[destination_lane];
                    for q in 0..Q3 {
                        let bit = 1u32 << q;
                        let value = if wall_mask & bit != 0 {
                            let velocity = if stationary_mask & bit != 0 {
                                [0.0; 3]
                            } else {
                                self.effective_wall_velocity(x, y, z)
                            };
                            let e = E3[q];
                            let eu_wall = f64::from(e.0).mul_add(
                                velocity[0],
                                f64::from(e.1).mul_add(velocity[1], f64::from(e.2) * velocity[2]),
                            );
                            self.post[OPP3[q]][destination_tile].0[destination_lane]
                                + 2.0 * W3[q] * rho_post * eu_wall / CS2
                        } else if open_mask & bit != 0 {
                            // Every population on an open face-interior cell
                            // is replaced by the regularized pass below.
                            self.post[OPP3[q]][destination_tile].0[destination_lane]
                        } else {
                            let PullSource::Fluid { tile, lane } =
                                self.classify_pull(&self.solid, x, y, z, q)
                            else {
                                unreachable!(
                                    "compiled boundary masks must classify every non-fluid link"
                                )
                            };
                            self.post[q][tile].0[lane]
                        };
                        self.f[q][destination_tile].0[destination_lane] = value;
                    }
                }
            }
        }
    }

    fn apply_open_boundaries(&mut self) {
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
                        let (neighbor_rho, neighbor_velocity, stress) =
                            self.regularized_source(interior[0], interior[1], interior[2]);
                        let (rho, velocity) = match self.boundaries.face(face) {
                            FaceBoundary3::Velocity { velocity } => (neighbor_rho, velocity),
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
                }
            }
        }
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

    /// One velocity owns every exterior wall link of a boundary cell. At a
    /// seam where incident wall faces disagree (for example, moving lid meets
    /// stationary side wall), the entire cell is stationary. This keeps the
    /// moving-wall corrections pairwise balanced instead of injecting mass
    /// through an arbitrarily owned diagonal link.
    fn effective_wall_velocity(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        let mut velocity: Option<[f64; 3]> = None;
        for face in Face3::ALL {
            if !self.on_face(x, y, z, face) {
                continue;
            }
            let FaceBoundary3::Wall {
                velocity: face_velocity,
            } = self.boundaries.face(face)
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
            .field("tau", &self.tau)
            .field("force", &self.force)
            .field("boundaries", &self.boundaries)
            .field(
                "solid_cells",
                &self.solid.iter().map(|bits| bits.count_ones()).sum::<u32>(),
            )
            .finish_non_exhaustive()
    }
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
