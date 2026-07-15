//! The D3Q19 core (bead 84hv): 3-D BGK stream-and-collide with Guo body
//! forcing and halfway bounce-back walls, over 128-byte-aligned SoA
//! distributions in tile-major layout (plan §5) with the PULL scheme —
//! the 3-D sibling of the D2Q9 module in `lib.rs`, which it deliberately
//! mirrors (and leaves untouched).
//!
//! DETERMINISTIC TRAVERSAL (pinned): both the collide and stream loops
//! visit tiles in ascending tile index (x-fastest, then y, then z over
//! the tile grid) and cells within a tile in ascending local index
//! (x-fastest, then y, then z). Single-threaded in this bead; when WS1-D
//! parallelizes the sweep, per-cell writes stay slot-exclusive (collide
//! is pointwise; pull-streaming writes only the destination cell), so
//! this order is documentation of the CANONICAL schedule, not a hidden
//! result dependence — results must stay bit-identical across worker
//! counts by construction.
//!
//! Layout: the domain is split into 4×4×4 tiles (64 cells, 512 B per
//! field per tile, 128-byte aligned). Each of the 19 distribution
//! fields is its own `Vec<Tile>` (structure of arrays); cell `(x,y,z)`
//! lives in tile `(x/4, y/4, z/4)` at local `(x%4, y%4, z%4)`.
//! Dimensions must be multiples of 4 (asserted): the fixture scales in
//! WS1 are, and padding is a later concern, not silent slop.

mod boundary;

pub use boundary::{
    BoundaryGrid3, BoundaryLink3, BoundarySpec3, D3Q19_BOUNDARY_BIT_SEMANTICS_VERSION, Face3,
    FaceBoundary3, LinkMaskTile3,
};

/// Bit-semantics version of the D3Q19 surface (golden-couplings.json):
/// covers the velocity/weight/opposite tables and ordering, the
/// equilibrium form, the Guo forcing form, the pull-stream + halfway
/// bounce-back rules, and the pinned traversal order. Bump on ANY
/// change that can move result bits.
pub const D3Q19_BIT_SEMANTICS_VERSION: u32 = 1;

/// The D3Q19 population count.
pub const Q3: usize = 19;

/// The D3Q19 lattice velocities: rest, 6 face neighbors, 12 edge
/// neighbors — opposites adjacent (`2k-1` ↔ `2k`), so the opposite
/// table is verifiable at a glance.
pub const E3: [(i32, i32, i32); Q3] = [
    (0, 0, 0),
    (1, 0, 0),
    (-1, 0, 0),
    (0, 1, 0),
    (0, -1, 0),
    (0, 0, 1),
    (0, 0, -1),
    (1, 1, 0),
    (-1, -1, 0),
    (1, -1, 0),
    (-1, 1, 0),
    (1, 0, 1),
    (-1, 0, -1),
    (1, 0, -1),
    (-1, 0, 1),
    (0, 1, 1),
    (0, -1, -1),
    (0, 1, -1),
    (0, -1, 1),
];

/// The D3Q19 weights: 1/3 rest, 1/18 per face, 1/36 per edge.
pub const W3: [f64; Q3] = [
    1.0 / 3.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 18.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Integer weights ×36 — the EXACT arithmetic the lattice-invariant
/// tests use (Σ = 36, moments in integers, no float tolerance).
pub const W36: [i64; Q3] = [12, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

/// Opposite-direction indices (for bounce-back): rest is self-opposite,
/// then adjacent pairs.
pub const OPP3: [usize; Q3] = [
    0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17,
];

/// Tile edge in cells (4×4×4 = 64 cells per tile).
pub const TILE: usize = 4;
const TILE_CELLS: usize = TILE * TILE * TILE;

/// One 4×4×4 tile of one scalar field: 512 B, 128-byte aligned — the
/// plan §5 SoA tile-major unit.
#[derive(Clone)]
#[repr(align(128))]
struct Tile([f64; TILE_CELLS]);

impl Tile {
    fn filled(value: f64) -> Tile {
        Tile([value; TILE_CELLS])
    }
}

/// The D3Q19 equilibrium distribution at density `rho` and velocity `u`.
#[must_use]
pub fn equilibrium3(rho: f64, u: [f64; 3]) -> [f64; Q3] {
    let usq = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    let mut f = [0.0; Q3];
    for i in 0..Q3 {
        let (ex, ey, ez) = (f64::from(E3[i].0), f64::from(E3[i].1), f64::from(E3[i].2));
        let eu = ex * u[0] + ey * u[1] + ez * u[2];
        f[i] = W3[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * usq);
    }
    f
}

/// Selectable D3Q19 collision law for the shared per-cell kernel.
///
/// The first rung contains only the frozen BGK law. Central-moment and
/// reduced-cumulant rungs extend this enum without duplicating the cell
/// authority path in [`Duct`] and [`BoundaryGrid3`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CollisionModel3 {
    /// Single-relaxation-time BGK collision.
    Bgk {
        /// Relaxation time; must be finite and greater than one half.
        tau: f64,
    },
}

impl CollisionModel3 {
    /// Validate all declared relaxation parameters.
    ///
    /// # Errors
    /// [`CollisionError3::InvalidRelaxationTime`] when a parameter is outside
    /// the finite positive-viscosity window.
    pub fn validate(self) -> Result<(), CollisionError3> {
        match self {
            Self::Bgk { tau } if tau.is_finite() && tau > 0.5 => Ok(()),
            Self::Bgk { tau } => Err(CollisionError3::InvalidRelaxationTime { tau }),
        }
    }

    fn tau(self) -> f64 {
        match self {
            Self::Bgk { tau } => tau,
        }
    }
}

/// Fail-closed D3Q19 cell-collision diagnostics.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CollisionError3 {
    /// Relaxation time is non-finite or does not yield positive viscosity.
    InvalidRelaxationTime {
        /// Rejected relaxation time.
        tau: f64,
    },
    /// A body-force component is non-finite.
    NonFiniteForce {
        /// Cartesian axis, `0..3`.
        axis: usize,
        /// Rejected component.
        value: f64,
    },
    /// An incoming population is non-finite.
    NonFinitePopulation {
        /// D3Q19 direction index.
        direction: usize,
        /// Rejected population.
        value: f64,
    },
    /// Incoming populations do not define positive finite mass.
    NonPositiveDensity {
        /// Computed density.
        rho: f64,
    },
    /// A force-corrected velocity component is non-finite.
    NonFiniteVelocity {
        /// Cartesian axis, `0..3`.
        axis: usize,
        /// Computed component.
        value: f64,
    },
    /// Collision produced a non-finite outgoing population.
    NonFiniteOutput {
        /// D3Q19 direction index.
        direction: usize,
        /// Computed population.
        value: f64,
    },
}

impl core::fmt::Display for CollisionError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidRelaxationTime { tau } => {
                write!(f, "D3Q19 relaxation time {tau} must be finite and > 0.5")
            }
            Self::NonFiniteForce { axis, value } => {
                write!(f, "D3Q19 force axis {axis} is non-finite ({value})")
            }
            Self::NonFinitePopulation { direction, value } => write!(
                f,
                "D3Q19 incoming population {direction} is non-finite ({value})"
            ),
            Self::NonPositiveDensity { rho } => {
                write!(f, "D3Q19 density must be positive and finite (got {rho})")
            }
            Self::NonFiniteVelocity { axis, value } => {
                write!(f, "D3Q19 velocity axis {axis} is non-finite ({value})")
            }
            Self::NonFiniteOutput { direction, value } => write!(
                f,
                "D3Q19 outgoing population {direction} is non-finite ({value})"
            ),
        }
    }
}

impl core::error::Error for CollisionError3 {}

/// Collide one D3Q19 population vector under an explicit checked model and
/// Guo body force.
///
/// The BGK arithmetic deliberately retains the two existing evaluation paths:
/// this public entry point uses the frozen boundary-grid expression, while
/// [`Duct`] delegates through its private frozen axial-force expression. This
/// extraction therefore centralizes authority without silently refreezing
/// either golden surface.
///
/// # Errors
/// [`CollisionError3`] for inadmissible model parameters, force, incoming
/// state, macroscopic state, or outgoing populations.
pub fn collide_cell3(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    collide_cell3_with_projection(populations, model, force, false)
}

fn collide_axial_z_cell3(
    populations: [f64; Q3],
    model: CollisionModel3,
    gz: f64,
) -> Result<[f64; Q3], CollisionError3> {
    collide_cell3_with_projection(populations, model, [0.0, 0.0, gz], true)
}

fn collide_cell3_with_projection(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
    frozen_axial_z_projection: bool,
) -> Result<[f64; Q3], CollisionError3> {
    model.validate()?;
    for (axis, value) in force.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFiniteForce { axis, value });
        }
    }

    let mut rho = 0.0;
    let mut momentum = [0.0; 3];
    for (direction, value) in populations.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFinitePopulation { direction, value });
        }
        rho += value;
        momentum[0] += f64::from(E3[direction].0) * value;
        momentum[1] += f64::from(E3[direction].1) * value;
        momentum[2] += f64::from(E3[direction].2) * value;
    }
    if !(rho.is_finite() && rho > 0.0) {
        return Err(CollisionError3::NonPositiveDensity { rho });
    }
    let velocity = core::array::from_fn(|axis| (momentum[axis] + 0.5 * force[axis]) / rho);
    for (axis, value) in velocity.into_iter().enumerate() {
        if !value.is_finite() {
            return Err(CollisionError3::NonFiniteVelocity { axis, value });
        }
    }

    let tau = model.tau();
    let equilibrium = equilibrium3(rho, velocity);
    let coefficient = 1.0 - 0.5 / tau;
    let cs2 = 1.0 / 3.0;
    let cs4 = cs2 * cs2;
    let mut post = [0.0; Q3];
    for direction in 0..Q3 {
        let e = [
            f64::from(E3[direction].0),
            f64::from(E3[direction].1),
            f64::from(E3[direction].2),
        ];
        let eu = if frozen_axial_z_projection {
            e[0] * velocity[0] + e[1] * velocity[1] + e[2] * velocity[2]
        } else {
            e.iter()
                .zip(velocity)
                .map(|(component, u)| *component * u)
                .sum::<f64>()
        };
        let forcing = if frozen_axial_z_projection {
            coefficient * W3[direction] * (3.0 * (e[2] - velocity[2]) + 9.0 * eu * e[2]) * force[2]
        } else {
            let projection = (0..3)
                .map(|axis| ((e[axis] - velocity[axis]) / cs2 + eu * e[axis] / cs4) * force[axis])
                .sum::<f64>();
            coefficient * W3[direction] * projection
        };
        let value = populations[direction]
            + (equilibrium[direction] - populations[direction]) / tau
            + forcing;
        if !value.is_finite() {
            return Err(CollisionError3::NonFiniteOutput { direction, value });
        }
        post[direction] = value;
    }
    Ok(post)
}

/// A D3Q19 duct: halfway bounce-back walls on the x and y boundaries,
/// periodic in z, driven by a body force `gz` along z — the 3-D
/// Poiseuille fixture (plan §14.1). Densities start at 1, velocities at
/// rest, unless seeded through [`Duct::perturb`].
pub struct Duct {
    nx: usize,
    ny: usize,
    nz: usize,
    tau: f64,
    gz: f64,
    /// SoA distributions: one tile-major field per population.
    f: [Vec<Tile>; Q3],
    /// Post-collision scratch (pull streaming reads this).
    post: [Vec<Tile>; Q3],
}

impl Duct {
    /// A duct at rest (unit density) with relaxation time `tau` and body
    /// force `gz`. Every dimension must be a positive multiple of
    /// [`TILE`].
    ///
    /// # Panics
    /// If any dimension is zero or not a multiple of [`TILE`], if `tau` is not
    /// finite and greater than one half, or if `gz` is not finite.
    #[must_use]
    pub fn new(nx: usize, ny: usize, nz: usize, tau: f64, gz: f64) -> Duct {
        assert!(
            nx > 0
                && ny > 0
                && nz > 0
                && nx.is_multiple_of(TILE)
                && ny.is_multiple_of(TILE)
                && nz.is_multiple_of(TILE),
            "duct dimensions must be positive multiples of {TILE} (got {nx}x{ny}x{nz})"
        );
        CollisionModel3::Bgk { tau }
            .validate()
            .expect("duct relaxation time must be finite and greater than 0.5");
        assert!(gz.is_finite(), "duct body force must be finite");
        let tiles = (nx / TILE) * (ny / TILE) * (nz / TILE);
        let f0 = equilibrium3(1.0, [0.0; 3]);
        let f = core::array::from_fn(|i| vec![Tile::filled(f0[i]); tiles]);
        let post = core::array::from_fn(|i| vec![Tile::filled(f0[i]); tiles]);
        Duct {
            nx,
            ny,
            nz,
            tau,
            gz,
            f,
            post,
        }
    }

    /// Tile index and local lane of cell `(x, y, z)` — the tile-major
    /// address map (x-fastest at both levels; the pinned order).
    #[inline]
    fn addr(&self, x: usize, y: usize, z: usize) -> (usize, usize) {
        let (ntx, nty) = (self.nx / TILE, self.ny / TILE);
        let tile = (z / TILE * nty + y / TILE) * ntx + x / TILE;
        let lane = (z % TILE * TILE + y % TILE) * TILE + x % TILE;
        (tile, lane)
    }

    /// The kinematic viscosity `ν = (τ − ½)/3`.
    #[must_use]
    pub fn viscosity(&self) -> f64 {
        (self.tau - 0.5) / 3.0
    }

    /// The macroscopic density at `(x, y, z)`.
    #[must_use]
    pub fn density(&self, x: usize, y: usize, z: usize) -> f64 {
        let (tile, lane) = self.addr(x, y, z);
        (0..Q3).map(|i| self.f[i][tile].0[lane]).sum()
    }

    /// The macroscopic velocity at `(x, y, z)` (with the Guo half-force
    /// momentum correction).
    #[must_use]
    pub fn velocity(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        let (tile, lane) = self.addr(x, y, z);
        let mut rho = 0.0;
        let mut m = [0.0; 3];
        for (i, e) in E3.iter().enumerate() {
            let fi = self.f[i][tile].0[lane];
            rho += fi;
            m[0] += f64::from(e.0) * fi;
            m[1] += f64::from(e.1) * fi;
            m[2] += f64::from(e.2) * fi;
        }
        // Guo: half the force is added to the momentum.
        [m[0] / rho, m[1] / rho, (m[2] + 0.5 * self.gz) / rho]
    }

    /// Total mass (conserved by construction; drift is roundoff only).
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        self.f
            .iter()
            .flat_map(|field| field.iter())
            .map(|tile| tile.0.iter().sum::<f64>())
            .sum()
    }

    /// Deterministically perturb the resting state: cell `(x, y, z)`
    /// gets its density shifted by `amplitude · h(x, y, z)` where `h`
    /// is a fixed integer hash mapped to `[-1, 1)` — the seeded golden
    /// fixture's initial condition (no RNG dependency; the hash IS the
    /// seed schedule).
    pub fn perturb(&mut self, seed: u64, amplitude: f64) {
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    let mut h = seed
                        ^ (x as u64)
                            .wrapping_mul(0x9e37_79b9_7f4a_7c15)
                            .wrapping_add((y as u64).wrapping_mul(0xbf58_476d_1ce4_e5b9))
                            .wrapping_add((z as u64).wrapping_mul(0x94d0_49bb_1331_11eb));
                    h ^= h >> 30;
                    h = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
                    h ^= h >> 27;
                    // map to [-1, 1)
                    let unit = (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0;
                    let rho = 1.0 + amplitude * unit;
                    let feq = equilibrium3(rho, [0.0; 3]);
                    let (tile, lane) = self.addr(x, y, z);
                    for (field, feq_i) in self.f.iter_mut().zip(feq) {
                        field[tile].0[lane] = feq_i;
                    }
                }
            }
        }
    }

    /// One collide-force-stream step (BGK + Guo forcing, pull scheme,
    /// halfway bounce-back x/y walls, periodic z). Traversal follows the
    /// pinned tile-major order documented at module level.
    pub fn step(&mut self) {
        // Collide + Guo forcing (pointwise): write post-collision
        // populations into `post`, visiting tiles then lanes ascending.
        let tiles = self.f[0].len();
        for tile in 0..tiles {
            for lane in 0..TILE_CELLS {
                let populations = core::array::from_fn(|direction| self.f[direction][tile].0[lane]);
                let post = collide_axial_z_cell3(
                    populations,
                    CollisionModel3::Bgk { tau: self.tau },
                    self.gz,
                )
                .expect("Duct constructor and prior collision state admit BGK/Guo collision");
                for (field, value) in self.post.iter_mut().zip(post) {
                    field[tile].0[lane] = value;
                }
            }
        }
        // Pull streaming: destination (x,y,z) reads post[source] where
        // source = destination − e_i; crossing an x or y boundary means
        // the population came off a wall — halfway bounce-back reflects
        // the OPPOSITE post-collision population of the destination cell
        // itself; z wraps (periodic). Offsets are ±1 only.
        for z in 0..self.nz {
            for y in 0..self.ny {
                for x in 0..self.nx {
                    let (dtile, dlane) = self.addr(x, y, z);
                    for i in 0..Q3 {
                        let (ex, ey, ez) = E3[i];
                        let sx = match ex {
                            1 if x == 0 => None,
                            1 => Some(x - 1),
                            -1 if x + 1 == self.nx => None,
                            -1 => Some(x + 1),
                            _ => Some(x),
                        };
                        let sy = match ey {
                            1 if y == 0 => None,
                            1 => Some(y - 1),
                            -1 if y + 1 == self.ny => None,
                            -1 => Some(y + 1),
                            _ => Some(y),
                        };
                        let sz = match ez {
                            1 => (z + self.nz - 1) % self.nz,
                            -1 => (z + 1) % self.nz,
                            _ => z,
                        };
                        self.f[i][dtile].0[dlane] = match (sx, sy) {
                            (Some(sx), Some(sy)) => {
                                let (stile, slane) = self.addr(sx, sy, sz);
                                self.post[i][stile].0[slane]
                            }
                            // wall crossing: bounce back in place.
                            _ => self.post[OPP3[i]][dtile].0[dlane],
                        };
                    }
                }
            }
        }
    }

    /// Run `steps` steps.
    pub fn run(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    /// The `z`-velocity over the duct cross-section at `z = 0`, row-major
    /// in `(y, x)` — the profile the analytic duct series certifies.
    #[must_use]
    pub fn z_velocity_section(&self) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.nx * self.ny);
        for y in 0..self.ny {
            for x in 0..self.nx {
                out.push(self.velocity(x, y, 0)[2]);
            }
        }
        out
    }
}

impl core::fmt::Debug for Duct {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Duct")
            .field("nx", &self.nx)
            .field("ny", &self.ny)
            .field("nz", &self.nz)
            .field("tau", &self.tau)
            .field("gz", &self.gz)
            .finish_non_exhaustive()
    }
}

/// The analytic steady rectangular-duct `z`-velocity at lattice cell
/// `(x, y)` for a duct of `nx × ny` cells under body acceleration `gz`,
/// with halfway bounce-back walls at `x = −½`, `x = nx − ½`, `y = −½`,
/// `y = ny − ½` (so the full widths are exactly `nx` and `ny`):
///
/// `u(X,Y) = (16 gz a²)/(ν π³) Σ_{n odd} (−1)^((n−1)/2) n⁻³
///           [1 − cosh(nπY/2a)/cosh(nπb/2a)] cos(nπX/2a)`
///
/// with `a = nx/2`, `b = ny/2`, and `(X, Y)` the cell center relative to
/// the duct center. The series is truncated at `n = 99` (the 1/n³ decay
/// puts the tail below 1e-6 relative — far under the 3% acceptance bar);
/// aspect ratios up to ~4 stay clear of `cosh` overflow.
#[must_use]
pub fn duct_analytic(gz: f64, viscosity: f64, nx: usize, ny: usize, x: usize, y: usize) -> f64 {
    let a = nx as f64 / 2.0;
    let b = ny as f64 / 2.0;
    let cx = x as f64 - (nx as f64 - 1.0) / 2.0;
    let cy = y as f64 - (ny as f64 - 1.0) / 2.0;
    let mut sum = 0.0;
    let mut sign = 1.0;
    let mut n = 1u32;
    while n <= 99 {
        let nf = f64::from(n);
        // n³ and π³ as explicit products — `powi` is the build-mode
        // determinism hazard class (det:: doctrine / check-powi lint).
        let k = nf * core::f64::consts::PI / (2.0 * a);
        let term = (1.0 - (k * cy).cosh() / (k * b).cosh()) * (k * cx).cos() / (nf * nf * nf);
        sum += sign * term;
        sign = -sign;
        n += 2;
    }
    let pi = core::f64::consts::PI;
    16.0 * gz * a * a / (viscosity * (pi * pi * pi)) * sum
}
