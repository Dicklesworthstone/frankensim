//! Free-surface 3-D on sparse active tiles (bead sxnm, plan §15.3): the
//! pour primitive's physics layer.
//!
//! This is a FAITHFUL PORT of the 2-D Körner-style VOF machinery in
//! `crate::freesurface` onto the D3Q19 sparse active-tile grid — the same
//! cell taxonomy (fluid / interface / gas), the same pairwise-antisymmetric
//! mass exchange, the same gas-side population reconstruction, the same
//! conversion cascade with a conservation carry, and the same
//! surface-tension reference-density port. Resolved interfaces use volume-
//! fraction height functions; incomplete stencils retain smoothed-fill
//! curvature with the contact-model wall ghost. An opt-in planar equilibrium
//! contact-angle boundary supplies a geometric wall ghost to that same
//! curvature and pressure-reconstruction operator.
//!
//! Sparse specifics owned by this layer:
//! - Inactive tiles and out-of-domain space are WALL (exactly the sparse
//!   sweep's convention); interface cells are inserted only between fluid
//!   and gas.
//! - Tile activation as fluid advances: when the conversion cascade must
//!   promote a gas cell that lives in an inactive in-domain tile, that tile
//!   is activated first (equilibrium-initialized, all lanes gas, zero
//!   tracked mass) — the WS1-D hook this bead exists to consume.
//! - Per-lane classification and tracked mass are keyed by Morton key, so
//!   slot reshuffles from activation cannot move free-surface state.
//!
//! Determinism: every sweep runs in ascending Morton-key order and
//! ascending lane order; conversion cascades process cells in collection
//! order; ties cannot depend on activation history (the sparse layer
//! guarantees set-only ordering).
//!
//! No-claims: this layer claims mass-ledger conservation (tracked per step,
//! worst violation exposed) and deterministic replay — not throughput, not
//! turbulence physics or calibrated contact-line dynamics. Height functions
//! require a resolved single interface; near walls and unresolved topology
//! the smoothed-fill closure remains. The imposed equilibrium-angle
//! closure is limited to one planar domain wall; it supplies no hysteresis,
//! thermocapillary stress, moving substrate or surface-energy balance.

use std::collections::BTreeMap;

use fs_exec::CancelGate;

use super::sparse::{SparseError3, SparseGrid3, demorton3, morton3};
use super::{CollisionError3, E3, OPP3, Q3, TILE, equilibrium3};

mod initial;
pub use initial::FreeSurfaceInitialCell3;
mod curvature;
mod forcing;
pub(crate) mod thermocapillary;
pub(crate) mod transport;
pub(crate) mod wetting;
use transport::{MassLocation3, MassTransport3};

/// Free-surface mass/geometry semantics, independent of the D3Q19 core.
/// Version 2 replaces the legacy global conversion pool with neighbor-local
/// transfers. This intentionally changes evolving fill, pressure and pour bits;
/// the core collision/streaming surface and its goldens are unchanged.
pub const FREE_SURFACE3_BIT_SEMANTICS_VERSION: u32 = 2;

/// Conversion semantics for new simulations versus historical trajectory replay.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum ConversionRedistribution3 {
    /// Transfer excess and debt only to surviving neighboring interface cells.
    #[default]
    Local,
    /// Reproduce version-1 global pooling, including its nonlocal mass and heat
    /// mixing between disconnected drops. Intended only for historical replay.
    LegacyGlobal,
}

/// D3Q19 lattice speed of sound squared.
const CS2_3: f64 = 1.0 / 3.0;

/// Cells per tile (4×4×4).
const TILE_CELLS: usize = TILE * TILE * TILE;

/// Conversion hysteresis band, identical to the 2-D crate's `EPS`.
const EPS: f64 = 1e-3;

/// Per-cell free-surface classification (walls are NOT a cell state here:
/// inactive tiles and out-of-domain space are wall by sparse convention).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Cell3 {
    /// Fully wetted; mass is `Σf` implicitly.
    Fluid,
    /// Partially filled; mass tracked explicitly.
    Interface,
    /// Empty; carries no mass and its populations are never read.
    Gas,
}

/// Contact-line model for the wall ghost of the fill field — the 2-D
/// crate's `ContactModel`, re-declared here so the sparse layer does not
/// depend on the dense 2-D module.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactModel3 {
    /// Wall mirrors the querying cell's fill (≈ 90° contact angle).
    Neutral,
    /// Wall reads as fully wet (spreading).
    Wetting,
}

/// One planar domain wall carrying an imposed equilibrium contact angle.
/// Internal solid tiles and the other five domain walls retain their existing
/// neutral ghost; this does not infer normals for a staircase substrate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlanarWall3 {
    /// Minimum x face, inward wall normal +x.
    XMin,
    /// Maximum x face, inward wall normal -x.
    XMax,
    /// Minimum y face, inward wall normal +y.
    YMin,
    /// Maximum y face, inward wall normal -y.
    YMax,
    /// Minimum z face, inward wall normal +z.
    ZMin,
    /// Maximum z face, inward wall normal -z.
    ZMax,
}

impl PlanarWall3 {
    fn axis(self) -> usize {
        match self {
            Self::XMin | Self::XMax => 0,
            Self::YMin | Self::YMax => 1,
            Self::ZMin | Self::ZMax => 2,
        }
    }

    fn upper(self) -> bool {
        matches!(self, Self::XMax | Self::YMax | Self::ZMax)
    }
}

#[derive(Debug, Clone, Copy)]
struct PlanarContactAngle3 {
    wall: PlanarWall3,
    angle_radians: f64,
    cotangent: f64,
}

/// Cumulative conversion ledger (port of the 2-D `ConversionStats`).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ConversionStats3 {
    /// Interface → fluid events.
    pub to_fluid: u64,
    /// Interface → gas events.
    pub to_gas: u64,
    /// Gas → interface events (closure repair).
    pub gas_to_interface: u64,
    /// Fluid → interface events (closure repair).
    pub fluid_to_interface: u64,
    /// Tiles activated because fluid advanced into them.
    pub tiles_activated: u64,
}

/// Typed refusal from free-surface construction or a step.
#[derive(Debug, Clone, PartialEq)]
pub enum FreeSurfaceError3 {
    /// The underlying sparse grid refused.
    Grid(SparseError3),
    /// The shared per-cell collision kernel refused at a wet cell.
    Collision {
        /// Morton key of the refusing tile.
        tile_key: u64,
        /// Lane (0..64) of the refusing cell.
        lane: usize,
        /// The underlying refusal.
        source: CollisionError3,
    },
    /// Construction found no fluid cell (an all-gas free surface is a
    /// degenerate fixture, refused rather than silently simulated).
    NoFluid,
    /// The imposed angle or its finite-difference ghost slope is inadmissible.
    InvalidContactAngle,
    /// Initial liquid state cannot be represented without changing its mass
    /// or physical velocity. No candidate is published on this refusal.
    InvalidInitialState {
        /// Offending cell, if the refusal concerns one supplied cell.
        coordinates: Option<[usize; 3]>,
        /// Failed admission condition.
        reason: &'static str,
    },
    /// A bounded companion transport could not retain the native mass routes.
    Transport(&'static str),
    /// Cancellation was observed before physical state publication.
    Cancelled,
}

impl core::fmt::Display for FreeSurfaceError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            FreeSurfaceError3::Grid(err) => write!(f, "sparse grid refused: {err}"),
            FreeSurfaceError3::Collision {
                tile_key,
                lane,
                source,
            } => write!(
                f,
                "free-surface collision refused at tile key {tile_key:#x} lane {lane}: {source}"
            ),
            FreeSurfaceError3::NoFluid => f.write_str("free-surface fixture contains no fluid"),
            FreeSurfaceError3::InvalidContactAngle => f.write_str(
                "contact angle must lie in (0, pi) radians with a representable wall slope",
            ),
            FreeSurfaceError3::InvalidInitialState {
                coordinates,
                reason,
            } => {
                write!(
                    f,
                    "free-surface initial state at {coordinates:?} refused: {reason}"
                )
            }
            FreeSurfaceError3::Transport(reason) => {
                write!(f, "free-surface companion transport refused: {reason}")
            }
            FreeSurfaceError3::Cancelled => {
                f.write_str("free-surface step cancelled before publication")
            }
        }
    }
}

impl std::error::Error for FreeSurfaceError3 {}

impl From<SparseError3> for FreeSurfaceError3 {
    fn from(err: SparseError3) -> Self {
        FreeSurfaceError3::Grid(err)
    }
}

/// Per-tile free-surface state, keyed by Morton key (slot-reshuffle proof).
#[derive(Clone)]
struct FsTile {
    cells: [Cell3; TILE_CELLS],
    mass: [f64; TILE_CELLS],
}

impl FsTile {
    fn all_gas() -> FsTile {
        FsTile {
            cells: [Cell3::Gas; TILE_CELLS],
            mass: [0.0; TILE_CELLS],
        }
    }
}

/// Wet-cell state from the native sparse populations and tracked fill.
///
/// Velocity includes the same Guo half-force correction used by collision.
/// Density and velocity are numerical lattice quantities, not a gas EOS.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeSurfaceSample3 {
    /// Fluid or interface; gas cells have no liquid sample.
    pub cell: Cell3,
    /// Native tracked liquid fill fraction.
    pub fill_fraction: f64,
    /// Sum of the wet cell's populations (reference density one).
    pub density: f64,
    /// Lattice velocity including half the actual retained Guo force density.
    pub velocity: [f64; 3],
}

/// 3-D free surface over a sparse active-tile D3Q19 grid.
#[derive(Clone)]
pub struct FreeSurface3 {
    grid: SparseGrid3,
    tiles: BTreeMap<u64, FsTile>,
    /// Permanently solid tiles (tank walls, fixtures): advance-driven
    /// activation refuses to open these — they are wall FOREVER, unlike
    /// plain inactive tiles which are merely not-yet-reached gas.
    solid: std::collections::BTreeSet<u64>,
    /// Surface tension coefficient (0 = off), 2-D port.
    sigma: f64,
    /// Source-evolved coefficients on the current wet cells. The original
    /// scalar path remains unchanged when all coefficients equal `sigma`.
    local_sigma: Option<BTreeMap<(u64, usize), f64>>,
    redistribution: ConversionRedistribution3,
    contact: ContactModel3,
    planar_contact: Option<PlanarContactAngle3>,
    planar_contact_field: Option<wetting::PlanarContactField3>,
    /// In Local mode, unresolved mass with no adjacent interface recipient.
    /// It remains in the ledger; no location is available to redistribute it.
    carry: f64,
    conversions: ConversionStats3,
    /// Worst single-step relative ledger violation observed so far —
    /// the bead demands this be logged, not assumed.
    worst_step_violation: f64,
    /// Smoothed fill field, refreshed at the top of each step.
    fill_smooth: BTreeMap<u64, [f64; TILE_CELLS]>,
    /// Actual force frame of each wet population, when it differs from the
    /// uniform reference gravity. This history travels with population state.
    local_force: Option<BTreeMap<(u64, usize), [f64; 3]>>,
}

impl core::fmt::Debug for FreeSurface3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("FreeSurface3")
            .field("active_tiles", &self.tiles.len())
            .field("sigma", &self.sigma)
            .field("redistribution", &self.redistribution)
            .field("contact", &self.contact)
            .field("planar_contact", &self.planar_contact)
            .field("carry", &self.carry)
            .field("conversions", &self.conversions)
            .field("worst_step_violation", &self.worst_step_violation)
            .finish_non_exhaustive()
    }
}

/// Local (x, y, z) of a lane, inverse of the sparse layer's lane packing.
fn lane_coords(lane: usize) -> (usize, usize, usize) {
    (lane % TILE, (lane / TILE) % TILE, lane / (TILE * TILE))
}

fn lane_of(lx: usize, ly: usize, lz: usize) -> usize {
    (lz * TILE + ly) * TILE + lx
}

impl FreeSurface3 {
    /// Bound shared by the height-function and fallback curvature estimates.
    /// This numerical clamp is not an assertion about physical interface shape.
    pub const CURVATURE_LIMIT: f64 = 1.0;

    /// Sample a wet cell by global coordinates, using the native population
    /// moments and Guo force convention. Gas, wall and out-of-domain return None.
    /// This is a read-only diagnostic; it cannot change the solver trajectory.
    #[must_use]
    pub fn sample_at(&self, x: usize, y: usize, z: usize) -> Option<FreeSurfaceSample3> {
        let (slot, lane) = self.grid.resolve_source(
            i64::try_from(x).ok()?,
            i64::try_from(y).ok()?,
            i64::try_from(z).ok()?,
        )?;
        let key = self.grid.active_keys()[slot];
        let cell = self.tiles.get(&key)?.cells[lane];
        if cell == Cell3::Gas {
            return None;
        }
        let (density, raw_velocity) = self.grid.cell_macros(slot, lane);
        let force = self.force_at_lane(key, lane);
        Some(FreeSurfaceSample3 {
            cell,
            fill_fraction: self.fill(key, lane),
            density,
            velocity: core::array::from_fn(|axis| raw_velocity[axis] + 0.5 * force[axis] / density),
        })
    }

    /// Build a free surface over an `nx × ny × nz` domain: `fluid`
    /// classifies each global cell; tiles containing fluid plus their
    /// in-domain neighbor margin are activated (the margin holds the gas
    /// the interface will advance into); everything else stays inactive
    /// (wall). Interface cells are inserted between fluid and gas exactly
    /// as in 2-D (fluid cell with ≥1 gas neighbor), masses initialized
    /// (fluid `Σf`, interface `Σf/2`), and fluid/gas closure is asserted.
    ///
    /// # Errors
    /// [`FreeSurfaceError3::Grid`] for inadmissible dims,
    /// [`FreeSurfaceError3::NoFluid`] for an all-gas fixture.
    ///
    /// # Panics
    /// If interface insertion fails to separate fluid from gas
    /// (impossible by construction).
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        tau: f64,
        gravity: [f64; 3],
        sigma: f64,
        contact: ContactModel3,
        fluid: impl Fn(usize, usize, usize) -> bool,
    ) -> Result<FreeSurface3, FreeSurfaceError3> {
        Self::with_solid_tiles(nx, ny, nz, tau, gravity, sigma, contact, &[], fluid)
    }

    /// [`FreeSurface3::new`] with permanently solid tiles (tank walls and
    /// fixtures at tile granularity): solid tiles are never activated —
    /// not at construction, not by fluid advance — so they behave as
    /// walls for streaming, exchange, and the fill-field ghost alike.
    ///
    /// # Errors
    /// As [`FreeSurface3::new`]; also refuses nonfinite/negative tension or a
    /// coefficient that can make the clamped capillary boundary density
    /// nonpositive via [`FreeSurfaceError3::InvalidInitialState`].
    ///
    /// # Panics
    /// If interface insertion fails to separate fluid from gas, or if a
    /// fluid cell is declared inside a solid tile (contradictory fixture).
    #[allow(clippy::too_many_arguments)]
    pub fn with_solid_tiles(
        nx: usize,
        ny: usize,
        nz: usize,
        tau: f64,
        gravity: [f64; 3],
        sigma: f64,
        contact: ContactModel3,
        solid_tiles: &[(u32, u32, u32)],
        fluid: impl Fn(usize, usize, usize) -> bool,
    ) -> Result<FreeSurface3, FreeSurfaceError3> {
        if !sigma.is_finite() || sigma < 0.0 || sigma * Self::CURVATURE_LIMIT / CS2_3 >= 1.0 {
            return Err(FreeSurfaceError3::InvalidInitialState {
                coordinates: None,
                reason: "capillary boundary density can be nonpositive",
            });
        }
        let mut grid = SparseGrid3::new(nx, ny, nz, tau, gravity)?;
        let solid: std::collections::BTreeSet<u64> = solid_tiles
            .iter()
            .map(|&(tx, ty, tz)| morton3(tx, ty, tz))
            .collect();
        let (ntx, nty, ntz) = grid.tile_dims();

        // Fluid tiles + in-domain neighbor margin.
        let mut fluid_tiles: Vec<(u32, u32, u32)> = Vec::new();
        for tz in 0..ntz {
            for ty in 0..nty {
                for tx in 0..ntx {
                    let has_fluid = (0..TILE_CELLS).any(|lane| {
                        let (lx, ly, lz) = lane_coords(lane);
                        fluid(tx * TILE + lx, ty * TILE + ly, tz * TILE + lz)
                    });
                    if has_fluid {
                        #[allow(clippy::cast_possible_truncation)]
                        let coords = (tx as u32, ty as u32, tz as u32);
                        assert!(
                            !solid.contains(&morton3(coords.0, coords.1, coords.2)),
                            "fixture declares fluid inside solid tile {coords:?}"
                        );
                        fluid_tiles.push(coords);
                    }
                }
            }
        }
        if fluid_tiles.is_empty() {
            return Err(FreeSurfaceError3::NoFluid);
        }
        let mut active: Vec<(u32, u32, u32)> = Vec::new();
        for &(tx, ty, tz) in &fluid_tiles {
            for dz in -1i64..=1 {
                for dy in -1i64..=1 {
                    for dx in -1i64..=1 {
                        let (ax, ay, az) =
                            (i64::from(tx) + dx, i64::from(ty) + dy, i64::from(tz) + dz);
                        if ax >= 0
                            && ay >= 0
                            && az >= 0
                            && (ax as usize) < ntx
                            && (ay as usize) < nty
                            && (az as usize) < ntz
                        {
                            #[allow(clippy::cast_possible_truncation)]
                            let coords = (ax as u32, ay as u32, az as u32);
                            if !solid.contains(&morton3(coords.0, coords.1, coords.2)) {
                                active.push(coords);
                            }
                        }
                    }
                }
            }
        }
        active.sort_unstable();
        active.dedup();
        grid.activate_tiles(&active)?;

        // Classify lanes; populations are already at rest equilibrium.
        let mut tiles: BTreeMap<u64, FsTile> = BTreeMap::new();
        for &key in grid.active_keys() {
            let (tx, ty, tz) = demorton3(key);
            let mut tile = FsTile::all_gas();
            for lane in 0..TILE_CELLS {
                let (lx, ly, lz) = lane_coords(lane);
                if fluid(
                    tx as usize * TILE + lx,
                    ty as usize * TILE + ly,
                    tz as usize * TILE + lz,
                ) {
                    tile.cells[lane] = Cell3::Fluid;
                }
            }
            tiles.insert(key, tile);
        }

        let mut fs = FreeSurface3 {
            grid,
            tiles,
            solid,
            sigma,
            local_sigma: None,
            redistribution: ConversionRedistribution3::Local,
            contact,
            planar_contact: None,
            planar_contact_field: None,
            carry: 0.0,
            conversions: ConversionStats3::default(),
            worst_step_violation: 0.0,
            fill_smooth: BTreeMap::new(),
            local_force: None,
        };

        // Interface insertion: any fluid cell with a gas neighbor.
        let mut promote: Vec<(u64, usize)> = Vec::new();
        for (&key, tile) in &fs.tiles {
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] != Cell3::Fluid {
                    continue;
                }
                for q in 1..Q3 {
                    if matches!(
                        fs.neighbor_state(key, lane, q),
                        Neighbor::Active(_, _, Cell3::Gas)
                    ) {
                        promote.push((key, lane));
                        break;
                    }
                }
            }
        }
        for &(key, lane) in &promote {
            fs.tiles.get_mut(&key).expect("promoted key active").cells[lane] = Cell3::Interface;
        }

        // Mass init: fluid = Σf, interface = Σf/2.
        for &key in &fs.tiles.keys().copied().collect::<Vec<_>>() {
            let slot = fs.grid.slot_of(key).expect("active");
            let tile = fs.tiles.get_mut(&key).expect("active");
            for lane in 0..TILE_CELLS {
                let rho: f64 = fs.grid.populations(slot, lane).iter().sum();
                tile.mass[lane] = match tile.cells[lane] {
                    Cell3::Fluid => rho,
                    Cell3::Interface => 0.5 * rho,
                    Cell3::Gas => 0.0,
                };
            }
        }
        fs.assert_closure();
        Ok(fs)
    }

    fn assert_closure(&self) {
        for (&key, tile) in &self.tiles {
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] != Cell3::Fluid {
                    continue;
                }
                for q in 1..Q3 {
                    assert!(
                        !matches!(
                            self.neighbor_state(key, lane, q),
                            Neighbor::Active(_, _, Cell3::Gas)
                        ),
                        "closure violated: fluid touches gas at key {key:#x} lane {lane}"
                    );
                }
            }
        }
    }

    /// Fill fraction (fluid 1, gas 0, interface `m/ρ` clamped) — 2-D port.
    #[must_use]
    pub fn fill(&self, key: u64, lane: usize) -> f64 {
        let tile = &self.tiles[&key];
        match tile.cells[lane] {
            Cell3::Fluid => 1.0,
            Cell3::Interface => {
                let slot = self.grid.slot_of(key).expect("active");
                let rho: f64 = self.grid.populations(slot, lane).iter().sum();
                (tile.mass[lane] / rho.max(1e-12)).clamp(0.0, 1.0)
            }
            Cell3::Gas => 0.0,
        }
    }

    /// The strict ledger: `Σ_fluid Σf + Σ_interface m + carry` — 2-D port.
    #[must_use]
    pub fn ledger_mass(&self) -> f64 {
        let mut total = self.carry;
        for (&key, tile) in &self.tiles {
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                match tile.cells[lane] {
                    Cell3::Fluid => {
                        total += self.grid.populations(slot, lane).iter().sum::<f64>();
                    }
                    Cell3::Interface => total += tile.mass[lane],
                    Cell3::Gas => {}
                }
            }
        }
        total
    }

    /// Ledger mass restricted to cells whose GLOBAL coordinates satisfy
    /// `region` (fluid Σf + interface m; the global carry is NOT included
    /// — read it via [`FreeSurface3::carry`] when partitioning the full
    /// ledger, since redistribution mass belongs to no one region).
    #[must_use]
    pub fn region_mass(&self, region: impl Fn(usize, usize, usize) -> bool) -> f64 {
        let mut total = 0.0;
        for (&key, tile) in &self.tiles {
            let slot = self.grid.slot_of(key).expect("active");
            let (tx, ty, tz) = demorton3(key);
            for lane in 0..TILE_CELLS {
                let (lx, ly, lz) = lane_coords(lane);
                if !region(
                    tx as usize * TILE + lx,
                    ty as usize * TILE + ly,
                    tz as usize * TILE + lz,
                ) {
                    continue;
                }
                match tile.cells[lane] {
                    Cell3::Fluid => {
                        total += self.grid.populations(slot, lane).iter().sum::<f64>();
                    }
                    Cell3::Interface => total += tile.mass[lane],
                    Cell3::Gas => {}
                }
            }
        }
        total
    }

    /// The conversion-conservation carry (global, unregioned). Orphan mass or
    /// debt stays here in Local mode, without redistribution along unknown routes.
    #[must_use]
    pub fn carry(&self) -> f64 {
        self.carry
    }

    /// Select conversion semantics explicitly. Constructors use [`ConversionRedistribution3::Local`].
    /// The legacy choice exists to replay previously frozen version-1 results;
    /// it does not provide spatially local drop or thermal transport.
    #[must_use]
    pub fn with_conversion_redistribution(mut self, model: ConversionRedistribution3) -> Self {
        self.redistribution = model;
        self
    }

    /// Cumulative conversion statistics.
    #[must_use]
    pub fn conversions(&self) -> ConversionStats3 {
        self.conversions
    }

    /// Worst single-step relative ledger violation observed so far.
    #[must_use]
    pub fn worst_step_violation(&self) -> f64 {
        self.worst_step_violation
    }

    /// Cell state at (key, lane).
    #[must_use]
    pub fn cell(&self, key: u64, lane: usize) -> Cell3 {
        self.tiles[&key].cells[lane]
    }

    /// The wrapped sparse grid (read-only view).
    #[must_use]
    pub fn grid(&self) -> &SparseGrid3 {
        &self.grid
    }

    /// Impose a homogeneous equilibrium angle through the liquid at one planar
    /// domain face. Other walls use the neutral closure. The angle is not a
    /// dynamic contact-line law or a claim of calibrated sessile-drop accuracy.
    /// Exactly pi/2 preserves the existing neutral trajectory bit for bit.
    ///
    /// # Errors
    /// Refuses invalid angles before changing the boundary condition.
    pub fn set_planar_contact_angle(
        &mut self,
        wall: PlanarWall3,
        angle_radians: f64,
    ) -> Result<(), FreeSurfaceError3> {
        if !angle_radians.is_finite()
            || angle_radians <= 0.0
            || angle_radians >= core::f64::consts::PI
        {
            return Err(FreeSurfaceError3::InvalidContactAngle);
        }
        let cotangent = if angle_radians == core::f64::consts::FRAC_PI_2 {
            0.0
        } else {
            let (sin, cos) = fs_math::det::sin_cos(angle_radians);
            cos / sin
        };
        let stencil_bound = 4.0 * cotangent;
        if !stencil_bound.is_finite() || !(stencil_bound * stencil_bound).is_finite() {
            return Err(FreeSurfaceError3::InvalidContactAngle);
        }
        self.planar_contact = Some(PlanarContactAngle3 {
            wall,
            angle_radians,
            cotangent,
        });
        self.planar_contact_field = None;
        self.contact = ContactModel3::Neutral;
        Ok(())
    }

    /// The explicitly assigned domain wall and equilibrium angle [rad].
    #[must_use]
    pub fn planar_contact_angle(&self) -> Option<(PlanarWall3, f64)> {
        self.planar_contact.map(|c| (c.wall, c.angle_radians))
    }

    /// Highest global z coordinate carrying fluid or interface — the 3-D
    /// front-position probe for Martin-Moyce style fixtures reads an axis
    /// extreme; callers pick the axis by fixture orientation.
    #[must_use]
    pub fn wet_extent(&self) -> Option<(i64, i64, i64)> {
        let mut extent: Option<(i64, i64, i64)> = None;
        for (&key, tile) in &self.tiles {
            let (tx, ty, tz) = demorton3(key);
            for lane in 0..TILE_CELLS {
                if matches!(tile.cells[lane], Cell3::Fluid | Cell3::Interface) {
                    let (lx, ly, lz) = lane_coords(lane);
                    let g = (
                        i64::from(tx) * TILE as i64 + lx as i64,
                        i64::from(ty) * TILE as i64 + ly as i64,
                        i64::from(tz) * TILE as i64 + lz as i64,
                    );
                    extent = Some(match extent {
                        None => g,
                        Some((mx, my, mz)) => (mx.max(g.0), my.max(g.1), mz.max(g.2)),
                    });
                }
            }
        }
        extent
    }

    /// Surge-front x position: the maximum global x among fluid/interface
    /// cells with `z < z_below` — the 3-D analog of the 2-D battery's
    /// bottom-row `surge_front` probe (restricting to the bottom slab
    /// keeps splashes out of the front measurement).
    #[must_use]
    pub fn surge_front_x(&self, z_below: usize) -> Option<i64> {
        let mut front: Option<i64> = None;
        for (&key, tile) in &self.tiles {
            let (tx, ty, tz) = demorton3(key);
            let _ = ty;
            for lane in 0..TILE_CELLS {
                if matches!(tile.cells[lane], Cell3::Fluid | Cell3::Interface) {
                    let (lx, _ly, lz) = lane_coords(lane);
                    let gz = tz as usize * TILE + lz;
                    if gz < z_below {
                        let gx = i64::from(tx) * TILE as i64 + lx as i64;
                        front = Some(front.map_or(gx, |f| f.max(gx)));
                    }
                }
            }
        }
        front
    }

    /// Neighbor of (key, lane) in direction `q` (see [`Neighbor`]).
    fn neighbor_state(&self, key: u64, lane: usize, q: usize) -> Neighbor {
        let (tx, ty, tz) = demorton3(key);
        let (lx, ly, lz) = lane_coords(lane);
        let gx = i64::from(tx) * TILE as i64 + lx as i64 + i64::from(E3[q].0);
        let gy = i64::from(ty) * TILE as i64 + ly as i64 + i64::from(E3[q].1);
        let gz = i64::from(tz) * TILE as i64 + lz as i64 + i64::from(E3[q].2);
        match self.grid.resolve_source(gx, gy, gz) {
            None => Neighbor::Wall,
            Some((slot, nb_lane)) => {
                let nb_key = self.grid.active_keys()[slot];
                Neighbor::Active(nb_key, nb_lane, self.tiles[&nb_key].cells[nb_lane])
            }
        }
    }

    /// Global coordinates of the tile+lane cell.
    fn global_of(key: u64, lane: usize) -> (i64, i64, i64) {
        let (tx, ty, tz) = demorton3(key);
        let (lx, ly, lz) = lane_coords(lane);
        (
            i64::from(tx) * TILE as i64 + lx as i64,
            i64::from(ty) * TILE as i64 + ly as i64,
            i64::from(tz) * TILE as i64 + lz as i64,
        )
    }

    /// Smoothed-fill value with the wall ghost per contact model (2-D
    /// `phi_at` port): `from` is the querying cell.
    fn phi_or_ghost(&self, gx: i64, gy: i64, gz: i64, from: (u64, usize)) -> f64 {
        match self.grid.resolve_source(gx, gy, gz) {
            None => self
                .planar_angle_ghost([gx, gy, gz])
                .unwrap_or_else(|| match self.contact {
                    ContactModel3::Neutral => self.fill(from.0, from.1),
                    ContactModel3::Wetting => 1.0,
                }),
            Some((slot, lane)) => {
                let key = self.grid.active_keys()[slot];
                self.fill_smooth[&key][lane]
            }
        }
    }

    // Geometric equilibrium closure: with n_w pointing from wall into fluid,
    // outward liquid normal n = -grad(phi)/|grad(phi)| must satisfy
    // n.n_w = cos(theta). Thus d(phi)/d(n_w) = -|grad_t(phi)| cot(theta).
    // Extend the smoothed fill into the wall without clipping its derivative.
    // This supplies the normal in the existing sigma*curvature pressure
    // reconstruction; tracked liquid mass is never edited here.
    fn planar_angle_ghost(&self, position: [i64; 3]) -> Option<f64> {
        let contact = self.planar_contact_at(position)?;
        if contact.cotangent == 0.0 {
            return None;
        }
        let axis = contact.wall.axis();
        let dims = self.grid.tile_dims();
        let extents = [dims.0 * TILE, dims.1 * TILE, dims.2 * TILE];
        let edge = if contact.wall.upper() {
            extents[axis] as i64 - 1
        } else {
            0
        };
        let distance = if contact.wall.upper() {
            position[axis] - edge
        } else {
            -position[axis]
        };
        if distance <= 0 {
            return None;
        }
        let mut base = position;
        base[axis] = edge;
        let sample = |p: [i64; 3]| -> Option<f64> {
            let (slot, lane) = self.grid.resolve_source(p[0], p[1], p[2])?;
            Some(self.fill_smooth[&self.grid.active_keys()[slot]][lane])
        };
        let value = sample(base)?;
        let mut tangent_squared = 0.0;
        for tangent in 0..3 {
            if tangent == axis {
                continue;
            }
            let mut plus = base;
            let mut minus = base;
            plus[tangent] += 1;
            minus[tangent] -= 1;
            let derivative = match (sample(plus), sample(minus)) {
                (Some(p), Some(m)) => (p - m) / 2.0,
                (Some(p), None) => p - value,
                (None, Some(m)) => value - m,
                (None, None) => 0.0,
            };
            tangent_squared += derivative * derivative;
        }
        Some(value + distance as f64 * tangent_squared.sqrt() * contact.cotangent)
    }

    /// Refresh the smoothed fill field (average over self + non-wall
    /// neighbors, 2-D port with the 18-neighbor stencil).
    fn refresh_fill(&mut self, gate: &CancelGate) -> Result<(), FreeSurfaceError3> {
        let mut raw: BTreeMap<u64, [f64; TILE_CELLS]> = BTreeMap::new();
        for &key in self.tiles.keys() {
            check_cancel(gate)?;
            raw.insert(key, core::array::from_fn(|lane| self.fill(key, lane)));
        }
        let mut smooth = BTreeMap::new();
        for &key in self.tiles.keys() {
            check_cancel(gate)?;
            let mut vals = [0.0; TILE_CELLS];
            for (lane, out) in vals.iter_mut().enumerate() {
                let mut acc = raw[&key][lane];
                let mut count = 1.0;
                for q in 1..Q3 {
                    if let Neighbor::Active(nb_key, nb_lane, _) = self.neighbor_state(key, lane, q)
                    {
                        acc += raw[&nb_key][nb_lane];
                        count += 1.0;
                    }
                }
                *out = acc / count;
            }
            smooth.insert(key, vals);
        }
        self.fill_smooth = smooth;
        Ok(())
    }

    /// Reference density for gas reconstruction: `1 + σκ/cs²` (2-D port;
    /// σ = 0 short-circuits to exactly 1).
    fn reference_density(&self, key: u64, lane: usize) -> f64 {
        let sigma = self
            .local_sigma
            .as_ref()
            .map_or(self.sigma, |field| field[&(key, lane)]);
        if sigma == 0.0 {
            return 1.0;
        }
        let kappa = self.curvature(key, lane);
        sigma.mul_add(kappa / CS2_3, 1.0)
    }

    /// Sum of principal curvatures from raw VOF height functions where a
    /// complete single-interface stencil exists. Walls and underresolved
    /// geometry retain `div(n̂)`, `n̂ = −∇φ/|∇φ|` of the smoothed fill,
    /// so the contact-angle ghost remains the boundary-condition owner.
    fn curvature(&self, key: u64, lane: usize) -> f64 {
        let (cx, cy, cz) = Self::global_of(key, lane);
        let phi = |dx: i64, dy: i64, dz: i64| -> f64 {
            self.phi_or_ghost(cx + dx, cy + dy, cz + dz, (key, lane))
        };
        // Preserve the original ghost's two-cell support. The height stencil
        // includes its diagonal samples but spans only one transverse cell.
        let wall_free = (0..3).all(|axis| {
            [-2, 2].into_iter().all(|distance| {
                let mut p = [cx, cy, cz];
                p[axis] += distance;
                self.grid.resolve_source(p[0], p[1], p[2]).is_some()
            })
        });
        if wall_free {
            let normal = [
                phi(-1, 0, 0) - phi(1, 0, 0),
                phi(0, -1, 0) - phi(0, 1, 0),
                phi(0, 0, -1) - phi(0, 0, 1),
            ];
            if let Some(kappa) = curvature::height_function(normal, |offset| {
                let (slot, sample_lane) =
                    self.grid
                        .resolve_source(cx + offset[0], cy + offset[1], cz + offset[2])?;
                Some(self.fill(self.grid.active_keys()[slot], sample_lane))
            }) {
                return kappa.clamp(-Self::CURVATURE_LIMIT, Self::CURVATURE_LIMIT);
            }
        }
        let nhat = |dx: i64, dy: i64, dz: i64| -> [f64; 3] {
            let gx = (phi(dx + 1, dy, dz) - phi(dx - 1, dy, dz)) / 2.0;
            let gy = (phi(dx, dy + 1, dz) - phi(dx, dy - 1, dz)) / 2.0;
            let gz = (phi(dx, dy, dz + 1) - phi(dx, dy, dz - 1)) / 2.0;
            let m = (gx * gx + gy * gy + gz * gz).sqrt().max(1e-9);
            [-gx / m, -gy / m, -gz / m]
        };
        let div = (nhat(1, 0, 0)[0] - nhat(-1, 0, 0)[0]) / 2.0
            + (nhat(0, 1, 0)[1] - nhat(0, -1, 0)[1]) / 2.0
            + (nhat(0, 0, 1)[2] - nhat(0, 0, -1)[2]) / 2.0;
        div.clamp(-Self::CURVATURE_LIMIT, Self::CURVATURE_LIMIT)
    }

    /// One free-surface step (2-D `step` port): refresh fill, collide wet
    /// cells, pairwise-antisymmetric mass exchange + pull-streaming with
    /// gas reconstruction, commit, conversion cascade with conservative
    /// redistribution. Tracks the per-step ledger violation.
    ///
    /// # Errors
    /// [`FreeSurfaceError3::Collision`] fail-closed on the first refusing
    /// wet cell in canonical order (state unchanged in that case).
    pub fn step(&mut self) -> Result<(), FreeSurfaceError3> {
        self.step_with_gate(&CancelGate::new_clock_free())
    }

    /// The same step with cancellation polls between active tiles and before
    /// publication. Cancellation preserves populations, tracked mass and all
    /// conversion counters. Once publication begins its conversion cascade
    /// finishes; a later cancellation request is observed by the next step.
    pub fn step_with_gate(&mut self, gate: &CancelGate) -> Result<(), FreeSurfaceError3> {
        self.step_inner(gate, None, None)
    }

    fn step_inner(
        &mut self,
        gate: &CancelGate,
        mut transport: Option<&mut MassTransport3>,
        local_tau: Option<&BTreeMap<(u64, usize), f64>>,
    ) -> Result<(), FreeSurfaceError3> {
        check_cancel(gate)?;
        if let Some(field) = &self.local_sigma {
            self.admit_surface_tensions(field, gate)?;
        }
        let ledger_before = self.ledger_mass();
        self.refresh_fill(gate)?;
        let (model, _) = self.grid.collision();

        // Collide wet lanes into a post map (canonical order). Gas lanes
        // keep their raw populations (never read back, kept defined).
        let keys: Vec<u64> = self.tiles.keys().copied().collect();
        let mut post: BTreeMap<u64, Vec<[f64; Q3]>> = BTreeMap::new();
        for &key in &keys {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active");
            let mut tile_post = vec![[0.0; Q3]; TILE_CELLS];
            for (lane, out) in tile_post.iter_mut().enumerate() {
                let populations = self.grid.populations(slot, lane);
                *out = match self.tiles[&key].cells[lane] {
                    Cell3::Gas => populations,
                    Cell3::Fluid | Cell3::Interface => {
                        let cell_model = match local_tau {
                            None => model,
                            Some(values) => super::CollisionModel3::Bgk {
                                tau: *values.get(&(key, lane)).ok_or(
                                    FreeSurfaceError3::Transport(
                                        "missing source viscosity at a wet cell",
                                    ),
                                )?,
                            },
                        };
                        super::collide_cell3(populations, cell_model, self.force_at_lane(key, lane))
                            .map_err(|source| FreeSurfaceError3::Collision {
                                tile_key: key,
                                lane,
                                source,
                            })?
                    }
                };
            }
            post.insert(key, tile_post);
        }

        // Pre-collision fills (2-D uses fills computed before streaming).
        let fills: BTreeMap<u64, [f64; TILE_CELLS]> = keys
            .iter()
            .map(|&key| {
                let mut vals = [0.0; TILE_CELLS];
                for (lane, v) in vals.iter_mut().enumerate() {
                    *v = self.fill(key, lane);
                }
                (key, vals)
            })
            .collect();

        // Exchange + stream (canonical order), committing per wet lane.
        let mut new_populations: Vec<(u64, usize, [f64; Q3])> = Vec::new();
        let mut mass_delta: Vec<(u64, usize, f64)> = Vec::new();
        for &key in &keys {
            check_cancel(gate)?;
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                let flag = self.tiles[&key].cells[lane];
                if flag == Cell3::Gas {
                    continue;
                }
                // The optional companion records these SAME pair exchanges.
                // Fluid Σf tracks its mass through streaming, not another update.
                if flag == Cell3::Interface || transport.is_some() {
                    let mut dm = 0.0f64;
                    for q in 1..Q3 {
                        if let Neighbor::Active(nb_key, nb_lane, nb_cell) =
                            self.neighbor_state(key, lane, q)
                        {
                            let w = match (flag, nb_cell) {
                                (_, Cell3::Gas) => 0.0,
                                (Cell3::Interface, Cell3::Interface) => {
                                    f64::midpoint(fills[&key][lane], fills[&nb_key][nb_lane])
                                }
                                _ => 1.0,
                            };
                            if w > 0.0 {
                                let incoming =
                                    w * (post[&nb_key][nb_lane][OPP3[q]] - post[&key][lane][q]);
                                dm += incoming;
                                if (key, lane) < (nb_key, nb_lane)
                                    && let Some(trace) = transport.as_deref_mut()
                                {
                                    trace.exchange(
                                        MassLocation3::Cell(nb_key, nb_lane),
                                        MassLocation3::Cell(key, lane),
                                        incoming,
                                    )?;
                                }
                            }
                        }
                    }
                    if flag == Cell3::Interface {
                        mass_delta.push((key, lane, dm));
                    }
                }
                // Pull-stream with gas reconstruction (2-D port).
                let populations = self.grid.populations(slot, lane);
                let rho_pre: f64 = populations.iter().sum();
                let mut u_pre = [
                    populations
                        .iter()
                        .enumerate()
                        .map(|(q, f)| f * f64::from(E3[q].0))
                        .sum::<f64>()
                        / rho_pre,
                    populations
                        .iter()
                        .enumerate()
                        .map(|(q, f)| f * f64::from(E3[q].1))
                        .sum::<f64>()
                        / rho_pre,
                    populations
                        .iter()
                        .enumerate()
                        .map(|(q, f)| f * f64::from(E3[q].2))
                        .sum::<f64>()
                        / rho_pre,
                ];
                let local_force = self.force_at_lane(key, lane);
                let base_force = self.grid.collision().1;
                let perturbation = core::array::from_fn(|d| local_force[d] - base_force[d]);
                let physical_velocity =
                    core::array::from_fn(|d| u_pre[d] + 0.5 * local_force[d] / rho_pre);
                if self.local_force.is_some() {
                    for axis in 0..3 {
                        u_pre[axis] += 0.5 * perturbation[axis] / rho_pre;
                    }
                }
                let mut new_f = [0.0f64; Q3];
                // One geometry and boundary velocity serve every missing gas
                // population. Interior cells never compute capillary curvature.
                let mut gas_equilibrium = None;
                let (cx, cy, cz) = Self::global_of(key, lane);
                for (q, out) in new_f.iter_mut().enumerate() {
                    // Pull source lies opposite the direction of travel.
                    let sx = cx - i64::from(E3[q].0);
                    let sy = cy - i64::from(E3[q].1);
                    let sz = cz - i64::from(E3[q].2);
                    *out = match self.grid.resolve_source(sx, sy, sz) {
                        None => post[&key][lane][OPP3[q]],
                        Some((src_slot, src_lane)) => {
                            let src_key = self.grid.active_keys()[src_slot];
                            match self.tiles[&src_key].cells[src_lane] {
                                Cell3::Gas => {
                                    let eq = gas_equilibrium.get_or_insert_with(|| {
                                        equilibrium3(self.reference_density(key, lane), u_pre)
                                    });
                                    let value = eq[q] + eq[OPP3[q]] - post[&key][lane][OPP3[q]];
                                    if self.local_force.is_some() {
                                        value
                                            - 0.5
                                                * (forcing::source(
                                                    q,
                                                    physical_velocity,
                                                    perturbation,
                                                ) + forcing::source(
                                                    OPP3[q],
                                                    physical_velocity,
                                                    perturbation,
                                                ))
                                    } else {
                                        value
                                    }
                                }
                                Cell3::Fluid | Cell3::Interface => post[&src_key][src_lane][q],
                            }
                        }
                    };
                }
                new_populations.push((key, lane, new_f));
            }
        }
        check_cancel(gate)?;
        for (key, lane, f) in new_populations {
            let slot = self.grid.slot_of(key).expect("active");
            self.grid.set_populations(slot, lane, f);
        }
        for (key, lane, dm) in mass_delta {
            self.tiles.get_mut(&key).expect("active").mass[lane] += dm;
        }

        self.apply_conversions(transport)?;

        let ledger_after = self.ledger_mass();
        let violation = ((ledger_after - ledger_before) / ledger_before.max(1e-12)).abs();
        if violation > self.worst_step_violation {
            self.worst_step_violation = violation;
        }
        Ok(())
    }

    /// Conversion cascade with conservative neighbor-local redistribution,
    /// activating inactive in-domain tiles when fluid advances into them.
    fn apply_conversions(
        &mut self,
        mut transport: Option<&mut MassTransport3>,
    ) -> Result<(), FreeSurfaceError3> {
        let mut excess = Vec::new();
        let keys: Vec<u64> = self.tiles.keys().copied().collect();
        let mut to_fluid: Vec<(u64, usize)> = Vec::new();
        let mut to_gas: Vec<(u64, usize)> = Vec::new();
        for &key in &keys {
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                if self.tiles[&key].cells[lane] != Cell3::Interface {
                    continue;
                }
                let rho: f64 = self.grid.populations(slot, lane).iter().sum();
                let mass = self.tiles[&key].mass[lane];
                if mass > (1.0 + EPS) * rho {
                    to_fluid.push((key, lane));
                } else if mass < -EPS * rho {
                    to_gas.push((key, lane));
                }
            }
        }

        // Interface → fluid; gas neighbors become interface (activating
        // their tile first when fluid advances into inactive space).
        for &(key, lane) in &to_fluid {
            let slot = self.grid.slot_of(key).expect("active");
            let rho: f64 = self.grid.populations(slot, lane).iter().sum();
            {
                let tile = self.tiles.get_mut(&key).expect("active");
                excess.push((key, lane, tile.mass[lane] - rho));
                tile.cells[lane] = Cell3::Fluid;
                tile.mass[lane] = rho;
            }
            self.conversions.to_fluid += 1;
            for q in 1..Q3 {
                let (cx, cy, cz) = Self::global_of(key, lane);
                let (gx, gy, gz) = (
                    cx + i64::from(E3[q].0),
                    cy + i64::from(E3[q].1),
                    cz + i64::from(E3[q].2),
                );
                self.ensure_active_for_advance(gx, gy, gz)?;
                let Some((nb_slot, nb_lane)) = self.grid.resolve_source(gx, gy, gz) else {
                    continue; // out of domain: wall
                };
                let nb_key = self.grid.active_keys()[nb_slot];
                if self.tiles[&nb_key].cells[nb_lane] != Cell3::Gas {
                    continue;
                }
                // Initialize from the average of wet neighbors (2-D port).
                let mut rho_avg = 0.0;
                let mut u_avg = [0.0f64; 3];
                let mut force_avg = [0.0f64; 3];
                let mut cnt = 0.0;
                for q2 in 1..Q3 {
                    if let Neighbor::Active(nn_key, nn_lane, nn_cell) =
                        self.neighbor_state(nb_key, nb_lane, q2)
                        && matches!(nn_cell, Cell3::Fluid | Cell3::Interface)
                    {
                        let nn_slot = self.grid.slot_of(nn_key).expect("active");
                        let f = self.grid.populations(nn_slot, nn_lane);
                        let r: f64 = f.iter().sum();
                        rho_avg += r;
                        for (q3, fq) in f.iter().enumerate() {
                            u_avg[0] += fq * f64::from(E3[q3].0) / r;
                            u_avg[1] += fq * f64::from(E3[q3].1) / r;
                            u_avg[2] += fq * f64::from(E3[q3].2) / r;
                        }
                        if self.local_force.is_some() {
                            let force = self.force_at_lane(nn_key, nn_lane);
                            let base = self.grid.collision().1;
                            for axis in 0..3 {
                                u_avg[axis] += 0.5 * (force[axis] - base[axis]) / r;
                                force_avg[axis] += force[axis];
                            }
                        }
                        cnt += 1.0;
                    }
                }
                if cnt > 0.0 {
                    rho_avg /= cnt;
                    u_avg = [u_avg[0] / cnt, u_avg[1] / cnt, u_avg[2] / cnt];
                } else {
                    rho_avg = 1.0;
                }
                let base_force = self.grid.collision().1;
                let mut initial = equilibrium3(rho_avg, u_avg);
                if let Some(forces) = &mut self.local_force {
                    if cnt > 0.0 {
                        for axis in 0..3 {
                            force_avg[axis] /= cnt;
                        }
                    } else {
                        return Err(FreeSurfaceError3::Transport(
                            "new forced interface needs a wet source neighbor",
                        ));
                    }
                    let physical_velocity =
                        core::array::from_fn(|d| u_avg[d] + 0.5 * base_force[d] / rho_avg);
                    initial = forcing::reframe(initial, physical_velocity, base_force, force_avg)?;
                    forces.insert((nb_key, nb_lane), force_avg);
                }
                let nb_slot = self.grid.slot_of(nb_key).expect("active");
                self.grid.set_populations(nb_slot, nb_lane, initial);
                let nb_tile = self.tiles.get_mut(&nb_key).expect("active");
                nb_tile.cells[nb_lane] = Cell3::Interface;
                nb_tile.mass[nb_lane] = 0.0;
                self.conversions.gas_to_interface += 1;
                // INVARIANT: every interface cell has a fully active
                // neighborhood. Without this, an interface reaching the
                // activation frontier is backed by artificial wall,
                // behaves like fluid-against-wall, never crosses the
                // to-fluid threshold, and the front pins at the frontier
                // forever (observed: extent stuck at the margin edge).
                self.ensure_interface_neighborhood(nb_key, nb_lane)?;
            }
        }

        // Interface → gas; fluid neighbors become interface (Σf IS their
        // mass — ledger unchanged), 2-D port.
        for &(key, lane) in &to_gas {
            if self.tiles[&key].cells[lane] != Cell3::Interface {
                continue; // re-flagged by the cascade
            }
            {
                let tile = self.tiles.get_mut(&key).expect("active");
                excess.push((key, lane, tile.mass[lane]));
                tile.cells[lane] = Cell3::Gas;
                tile.mass[lane] = 0.0;
            }
            self.conversions.to_gas += 1;
            if let Some(forces) = &mut self.local_force {
                forces.remove(&(key, lane));
            }
            for q in 1..Q3 {
                if let Neighbor::Active(nb_key, nb_lane, Cell3::Fluid) =
                    self.neighbor_state(key, lane, q)
                {
                    let nb_slot = self.grid.slot_of(nb_key).expect("active");
                    let rho: f64 = self.grid.populations(nb_slot, nb_lane).iter().sum();
                    let nb_tile = self.tiles.get_mut(&nb_key).expect("active");
                    nb_tile.cells[nb_lane] = Cell3::Interface;
                    nb_tile.mass[nb_lane] = rho;
                    self.conversions.fluid_to_interface += 1;
                    self.ensure_interface_neighborhood(nb_key, nb_lane)?;
                }
            }
        }

        if self.redistribution == ConversionRedistribution3::LegacyGlobal {
            // Keep the historical operation tree only for explicit version-1
            // replay. Its global pooling is intentionally absent by default.
            let mut pool = std::mem::take(&mut self.carry);
            for (key, lane, amount) in excess {
                if let Some(trace) = transport.as_deref_mut() {
                    trace.exchange(MassLocation3::Cell(key, lane), MassLocation3::Carry, amount)?;
                }
                pool += amount;
            }
            let count = self
                .tiles
                .values()
                .map(|tile| {
                    tile.cells
                        .iter()
                        .filter(|&&cell| cell == Cell3::Interface)
                        .count()
                })
                .sum::<usize>();
            if count == 0 {
                self.carry = pool;
            } else {
                let share = pool / count as f64;
                for (&key, tile) in &mut self.tiles {
                    for lane in 0..TILE_CELLS {
                        if tile.cells[lane] == Cell3::Interface {
                            if let Some(trace) = transport.as_deref_mut() {
                                trace.exchange(
                                    MassLocation3::Carry,
                                    MassLocation3::Cell(key, lane),
                                    share,
                                )?;
                            }
                            tile.mass[lane] += share;
                        }
                    }
                }
            }
            return Ok(());
        }

        // Finish the complete flag cascade before finding local recipients.
        // Direct routes preserve each drop's mass AND extensive companions;
        // a shared zero-storage Carry node would spuriously mix their heat.
        // Fixed-size neighbor storage avoids a whole-domain redistribution pass.
        for (key, lane, amount) in excess {
            if amount == 0.0 {
                continue;
            }
            let mut recipients = [(0, 0); Q3 - 1];
            let mut count = 0;
            for q in 1..Q3 {
                if let Neighbor::Active(nb_key, nb_lane, Cell3::Interface) =
                    self.neighbor_state(key, lane, q)
                {
                    recipients[count] = (nb_key, nb_lane);
                    count += 1;
                }
            }
            if count == 0 {
                // A scalar carry has no position. Keep orphan excess/debt
                // explicitly unresolved rather than injecting it elsewhere.
                if let Some(trace) = transport.as_deref_mut() {
                    trace.exchange(MassLocation3::Cell(key, lane), MassLocation3::Carry, amount)?;
                }
                self.carry += amount;
                continue;
            }
            let share = amount / count as f64;
            let mut remaining = amount;
            for (j, &(nb_key, nb_lane)) in recipients[..count].iter().enumerate() {
                let transfer = if j + 1 == count { remaining } else { share };
                if let Some(trace) = transport.as_deref_mut() {
                    trace.exchange(
                        MassLocation3::Cell(key, lane),
                        MassLocation3::Cell(nb_key, nb_lane),
                        transfer,
                    )?;
                }
                self.tiles.get_mut(&nb_key).expect("active recipient").mass[nb_lane] += transfer;
                remaining -= transfer;
            }
        }
        Ok(())
    }

    /// Activate every inactive in-domain tile touching the 18-neighborhood
    /// of a newly created interface cell, upholding the invariant that an
    /// interface cell never sees artificial wall where gas should be.
    fn ensure_interface_neighborhood(
        &mut self,
        key: u64,
        lane: usize,
    ) -> Result<(), FreeSurfaceError3> {
        let (cx, cy, cz) = Self::global_of(key, lane);
        for q in 1..Q3 {
            self.ensure_active_for_advance(
                cx + i64::from(E3[q].0),
                cy + i64::from(E3[q].1),
                cz + i64::from(E3[q].2),
            )?;
        }
        Ok(())
    }

    /// Activate the tile containing an in-domain global cell if inactive —
    /// fluid is advancing into it (equilibrium-initialized, all lanes gas,
    /// zero tracked mass). Out-of-domain coordinates are ignored (wall).
    fn ensure_active_for_advance(
        &mut self,
        gx: i64,
        gy: i64,
        gz: i64,
    ) -> Result<(), FreeSurfaceError3> {
        if gx < 0 || gy < 0 || gz < 0 {
            return Ok(());
        }
        let (ntx, nty, ntz) = self.grid.tile_dims();
        let (tx, ty, tz) = (gx as usize / TILE, gy as usize / TILE, gz as usize / TILE);
        if tx >= ntx || ty >= nty || tz >= ntz {
            return Ok(());
        }
        #[allow(clippy::cast_possible_truncation)]
        let (tx, ty, tz) = (tx as u32, ty as u32, tz as u32);
        let key = morton3(tx, ty, tz);
        if self.tiles.contains_key(&key) || self.solid.contains(&key) {
            return Ok(());
        }
        self.grid.activate_tiles(&[(tx, ty, tz)])?;
        self.tiles.insert(key, FsTile::all_gas());
        self.conversions.tiles_activated += 1;
        Ok(())
    }

    /// Retire tiles whose lanes are ALL gas (memory maintenance; the
    /// free-surface ledger is unchanged because gas lanes carry no ledger
    /// mass — the sparse-layer mass returned by retirement is stale gas
    /// population data, deliberately dropped). Returns retired tile count.
    ///
    /// # Errors
    /// [`FreeSurfaceError3::Grid`] if retirement refuses (not expected for
    /// keys taken from the active set).
    pub fn retire_gas_tiles(&mut self) -> Result<usize, FreeSurfaceError3> {
        let mut retire: Vec<(u32, u32, u32)> = Vec::new();
        let mut retire_keys: Vec<u64> = Vec::new();
        for (&key, tile) in &self.tiles {
            if tile.cells.iter().all(|c| *c == Cell3::Gas) {
                let (tx, ty, tz) = demorton3(key);
                retire.push((tx, ty, tz));
                retire_keys.push(key);
            }
        }
        if retire.is_empty() {
            return Ok(0);
        }
        self.grid.deactivate_tiles(&retire)?;
        for key in retire_keys {
            self.tiles.remove(&key);
            self.fill_smooth.remove(&key);
        }
        Ok(retire.len())
    }
}

fn check_cancel(gate: &CancelGate) -> Result<(), FreeSurfaceError3> {
    if gate.is_requested() {
        Err(FreeSurfaceError3::Cancelled)
    } else {
        Ok(())
    }
}

/// Neighbor resolution outcome.
enum Neighbor {
    /// Out-of-domain or inactive-tile space: wall by sparse convention.
    Wall,
    /// A cell in an active tile.
    Active(u64, usize, Cell3),
}

#[cfg(test)]
mod capillary_boundary_tests {
    use super::*;

    fn separated_drops() -> FreeSurface3 {
        FreeSurface3::new(
            24,
            12,
            12,
            0.8,
            [0.0; 3],
            0.002,
            ContactModel3::Neutral,
            |x, y, z| {
                (2..7).contains(&y)
                    && (2..7).contains(&z)
                    && ((2..7).contains(&x) || (17..22).contains(&x))
            },
        )
        .unwrap()
    }

    fn assert_route_balance(trace: &MassTransport3) {
        let mut balance = trace.before.clone();
        for route in &trace.routes {
            *balance.entry(route.from).or_default() -= route.mass;
            *balance.entry(route.to).or_default() += route.mass;
        }
        for (&key, &mass) in &trace.after {
            *balance.entry(key).or_default() -= mass;
        }
        for (key, residual) in balance {
            assert!(residual.abs() < 2e-14, "{key:?}: {residual}");
        }
    }

    #[test]
    fn g0_conversion_excess_and_debt_use_local_mass_routes() {
        for amount in [0.25, -0.25] {
            let mut sim = separated_drops();
            let (slot, lane) = sim.grid.resolve_source(6, 4, 4).unwrap();
            let key = sim.grid.active_keys()[slot];
            assert_eq!(sim.tiles[&key].cells[lane], Cell3::Interface);
            let rho = sim.grid.populations(slot, lane).iter().sum::<f64>();
            sim.tiles.get_mut(&key).unwrap().mass[lane] =
                if amount > 0.0 { rho + amount } else { amount };
            let before = sim.mass_snapshot();
            let mass = sim.ledger_mass();
            let distant: Vec<_> = before
                .keys()
                .filter_map(|&location| {
                    let p = FreeSurface3::transport_coordinates(location)?;
                    (p[0] >= 14).then_some((location, before[&location].to_bits()))
                })
                .collect();
            let mut trace = MassTransport3::new(before, 1000);
            sim.apply_conversions(Some(&mut trace)).unwrap();
            trace.after = sim.mass_snapshot();
            assert_route_balance(&trace);
            assert!((sim.ledger_mass() - mass).abs() < 1e-11);
            assert_eq!(sim.carry, 0.0);
            assert!(!trace.routes.is_empty());
            for route in &trace.routes {
                let from = FreeSurface3::transport_coordinates(route.from).unwrap();
                let to = FreeSurface3::transport_coordinates(route.to).unwrap();
                assert!(from[0] < 14 && to[0] < 14);
                assert!((0..3).all(|d| from[d].abs_diff(to[d]) <= 1));
                assert_eq!(
                    if amount > 0.0 { route.from } else { route.to },
                    MassLocation3::Cell(key, lane)
                );
            }
            for (location, bits) in distant {
                assert_eq!(trace.after[&location].to_bits(), bits);
            }
        }
    }

    #[test]
    fn g0_orphan_conversion_debt_stays_in_carry_without_remote_injection() {
        let mut sim = separated_drops();
        let (slot, lane) = sim.grid.resolve_source(11, 4, 4).unwrap();
        let key = sim.grid.active_keys()[slot];
        sim.tiles.get_mut(&key).unwrap().cells[lane] = Cell3::Interface;
        sim.tiles.get_mut(&key).unwrap().mass[lane] = -0.25;
        let mut trace = MassTransport3::new(sim.mass_snapshot(), 1000);
        sim.apply_conversions(Some(&mut trace)).unwrap();
        trace.after = sim.mass_snapshot();
        assert_route_balance(&trace);
        assert_eq!(trace.routes.len(), 1);
        assert_eq!(trace.routes[0].from, MassLocation3::Carry);
        assert_eq!(trace.routes[0].to, MassLocation3::Cell(key, lane));
        assert_eq!(sim.carry, -0.25);
        let retained = trace.after;
        sim.apply_conversions(None).unwrap();
        assert_eq!(sim.mass_snapshot(), retained);
    }

    #[test]
    fn g0_scalar_capillarity_refuses_nonphysical_boundary_densities() {
        for sigma in [-0.01, f64::NAN, f64::INFINITY, CS2_3, 1.0] {
            assert!(matches!(
                FreeSurface3::new(
                    8,
                    8,
                    8,
                    0.8,
                    [0.0; 3],
                    sigma,
                    ContactModel3::Neutral,
                    |_, _, _| panic!("invalid sigma must refuse before geometry"),
                ),
                Err(FreeSurfaceError3::InvalidInitialState { .. })
            ));
        }
    }

    #[test]
    fn g1_resolved_vof_cap_drives_native_laplace_pressure_and_inward_motion() {
        let radius = 16.0_f64;
        let sigma = 0.001;
        let mut surface = FreeSurface3::new(
            16,
            16,
            16,
            0.8,
            [0.0; 3],
            sigma,
            ContactModel3::Neutral,
            |_, _, _| true,
        )
        .unwrap();
        // Independent cell-volume quadrature of a resolved spherical cap.
        // The populations begin at zero pressure jump, so the first response
        // should accelerate inward; this is not a static-equilibrium claim.
        let node = 0.5 / 3.0_f64.sqrt();
        for (&key, tile) in &mut surface.tiles {
            for lane in 0..TILE_CELLS {
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                let mut fill = 0.0;
                for du in [-node, node] {
                    for dv in [-node, node] {
                        let u = x as f64 - 8.0 + du;
                        let v = y as f64 - 8.0 + dv;
                        let height = 7.2 + (radius * radius - u * u - v * v).sqrt() - radius;
                        fill += 0.25 * (height - z as f64 + 0.5).clamp(0.0, 1.0);
                    }
                }
                tile.cells[lane] = if fill == 0.0 {
                    Cell3::Gas
                } else if fill == 1.0 {
                    Cell3::Fluid
                } else {
                    Cell3::Interface
                };
                tile.mass[lane] = fill;
            }
        }
        // Full cells bordering gas still need the original interface closure.
        let mut promote = Vec::new();
        for (&key, tile) in &surface.tiles {
            for lane in 0..TILE_CELLS {
                if tile.cells[lane] == Cell3::Fluid
                    && (1..Q3).any(|q| {
                        matches!(
                            surface.neighbor_state(key, lane, q),
                            Neighbor::Active(_, _, Cell3::Gas)
                        )
                    })
                {
                    promote.push((key, lane));
                }
            }
        }
        for (key, lane) in promote {
            surface.tiles.get_mut(&key).unwrap().cells[lane] = Cell3::Interface;
        }
        surface.refresh_fill(&CancelGate::new_clock_free()).unwrap();
        let (slot, lane) = surface.grid.resolve_source(8, 8, 7).unwrap();
        let key = surface.grid.active_keys()[slot];
        let pressure = (surface.reference_density(key, lane) - 1.0) * CS2_3;
        assert!((pressure / (2.0 * sigma / radius) - 1.0).abs() < 0.01);
        let mass = surface.ledger_mass();
        surface.step().unwrap();
        assert!((surface.ledger_mass() - mass).abs() < 1e-10);
        assert!(surface.sample_at(8, 8, 7).unwrap().velocity[2] < 0.0);
    }

    #[test]
    fn g1_planar_wall_ghost_has_the_prescribed_angle_and_zero_planar_curvature() {
        for wall in [
            PlanarWall3::XMin,
            PlanarWall3::XMax,
            PlanarWall3::YMin,
            PlanarWall3::YMax,
            PlanarWall3::ZMin,
            PlanarWall3::ZMax,
        ] {
            for angle in [
                core::f64::consts::PI / 3.0,
                2.0 * core::f64::consts::PI / 3.0,
            ] {
                let mut surface = FreeSurface3::new(
                    8,
                    8,
                    8,
                    0.8,
                    [0.0; 3],
                    0.001,
                    ContactModel3::Neutral,
                    |_, _, _| true,
                )
                .unwrap();
                surface.set_planar_contact_angle(wall, angle).unwrap();
                let axis = wall.axis();
                let tangent = (axis + 1) % 3;
                // Independent analytical plane: its normal makes angle theta
                // with the inward wall normal. No liquid-motion solver is
                // involved in assigning this manufactured smoothed fill.
                for &key in surface.grid.active_keys() {
                    let values = core::array::from_fn(|lane| {
                        let (x, y, z) = FreeSurface3::global_of(key, lane);
                        let p = [x, y, z];
                        let inward = if wall.upper() { 7 - p[axis] } else { p[axis] };
                        0.5 - 0.03 * p[tangent] as f64
                            - 0.03 * (angle.cos() / angle.sin()) * inward as f64
                    });
                    surface.fill_smooth.insert(key, values);
                }
                let mut p = [3i64; 3];
                p[axis] = if wall.upper() { 7 } else { 0 };
                let (slot, lane) = surface.grid.resolve_source(p[0], p[1], p[2]).unwrap();
                let key = surface.grid.active_keys()[slot];
                let mut gradient = [0.0; 3];
                for component in 0..3 {
                    let mut plus = p;
                    let mut minus = p;
                    plus[component] += 1;
                    minus[component] -= 1;
                    gradient[component] =
                        (surface.phi_or_ghost(plus[0], plus[1], plus[2], (key, lane))
                            - surface.phi_or_ghost(minus[0], minus[1], minus[2], (key, lane)))
                            / 2.0;
                }
                let magnitude = gradient.iter().map(|x| x * x).sum::<f64>().sqrt();
                let inward_sign = if wall.upper() { -1.0 } else { 1.0 };
                let normal_dot_wall = -gradient[axis] * inward_sign / magnitude;
                assert!(
                    (normal_dot_wall - angle.cos()).abs() < 2e-12,
                    "{wall:?}, angle {angle}: normal dot wall {normal_dot_wall}"
                );
                assert!(surface.curvature(key, lane).abs() < 2e-12);
            }
        }
    }

    #[test]
    fn g1_spherical_laplace_pressure_refines_with_the_existing_curvature_operator() {
        let relative_error = |radius: usize| {
            let cells = 4 * radius;
            let center = (2 * radius) as f64;
            let mut surface = FreeSurface3::new(
                cells,
                cells,
                cells,
                0.8,
                [0.0; 3],
                0.001,
                ContactModel3::Neutral,
                |_, _, _| true,
            )
            .unwrap();
            for &key in surface.grid.active_keys() {
                let values = core::array::from_fn(|lane| {
                    let (x, y, z) = FreeSurface3::global_of(key, lane);
                    let d = [x as f64 - center, y as f64 - center, z as f64 - center];
                    let r = d.iter().map(|x| x * x).sum::<f64>().sqrt();
                    0.5 + 0.01 * (radius as f64 - r)
                });
                surface.fill_smooth.insert(key, values);
            }
            let (slot, lane) = surface
                .grid
                .resolve_source(
                    (3 * radius) as i64,
                    (2 * radius) as i64,
                    (2 * radius) as i64,
                )
                .unwrap();
            let key = surface.grid.active_keys()[slot];
            let curvature = surface.curvature(key, lane);
            let density = surface.reference_density(key, lane);
            assert!(((density - 1.0) * CS2_3 - 0.001 * curvature).abs() < 1e-16);
            (curvature * radius as f64 / 2.0 - 1.0).abs()
        };
        let coarse = relative_error(4);
        let fine = relative_error(8);
        assert!(coarse < 0.06);
        assert!(
            fine < 0.016 && fine < coarse / 3.0,
            "sphere curvature errors coarse={coarse}, fine={fine}"
        );
    }
}
