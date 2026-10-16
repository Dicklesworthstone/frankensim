//! Free-surface 3-D on sparse active tiles (bead sxnm, plan §15.3): the
//! pour primitive's physics layer.
//!
//! This is a FAITHFUL PORT of the 2-D Körner-style VOF machinery in
//! `crate::freesurface` onto the D3Q19 sparse active-tile grid — the same
//! cell taxonomy (fluid / interface / gas), the same pairwise-antisymmetric
//! mass exchange, the same gas-side population reconstruction, the same
//! conversion cascade with a conservation carry, and the same
//! surface-tension reference-density port (smoothed-fill curvature with the
//! contact-model wall ghost). No new physics is introduced here by design
//! (the bead's own rule); every formula is the 2-D crate's, generalized to
//! three dimensions and Morton/lane canonical order.
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
//! turbulence physics, not contact-angle physics beyond the 2-D crate's
//! neutral/wetting ghost, and no interface-reconstruction geometry beyond
//! the smoothed-fill normal the 2-D crate uses.

use std::collections::BTreeMap;

use super::sparse::{SparseError3, SparseGrid3, demorton3, morton3};
use super::{CollisionError3, E3, OPP3, Q3, TILE, equilibrium3};

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

/// 3-D free surface over a sparse active-tile D3Q19 grid.
pub struct FreeSurface3 {
    grid: SparseGrid3,
    tiles: BTreeMap<u64, FsTile>,
    /// Permanently solid tiles (tank walls, fixtures): advance-driven
    /// activation refuses to open these — they are wall FOREVER, unlike
    /// plain inactive tiles which are merely not-yet-reached gas.
    solid: std::collections::BTreeSet<u64>,
    /// Surface tension coefficient (0 = off), 2-D port.
    sigma: f64,
    contact: ContactModel3,
    /// Conversion-conservation carry (mass awaiting redistribution).
    carry: f64,
    conversions: ConversionStats3,
    /// Worst single-step relative ledger violation observed so far —
    /// the bead demands this be logged, not assumed.
    worst_step_violation: f64,
    /// Smoothed fill field, refreshed at the top of each step.
    fill_smooth: BTreeMap<u64, [f64; TILE_CELLS]>,
}

impl core::fmt::Debug for FreeSurface3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("FreeSurface3")
            .field("active_tiles", &self.tiles.len())
            .field("sigma", &self.sigma)
            .field("contact", &self.contact)
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
    /// As [`FreeSurface3::new`]; additionally refuses via
    /// [`FreeSurfaceError3::Grid`] if a solid tile lies out of domain.
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
            contact,
            carry: 0.0,
            conversions: ConversionStats3::default(),
            worst_step_violation: 0.0,
            fill_smooth: BTreeMap::new(),
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

    /// The conversion-conservation carry (global, unregioned).
    #[must_use]
    pub fn carry(&self) -> f64 {
        self.carry
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
            None => match self.contact {
                ContactModel3::Neutral => self.fill(from.0, from.1),
                ContactModel3::Wetting => 1.0,
            },
            Some((slot, lane)) => {
                let key = self.grid.active_keys()[slot];
                self.fill_smooth[&key][lane]
            }
        }
    }

    /// Refresh the smoothed fill field (average over self + non-wall
    /// neighbors, 2-D port with the 18-neighbor stencil).
    fn refresh_fill(&mut self) {
        let raw: BTreeMap<u64, [f64; TILE_CELLS]> = self
            .tiles
            .keys()
            .map(|&key| {
                let mut vals = [0.0; TILE_CELLS];
                for (lane, v) in vals.iter_mut().enumerate() {
                    *v = self.fill(key, lane);
                }
                (key, vals)
            })
            .collect();
        let mut smooth = BTreeMap::new();
        for &key in self.tiles.keys() {
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
    }

    /// Reference density for gas reconstruction: `1 + σκ/cs²` (2-D port;
    /// σ = 0 short-circuits to exactly 1).
    fn reference_density(&self, key: u64, lane: usize) -> f64 {
        if self.sigma == 0.0 {
            return 1.0;
        }
        let kappa = self.curvature(key, lane);
        self.sigma.mul_add(kappa / CS2_3, 1.0)
    }

    /// Curvature of the smoothed fill at the cell: `div(n̂)` with
    /// `n̂ = −∇φ/|∇φ|`, central differences along the three axes (the
    /// exact 3-D analog of the 2-D construction).
    fn curvature(&self, key: u64, lane: usize) -> f64 {
        let (cx, cy, cz) = Self::global_of(key, lane);
        let phi = |dx: i64, dy: i64, dz: i64| -> f64 {
            self.phi_or_ghost(cx + dx, cy + dy, cz + dz, (key, lane))
        };
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
        div.clamp(-1.0, 1.0)
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
        let ledger_before = self.ledger_mass();
        self.refresh_fill();
        let (model, force) = self.grid.collision();

        // Collide wet lanes into a post map (canonical order). Gas lanes
        // keep their raw populations (never read back, kept defined).
        let keys: Vec<u64> = self.tiles.keys().copied().collect();
        let mut post: BTreeMap<u64, Vec<[f64; Q3]>> = BTreeMap::new();
        for &key in &keys {
            let slot = self.grid.slot_of(key).expect("active");
            let mut tile_post = vec![[0.0; Q3]; TILE_CELLS];
            for (lane, out) in tile_post.iter_mut().enumerate() {
                let populations = self.grid.populations(slot, lane);
                *out = match self.tiles[&key].cells[lane] {
                    Cell3::Gas => populations,
                    Cell3::Fluid | Cell3::Interface => {
                        super::collide_cell3(populations, model, force).map_err(|source| {
                            FreeSurfaceError3::Collision {
                                tile_key: key,
                                lane,
                                source,
                            }
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
            let slot = self.grid.slot_of(key).expect("active");
            for lane in 0..TILE_CELLS {
                let flag = self.tiles[&key].cells[lane];
                if flag == Cell3::Gas {
                    continue;
                }
                // Mass exchange (interface only; fluid Σf tracks itself).
                if flag == Cell3::Interface {
                    let mut dm = 0.0f64;
                    for q in 1..Q3 {
                        if let Neighbor::Active(nb_key, nb_lane, nb_cell) =
                            self.neighbor_state(key, lane, q)
                        {
                            let w = match nb_cell {
                                Cell3::Fluid => 1.0,
                                Cell3::Interface => {
                                    f64::midpoint(fills[&key][lane], fills[&nb_key][nb_lane])
                                }
                                Cell3::Gas => 0.0,
                            };
                            if w > 0.0 {
                                dm += w * (post[&nb_key][nb_lane][OPP3[q]] - post[&key][lane][q]);
                            }
                        }
                    }
                    mass_delta.push((key, lane, dm));
                }
                // Pull-stream with gas reconstruction (2-D port).
                let populations = self.grid.populations(slot, lane);
                let rho_pre: f64 = populations.iter().sum();
                let u_pre = [
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
                let mut new_f = [0.0f64; Q3];
                for (q, out) in new_f.iter_mut().enumerate() {
                    // Pull source lies opposite the direction of travel.
                    let (cx, cy, cz) = Self::global_of(key, lane);
                    let sx = cx - i64::from(E3[q].0);
                    let sy = cy - i64::from(E3[q].1);
                    let sz = cz - i64::from(E3[q].2);
                    *out = match self.grid.resolve_source(sx, sy, sz) {
                        None => post[&key][lane][OPP3[q]],
                        Some((src_slot, src_lane)) => {
                            let src_key = self.grid.active_keys()[src_slot];
                            match self.tiles[&src_key].cells[src_lane] {
                                Cell3::Gas => {
                                    let rho_ref = self.reference_density(key, lane);
                                    let eq = equilibrium3(rho_ref, u_pre);
                                    eq[q] + eq[OPP3[q]] - post[&key][lane][OPP3[q]]
                                }
                                Cell3::Fluid | Cell3::Interface => post[&src_key][src_lane][q],
                            }
                        }
                    };
                }
                new_populations.push((key, lane, new_f));
            }
        }
        for (key, lane, f) in new_populations {
            let slot = self.grid.slot_of(key).expect("active");
            self.grid.set_populations(slot, lane, f);
        }
        for (key, lane, dm) in mass_delta {
            self.tiles.get_mut(&key).expect("active").mass[lane] += dm;
        }

        self.apply_conversions()?;

        let ledger_after = self.ledger_mass();
        let violation = ((ledger_after - ledger_before) / ledger_before.max(1e-12)).abs();
        if violation > self.worst_step_violation {
            self.worst_step_violation = violation;
        }
        Ok(())
    }

    /// Conversion cascade with conservative redistribution (2-D port),
    /// activating inactive in-domain tiles when fluid advances into them.
    fn apply_conversions(&mut self) -> Result<(), FreeSurfaceError3> {
        let mut excess_pool = std::mem::take(&mut self.carry);
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
                excess_pool += tile.mass[lane] - rho;
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
                        cnt += 1.0;
                    }
                }
                if cnt > 0.0 {
                    rho_avg /= cnt;
                    u_avg = [u_avg[0] / cnt, u_avg[1] / cnt, u_avg[2] / cnt];
                } else {
                    rho_avg = 1.0;
                }
                let nb_slot = self.grid.slot_of(nb_key).expect("active");
                self.grid
                    .set_populations(nb_slot, nb_lane, equilibrium3(rho_avg, u_avg));
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
                excess_pool += tile.mass[lane];
                tile.cells[lane] = Cell3::Gas;
                tile.mass[lane] = 0.0;
            }
            self.conversions.to_gas += 1;
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

        // Conservative redistribution over interface cells (2-D port).
        let mut interface_count = 0usize;
        for tile in self.tiles.values() {
            interface_count += tile
                .cells
                .iter()
                .filter(|c| **c == Cell3::Interface)
                .count();
        }
        if interface_count == 0 {
            self.carry = excess_pool;
        } else {
            #[allow(clippy::cast_precision_loss)]
            let share = excess_pool / interface_count as f64;
            for tile in self.tiles.values_mut() {
                for lane in 0..TILE_CELLS {
                    if tile.cells[lane] == Cell3::Interface {
                        tile.mass[lane] += share;
                    }
                }
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

/// Neighbor resolution outcome.
enum Neighbor {
    /// Out-of-domain or inactive-tile space: wall by sparse convention.
    Wall,
    /// A cell in an active tile.
    Active(u64, usize, Cell3),
}
