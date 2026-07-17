//! Sparse active-tile D3Q19 infrastructure (bead sjro, plan §5 + §14.1).
//!
//! Most of a pour domain is empty: the free-surface flagship and the scale
//! workstream need stream+collide over only the ACTIVE 4×4×4 tiles. This
//! module provides the Morton-ordered active-tile grid, a serial two-pass
//! sweep (collide, then pull-stream), and a TilePool-parallel sweep through
//! the [`fs_exec::KernelRunner`] seam with the same
//! bitwise-across-worker-counts guarantee the pooled FFT path proved:
//! every output tile is written by exactly one kernel tile from read-only
//! inputs, kernel tiles partition the Morton-sorted slot vector into
//! worker-count-independent groups, and the transactional destination becomes
//! published state only after both passes finish — so scheduling and a drained
//! cancellation can influence timing but never published bytes.
//!
//! Active-tile ordering is the sorted Morton key order, NEVER insertion
//! order: activation history cannot leak into iteration order, sweep order,
//! or reduction order.
//!
//! Occupancy is whole-tile in this increment: an active tile is 64 fluid
//! cells; an inactive tile behaves as solid wall (halfway bounce-back on
//! every link into it), and the domain boundary behaves the same way, so a
//! fully-enclosed seeded box conserves mass. Per-cell occupancy masks and
//! the free-surface activation/deactivation rules (WS1-E) build on top of
//! this layer in a later increment and are explicitly out of scope here.
//!
//! No-claims: this increment claims determinism and memory proportionality,
//! not throughput (the GLUP/s perf lane is bead 712t); collision is the
//! shared per-cell kernel (scalar, not the SIMD duct path); tile size is
//! pinned at the dense core's 4³ (the 4³→8³ autotune sweep is follow-on
//! work recorded through the fs-exec tune path when the perf lane lands).

use std::collections::{BTreeMap, BTreeSet};
use std::ops::ControlFlow;
use std::sync::Mutex;

use fs_exec::{CancelGate, Cancelled, Cx, KernelRunner, Reduce, RunReport, TileKernel, TilePlan};

use super::{CollisionError3, CollisionModel3, E3, OPP3, Q3, TILE, Tile, equilibrium3};

/// Number of cells in one tile (4×4×4).
const TILE_CELLS: usize = TILE * TILE * TILE;

/// Morton coordinates are interleaved from 21 bits per axis (3×21 = 63).
pub const MORTON_COORD_BITS: u32 = 21;
const MORTON_COORD_LIMIT: u32 = 1 << MORTON_COORD_BITS;

/// Consecutive Morton-sorted tile slots processed by one kernel tile. A
/// worker-count-INDEPENDENT constant: the tile plan (and therefore the
/// reduction shape) is a function of the active set alone, so worker count
/// can never change how work is grouped, only who executes it.
pub const SPARSE_SWEEP_GROUP_TILES: usize = 8;

/// Spread the low 21 bits of `v` so consecutive bits land 3 apart.
const fn spread3(v: u64) -> u64 {
    let mut x = v & 0x1f_ffff;
    x = (x | (x << 32)) & 0x1f_0000_0000_ffff;
    x = (x | (x << 16)) & 0x1f_0000_ff00_00ff;
    x = (x | (x << 8)) & 0x100f_00f0_0f00_f00f;
    x = (x | (x << 4)) & 0x10c3_0c30_c30c_30c3;
    (x | (x << 2)) & 0x1249_2492_4924_9249
}

/// Collapse every third bit of `v` back into the low 21 bits.
const fn compact3(v: u64) -> u64 {
    let mut x = v & 0x1249_2492_4924_9249;
    x = (x | (x >> 2)) & 0x10c3_0c30_c30c_30c3;
    x = (x | (x >> 4)) & 0x100f_00f0_0f00_f00f;
    x = (x | (x >> 8)) & 0x1f_0000_ff00_00ff;
    x = (x | (x >> 16)) & 0x1f_0000_0000_ffff;
    (x | (x >> 32)) & 0x1f_ffff
}

/// Morton (Z-order) key of tile coordinates: bit `3k` of the key is bit `k`
/// of `tx`, bit `3k+1` is bit `k` of `ty`, bit `3k+2` is bit `k` of `tz`.
///
/// # Panics
/// If any coordinate needs more than 21 bits.
#[must_use]
pub fn morton3(tx: u32, ty: u32, tz: u32) -> u64 {
    assert!(
        tx < MORTON_COORD_LIMIT && ty < MORTON_COORD_LIMIT && tz < MORTON_COORD_LIMIT,
        "tile coordinate exceeds the 21-bit Morton range"
    );
    spread3(u64::from(tx)) | (spread3(u64::from(ty)) << 1) | (spread3(u64::from(tz)) << 2)
}

/// Inverse of [`morton3`].
#[must_use]
#[allow(clippy::cast_possible_truncation)]
pub fn demorton3(key: u64) -> (u32, u32, u32) {
    (
        compact3(key) as u32,
        compact3(key >> 1) as u32,
        compact3(key >> 2) as u32,
    )
}

/// One active tile's populations: 19 SoA lanes of one 4×4×4 tile.
#[derive(Clone)]
struct TileBlock {
    f: [Tile; Q3],
}

impl TileBlock {
    fn equilibrium(rho: f64) -> TileBlock {
        let eq = equilibrium3(rho, [0.0; 3]);
        TileBlock {
            f: core::array::from_fn(|q| Tile([eq[q]; TILE_CELLS])),
        }
    }
}

/// Bytes of population state held per active tile (published, collided, and
/// transactional-destination buffers).
#[must_use]
pub const fn state_bytes_per_tile() -> usize {
    3 * Q3 * core::mem::size_of::<Tile>()
}

/// Wall-clock envelope and executor facts for one pooled sparse pass.
///
/// `wall_ns` includes the serving runner's scheduling and drain overhead. It
/// is measurement-only and cannot affect the published lattice state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparsePassObservation {
    /// Executor-owned report for the exact pass invocation.
    pub executor: RunReport,
    /// Saturating positive wall time around the complete runner call.
    pub wall_ns: u64,
}

/// Per-group local-stream and cross-tile halo wall envelopes.
///
/// Groups retain canonical Morton-slot identity. Their wall values may overlap
/// across workers and therefore must be composed as a task DAG, not summed as
/// if they were serial stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SparseStreamGroupObservation {
    /// Worker-count-independent group identity.
    pub group: u64,
    /// First Morton-sorted active-tile slot in this group.
    pub first_tile_slot: usize,
    /// Active tiles covered by this group.
    pub tiles: usize,
    /// Saturating positive wall time for same-tile population pulls.
    pub local_stream_wall_ns: u64,
    /// Saturating positive wall time for cross-tile/domain lookup and pulls.
    pub halo_wall_ns: u64,
}

/// Opt-in timing evidence from one successfully published pooled sweep.
///
/// The ordinary [`SparseGrid3::step_pooled`] path creates no timers or timing
/// buffers. The observed path retains the same two executor passes and the
/// same failure-atomic publication boundary, while conservatively measuring
/// the split local-stream/halo loops used by both paths.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseSweepObservation {
    /// Active tile count used by both passes.
    pub active_tiles: usize,
    /// Normalized serving-worker count.
    pub workers: usize,
    /// Collide-pass executor and wall envelope.
    pub collide: SparsePassObservation,
    /// Combined local-stream/halo pass executor and wall envelope.
    pub stream: SparsePassObservation,
    /// Canonical per-group local-stream/halo observations.
    pub stream_groups: Vec<SparseStreamGroupObservation>,
    /// Saturating positive wall time for the final publication swap/counter.
    pub publication_wall_ns: u64,
}

/// Typed refusal from sparse-grid construction or a sweep.
#[derive(Debug, Clone, PartialEq)]
pub enum SparseError3 {
    /// Domain dimensions must be positive multiples of the tile edge.
    Dims { nx: usize, ny: usize, nz: usize },
    /// A tile coordinate lies outside the declared domain.
    TileOutOfDomain { tx: u32, ty: u32, tz: u32 },
    /// The per-cell collision kernel refused; the cell is identified by
    /// Morton tile key and lane for forensics.
    Collision {
        /// Morton key of the refusing tile.
        tile_key: u64,
        /// Lane (0..64) of the refusing cell inside the tile.
        lane: usize,
        /// The underlying collision refusal.
        source: CollisionError3,
    },
    /// The pooled sweep was cancelled before completion; the pre-sweep
    /// state is intact and the step may be re-issued deterministically.
    Cancelled,
    /// The kernel pool refused the run.
    Pool(String),
    /// A successful observed pass did not retain every canonical group row.
    IncompleteObservation {
        /// Groups required by the immutable active-set plan.
        expected_groups: usize,
        /// Complete group rows recovered before publication.
        observed_groups: usize,
    },
}

impl core::fmt::Display for SparseError3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SparseError3::Dims { nx, ny, nz } => write!(
                f,
                "sparse domain dims ({nx},{ny},{nz}) must be positive multiples of {TILE}"
            ),
            SparseError3::TileOutOfDomain { tx, ty, tz } => {
                write!(f, "tile ({tx},{ty},{tz}) lies outside the declared domain")
            }
            SparseError3::Collision {
                tile_key,
                lane,
                source,
            } => write!(
                f,
                "collision refused at tile key {tile_key:#x} lane {lane}: {source}"
            ),
            SparseError3::Cancelled => f.write_str("sparse sweep cancelled before completion"),
            SparseError3::Pool(detail) => write!(f, "kernel pool refused the sweep: {detail}"),
            SparseError3::IncompleteObservation {
                expected_groups,
                observed_groups,
            } => write!(
                f,
                "sparse sweep observed {observed_groups}/{expected_groups} stream groups"
            ),
        }
    }
}

impl std::error::Error for SparseError3 {}

/// First refusal in ascending kernel-tile order — the deterministic error
/// reduction for pooled sweeps ([`Reduce::merge`] is always applied in
/// ascending tile order, so "first" is well-defined and scheduling-free).
struct FirstRefusal(Option<SparseError3>);

impl Reduce for FirstRefusal {
    fn identity() -> Self {
        FirstRefusal(None)
    }

    fn merge(self, other: Self) -> Self {
        FirstRefusal(self.0.or(other.0))
    }
}

/// Morton-ordered sparse active-tile D3Q19 grid.
///
/// State exists only for active tiles; memory is proportional to the active
/// set ([`SparseGrid3::allocated_state_bytes`]). Sweeps are two passes —
/// collide (own-tile only), then pull-stream (own tile + face/edge
/// neighbors' post-collision state) — followed by one failure-atomic
/// publication swap, so every pass writes each output exactly once from
/// read-only inputs.
pub struct SparseGrid3 {
    ntx: usize,
    nty: usize,
    ntz: usize,
    model: CollisionModel3,
    force: [f64; 3],
    /// Ascending Morton keys of the active tiles; parallel to the state
    /// vectors. THE canonical iteration order.
    keys: Vec<u64>,
    /// Morton key → slot index into `keys`/`pre`/`post`/`next`.
    index: BTreeMap<u64, usize>,
    pre: Vec<TileBlock>,
    post: Vec<TileBlock>,
    /// Stream destination kept private until a complete pass commits.
    next: Vec<TileBlock>,
    steps: u64,
}

impl core::fmt::Debug for SparseGrid3 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("SparseGrid3")
            .field("tile_dims", &(self.ntx, self.nty, self.ntz))
            .field("active_tiles", &self.keys.len())
            .field("steps", &self.steps)
            .field("model", &self.model)
            .field("force", &self.force)
            .finish_non_exhaustive()
    }
}

impl SparseGrid3 {
    /// An empty sparse grid over an `nx × ny × nz` cell domain (each a
    /// positive multiple of the 4-cell tile edge) with BGK relaxation time
    /// `tau` and a uniform body force.
    ///
    /// # Errors
    /// [`SparseError3::Dims`] for inadmissible dimensions.
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        tau: f64,
        force: [f64; 3],
    ) -> Result<SparseGrid3, SparseError3> {
        if nx == 0 || ny == 0 || nz == 0 || nx % TILE != 0 || ny % TILE != 0 || nz % TILE != 0 {
            return Err(SparseError3::Dims { nx, ny, nz });
        }
        Ok(SparseGrid3 {
            ntx: nx / TILE,
            nty: ny / TILE,
            ntz: nz / TILE,
            model: CollisionModel3::Bgk { tau },
            force,
            keys: Vec::new(),
            index: BTreeMap::new(),
            pre: Vec::new(),
            post: Vec::new(),
            next: Vec::new(),
            steps: 0,
        })
    }

    /// Activate tiles (idempotent), initializing new tiles at rest-density
    /// equilibrium. Existing tile state is preserved; the active list is
    /// re-sorted into Morton order, so iteration order depends only on the
    /// final active SET, never on activation history.
    ///
    /// # Errors
    /// [`SparseError3::TileOutOfDomain`] if any coordinate lies outside the
    /// domain (no partial activation is applied).
    pub fn activate_tiles(&mut self, tiles: &[(u32, u32, u32)]) -> Result<(), SparseError3> {
        for &(tx, ty, tz) in tiles {
            if tx as usize >= self.ntx || ty as usize >= self.nty || tz as usize >= self.ntz {
                return Err(SparseError3::TileOutOfDomain { tx, ty, tz });
            }
        }
        let mut merged: BTreeMap<u64, Option<usize>> = self
            .index
            .iter()
            .map(|(&key, &slot)| (key, Some(slot)))
            .collect();
        for &(tx, ty, tz) in tiles {
            merged.entry(morton3(tx, ty, tz)).or_insert(None);
        }
        let mut keys = Vec::with_capacity(merged.len());
        let mut pre = Vec::with_capacity(merged.len());
        let mut post = Vec::with_capacity(merged.len());
        let mut next = Vec::with_capacity(merged.len());
        let mut index = BTreeMap::new();
        for (slot, (key, old)) in merged.into_iter().enumerate() {
            keys.push(key);
            match old {
                Some(old_slot) => {
                    pre.push(self.pre[old_slot].clone());
                    post.push(self.post[old_slot].clone());
                    next.push(self.next[old_slot].clone());
                }
                None => {
                    pre.push(TileBlock::equilibrium(1.0));
                    post.push(TileBlock::equilibrium(1.0));
                    next.push(TileBlock::equilibrium(1.0));
                }
            }
            index.insert(key, slot);
        }
        self.keys = keys;
        self.pre = pre;
        self.post = post;
        self.next = next;
        self.index = index;
        Ok(())
    }

    /// Deactivate tiles, dropping their state and returning the exact mass
    /// removed (summed in canonical Morton/q/lane order) — the WS1-E
    /// mass-ledger hook: a free-surface activation rule that retires tiles
    /// must account for this mass, never silently drop it. Deactivating a
    /// tile that is not active is a no-op contributing zero mass. Surviving
    /// tiles keep their state bit-exactly; the active list stays
    /// Morton-sorted.
    ///
    /// # Errors
    /// [`SparseError3::TileOutOfDomain`] if any coordinate lies outside the
    /// domain (no partial deactivation is applied).
    pub fn deactivate_tiles(&mut self, tiles: &[(u32, u32, u32)]) -> Result<f64, SparseError3> {
        for &(tx, ty, tz) in tiles {
            if tx as usize >= self.ntx || ty as usize >= self.nty || tz as usize >= self.ntz {
                return Err(SparseError3::TileOutOfDomain { tx, ty, tz });
            }
        }
        let mut retire: BTreeSet<u64> = BTreeSet::new();
        for &(tx, ty, tz) in tiles {
            retire.insert(morton3(tx, ty, tz));
        }
        let mut removed_mass = 0.0;
        for &key in &retire {
            if let Some(&slot) = self.index.get(&key) {
                let block = &self.pre[slot];
                for q in 0..Q3 {
                    for lane in 0..TILE_CELLS {
                        removed_mass += block.f[q].0[lane];
                    }
                }
            }
        }
        let survivors: Vec<usize> = (0..self.keys.len())
            .filter(|&slot| !retire.contains(&self.keys[slot]))
            .collect();
        let mut keys = Vec::with_capacity(survivors.len());
        let mut pre = Vec::with_capacity(survivors.len());
        let mut post = Vec::with_capacity(survivors.len());
        let mut next = Vec::with_capacity(survivors.len());
        let mut index = BTreeMap::new();
        for (new_slot, &old_slot) in survivors.iter().enumerate() {
            keys.push(self.keys[old_slot]);
            pre.push(self.pre[old_slot].clone());
            post.push(self.post[old_slot].clone());
            next.push(self.next[old_slot].clone());
            index.insert(self.keys[old_slot], new_slot);
        }
        self.keys = keys;
        self.pre = pre;
        self.post = post;
        self.next = next;
        self.index = index;
        Ok(removed_mass)
    }

    /// Number of active tiles.
    #[must_use]
    pub fn active_tiles(&self) -> usize {
        self.keys.len()
    }

    /// Completed sweep steps.
    #[must_use]
    pub fn steps(&self) -> u64 {
        self.steps
    }

    /// Bytes of population state currently allocated — exactly
    /// proportional to the active set.
    #[must_use]
    pub fn allocated_state_bytes(&self) -> usize {
        (self.pre.len() + self.post.len() + self.next.len()) * Q3 * core::mem::size_of::<Tile>()
    }

    /// Whether the tile at Morton key `key` is active.
    #[must_use]
    pub fn is_active(&self, tx: u32, ty: u32, tz: u32) -> bool {
        self.index.contains_key(&morton3(tx, ty, tz))
    }

    /// Deterministic seeded perturbation: every cell's density is offset by
    /// a splitmix64 hash of `(seed, morton key, lane)` scaled to
    /// `±amplitude`, and populations reset to equilibrium at that density.
    /// The hash IS the seed schedule — no RNG state, bitwise reproducible.
    pub fn perturb(&mut self, seed: u64, amplitude: f64) {
        fn splitmix64(mut z: u64) -> u64 {
            z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        }
        for (slot, &key) in self.keys.iter().enumerate() {
            for lane in 0..TILE_CELLS {
                let h = splitmix64(seed ^ key.wrapping_mul(0x100_0000_01b3) ^ (lane as u64));
                #[allow(clippy::cast_precision_loss)]
                let unit = (h >> 11) as f64 / (1u64 << 53) as f64;
                let rho = 1.0 + amplitude * (2.0 * unit - 1.0);
                let eq = equilibrium3(rho, [0.0; 3]);
                for q in 0..Q3 {
                    self.pre[slot].f[q].0[lane] = eq[q];
                }
            }
        }
    }

    /// Total mass over the active set, summed in canonical Morton/lane/q
    /// order (deterministic bytes).
    #[must_use]
    pub fn total_mass(&self) -> f64 {
        let mut sum = 0.0;
        for block in &self.pre {
            for q in 0..Q3 {
                for lane in 0..TILE_CELLS {
                    sum += block.f[q].0[lane];
                }
            }
        }
        sum
    }

    /// Exact bits of the full population state in canonical order, for
    /// golden hashing and bitwise comparison across worker counts.
    #[must_use]
    pub fn state_bits(&self) -> Vec<u64> {
        let mut bits = Vec::with_capacity(self.pre.len() * Q3 * TILE_CELLS);
        for block in &self.pre {
            for q in 0..Q3 {
                for lane in 0..TILE_CELLS {
                    bits.push(block.f[q].0[lane].to_bits());
                }
            }
        }
        bits
    }

    /// Density and velocity of the cell at `lane` of active slot `slot`
    /// (test/diagnostic surface; canonical slot order is Morton order).
    #[must_use]
    pub fn cell_macros(&self, slot: usize, lane: usize) -> (f64, [f64; 3]) {
        let mut rho = 0.0;
        let mut mom = [0.0; 3];
        for q in 0..Q3 {
            let f = self.pre[slot].f[q].0[lane];
            rho += f;
            mom[0] += f * f64::from(E3[q].0);
            mom[1] += f * f64::from(E3[q].1);
            mom[2] += f * f64::from(E3[q].2);
        }
        (rho, [mom[0] / rho, mom[1] / rho, mom[2] / rho])
    }

    /// Ascending Morton keys of the active tiles (the canonical slot
    /// order) — the free-surface layer's topology view.
    pub(super) fn active_keys(&self) -> &[u64] {
        &self.keys
    }

    /// Slot index of an active tile by Morton key.
    pub(super) fn slot_of(&self, key: u64) -> Option<usize> {
        self.index.get(&key).copied()
    }

    /// Domain extent in tiles.
    pub(super) fn tile_dims(&self) -> (usize, usize, usize) {
        (self.ntx, self.nty, self.ntz)
    }

    /// The grid's collision model and body force (the free-surface step
    /// shares the exact per-cell collision authority).
    pub(super) fn collision(&self) -> (CollisionModel3, [f64; 3]) {
        (self.model, self.force)
    }

    /// Published populations of one cell.
    pub(super) fn populations(&self, slot: usize, lane: usize) -> [f64; Q3] {
        core::array::from_fn(|q| self.pre[slot].f[q].0[lane])
    }

    /// Overwrite the published populations of one cell (free-surface
    /// streaming commits through this; the sparse sweep never does).
    pub(super) fn set_populations(&mut self, slot: usize, lane: usize, f: [f64; Q3]) {
        for q in 0..Q3 {
            self.pre[slot].f[q].0[lane] = f[q];
        }
    }

    /// Resolve a global cell coordinate to `(active slot, lane)`, or
    /// `None` for out-of-domain / inactive-tile space (wall semantics).
    pub(super) fn resolve_source(&self, sx: i64, sy: i64, sz: i64) -> Option<(usize, usize)> {
        source_slot(sx, sy, sz, &self.index, (self.ntx, self.nty, self.ntz))
    }

    /// One serial sweep step: collide every active cell, pull-stream into a
    /// private transactional destination, then publish by swapping buffers.
    /// The serial path is the bitwise reference the pooled path must reproduce.
    ///
    /// # Errors
    /// [`SparseError3::Collision`] fail-closed on the first refusing cell
    /// in canonical order; the pre-step state is left intact.
    pub fn step_serial(&mut self) -> Result<(), SparseError3> {
        for slot in 0..self.keys.len() {
            collide_block(
                &self.pre[slot],
                &mut self.post[slot],
                self.keys[slot],
                self.model,
                self.force,
            )?;
        }
        for slot in 0..self.keys.len() {
            stream_local_block(&mut self.next[slot], slot, &self.post);
            stream_halo_block(
                &mut self.next[slot],
                slot,
                &self.post,
                &self.keys,
                &self.index,
                (self.ntx, self.nty, self.ntz),
            );
        }
        std::mem::swap(&mut self.pre, &mut self.next);
        self.steps += 1;
        Ok(())
    }

    /// One pooled sweep step through any [`KernelRunner`]: identical bytes
    /// to [`SparseGrid3::step_serial`] for every worker count, because each
    /// kernel tile writes a disjoint Morton-contiguous group of slots from
    /// read-only inputs and the group size is worker-count-independent.
    ///
    /// # Errors
    /// [`SparseError3::Collision`] (first refusal in canonical order),
    /// [`SparseError3::Cancelled`] if the gate trips (pre-step state
    /// intact; re-issuing the step after cancellation is deterministic),
    /// [`SparseError3::Pool`] for pool-level refusals.
    pub fn step_pooled<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
    ) -> Result<(), SparseError3> {
        self.step_pooled_inner(pool, gate, false).map(|_| ())
    }

    /// One pooled sweep with opt-in measurement-only pass/group observation.
    ///
    /// This uses the same worker-count-independent group plan, collide pass,
    /// combined stream/halo pass, and failure-atomic publication boundary as
    /// [`SparseGrid3::step_pooled`]. Timing adds one pair of clock reads per
    /// completed stream group and around each executor pass. The returned
    /// values are envelope-class telemetry; no numerical result depends on
    /// them.
    ///
    /// # Errors
    /// As [`SparseGrid3::step_pooled`], plus
    /// [`SparseError3::IncompleteObservation`] if a successful executor pass
    /// fails to retain every canonical group row. In every refusal case the
    /// pre-step published state remains intact.
    pub fn step_pooled_observed<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
    ) -> Result<SparseSweepObservation, SparseError3> {
        self.step_pooled_inner(pool, gate, true)?
            .ok_or(SparseError3::IncompleteObservation {
                expected_groups: self.keys.len().div_ceil(SPARSE_SWEEP_GROUP_TILES).max(1),
                observed_groups: 0,
            })
    }

    fn step_pooled_inner<P: KernelRunner>(
        &mut self,
        pool: &P,
        gate: &CancelGate,
        observed: bool,
    ) -> Result<Option<SparseSweepObservation>, SparseError3> {
        let group_count = self.keys.len().div_ceil(SPARSE_SWEEP_GROUP_TILES).max(1);
        let groups = group_count as u64;
        let stream_timing = observed.then(|| {
            (0..group_count)
                .map(|_| Mutex::new(None::<SparseStreamGroupObservation>))
                .collect::<Vec<_>>()
        });
        let active_tiles = self.keys.len();
        let workers = pool.workers();

        let collide_started = observed.then(std::time::Instant::now);
        let (collide_outcome, collide_executor) = {
            let kernel = CollideKernel {
                pre: &self.pre,
                chunks: self
                    .post
                    .chunks_mut(SPARSE_SWEEP_GROUP_TILES)
                    .map(Mutex::new)
                    .collect(),
                keys: &self.keys,
                model: self.model,
                force: self.force,
                groups,
            };
            run_sweep(pool, gate, &kernel)
        };
        let collide_wall_ns = collide_started.map(|start| nonzero_wall_ns(start.elapsed()));
        collide_outcome?;

        let stream_started = observed.then(std::time::Instant::now);
        let (stream_outcome, stream_executor) = {
            let kernel = StreamKernel {
                src: &self.post,
                chunks: self
                    .next
                    .chunks_mut(SPARSE_SWEEP_GROUP_TILES)
                    .map(Mutex::new)
                    .collect(),
                keys: &self.keys,
                index: &self.index,
                tile_dims: (self.ntx, self.nty, self.ntz),
                groups,
                timings: stream_timing.as_deref(),
            };
            run_sweep(pool, gate, &kernel)
        };
        let stream_wall_ns = stream_started.map(|start| nonzero_wall_ns(start.elapsed()));
        stream_outcome?;

        let mut stream_groups = Vec::new();
        if let Some(timing) = stream_timing {
            stream_groups.reserve(timing.len());
            for slot in timing {
                if let Ok(Some(group)) = slot.into_inner() {
                    stream_groups.push(group);
                }
            }
            if stream_groups.len() != group_count {
                return Err(SparseError3::IncompleteObservation {
                    expected_groups: group_count,
                    observed_groups: stream_groups.len(),
                });
            }
        }

        let publication_started = observed.then(std::time::Instant::now);
        std::mem::swap(&mut self.pre, &mut self.next);
        self.steps += 1;
        let publication_wall_ns = publication_started.map(|start| nonzero_wall_ns(start.elapsed()));

        Ok(
            match (collide_wall_ns, stream_wall_ns, publication_wall_ns) {
                (Some(collide_wall_ns), Some(stream_wall_ns), Some(publication_wall_ns)) => {
                    Some(SparseSweepObservation {
                        active_tiles,
                        workers,
                        collide: SparsePassObservation {
                            executor: collide_executor,
                            wall_ns: collide_wall_ns,
                        },
                        stream: SparsePassObservation {
                            executor: stream_executor,
                            wall_ns: stream_wall_ns,
                        },
                        stream_groups,
                        publication_wall_ns,
                    })
                }
                _ => None,
            },
        )
    }
}

/// Drive one pass and normalize pool/cancel/collision refusals.
fn run_sweep<P: KernelRunner, K: TileKernel<Out = FirstRefusal>>(
    pool: &P,
    gate: &CancelGate,
    kernel: &K,
) -> (Result<(), SparseError3>, RunReport) {
    let (outcome, report) = pool.run_with_gate(kernel, gate);
    let outcome = match outcome {
        Ok(FirstRefusal(None)) => Ok(()),
        Ok(FirstRefusal(Some(refusal))) => Err(refusal),
        Err(err) => {
            if gate.is_requested() {
                Err(SparseError3::Cancelled)
            } else {
                Err(SparseError3::Pool(format!("{err:?}")))
            }
        }
    };
    (outcome, report)
}

/// Collide every cell of `pre` into `post` (own-tile data only).
fn collide_block(
    pre: &TileBlock,
    post: &mut TileBlock,
    tile_key: u64,
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<(), SparseError3> {
    for lane in 0..TILE_CELLS {
        let populations: [f64; Q3] = core::array::from_fn(|q| pre.f[q].0[lane]);
        let collided = collide_or_refuse(populations, model, force).map_err(|source| {
            SparseError3::Collision {
                tile_key,
                lane,
                source,
            }
        })?;
        for q in 0..Q3 {
            post.f[q].0[lane] = collided[q];
        }
    }
    Ok(())
}

fn collide_or_refuse(
    populations: [f64; Q3],
    model: CollisionModel3,
    force: [f64; 3],
) -> Result<[f64; Q3], CollisionError3> {
    super::collide_cell3(populations, model, force)
}

/// Same-tile half of one pull-stream group. This loop has no sparse-map
/// lookup: every selected source lane is known to live in `slot`.
fn stream_local_block(dst: &mut TileBlock, slot: usize, src: &[TileBlock]) {
    for lz in 0..TILE {
        for ly in 0..TILE {
            for lx in 0..TILE {
                let lane = (lz * TILE + ly) * TILE + lx;
                for q in 0..Q3 {
                    if let Some(src_lane) = local_source_lane(lx, ly, lz, q) {
                        dst.f[q].0[lane] = src[slot].f[q].0[src_lane];
                    }
                }
            }
        }
    }
}

/// Cross-tile/domain half of one pull-stream group. A link whose source cell
/// lies in an inactive tile (or outside the domain) bounces back: the cell
/// keeps its own post-collision opposite population — inactive space is solid
/// wall in this increment. These writes are disjoint from
/// [`stream_local_block`].
fn stream_halo_block(
    dst: &mut TileBlock,
    slot: usize,
    src: &[TileBlock],
    keys: &[u64],
    index: &BTreeMap<u64, usize>,
    tile_dims: (usize, usize, usize),
) {
    let (tx, ty, tz) = demorton3(keys[slot]);
    for lz in 0..TILE {
        for ly in 0..TILE {
            for lx in 0..TILE {
                let lane = (lz * TILE + ly) * TILE + lx;
                let gx = tx as i64 * TILE as i64 + lx as i64;
                let gy = ty as i64 * TILE as i64 + ly as i64;
                let gz = tz as i64 * TILE as i64 + lz as i64;
                for q in 0..Q3 {
                    if local_source_lane(lx, ly, lz, q).is_some() {
                        continue;
                    }
                    let sx = gx - i64::from(E3[q].0);
                    let sy = gy - i64::from(E3[q].1);
                    let sz = gz - i64::from(E3[q].2);
                    dst.f[q].0[lane] = match source_slot(sx, sy, sz, index, tile_dims) {
                        Some((src_slot, src_lane)) => src[src_slot].f[q].0[src_lane],
                        None => src[slot].f[OPP3[q]].0[lane],
                    };
                }
            }
        }
    }
}

/// Source lane when a pull remains inside the destination tile.
fn local_source_lane(lx: usize, ly: usize, lz: usize, q: usize) -> Option<usize> {
    let sx = lx as i64 - i64::from(E3[q].0);
    let sy = ly as i64 - i64::from(E3[q].1);
    let sz = lz as i64 - i64::from(E3[q].2);
    let tile = TILE as i64;
    if !(0..tile).contains(&sx) || !(0..tile).contains(&sy) || !(0..tile).contains(&sz) {
        return None;
    }
    Some(((sz as usize * TILE + sy as usize) * TILE) + sx as usize)
}

/// Resolve a global cell coordinate to `(active slot, lane)`, or `None` if
/// it lies outside the domain or in an inactive tile.
fn source_slot(
    sx: i64,
    sy: i64,
    sz: i64,
    index: &BTreeMap<u64, usize>,
    (ntx, nty, ntz): (usize, usize, usize),
) -> Option<(usize, usize)> {
    if sx < 0 || sy < 0 || sz < 0 {
        return None;
    }
    let (sx, sy, sz) = (sx as usize, sy as usize, sz as usize);
    let (stx, sty, stz) = (sx / TILE, sy / TILE, sz / TILE);
    if stx >= ntx || sty >= nty || stz >= ntz {
        return None;
    }
    #[allow(clippy::cast_possible_truncation)]
    let key = morton3(stx as u32, sty as u32, stz as u32);
    let slot = *index.get(&key)?;
    let lane = ((sz % TILE) * TILE + (sy % TILE)) * TILE + (sx % TILE);
    Some((slot, lane))
}

/// Collide pass: kernel tile `g` collides the Morton-contiguous slot group
/// `g` from the read-only pre buffer into its exclusively-owned post chunk.
struct CollideKernel<'a> {
    pre: &'a [TileBlock],
    chunks: Vec<Mutex<&'a mut [TileBlock]>>,
    keys: &'a [u64],
    model: CollisionModel3,
    force: [f64; 3],
    groups: u64,
}

impl TileKernel for CollideKernel<'_> {
    type Out = FirstRefusal;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("fs-lbm/d3q19-sparse-collide", self.groups)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, FirstRefusal> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let group = tile as usize;
        let base = group * SPARSE_SWEEP_GROUP_TILES;
        let mut chunk = self.chunks[group].lock().expect("collide chunk poisoned");
        for (offset, post) in chunk.iter_mut().enumerate() {
            let slot = base + offset;
            if let Err(refusal) = collide_block(
                &self.pre[slot],
                post,
                self.keys[slot],
                self.model,
                self.force,
            ) {
                return ControlFlow::Continue(FirstRefusal(Some(refusal)));
            }
        }
        ControlFlow::Continue(FirstRefusal(None))
    }
}

/// Stream pass: kernel tile `g` pull-streams slot group `g` from the
/// read-only post-collision buffer into an exclusively-owned transactional
/// destination chunk. The caller publishes that buffer only after every
/// kernel tile completes.
struct StreamKernel<'a> {
    src: &'a [TileBlock],
    chunks: Vec<Mutex<&'a mut [TileBlock]>>,
    keys: &'a [u64],
    index: &'a BTreeMap<u64, usize>,
    tile_dims: (usize, usize, usize),
    groups: u64,
    timings: Option<&'a [Mutex<Option<SparseStreamGroupObservation>>]>,
}

impl TileKernel for StreamKernel<'_> {
    type Out = FirstRefusal;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("fs-lbm/d3q19-sparse-stream", self.groups)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, FirstRefusal> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let group = tile as usize;
        let base = group * SPARSE_SWEEP_GROUP_TILES;
        let mut chunk = self.chunks[group].lock().expect("stream chunk poisoned");
        if let Some(timings) = self.timings {
            let local_started = std::time::Instant::now();
            for (offset, dst) in chunk.iter_mut().enumerate() {
                stream_local_block(dst, base + offset, self.src);
            }
            let local_stream_wall_ns = nonzero_wall_ns(local_started.elapsed());

            let halo_started = std::time::Instant::now();
            for (offset, dst) in chunk.iter_mut().enumerate() {
                stream_halo_block(
                    dst,
                    base + offset,
                    self.src,
                    self.keys,
                    self.index,
                    self.tile_dims,
                );
            }
            let halo_wall_ns = nonzero_wall_ns(halo_started.elapsed());
            *timings[group].lock().expect("stream timing slot poisoned") =
                Some(SparseStreamGroupObservation {
                    group: tile,
                    first_tile_slot: base,
                    tiles: chunk.len(),
                    local_stream_wall_ns,
                    halo_wall_ns,
                });
        } else {
            for (offset, dst) in chunk.iter_mut().enumerate() {
                stream_local_block(dst, base + offset, self.src);
            }
            for (offset, dst) in chunk.iter_mut().enumerate() {
                stream_halo_block(
                    dst,
                    base + offset,
                    self.src,
                    self.keys,
                    self.index,
                    self.tile_dims,
                );
            }
        }
        ControlFlow::Continue(FirstRefusal(None))
    }
}

/// Envelope-class ns from a duration. A completed non-empty phase is clamped
/// to one nanosecond when the host clock's resolution reports zero.
fn nonzero_wall_ns(elapsed: std::time::Duration) -> u64 {
    u64::try_from(elapsed.as_nanos()).unwrap_or(u64::MAX).max(1)
}
