//! Tiled field container: dense 3D data in Morton/tile-major layout with
//! per-tile metadata, halo gathers for stencil kernels, and shard views for
//! first-touch initialization (plan §5.3).
//!
//! Layout: storage is tile-major with tiles ordered by Z-ORDER RANK and
//! cells within a tile in x-fastest row order. Tile bases are 128-byte
//! aligned whenever `edge³ * size_of::<T>()` is a multiple of 128 (all f32/
//! f64 fields at every supported edge) — the unconditional alignment policy
//! (fs-alloc's `ALLOC_ALIGN`).

use crate::tile::{TileCoord, TileError, TileGrid};

/// Boundary handling for halo gathers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary<T> {
    /// Out-of-domain cells clamp to the nearest in-domain cell.
    Clamp,
    /// Out-of-domain cells read as a constant.
    Constant(T),
    /// Coordinates wrap (periodic domain).
    Periodic,
}

/// Per-tile metadata: first-touch ownership tag, dirty flag for incremental
/// algorithms, and an occupancy hint (the hook FrankenVDB-style sparse
/// charts key off later — full per-cell masks are that bead's business).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TileMeta {
    /// Shard (CCD) that first-touch-initialized this tile; `u16::MAX` until
    /// an owner is recorded.
    pub owner_shard: u16,
    /// Set by writers that want incremental consumers to revisit this tile.
    pub dirty: bool,
    /// Occupancy hint for sparse-use consumers.
    pub occupied: bool,
}

impl Default for TileMeta {
    fn default() -> TileMeta {
        // owner_shard MUST default to the `u16::MAX` "unowned" sentinel, not the
        // derived `0` — shard 0 is a valid owner, so a derived Default would
        // silently mark every fresh TileMeta as owned by shard 0, contradicting
        // the field's documented invariant.
        TileMeta {
            owner_shard: u16::MAX,
            dirty: false,
            occupied: false,
        }
    }
}

/// 128-byte-aligned owned buffer (safe over-allocate + `align_offset`; no
/// unsafe). Alignment is best-effort for exotic element sizes and exact for
/// power-of-two sizes — asserted in tests for the supported field types.
#[derive(Debug)]
struct AlignedBuf<T> {
    raw: Vec<T>,
    off: usize,
    len: usize,
}

const BUF_ALIGN: usize = 128;

impl<T: Copy> AlignedBuf<T> {
    fn new(len: usize, fill: T) -> AlignedBuf<T> {
        let pad = if size_of::<T>() == 0 {
            0
        } else {
            BUF_ALIGN.div_ceil(size_of::<T>())
        };
        let raw = vec![fill; len + pad];
        let off = match raw.as_ptr().align_offset(BUF_ALIGN) {
            usize::MAX => 0, // impossible alignment for T: best-effort, documented
            o if o <= pad => o,
            _ => 0,
        };
        AlignedBuf { raw, off, len }
    }

    fn as_slice(&self) -> &[T] {
        &self.raw[self.off..self.off + self.len]
    }

    fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.raw[self.off..self.off + self.len]
    }
}

/// Dense tiled field over a [`TileGrid`]. `T` is the storage type; the
/// f32-storage/f64-evaluate SDF convention is [`TiledField::<f32>::eval_f64`].
#[derive(Debug)]
pub struct TiledField<T> {
    grid: TileGrid,
    data: AlignedBuf<T>,
    meta: Vec<TileMeta>,
    /// Z-order rank per linear tile index (storage slot lookup).
    rank_of_linear: Vec<u32>,
}

impl<T: Copy> TiledField<T> {
    /// Build a field over `grid`, every cell set to `fill` (including the
    /// padding cells of partial edge tiles — they are never observable
    /// through cell APIs but keep tile slices uniform).
    #[must_use]
    pub fn new(grid: TileGrid, fill: T) -> TiledField<T> {
        let vol = grid.edge().volume() as usize;
        let tiles = grid.tile_count() as usize;
        TiledField {
            rank_of_linear: grid.zorder_ranks(),
            data: AlignedBuf::new(tiles * vol, fill),
            meta: vec![
                TileMeta {
                    owner_shard: u16::MAX,
                    ..TileMeta::default()
                };
                tiles
            ],
            grid,
        }
    }

    /// Convenience: build the grid and field together.
    ///
    /// # Errors
    /// Propagates [`TileError`] from [`TileGrid::new`].
    pub fn with_dims(
        cells: [u32; 3],
        edge: crate::tile::TileEdge,
        fill: T,
    ) -> Result<TiledField<T>, TileError> {
        Ok(TiledField::new(TileGrid::new(cells, edge)?, fill))
    }

    /// The grid this field lives on.
    #[must_use]
    pub fn grid(&self) -> &TileGrid {
        &self.grid
    }

    fn slot(&self, t: TileCoord) -> usize {
        let [dx, dy, _] = self.grid.tile_dims();
        let linear = (t.z * dy + t.y) * dx + t.x;
        self.rank_of_linear[linear as usize] as usize
    }

    /// The `edge³` cells of tile `t` (x-fastest row order; partial tiles
    /// include padding cells).
    #[must_use]
    pub fn tile_slice(&self, t: TileCoord) -> &[T] {
        let vol = self.grid.edge().volume() as usize;
        let s = self.slot(t) * vol;
        &self.data.as_slice()[s..s + vol]
    }

    /// Mutable tile access; sets the dirty flag.
    pub fn tile_slice_mut(&mut self, t: TileCoord) -> &mut [T] {
        let vol = self.grid.edge().volume() as usize;
        let s = self.slot(t) * vol;
        let slot = self.slot(t);
        self.meta[slot].dirty = true;
        &mut self.data.as_mut_slice()[s..s + vol]
    }

    /// Read one cell (must be in-domain; panics on out-of-domain cells —
    /// programmer-error contract, checked before any arithmetic).
    #[must_use]
    pub fn get(&self, cell: [u32; 3]) -> T {
        let dims = self.grid.cell_dims();
        assert!(
            cell[0] < dims[0] && cell[1] < dims[1] && cell[2] < dims[2],
            "cell {cell:?} outside domain {dims:?}"
        );
        let (t, w) = self.grid.tile_of_cell(cell);
        self.tile_slice(t)[w as usize]
    }

    /// Write one cell (same domain contract as [`TiledField::get`]).
    pub fn set(&mut self, cell: [u32; 3], value: T) {
        let dims = self.grid.cell_dims();
        assert!(
            cell[0] < dims[0] && cell[1] < dims[1] && cell[2] < dims[2],
            "cell {cell:?} outside domain {dims:?}"
        );
        let (t, w) = self.grid.tile_of_cell(cell);
        self.tile_slice_mut(t)[w as usize] = value;
    }

    /// Per-tile metadata.
    #[must_use]
    pub fn meta(&self, t: TileCoord) -> TileMeta {
        self.meta[self.slot(t)]
    }

    /// Mutable per-tile metadata (ownership tags, dirty/occupancy flags).
    pub fn meta_mut(&mut self, t: TileCoord) -> &mut TileMeta {
        let slot = self.slot(t);
        &mut self.meta[slot]
    }

    /// Read a cell through a boundary policy (in-domain cells read
    /// normally; out-of-domain per `bc`). Coordinates are signed to admit
    /// ghost offsets.
    #[must_use]
    pub fn get_bc(&self, cell: [i64; 3], bc: Boundary<T>) -> T {
        let dims = self.grid.cell_dims();
        let mut c = [0u32; 3];
        for a in 0..3 {
            let d = i64::from(dims[a]);
            c[a] = match bc {
                Boundary::Clamp => cell[a].clamp(0, d - 1) as u32,
                Boundary::Periodic => cell[a].rem_euclid(d) as u32,
                Boundary::Constant(v) => {
                    if cell[a] < 0 || cell[a] >= d {
                        return v;
                    }
                    cell[a] as u32
                }
            };
        }
        self.get(c)
    }

    /// Gather the `(edge+2)³` halo block of tile `t`: the tile's cells plus
    /// a one-cell ghost layer, in deterministic z-then-y-then-x order.
    /// Ghost cells (and partial-tile padding positions) read through `bc`.
    /// `out` is cleared and refilled.
    pub fn gather_halo(&self, t: TileCoord, bc: Boundary<T>, out: &mut Vec<T>) {
        let e = i64::from(self.grid.edge().cells());
        let base = [
            i64::from(t.x) * e - 1,
            i64::from(t.y) * e - 1,
            i64::from(t.z) * e - 1,
        ];
        let side = e + 2;
        out.clear();
        out.reserve((side * side * side) as usize);
        for lz in 0..side {
            for ly in 0..side {
                for lx in 0..side {
                    out.push(self.get_bc([base[0] + lx, base[1] + ly, base[2] + lz], bc));
                }
            }
        }
    }

    /// Optimized halo gather: for tiles whose cells are all in-domain, body
    /// rows copy directly from the tile slice and only the one-cell ghost
    /// shell goes through the boundary policy; partial edge tiles fall back
    /// to the reference path. Bit-identical to [`TiledField::gather_halo`]
    /// — a conformance law, not an assumption.
    pub fn gather_halo_fast(&self, t: TileCoord, bc: Boundary<T>, out: &mut Vec<T>) {
        let e = self.grid.edge().cells();
        let dims = self.grid.cell_dims();
        let full = (t.x + 1) * e <= dims[0] && (t.y + 1) * e <= dims[1] && (t.z + 1) * e <= dims[2];
        if !full {
            self.gather_halo(t, bc, out);
            return;
        }
        let eu = e as usize;
        let side = eu + 2;
        let side_i = i64::from(e) + 2;
        let tile = self.tile_slice(t);
        out.clear();
        out.resize(side * side * side, tile[0]);
        let base = [
            i64::from(t.x) * i64::from(e) - 1,
            i64::from(t.y) * i64::from(e) - 1,
            i64::from(t.z) * i64::from(e) - 1,
        ];
        for lz in 0..side_i {
            for ly in 0..side_i {
                let (lzu, lyu) = (lz as usize, ly as usize);
                let row_start = (lzu * side + lyu) * side;
                let ghost_row = lzu == 0 || lzu == side - 1 || lyu == 0 || lyu == side - 1;
                if ghost_row {
                    for lx in 0..side_i {
                        out[row_start + lx as usize] =
                            self.get_bc([base[0] + lx, base[1] + ly, base[2] + lz], bc);
                    }
                } else {
                    // Body row: copy the tile's x-run, ghost the two ends.
                    out[row_start] = self.get_bc([base[0], base[1] + ly, base[2] + lz], bc);
                    let trow = ((lzu - 1) * eu + (lyu - 1)) * eu;
                    out[row_start + 1..row_start + 1 + eu].copy_from_slice(&tile[trow..trow + eu]);
                    out[row_start + side - 1] =
                        self.get_bc([base[0] + side_i - 1, base[1] + ly, base[2] + lz], bc);
                }
            }
        }
    }

    /// Split the field into per-shard disjoint mutable views for first-touch
    /// initialization: shard `s` receives the contiguous z-order slot range
    /// `map.slots_of(s)` and its owner tags are stamped. The caller (the
    /// executor) runs each view's initializer on the worker that owns the
    /// shard — that placement is fs-exec's contract, not this crate's.
    pub fn shard_views_mut(
        &mut self,
        map: &crate::affinity::AffinityMap,
    ) -> Vec<ShardViewMut<'_, T>> {
        assert_eq!(
            u64::from(map.tile_count()),
            self.grid.tile_count(),
            "affinity map and field must share the tile count"
        );
        let vol = self.grid.edge().volume() as usize;
        for (slot, m) in self.meta.iter_mut().enumerate() {
            m.owner_shard = map.shard_of_slot(slot as u32);
        }
        let mut views = Vec::with_capacity(map.shard_count() as usize);
        let mut rest = self.data.as_mut_slice();
        let mut consumed = 0usize;
        for shard in 0..map.shard_count() {
            let slots = map.slots_of(shard);
            let bytes = (slots.end - slots.start) as usize * vol;
            debug_assert_eq!(slots.start as usize * vol, consumed);
            let (head, tail) = rest.split_at_mut(bytes);
            views.push(ShardViewMut {
                shard,
                first_slot: slots.start,
                cells_per_tile: vol,
                data: head,
            });
            rest = tail;
            consumed += bytes;
        }
        views
    }
}

impl TiledField<f32> {
    /// The f32-storage / f64-evaluate convention of SDF charts: store
    /// narrow, evaluate wide (plan §7.2).
    #[must_use]
    pub fn eval_f64(&self, cell: [u32; 3]) -> f64 {
        f64::from(self.get(cell))
    }
}

/// One shard's disjoint mutable window of a tiled field (contiguous z-order
/// slots). `Send` so the executor can hand each view to its owning worker.
pub struct ShardViewMut<'f, T> {
    shard: u16,
    first_slot: u32,
    cells_per_tile: usize,
    data: &'f mut [T],
}

impl<T> ShardViewMut<'_, T> {
    /// The shard this view belongs to.
    #[must_use]
    pub fn shard(&self) -> u16 {
        self.shard
    }

    /// Number of tiles in this view.
    #[must_use]
    pub fn tile_count(&self) -> usize {
        self.data.len() / self.cells_per_tile
    }

    /// Mutable cells of the i-th tile in this view (its global z-order slot
    /// is `first_slot + i`).
    pub fn tile_mut(&mut self, i: usize) -> &mut [T] {
        &mut self.data[i * self.cells_per_tile..(i + 1) * self.cells_per_tile]
    }

    /// Global z-order slot of the i-th tile in this view.
    #[must_use]
    pub fn global_slot(&self, i: usize) -> u32 {
        self.first_slot + i as u32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tile::TileEdge;

    #[test]
    fn cell_get_set_roundtrip_across_tiles_and_partial_edges() {
        let mut f = TiledField::with_dims([20, 12, 9], TileEdge::E8, 0i32).expect("field");
        let dims = f.grid().cell_dims();
        let mut v = 1;
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    f.set([x, y, z], v);
                    v += 1;
                }
            }
        }
        let mut want = 1;
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    assert_eq!(f.get([x, y, z]), want);
                    want += 1;
                }
            }
        }
    }

    #[test]
    fn tile_bases_are_128_byte_aligned_for_field_types() {
        let f = TiledField::with_dims([64, 64, 64], TileEdge::E8, 0.0f32).expect("field");
        let t0 = f.tile_slice(TileCoord { x: 0, y: 0, z: 0 });
        assert_eq!(t0.as_ptr() as usize % 128, 0, "f32 E8 tile base");
        let f = TiledField::with_dims([32, 32, 32], TileEdge::E4, 0.0f64).expect("field");
        let t0 = f.tile_slice(TileCoord { x: 0, y: 0, z: 0 });
        assert_eq!(t0.as_ptr() as usize % 128, 0, "f64 E4 tile base");
    }

    #[test]
    fn dirty_and_meta_flags_track_writes() {
        let mut f = TiledField::with_dims([16, 16, 16], TileEdge::E8, 0u8).expect("field");
        let t = TileCoord { x: 1, y: 0, z: 0 };
        assert!(!f.meta(t).dirty);
        f.set([9, 0, 0], 5);
        assert!(f.meta(t).dirty);
        assert!(!f.meta(TileCoord { x: 0, y: 0, z: 0 }).dirty);
        f.meta_mut(t).occupied = true;
        assert!(f.meta(t).occupied);
    }

    #[test]
    fn f32_storage_evaluates_as_f64() {
        let mut f = TiledField::with_dims([8, 8, 8], TileEdge::E8, 0.0f32).expect("field");
        f.set([1, 2, 3], 0.5f32);
        let v: f64 = f.eval_f64([1, 2, 3]);
        assert!((v - 0.5).abs() < 1e-12);
    }

    #[test]
    fn halo_boundary_policies_agree_with_reference_semantics() {
        let mut f = TiledField::with_dims([8, 8, 8], TileEdge::E8, 0i64).expect("field");
        let dims = f.grid().cell_dims();
        for z in 0..8 {
            for y in 0..8 {
                for x in 0..8 {
                    f.set([x, y, z], i64::from((z * 64 + y * 8 + x) + 1));
                }
            }
        }
        // Clamp: ghost (-1,0,0) reads cell (0,0,0).
        assert_eq!(f.get_bc([-1, 0, 0], Boundary::Clamp), f.get([0, 0, 0]));
        // Periodic: ghost (-1,0,0) reads cell (7,0,0).
        assert_eq!(f.get_bc([-1, 0, 0], Boundary::Periodic), f.get([7, 0, 0]));
        // Constant: ghost reads the constant.
        assert_eq!(f.get_bc([-1, 0, 0], Boundary::Constant(-9)), -9);
        // Full halo has (8+2)³ cells in deterministic order; spot-check the
        // first (ghost corner) and center entries.
        let mut halo = Vec::new();
        f.gather_halo(
            TileCoord { x: 0, y: 0, z: 0 },
            Boundary::Constant(-1),
            &mut halo,
        );
        assert_eq!(halo.len(), 1000);
        assert_eq!(halo[0], -1, "ghost corner");
        let side = 10;
        let center_idx = (5 * side + 5) * side + 5; // halo-local (5,5,5) = cell (4,4,4)
        assert_eq!(halo[center_idx], f.get([4, 4, 4]));
        let _ = dims;
    }
}
