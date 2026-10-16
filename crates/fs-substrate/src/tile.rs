//! Tile geometry: the logical identity layer of the Morton/tile-major
//! layout (plan §5.3). Tiles are THE unit of scheduling, cancellation,
//! determinism, and NUMA placement; a [`TileId`] is the stable
//! `(grid-hash, z-order code)` identity that keys RNG streams,
//! deterministic reductions, and ledger events, which is what makes results
//! independent of which worker ran which tile (Decalogue P2 foundation).

use crate::morton::{MORTON_COORD_LIMIT, morton3_encode};
use core::fmt;

/// Tile edge length in cells (8³ default; 4³/16³ are the autotuner's other
/// candidates — plan §5.3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TileEdge {
    /// 4³ = 64 cells.
    E4,
    /// 8³ = 512 cells (default).
    #[default]
    E8,
    /// 16³ = 4096 cells.
    E16,
}

impl TileEdge {
    /// Edge length in cells.
    #[must_use]
    pub const fn cells(self) -> u32 {
        match self {
            TileEdge::E4 => 4,
            TileEdge::E8 => 8,
            TileEdge::E16 => 16,
        }
    }

    /// Cells per tile (edge³).
    #[must_use]
    pub const fn volume(self) -> u32 {
        let e = self.cells();
        e * e * e
    }
}

/// Tile coordinate within a [`TileGrid`] (units: tiles, not cells).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TileCoord {
    /// Tile x index.
    pub x: u32,
    /// Tile y index.
    pub y: u32,
    /// Tile z index.
    pub z: u32,
}

/// Stable logical tile identity: `(grid, code)` where `grid` hashes the
/// grid geometry and `code` is the tile coordinate's Morton code. The SAME
/// identity keys RNG streams, reduction slots, and ledger events.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct TileId {
    /// Grid-geometry hash (see [`TileGrid::grid_hash`]).
    pub grid: u64,
    /// Morton code of the tile coordinate.
    pub code: u64,
}

impl TileId {
    /// Single 128-bit key for counter-based RNG streams (Philox keying by
    /// logical identity — plan §5.2).
    #[must_use]
    pub const fn stream_key(self) -> u128 {
        ((self.grid as u128) << 64) | self.code as u128
    }
}

impl fmt::Display for TileId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:016x}/{:012x}", self.grid, self.code)
    }
}

/// Structured tile-geometry error (Decalogue P10: teaching refusals).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TileError {
    /// A cell dimension is zero.
    ZeroDim {
        /// The offending `[x, y, z]` cell dimensions.
        cells: [u32; 3],
    },
    /// The tile count on an axis exceeds the 21-bit Morton domain.
    DomainTooLarge {
        /// The offending `[x, y, z]` cell dimensions.
        cells: [u32; 3],
        /// The tile edge that was requested.
        edge: u32,
        /// Maximum cells per axis at this edge.
        max_cells_per_axis: u64,
    },
    /// The TOTAL tile count (product of the per-axis tile dimensions) exceeds
    /// the `u32` tile-index space used by `AffinityMap`/`zorder_ranks`/`slot`.
    TileCountOverflow {
        /// The offending `[x, y, z]` cell dimensions.
        cells: [u32; 3],
        /// The tile edge that was requested.
        edge: u32,
        /// The total tile count that overflowed the `u32` index space.
        tile_count: u64,
    },
}

impl fmt::Display for TileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TileError::ZeroDim { cells } => write!(
                f,
                "tile grid dimensions must be positive, got {cells:?}; empty domains have no \
                 tiles — construct nothing instead"
            ),
            TileError::DomainTooLarge {
                cells,
                edge,
                max_cells_per_axis,
            } => write!(
                f,
                "domain {cells:?} at tile edge {edge} exceeds the 21-bit Morton axis budget \
                 ({max_cells_per_axis} cells/axis); raise the tile edge or split the domain"
            ),
            TileError::TileCountOverflow {
                cells,
                edge,
                tile_count,
            } => write!(
                f,
                "domain {cells:?} at tile edge {edge} has {tile_count} tiles, which exceeds the \
                 u32 tile-index space (max {}); raise the tile edge or split the domain",
                u32::MAX
            ),
        }
    }
}

impl core::error::Error for TileError {}

/// A dense tile grid covering a 3D cell domain (ceil-division per axis; the
/// last tile on an axis may be partial — cell queries respect the true cell
/// dimensions).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TileGrid {
    cell_dims: [u32; 3],
    tile_dims: [u32; 3],
    edge: TileEdge,
    grid_hash: u64,
}

impl TileGrid {
    /// Build a grid over `cells` with the given tile edge.
    ///
    /// # Errors
    /// [`TileError::ZeroDim`] for empty domains; [`TileError::DomainTooLarge`]
    /// when the per-axis tile count would overflow the Morton domain.
    pub fn new(cells: [u32; 3], edge: TileEdge) -> Result<TileGrid, TileError> {
        if cells.contains(&0) {
            return Err(TileError::ZeroDim { cells });
        }
        let e = edge.cells();
        let tile_dims = [
            cells[0].div_ceil(e),
            cells[1].div_ceil(e),
            cells[2].div_ceil(e),
        ];
        if tile_dims.iter().any(|&t| t >= MORTON_COORD_LIMIT) {
            return Err(TileError::DomainTooLarge {
                cells,
                edge: e,
                max_cells_per_axis: u64::from(MORTON_COORD_LIMIT - 1) * u64::from(e),
            });
        }
        // The per-axis Morton budget (< 2^21 each) still allows a PRODUCT up to
        // ~2^63, but `AffinityMap`/`zorder_ranks`/`slot` index tiles with u32
        // (`grid.tile_count() as u32`), silently TRUNCATING a larger count into
        // a wrong partition. Fail closed at the u32 index-space limit. (Each
        // factor < 2^21, so the u64 product cannot itself overflow.)
        let total_tiles =
            u64::from(tile_dims[0]) * u64::from(tile_dims[1]) * u64::from(tile_dims[2]);
        if total_tiles > u64::from(u32::MAX) {
            return Err(TileError::TileCountOverflow {
                cells,
                edge: e,
                tile_count: total_tiles,
            });
        }
        let stable = format!(
            "tile-grid-v1|{}x{}x{}|edge={e}",
            cells[0], cells[1], cells[2]
        );
        Ok(TileGrid {
            cell_dims: cells,
            tile_dims,
            edge,
            grid_hash: fs_obs::fnv1a64(stable.as_bytes()),
        })
    }

    /// Cell dimensions `[x, y, z]`.
    #[must_use]
    pub const fn cell_dims(&self) -> [u32; 3] {
        self.cell_dims
    }

    /// Tile dimensions `[x, y, z]`.
    #[must_use]
    pub const fn tile_dims(&self) -> [u32; 3] {
        self.tile_dims
    }

    /// The tile edge.
    #[must_use]
    pub const fn edge(&self) -> TileEdge {
        self.edge
    }

    /// Total tile count.
    #[must_use]
    pub fn tile_count(&self) -> u64 {
        u64::from(self.tile_dims[0]) * u64::from(self.tile_dims[1]) * u64::from(self.tile_dims[2])
    }

    /// Stable hash of the grid geometry (cells + edge). Two fields over the
    /// same domain share tile identities; different geometry never collides
    /// logically (FNV collisions are a ledger-display concern only).
    #[must_use]
    pub const fn grid_hash(&self) -> u64 {
        self.grid_hash
    }

    /// The stable identity of tile `t`.
    #[must_use]
    pub fn tile_id(&self, t: TileCoord) -> TileId {
        TileId {
            grid: self.grid_hash,
            code: morton3_encode(t.x, t.y, t.z),
        }
    }

    /// Map a world cell to its tile and within-tile linear offset
    /// (x fastest, then y, then z — the SIMD-friendly row order).
    #[must_use]
    pub fn tile_of_cell(&self, cell: [u32; 3]) -> (TileCoord, u32) {
        let e = self.edge.cells();
        let t = TileCoord {
            x: cell[0] / e,
            y: cell[1] / e,
            z: cell[2] / e,
        };
        let (lx, ly, lz) = (cell[0] % e, cell[1] % e, cell[2] % e);
        (t, (lz * e + ly) * e + lx)
    }

    /// Inverse of [`TileGrid::tile_of_cell`].
    #[must_use]
    pub fn cell_of(&self, t: TileCoord, within: u32) -> [u32; 3] {
        let e = self.edge.cells();
        let lx = within % e;
        let ly = (within / e) % e;
        let lz = within / (e * e);
        [t.x * e + lx, t.y * e + ly, t.z * e + lz]
    }

    /// True when `t` touches the domain boundary.
    #[must_use]
    pub fn is_boundary(&self, t: TileCoord) -> bool {
        t.x == 0
            || t.y == 0
            || t.z == 0
            || t.x == self.tile_dims[0] - 1
            || t.y == self.tile_dims[1] - 1
            || t.z == self.tile_dims[2] - 1
    }

    /// Linear iteration order (x fastest). Deterministic.
    pub fn iter_linear(&self) -> impl Iterator<Item = TileCoord> + '_ {
        let [dx, dy, dz] = self.tile_dims;
        (0..dz).flat_map(move |z| {
            (0..dy).flat_map(move |y| (0..dx).map(move |x| TileCoord { x, y, z }))
        })
    }

    /// Z-order iteration: all tiles sorted by Morton code. Deterministic;
    /// allocates (non-power-of-two grids have gaps in code space, so
    /// enumerate-and-sort is the honest total order).
    #[must_use]
    pub fn iter_zorder(&self) -> Vec<TileCoord> {
        let mut tiles: Vec<TileCoord> = self.iter_linear().collect();
        tiles.sort_unstable_by_key(|t| morton3_encode(t.x, t.y, t.z));
        tiles
    }

    /// Boundary-first iteration: domain-boundary tiles (in z-order), then
    /// interior tiles (in z-order) — the order that lets halo/communication
    /// work start before interior compute. Deterministic.
    #[must_use]
    pub fn iter_boundary_first(&self) -> Vec<TileCoord> {
        let z = self.iter_zorder();
        let mut out = Vec::with_capacity(z.len());
        out.extend(z.iter().copied().filter(|&t| self.is_boundary(t)));
        out.extend(z.iter().copied().filter(|&t| !self.is_boundary(t)));
        out
    }

    /// The z-order rank of every tile, indexed by linear tile index
    /// (`(z * dy + y) * dx + x`). Ranks are the storage slots of the
    /// tile-major layout and the units the affinity partitioner splits.
    #[must_use]
    pub fn zorder_ranks(&self) -> Vec<u32> {
        let tiles = self.iter_zorder();
        let [dx, dy, _] = self.tile_dims;
        let mut rank_of_linear = vec![0u32; tiles.len()];
        for (rank, t) in tiles.iter().enumerate() {
            let linear = (t.z * dy + t.y) * dx + t.x;
            rank_of_linear[linear as usize] = rank as u32;
        }
        rank_of_linear
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_empty_and_oversized_domains_with_teaching_errors() {
        let e = TileGrid::new([0, 8, 8], TileEdge::E8).unwrap_err();
        assert!(matches!(e, TileError::ZeroDim { .. }), "{e}");
        let e = TileGrid::new([u32::MAX, 8, 8], TileEdge::E4).unwrap_err();
        assert!(matches!(e, TileError::DomainTooLarge { .. }), "{e}");
        assert!(e.to_string().contains("21-bit"), "{e}");
        // Each axis is within the 21-bit Morton budget (1626 < 2^21), but the
        // PRODUCT (~4.30e9) exceeds the u32 tile-index space AffinityMap/slot
        // use — reject rather than silently truncating `tile_count() as u32`.
        let e = TileGrid::new([13008, 13008, 13008], TileEdge::E8).unwrap_err();
        assert!(matches!(e, TileError::TileCountOverflow { .. }), "{e}");
        assert!(e.to_string().contains("u32 tile-index space"), "{e}");
        assert_eq!(u64::from(u32::MAX), 4_294_967_295);
    }

    #[test]
    fn cell_tile_roundtrip_including_partial_edge_tiles() {
        // 20 cells at edge 8 → 3 tiles per axis, last tile partial.
        let g = TileGrid::new([20, 20, 20], TileEdge::E8).expect("grid");
        assert_eq!(g.tile_dims(), [3, 3, 3]);
        for cell in [[0, 0, 0], [7, 7, 7], [8, 0, 19], [19, 19, 19]] {
            let (t, w) = g.tile_of_cell(cell);
            assert_eq!(g.cell_of(t, w), cell);
            assert!(w < g.edge().volume());
        }
    }

    #[test]
    fn iteration_orders_are_permutations_of_each_other() {
        let g = TileGrid::new([24, 16, 8], TileEdge::E8).expect("grid");
        let mut lin: Vec<TileCoord> = g.iter_linear().collect();
        let mut zor = g.iter_zorder();
        let mut bf = g.iter_boundary_first();
        assert_eq!(lin.len() as u64, g.tile_count());
        lin.sort_unstable();
        zor.sort_unstable();
        bf.sort_unstable();
        assert_eq!(lin, zor);
        assert_eq!(lin, bf);
    }

    #[test]
    fn boundary_first_puts_every_boundary_tile_ahead() {
        let g = TileGrid::new([32, 32, 32], TileEdge::E8).expect("grid");
        let bf = g.iter_boundary_first();
        let first_interior = bf
            .iter()
            .position(|&t| !g.is_boundary(t))
            .expect("interior");
        assert!(bf[..first_interior].iter().all(|&t| g.is_boundary(t)));
        assert!(bf[first_interior..].iter().all(|&t| !g.is_boundary(t)));
        // 4³ tile grid: 64 tiles, 8 interior.
        assert_eq!(bf.len() - first_interior, 8);
    }

    #[test]
    fn tile_ids_are_stable_and_geometry_sensitive() {
        let g1 = TileGrid::new([64, 64, 64], TileEdge::E8).expect("grid");
        let g2 = TileGrid::new([64, 64, 64], TileEdge::E8).expect("grid");
        let g3 = TileGrid::new([64, 64, 64], TileEdge::E4).expect("grid");
        let t = TileCoord { x: 1, y: 2, z: 3 };
        assert_eq!(g1.tile_id(t), g2.tile_id(t), "same geometry, same identity");
        assert_ne!(
            g1.tile_id(t).grid,
            g3.tile_id(t).grid,
            "different edge, different grid identity"
        );
        let key = g1.tile_id(t).stream_key();
        assert_eq!((key >> 64) as u64, g1.grid_hash());
        assert_eq!(key as u64, crate::morton::morton3_encode(1, 2, 3));
    }

    #[test]
    fn zorder_ranks_invert_the_zorder_iteration() {
        let g = TileGrid::new([24, 16, 8], TileEdge::E8).expect("grid");
        let ranks = g.zorder_ranks();
        let [dx, dy, _] = g.tile_dims();
        for (rank, t) in g.iter_zorder().iter().enumerate() {
            let linear = (t.z * dy + t.y) * dx + t.x;
            assert_eq!(ranks[linear as usize] as usize, rank);
        }
    }
}
