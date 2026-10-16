//! The quadtree background grid — the 2D restriction of the plan's
//! octree, sharing FrankenVDB's dyadic tile alignment: a cell at level
//! ℓ is the dyadic box `[i·2⁻ℓ, (i+1)·2⁻ℓ] × [j·2⁻ℓ, (j+1)·2⁻ℓ]` over
//! the unit square, so a leaf at tile depth is exactly one FrankenVDB
//! leaf face and one tree geometry serves SDF storage, LBM lattices,
//! and CutFEM cells.
//!
//! Non-uniform resolution comes free: leaves live at mixed levels
//! under a 2:1 edge balance, and the resulting hanging nodes are
//! constrained IN THE ELEMENT SPACE (midpoint = average of the coarse
//! edge's endpoints), not by mesh surgery. Node identities are lattice
//! coordinates at the finest level — dyadic, hence exactly
//! representable as f64 positions.

use crate::sdf::CutSdf;
use std::collections::BTreeSet;

/// A cell key: (level, i, j) with `0 ≤ i, j < 2^level`.
pub type CellKey = (u32, u32, u32);

/// A node key: lattice coordinates at the tree's max level,
/// `0 ≤ g ≤ 2^max_level`.
pub type NodeKey = (u32, u32);

/// Cartesian axis normal to a shared cell-face patch.
///
/// Patch orientation is always positive: [`Self::X`] points from the
/// negative-x cell to the positive-x cell, and [`Self::Y`] from the
/// negative-y cell to the positive-y cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaceAxis {
    /// Positive x-axis normal.
    X,
    /// Positive y-axis normal.
    Y,
}

impl FaceAxis {
    /// Zero-based Cartesian component index.
    #[must_use]
    pub const fn index(self) -> usize {
        match self {
            Self::X => 0,
            Self::Y => 1,
        }
    }

    /// Positive-axis unit normal.
    #[must_use]
    pub const fn normal(self) -> [f64; 2] {
        match self {
            Self::X => [1.0, 0.0],
            Self::Y => [0.0, 1.0],
        }
    }
}

/// One oriented face of a cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FaceDirection {
    /// Face with outward normal `(-1, 0)`.
    NegativeX,
    /// Face with outward normal `(1, 0)`.
    PositiveX,
    /// Face with outward normal `(0, -1)`.
    NegativeY,
    /// Face with outward normal `(0, 1)`.
    PositiveY,
}

impl FaceDirection {
    /// Deterministic legacy direction order: `-x, +x, -y, +y`.
    pub const ALL: [Self; 4] = [
        Self::NegativeX,
        Self::PositiveX,
        Self::NegativeY,
        Self::PositiveY,
    ];

    const fn axis(self) -> FaceAxis {
        match self {
            Self::NegativeX | Self::PositiveX => FaceAxis::X,
            Self::NegativeY | Self::PositiveY => FaceAxis::Y,
        }
    }

    const fn is_positive(self) -> bool {
        matches!(self, Self::PositiveX | Self::PositiveY)
    }
}

/// Exact dyadic overlap of two 2:1 leaf-cell faces.
///
/// Geometry is stored as integer coordinates on the quadtree's finest node
/// lattice. Floating-point accessors therefore perform only exact dyadic
/// scaling. `negative` and `positive` follow [`FaceAxis`]'s positive-axis
/// orientation even when their tree levels differ.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SharedFacePatch {
    negative: CellKey,
    positive: CellKey,
    axis: FaceAxis,
    lattice_level: u32,
    coordinate: u32,
    tangent_start: u32,
    tangent_end: u32,
}

impl SharedFacePatch {
    /// Cell on the negative side of the oriented face.
    #[must_use]
    pub const fn negative_cell(self) -> CellKey {
        self.negative
    }

    /// Cell on the positive side of the oriented face.
    #[must_use]
    pub const fn positive_cell(self) -> CellKey {
        self.positive
    }

    /// Positive-axis canonical cell pair.
    #[must_use]
    pub const fn oriented_cells(self) -> (CellKey, CellKey) {
        (self.negative, self.positive)
    }

    /// Lexicographically canonical cell pair used by retained evidence maps.
    #[must_use]
    pub fn canonical_cells(self) -> (CellKey, CellKey) {
        if self.negative < self.positive {
            (self.negative, self.positive)
        } else {
            (self.positive, self.negative)
        }
    }

    /// Face-normal axis. Its normal always points in the positive direction.
    #[must_use]
    pub const fn axis(self) -> FaceAxis {
        self.axis
    }

    /// Exact dyadic face coordinate.
    #[must_use]
    pub fn coordinate(self) -> f64 {
        f64::from(self.coordinate) / f64::from(1u32 << self.lattice_level)
    }

    /// Exact dyadic tangent interval, increasing along the other axis.
    #[must_use]
    pub fn tangent_interval(self) -> (f64, f64) {
        let scale = f64::from(1u32 << self.lattice_level);
        (
            f64::from(self.tangent_start) / scale,
            f64::from(self.tangent_end) / scale,
        )
    }

    /// Exact patch length.
    #[must_use]
    pub fn length(self) -> f64 {
        let (start, end) = self.tangent_interval();
        end - start
    }

    /// Ghost face size `h_F = min(h_negative, h_positive)`.
    #[must_use]
    pub fn h_f(self) -> f64 {
        let level = self.negative.0.max(self.positive.0);
        1.0 / f64::from(1u32 << level)
    }

    /// The other cell sharing this patch.
    #[must_use]
    pub fn other(self, cell: CellKey) -> Option<CellKey> {
        if cell == self.negative {
            Some(self.positive)
        } else if cell == self.positive {
            Some(self.negative)
        } else {
            None
        }
    }
}

/// Fail-closed diagnostics for malformed shared-face topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceTopologyError {
    /// The requested cell is not a leaf of this quadtree.
    NotLeaf {
        /// Requested cell.
        cell: CellKey,
    },
    /// A retained leaf has an invalid level or index for this tree.
    InvalidCell {
        /// Malformed leaf.
        cell: CellKey,
    },
    /// Two face neighbors violate the 2:1 level-difference bound.
    Unbalanced {
        /// First cell.
        cell: CellKey,
        /// Face neighbor more than one level away.
        neighbor: CellKey,
    },
    /// Two requested leaves do not share a positive-length face patch.
    NotFaceNeighbors {
        /// First cell.
        cell: CellKey,
        /// Second cell.
        neighbor: CellKey,
    },
    /// Neighbor patches overlap or leave a gap on an interior cell face.
    InvalidCoverage {
        /// Cell whose face is malformed.
        cell: CellKey,
        /// Requested face direction.
        direction: FaceDirection,
        /// First uncovered or multiply-covered finest-lattice coordinate.
        at: u32,
    },
}

impl core::fmt::Display for FaceTopologyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NotLeaf { cell } => write!(f, "cell {cell:?} is not a quadtree leaf"),
            Self::InvalidCell { cell } => write!(
                f,
                "cell {cell:?} is outside the quadtree's level/index bounds"
            ),
            Self::Unbalanced { cell, neighbor } => write!(
                f,
                "face neighbors {cell:?} and {neighbor:?} violate 2:1 balance"
            ),
            Self::NotFaceNeighbors { cell, neighbor } => write!(
                f,
                "cells {cell:?} and {neighbor:?} do not share a positive-length face"
            ),
            Self::InvalidCoverage {
                cell,
                direction,
                at,
            } => write!(
                f,
                "cell {cell:?} has overlapping or incomplete {direction:?} face coverage at finest-lattice coordinate {at}"
            ),
        }
    }
}

impl std::error::Error for FaceTopologyError {}

#[derive(Debug, Clone, Copy)]
struct LatticeRect {
    x0: u32,
    x1: u32,
    y0: u32,
    y1: u32,
}

fn same_level_face_cell(cell: CellKey, direction: FaceDirection) -> Option<CellKey> {
    let (level, i, j) = cell;
    let count = 1u32 << level;
    match direction {
        FaceDirection::NegativeX => i.checked_sub(1).map(|neighbor_i| (level, neighbor_i, j)),
        FaceDirection::PositiveX => (i + 1 < count).then_some((level, i + 1, j)),
        FaceDirection::NegativeY => j.checked_sub(1).map(|neighbor_j| (level, i, neighbor_j)),
        FaceDirection::PositiveY => (j + 1 < count).then_some((level, i, j + 1)),
    }
}

fn validate_face_coverage(
    cell: CellKey,
    direction: FaceDirection,
    tangent_start: u32,
    tangent_end: u32,
    patches: &[SharedFacePatch],
) -> Result<(), FaceTopologyError> {
    let mut cursor = tangent_start;
    for patch in patches {
        if patch.tangent_start != cursor || patch.tangent_end <= patch.tangent_start {
            return Err(FaceTopologyError::InvalidCoverage {
                cell,
                direction,
                at: cursor.min(patch.tangent_start),
            });
        }
        cursor = patch.tangent_end;
    }
    if cursor == tangent_end {
        Ok(())
    } else {
        Err(FaceTopologyError::InvalidCoverage {
            cell,
            direction,
            at: cursor,
        })
    }
}

/// A 2:1-balanced quadtree over the unit square.
#[derive(Debug, Clone)]
pub struct Quadtree {
    max_level: u32,
    leaves: BTreeSet<CellKey>,
}

impl Quadtree {
    /// A uniform grid: every leaf at `level` (which is also the node
    /// lattice level).
    #[must_use]
    pub fn uniform(level: u32) -> Quadtree {
        Quadtree::with_room(level, level)
    }

    /// A uniform grid at `base` with refinement headroom to
    /// `max_level` (the node lattice level).
    ///
    /// # Panics
    /// If `base > max_level` or `max_level > 16` (the dyadic lattice
    /// must stay comfortably inside u32/f64 exactness).
    #[must_use]
    pub fn with_room(base: u32, max_level: u32) -> Quadtree {
        assert!(base <= max_level, "base level exceeds max level");
        assert!(max_level <= 16, "max level capped at 16");
        let n = 1u32 << base;
        let mut leaves = BTreeSet::new();
        for i in 0..n {
            for j in 0..n {
                leaves.insert((base, i, j));
            }
        }
        Quadtree { max_level, leaves }
    }

    /// The node lattice level.
    #[must_use]
    pub fn max_level(&self) -> u32 {
        self.max_level
    }

    /// Leaf cells in deterministic (BTree) order.
    pub fn leaves(&self) -> impl Iterator<Item = CellKey> + '_ {
        self.leaves.iter().copied()
    }

    /// Number of leaves.
    #[must_use]
    pub fn leaf_count(&self) -> usize {
        self.leaves.len()
    }

    /// Is this cell a leaf?
    #[must_use]
    pub fn is_leaf(&self, c: CellKey) -> bool {
        self.leaves.contains(&c)
    }

    /// The cell's box `(lo, hi)` — exact dyadic f64.
    #[must_use]
    pub fn rect(&self, c: CellKey) -> ([f64; 2], [f64; 2]) {
        let (lv, i, j) = c;
        let h = 1.0 / f64::from(1u32 << lv);
        let lo = [f64::from(i) * h, f64::from(j) * h];
        (lo, [lo[0] + h, lo[1] + h])
    }

    /// Cell side length.
    #[must_use]
    pub fn cell_h(&self, c: CellKey) -> f64 {
        1.0 / f64::from(1u32 << c.0)
    }

    /// The cell's four corner nodes, counterclockwise from the lower
    /// left: (0,0), (1,0), (1,1), (0,1) in cell-local coordinates.
    #[must_use]
    pub fn corner_nodes(&self, c: CellKey) -> [NodeKey; 4] {
        let (lv, i, j) = c;
        let s = 1u32 << (self.max_level - lv);
        [
            (i * s, j * s),
            ((i + 1) * s, j * s),
            ((i + 1) * s, (j + 1) * s),
            (i * s, (j + 1) * s),
        ]
    }

    /// A node's position — exact dyadic f64.
    #[must_use]
    pub fn node_pos(&self, n: NodeKey) -> [f64; 2] {
        let h = 1.0 / f64::from(1u32 << self.max_level);
        [f64::from(n.0) * h, f64::from(n.1) * h]
    }

    /// The node lattice extent (`2^max_level`, the largest coordinate).
    #[must_use]
    pub fn node_extent(&self) -> u32 {
        1u32 << self.max_level
    }

    /// Split a leaf into its four children.
    ///
    /// # Panics
    /// If `c` is not a leaf or is already at max level.
    pub fn split(&mut self, c: CellKey) {
        assert!(self.leaves.remove(&c), "split of a non-leaf {c:?}");
        let (lv, i, j) = c;
        assert!(lv < self.max_level, "split below max level {c:?}");
        for di in 0..2u32 {
            for dj in 0..2u32 {
                self.leaves.insert((lv + 1, 2 * i + di, 2 * j + dj));
            }
        }
    }

    /// The leaf containing a strictly interior point (cell-boundary
    /// points resolve to an arbitrary abutting leaf; callers pass cell
    /// centers). `None` outside `(0,1)²`.
    #[must_use]
    pub fn find_leaf_at(&self, x: f64, y: f64) -> Option<CellKey> {
        if !(0.0..1.0).contains(&x) || !(0.0..1.0).contains(&y) {
            return None;
        }
        for lv in 0..=self.max_level {
            let n = f64::from(1u32 << lv);
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let key = (lv, (x * n).floor() as u32, (y * n).floor() as u32);
            if self.leaves.contains(&key) {
                return Some(key);
            }
        }
        None
    }

    /// The leaf covering the same-level neighbor of `c` in direction
    /// `dir` (0 −x, 1 +x, 2 −y, 3 +y); `None` at the domain boundary.
    #[must_use]
    pub fn covering_neighbor(&self, c: CellKey, dir: u8) -> Option<CellKey> {
        let (lv, i, j) = c;
        let n = i64::from(1u32 << lv);
        let (ni, nj) = match dir {
            0 => (i64::from(i) - 1, i64::from(j)),
            1 => (i64::from(i) + 1, i64::from(j)),
            2 => (i64::from(i), i64::from(j) - 1),
            _ => (i64::from(i), i64::from(j) + 1),
        };
        if ni < 0 || nj < 0 || ni >= n || nj >= n {
            return None;
        }
        let h = 1.0 / f64::from(1u32 << lv);
        #[allow(clippy::cast_precision_loss)]
        let cx = (ni as f64 + 0.5) * h;
        #[allow(clippy::cast_precision_loss)]
        let cy = (nj as f64 + 0.5) * h;
        self.find_leaf_at(cx, cy)
    }

    /// Every leaf neighbor sharing an exact positive-length patch of one cell
    /// face, ordered by increasing tangent coordinate.
    ///
    /// A coarse face has two fine neighbors under 2:1 balance; the reverse
    /// query on either fine cell returns the same single patch. Interior faces
    /// must be covered exactly once with no gaps. Domain-boundary faces return
    /// an empty vector.
    ///
    /// # Errors
    /// Refuses a non-leaf request, malformed retained cell, level jump greater
    /// than one, overlapping neighbor patches, or incomplete interior-face
    /// coverage.
    pub fn face_neighbors(
        &self,
        cell: CellKey,
        direction: FaceDirection,
    ) -> Result<Vec<SharedFacePatch>, FaceTopologyError> {
        if !self.is_leaf(cell) {
            return Err(FaceTopologyError::NotLeaf { cell });
        }
        let rect = self.lattice_rect(cell)?;
        let (face_coordinate, tangent_start, tangent_end) = match direction {
            FaceDirection::NegativeX => (rect.x0, rect.y0, rect.y1),
            FaceDirection::PositiveX => (rect.x1, rect.y0, rect.y1),
            FaceDirection::NegativeY => (rect.y0, rect.x0, rect.x1),
            FaceDirection::PositiveY => (rect.y1, rect.x0, rect.x1),
        };
        let extent = self.node_extent();
        let on_boundary = if direction.is_positive() {
            face_coordinate == extent
        } else {
            face_coordinate == 0
        };
        if on_boundary {
            return Ok(Vec::new());
        }

        let same_level =
            same_level_face_cell(cell, direction).ok_or(FaceTopologyError::InvalidCoverage {
                cell,
                direction,
                at: tangent_start,
            })?;
        let mut candidates = Vec::with_capacity(4);
        candidates.push(same_level);
        if cell.0 > 0 {
            candidates.push((cell.0 - 1, same_level.1 / 2, same_level.2 / 2));
        }
        if cell.0 < self.max_level {
            let fine_level = cell.0 + 1;
            let fine_i = 2 * same_level.1;
            let fine_j = 2 * same_level.2;
            match direction {
                FaceDirection::NegativeX => {
                    candidates.push((fine_level, fine_i + 1, fine_j));
                    candidates.push((fine_level, fine_i + 1, fine_j + 1));
                }
                FaceDirection::PositiveX => {
                    candidates.push((fine_level, fine_i, fine_j));
                    candidates.push((fine_level, fine_i, fine_j + 1));
                }
                FaceDirection::NegativeY => {
                    candidates.push((fine_level, fine_i, fine_j + 1));
                    candidates.push((fine_level, fine_i + 1, fine_j + 1));
                }
                FaceDirection::PositiveY => {
                    candidates.push((fine_level, fine_i, fine_j));
                    candidates.push((fine_level, fine_i + 1, fine_j));
                }
            }
        }

        let mut patches = Vec::with_capacity(2);
        for neighbor in candidates {
            if !self.leaves.contains(&neighbor) {
                continue;
            }
            let patch = self.shared_face_patch_inner(cell, neighbor)?.ok_or(
                FaceTopologyError::InvalidCoverage {
                    cell,
                    direction,
                    at: tangent_start,
                },
            )?;
            let cell_is_negative = patch.negative == cell;
            if patch.axis != direction.axis() || cell_is_negative != direction.is_positive() {
                return Err(FaceTopologyError::InvalidCoverage {
                    cell,
                    direction,
                    at: tangent_start,
                });
            }
            patches.push(patch);
        }
        patches.sort_unstable_by_key(|patch| {
            (
                patch.tangent_start,
                patch.tangent_end,
                patch.canonical_cells(),
            )
        });
        if let Err(coverage_error) =
            validate_face_coverage(cell, direction, tangent_start, tangent_end, &patches)
        {
            // Valid 2:1 trees never reach this scan. It exists only to retain
            // a precise `Unbalanced`/`InvalidCell` refusal for corrupted
            // private topology instead of misreporting every defect as a gap.
            self.diagnose_malformed_face(cell)?;
            return Err(coverage_error);
        }
        Ok(patches)
    }

    fn diagnose_malformed_face(&self, cell: CellKey) -> Result<(), FaceTopologyError> {
        for &neighbor in &self.leaves {
            if neighbor != cell {
                let _ = self.shared_face_patch_inner(cell, neighbor)?;
            }
        }
        Ok(())
    }

    /// Exact positive-axis overlap patch shared by two leaves.
    ///
    /// # Errors
    /// Refuses non-leaves, malformed/unbalanced cells, or cells that do not
    /// share a positive-length face.
    pub fn shared_face_patch(
        &self,
        cell: CellKey,
        neighbor: CellKey,
    ) -> Result<SharedFacePatch, FaceTopologyError> {
        if !self.is_leaf(cell) {
            return Err(FaceTopologyError::NotLeaf { cell });
        }
        if !self.is_leaf(neighbor) {
            return Err(FaceTopologyError::NotLeaf { cell: neighbor });
        }
        self.shared_face_patch_inner(cell, neighbor)?
            .ok_or(FaceTopologyError::NotFaceNeighbors { cell, neighbor })
    }

    fn shared_face_patch_inner(
        &self,
        cell: CellKey,
        neighbor: CellKey,
    ) -> Result<Option<SharedFacePatch>, FaceTopologyError> {
        let a = self.lattice_rect(cell)?;
        let b = self.lattice_rect(neighbor)?;
        let level_gap = cell.0.abs_diff(neighbor.0);

        let x_overlap_start = a.x0.max(b.x0);
        let x_overlap_end = a.x1.min(b.x1);
        let y_overlap_start = a.y0.max(b.y0);
        let y_overlap_end = a.y1.min(b.y1);
        let patch = if a.x1 == b.x0 && y_overlap_start < y_overlap_end {
            Some(SharedFacePatch {
                negative: cell,
                positive: neighbor,
                axis: FaceAxis::X,
                lattice_level: self.max_level,
                coordinate: a.x1,
                tangent_start: y_overlap_start,
                tangent_end: y_overlap_end,
            })
        } else if b.x1 == a.x0 && y_overlap_start < y_overlap_end {
            Some(SharedFacePatch {
                negative: neighbor,
                positive: cell,
                axis: FaceAxis::X,
                lattice_level: self.max_level,
                coordinate: b.x1,
                tangent_start: y_overlap_start,
                tangent_end: y_overlap_end,
            })
        } else if a.y1 == b.y0 && x_overlap_start < x_overlap_end {
            Some(SharedFacePatch {
                negative: cell,
                positive: neighbor,
                axis: FaceAxis::Y,
                lattice_level: self.max_level,
                coordinate: a.y1,
                tangent_start: x_overlap_start,
                tangent_end: x_overlap_end,
            })
        } else if b.y1 == a.y0 && x_overlap_start < x_overlap_end {
            Some(SharedFacePatch {
                negative: neighbor,
                positive: cell,
                axis: FaceAxis::Y,
                lattice_level: self.max_level,
                coordinate: b.y1,
                tangent_start: x_overlap_start,
                tangent_end: x_overlap_end,
            })
        } else {
            None
        };

        if patch.is_some() && level_gap > 1 {
            return Err(FaceTopologyError::Unbalanced { cell, neighbor });
        }
        Ok(patch)
    }

    fn lattice_rect(&self, cell: CellKey) -> Result<LatticeRect, FaceTopologyError> {
        let (level, i, j) = cell;
        if level > self.max_level {
            return Err(FaceTopologyError::InvalidCell { cell });
        }
        let count = 1u32 << level;
        if i >= count || j >= count {
            return Err(FaceTopologyError::InvalidCell { cell });
        }
        let side = 1u32 << (self.max_level - level);
        Ok(LatticeRect {
            x0: i * side,
            x1: (i + 1) * side,
            y0: j * side,
            y1: (j + 1) * side,
        })
    }

    /// Refine every leaf below `target` whose box satisfies `pred`
    /// until no such leaf remains, then restore the 2:1 balance.
    pub fn refine_where(&mut self, target: u32, pred: &dyn Fn([f64; 2], [f64; 2]) -> bool) {
        assert!(target <= self.max_level, "target beyond max level");
        loop {
            let due: Vec<CellKey> = self
                .leaves
                .iter()
                .copied()
                .filter(|&c| {
                    let (lo, hi) = self.rect(c);
                    c.0 < target && pred(lo, hi)
                })
                .collect();
            if due.is_empty() {
                break;
            }
            for c in due {
                self.split(c);
            }
        }
        self.balance();
    }

    /// Refine the interface band to `target`: any leaf whose box,
    /// INFLATED by one cell width on every side, has a certified
    /// enclosure straddling zero. The inflation guarantees that cut
    /// cells AND their face-neighbors reach `target`. A bounded closure pass
    /// then splits the coarser side of any remaining cut-cell face mismatch,
    /// including finer neighbors introduced by unrelated adaptive refinement.
    /// This re-establishes the scalar frontend's uniform-band precondition;
    /// vector elasticity can also integrate balanced coarse/fine shared
    /// patches directly.
    pub fn refine_toward_interface(&mut self, sdf: &dyn CutSdf, target: u32) {
        self.refine_where(target, &|lo: [f64; 2], hi: [f64; 2]| {
            let h = hi[0] - lo[0];
            let ilo = [(lo[0] - h).max(0.0), (lo[1] - h).max(0.0)];
            let ihi = [(hi[0] + h).min(1.0), (hi[1] + h).min(1.0)];
            sdf.enclose(ilo, ihi).contains_zero()
        });
        loop {
            let mut due = BTreeSet::new();
            for cell in self.leaves.iter().copied() {
                let (lo, hi) = self.rect(cell);
                if !sdf.enclose(lo, hi).contains_zero() {
                    continue;
                }
                for direction in 0..4 {
                    let Some(neighbor) = self.covering_neighbor(cell, direction) else {
                        continue;
                    };
                    if cell.0 < neighbor.0 {
                        due.insert(cell);
                    } else if neighbor.0 < cell.0 {
                        due.insert(neighbor);
                    }
                }
            }
            if due.is_empty() {
                break;
            }
            for cell in due {
                if self.is_leaf(cell) {
                    self.split(cell);
                }
            }
            self.balance();
        }
    }

    /// A copy with every leaf split once — the enriched-space grid the
    /// DWR estimator solves its adjoint on (the "higher-resolution"
    /// enrichment option). Headroom grows by one level.
    ///
    /// # Panics
    /// If the tree is already at the level-16 lattice cap.
    #[must_use]
    pub fn refined_once(&self) -> Quadtree {
        assert!(self.max_level < 16, "refined_once at the lattice cap");
        let mut out = Quadtree {
            max_level: self.max_level + 1,
            leaves: BTreeSet::new(),
        };
        for &(lv, i, j) in &self.leaves {
            for di in 0..2u32 {
                for dj in 0..2u32 {
                    out.leaves.insert((lv + 1, 2 * i + di, 2 * j + dj));
                }
            }
        }
        out
    }

    /// Restore the 2:1 edge balance by splitting too-coarse neighbors.
    pub fn balance(&mut self) {
        loop {
            let mut due: BTreeSet<CellKey> = BTreeSet::new();
            for &c in &self.leaves {
                if c.0 == 0 {
                    continue;
                }
                for dir in 0..4u8 {
                    if let Some(nb) = self.covering_neighbor(c, dir)
                        && nb.0 + 1 < c.0
                    {
                        due.insert(nb);
                    }
                }
            }
            if due.is_empty() {
                break;
            }
            for c in due {
                if self.is_leaf(c) {
                    self.split(c);
                }
            }
        }
    }

    /// Hanging-node constraints among `nodes` (the mesh nodes of the
    /// `active` leaves): for every active leaf edge whose integer
    /// midpoint is itself a mesh node, that midpoint hangs and is
    /// constrained to the average of the edge's endpoints. With the
    /// 2:1 balance, chains through such constraints strictly coarsen
    /// and terminate.
    #[must_use]
    pub fn hanging_constraints(
        &self,
        active: &BTreeSet<CellKey>,
        nodes: &BTreeSet<NodeKey>,
    ) -> Vec<(NodeKey, [(NodeKey, f64); 2])> {
        let mut out = Vec::new();
        for &c in active {
            let s = 1u32 << (self.max_level - c.0);
            if s < 2 {
                continue;
            }
            let corners = self.corner_nodes(c);
            for e in 0..4 {
                let a = corners[e];
                let b = corners[(e + 1) % 4];
                let m = (u32::midpoint(a.0, b.0), u32::midpoint(a.1, b.1));
                if nodes.contains(&m) {
                    out.push((m, [(a, 0.5), (b, 0.5)]));
                }
            }
        }
        out.sort_unstable_by_key(|(m, _)| *m);
        out.dedup_by_key(|(m, _)| *m);
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Circle, CutFemError, FemParams, Space};

    #[test]
    fn balance_holds_two_to_one() {
        let mut q = Quadtree::with_room(2, 5);
        // Refine one corner hard; balance must grade the rest.
        q.refine_where(5, &|lo: [f64; 2], _hi: [f64; 2]| {
            lo[0] < 0.13 && lo[1] < 0.13
        });
        for c in q.leaves().collect::<Vec<_>>() {
            for dir in 0..4u8 {
                if let Some(nb) = q.covering_neighbor(c, dir) {
                    assert!(
                        nb.0 + 1 >= c.0 && c.0 + 1 >= nb.0,
                        "2:1 violated: {c:?} vs {nb:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn leaves_partition_unit_area() {
        let mut q = Quadtree::with_room(1, 4);
        q.refine_where(4, &|lo: [f64; 2], _| lo[0] < 0.3);
        let area: f64 = q
            .leaves()
            .map(|c| {
                let h = q.cell_h(c);
                h * h
            })
            .sum();
        assert!((area - 1.0).abs() < 1e-12, "leaf area sum {area}");
    }

    #[test]
    fn interface_refinement_recloses_scalar_band_after_adjacent_split() {
        let disk = Circle {
            center: [0.5, 0.5],
            radius: 0.42,
        };
        let mut grid = Quadtree::with_room(4, 5);
        grid.split((4, 13, 7));

        let cut = (4, 14, 7);
        let finer_neighbor = (5, 27, 15);
        let (cut_lo, cut_hi) = grid.rect(cut);
        assert!(disk.enclose(cut_lo, cut_hi).contains_zero());
        assert_eq!(grid.covering_neighbor(cut, 0), Some(finer_neighbor));
        assert!(matches!(
            Space::build(&grid, &disk, FemParams::default()),
            Err(CutFemError::CutBandNotUniform { cell, neighbor })
                if cell == cut && neighbor == finer_neighbor
        ));

        // The target alone cannot repair this shape: the cut cell is already
        // at level four, while an off-band adaptive split produced level five.
        grid.refine_toward_interface(&disk, 4);
        assert!(!grid.is_leaf(cut), "the coarse cut leaf must be split");
        for cell in grid.leaves().collect::<Vec<_>>() {
            let (lo, hi) = grid.rect(cell);
            if !disk.enclose(lo, hi).contains_zero() {
                continue;
            }
            for direction in 0..4 {
                if let Some(neighbor) = grid.covering_neighbor(cell, direction) {
                    assert_eq!(
                        cell.0, neighbor.0,
                        "scalar cut band remained graded at {cell:?}/{neighbor:?}"
                    );
                }
            }
        }
        Space::build(&grid, &disk, FemParams::default())
            .expect("reclosed scalar interface band must admit ghost assembly");
    }

    #[test]
    fn coarse_face_has_two_exact_fine_reverse_patches() {
        let mut q = Quadtree::with_room(1, 2);
        q.split((1, 1, 0));

        let coarse = (1, 0, 0);
        let patches = q
            .face_neighbors(coarse, FaceDirection::PositiveX)
            .expect("balanced coarse face topology");
        assert_eq!(patches.len(), 2);
        assert_eq!(patches[0].oriented_cells(), (coarse, (2, 2, 0)));
        assert_eq!(patches[1].oriented_cells(), (coarse, (2, 2, 1)));
        assert_eq!(patches[0].axis(), FaceAxis::X);
        assert_eq!(patches[0].coordinate().to_bits(), 0.5f64.to_bits());
        let first_interval = patches[0].tangent_interval();
        let second_interval = patches[1].tangent_interval();
        assert_eq!(
            (first_interval.0.to_bits(), first_interval.1.to_bits()),
            (0.0f64.to_bits(), 0.25f64.to_bits())
        );
        assert_eq!(
            (second_interval.0.to_bits(), second_interval.1.to_bits()),
            (0.25f64.to_bits(), 0.5f64.to_bits())
        );
        assert_eq!(patches[0].h_f().to_bits(), 0.25f64.to_bits());
        assert_eq!(patches[1].h_f().to_bits(), 0.25f64.to_bits());

        for patch in patches {
            let fine = patch.positive_cell();
            let reverse = q
                .face_neighbors(fine, FaceDirection::NegativeX)
                .expect("fine reverse face topology");
            assert_eq!(reverse, vec![patch]);
            assert_eq!(q.shared_face_patch(fine, coarse), Ok(patch));
        }
    }

    #[test]
    fn coarse_y_face_has_two_exact_fine_reverse_patches() {
        let mut q = Quadtree::with_room(1, 2);
        q.split((1, 0, 1));

        let coarse = (1, 0, 0);
        let patches = q
            .face_neighbors(coarse, FaceDirection::PositiveY)
            .expect("balanced coarse y-face topology");
        assert_eq!(patches.len(), 2);
        assert_eq!(patches[0].oriented_cells(), (coarse, (2, 0, 2)));
        assert_eq!(patches[1].oriented_cells(), (coarse, (2, 1, 2)));
        assert_eq!(patches[0].axis(), FaceAxis::Y);
        assert_eq!(patches[0].coordinate().to_bits(), 0.5f64.to_bits());
        let first_interval = patches[0].tangent_interval();
        let second_interval = patches[1].tangent_interval();
        assert_eq!(
            (first_interval.0.to_bits(), first_interval.1.to_bits()),
            (0.0f64.to_bits(), 0.25f64.to_bits())
        );
        assert_eq!(
            (second_interval.0.to_bits(), second_interval.1.to_bits()),
            (0.25f64.to_bits(), 0.5f64.to_bits())
        );
        assert_eq!(patches[0].h_f().to_bits(), 0.25f64.to_bits());
        assert_eq!(patches[1].h_f().to_bits(), 0.25f64.to_bits());

        for patch in patches {
            let fine = patch.positive_cell();
            let reverse = q
                .face_neighbors(fine, FaceDirection::NegativeY)
                .expect("fine reverse y-face topology");
            assert_eq!(reverse, vec![patch]);
            assert_eq!(q.shared_face_patch(fine, coarse), Ok(patch));
        }
    }

    #[test]
    fn face_topology_refuses_unbalanced_and_incomplete_partitions() {
        let unbalanced = Quadtree {
            max_level: 3,
            leaves: BTreeSet::from([(1, 0, 0), (3, 4, 0)]),
        };
        assert!(matches!(
            unbalanced.face_neighbors((1, 0, 0), FaceDirection::PositiveX),
            Err(FaceTopologyError::Unbalanced { .. })
        ));

        let incomplete = Quadtree {
            max_level: 2,
            leaves: BTreeSet::from([(1, 0, 0)]),
        };
        assert!(matches!(
            incomplete.face_neighbors((1, 0, 0), FaceDirection::PositiveX),
            Err(FaceTopologyError::InvalidCoverage { .. })
        ));
        assert!(matches!(
            incomplete.shared_face_patch((1, 0, 0), (1, 1, 0)),
            Err(FaceTopologyError::NotLeaf { .. })
        ));
    }
}
