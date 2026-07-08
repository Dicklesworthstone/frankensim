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
    /// cells AND their face-neighbors reach `target`, which is exactly
    /// the uniform-band precondition ghost-penalty assembly requires.
    pub fn refine_toward_interface(&mut self, sdf: &dyn CutSdf, target: u32) {
        self.refine_where(target, &|lo: [f64; 2], hi: [f64; 2]| {
            let h = hi[0] - lo[0];
            let ilo = [(lo[0] - h).max(0.0), (lo[1] - h).max(0.0)];
            let ihi = [(hi[0] + h).min(1.0), (hi[1] + h).min(1.0)];
            sdf.enclose(ilo, ihi).contains_zero()
        });
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
}
