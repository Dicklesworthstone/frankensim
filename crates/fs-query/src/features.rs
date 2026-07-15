//! Feature complexes and conservative CCD candidate enumeration
//! (bead rjnd, E1 query upgrades, part 3).
//!
//! A [`FeatureComplex`] is the typed vertex/edge/face decomposition of
//! a triangle boundary, each feature carrying an outward-rounded AABB.
//! [`ccd_candidates`] enumerates feature pairs from two complexes
//! whose motion-inflated boxes overlap, through a deterministic
//! median-split BVH.
//!
//! The guarantee is one-sided and conservative: any feature pair whose
//! true swept geometry could come within the caller's declared motion
//! inflations is INCLUDED (boxes are outward-rounded and inflation is
//! added with upward rounding). Nothing narrower is claimed — the
//! output is a candidate SUPERSET for a downstream narrow phase, not a
//! contact result. Output order is deterministic (lexicographic by
//! feature identifier), and the pair count is bounded by an explicit
//! refusal, never truncated silently.

use crate::{ContactInflation, QueryError};
use fs_exec::Cx;

/// Hard bound on features per complex.
pub const MAX_COMPLEX_FEATURES: usize = 1 << 20;

/// Cancellation-poll stride in visited BVH node pairs.
const CHECKPOINT_STRIDE: usize = 1024;

/// One boundary feature of a complex, indexing its vertex array.
/// Edges store `a < b`; faces store their triangle's vertex triple in
/// input order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Feature {
    /// A mesh vertex.
    Vertex(u32),
    /// An undirected edge (canonical `a < b`).
    Edge(u32, u32),
    /// A triangle, by its deterministic input index.
    Face(u32),
}

#[derive(Debug, Clone, Copy)]
struct FeatureBox {
    min: [f64; 3],
    max: [f64; 3],
}

impl FeatureBox {
    fn of_points(points: &[[f64; 3]]) -> FeatureBox {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for p in points {
            for a in 0..3 {
                min[a] = min[a].min(p[a]);
                max[a] = max[a].max(p[a]);
            }
        }
        // Outward rounding keeps the box conservative against the
        // exact feature even after later arithmetic on the bounds.
        FeatureBox {
            min: [min[0].next_down(), min[1].next_down(), min[2].next_down()],
            max: [max[0].next_up(), max[1].next_up(), max[2].next_up()],
        }
    }

    fn inflated(&self, pad: f64) -> FeatureBox {
        FeatureBox {
            min: [
                (self.min[0] - pad).next_down(),
                (self.min[1] - pad).next_down(),
                (self.min[2] - pad).next_down(),
            ],
            max: [
                (self.max[0] + pad).next_up(),
                (self.max[1] + pad).next_up(),
                (self.max[2] + pad).next_up(),
            ],
        }
    }

    fn overlaps(&self, other: &FeatureBox) -> bool {
        (0..3).all(|a| self.min[a] <= other.max[a] && other.min[a] <= self.max[a])
    }

    fn merged(&self, other: &FeatureBox) -> FeatureBox {
        FeatureBox {
            min: [
                self.min[0].min(other.min[0]),
                self.min[1].min(other.min[1]),
                self.min[2].min(other.min[2]),
            ],
            max: [
                self.max[0].max(other.max[0]),
                self.max[1].max(other.max[1]),
                self.max[2].max(other.max[2]),
            ],
        }
    }

    fn center(&self, axis: usize) -> f64 {
        f64::midpoint(self.min[axis], self.max[axis])
    }

    fn widest_axis(&self) -> usize {
        let spans = [
            self.max[0] - self.min[0],
            self.max[1] - self.min[1],
            self.max[2] - self.min[2],
        ];
        let mut axis = 0;
        if spans[1] > spans[axis] {
            axis = 1;
        }
        if spans[2] > spans[axis] {
            axis = 2;
        }
        axis
    }
}

/// The typed vertex/edge/face decomposition of a triangle boundary,
/// with deterministic feature order and per-feature conservative boxes.
#[derive(Debug, Clone)]
pub struct FeatureComplex {
    features: Vec<Feature>,
    boxes: Vec<FeatureBox>,
}

impl FeatureComplex {
    /// Build from vertex positions and triangle index triples.
    ///
    /// Vertices become features in index order, undirected edges are
    /// deduplicated canonically (`a < b`, sorted), and each triangle
    /// becomes a face feature in input order — so the complex is a
    /// pure function of the input arrays.
    ///
    /// # Errors
    /// [`QueryError::InvalidBoundaryIndex`] for an out-of-range or
    /// degenerate (repeated-vertex) triangle;
    /// [`QueryError::InvalidPointSample`] for a non-finite position;
    /// [`QueryError::FeatureComplexTooLarge`] beyond
    /// [`MAX_COMPLEX_FEATURES`] total features.
    pub fn from_triangles(
        positions: &[[f64; 3]],
        triangles: &[[u32; 3]],
    ) -> Result<FeatureComplex, QueryError> {
        for p in positions {
            if !(p[0].is_finite() && p[1].is_finite() && p[2].is_finite()) {
                return Err(QueryError::InvalidPointSample { at: *p });
            }
        }
        let mut edges: Vec<(u32, u32)> = Vec::with_capacity(triangles.len() * 3);
        for (triangle, tri) in triangles.iter().enumerate() {
            for (corner, &index) in tri.iter().enumerate() {
                if index as usize >= positions.len() {
                    return Err(QueryError::InvalidBoundaryIndex {
                        triangle,
                        corner,
                        index,
                        positions: positions.len(),
                    });
                }
            }
            if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                return Err(QueryError::InvalidBoundaryIndex {
                    triangle,
                    corner: 0,
                    index: tri[0],
                    positions: positions.len(),
                });
            }
            for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[0], tri[2])] {
                edges.push((a.min(b), a.max(b)));
            }
        }
        edges.sort_unstable();
        edges.dedup();
        let total = positions
            .len()
            .saturating_add(edges.len())
            .saturating_add(triangles.len());
        if total > MAX_COMPLEX_FEATURES {
            return Err(QueryError::FeatureComplexTooLarge {
                features: total,
                max: MAX_COMPLEX_FEATURES,
            });
        }
        let mut features = Vec::with_capacity(total);
        let mut boxes = Vec::with_capacity(total);
        for (i, p) in positions.iter().enumerate() {
            #[allow(clippy::cast_possible_truncation)]
            features.push(Feature::Vertex(i as u32));
            boxes.push(FeatureBox::of_points(std::slice::from_ref(p)));
        }
        for &(a, b) in &edges {
            features.push(Feature::Edge(a, b));
            boxes.push(FeatureBox::of_points(&[
                positions[a as usize],
                positions[b as usize],
            ]));
        }
        for (t, tri) in triangles.iter().enumerate() {
            #[allow(clippy::cast_possible_truncation)]
            features.push(Feature::Face(t as u32));
            boxes.push(FeatureBox::of_points(&[
                positions[tri[0] as usize],
                positions[tri[1] as usize],
                positions[tri[2] as usize],
            ]));
        }
        Ok(FeatureComplex { features, boxes })
    }

    /// Number of features (vertices + unique edges + faces).
    #[must_use]
    pub fn len(&self) -> usize {
        self.features.len()
    }

    /// Whether the complex has no features.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.features.is_empty()
    }

    /// The feature at a deterministic complex index.
    #[must_use]
    pub fn feature(&self, index: usize) -> Option<Feature> {
        self.features.get(index).copied()
    }
}

/// A deterministic median-split BVH over one complex's inflated boxes.
struct Bvh {
    nodes: Vec<BvhNode>,
    /// Feature indices, permuted so each leaf owns a contiguous range.
    order: Vec<u32>,
    boxes: Vec<FeatureBox>,
}

struct BvhNode {
    bounds: FeatureBox,
    /// Leaf: `(start, len)` into `order`; internal: child node ids.
    kind: BvhKind,
}

enum BvhKind {
    Leaf { start: u32, len: u32 },
    Internal { left: u32, right: u32 },
}

const LEAF_SIZE: usize = 8;

impl Bvh {
    fn build(complex: &FeatureComplex, pad: f64) -> Bvh {
        let boxes: Vec<FeatureBox> = complex.boxes.iter().map(|b| b.inflated(pad)).collect();
        #[allow(clippy::cast_possible_truncation)]
        let mut order: Vec<u32> = (0..boxes.len() as u32).collect();
        let mut nodes = Vec::with_capacity(boxes.len() / LEAF_SIZE * 2 + 1);
        if !order.is_empty() {
            build_node(&boxes, &mut order, 0, &mut nodes);
        }
        Bvh {
            nodes,
            order,
            boxes,
        }
    }
}

fn build_node(boxes: &[FeatureBox], order: &mut [u32], base: u32, nodes: &mut Vec<BvhNode>) -> u32 {
    let mut bounds = boxes[order[0] as usize];
    for &i in order.iter().skip(1) {
        bounds = bounds.merged(&boxes[i as usize]);
    }
    #[allow(clippy::cast_possible_truncation)]
    let id = nodes.len() as u32;
    nodes.push(BvhNode {
        bounds,
        kind: BvhKind::Leaf {
            start: base,
            len: 0,
        },
    });
    if order.len() <= LEAF_SIZE {
        #[allow(clippy::cast_possible_truncation)]
        let len = order.len() as u32;
        nodes[id as usize].kind = BvhKind::Leaf { start: base, len };
        return id;
    }
    // Median split on the widest axis of the node bounds; ties in the
    // comparator fall back to the feature index so the permutation is
    // a pure function of the input.
    let axis = bounds.widest_axis();
    let mid = order.len() / 2;
    order.select_nth_unstable_by(mid, |&a, &b| {
        boxes[a as usize]
            .center(axis)
            .total_cmp(&boxes[b as usize].center(axis))
            .then_with(|| a.cmp(&b))
    });
    let (lo, hi) = order.split_at_mut(mid);
    #[allow(clippy::cast_possible_truncation)]
    let left = build_node(boxes, lo, base, nodes);
    #[allow(clippy::cast_possible_truncation)]
    let right = build_node(boxes, hi, base + mid as u32, nodes);
    nodes[id as usize].kind = BvhKind::Internal { left, right };
    id
}

/// Conservative CCD candidate enumeration between two complexes.
///
/// Each side's feature boxes are inflated by its own declared motion
/// bound (a certified radius the caller derives from velocities and
/// the CCD window); a feature pair is emitted exactly when the
/// inflated boxes overlap. The result is sorted lexicographically by
/// `(index_a, index_b)` and is therefore a pure function of the
/// inputs.
///
/// # Errors
/// [`QueryError::FeatureInvalidInflation`] for non-finite or negative
/// motion bounds; [`QueryError::FeatureTooManyPairs`] when the
/// candidate count exceeds `max_pairs` (refusal, not truncation);
/// [`QueryError::Cancelled`] on cancellation.
pub fn ccd_candidates(
    a: &FeatureComplex,
    b: &FeatureComplex,
    motion_a: f64,
    motion_b: f64,
    max_pairs: usize,
    cx: &Cx<'_>,
) -> Result<Vec<(usize, usize)>, QueryError> {
    for motion in [motion_a, motion_b] {
        if !motion.is_finite() || motion < 0.0 {
            return Err(QueryError::FeatureInvalidInflation {
                inflation_bits: motion.to_bits(),
            });
        }
    }
    if a.is_empty() || b.is_empty() {
        return Ok(Vec::new());
    }
    let bvh_a = Bvh::build(a, motion_a);
    let bvh_b = Bvh::build(b, motion_b);
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    let mut stack: Vec<(u32, u32)> = vec![(0, 0)];
    let mut visited = 0usize;
    while let Some((na, nb)) = stack.pop() {
        if visited.is_multiple_of(CHECKPOINT_STRIDE) && cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        visited += 1;
        let node_a = &bvh_a.nodes[na as usize];
        let node_b = &bvh_b.nodes[nb as usize];
        if !node_a.bounds.overlaps(&node_b.bounds) {
            continue;
        }
        match (&node_a.kind, &node_b.kind) {
            (&BvhKind::Leaf { start: sa, len: la }, &BvhKind::Leaf { start: sb, len: lb }) => {
                for &fa in &bvh_a.order[sa as usize..(sa + la) as usize] {
                    for &fb in &bvh_b.order[sb as usize..(sb + lb) as usize] {
                        if bvh_a.boxes[fa as usize].overlaps(&bvh_b.boxes[fb as usize]) {
                            if pairs.len() >= max_pairs {
                                return Err(QueryError::FeatureTooManyPairs { max: max_pairs });
                            }
                            pairs.push((fa as usize, fb as usize));
                        }
                    }
                }
            }
            (&BvhKind::Internal { left, right }, _) => {
                stack.push((left, nb));
                stack.push((right, nb));
            }
            (&BvhKind::Leaf { .. }, &BvhKind::Internal { left, right }) => {
                stack.push((na, left));
                stack.push((na, right));
            }
        }
    }
    pairs.sort_unstable();
    Ok(pairs)
}

/// Conservative CCD candidates with each side's representation/motion error
/// added to its declared swept radius.
///
/// The added radius is proof-bearing and rounded outward before BVH
/// construction. Increasing either inflation can therefore only retain or add
/// candidate pairs. Exact-zero inflations delegate with the original motion
/// bits unchanged.
///
/// # Errors
/// The refusals from [`ccd_candidates`], plus
/// [`QueryError::InvalidContactInflation`] if a finite motion radius cannot be
/// inflated without overflow.
pub fn ccd_candidates_with_inflation(
    a: &FeatureComplex,
    b: &FeatureComplex,
    motion_a: f64,
    motion_b: f64,
    max_pairs: usize,
    cx: &Cx<'_>,
    inflation_a: ContactInflation,
    inflation_b: ContactInflation,
) -> Result<Vec<(usize, usize)>, QueryError> {
    if inflation_a.radius() == 0.0 && inflation_b.radius() == 0.0 {
        return ccd_candidates(a, b, motion_a, motion_b, max_pairs, cx);
    }
    for motion in [motion_a, motion_b] {
        if !motion.is_finite() || motion < 0.0 {
            return Err(QueryError::FeatureInvalidInflation {
                inflation_bits: motion.to_bits(),
            });
        }
    }
    let motion_a = inflation_a.inflate_nonnegative(motion_a)?;
    let motion_b = inflation_b.inflate_nonnegative(motion_b)?;
    ccd_candidates(a, b, motion_a, motion_b, max_pairs, cx)
}
