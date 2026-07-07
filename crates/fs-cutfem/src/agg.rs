//! Aggregated-element fallback (belt + suspenders with ghost penalty):
//! a cell whose inside fraction is pathologically small cannot by
//! itself support its DOFs — the stiffness it contributes is O(cut
//! fraction) and the conditioning collapses as the cut degenerates.
//! Aggregation constrains every node supported ONLY by small-cut cells
//! to the Q1 polynomial EXTENSION from a nearby well-cut anchor cell
//! (bilinear extrapolation of the anchor's nodal values), removing the
//! degenerate DOFs from the free set entirely.
//!
//! Policy (documented, conformance-logged): a cut cell is SMALL below
//! `small_fraction`, a valid anchor is Inside or a cut cell at or
//! above `good_fraction`; anchors are found by deterministic
//! breadth-first search over face-adjacent active cells (≤ 4 rings);
//! failure to find one is a structured refusal
//! ([`crate::CutFemError::AggregationNoAnchor`]), never a silent
//! degradation.

use crate::CutFemError;
use crate::fem::CellClass;
use crate::grid::{CellKey, NodeKey, Quadtree};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

/// Aggregation thresholds (inside-area fractions in [0, 1]).
#[derive(Debug, Clone, Copy)]
pub struct AggPolicy {
    /// Cells strictly below this fraction are "small" (default 0.1).
    pub small_fraction: f64,
    /// Anchor cells need at least this fraction (default 0.5).
    pub good_fraction: f64,
}

impl Default for AggPolicy {
    fn default() -> Self {
        AggPolicy {
            small_fraction: 0.1,
            good_fraction: 0.5,
        }
    }
}

/// The computed aggregation: extrapolation constraints plus the
/// policy log rows (one JSON object per aggregated node).
pub(crate) struct AggOutcome {
    pub constraints: Vec<(NodeKey, Vec<(NodeKey, f64)>)>,
    pub log: Vec<String>,
}

pub(crate) fn aggregate(
    grid: &Quadtree,
    class: &BTreeMap<CellKey, CellClass>,
    frac: &BTreeMap<CellKey, f64>,
    active: &BTreeSet<CellKey>,
    already_constrained: &BTreeSet<NodeKey>,
    policy: AggPolicy,
) -> Result<AggOutcome, CutFemError> {
    let small = |c: &CellKey| {
        class.get(c) == Some(&CellClass::Cut)
            && frac.get(c).copied().unwrap_or(0.0) < policy.small_fraction
    };
    let good = |c: &CellKey| {
        class.get(c) == Some(&CellClass::Inside)
            || frac.get(c).copied().unwrap_or(0.0) >= policy.good_fraction
    };

    let mut node_leaves: BTreeMap<NodeKey, Vec<CellKey>> = BTreeMap::new();
    for &c in active {
        for n in grid.corner_nodes(c) {
            node_leaves.entry(n).or_default().push(c);
        }
    }

    let mut out = AggOutcome {
        constraints: Vec::new(),
        log: Vec::new(),
    };
    for (&node, leaves) in &node_leaves {
        if already_constrained.contains(&node) || !leaves.iter().all(small) {
            continue;
        }
        let anchor = find_anchor(grid, active, leaves, &good)
            .ok_or(CutFemError::AggregationNoAnchor { node })?;
        let (lo, hi) = grid.rect(anchor);
        let p = grid.node_pos(node);
        let xi = (p[0] - lo[0]) / (hi[0] - lo[0]);
        let et = (p[1] - lo[1]) / (hi[1] - lo[1]);
        let w = [
            (1.0 - xi) * (1.0 - et),
            xi * (1.0 - et),
            xi * et,
            (1.0 - xi) * et,
        ];
        let corners = grid.corner_nodes(anchor);
        let terms: Vec<(NodeKey, f64)> = corners.iter().copied().zip(w).collect();
        let mut row = String::new();
        let _ = write!(
            row,
            "{{\"node\":[{},{}],\"anchor\":[{},{},{}],\"xi\":{xi:.4},\"eta\":{et:.4}}}",
            node.0, node.1, anchor.0, anchor.1, anchor.2
        );
        out.log.push(row);
        out.constraints.push((node, terms));
    }
    Ok(out)
}

/// Deterministic ring-by-ring BFS over face-adjacent active cells.
fn find_anchor(
    grid: &Quadtree,
    active: &BTreeSet<CellKey>,
    start: &[CellKey],
    good: &dyn Fn(&CellKey) -> bool,
) -> Option<CellKey> {
    let mut visited: BTreeSet<CellKey> = start.iter().copied().collect();
    let mut ring: BTreeSet<CellKey> = visited.clone();
    for _ in 0..4 {
        let mut next: BTreeSet<CellKey> = BTreeSet::new();
        for &c in &ring {
            for dir in 0..4u8 {
                if let Some(nb) = grid.covering_neighbor(c, dir) {
                    if active.contains(&nb) && !visited.contains(&nb) {
                        next.insert(nb);
                    }
                }
            }
        }
        for &c in &next {
            if good(&c) {
                return Some(c);
            }
        }
        visited.extend(next.iter().copied());
        if next.is_empty() {
            return None;
        }
        ring = next;
    }
    None
}
