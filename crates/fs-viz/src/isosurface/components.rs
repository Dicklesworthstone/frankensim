//! Connected strict-sublevel support on the native marching-tetrahedra grid.

use super::{Grid3, IsoSurfaceError, TETRAHEDRA3};

/// Node labels of the piecewise-linear field's strict-sublevel components.
///
/// Inside nodes are joined only through edges of the original six tetrahedra
/// per cube. Nodes exactly on the level are excluded. Labels are contiguous,
/// ordered by each component's first x-fastest node, and repeatable for an
/// unchanged field. They are snapshot labels, not tracked physical identities.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SublevelComponents3 {
    labels: Vec<usize>,
    count: usize,
}

impl SublevelComponents3 {
    /// Number of nonempty components.
    #[must_use]
    pub const fn component_count(&self) -> usize {
        self.count
    }

    /// Number of nodes in the original grid, including outside nodes.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.labels.len()
    }

    /// Component of an x-fastest node, or none for outside/out-of-range nodes.
    #[must_use]
    pub fn component_of_node(&self, index: usize) -> Option<usize> {
        self.labels
            .get(index)
            .copied()
            .filter(|&label| label != usize::MAX)
    }
}

impl Grid3 {
    /// Label connected components of the same strict-sublevel field used by
    /// the native volume query and mesher. No surface extraction is required.
    ///
    /// The admitted grid count bounds storage and work. Union by rank with
    /// path compression uses one node-sized parent array plus one byte per
    /// node; the parent array becomes the returned labels. No adjacency graph
    /// or per-edge allocation is built. The final count must fit the explicit
    /// positive component limit, otherwise no partial labels are returned.
    pub fn sublevel_components(
        &self,
        iso: f64,
        component_limit: usize,
    ) -> Result<SublevelComponents3, IsoSurfaceError> {
        self.sublevel_components_with_control(iso, component_limit, || true)
    }

    /// Cancellable support labeling, polled per node/cube and before return.
    pub fn sublevel_components_with_control(
        &self,
        iso: f64,
        component_limit: usize,
        mut control: impl FnMut() -> bool,
    ) -> Result<SublevelComponents3, IsoSurfaceError> {
        if !iso.is_finite() {
            return Err(IsoSurfaceError::NonFiniteIso { iso });
        }
        if component_limit == 0 {
            return Err(IsoSurfaceError::ZeroComponentLimit);
        }
        if !control() {
            return Err(IsoSurfaceError::Cancelled);
        }
        let nodes = self.values.len();
        let mut parents = Vec::new();
        let mut ranks = Vec::new();
        parents
            .try_reserve_exact(nodes)
            .map_err(|_| IsoSurfaceError::ComponentAllocationFailed { nodes })?;
        ranks
            .try_reserve_exact(nodes)
            .map_err(|_| IsoSurfaceError::ComponentAllocationFailed { nodes })?;
        for (node, &value) in self.values.iter().enumerate() {
            if !control() {
                return Err(IsoSurfaceError::Cancelled);
            }
            parents.push(if value < iso { node } else { usize::MAX });
            ranks.push(0u8);
        }
        let [nx, ny, nz] = self.dimensions;
        for k in 0..(nz - 1) {
            for j in 0..(ny - 1) {
                for i in 0..(nx - 1) {
                    if !control() {
                        return Err(IsoSurfaceError::Cancelled);
                    }
                    let cube = self.cube_nodes(i, j, k);
                    for tetrahedron in TETRAHEDRA3 {
                        let mut first = None;
                        for local in tetrahedron {
                            let node = cube[local];
                            if parents[node] == usize::MAX {
                                continue;
                            }
                            if let Some(other) = first {
                                join(&mut parents, &mut ranks, other, node);
                            } else {
                                first = Some(node);
                            }
                        }
                    }
                }
            }
        }
        // Flatten before reusing root entries as encoded component labels.
        for node in 0..nodes {
            if !control() {
                return Err(IsoSurfaceError::Cancelled);
            }
            if parents[node] != usize::MAX {
                parents[node] = root(&mut parents, node);
            }
        }
        drop(ranks);
        let mut count = 0;
        for node in 0..nodes {
            if !control() {
                return Err(IsoSurfaceError::Cancelled);
            }
            let parent = parents[node];
            if parent == usize::MAX || parent >= nodes {
                continue;
            }
            if parents[parent] < nodes {
                if count == component_limit {
                    return Err(IsoSurfaceError::ComponentBudgetExceeded {
                        limit: component_limit,
                    });
                }
                parents[parent] = nodes
                    .checked_add(count)
                    .filter(|&label| label < usize::MAX)
                    .ok_or(IsoSurfaceError::ComponentIndexOverflow)?;
                count += 1;
            }
            parents[node] = parents[parent];
        }
        for label in &mut parents {
            if !control() {
                return Err(IsoSurfaceError::Cancelled);
            }
            if *label != usize::MAX {
                *label -= nodes;
            }
        }
        if !control() {
            return Err(IsoSurfaceError::Cancelled);
        }
        Ok(SublevelComponents3 {
            labels: parents,
            count,
        })
    }
}

fn root(parents: &mut [usize], mut node: usize) -> usize {
    let mut ancestor = node;
    while parents[ancestor] != ancestor {
        ancestor = parents[ancestor];
    }
    while parents[node] != node {
        let next = parents[node];
        parents[node] = ancestor;
        node = next;
    }
    ancestor
}

fn join(parents: &mut [usize], ranks: &mut [u8], first: usize, second: usize) {
    let (mut a, mut b) = (root(parents, first), root(parents, second));
    if a == b {
        return;
    }
    if ranks[a] < ranks[b] || (ranks[a] == ranks[b] && a > b) {
        core::mem::swap(&mut a, &mut b);
    }
    parents[b] = a;
    if ranks[a] == ranks[b] {
        ranks[a] += 1;
    }
}
