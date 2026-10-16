//! Deterministic nested dissection (George1973): interiors before separators.
//! Graph path: index-ordered heavy-edge coarsening, balanced bisection, FM
//! refinement on every uncoarsening level, then an explicit vertex cover of
//! the cut. Optional coordinates select geometric bisection instead. AMD
//! orders leaves and separators; no second symbolic/factorization owner.

use crate::{
    Coo, Csr,
    direct::{LdltError, amd_order},
};
use std::collections::{BTreeMap, BinaryHeap, VecDeque};

/// Copyable ordering policy; coordinates are borrowed separately so choosing
/// ND does not make every existing factor/slice option own a mesh copy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NdOptions {
    /// Maximum subproblem size directly handed to AMD; must be nonzero.
    pub leaf_size: usize,
    /// FM passes per coarsening level, in 0..=64.
    pub refinement_passes: usize,
}

/// Original vertex ids. Every branch certifies that there is no edge between
/// its left and right descendants after removing `separator`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SeparatorTree {
    /// AMD-ordered leaf (or a subproblem without a useful split).
    Leaf(Vec<usize>),
    /// Interiors are eliminated before all separator vertices.
    Branch {
        /// First interior subtree.
        left: Box<SeparatorTree>,
        /// Second interior subtree, with no edges to the first.
        right: Box<SeparatorTree>,
        /// Original ids covering every edge of the cut.
        separator: Vec<usize>,
    },
}

#[derive(Clone, Debug, PartialEq, Eq)]
/// Permutation plus the explicit disconnection witnesses.
pub struct NdOrdering {
    /// perm[k] is the original vertex eliminated kth.
    pub permutation: Vec<usize>,
    /// Hierarchical separator witness in original vertex ids.
    pub tree: SeparatorTree,
}

type Graph = Vec<Vec<(usize, usize)>>;

/// Topological ND, with no inferred geometric coordinates.
pub fn nd_order(a: &Csr, opts: &NdOptions) -> Result<NdOrdering, LdltError> {
    nd_order_with_coordinates(a, opts, None)
}

/// Geometric ND when coordinates are supplied, multilevel graph ND otherwise.
/// Stored entries, including stored zeros, are structural edges. The input
/// must be square and structurally symmetric, exactly like the LDLT owner.
pub fn nd_order_with_coordinates(
    a: &Csr,
    opts: &NdOptions,
    coordinates: Option<&[[f64; 3]]>,
) -> Result<NdOrdering, LdltError> {
    if a.nrows() != a.ncols() {
        return Err(LdltError::NotSquare {
            nrows: a.nrows(),
            ncols: a.ncols(),
        });
    }
    if opts.leaf_size == 0 || opts.refinement_passes > 64 {
        return Err(LdltError::InvalidOptions {
            field: "ND leaf_size>0 and refinement_passes<=64",
        });
    }
    if let Some(x) = coordinates {
        if x.len() != a.nrows() || x.iter().flatten().any(|v| !v.is_finite()) {
            return Err(LdltError::InvalidOptions {
                field: "ND coordinates: one finite xyz per row",
            });
        }
    }
    let mut graph = vec![Vec::new(); a.nrows()];
    for (r, neighbors) in graph.iter_mut().enumerate() {
        for &c in a.row(r).0 {
            if a.row(c).0.binary_search(&r).is_err() {
                return Err(LdltError::StructurallyAsymmetric { row: r, col: c });
            }
            if c != r {
                neighbors.push((c, 1));
            }
        }
    }
    let vertices: Vec<_> = (0..graph.len()).collect();
    if graph.iter().all(Vec::is_empty) {
        // AMD on an edgeless graph is exactly ascending ids. Avoid recursively
        // partitioning a million-entry diagonal mass matrix for no fill gain.
        return Ok(NdOrdering {
            permutation: vertices.clone(),
            tree: SeparatorTree::Leaf(vertices),
        });
    }
    let mut permutation = Vec::with_capacity(graph.len());
    let tree = dissect(&graph, &vertices, coordinates, opts, &mut permutation);
    Ok(NdOrdering { permutation, tree })
}

fn local_graph(graph: &Graph, vertices: &[usize]) -> Graph {
    vertices
        .iter()
        .map(|&v| {
            graph[v]
                .iter()
                .filter_map(|&(u, w)| vertices.binary_search(&u).ok().map(|i| (i, w)))
                .collect()
        })
        .collect()
}

fn amd_vertices(graph: &Graph, vertices: &[usize], out: &mut Vec<usize>) {
    let mut coo = Coo::new(vertices.len(), vertices.len());
    for (i, &v) in vertices.iter().enumerate() {
        coo.push(i, i, 1.0);
        for &(u, _) in &graph[v] {
            if let Ok(j) = vertices.binary_search(&u) {
                coo.push(i, j, 1.0);
            }
        }
    }
    out.extend(amd_order(&coo.assemble()).into_iter().map(|i| vertices[i]));
}

fn dissect(
    graph: &Graph,
    vertices: &[usize],
    xyz: Option<&[[f64; 3]]>,
    opts: &NdOptions,
    out: &mut Vec<usize>,
) -> SeparatorTree {
    if vertices.len() <= opts.leaf_size {
        amd_vertices(graph, vertices, out);
        return SeparatorTree::Leaf(vertices.to_vec());
    }
    let local = local_graph(graph, vertices);
    let side = if let Some(xyz) = xyz {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for &v in vertices {
            for a in 0..3 {
                lo[a] = lo[a].min(xyz[v][a]);
                hi[a] = hi[a].max(xyz[v][a]);
            }
        }
        let axis = (0..3)
            .max_by(|&a, &b| (hi[a] - lo[a]).total_cmp(&(hi[b] - lo[b])).then(b.cmp(&a)))
            .unwrap_or(0);
        let mut order: Vec<_> = (0..vertices.len()).collect();
        order.sort_by(|&a, &b| {
            xyz[vertices[a]][axis]
                .total_cmp(&xyz[vertices[b]][axis])
                .then(vertices[a].cmp(&vertices[b]))
        });
        let pivot = xyz[vertices[order[(vertices.len() - 1) / 2]]][axis];
        // Keep coincident coordinate planes together (including the three
        // displacement components of a mesh vertex). Splitting a median
        // plane by row id needlessly turns a flat grid separator into stairs.
        let side = vertices.iter().map(|&v| xyz[v][axis] > pivot).collect();
        side
    } else {
        multilevel(&local, opts.refinement_passes)
    };
    // Greedy vertex cover of the bipartite cut: every crossing edge has a
    // separator endpoint. Highest live cut degree, then lowest original id.
    let mut degree: Vec<_> = local
        .iter()
        .enumerate()
        .map(|(v, adj)| adj.iter().filter(|&&(u, _)| side[u] != side[v]).count())
        .collect();
    let mut heap: BinaryHeap<_> = degree
        .iter()
        .enumerate()
        .filter(|&(_, d)| *d > 0)
        .map(|(i, &d)| (d, core::cmp::Reverse(i)))
        .collect();
    let mut separator = vec![false; vertices.len()];
    while let Some((d, core::cmp::Reverse(v))) = heap.pop() {
        if separator[v] || degree[v] != d {
            continue;
        }
        separator[v] = true;
        for &(u, _) in &local[v] {
            if !separator[u] && side[u] != side[v] {
                degree[u] -= 1;
                if degree[u] > 0 {
                    heap.push((degree[u], core::cmp::Reverse(u)));
                }
            }
        }
        degree[v] = 0;
    }
    let mut left = Vec::new();
    let mut right = Vec::new();
    let mut sep = Vec::new();
    for (i, &v) in vertices.iter().enumerate() {
        if separator[i] {
            sep.push(v);
        } else if side[i] {
            right.push(v);
        } else {
            left.push(v);
        }
    }
    // Complete/expander graphs may have no useful balanced vertex separator.
    // AMD the whole subproblem; never claim an empty half separated anything.
    if left.is_empty() || right.is_empty() {
        amd_vertices(graph, vertices, out);
        return SeparatorTree::Leaf(vertices.to_vec());
    }
    let left = Box::new(dissect(graph, &left, xyz, opts, out));
    let right = Box::new(dissect(graph, &right, xyz, opts, out));
    amd_vertices(graph, &sep, out);
    SeparatorTree::Branch {
        left,
        right,
        separator: sep,
    }
}

fn multilevel(graph: &Graph, passes: usize) -> Vec<bool> {
    let mut levels = vec![(graph.clone(), vec![1usize; graph.len()])];
    let mut maps = Vec::new();
    while levels.last().map_or(0, |x| x.0.len()) > 64 {
        let (g, weights) = levels.last().expect("root level exists");
        let mut map = vec![usize::MAX; g.len()];
        let mut cw = Vec::new();
        for v in 0..g.len() {
            if map[v] != usize::MAX {
                continue;
            }
            let partner = g[v]
                .iter()
                .filter(|&&(u, _)| map[u] == usize::MAX)
                .max_by(|&&(u, w), &&(v, z)| w.cmp(&z).then(v.cmp(&u)))
                .map(|&(u, _)| u);
            let c = cw.len();
            map[v] = c;
            let mut weight = weights[v];
            if let Some(u) = partner {
                map[u] = c;
                weight += weights[u];
            }
            cw.push(weight);
        }
        if cw.len() * 4 > g.len() * 3 {
            break;
        }
        let mut edges = vec![BTreeMap::<usize, usize>::new(); cw.len()];
        for v in 0..g.len() {
            for &(u, w) in &g[v] {
                if map[v] != map[u] {
                    *edges[map[v]].entry(map[u]).or_default() += w;
                }
            }
        }
        let coarse = edges.into_iter().map(|e| e.into_iter().collect()).collect();
        maps.push(map);
        levels.push((coarse, cw));
    }
    let (g, w) = levels.last().expect("root level exists");
    let mut side = initial_partition(g, w);
    refine(g, w, &mut side, passes);
    while let Some(map) = maps.pop() {
        levels.pop();
        let (g, w) = levels.last().expect("finer level exists");
        side = map.into_iter().map(|v| side[v]).collect();
        refine(g, w, &mut side, passes);
    }
    side
}

fn distances(graph: &Graph, start: usize) -> Vec<usize> {
    let mut d = vec![usize::MAX; graph.len()];
    let mut q = VecDeque::new();
    d[start] = 0;
    q.push_back(start);
    while let Some(v) = q.pop_front() {
        for &(u, _) in &graph[v] {
            if d[u] == usize::MAX {
                d[u] = d[v] + 1;
                q.push_back(u);
            }
        }
    }
    d
}
fn initial_partition(graph: &Graph, weights: &[usize]) -> Vec<bool> {
    if graph.is_empty() {
        return Vec::new();
    }
    let d = distances(graph, 0);
    let seed = (0..graph.len())
        .filter(|&v| d[v] != usize::MAX)
        .max_by_key(|&v| (d[v], core::cmp::Reverse(v)))
        .unwrap_or(0);
    let d = distances(graph, seed);
    let mut order: Vec<_> = (0..graph.len()).collect();
    order.sort_by_key(|&v| (d[v], v));
    let total: usize = weights.iter().sum();
    let mut sum = 0;
    let mut side = vec![true; graph.len()];
    for v in order {
        if sum >= total / 2 {
            break;
        }
        side[v] = false;
        sum += weights[v];
    }
    side
}

// Fiduccia--Mattheyses edge-cut passes: locked single moves, incremental
// gains, bounded imbalance, best-prefix rollback. No floating-point scores.
fn refine(graph: &Graph, weights: &[usize], side: &mut Vec<bool>, passes: usize) {
    let total: usize = weights.iter().sum();
    let cap = total.saturating_mul(3) / 5 + weights.iter().copied().max().unwrap_or(0);
    for _ in 0..passes {
        let initial = side.clone();
        let mut moves = Vec::new();
        let mut best_len = 0;
        let mut best_gain = 0isize;
        let mut cumulative = 0isize;
        let mut locked = vec![false; graph.len()];
        let mut gain: Vec<isize> = graph
            .iter()
            .enumerate()
            .map(|(v, adj)| {
                adj.iter()
                    .map(|&(u, w)| {
                        if side[u] != side[v] {
                            w as isize
                        } else {
                            -(w as isize)
                        }
                    })
                    .sum()
            })
            .collect();
        let mut heap: BinaryHeap<_> = gain
            .iter()
            .enumerate()
            .map(|(v, &g)| (g, core::cmp::Reverse(v)))
            .collect();
        let mut sizes = [0usize; 2];
        for v in 0..graph.len() {
            sizes[usize::from(side[v])] += weights[v];
        }
        while let Some((g, core::cmp::Reverse(v))) = heap.pop() {
            if locked[v] || g != gain[v] {
                continue;
            }
            let from = usize::from(side[v]);
            let to = 1 - from;
            if sizes[to] + weights[v] > cap || sizes[from] <= weights[v] {
                continue;
            }
            locked[v] = true;
            sizes[from] -= weights[v];
            sizes[to] += weights[v];
            side[v] = !side[v];
            cumulative += g;
            moves.push(v);
            if cumulative > best_gain {
                best_gain = cumulative;
                best_len = moves.len();
            }
            for &(u, w) in &graph[v] {
                if !locked[u] {
                    gain[u] += if side[u] == side[v] {
                        -(2 * w as isize)
                    } else {
                        2 * w as isize
                    };
                    heap.push((gain[u], core::cmp::Reverse(u)));
                }
            }
        }
        *side = initial;
        for &v in &moves[..best_len] {
            side[v] = !side[v];
        }
        if best_gain == 0 {
            break;
        }
    }
}
