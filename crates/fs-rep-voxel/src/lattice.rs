//! Explicit lattice/strut graphs: nodes + struts with cross-section
//! attributes, FrankenNetworkx round-trip (ground structures, infill
//! graphs), and REALIZATION — strut graph → smooth-min capsule SDF, a
//! watertight solid by construction (a continuous level set has no open
//! boundary). Degenerate graphs (coincident nodes, zero-length struts)
//! are structured refusals.

use crate::VoxelError;
use fnx_classes::Graph;
use fnx_runtime::{CgseValue, CompatibilityMode};
use fs_math::det;

/// One lattice node.
#[derive(Debug, Clone, PartialEq)]
pub struct LatticeNode {
    /// Position (m).
    pub pos: [f64; 3],
    /// Junction blend radius (m).
    pub radius: f64,
}

/// One strut connecting two node indices.
#[derive(Debug, Clone, PartialEq)]
pub struct Strut {
    /// First node index.
    pub a: usize,
    /// Second node index.
    pub b: usize,
    /// Cross-section radius (m).
    pub radius: f64,
}

/// A validated lattice graph.
#[derive(Debug, Clone, PartialEq)]
pub struct LatticeGraph {
    /// Nodes.
    pub nodes: Vec<LatticeNode>,
    /// Struts.
    pub struts: Vec<Strut>,
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl LatticeGraph {
    /// Validate and construct. Coincident nodes (within `1e-12` m),
    /// zero-length struts, out-of-range indices, and non-positive radii
    /// are STRUCTURED errors, not downstream NaNs.
    ///
    /// # Errors
    /// [`VoxelError::Lattice`] with the offending element named.
    pub fn new(nodes: Vec<LatticeNode>, struts: Vec<Strut>) -> Result<Self, VoxelError> {
        for (i, n) in nodes.iter().enumerate() {
            if n.pos.iter().any(|v| !v.is_finite()) || !(n.radius > 0.0) {
                return Err(VoxelError::Lattice {
                    what: format!("node {i}: non-finite position or non-positive radius"),
                });
            }
        }
        for i in 0..nodes.len() {
            for j in (i + 1)..nodes.len() {
                let d = sub(nodes[i].pos, nodes[j].pos);
                if dot(d, d) < 1e-24 {
                    return Err(VoxelError::Lattice {
                        what: format!("nodes {i} and {j} are coincident"),
                    });
                }
            }
        }
        for (k, s) in struts.iter().enumerate() {
            if s.a >= nodes.len() || s.b >= nodes.len() {
                return Err(VoxelError::Lattice {
                    what: format!("strut {k} references a missing node"),
                });
            }
            if s.a == s.b {
                return Err(VoxelError::Lattice {
                    what: format!("strut {k} is zero-length (both ends at node {})", s.a),
                });
            }
            if !(s.radius > 0.0 && s.radius.is_finite()) {
                return Err(VoxelError::Lattice {
                    what: format!("strut {k}: non-positive radius"),
                });
            }
        }
        Ok(LatticeGraph { nodes, struts })
    }

    /// Export to a FrankenNetworkx graph: node ids `"n<i>"` with
    /// position/radius attributes; edges carry the cross-section radius.
    #[must_use]
    pub fn to_fnx(&self) -> Graph {
        let mut g = Graph::new(CompatibilityMode::Strict);
        for (i, n) in self.nodes.iter().enumerate() {
            let mut attrs = fnx_classes::AttrMap::new();
            attrs.insert("x".to_string(), CgseValue::Float(n.pos[0]));
            attrs.insert("y".to_string(), CgseValue::Float(n.pos[1]));
            attrs.insert("z".to_string(), CgseValue::Float(n.pos[2]));
            attrs.insert("radius".to_string(), CgseValue::Float(n.radius));
            g.add_node_with_attrs(format!("n{i}"), attrs);
        }
        for s in &self.struts {
            let mut attrs = fnx_classes::AttrMap::new();
            attrs.insert("radius".to_string(), CgseValue::Float(s.radius));
            let _ = g.add_edge_with_attrs(format!("n{}", s.a), format!("n{}", s.b), attrs);
        }
        g
    }

    /// Import back from FrankenNetworkx (inverse of [`Self::to_fnx`]).
    ///
    /// # Errors
    /// [`VoxelError::Graph`] on missing/mistyped attributes;
    /// [`VoxelError::Lattice`] if the recovered lattice is degenerate.
    pub fn from_fnx(g: &Graph) -> Result<Self, VoxelError> {
        let float_attr = |attrs: &fnx_classes::AttrMap, key: &str| -> Result<f64, VoxelError> {
            match attrs.get(key) {
                Some(CgseValue::Float(v)) => Ok(*v),
                Some(CgseValue::Int(v)) => {
                    #[allow(clippy::cast_precision_loss)]
                    Ok(*v as f64)
                }
                _ => Err(VoxelError::Graph {
                    what: format!("missing/mistyped float attribute {key:?}"),
                }),
            }
        };
        let mut ids: Vec<String> = g.nodes_ordered().iter().map(|s| (*s).to_string()).collect();
        ids.sort_by_key(|id| {
            id.strip_prefix('n')
                .and_then(|t| t.parse::<usize>().ok())
                .unwrap_or(usize::MAX)
        });
        let mut nodes = Vec::with_capacity(ids.len());
        for id in &ids {
            let attrs = g.node_attrs(id).ok_or_else(|| VoxelError::Graph {
                what: format!("node {id} has no attributes"),
            })?;
            nodes.push(LatticeNode {
                pos: [
                    float_attr(attrs, "x")?,
                    float_attr(attrs, "y")?,
                    float_attr(attrs, "z")?,
                ],
                radius: float_attr(attrs, "radius")?,
            });
        }
        let index_of = |id: &str| -> Result<usize, VoxelError> {
            ids.iter().position(|x| x == id).ok_or_else(|| VoxelError::Graph {
                what: format!("edge references unknown node {id}"),
            })
        };
        let mut struts = Vec::new();
        for (i, id) in ids.iter().enumerate() {
            if let Some(neighbors) = g.neighbors(id) {
                for other in neighbors {
                    let j = index_of(other)?;
                    if j <= i {
                        continue; // undirected: keep one direction
                    }
                    let attrs = g.edge_attrs(id, other).ok_or_else(|| VoxelError::Graph {
                        what: format!("edge {id}-{other} has no attributes"),
                    })?;
                    struts.push(Strut {
                        a: i,
                        b: j,
                        radius: float_attr(attrs, "radius")?,
                    });
                }
            }
        }
        LatticeGraph::new(nodes, struts)
    }

    /// The realized solid's signed distance at a point: smooth-min of
    /// capsule SDFs (junction blending radius = the larger node blend).
    /// Negative inside — the level set of a continuous function, hence
    /// watertight by construction.
    #[must_use]
    pub fn sdf(&self, p: [f64; 3]) -> f64 {
        let mut d = f64::INFINITY;
        for s in &self.struts {
            let cap = capsule_distance(p, self.nodes[s.a].pos, self.nodes[s.b].pos, s.radius);
            let blend = self.nodes[s.a].radius.max(self.nodes[s.b].radius);
            d = smooth_min(d, cap, blend);
        }
        d
    }

    /// Per-strut realization receipts (ledger logs): JSON lines with
    /// endpoints, length, and radius.
    #[must_use]
    pub fn realization_receipts(&self) -> Vec<String> {
        self.struts
            .iter()
            .enumerate()
            .map(|(k, s)| {
                let d = sub(self.nodes[s.b].pos, self.nodes[s.a].pos);
                format!(
                    "{{\"kind\":\"strut-realization\",\"strut\":{k},\"a\":{},\"b\":{},\
                     \"length\":{},\"radius\":{}}}",
                    s.a,
                    s.b,
                    det::sqrt(dot(d, d)),
                    s.radius
                )
            })
            .collect()
    }
}

/// Exact distance to a capsule (segment with radius).
fn capsule_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3], radius: f64) -> f64 {
    let ab = sub(b, a);
    let ap = sub(p, a);
    let t = (dot(ap, ab) / dot(ab, ab)).clamp(0.0, 1.0);
    let closest = [a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2]];
    let d = sub(p, closest);
    det::sqrt(dot(d, d)) - radius
}

/// Polynomial smooth minimum (junction blending); reduces to `min` as
/// `k → 0`, never exceeds either operand by more than the blend.
fn smooth_min(a: f64, b: f64, k: f64) -> f64 {
    if !a.is_finite() {
        return b;
    }
    let h = (k - (a - b).abs()).max(0.0) / k.max(f64::MIN_POSITIVE);
    a.min(b) - h * h * k * 0.25
}
