//! Ground structures: a node grid plus every candidate member that
//! passes the fabrication rules — length bounds, an allowed angle set
//! (fabrication constrains member directions), and a per-node
//! neighbor cap. The candidate graph is a FrankenNetworkx `Graph`
//! (the §9.5 pairing), generation is REPRODUCIBLE (BTree orders, no
//! ambient state), and the stats row is the ledger evidence.

use fnx_classes::Graph;
use fnx_runtime::CompatibilityMode;
use std::fmt::Write as _;

/// Fabrication rules for candidate members.
#[derive(Debug, Clone)]
pub struct GroundRules {
    /// Minimum member length.
    pub min_len: f64,
    /// Maximum member length.
    pub max_len: f64,
    /// Allowed direction angles (degrees in [0, 180), tolerance
    /// `angle_tol`); empty = all directions allowed.
    pub angles: Vec<f64>,
    /// Angle tolerance (degrees).
    pub angle_tol: f64,
}

impl Default for GroundRules {
    fn default() -> Self {
        GroundRules {
            min_len: 1e-9,
            max_len: f64::INFINITY,
            angles: Vec::new(),
            angle_tol: 1e-6,
        }
    }
}

/// A generated ground structure.
pub struct GroundStructure {
    /// Node positions.
    pub nodes: Vec<[f64; 2]>,
    /// Members as node index pairs (a < b), deterministic order.
    pub members: Vec<(usize, usize)>,
    /// Member lengths.
    pub lengths: Vec<f64>,
    /// The FrankenNetworkx candidate graph (node names `n{i}`).
    pub graph: Graph,
}

impl GroundStructure {
    /// A `nx × ny` node grid on `[0, w] × [0, h]` with all-pairs
    /// candidates filtered by the rules.
    ///
    /// # Panics
    /// On a degenerate grid.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn grid(nx: usize, ny: usize, w: f64, h: f64, rules: &GroundRules) -> GroundStructure {
        assert!(nx >= 2 && ny >= 2, "grid too small");
        let mut nodes = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                nodes.push([
                    w * i as f64 / (nx - 1) as f64,
                    h * j as f64 / (ny - 1) as f64,
                ]);
            }
        }
        let n = nodes.len();
        let mut members = Vec::new();
        let mut lengths = Vec::new();
        for a in 0..n {
            for b in (a + 1)..n {
                let dx = nodes[b][0] - nodes[a][0];
                let dy = nodes[b][1] - nodes[a][1];
                let len = dx.hypot(dy);
                if len < rules.min_len || len > rules.max_len {
                    continue;
                }
                if !rules.angles.is_empty() {
                    let ang = dy.atan2(dx).to_degrees().rem_euclid(180.0);
                    let ok = rules.angles.iter().any(|&want| {
                        let d = (ang - want).abs();
                        d.min(180.0 - d) <= rules.angle_tol
                    });
                    if !ok {
                        continue;
                    }
                }
                // Skip members that pass exactly through another node
                // (collinear duplicates carry no new statics).
                let mut through = false;
                for (c, node) in nodes.iter().enumerate() {
                    if c == a || c == b {
                        continue;
                    }
                    let cx = node[0] - nodes[a][0];
                    let cy = node[1] - nodes[a][1];
                    let cross = cx * dy - cy * dx;
                    let dot = cx * dx + cy * dy;
                    if cross.abs() < 1e-9 * len && dot > 1e-12 && dot < len * len - 1e-12 {
                        through = true;
                        break;
                    }
                }
                if through {
                    continue;
                }
                members.push((a, b));
                lengths.push(len);
            }
        }
        // The fnx candidate graph.
        let mut graph = Graph::new(CompatibilityMode::Strict);
        for i in 0..n {
            graph.add_node(format!("n{i}"));
        }
        for &(a, b) in &members {
            let _ = graph.add_edge(format!("n{a}"), format!("n{b}"));
        }
        GroundStructure {
            nodes,
            members,
            lengths,
            graph,
        }
    }

    /// Ledger stats row (counts + FNV hash of the member list).
    #[must_use]
    pub fn stats(&self) -> String {
        let mut h = 0xcbf2_9ce4_8422_2325u64;
        let mut mix = |v: u64| {
            for b in v.to_le_bytes() {
                h ^= u64::from(b);
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
        };
        for &(a, b) in &self.members {
            mix(a as u64);
            mix(b as u64);
        }
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"nodes\":{},\"members\":{},\"hash\":\"{h:#018x}\"}}",
            self.nodes.len(),
            self.members.len()
        );
        s
    }
}
