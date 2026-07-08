//! The RT0 mesh substrate: triangles plus a GLOBALLY ORIENTED edge
//! table. Every interior edge carries one fixed unit normal (owned by
//! the lower-index adjacent triangle); each triangle records a sign
//! per local edge telling whether its outward normal agrees with the
//! global one — the whole H(div) bookkeeping in two arrays.

use fs_solid::Mesh2;

/// One edge of the triangulation.
#[derive(Debug, Clone)]
pub struct Edge {
    /// Endpoint vertex indices (a < b).
    pub verts: (usize, usize),
    /// Adjacent triangles (t1 = usize::MAX on the boundary).
    pub tris: (usize, usize),
    /// Length.
    pub len: f64,
    /// Global unit normal (outward from `tris.0`).
    pub normal: [f64; 2],
    /// Midpoint.
    pub mid: [f64; 2],
}

/// Triangle mesh with oriented edges (the RT0 substrate).
pub struct TriMesh {
    /// Vertex positions.
    pub verts: Vec<[f64; 2]>,
    /// Triangles (vertex triples, counterclockwise).
    pub tris: Vec<[usize; 3]>,
    /// Edges.
    pub edges: Vec<Edge>,
    /// Per-triangle local-edge → (global edge, sign): +1 when the
    /// triangle's outward normal equals the global normal.
    pub tri_edges: Vec<[(usize, f64); 3]>,
    /// Triangle areas.
    pub areas: Vec<f64>,
    /// Triangle centroids.
    pub centroids: Vec<[f64; 2]>,
}

impl TriMesh {
    /// Build from an fs-solid triangle mesh (CCW triangles).
    ///
    /// # Panics
    /// On non-triangle elements or degenerate geometry.
    #[must_use]
    pub fn from_mesh2(m: &Mesh2) -> TriMesh {
        let verts: Vec<[f64; 2]> = m.nodes.clone();
        let tris: Vec<[usize; 3]> = m
            .elems
            .iter()
            .map(|e| {
                assert_eq!(e.len(), 3, "triangles only");
                [e[0], e[1], e[2]]
            })
            .collect();
        let mut edge_map: std::collections::BTreeMap<(usize, usize), usize> =
            std::collections::BTreeMap::new();
        let mut edges: Vec<Edge> = Vec::new();
        let mut tri_edges = vec![[(0usize, 0.0f64); 3]; tris.len()];
        let mut areas = Vec::with_capacity(tris.len());
        let mut centroids = Vec::with_capacity(tris.len());
        for (t, tri) in tris.iter().enumerate() {
            let p: [[f64; 2]; 3] = core::array::from_fn(|k| verts[tri[k]]);
            let area = 0.5
                * ((p[1][0] - p[0][0]) * (p[2][1] - p[0][1])
                    - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]));
            assert!(area > 1e-14, "CCW nondegenerate triangles required");
            areas.push(area);
            centroids.push([
                (p[0][0] + p[1][0] + p[2][0]) / 3.0,
                (p[0][1] + p[1][1] + p[2][1]) / 3.0,
            ]);
            // Local edge k is OPPOSITE vertex k: (k+1, k+2).
            for k in 0..3 {
                let (a, b) = (tri[(k + 1) % 3], tri[(k + 2) % 3]);
                let key = (a.min(b), a.max(b));
                let (pa, pb) = (verts[key.0], verts[key.1]);
                let dx = pb[0] - pa[0];
                let dy = pb[1] - pa[1];
                let len = dx.hypot(dy);
                // Outward normal of THIS triangle on this edge: rotate
                // the CCW edge direction (from tri[(k+1)] to tri[(k+2)])
                // by −90°.
                let (ex, ey) = (
                    verts[tri[(k + 2) % 3]][0] - verts[tri[(k + 1) % 3]][0],
                    verts[tri[(k + 2) % 3]][1] - verts[tri[(k + 1) % 3]][1],
                );
                let outward = [ey / len, -ex / len];
                let idx = *edge_map.entry(key).or_insert_with(|| {
                    edges.push(Edge {
                        verts: key,
                        tris: (t, usize::MAX),
                        len,
                        normal: outward,
                        mid: [
                            f64::midpoint(pa[0], pb[0]),
                            f64::midpoint(pa[1], pb[1]),
                        ],
                    });
                    edges.len() - 1
                });
                let e = &mut edges[idx];
                let sign = if e.tris.0 == t {
                    1.0
                } else {
                    e.tris.1 = t;
                    -1.0
                };
                tri_edges[t][k] = (idx, sign);
            }
        }
        TriMesh {
            verts,
            tris,
            edges,
            tri_edges,
            areas,
            centroids,
        }
    }

    /// Is this edge on the boundary?
    #[must_use]
    pub fn is_boundary(&self, e: usize) -> bool {
        self.edges[e].tris.1 == usize::MAX
    }

    /// RT0 basis value of local edge `k` of triangle `t` at point `x`
    /// (GLOBAL orientation folded in).
    #[must_use]
    pub fn rt0(&self, t: usize, k: usize, x: [f64; 2]) -> [f64; 2] {
        let (e, sign) = self.tri_edges[t][k];
        let opp = self.verts[self.tris[t][k]];
        let c = sign * self.edges[e].len / (2.0 * self.areas[t]);
        [c * (x[0] - opp[0]), c * (x[1] - opp[1])]
    }

    /// Divergence of the RT0 basis of local edge `k` on triangle `t`
    /// (constant).
    #[must_use]
    pub fn rt0_div(&self, t: usize, k: usize) -> f64 {
        let (e, sign) = self.tri_edges[t][k];
        sign * self.edges[e].len / self.areas[t]
    }

    /// Gradient (2×2, row = component) of the RT0 basis — constant per
    /// cell: ∇φ = c·I.
    #[must_use]
    pub fn rt0_grad(&self, t: usize, k: usize) -> f64 {
        let (e, sign) = self.tri_edges[t][k];
        sign * self.edges[e].len / (2.0 * self.areas[t])
    }
}
