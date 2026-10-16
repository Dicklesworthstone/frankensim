//! Half-edge surface meshes (plan §7.2): manifold connectivity with
//! boundary support and an Euler-operator editing core. The connectivity
//! INVARIANTS are the contract — `twin(twin(h)) == h`, `next` cycles
//! faces, boundary half-edges carry no face — and the property battery
//! maintains them under random edit sequences.

use fs_geom::Point3;

/// Sentinel for "no face" (boundary half-edges).
pub const NO_FACE: u32 = u32::MAX;

/// One half-edge record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HalfEdge {
    /// Origin vertex index.
    pub origin: u32,
    /// Opposite half-edge index.
    pub twin: u32,
    /// Next half-edge around the face (or around the boundary loop).
    pub next: u32,
    /// Incident face, or [`NO_FACE`] on the boundary.
    pub face: u32,
}

/// A manifold triangle mesh in half-edge form.
#[derive(Debug, Clone)]
pub struct HalfEdgeMesh {
    /// Vertex positions.
    pub positions: Vec<Point3>,
    /// Half-edge records.
    pub half_edges: Vec<HalfEdge>,
    /// One half-edge per face.
    pub face_edge: Vec<u32>,
}

/// Structured build failure (Decalogue P10).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MeshBuildError {
    /// An edge is shared by more than two faces (non-manifold).
    NonManifoldEdge {
        /// The offending vertex pair.
        edge: (u32, u32),
    },
    /// A face references a vertex that does not exist.
    VertexOutOfRange {
        /// The offending face index.
        face: u32,
    },
    /// A boundary vertex belongs to more than one boundary loop.
    NonManifoldVertex {
        /// The offending vertex.
        vertex: u32,
    },
}

impl core::fmt::Display for MeshBuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            MeshBuildError::NonManifoldEdge { edge } => write!(
                f,
                "half-edge build refused: edge ({}, {}) is shared by more than two faces — \
                 run the repair suite's soup path first (winding-number classification does \
                 not need manifoldness; half-edge connectivity does)",
                edge.0, edge.1
            ),
            MeshBuildError::VertexOutOfRange { face } => {
                write!(
                    f,
                    "half-edge build refused: face {face} references a missing vertex"
                )
            }
            MeshBuildError::NonManifoldVertex { vertex } => write!(
                f,
                "half-edge build refused: boundary vertex {vertex} sits on more than one \
                 boundary loop (non-manifold pinch) — repair the soup first"
            ),
        }
    }
}

impl core::error::Error for MeshBuildError {}

impl HalfEdgeMesh {
    /// Build from indexed triangles. Boundary edges get twin half-edges
    /// with [`NO_FACE`]; boundary `next` links walk each boundary loop.
    ///
    /// # Errors
    /// [`MeshBuildError`] on non-manifold or out-of-range input (the soup
    /// path is the repair suite, not this constructor).
    pub fn from_triangles(
        positions: Vec<Point3>,
        triangles: &[[u32; 3]],
    ) -> Result<HalfEdgeMesh, MeshBuildError> {
        use std::collections::BTreeMap;
        let nv = positions.len() as u32;
        let mut half_edges: Vec<HalfEdge> = Vec::with_capacity(triangles.len() * 3);
        let mut face_edge = Vec::with_capacity(triangles.len());
        let mut edge_map: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        for (fi, tri) in triangles.iter().enumerate() {
            if tri.iter().any(|&v| v >= nv) {
                return Err(MeshBuildError::VertexOutOfRange { face: fi as u32 });
            }
            let base = half_edges.len() as u32;
            face_edge.push(base);
            for c in 0..3 {
                let (a, b) = (tri[c], tri[(c + 1) % 3]);
                half_edges.push(HalfEdge {
                    origin: a,
                    twin: u32::MAX,
                    next: base + ((c as u32 + 1) % 3),
                    face: fi as u32,
                });
                let he = base + c as u32;
                if let Some(&other) = edge_map.get(&(b, a)) {
                    if half_edges[other as usize].twin != u32::MAX {
                        return Err(MeshBuildError::NonManifoldEdge { edge: (a, b) });
                    }
                    half_edges[he as usize].twin = other;
                    half_edges[other as usize].twin = he;
                } else if edge_map.insert((a, b), he).is_some() {
                    return Err(MeshBuildError::NonManifoldEdge { edge: (a, b) });
                }
            }
        }
        // Boundary loops: create NO_FACE twins for unmatched half-edges,
        // then link their `next` pointers by walking around each loop.
        let mut boundary_of: BTreeMap<u32, u32> = BTreeMap::new(); // origin -> boundary he
        let unmatched: Vec<u32> = (0..half_edges.len() as u32)
            .filter(|&h| half_edges[h as usize].twin == u32::MAX)
            .collect();
        for &h in &unmatched {
            let inner = half_edges[h as usize];
            let dest = half_edges[inner.next as usize].origin;
            let b = half_edges.len() as u32;
            half_edges.push(HalfEdge {
                origin: dest,
                twin: h,
                next: u32::MAX,
                face: NO_FACE,
            });
            half_edges[h as usize].twin = b;
            if boundary_of.insert(dest, b).is_some() {
                return Err(MeshBuildError::NonManifoldVertex { vertex: dest });
            }
        }
        let boundary: Vec<u32> = boundary_of.values().copied().collect();
        for &b in &boundary {
            // The next boundary half-edge starts where this one ends: the
            // origin of this one's twin.
            let dest = half_edges[half_edges[b as usize].twin as usize].origin;
            let next = boundary_of
                .get(&dest)
                .copied()
                .expect("manifold boundary loops close");
            half_edges[b as usize].next = next;
        }
        Ok(HalfEdgeMesh {
            positions,
            half_edges,
            face_edge,
        })
    }

    /// Face count.
    #[must_use]
    pub fn face_count(&self) -> usize {
        self.face_edge.len()
    }

    /// Undirected edge count.
    #[must_use]
    pub fn edge_count(&self) -> usize {
        self.half_edges.len() / 2
    }

    /// The three vertex indices of a face.
    #[must_use]
    pub fn face_vertices(&self, face: u32) -> [u32; 3] {
        let h0 = self.face_edge[face as usize];
        let h1 = self.half_edges[h0 as usize].next;
        let h2 = self.half_edges[h1 as usize].next;
        [
            self.half_edges[h0 as usize].origin,
            self.half_edges[h1 as usize].origin,
            self.half_edges[h2 as usize].origin,
        ]
    }

    /// Check every structural invariant; returns the first violation as a
    /// teaching string (`None` = healthy). The property battery calls this
    /// after every random edit.
    #[must_use]
    pub fn check_invariants(&self) -> Option<String> {
        for (i, he) in self.half_edges.iter().enumerate() {
            let i = i as u32;
            let Some(twin) = self.half_edges.get(he.twin as usize) else {
                return Some(format!("half-edge {i} has an out-of-range twin"));
            };
            if twin.twin != i {
                return Some(format!("twin(twin({i})) = {} != {i}", twin.twin));
            }
            if he.face != NO_FACE {
                // Interior: next³ must return home (triangles).
                let n1 = self.half_edges[he.next as usize];
                let n2 = self.half_edges[n1.next as usize];
                if n2.next != i {
                    return Some(format!("face cycle at {i} does not close in 3 steps"));
                }
                // twin's origin is this half-edge's destination.
                if twin.origin != n1.origin {
                    return Some(format!("twin origin mismatch at {i}"));
                }
            }
        }
        for (fi, &h) in self.face_edge.iter().enumerate() {
            if self.half_edges[h as usize].face != fi as u32 {
                return Some(format!("face_edge[{fi}] points at a foreign half-edge"));
            }
        }
        None
    }

    /// Euler characteristic V − E + F (boundary loops not added back —
    /// callers with genus/boundary knowledge interpret it).
    #[must_use]
    pub fn euler_characteristic(&self) -> i64 {
        i64::try_from(self.positions.len()).expect("mesh sizes fit i64")
            - i64::try_from(self.edge_count()).expect("mesh sizes fit i64")
            + i64::try_from(self.face_count()).expect("mesh sizes fit i64")
    }

    /// Flip an interior edge (the Euler-operator editing core's first
    /// citizen). Returns false (a non-event) when the half-edge is on the
    /// boundary or the flip would create a duplicate edge.
    pub fn flip_edge(&mut self, h: u32) -> bool {
        let t = self.half_edges[h as usize].twin;
        let (hf, tf) = (
            self.half_edges[h as usize].face,
            self.half_edges[t as usize].face,
        );
        if hf == NO_FACE || tf == NO_FACE {
            return false;
        }
        // Gather the quad: h runs a->b with face (a,b,c); t runs b->a
        // with face (b,a,d).
        let hn = self.half_edges[h as usize].next;
        let hp = self.half_edges[hn as usize].next;
        let tn = self.half_edges[t as usize].next;
        let tp = self.half_edges[tn as usize].next;
        let c = self.half_edges[hp as usize].origin;
        let d = self.half_edges[tp as usize].origin;
        // Refuse flips that would duplicate an existing edge (c, d).
        // Linear scan: correctness over speed (remeshing-scale edit
        // throughput is the anisotropic-remesh bead's problem).
        for he in &self.half_edges {
            if he.origin == c && self.half_edges[he.twin as usize].origin == d {
                return false;
            }
        }
        // Relink (Botsch et al.): h/t become the diagonal c->d / d->c;
        // face hf becomes (c, d, b) with cycle h -> tp -> hn, face tf
        // becomes (d, c, a) with cycle t -> hp -> tn. The invariant checker
        // in the property battery audits every flip.
        self.half_edges[h as usize].origin = c;
        self.half_edges[t as usize].origin = d;
        self.half_edges[h as usize].next = tp;
        self.half_edges[tp as usize].next = hn;
        self.half_edges[hn as usize].next = h;
        self.half_edges[t as usize].next = hp;
        self.half_edges[hp as usize].next = tn;
        self.half_edges[tn as usize].next = t;
        self.half_edges[tp as usize].face = hf;
        self.half_edges[hn as usize].face = hf;
        self.half_edges[hp as usize].face = tf;
        self.half_edges[tn as usize].face = tf;
        self.face_edge[hf as usize] = h;
        self.face_edge[tf as usize] = t;
        true
    }
}
