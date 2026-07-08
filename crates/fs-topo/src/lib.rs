//! fs-topo — validity and topology certificates (plan §7.8). Layer: L2.
//!
//! Three certificate families, none of them sampling heuristics:
//!
//! - [`manifold_certificate`] — combinatorial manifoldness (edge-use,
//!   half-edge round-trip, orientability, closedness) plus geometric
//!   red flags (degenerate faces, fold-overs), with every defect
//!   LOCALIZED to its faces/edges;
//! - [`self_intersection_certificate`] — non-intersection as a PROOF:
//!   sweep-and-prune broad phase, then an EXACT narrow phase built on
//!   fs-ivl's exact `orient3d`/`orient2d` (Guigue–Devillers with exact
//!   signs). A PASS cannot be falsely claimed — the arithmetic is
//!   exact; exact-contact configurations are reported CONSERVATIVELY
//!   as touching (bounded, listed false-FAILs, per the acceptance
//!   contract);
//! - [`crate::cubical`] — Betti numbers of voxel solids by union-find
//!   plus exact Euler characteristic duality, true 0-dimensional
//!   persistence (elder rule over the filtration), persistence-aware
//!   feature counting, and chart-level topology verification with
//!   HONEST resolution caveats.

pub mod cubical;
mod intersect;

#[cfg(feature = "moonshot-topo-persistence")]
pub mod penalty;

pub use intersect::{
    IntersectKind, SelfIntersectReport, self_intersection_certificate, tri_tri_intersect,
};

use fs_geom::{Point3, Vec3};
use fs_rep_mesh::{HalfEdgeMesh, Soup, winding_exact};
use std::collections::BTreeMap;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// One localized manifoldness defect.
#[derive(Debug, Clone, PartialEq)]
pub enum ManifoldDefect {
    /// An edge used by only one face (open boundary).
    BoundaryEdge {
        /// Vertex pair.
        edge: [u32; 2],
    },
    /// An edge used by more than two faces (fin).
    NonManifoldEdge {
        /// Vertex pair.
        edge: [u32; 2],
        /// How many faces use it.
        uses: u32,
    },
    /// Two faces traverse a shared edge in the SAME direction
    /// (inconsistent orientation).
    MisorientedEdge {
        /// Vertex pair.
        edge: [u32; 2],
    },
    /// A zero-area or repeated-vertex face.
    DegenerateFace {
        /// Face index.
        face: usize,
    },
    /// Adjacent faces folded back onto each other (dihedral ≈ π).
    FoldedEdge {
        /// Vertex pair.
        edge: [u32; 2],
    },
    /// The half-edge builder refused outright (its teaching text).
    BuildRefusal {
        /// The builder's message.
        message: String,
    },
}

/// The manifoldness certificate.
#[derive(Debug, Clone)]
pub struct ManifoldReport {
    /// Combinatorially manifold (every edge used exactly twice, half-
    /// edge structure builds and round-trips).
    pub manifold: bool,
    /// Closed (no boundary edges).
    pub closed: bool,
    /// Consistently oriented (every shared edge traversed both ways),
    /// and outward (winding +1 at the interior probe) when closed.
    pub oriented: bool,
    /// Localized defects (empty ⟺ all three flags hold).
    pub defects: Vec<ManifoldDefect>,
}

impl ManifoldReport {
    /// True when the soup is a closed, oriented, manifold surface.
    #[must_use]
    pub fn certified(&self) -> bool {
        self.manifold && self.closed && self.oriented && self.defects.is_empty()
    }
}

/// Combinatorial + geometric manifoldness with defect localization.
/// `interior_probe` is a point expected inside (outwardness check);
/// pass `None` to skip the orientation-sign check.
#[must_use]
pub fn manifold_certificate(soup: &Soup, interior_probe: Option<Point3>) -> ManifoldReport {
    let mut defects = Vec::new();
    // Edge-use census with direction bookkeeping.
    let mut uses: BTreeMap<[u32; 2], (u32, i32)> = BTreeMap::new();
    for (fi, t) in soup.triangles.iter().enumerate() {
        // Degeneracy: repeated indices or zero area.
        let [a, b, c] = *t;
        let repeated = a == b || b == c || a == c;
        let pa = soup.positions[a as usize];
        let n = cross(
            soup.positions[b as usize].delta_from(pa),
            soup.positions[c as usize].delta_from(pa),
        );
        if repeated || n.norm() < 1e-30 {
            defects.push(ManifoldDefect::DegenerateFace { face: fi });
        }
        for k in 0..3 {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            let key = if u < v { [u, v] } else { [v, u] };
            let dir = if u < v { 1 } else { -1 };
            let e = uses.entry(key).or_insert((0, 0));
            e.0 += 1;
            e.1 += dir;
        }
    }
    let mut closed = true;
    let mut manifold = true;
    let mut oriented = true;
    for (&edge, &(count, dir_sum)) in &uses {
        match count {
            1 => {
                closed = false;
                defects.push(ManifoldDefect::BoundaryEdge { edge });
            }
            2 => {
                // Two uses must traverse in OPPOSITE directions.
                if dir_sum != 0 {
                    oriented = false;
                    defects.push(ManifoldDefect::MisorientedEdge { edge });
                }
            }
            n => {
                manifold = false;
                defects.push(ManifoldDefect::NonManifoldEdge { edge, uses: n });
            }
        }
    }
    // Half-edge round-trip (vertex-link conditions live in the builder).
    if manifold && closed && oriented {
        match HalfEdgeMesh::from_triangles(soup.positions.clone(), &soup.triangles) {
            Ok(he) => {
                if let Some(v) = he.check_invariants() {
                    manifold = false;
                    defects.push(ManifoldDefect::BuildRefusal { message: v });
                }
            }
            Err(e) => {
                manifold = false;
                defects.push(ManifoldDefect::BuildRefusal {
                    message: e.to_string(),
                });
            }
        }
    }
    // Fold-over red flags: adjacent faces with near-antiparallel normals.
    let mut face_of: BTreeMap<[u32; 2], Vec<usize>> = BTreeMap::new();
    for (fi, t) in soup.triangles.iter().enumerate() {
        for k in 0..3 {
            let (u, v) = (t[k], t[(k + 1) % 3]);
            let key = if u < v { [u, v] } else { [v, u] };
            face_of.entry(key).or_default().push(fi);
        }
    }
    for (&edge, fs) in &face_of {
        if let [f1, f2] = fs.as_slice() {
            let n1 = face_normal(soup, *f1);
            let n2 = face_normal(soup, *f2);
            let den = n1.norm() * n2.norm();
            if den > 1e-30 && n1.dot(n2) / den < -0.999 {
                defects.push(ManifoldDefect::FoldedEdge { edge });
            }
        }
    }
    // Outwardness (only meaningful when closed + consistent).
    if closed
        && oriented
        && manifold
        && let Some(p) = interior_probe
    {
        let w = winding_exact(soup, p);
        if (w - 1.0).abs() > 0.5 {
            oriented = false;
            defects.push(ManifoldDefect::MisorientedEdge { edge: [0, 0] });
        }
    }
    ManifoldReport {
        manifold,
        closed,
        oriented,
        defects,
    }
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn face_normal(soup: &Soup, f: usize) -> Vec3 {
    let [a, b, c] = soup.triangles[f].map(|v| soup.positions[v as usize]);
    cross(b.delta_from(a), c.delta_from(a))
}

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
