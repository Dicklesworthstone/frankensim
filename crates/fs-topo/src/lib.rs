//! fs-topo — validity and topology certificates (plan §7.8). Layer: L2.
//!
//! Three certificate families, none of them sampling heuristics:
//!
//! - [`manifold_certificate`] — combinatorial manifoldness (edge-use,
//!   half-edge round-trip, orientability, closedness) plus geometric
//!   red flags (degenerate faces, fold-overs) and an OPTIONAL outward-
//!   orientation probe. Combinatorial defects are LOCALIZED to their
//!   faces/edges; the global inward-orientation verdict is reported as
//!   the global [`ManifoldDefect::InwardOrientation`], never disguised
//!   as a local edge. Geometric checks run only on an ADMITTED soup
//!   (in-range indices, finite coordinates) and their refusals are
//!   recorded rather than absorbed;
//! - [`self_intersection_certificate`] — non-intersection as a PROOF
//!   over an ADMITTED soup: sweep-and-prune broad phase, then an EXACT
//!   narrow phase built on fs-ivl's exact `orient3d`/`orient2d`. Faces
//!   sharing a vertex or an edge are decided against their SHARED
//!   FEATURE instead of being skipped, so a face that pierces a
//!   neighbour it touches is caught. A PASS on an admitted soup cannot
//!   be falsely claimed — the arithmetic is exact; exact-contact
//!   configurations are reported CONSERVATIVELY as touching (bounded,
//!   listed false-FAILs, per the acceptance contract);
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
    IntersectKind, SelfIntersectRefusal, SelfIntersectReport, self_intersection_certificate,
    tri_tri_intersect,
};

use fs_geom::{Point3, Vec3};
use fs_rep_mesh::{HalfEdgeMesh, Soup, winding_exact};
use std::collections::{BTreeMap, BTreeSet};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// One manifoldness defect. The combinatorial and geometric variants are
/// LOCALIZED to the face or edge that witnesses them; the two global
/// variants ([`Self::InwardOrientation`], [`Self::IndeterminateWinding`])
/// and the admission refusals say so in their own shape rather than
/// borrowing a sentinel edge.
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
    /// GLOBAL, not localized: every shared edge is traversed both ways —
    /// the surface is consistently oriented — but it faces INWARD at the
    /// probe, where the exact winding number is `winding` instead of +1.
    /// No individual edge is at fault, so none is named: repairing this
    /// means reversing the whole surface, not editing one edge.
    InwardOrientation {
        /// The interior probe the caller supplied.
        probe: Point3,
        /// The exact winding number observed there.
        winding: f64,
    },
    /// The outwardness probe produced no verdict: the probe itself was
    /// non-finite, or the winding sum was. Nothing is claimed about the
    /// surface's facing.
    IndeterminateWinding {
        /// The interior probe the caller supplied.
        probe: Point3,
    },
    /// A referenced vertex position carries a non-finite coordinate.
    /// EVERY position-dependent check (degenerate faces, fold-overs,
    /// outwardness) was skipped; only the combinatorial verdicts stand.
    NonFiniteVertex {
        /// The offending vertex index.
        vertex: u32,
    },
    /// A triangle references a vertex index outside `positions`. Every
    /// position-dependent check was skipped.
    VertexIndexOutOfRange {
        /// Face index.
        face: usize,
        /// The out-of-range vertex index.
        vertex: u32,
    },
}

/// The manifoldness certificate.
///
/// The three `bool` flags are purely COMBINATORIAL: they read vertex
/// indices only and are therefore meaningful even when the geometry was
/// refused. Everything that depends on coordinates lives in `outward`
/// (an explicit three-state answer) or in `defects`.
#[derive(Debug, Clone)]
pub struct ManifoldReport {
    /// Combinatorially manifold (every edge used exactly twice, half-
    /// edge structure builds and round-trips).
    pub manifold: bool,
    /// Closed (no boundary edges).
    pub closed: bool,
    /// Consistently oriented: every shared edge is traversed in BOTH
    /// directions by its two incident faces.
    ///
    /// This is a COMBINATORIAL claim only. It says nothing about which
    /// side of the surface the face normals point at — a mesh with every
    /// triangle reversed is still consistently oriented. Outwardness is
    /// reported separately in [`Self::outward`].
    pub consistently_oriented: bool,
    /// Outwardness verdict from the interior probe:
    ///
    /// - `Some(true)` — the exact winding number at the probe is +1, so
    ///   the (closed, consistently oriented, manifold) surface faces
    ///   OUTWARD;
    /// - `Some(false)` — a probe ran and the surface faces inward
    ///   (see [`ManifoldDefect::InwardOrientation`]);
    /// - `None` — NO outwardness check ran. Either no probe was
    ///   supplied, or the combinatorial checks already failed, or the
    ///   geometry was refused, or the winding was indeterminate. A
    ///   `None` here is a no-claim, not a pass.
    pub outward: Option<bool>,
    /// Defects and admission refusals (empty ⟺ nothing was found and
    /// nothing was refused).
    pub defects: Vec<ManifoldDefect>,
}

impl ManifoldReport {
    /// True when the soup is a closed, consistently oriented, manifold
    /// surface with no defect AND an interior probe CONFIRMED that it
    /// faces outward.
    ///
    /// Outwardness is required, not assumed: without a probe the
    /// certificate has no evidence about facing, so `certified()` is
    /// false. Use [`Self::combinatorially_certified`] for the weaker,
    /// probe-free claim.
    #[must_use]
    pub fn certified(&self) -> bool {
        self.combinatorially_certified() && self.outward == Some(true)
    }

    /// The weaker COMBINATORIAL certificate: closed, manifold, and
    /// consistently oriented, with no defect and no refusal — and NO
    /// claim about which way the surface faces.
    #[must_use]
    pub fn combinatorially_certified(&self) -> bool {
        self.manifold && self.closed && self.consistently_oriented && self.defects.is_empty()
    }
}

fn finite_point(p: Point3) -> bool {
    p.x.is_finite() && p.y.is_finite() && p.z.is_finite()
}

/// Admit the soup's geometry before any position-dependent check reads
/// it. Returns `(indices_in_range, geometry_admitted, refusal defects)`;
/// `geometry_admitted` implies `indices_in_range`.
fn admit(soup: &Soup) -> (bool, bool, Vec<ManifoldDefect>) {
    let mut defects = Vec::new();
    let vertex_count = soup.positions.len();
    let mut indices_in_range = true;
    for (face, t) in soup.triangles.iter().enumerate() {
        for &vertex in t {
            if vertex as usize >= vertex_count {
                indices_in_range = false;
                defects.push(ManifoldDefect::VertexIndexOutOfRange { face, vertex });
            }
        }
    }
    if !indices_in_range {
        return (false, false, defects);
    }
    let mut non_finite: BTreeSet<u32> = BTreeSet::new();
    for t in &soup.triangles {
        for &vertex in t {
            if !finite_point(soup.positions[vertex as usize]) {
                non_finite.insert(vertex);
            }
        }
    }
    let geometry_admitted = non_finite.is_empty();
    defects.extend(
        non_finite
            .into_iter()
            .map(|vertex| ManifoldDefect::NonFiniteVertex { vertex }),
    );
    (true, geometry_admitted, defects)
}

/// Combinatorial + geometric manifoldness with defect localization.
/// `interior_probe` is a point expected inside (outwardness check);
/// pass `None` to skip it — but note that skipping it leaves
/// [`ManifoldReport::outward`] at `None` and makes
/// [`ManifoldReport::certified`] false, because no evidence about the
/// surface's facing was gathered.
///
/// The geometry is ADMITTED before it is read: out-of-range vertex
/// indices and non-finite coordinates are recorded as defects and every
/// position-dependent check (degenerate faces, fold-overs, outwardness)
/// is SKIPPED. Those checks are sign/threshold comparisons that a NaN
/// quietly falsifies, so a NaN soup would otherwise certify.
#[must_use]
pub fn manifold_certificate(soup: &Soup, interior_probe: Option<Point3>) -> ManifoldReport {
    // ---- Admission first: everything below indexes positions.
    let (indices_in_range, geometry_admitted, mut defects) = admit(soup);
    // Edge-use census with direction bookkeeping.
    let mut uses: BTreeMap<[u32; 2], (u32, i32)> = BTreeMap::new();
    for (fi, t) in soup.triangles.iter().enumerate() {
        // Degeneracy: repeated indices always; zero area only when the
        // coordinates were admitted (`n.norm() < 1e-30` is false for NaN).
        let [a, b, c] = *t;
        let repeated = a == b || b == c || a == c;
        let zero_area = geometry_admitted && {
            let pa = soup.positions[a as usize];
            let n = cross(
                soup.positions[b as usize].delta_from(pa),
                soup.positions[c as usize].delta_from(pa),
            );
            n.norm() < 1e-30
        };
        if repeated || zero_area {
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
    let mut consistently_oriented = true;
    for (&edge, &(count, dir_sum)) in &uses {
        match count {
            1 => {
                closed = false;
                defects.push(ManifoldDefect::BoundaryEdge { edge });
            }
            2 => {
                // Two uses must traverse in OPPOSITE directions.
                if dir_sum != 0 {
                    consistently_oriented = false;
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
    // Combinatorial, but the builder needs in-range indices.
    if manifold && closed && consistently_oriented && indices_in_range {
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
    // A NaN coordinate makes every comparison there false, so this runs
    // only on admitted geometry.
    if geometry_admitted {
        push_folded_edges(soup, &mut defects);
    }
    // Outwardness (only meaningful when closed + consistent + admitted).
    let mut outward = None;
    if geometry_admitted
        && closed
        && consistently_oriented
        && manifold
        && let Some(p) = interior_probe
    {
        outward = probe_outwardness(soup, p, &mut defects);
    }
    ManifoldReport {
        manifold,
        closed,
        consistently_oriented,
        outward,
        defects,
    }
}

/// Fold-over red flags: two faces sharing an edge with near-antiparallel
/// normals. Requires ADMITTED geometry — every comparison here is false
/// for a NaN coordinate.
fn push_folded_edges(soup: &Soup, defects: &mut Vec<ManifoldDefect>) {
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
}

/// The outwardness verdict from one interior probe. `None` is a NO-CLAIM:
/// the probe or the winding sum was non-finite, so nothing about the
/// surface's facing was established.
fn probe_outwardness(soup: &Soup, p: Point3, defects: &mut Vec<ManifoldDefect>) -> Option<bool> {
    let w = if finite_point(p) {
        winding_exact(soup, p)
    } else {
        f64::NAN
    };
    if !w.is_finite() {
        defects.push(ManifoldDefect::IndeterminateWinding { probe: p });
        return None;
    }
    let is_outward = (w - 1.0).abs() <= 0.5;
    if !is_outward {
        defects.push(ManifoldDefect::InwardOrientation {
            probe: p,
            winding: w,
        });
    }
    Some(is_outward)
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
