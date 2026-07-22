//! The conduction mesh: an `fs-feec` tet complex with its element
//! geometry and its EXTRACTED boundary — the set of faces incident to
//! exactly one tet, each with an area, an outward unit normal, and a
//! centroid.
//!
//! Nothing here re-derives element geometry: signed volumes and the
//! constant barycentric gradients ∇λ_a come from
//! [`fs_feec::element_geometry`], which routes them through fs-la's
//! batched small-dense kernels. This module adds only the boundary
//! extraction and the degeneracy pre-check that turns fs-feec's
//! (correct, but panicking) degeneracy assertion into a typed refusal.

use fs_feec::{ElementGeometry, element_geometry};
use fs_rep_mesh::TetComplex;

use crate::ConductionError;

/// One face of the mesh boundary: incident to exactly one tetrahedron.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryFace {
    /// Index into the complex's canonical face table.
    pub face: usize,
    /// The face's SORTED vertex triple (the complex's canonical order).
    pub vertices: [u32; 3],
    /// Index of the owning tetrahedron.
    pub element: usize,
    /// Area in m².
    pub area: f64,
    /// Unit normal pointing OUT of the domain.
    pub outward_normal: [f64; 3],
    /// Face centroid (the arithmetic mean of the three vertices).
    pub centroid: [f64; 3],
}

/// A tet complex prepared for conduction assembly.
///
/// Not `Clone`/`Debug`: `fs_feec::ElementGeometry` is neither, and
/// re-deriving those here would either copy fs-feec's data model or
/// print a megabyte of gradients. Share a mesh by reference.
pub struct ConductionMesh {
    complex: TetComplex,
    positions: Vec<[f64; 3]>,
    geometry: ElementGeometry,
    boundary: Vec<BoundaryFace>,
    total_volume: f64,
}

impl core::fmt::Debug for ConductionMesh {
    /// Shape only. A conduction mesh holds one gradient triple per
    /// element; printing them would bury the diagnosis it was meant to
    /// support.
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ConductionMesh")
            .field("vertices", &self.complex.vertex_count)
            .field("elements", &self.complex.tets.len())
            .field("boundary_faces", &self.boundary.len())
            .field("total_volume", &self.total_volume)
            .finish_non_exhaustive()
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -(a[2] * b[1])),
        a[2].mul_add(b[0], -(a[0] * b[2])),
        a[0].mul_add(b[1], -(a[1] * b[0])),
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Relative degeneracy floor: a tet whose |signed volume| is below this
/// fraction of the mean element volume cannot carry a usable gradient
/// operator and is refused rather than inverted.
const DEGENERACY_FLOOR: f64 = 1e-12;

impl ConductionMesh {
    /// Prepare a complex for assembly.
    ///
    /// # Errors
    /// [`ConductionError::Mesh`] for a positions/complex size mismatch,
    /// an empty complex, or a non-finite coordinate;
    /// [`ConductionError::DegenerateElement`] when a tet's signed volume
    /// is below the relative degeneracy floor.
    pub fn new(
        complex: TetComplex,
        positions: Vec<[f64; 3]>,
    ) -> Result<ConductionMesh, ConductionError> {
        if positions.len() != complex.vertex_count {
            return Err(ConductionError::Mesh {
                what: format!(
                    "{} positions for {} vertices",
                    positions.len(),
                    complex.vertex_count
                ),
                fix: "supply exactly one position per complex vertex".to_string(),
            });
        }
        if complex.tets.is_empty() {
            return Err(ConductionError::Mesh {
                what: "the complex has no tetrahedra".to_string(),
                fix: "supply a complex with at least one element".to_string(),
            });
        }
        for (v, p) in positions.iter().enumerate() {
            for (k, &c) in p.iter().enumerate() {
                if !c.is_finite() {
                    return Err(ConductionError::Mesh {
                        what: format!("vertex {v} coordinate {k} is not finite"),
                        fix: "supply finite coordinates".to_string(),
                    });
                }
            }
        }

        // Degeneracy pre-check: fs-feec's element_geometry ASSERTS on a
        // singular Jacobian. Refuse first so a bad mesh becomes a typed
        // value instead of a panic crossing the crate boundary.
        let mut raw = Vec::with_capacity(complex.tets.len());
        let mut total = 0.0f64;
        for tet in &complex.tets {
            let p0 = positions[tet[0] as usize];
            let e1 = sub(positions[tet[1] as usize], p0);
            let e2 = sub(positions[tet[2] as usize], p0);
            let e3 = sub(positions[tet[3] as usize], p0);
            let vol = dot(cross(e1, e2), e3) / 6.0;
            total += vol.abs();
            raw.push(vol);
        }
        let mean = total / raw.len() as f64;
        for (e, &vol) in raw.iter().enumerate() {
            if !vol.is_finite() || vol.abs() <= DEGENERACY_FLOOR * mean {
                return Err(ConductionError::DegenerateElement {
                    element: e,
                    signed_volume: vol,
                });
            }
        }

        let geometry = element_geometry(&complex, &positions);
        let boundary = extract_boundary(&complex, &positions)?;
        Ok(ConductionMesh {
            complex,
            positions,
            geometry,
            boundary,
            total_volume: total,
        })
    }

    /// The underlying complex.
    #[must_use]
    pub const fn complex(&self) -> &TetComplex {
        &self.complex
    }

    /// Vertex positions in canonical vertex order.
    #[must_use]
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// The fs-feec element geometry (signed volumes, ∇λ, gradient Gram).
    #[must_use]
    pub const fn geometry(&self) -> &ElementGeometry {
        &self.geometry
    }

    /// The extracted boundary faces, in ascending canonical face order.
    #[must_use]
    pub fn boundary(&self) -> &[BoundaryFace] {
        &self.boundary
    }

    /// Vertex count (the P₁ degree-of-freedom count before Dirichlet
    /// elimination).
    #[must_use]
    pub const fn vertex_count(&self) -> usize {
        self.complex.vertex_count
    }

    /// Element (tetrahedron) count.
    #[must_use]
    pub fn element_count(&self) -> usize {
        self.complex.tets.len()
    }

    /// Unsigned volume of element `e` in m³.
    #[must_use]
    pub fn element_volume(&self, e: usize) -> f64 {
        self.geometry.vol_signed[e].abs()
    }

    /// Total unsigned mesh volume in m³.
    #[must_use]
    pub const fn total_volume(&self) -> f64 {
        self.total_volume
    }

    /// Total boundary area in m².
    #[must_use]
    pub fn boundary_area(&self) -> f64 {
        self.boundary.iter().map(|f| f.area).sum()
    }

    /// The mesh size parameter used by the G1 ladders: the longest edge
    /// of any element (a deterministic max over a fixed traversal).
    #[must_use]
    pub fn max_edge_length(&self) -> f64 {
        let mut worst = 0.0f64;
        for &[a, b] in &self.complex.edges {
            let d = sub(self.positions[b as usize], self.positions[a as usize]);
            worst = worst.max(fs_math::det::sqrt(dot(d, d)));
        }
        worst
    }
}

fn extract_boundary(
    complex: &TetComplex,
    positions: &[[f64; 3]],
) -> Result<Vec<BoundaryFace>, ConductionError> {
    let nf = complex.faces.len();
    let mut uses = vec![0u32; nf];
    // (owning tet, the tet vertex NOT on the face) — only meaningful for
    // faces used once, which is exactly the boundary set.
    let mut owner = vec![(usize::MAX, 0u32); nf];
    for (t, tet) in complex.tets.iter().enumerate() {
        for skip in 0..4 {
            let mut tri = [0u32; 3];
            let mut k = 0;
            for (i, &v) in tet.iter().enumerate() {
                if i != skip {
                    tri[k] = v;
                    k += 1;
                }
            }
            tri.sort_unstable();
            let fid = complex
                .faces
                .binary_search(&tri)
                .map_err(|_| ConductionError::Mesh {
                    what: format!("tet {t} face {tri:?} is missing from the complex face table"),
                    fix: "rebuild the complex with TetComplex::from_tets".to_string(),
                })?;
            uses[fid] += 1;
            owner[fid] = (t, tet[skip]);
        }
    }

    let mut out = Vec::new();
    for (fid, &count) in uses.iter().enumerate() {
        if count != 1 {
            continue;
        }
        let tri = complex.faces[fid];
        let p0 = positions[tri[0] as usize];
        let p1 = positions[tri[1] as usize];
        let p2 = positions[tri[2] as usize];
        let c = cross(sub(p1, p0), sub(p2, p0));
        let twice_area = fs_math::det::sqrt(dot(c, c));
        if !(twice_area.is_finite() && twice_area > 0.0) {
            return Err(ConductionError::Mesh {
                what: format!("boundary face {fid} has zero area"),
                fix: "remove the sliver element that produced it".to_string(),
            });
        }
        let mut n = [c[0] / twice_area, c[1] / twice_area, c[2] / twice_area];
        let (element, opposite) = owner[fid];
        // Outward = away from the tet's fourth vertex.
        if dot(n, sub(positions[opposite as usize], p0)) > 0.0 {
            n = [-n[0], -n[1], -n[2]];
        }
        let centroid = [
            (p0[0] + p1[0] + p2[0]) / 3.0,
            (p0[1] + p1[1] + p2[1]) / 3.0,
            (p0[2] + p1[2] + p2[2]) / 3.0,
        ];
        out.push(BoundaryFace {
            face: fid,
            vertices: tri,
            element,
            area: 0.5 * twice_area,
            outward_normal: n,
            centroid,
        });
    }
    Ok(out)
}
