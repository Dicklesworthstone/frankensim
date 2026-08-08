//! fs-plate — orthotropic thin-plate bending for instrument bodies (bead
//! frankensim-fsim-plates-shells-kj3s0, musical-acoustics program): DKT
//! (discrete Kirchhoff triangle) bending elements with an orthotropic
//! D-matrix from `fs_material::OrthotropicElastic`, membrane-prestress
//! geometric stiffness (drumheads, soundboard downbearing), offset
//! stiffener beams (bracing eccentricity), and modal analysis through
//! `fs_modal::slice_window` so every frequency arrives with a certified
//! interval and an inertia-certified count.
//!
//! DOF CONVENTION (pinned by the constant-curvature patch test, which fails
//! under any other mapping): each node carries `(w, wx, wy)` where `wx =
//! ∂w/∂x` and `wy = ∂w/∂y` are SLOPES, not signed rotations — no
//! rotation-sign zoo. Curvatures are `κ = [∂βx/∂x, ∂βy/∂y, ∂βx/∂y +
//! ∂βy/∂x]` with `β = ∇w` on the discrete Kirchhoff constraint lines, and
//! the bending energy density is `½ κᵀ D κ`.
//!
//! ELEMENT — Batoz's DKT in the compact area-coordinate form: quadratic
//! rotation fields `βx = Hxᵀu`, `βy = Hyᵀu` built from the six quadratic
//! shape functions and the per-edge geometry coefficients (a..e), with
//! Kirchhoff constraints enforced discretely at corners and mid-sides.
//! B-matrices are evaluated at the three mid-side points (degree-2 exact),
//! so `Kₑ = (A/3)·Σ_g BᵀDB` is the exact integral of the quadratic
//! integrand. Verification: rigid-body motions produce zero strain energy
//! to roundoff and a constant-curvature field reproduces `A·κᵀDκ` exactly
//! (both tested per element), the standard FEM patch certificates.
//!
//! Mass is LUMPED (translational `ρhA/3` per node, rotary `ρh³/12·A/3` on
//! slope DOFs) — SPD by construction, which the fs-modal pencil contract
//! requires; a consistent mass matrix is a recorded follow-up. Membrane
//! prestress enters as the standard P1 geometric stiffness `K_G = T·∫∇w·∇w`
//! (drumheads are the K_G-dominated, D→0 limit — tested against continuum
//! membrane frequencies). Stiffeners are 2-node Hermite beams on plate node
//! lines with parallel-axis eccentricity (`EI_eff = EI + EAe²`) plus GJ
//! torsion on the cross-slope.
//!
//! Determinism: assembly through `fs_sparse::Coo` canonical accumulation;
//! everything sequential; repeat runs bitwise (inherited and tested
//! downstream through fs-modal).

use fs_material::elastic::OrthotropicElastic;
pub use fs_modal::{ModalError, ModePair, SliceOptions, SliceReport};
use fs_sparse::{Coo, Csr};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Typed refusals. Display strings carry stable `FS-PLATE-*` codes.
#[derive(Debug, Clone, PartialEq)]
pub enum PlateError {
    /// A material/section parameter is non-positive or non-finite.
    BadSection {
        /// Which parameter refused.
        what: &'static str,
    },
    /// Mesh geometry refused: a triangle is degenerate (near-zero area) or
    /// inverted beyond the quality floor.
    DegenerateElement {
        /// Element index.
        element: usize,
        /// Its doubled signed area.
        twice_area: f64,
    },
    /// A stiffener path references nodes outside the mesh or is too short.
    BadStiffener {
        /// Diagnostic.
        what: &'static str,
    },
    /// Modal stage refusal, forwarded from fs-modal.
    Modal(ModalError),
}

impl core::fmt::Display for PlateError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PlateError::BadSection { what } => write!(f, "FS-PLATE-BAD-SECTION: {what}"),
            PlateError::DegenerateElement {
                element,
                twice_area,
            } => write!(
                f,
                "FS-PLATE-DEGENERATE-ELEMENT: element {element} has 2A = {twice_area:.3e}"
            ),
            PlateError::BadStiffener { what } => write!(f, "FS-PLATE-BAD-STIFFENER: {what}"),
            PlateError::Modal(e) => write!(f, "FS-PLATE-MODAL: {e}"),
        }
    }
}

impl std::error::Error for PlateError {}

impl From<ModalError> for PlateError {
    fn from(e: ModalError) -> PlateError {
        PlateError::Modal(e)
    }
}

// ---------------------------------------------------------------------------
// Section: orthotropic bending rigidity
// ---------------------------------------------------------------------------

/// Plate section: bending rigidity matrix D (3×3, Voigt `[κxx, κyy, 2κ...]`
/// row order `[Mxx, Myy, Mxy]`), thickness, and areal density. SI units
/// throughout (Pa, m, kg/m³); a dimensioned fs-qty front door is a recorded
/// follow-up on the bead.
#[derive(Debug, Clone, Copy)]
pub struct PlateSection {
    /// Bending rigidity `[[D11, D12, 0], [D12, D22, 0], [0, 0, D33]]`
    /// stored dense row-major 3×3 (a full anisotropic D is accepted).
    pub d: [f64; 9],
    /// Thickness h [m].
    pub thickness: f64,
    /// Volumetric density ρ [kg/m³].
    pub density: f64,
}

impl PlateSection {
    /// Orthotropic section with material axes 1,2 aligned to the plate x,y
    /// axes: `D11 = E1h³/12(1−ν12ν21)`, `D22 = E2h³/…`, `D12 = ν12E2h³/…`
    /// (symmetric since `ν21E1 = ν12E2`), `D33 = G12h³/12`. The grain axis
    /// (L) maps to material axis 1 — the matdb axis-convention contract.
    ///
    /// # Errors
    /// [`PlateError::BadSection`] on non-positive thickness/density or a
    /// non-positive-definite reduced stiffness.
    pub fn orthotropic(
        law: &OrthotropicElastic,
        thickness: f64,
        density: f64,
    ) -> Result<PlateSection, PlateError> {
        if !(thickness > 0.0 && thickness.is_finite()) {
            return Err(PlateError::BadSection {
                what: "thickness must be positive and finite",
            });
        }
        if !(density > 0.0 && density.is_finite()) {
            return Err(PlateError::BadSection {
                what: "density must be positive and finite",
            });
        }
        let e1 = law.e[0];
        let e2 = law.e[1];
        let nu12 = law.nu[0];
        let g12 = law.g[0];
        let nu21 = nu12 * e2 / e1;
        let denom = 1.0 - nu12 * nu21;
        if !denom.is_finite() || denom <= 0.0 {
            return Err(PlateError::BadSection {
                what: "1 - nu12*nu21 must be positive (reduced stiffness definiteness)",
            });
        }
        let h3 = thickness * thickness * thickness / 12.0;
        let d11 = e1 * h3 / denom;
        let d22 = e2 * h3 / denom;
        let d12 = nu12 * e2 * h3 / denom;
        let d33 = g12 * h3;
        Ok(PlateSection {
            d: [d11, d12, 0.0, d12, d22, 0.0, 0.0, 0.0, d33],
            thickness,
            density,
        })
    }

    /// Isotropic convenience: `E`, `ν` with `G = E/2(1+ν)`.
    ///
    /// # Errors
    /// [`PlateError::BadSection`] as for [`PlateSection::orthotropic`].
    pub fn isotropic(
        e: f64,
        nu: f64,
        thickness: f64,
        density: f64,
    ) -> Result<PlateSection, PlateError> {
        if !(e > 0.0 && e.is_finite() && nu > -1.0 && nu < 0.5) {
            return Err(PlateError::BadSection {
                what: "isotropic constants out of range",
            });
        }
        let g = e / (2.0 * (1.0 + nu));
        let law =
            OrthotropicElastic::new([e, e, e], [nu, nu, nu], [g, g, g], 1.0).map_err(|_| {
                PlateError::BadSection {
                    what: "isotropic constants rejected by compliance check",
                }
            })?;
        PlateSection::orthotropic(&law, thickness, density)
    }
}

// ---------------------------------------------------------------------------
// Mesh
// ---------------------------------------------------------------------------

/// Triangle mesh with 3 DOF per node `(w, wx, wy)`.
#[derive(Debug, Clone)]
pub struct PlateMesh {
    /// Node coordinates (x, y).
    pub nodes: Vec<(f64, f64)>,
    /// Triangles as CCW node triples.
    pub tris: Vec<[usize; 3]>,
}

impl PlateMesh {
    /// Structured right-triangle mesh of an a×b rectangle with
    /// `nx`×`ny` CELLS (so `(nx+1)(ny+1)` nodes, `2·nx·ny` triangles).
    /// Node (i, j) sits at `(i·a/nx, j·b/ny)`, index `j·(nx+1)+i`.
    #[must_use]
    pub fn rectangle(a: f64, b: f64, nx: usize, ny: usize) -> PlateMesh {
        let mut nodes = Vec::with_capacity((nx + 1) * (ny + 1));
        for j in 0..=ny {
            for i in 0..=nx {
                nodes.push((a * i as f64 / nx as f64, b * j as f64 / ny as f64));
            }
        }
        let idx = |i: usize, j: usize| j * (nx + 1) + i;
        let mut tris = Vec::with_capacity(2 * nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                tris.push([idx(i, j), idx(i + 1, j), idx(i + 1, j + 1)]);
                tris.push([idx(i, j), idx(i + 1, j + 1), idx(i, j + 1)]);
            }
        }
        PlateMesh { nodes, tris }
    }

    /// Node count.
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Boundary node indices of a [`PlateMesh::rectangle`] mesh.
    #[must_use]
    pub fn rectangle_boundary(nx: usize, ny: usize) -> Vec<usize> {
        let idx = |i: usize, j: usize| j * (nx + 1) + i;
        let mut out = Vec::new();
        for j in 0..=ny {
            for i in 0..=nx {
                if i == 0 || j == 0 || i == nx || j == ny {
                    out.push(idx(i, j));
                }
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// DKT element
// ---------------------------------------------------------------------------

/// Per-edge DKT geometry coefficients for edges 4=(2,3), 5=(3,1), 6=(1,2).
struct EdgeCoeff {
    a: [f64; 3],
    b: [f64; 3],
    c: [f64; 3],
    d: [f64; 3],
    e: [f64; 3],
}

fn edge_coeffs(x: &[f64; 3], y: &[f64; 3]) -> EdgeCoeff {
    let pairs = [(1usize, 2usize), (2, 0), (0, 1)]; // edges 4,5,6 as (i,j)
    let mut c = EdgeCoeff {
        a: [0.0; 3],
        b: [0.0; 3],
        c: [0.0; 3],
        d: [0.0; 3],
        e: [0.0; 3],
    };
    for (k, &(i, j)) in pairs.iter().enumerate() {
        let xij = x[i] - x[j];
        let yij = y[i] - y[j];
        let l2 = xij * xij + yij * yij;
        c.a[k] = -xij / l2;
        c.b[k] = 0.75 * xij * yij / l2;
        c.c[k] = (0.25 * xij * xij - 0.5 * yij * yij) / l2;
        c.d[k] = -yij / l2;
        c.e[k] = (0.25 * yij * yij - 0.5 * xij * xij) / l2;
    }
    c
}

/// The six quadratic shape functions' (ξ, η) derivatives at a point.
fn quad_shape_derivs(xi: f64, eta: f64) -> ([f64; 6], [f64; 6]) {
    let dxi = [
        4.0 * (xi + eta) - 3.0, // N1,ξ
        4.0 * xi - 1.0,         // N2,ξ
        0.0,                    // N3,ξ
        4.0 * eta,              // N4,ξ
        -4.0 * eta,             // N5,ξ
        4.0 * (1.0 - 2.0 * xi - eta),
    ];
    let deta = [
        4.0 * (xi + eta) - 3.0, // N1,η
        0.0,                    // N2,η
        4.0 * eta - 1.0,        // N3,η
        4.0 * xi,               // N4,η
        4.0 * (1.0 - xi - 2.0 * eta),
        -4.0 * xi, // N6,η
    ];
    (dxi, deta)
}

/// Assemble the 9-vectors `Hx,•` and `Hy,•` from quadratic shape-function
/// derivatives `dn` (Batoz compact DKT form; per-node DOF order (w, wx, wy)).
fn h_vectors(cf: &EdgeCoeff, dn: &[f64; 6]) -> ([f64; 9], [f64; 9]) {
    // Corner→edge pairing: node 1 uses edges 5,6; node 2 edges 6,4;
    // node 3 edges 4,5 (edge indices 0,1,2 = DKT edges 4,5,6).
    let mut hx = [0.0f64; 9];
    let mut hy = [0.0f64; 9];
    let corner_edges = [(1usize, 2usize), (2, 0), (0, 1)]; // (m, k) per node: edges m,k
    // Batoz's published tables use DOFs (w, θx = ∂w/∂y, θy = −∂w/∂x); this
    // crate's DOFs are slopes (w, wx, wy) = (w, −θy, θx), so Batoz's
    // θy-columns land in the wx slots NEGATED and his θx-columns land in
    // the wy slots unchanged. The resulting internal rotation fields are
    // βx = −w,x and βy = −w,y — both negated, so the curvature quadratic
    // form is unchanged (verified by the constant-curvature patch test).
    for node in 0..3 {
        let (m, k) = corner_edges[node]; // node 1: (5,6) etc.
        let nm = dn[3 + m];
        let nk = dn[3 + k];
        let ni = dn[node];
        hx[3 * node] = 1.5 * (cf.a[k] * nk - cf.a[m] * nm);
        hx[3 * node + 1] = cf.c[m] * nm + cf.c[k] * nk - ni;
        hx[3 * node + 2] = cf.b[m] * nm + cf.b[k] * nk;
        hy[3 * node] = 1.5 * (cf.d[k] * nk - cf.d[m] * nm);
        hy[3 * node + 1] = cf.b[m] * nm + cf.b[k] * nk;
        hy[3 * node + 2] = cf.e[m] * nm + cf.e[k] * nk - ni;
    }
    (hx, hy)
}

/// DKT bending stiffness (9×9 row-major, DOFs `(w, wx, wy)` per node) plus
/// the doubled area. Exact three-point mid-side integration.
///
/// # Errors
/// [`PlateError::DegenerateElement`] when |2A| falls below `1e-14` of the
/// squared characteristic edge length.
pub fn dkt_stiffness(
    x: &[f64; 3],
    y: &[f64; 3],
    d: &[f64; 9],
    element: usize,
) -> Result<([f64; 81], f64), PlateError> {
    let x21 = x[1] - x[0];
    let x31 = x[2] - x[0];
    let y21 = y[1] - y[0];
    let y31 = y[2] - y[0];
    let twice_area = x21 * y31 - x31 * y21;
    let scale = (x21 * x21 + y21 * y21).max(x31 * x31 + y31 * y31);
    if !twice_area.is_finite() || twice_area.abs() <= 1e-14 * scale {
        return Err(PlateError::DegenerateElement {
            element,
            twice_area,
        });
    }
    let cf = edge_coeffs(x, y);
    let area = 0.5 * twice_area;
    let gauss = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)];
    let mut k = [0.0f64; 81];
    for &(xi, eta) in &gauss {
        let (dxi, deta) = quad_shape_derivs(xi, eta);
        let (hx_xi, hy_xi) = h_vectors(&cf, &dxi);
        let (hx_eta, hy_eta) = h_vectors(&cf, &deta);
        // Cartesian curvature rows via the affine map inverse.
        let mut b = [[0.0f64; 9]; 3];
        for c in 0..9 {
            let bx_x = (y31 * hx_xi[c] - y21 * hx_eta[c]) / twice_area;
            let by_y = (-x31 * hy_xi[c] + x21 * hy_eta[c]) / twice_area;
            let bx_y = (-x31 * hx_xi[c] + x21 * hx_eta[c]) / twice_area;
            let by_x = (y31 * hy_xi[c] - y21 * hy_eta[c]) / twice_area;
            b[0][c] = bx_x;
            b[1][c] = by_y;
            b[2][c] = bx_y + by_x;
        }
        let w = area / 3.0;
        for r in 0..9 {
            for c in 0..9 {
                let mut acc = 0.0;
                for p in 0..3 {
                    for q in 0..3 {
                        acc += b[p][r] * d[p * 3 + q] * b[q][c];
                    }
                }
                k[r * 9 + c] += w * acc;
            }
        }
    }
    Ok((k, twice_area))
}

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------

/// Boundary condition applied to every listed node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeSupport {
    /// Simply supported: `w = 0`, slopes free.
    SimplySupported,
    /// Clamped: `w = wx = wy = 0`.
    Clamped,
}

/// A straight stiffener along a path of mesh nodes.
#[derive(Debug, Clone)]
pub struct Stiffener {
    /// Node path (consecutive nodes are beam elements).
    pub nodes: Vec<usize>,
    /// Young's modulus [Pa].
    pub e: f64,
    /// Shear modulus [Pa] (torsion).
    pub g: f64,
    /// Cross-section area [m²].
    pub area: f64,
    /// Second moment about the beam's own bending axis [m⁴].
    pub inertia: f64,
    /// Torsion constant J [m⁴].
    pub torsion: f64,
    /// Eccentricity of the beam centroid from the plate midplane [m]
    /// (parallel-axis: effective `EI + EAe²`).
    pub eccentricity: f64,
    /// Density [kg/m³] (lumped translational mass contribution).
    pub density: f64,
}

/// Assembled reduced pencil (Dirichlet DOFs eliminated) plus bookkeeping.
#[derive(Debug, Clone)]
pub struct PlateModel {
    /// Reduced stiffness (bending + stiffeners + optional prestress K_G).
    pub k: Csr,
    /// Reduced lumped mass (SPD).
    pub m: Csr,
    /// For each full DOF (3·node + comp), its reduced index or `None`.
    pub dof_map: Vec<Option<usize>>,
    /// Number of free DOFs.
    pub free: usize,
}

/// Assembly options.
#[derive(Debug, Clone, Copy)]
pub struct AssemblyOptions {
    /// Uniform membrane pre-tension T [N/m] (0 disables `K_G`).
    pub pretension: f64,
    /// Boundary support applied to the caller-supplied boundary nodes.
    pub support: EdgeSupport,
}

/// Assemble the reduced (K, M) pencil for a plate mesh: DKT bending +
/// lumped mass + optional membrane-prestress `K_G` + stiffener beams, with
/// the given boundary nodes eliminated per `support`.
///
/// # Errors
/// [`PlateError`] on degenerate elements, bad sections, or bad stiffeners.
#[allow(clippy::too_many_lines)] // one coherent assembly pipeline
pub fn assemble(
    mesh: &PlateMesh,
    section: &PlateSection,
    boundary: &[usize],
    stiffeners: &[Stiffener],
    opts: &AssemblyOptions,
) -> Result<PlateModel, PlateError> {
    let nn = mesh.node_count();
    let ndof = 3 * nn;
    // DOF elimination map.
    let mut fixed = vec![false; ndof];
    for &node in boundary {
        match opts.support {
            EdgeSupport::SimplySupported => fixed[3 * node] = true,
            EdgeSupport::Clamped => {
                fixed[3 * node] = true;
                fixed[3 * node + 1] = true;
                fixed[3 * node + 2] = true;
            }
        }
    }
    let mut dof_map = vec![None; ndof];
    let mut free = 0usize;
    for (dof, slot) in dof_map.iter_mut().enumerate() {
        if !fixed[dof] {
            *slot = Some(free);
            free += 1;
        }
    }

    let mut kc = Coo::new(free, free);
    let mut mc = Coo::new(free, free);
    let push_sym = |coo: &mut Coo, gr: usize, gc: usize, v: f64| {
        if let (Some(r), Some(c)) = (dof_map[gr], dof_map[gc]) {
            coo.push(r, c, v);
        }
    };

    for (eidx, tri) in mesh.tris.iter().enumerate() {
        let x = [
            mesh.nodes[tri[0]].0,
            mesh.nodes[tri[1]].0,
            mesh.nodes[tri[2]].0,
        ];
        let y = [
            mesh.nodes[tri[0]].1,
            mesh.nodes[tri[1]].1,
            mesh.nodes[tri[2]].1,
        ];
        let (ke, twice_area) = dkt_stiffness(&x, &y, &section.d, eidx)?;
        let area = 0.5 * twice_area.abs();
        let gdof = |local: usize| 3 * tri[local / 3] + local % 3;
        for r in 0..9 {
            for c in 0..9 {
                let v = ke[r * 9 + c];
                if v != 0.0 {
                    push_sym(&mut kc, gdof(r), gdof(c), v);
                }
            }
        }
        // Lumped mass: translational ρhA/3, rotary ρh³/12·A/3 per node.
        let h = section.thickness;
        let mt = section.density * h * area / 3.0;
        let mr = section.density * h * h * h / 12.0 * area / 3.0;
        for node in 0..3 {
            push_sym(&mut mc, gdof(3 * node), gdof(3 * node), mt);
            push_sym(&mut mc, gdof(3 * node + 1), gdof(3 * node + 1), mr);
            push_sym(&mut mc, gdof(3 * node + 2), gdof(3 * node + 2), mr);
        }
        // Membrane prestress: P1 geometric stiffness on w DOFs,
        // K_G = T·A·(∇N)ᵀ(∇N) with the linear hat gradients.
        if opts.pretension != 0.0 {
            let t = opts.pretension;
            let bg = [
                [
                    (y[1] - y[2]) / twice_area,
                    (y[2] - y[0]) / twice_area,
                    (y[0] - y[1]) / twice_area,
                ],
                [
                    (x[2] - x[1]) / twice_area,
                    (x[0] - x[2]) / twice_area,
                    (x[1] - x[0]) / twice_area,
                ],
            ];
            for r in 0..3 {
                for c in 0..3 {
                    let v = t * area * (bg[0][r] * bg[0][c] + bg[1][r] * bg[1][c]);
                    push_sym(&mut kc, gdof(3 * r), gdof(3 * c), v);
                }
            }
        }
    }

    // Stiffeners: Hermite bending on (w, slope-along) with EI + EAe²,
    // torsion GJ on the cross-slope, lumped translational mass.
    for st in stiffeners {
        if st.nodes.len() < 2 {
            return Err(PlateError::BadStiffener {
                what: "a stiffener needs at least two nodes",
            });
        }
        if !(st.e > 0.0 && st.g > 0.0 && st.area > 0.0 && st.inertia >= 0.0 && st.torsion >= 0.0) {
            return Err(PlateError::BadStiffener {
                what: "stiffener constants must be positive",
            });
        }
        for pair in st.nodes.windows(2) {
            let (n1, n2) = (pair[0], pair[1]);
            if n1 >= nn || n2 >= nn {
                return Err(PlateError::BadStiffener {
                    what: "stiffener node outside the mesh",
                });
            }
            let (x1, y1) = mesh.nodes[n1];
            let (x2, y2) = mesh.nodes[n2];
            let l = ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)).sqrt();
            if !l.is_finite() || l <= 0.0 {
                return Err(PlateError::BadStiffener {
                    what: "zero-length stiffener segment",
                });
            }
            let (tx, ty) = ((x2 - x1) / l, (y2 - y1) / l);
            let ei = st.e * (st.inertia + st.area * st.eccentricity * st.eccentricity);
            // Hermite 4-DOF [w1, s1, w2, s2] with s = slope along the beam
            // (= tx·wx + ty·wy).
            let kb = [
                [12.0, 6.0 * l, -12.0, 6.0 * l],
                [6.0 * l, 4.0 * l * l, -6.0 * l, 2.0 * l * l],
                [-12.0, -6.0 * l, 12.0, -6.0 * l],
                [6.0 * l, 2.0 * l * l, -6.0 * l, 4.0 * l * l],
            ];
            let coef = ei / (l * l * l);
            // Map beam DOF b → (global dof, weight) pairs.
            let maps: [Vec<(usize, f64)>; 4] = [
                vec![(3 * n1, 1.0)],
                vec![(3 * n1 + 1, tx), (3 * n1 + 2, ty)],
                vec![(3 * n2, 1.0)],
                vec![(3 * n2 + 1, tx), (3 * n2 + 2, ty)],
            ];
            for r in 0..4 {
                for c in 0..4 {
                    let v = coef * kb[r][c];
                    for &(gr, wr) in &maps[r] {
                        for &(gc, wc) in &maps[c] {
                            push_sym(&mut kc, gr, gc, v * wr * wc);
                        }
                    }
                }
            }
            // Torsion on the cross-slope q = −ty·wx + tx·wy.
            let gj = st.g * st.torsion / l;
            let tmaps: [Vec<(usize, f64)>; 2] = [
                vec![(3 * n1 + 1, -ty), (3 * n1 + 2, tx)],
                vec![(3 * n2 + 1, -ty), (3 * n2 + 2, tx)],
            ];
            let kt = [[1.0, -1.0], [-1.0, 1.0]];
            for r in 0..2 {
                for c in 0..2 {
                    for &(gr, wr) in &tmaps[r] {
                        for &(gc, wc) in &tmaps[c] {
                            push_sym(&mut kc, gr, gc, gj * kt[r][c] * wr * wc);
                        }
                    }
                }
            }
            // Lumped translational beam mass.
            let mb = st.density * st.area * l / 2.0;
            push_sym(&mut mc, 3 * n1, 3 * n1, mb);
            push_sym(&mut mc, 3 * n2, 3 * n2, mb);
        }
    }

    Ok(PlateModel {
        k: kc.assemble(),
        m: mc.assemble(),
        dof_map,
        free,
    })
}

/// Certified modes of an assembled plate model in the angular-frequency-
/// squared window `(low, high]` — a thin front over
/// [`fs_modal::slice_window`]; every returned mode carries its certified
/// eigenvalue interval, and the count is inertia-certified.
///
/// # Errors
/// [`PlateError::Modal`] forwarding the fs-modal refusal.
pub fn modes(
    model: &PlateModel,
    window: (f64, f64),
    opts: &SliceOptions,
) -> Result<SliceReport, PlateError> {
    Ok(fs_modal::slice_window(&model.k, &model.m, window, opts)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn steel_section() -> PlateSection {
        PlateSection::isotropic(200e9, 0.3, 0.002, 7800.0).expect("steel")
    }

    /// Analytic SS isotropic plate: ω_mn = π²(m²/a² + n²/b²)·√(D/ρh).
    fn ss_omega(d: f64, rho_h: f64, a: f64, b: f64, m: usize, n: usize) -> f64 {
        let pi = std::f64::consts::PI;
        pi * pi * ((m * m) as f64 / (a * a) + (n * n) as f64 / (b * b)) * (d / rho_h).sqrt()
    }

    #[test]
    fn element_rigid_body_and_constant_curvature_patch_tests() {
        // Deliberately irregular triangle.
        let x = [0.1, 1.3, 0.4];
        let y = [-0.2, 0.1, 0.9];
        let sec = steel_section();
        let (ke, twice_area) = dkt_stiffness(&x, &y, &sec.d, 0).expect("element");
        let area = 0.5 * twice_area;
        // Symmetry.
        for r in 0..9 {
            for c in 0..9 {
                let s = (ke[r * 9 + c] - ke[c * 9 + r]).abs();
                assert!(s < 1e-6 * ke[0].abs().max(1.0), "K symmetry at ({r},{c})");
            }
        }
        // Rigid-body motions: w = 1; w = x (wx = 1); w = y (wy = 1) — zero
        // strain energy to roundoff.
        let dscale = sec.d[0];
        for rb in 0..3 {
            let mut u = [0.0f64; 9];
            for node in 0..3 {
                match rb {
                    0 => u[3 * node] = 1.0,
                    1 => {
                        u[3 * node] = x[node];
                        u[3 * node + 1] = 1.0;
                    }
                    _ => {
                        u[3 * node] = y[node];
                        u[3 * node + 2] = 1.0;
                    }
                }
            }
            let mut energy = 0.0;
            for r in 0..9 {
                for c in 0..9 {
                    energy += u[r] * ke[r * 9 + c] * u[c];
                }
            }
            assert!(
                energy.abs() < 1e-9 * dscale,
                "rigid-body mode {rb} stores energy {energy:.3e}"
            );
        }
        // Constant curvature: w = ½(κ1x² + κ2y²) + κ12·x·y with consistent
        // slope DOFs. Energy must equal A·κᵀDκ EXACTLY (quadratic field,
        // exact integration) — the patch certificate that pins every table
        // and the DOF convention.
        let (k1, k2, k12) = (0.7, -0.4, 0.3);
        let mut u = [0.0f64; 9];
        for node in 0..3 {
            let (xi, yi) = (x[node], y[node]);
            u[3 * node] = f64::midpoint(k1 * xi * xi, k2 * yi * yi) + k12 * xi * yi;
            u[3 * node + 1] = k1 * xi + k12 * yi;
            u[3 * node + 2] = k2 * yi + k12 * xi;
        }
        let mut energy = 0.0;
        for r in 0..9 {
            for c in 0..9 {
                energy += u[r] * ke[r * 9 + c] * u[c];
            }
        }
        let kappa = [k1, k2, 2.0 * k12];
        let mut want = 0.0;
        for p in 0..3 {
            for q in 0..3 {
                want += kappa[p] * sec.d[p * 3 + q] * kappa[q];
            }
        }
        want *= area;
        assert!(
            (energy - want).abs() < 1e-9 * want.abs(),
            "constant-curvature energy {energy} vs exact {want}"
        );
        println!("{{\"suite\":\"fs-plate\",\"case\":\"patch-tests\",\"verdict\":\"pass\"}}");
    }

    #[test]
    fn ss_square_frequencies_converge_to_navier() {
        let sec = steel_section();
        let d = sec.d[0];
        let rho_h = sec.density * sec.thickness;
        let a = 1.0;
        let mut errs = Vec::new();
        for nx in [8usize, 16] {
            let mesh = PlateMesh::rectangle(a, a, nx, nx);
            let boundary = PlateMesh::rectangle_boundary(nx, nx);
            let model = assemble(
                &mesh,
                &sec,
                &boundary,
                &[],
                &AssemblyOptions {
                    pretension: 0.0,
                    support: EdgeSupport::SimplySupported,
                },
            )
            .expect("assemble");
            let w11 = ss_omega(d, rho_h, a, a, 1, 1);
            let lam11 = w11 * w11;
            let rep =
                modes(&model, (0.5 * lam11, 1.5 * lam11), &SliceOptions::default()).expect("modes");
            assert_eq!(rep.expected, 1, "isolated fundamental at nx={nx}");
            let got = rep.modes[0].lambda.sqrt();
            errs.push((got - w11).abs() / w11);
        }
        // h² (or better) trend: halving h must cut the error by ≥ ~3×.
        assert!(
            errs[1] < errs[0] / 3.0,
            "convergence trend too slow: {errs:?}"
        );
        assert!(errs[1] < 0.01, "fine-mesh error {:.3e} above 1%", errs[1]);
        println!(
            "{{\"suite\":\"fs-plate\",\"case\":\"ss-navier-convergence\",\"errors\":[{:.3e},{:.3e}],\"verdict\":\"pass\"}}",
            errs[0], errs[1]
        );
    }

    #[test]
    fn orthotropic_navier_and_axis_swap_mutation() {
        // Spruce-like constants (order-of-magnitude): E_L ≫ E_R.
        let law = OrthotropicElastic::new(
            [12.0e9, 0.9e9, 0.5e9],
            [0.35, 0.4, 0.4],
            [0.75e9, 0.6e9, 0.7e9],
            1.0,
        )
        .expect("law");
        let sec = PlateSection::orthotropic(&law, 0.003, 450.0).expect("section");
        let rho_h = sec.density * sec.thickness;
        let (a, b) = (1.0, 0.6);
        let pi = std::f64::consts::PI;
        let navier = |sec: &PlateSection, m: usize, n: usize| -> f64 {
            let (d11, d12, d22, d33) = (sec.d[0], sec.d[1], sec.d[4], sec.d[8]);
            let ma2 = (m as f64) / a * ((m as f64) / a);
            let nb2 = (n as f64) / b * ((n as f64) / b);
            let num = d11 * ma2 * ma2 + 2.0 * (d12 + 2.0 * d33) * ma2 * nb2 + d22 * nb2 * nb2;
            let pi2 = pi * pi;
            (pi2 * pi2 * num / rho_h).sqrt()
        };
        let nx = 20;
        let ny = 12;
        let mesh = PlateMesh::rectangle(a, b, nx, ny);
        let boundary = PlateMesh::rectangle_boundary(nx, ny);
        let model = assemble(
            &mesh,
            &sec,
            &boundary,
            &[],
            &AssemblyOptions {
                pretension: 0.0,
                support: EdgeSupport::SimplySupported,
            },
        )
        .expect("assemble");
        let w11 = navier(&sec, 1, 1);
        let rep = modes(
            &model,
            (0.7 * w11 * w11, 1.4 * w11 * w11),
            &SliceOptions::default(),
        )
        .expect("modes");
        assert_eq!(rep.expected, 1);
        let got = rep.modes[0].lambda.sqrt();
        assert!(
            (got - w11).abs() / w11 < 0.02,
            "orthotropic fundamental {got} vs Navier {w11}"
        );
        // MUTATION: swapping E_L and E_R must shift the fundamental
        // detectably (the wrong-grain-axis error class).
        let swapped = OrthotropicElastic::new(
            [0.9e9, 12.0e9, 0.5e9],
            [0.35 * 0.9 / 12.0, 0.4, 0.4], // keep symmetry admissible
            [0.75e9, 0.6e9, 0.7e9],
            1.0,
        )
        .expect("swapped law");
        let sec_sw = PlateSection::orthotropic(&swapped, 0.003, 450.0).expect("swapped section");
        let w11_sw = navier(&sec_sw, 1, 1);
        assert!(
            (w11_sw - w11).abs() / w11 > 0.10,
            "axis swap must move the fundamental by >10%, got {w11} vs {w11_sw}"
        );
        println!(
            "{{\"suite\":\"fs-plate\",\"case\":\"orthotropic-navier\",\"rel_err\":{:.3e},\"swap_shift\":{:.3},\"verdict\":\"pass\"}}",
            (got - w11).abs() / w11,
            (w11_sw - w11).abs() / w11
        );
    }

    #[test]
    fn drumhead_prestress_limit_matches_membrane() {
        // K_G-dominated: near-zero D, tension T. Continuum membrane:
        // ω_mn = π√((m/a)² + (n/b)²)·√(T/ρh).
        let sec = PlateSection::isotropic(1.0e3, 0.3, 0.001, 1000.0).expect("floppy");
        let rho_h = sec.density * sec.thickness;
        let t = 2000.0; // N/m
        let a = 1.0;
        let pi = std::f64::consts::PI;
        let w11 = pi * (2.0f64).sqrt() * (t / rho_h).sqrt() / a;
        let mut errs = Vec::new();
        for nx in [10usize, 20] {
            let mesh = PlateMesh::rectangle(a, a, nx, nx);
            let boundary = PlateMesh::rectangle_boundary(nx, nx);
            let model = assemble(
                &mesh,
                &sec,
                &boundary,
                &[],
                &AssemblyOptions {
                    pretension: t,
                    support: EdgeSupport::SimplySupported,
                },
            )
            .expect("assemble");
            let lam = w11 * w11;
            let rep =
                modes(&model, (0.5 * lam, 1.5 * lam), &SliceOptions::default()).expect("modes");
            assert_eq!(rep.expected, 1);
            errs.push((rep.modes[0].lambda.sqrt() - w11).abs() / w11);
        }
        assert!(errs[1] < errs[0] / 2.5, "membrane trend: {errs:?}");
        assert!(errs[1] < 0.02, "fine membrane error {:.3e}", errs[1]);
        println!(
            "{{\"suite\":\"fs-plate\",\"case\":\"drumhead-membrane-limit\",\"errors\":[{:.3e},{:.3e}],\"verdict\":\"pass\"}}",
            errs[0], errs[1]
        );
    }

    #[test]
    fn stiffener_raises_frequencies_and_eccentricity_matters() {
        let sec = steel_section();
        let a = 1.0;
        let nx = 12;
        let mesh = PlateMesh::rectangle(a, a, nx, nx);
        let boundary = PlateMesh::rectangle_boundary(nx, nx);
        // QUARTER-line brace: a midline brace sits on the braced plate's
        // new fundamental nodal line (w = w'' = 0 along it), where neither
        // EI nor eccentricity can act — physically correct and useless as a
        // discriminator. On the quarter line the fundamental bends the
        // beam, so both effects are observable.
        let quarter = nx / 4;
        let path: Vec<usize> = (0..=nx).map(|i| quarter * (nx + 1) + i).collect();
        // A DELIBERATELY WEAK 5 mm square steel brace (EI ≈ 10 N·m² vs
        // plate D ≈ 147 N·m): a stiff brace saturates into a rigid line
        // support — the mode stops bending it and extra rigidity becomes
        // invisible (measured: 13× EI moved the fundamental by 2e-5).
        // In the weak regime the parallel-axis term EAe² ≈ 500 N·m²
        // multiplies the beam's rigidity ~50× and must be clearly visible.
        let brace = |ecc: f64| Stiffener {
            nodes: path.clone(),
            e: 200e9,
            g: 80e9,
            area: 2.5e-5,
            inertia: 5.2e-11,
            torsion: 8.8e-11,
            eccentricity: ecc,
            density: 7800.0,
        };
        let fundamental = |stiffeners: &[Stiffener]| -> f64 {
            let model = assemble(
                &mesh,
                &sec,
                &boundary,
                stiffeners,
                &AssemblyOptions {
                    pretension: 0.0,
                    support: EdgeSupport::SimplySupported,
                },
            )
            .expect("assemble");
            let d = sec.d[0];
            let rho_h = sec.density * sec.thickness;
            let w11 = ss_omega(d, rho_h, a, a, 1, 1);
            // Generous window: a strong midline brace can push the
            // fundamental toward the (1,2)-type mode (~2.5× in ω ⇒ ~6.3×
            // in λ), so search up to 25·λ11 (5× in ω) and take the LOWEST
            // harvested mode.
            let rep = modes(
                &model,
                (0.25 * w11 * w11, 25.0 * w11 * w11),
                &SliceOptions::default(),
            )
            .expect("modes");
            assert!(rep.expected >= 1, "window must contain the fundamental");
            rep.modes[0].lambda.sqrt()
        };
        let bare = fundamental(&[]);
        let flat = fundamental(&[brace(0.0)]);
        let offset = fundamental(&[brace(0.01)]);
        // Rayleigh: adding a PSD beam stiffness cannot lower the
        // fundamental (its tiny lumped mass can, marginally — hence the
        // 0.999 floor), and the weak flat brace stiffens only slightly.
        assert!(flat > bare * 0.999, "flat brace sanity: {bare} → {flat}");
        assert!(
            offset > flat * 1.05,
            "eccentricity must add parallel-axis stiffness: bare {bare}, flat {flat} → offset {offset}"
        );
        println!(
            "{{\"suite\":\"fs-plate\",\"case\":\"stiffener-eccentricity\",\"bare\":{bare:.1},\"flat\":{flat:.1},\"offset\":{offset:.1},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn stiffener_terms_isolated_match_rayleigh_scale() {
        // Term isolation: each stiffener ingredient must act at the scale
        // first-order Rayleigh analysis predicts for the (1,1) mode with a
        // quarter-line brace (bending ≤ +1%, torsion ≤ +1%, mass ≤ −1%).
        // This is the gate that caught a sign error in the Hermite matrix
        // (an asymmetric K inflated the fundamental 64%).
        let sec = steel_section();
        let a = 1.0;
        let nx = 12;
        let mesh = PlateMesh::rectangle(a, a, nx, nx);
        let boundary = PlateMesh::rectangle_boundary(nx, nx);
        let quarter = nx / 4;
        let path: Vec<usize> = (0..=nx).map(|i| quarter * (nx + 1) + i).collect();
        let mk = |inertia: f64, torsion: f64, area: f64| Stiffener {
            nodes: path.clone(),
            e: 200e9,
            g: 80e9,
            area,
            inertia,
            torsion,
            eccentricity: 0.0,
            density: 7800.0,
        };
        let fundamental = |st: &[Stiffener]| -> f64 {
            let model = assemble(
                &mesh,
                &sec,
                &boundary,
                st,
                &AssemblyOptions {
                    pretension: 0.0,
                    support: EdgeSupport::SimplySupported,
                },
            )
            .expect("assemble");
            let d = sec.d[0];
            let rho_h = sec.density * sec.thickness;
            let w11 = ss_omega(d, rho_h, a, a, 1, 1);
            let rep = modes(
                &model,
                (0.25 * w11 * w11, 25.0 * w11 * w11),
                &SliceOptions::default(),
            )
            .expect("modes");
            rep.modes[0].lambda.sqrt()
        };
        let bare = fundamental(&[]);
        let bend_only = fundamental(&[mk(5.2e-11, 0.0, 1e-12)]);
        let tors_only = fundamental(&[mk(0.0, 8.8e-11, 1e-12)]);
        let mass_only = fundamental(&[mk(0.0, 0.0, 2.5e-5)]);
        assert!(
            bend_only > bare && bend_only < bare * 1.012,
            "bending term outside Rayleigh scale: {bare} → {bend_only}"
        );
        assert!(
            tors_only > bare && tors_only < bare * 1.012,
            "torsion term outside Rayleigh scale: {bare} → {tors_only}"
        );
        assert!(
            mass_only < bare && mass_only > bare * 0.988,
            "mass term outside Rayleigh scale: {bare} → {mass_only}"
        );
        println!(
            "{{\"suite\":\"fs-plate\",\"case\":\"stiffener-term-isolation\",\"bare\":{bare:.3},\"bend\":{bend_only:.3},\"torsion\":{tors_only:.3},\"mass\":{mass_only:.3},\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn clamped_square_matches_leissa() {
        // Leissa (NASA SP-160) clamped square: λ = ω·a²·√(ρh/D) = 35.992
        // for the fundamental. Exercises the Clamped BC elimination path.
        let sec = steel_section();
        let d = sec.d[0];
        let rho_h = sec.density * sec.thickness;
        let a = 1.0;
        let w1 = 35.992 / (a * a) * (d / rho_h).sqrt();
        let mut errs = Vec::new();
        for nx in [10usize, 20] {
            let mesh = PlateMesh::rectangle(a, a, nx, nx);
            let boundary = PlateMesh::rectangle_boundary(nx, nx);
            let model = assemble(
                &mesh,
                &sec,
                &boundary,
                &[],
                &AssemblyOptions {
                    pretension: 0.0,
                    support: EdgeSupport::Clamped,
                },
            )
            .expect("assemble");
            // Next clamped mode sits at 73.41/35.99 ≈ 2.04× in ω, so a
            // (0.5, 2.5)·λ1 window isolates the fundamental... in λ that
            // next mode is 4.16×λ1, comfortably outside.
            let rep = modes(
                &model,
                (0.5 * w1 * w1, 2.5 * w1 * w1),
                &SliceOptions::default(),
            )
            .expect("modes");
            assert_eq!(rep.expected, 1, "isolated clamped fundamental at nx={nx}");
            errs.push((rep.modes[0].lambda.sqrt() - w1).abs() / w1);
        }
        assert!(errs[1] < errs[0] / 2.5, "clamped trend: {errs:?}");
        assert!(errs[1] < 0.02, "clamped fine-mesh error {:.3e}", errs[1]);
        println!(
            "{{\"suite\":\"fs-plate\",\"case\":\"clamped-leissa\",\"errors\":[{:.3e},{:.3e}],\"verdict\":\"pass\"}}",
            errs[0], errs[1]
        );
    }

    /// Nodal-w samples of a mode shape on the sub-grid of every `step`-th
    /// node of a `rectangle(nx, ny)` mesh (fixed DOFs sample as 0).
    fn sample_w(model: &PlateModel, phi: &[f64], nx: usize, ny: usize, step: usize) -> Vec<f64> {
        let mut out = Vec::new();
        for j in (0..=ny).step_by(step) {
            for i in (0..=nx).step_by(step) {
                let node = j * (nx + 1) + i;
                out.push(match model.dof_map[3 * node] {
                    Some(r) => phi[r],
                    None => 0.0,
                });
            }
        }
        out
    }

    /// Modal assurance criterion between two shape samples.
    fn mac(a: &[f64], b: &[f64]) -> f64 {
        let dot: f64 = a.iter().zip(b).map(|(x, y)| x * y).sum();
        let na: f64 = a.iter().map(|x| x * x).sum();
        let nb: f64 = b.iter().map(|x| x * x).sum();
        dot * dot / (na * nb)
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one coherent literature-case gate
    fn olson_hazell_stiffened_panel_matches_literature() {
        // LITERATURE CASE — Olson & Hazell (1977), the classical clamped
        // square integrally rib-stiffened plate: 203 mm × 203 mm × 1.37 mm,
        // one central rib 6.35 mm wide × 12.7 mm deep, E = 68.7 GPa,
        // ν = 0.3, ρ = 2820 kg/m³. Reference frequencies (Hz) as reproduced
        // in two independent secondary sources (Srivastava, Datta & Sheikh,
        // Shock and Vibration 11 (2004) 9–19, Table 4; Thinh, Binh & Tu,
        // Vietnam J. Mech. 35 (2013) 31–50, Table 1):
        let oh_theory = [718.1, 751.4, 997.4, 1007.1, 1419.8];
        let oh_experiment = [689.0, 725.0, 961.0, 986.0, 1376.0];
        let a = 0.203;
        let h = 1.37e-3;
        let e_al = 68.7e9;
        let rho = 2820.0;
        let sec = PlateSection::isotropic(e_al, 0.3, h, rho).expect("aluminium section");
        // Rib: rectangle 6.35 × 12.7 mm on one face ⇒ centroid offset from
        // the plate midplane e = h/2 + d/2. Torsion constant β·a·b³ with
        // a/b = 2 ⇒ β ≈ 0.229 (Roark). Full-composite parallel-axis model
        // (EI + EAe²) — the crate's stated eccentricity model.
        let (bw, bd) = (6.35e-3, 12.7e-3);
        let rib = |nodes: Vec<usize>| Stiffener {
            nodes,
            e: e_al,
            g: e_al / 2.6,
            area: bw * bd,
            inertia: bw * bd * bd * bd / 12.0,
            torsion: 0.229 * bd * bw * bw * bw,
            eccentricity: f64::midpoint(h, bd),
            density: rho,
        };
        // Two mesh densities; the coarse node grid is a subset of the fine
        // one (24 = 2·12), so MAC pairing works on shared nodal w samples.
        let (nx_c, nx_f) = (12usize, 24usize);
        let window = {
            let lo = 2.0 * std::f64::consts::PI * 400.0;
            let hi = 2.0 * std::f64::consts::PI * 1700.0;
            (lo * lo, hi * hi)
        };
        let run = |nx: usize| -> (Vec<f64>, Vec<Vec<f64>>) {
            let mesh = PlateMesh::rectangle(a, a, nx, nx);
            let boundary = PlateMesh::rectangle_boundary(nx, nx);
            let mid = nx / 2;
            let path: Vec<usize> = (0..=nx).map(|i| mid * (nx + 1) + i).collect();
            let model = assemble(
                &mesh,
                &sec,
                &boundary,
                &[rib(path)],
                &AssemblyOptions {
                    pretension: 0.0,
                    support: EdgeSupport::Clamped,
                },
            )
            .expect("assemble");
            let rep = modes(&model, window, &SliceOptions::default()).expect("modes");
            assert!(
                rep.expected >= 5,
                "window must certify at least five modes at nx={nx}, got {}",
                rep.expected
            );
            let freqs = rep
                .modes
                .iter()
                .map(|m| m.lambda.sqrt() / (2.0 * std::f64::consts::PI))
                .collect();
            let shapes = rep
                .modes
                .iter()
                .map(|m| sample_w(&model, &m.phi, nx, nx, nx / nx_c))
                .collect();
            (freqs, shapes)
        };
        let (freq_c, shape_c) = run(nx_c);
        let (freq_f, shape_f) = run(nx_f);
        // MAC pairing table: pair each of the first five fine-mesh modes
        // with its best coarse-mesh partner. The pairing must be
        // order-preserving (no mode crossing between densities) and strong
        // (MAC ≥ 0.90) — this is what certifies that the frequency rows
        // below compare like with like, not accidental neighbors.
        let mut rows = String::new();
        let mut last_pair = None;
        for (m, sf) in shape_f.iter().take(5).enumerate() {
            let (mut best, mut best_mac) = (0usize, -1.0f64);
            for (c, sc) in shape_c.iter().enumerate() {
                let v = mac(sf, sc);
                if v > best_mac {
                    best = c;
                    best_mac = v;
                }
            }
            assert!(
                best_mac >= 0.90,
                "mode {m}: best coarse MAC {best_mac:.3} below 0.90"
            );
            if let Some(prev) = last_pair {
                assert!(
                    best > prev,
                    "MAC pairing must be order-preserving: fine mode {m} pairs to \
                     coarse {best} after {prev}"
                );
            }
            last_pair = Some(best);
            let rel_theory = (freq_f[m] - oh_theory[m]).abs() / oh_theory[m];
            let rel_exp = (freq_f[m] - oh_experiment[m]).abs() / oh_experiment[m];
            // Authored acceptance: 8% of the Olson–Hazell THEORY column.
            // The model boundary earns the envelope: bending-only DKT with
            // a full-composite parallel-axis rib (no membrane DOFs) plus
            // lumped translational rib mass — the classical eccentric-beam
            // idealization, not the papers' in-plane-coupled elements.
            assert!(
                rel_theory < 0.08,
                "mode {m}: {:.1} Hz vs Olson–Hazell theory {:.1} Hz ({:.1}%)",
                freq_f[m],
                oh_theory[m],
                100.0 * rel_theory
            );
            {
                use core::fmt::Write as _;
                write!(
                    rows,
                    "{}{{\"mode\":{},\"freq_hz\":{:.1},\"coarse_hz\":{:.1},\"mac\":{:.4},\
                     \"oh_theory_hz\":{:.1},\"oh_experiment_hz\":{:.1},\
                     \"rel_err_theory\":{:.4},\"rel_err_experiment\":{:.4}}}",
                    if m == 0 { "" } else { "," },
                    m + 1,
                    freq_f[m],
                    freq_c[best],
                    best_mac,
                    oh_theory[m],
                    oh_experiment[m],
                    rel_theory,
                    rel_exp
                )
                .expect("write to String cannot fail");
            }
        }
        println!(
            "{{\"suite\":\"fs-plate\",\"case\":\"olson-hazell-stiffened-panel\",\
             \"source\":\"Olson & Hazell 1977 via Srivastava 2004 Table 4 + Thinh 2013 Table 1\",\
             \"plate\":\"203x203x1.37mm clamped\",\"rib\":\"6.35x12.7mm central, e=7.035mm\",\
             \"pairs\":[{rows}],\"verdict\":\"pass\"}}"
        );
    }

    #[test]
    fn named_refusals_fire() {
        // Degenerate element.
        let x = [0.0, 1.0, 2.0];
        let y = [0.0, 0.0, 0.0];
        let sec = steel_section();
        assert!(matches!(
            dkt_stiffness(&x, &y, &sec.d, 7),
            Err(PlateError::DegenerateElement { element: 7, .. })
        ));
        // Bad section.
        assert!(matches!(
            PlateSection::isotropic(200e9, 0.3, -1.0, 7800.0),
            Err(PlateError::BadSection { .. })
        ));
        // Bad stiffener.
        let mesh = PlateMesh::rectangle(1.0, 1.0, 2, 2);
        let boundary = PlateMesh::rectangle_boundary(2, 2);
        let bad = Stiffener {
            nodes: vec![0, 99],
            e: 1.0,
            g: 1.0,
            area: 1.0,
            inertia: 1.0,
            torsion: 1.0,
            eccentricity: 0.0,
            density: 1.0,
        };
        let err = assemble(
            &mesh,
            &sec,
            &boundary,
            &[bad],
            &AssemblyOptions {
                pretension: 0.0,
                support: EdgeSupport::SimplySupported,
            },
        )
        .unwrap_err();
        assert!(err.to_string().contains("FS-PLATE-BAD-STIFFENER"));
    }
}
