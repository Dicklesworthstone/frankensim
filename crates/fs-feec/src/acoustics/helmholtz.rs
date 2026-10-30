//! Interior Helmholtz finite elements on tetrahedral complexes (music bead `frankensim-irp-v5-instrument-realism-pexfp.93`,
//! M2.4-W1a, in the home M2.4-W0 chose: fs-feec).
//!
//! Equations (time convention `e^{+j w t}`, `k = w / c`, outward normal `n`). The pressure `p` solves
//! `laplacian p + k^2 p = 0` in the cavity. A boundary is one of
//! - rigid (`dp/dn = 0`), the natural condition, nothing to add;
//! - a specified inward particle velocity `v_in`: `dp/dn = +j w rho v_in`;
//! - a locally reacting impedance `Z`: `p = Z v_out`, i.e. `dp/dn = -j k y p` with the dimensionless specific admittance
//!   `y = rho c / Z`.
//!
//! Weak form on the hierarchical H1 space of [`SimplexSpace`] (order `r`, straight-sided tets):
//! `(K - k^2 M + j k sum_b y_b C_b) p = j k rho c v_in b_src`, with the stiffness `K`, the mass `M`, the boundary mass matrices
//! `C_b = int_b phi_i phi_j dS` and the load vectors `b_b = int_b phi_i dS`. Outward volume flow through a loaded boundary is
//! `U_b = (y_b / (rho c)) b_b . p`; the volume velocity transfer function between two boundaries is the ratio of an outward flow
//! to the inflow `v_in * area(src)`.
//!
//! Two paths. [`HelmholtzFe::modes`] gives the rigid-wall (Neumann) eigenpairs of `(K, M)` by `fs_modal::slice_window`, with the
//! inertia-certified count and per-mode residuals. [`ReducedHelmholtz`] is a Galerkin projection of the DRIVEN problem onto a basis made of
//! those modes plus one static (residual-flexibility) vector per boundary `(K + M)^-1 b_b`, M-orthonormalized: a frequency sweep is then
//! a dense `n x n` complex solve per frequency. The static vectors matter: a boundary-forced expansion in Neumann modes alone converges
//! only like `1/N`, because the loaded boundary is where the modes are least able to represent the field.
//!
//! Limits, stated: straight-sided elements (no isoparametric map, so a curved wall has an `O(h^2)` geometric floor), first-order
//! viscothermal boundary layers only inside their explicit thin-layer band, a uniform impedance law per boundary, one-way coupling (an admittance is an
//! input law, not a radiation solve). Refusals: order above [`MAX_ORDER`], a degenerate tet, a face that is not a boundary face,
//! a window holding more than [`MAX_MODES`] modes, a reduced basis above [`MAX_REDUCED_DIM`]; never a silent truncation.

use crate::highorder::quad1d::gauss_legendre;
use crate::highorder::simplex::SimplexSpace;
use fs_modal::{ModalError, SliceOptions, slice_window};
use fs_rep_mesh::TetComplex;
use fs_sparse::direct::{DirectOrdering, LdltError, LdltOptions, SymbolicLdlt};
use fs_sparse::{Coo, Csr};
use std::fmt;

/// Heat-diffusion modes on this same hierarchical volume space.
pub mod diffusion;

/// Viscous and thermal boundary layers on the original pressure traces.
pub mod boundary_layer;

/// Work-conjugate spatial port modes and multiport extraction.
pub mod multiport;

/// Highest element order accepted; `MAX_ORDER + 1` refuses.
pub const MAX_ORDER: usize = 4;
/// Most Neumann modes a window may hold; more refuses (the window is never truncated).
pub const MAX_MODES: usize = 512;
/// Largest reduced basis; more refuses.
pub const MAX_REDUCED_DIM: usize = 640;

/// A homogeneous gas, as read from a gas state by the caller (fs-feec stays below the material layer).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Medium {
    /// Density rho [kg/m^3].
    pub density_kg_m3: f64,
    /// Speed of sound c [m/s].
    pub sound_speed_m_s: f64,
}

/// Typed refusals of the interior Helmholtz FE.
#[derive(Debug, Clone, PartialEq)]
pub enum AcousticFeError {
    /// Density or sound speed is not finite and positive.
    InvalidMedium {
        /// Which field.
        what: &'static str,
    },
    /// Element order outside `1..=MAX_ORDER`.
    OrderOutOfRange {
        /// Requested order.
        order: usize,
        /// The cap.
        cap: usize,
    },
    /// The position table does not match the vertex count, or holds a non-finite coordinate.
    InvalidPositions {
        /// Vertex count of the complex.
        vertices: usize,
        /// Positions supplied.
        positions: usize,
    },
    /// A tet with zero or non-finite volume.
    DegenerateTet {
        /// Tet index.
        tet: usize,
    },
    /// A face given as a boundary is not a face of the complex.
    FaceNotInMesh {
        /// The face as given.
        face: [u32; 3],
    },
    /// A face given as a boundary is shared by two tets.
    FaceNotOnBoundary {
        /// The face as given.
        face: [u32; 3],
    },
    /// The window holds more modes than [`MAX_MODES`].
    TooManyModes {
        /// Modes found.
        found: usize,
        /// The cap.
        cap: usize,
    },
    /// The reduced basis would exceed [`MAX_REDUCED_DIM`].
    ReducedBasisTooLarge {
        /// Basis size.
        size: usize,
        /// The cap.
        cap: usize,
    },
    /// An input of the driven problem is inconsistent.
    InvalidInput {
        /// What is wrong.
        what: &'static str,
    },
    /// The reduced dense system is singular at this frequency.
    SingularAtFrequency {
        /// Frequency [Hz].
        frequency_hz: f64,
    },
    /// The modal facility refused.
    Modal(ModalError),
    /// The sparse direct solver refused.
    Factor(LdltError),
}

impl fmt::Display for AcousticFeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FS-FEEC-ACOUSTICS: {self:?}")
    }
}

impl std::error::Error for AcousticFeError {}

impl From<ModalError> for AcousticFeError {
    fn from(e: ModalError) -> Self {
        AcousticFeError::Modal(e)
    }
}

impl From<LdltError> for AcousticFeError {
    fn from(e: LdltError) -> Self {
        AcousticFeError::Factor(e)
    }
}

/// Rigid-wall eigenpairs of a cavity.
#[derive(Debug, Clone)]
pub struct HelmholtzModes {
    /// Eigenvalues `k^2` [1/m^2], ascending.
    pub lambda: Vec<f64>,
    /// Eigenfrequencies `c sqrt(max(lambda, 0)) / (2 pi)` [Hz], ascending.
    pub frequency_hz: Vec<f64>,
    /// Mass-normalized mode shapes (`phi^T M phi = 1`), one dof vector per mode.
    pub shapes: Vec<Vec<f64>>,
    /// Certified residual `||K phi - lambda M phi||_{M^-1}` of each mode.
    pub residual: Vec<f64>,
    /// Inertia-certified eigenvalue count strictly below the window's low end.
    pub below_low: usize,
    /// Inertia-certified eigenvalue count strictly below the window's high end.
    pub below_high: usize,
}

/// The interior Helmholtz finite-element space on one tet complex.
pub struct HelmholtzFe<'c> {
    space: SimplexSpace<'c>,
    positions: Vec<[f64; 3]>,
    medium: Medium,
    order: usize,
    stiffness: Csr,
    mass: Csr,
    /// Per canonical face: (first owning tet, its omitted local vertex, owner count).
    owners: Vec<(u32, u8, u8)>,
}

struct FaceRule {
    vertices: [u32; 3],
    dofs: Vec<usize>,
    phi: Vec<Vec<f64>>,
    weights: Vec<f64>,
    area: f64,
}

/// Triangle quadrature by the collapsed (Duffy) Gauss-Legendre map: barycentric points and weights summing to 1/2.
fn tri_quadrature(n: usize) -> Vec<([f64; 3], f64)> {
    let (nodes, weights) = gauss_legendre(n);
    let map = |x: f64| f64::midpoint(1.0, x);
    let mut out = Vec::with_capacity(n * n);
    for (&xu, &wu) in nodes.iter().zip(&weights) {
        let u = map(xu);
        for (&xv, &wv) in nodes.iter().zip(&weights) {
            let v = map(xv);
            let l1 = u;
            let l2 = v * (1.0 - u);
            out.push(([1.0 - l1 - l2, l1, l2], wu * wv * (1.0 - u) / 4.0));
        }
    }
    out
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn csr_sum(a: &Csr, b: &Csr) -> Csr {
    let mut coo = Coo::new(a.nrows(), a.ncols());
    for r in 0..a.nrows() {
        for m in [a, b] {
            let (cols, vals) = m.row(r);
            for (&c, &v) in cols.iter().zip(vals) {
                coo.push(r, c, v);
            }
        }
    }
    coo.assemble()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn spmv_new(a: &Csr, x: &[f64]) -> Vec<f64> {
    let mut y = vec![0.0; a.nrows()];
    a.spmv(x, &mut y);
    y
}

impl<'c> HelmholtzFe<'c> {
    /// Assemble the stiffness and mass of order `order` on `complex` with vertex `positions` [m].
    ///
    /// # Errors
    /// An invalid medium, an order outside `1..=MAX_ORDER`, a bad position table or a degenerate tet.
    pub fn new(
        complex: &'c TetComplex,
        positions: &[[f64; 3]],
        order: usize,
        medium: Medium,
    ) -> Result<Self, AcousticFeError> {
        if !(medium.density_kg_m3.is_finite() && medium.density_kg_m3 > 0.0) {
            return Err(AcousticFeError::InvalidMedium { what: "density" });
        }
        if !(medium.sound_speed_m_s.is_finite() && medium.sound_speed_m_s > 0.0) {
            return Err(AcousticFeError::InvalidMedium {
                what: "sound speed",
            });
        }
        if !(1..=MAX_ORDER).contains(&order) {
            return Err(AcousticFeError::OrderOutOfRange {
                order,
                cap: MAX_ORDER,
            });
        }
        if positions.len() != complex.vertex_count
            || positions.iter().any(|p| p.iter().any(|x| !x.is_finite()))
        {
            return Err(AcousticFeError::InvalidPositions {
                vertices: complex.vertex_count,
                positions: positions.len(),
            });
        }
        for (t, tet) in complex.tets.iter().enumerate() {
            let p0 = positions[tet[0] as usize];
            let six_v = dot3(
                sub3(positions[tet[1] as usize], p0),
                cross3(
                    sub3(positions[tet[2] as usize], p0),
                    sub3(positions[tet[3] as usize], p0),
                ),
            );
            if !(six_v.is_finite() && six_v != 0.0) {
                return Err(AcousticFeError::DegenerateTet { tet: t });
            }
        }
        let space = SimplexSpace::new(complex, order);
        let stiffness = space.stiffness(positions);
        let mass = space.mass(positions);
        let mut owners = vec![(0u32, 0u8, 0u8); complex.faces.len()];
        for (t, tet) in complex.tets.iter().enumerate() {
            for omit in 0..4usize {
                let mut key: Vec<u32> = (0..4).filter(|&l| l != omit).map(|l| tet[l]).collect();
                key.sort_unstable();
                let key = [key[0], key[1], key[2]];
                let fi = complex
                    .faces
                    .binary_search(&key)
                    .expect("every tet face is in the canonical table");
                if owners[fi].2 == 0 {
                    owners[fi] = (t as u32, omit as u8, 1);
                } else {
                    owners[fi].2 = owners[fi].2.saturating_add(1);
                }
            }
        }
        Ok(Self {
            space,
            positions: positions.to_vec(),
            medium,
            order,
            stiffness,
            mass,
            owners,
        })
    }

    /// Number of dofs.
    #[must_use]
    pub fn ndof(&self) -> usize {
        self.space.ndof
    }

    /// The stiffness `int grad phi_i . grad phi_j`.
    #[must_use]
    pub fn stiffness(&self) -> &Csr {
        &self.stiffness
    }

    /// The mass `int phi_i phi_j`.
    #[must_use]
    pub fn mass(&self) -> &Csr {
        &self.mass
    }

    /// The medium.
    #[must_use]
    pub const fn medium(&self) -> Medium {
        self.medium
    }

    /// Canonical faces owned by exactly one tet (the boundary), in table order.
    #[must_use]
    pub fn boundary_faces(&self) -> Vec<[u32; 3]> {
        self.space
            .complex
            .faces
            .iter()
            .zip(&self.owners)
            .filter(|(_, o)| o.2 == 1)
            .map(|(f, _)| *f)
            .collect()
    }

    fn face_rule(
        &self,
        face: [u32; 3],
        quad: &[([f64; 3], f64)],
    ) -> Result<FaceRule, AcousticFeError> {
        let mut key = face;
        key.sort_unstable();
        let fi = self
            .space
            .complex
            .faces
            .binary_search(&key)
            .map_err(|_| AcousticFeError::FaceNotInMesh { face })?;
        let (t, omit, count) = self.owners[fi];
        if count != 1 {
            return Err(AcousticFeError::FaceNotOnBoundary { face });
        }
        let tet = self.space.complex.tets[t as usize];
        let local: Vec<usize> = (0..4).filter(|&l| l != omit as usize).collect();
        let pa = self.positions[tet[local[0]] as usize];
        let pb = self.positions[tet[local[1]] as usize];
        let pc = self.positions[tet[local[2]] as usize];
        let cr = cross3(sub3(pb, pa), sub3(pc, pa));
        let area = 0.5 * dot3(cr, cr).sqrt();
        let mut dofs = Vec::new();
        let mut phi = Vec::new();
        for (global, f) in self.space.element_dofs(t as usize) {
            let vals: Vec<f64> = quad
                .iter()
                .map(|(b, _)| {
                    let mut lam = [0.0; 4];
                    lam[local[0]] = b[0];
                    lam[local[1]] = b[1];
                    lam[local[2]] = b[2];
                    f.eval(lam, self.order)
                })
                .collect();
            if vals.iter().any(|v| v.abs() > 1e-14) {
                dofs.push(global);
                phi.push(vals);
            }
        }
        let weights = quad.iter().map(|(_, w)| w * 2.0 * area).collect();
        Ok(FaceRule {
            vertices: [tet[local[0]], tet[local[1]], tet[local[2]]],
            dofs,
            phi,
            weights,
            area,
        })
    }

    /// Number of exterior boundary faces, without allocating a face table.
    #[must_use]
    pub fn boundary_face_count(&self) -> usize {
        self.owners.iter().filter(|owner| owner.2 == 1).count()
    }

    /// Original straight-sided volume geometry; no boundary remeshing is implied.
    #[must_use]
    pub fn positions(&self) -> &[[f64; 3]] { &self.positions }

    /// Original tetrahedra in their supplied vertex order.
    #[must_use]
    pub fn tetrahedra(&self) -> &[[u32; 4]] { &self.space.complex.tets }

    /// Original hierarchical polynomial order.
    #[must_use]
    pub const fn order(&self) -> usize { self.order }

    /// One boundary triangle in the caller's vertex order and its OUTWARD
    /// unit normal. The owning tetrahedron determines the sign, independently
    /// of how a caller ordered the three boundary indices.
    pub fn boundary_geometry(&self, face: [u32; 3])
        -> Result<([[f64; 3]; 3], [f64; 3], f64), AcousticFeError>
    {
        let mut key = face;
        key.sort_unstable();
        let fi = self.space.complex.faces.binary_search(&key)
            .map_err(|_| AcousticFeError::FaceNotInMesh { face })?;
        let (tet, omit, count) = self.owners[fi];
        if count != 1 { return Err(AcousticFeError::FaceNotOnBoundary { face }); }
        let vertices = face.map(|v| self.positions[v as usize]);
        let area_vector = cross3(sub3(vertices[1], vertices[0]), sub3(vertices[2], vertices[0]));
        let twice_area = fs_math::det::sqrt(dot3(area_vector, area_vector));
        let opposite = self.positions[self.space.complex.tets[tet as usize][omit as usize] as usize];
        let sign = if dot3(area_vector, sub3(opposite, vertices[0])) < 0.0 { 1.0 } else { -1.0 };
        if !twice_area.is_finite() || twice_area <= 0.0 {
            return Err(AcousticFeError::DegenerateTet { tet: tet as usize });
        }
        Ok((vertices, area_vector.map(|v| sign * v / twice_area), 0.5 * twice_area))
    }

    /// Work-conjugate boundary load b_i = integral(phi_i u_n dS) for
    /// piecewise-linear normal displacement/velocity, using the SAME pressure
    /// basis and face quadrature as the original Helmholtz operators.
    ///
    /// Each entry contains a boundary face and values at its three vertices in
    /// THAT order. This integrates the product, not the product of face means:
    /// nonuniform structural traces therefore couple to nonuniform pressure.
    /// No normal or sign is inferred here; the structural owner supplies an
    /// outward-from-gas trace. Repeated faces and nonfinite traces refuse.
    pub fn boundary_linear_load(&self, traces: &[([u32; 3], [f64; 3])])
        -> Result<Vec<f64>, AcousticFeError>
    {
        let quad = tri_quadrature(self.order + 1);
        let mut seen = std::collections::BTreeSet::new();
        let mut load = vec![0.0; self.ndof()];
        for &(face, values) in traces {
            let mut key = face;
            key.sort_unstable();
            if !seen.insert(key) || values.iter().any(|v| !v.is_finite()) {
                return Err(AcousticFeError::InvalidInput {
                    what: "boundary trace requires unique faces and finite nodal values",
                });
            }
            let rule = self.face_rule(face, &quad)?;
            let mut local_values = [0.0; 3];
            for (local, vertex) in rule.vertices.iter().enumerate() {
                let source = face.iter().position(|v| v == vertex)
                    .ok_or(AcousticFeError::InvalidInput { what: "boundary trace vertex mismatch" })?;
                local_values[local] = values[source];
            }
            for (local, &global) in rule.dofs.iter().enumerate() {
                let value: f64 = quad.iter().zip(&rule.weights).zip(&rule.phi[local])
                    .map(|(((barycentric, _), weight), pressure)| {
                        let normal: f64 = barycentric.iter().zip(local_values).map(|(b,u)| b*u).sum();
                        weight * pressure * normal
                    }).sum();
                load[global] += value;
            }
        }
        if load.iter().any(|v| !v.is_finite()) {
            return Err(AcousticFeError::InvalidInput { what: "boundary trace projection overflow" });
        }
        Ok(load)
    }

    /// The boundary mass matrix `C = int_faces phi_i phi_j dS` of the given boundary faces.
    ///
    /// # Errors
    /// A face that is not a boundary face of the complex.
    pub fn boundary_mass(&self, faces: &[[u32; 3]]) -> Result<Csr, AcousticFeError> {
        let quad = tri_quadrature(self.order + 1);
        let mut coo = Coo::new(self.ndof(), self.ndof());
        for &face in faces {
            let rule = self.face_rule(face, &quad)?;
            for (i, &gi) in rule.dofs.iter().enumerate() {
                for (j, &gj) in rule.dofs.iter().enumerate() {
                    let s: f64 = rule
                        .weights
                        .iter()
                        .zip(rule.phi[i].iter().zip(&rule.phi[j]))
                        .map(|(w, (a, b))| w * a * b)
                        .sum();
                    coo.push(gi, gj, s);
                }
            }
        }
        Ok(coo.assemble())
    }

    /// The load vector `b_i = int_faces phi_i dS` and the total area of the given boundary faces.
    ///
    /// # Errors
    /// A face that is not a boundary face of the complex.
    pub fn boundary_load(&self, faces: &[[u32; 3]]) -> Result<(Vec<f64>, f64), AcousticFeError> {
        let quad = tri_quadrature(self.order + 1);
        let mut b = vec![0.0; self.ndof()];
        let mut area = 0.0;
        for &face in faces {
            let rule = self.face_rule(face, &quad)?;
            area += rule.area;
            for (i, &gi) in rule.dofs.iter().enumerate() {
                let s: f64 = rule
                    .weights
                    .iter()
                    .zip(&rule.phi[i])
                    .map(|(w, a)| w * a)
                    .sum();
                b[gi] += s;
            }
        }
        Ok((b, area))
    }

    /// Rigid-wall eigenpairs with `window_hz.0 < f <= window_hz.1` (a negative low end includes the constant mode),
    /// with the inertia-certified count; a window holding more than [`MAX_MODES`] modes refuses.
    ///
    /// # Errors
    /// A malformed window, a modal-facility refusal, or too many modes.
    pub fn modes(
        &self,
        window_hz: (f64, f64),
        opts: &SliceOptions,
    ) -> Result<HelmholtzModes, AcousticFeError> {
        let c = self.medium.sound_speed_m_s;
        let to_lambda = |f: f64| {
            let k = 2.0 * core::f64::consts::PI * f / c;
            k * k.abs()
        };
        let report = slice_window(
            &self.stiffness,
            &self.mass,
            (to_lambda(window_hz.0), to_lambda(window_hz.1)),
            opts,
        )?;
        if report.modes.len() > MAX_MODES {
            return Err(AcousticFeError::TooManyModes {
                found: report.modes.len(),
                cap: MAX_MODES,
            });
        }
        let frequency_hz = report
            .modes
            .iter()
            .map(|m| c * m.lambda.max(0.0).sqrt() / (2.0 * core::f64::consts::PI))
            .collect();
        Ok(HelmholtzModes {
            lambda: report.modes.iter().map(|m| m.lambda).collect(),
            frequency_hz,
            residual: report.modes.iter().map(|m| m.residual).collect(),
            shapes: report.modes.into_iter().map(|m| m.phi).collect(),
            below_low: report.below_low,
            below_high: report.below_high,
        })
    }
}

/// A Galerkin reduction of the driven Helmholtz problem for a set of boundary patches (ports).
pub struct ReducedHelmholtz {
    medium: Medium,
    n: usize,
    k_r: Vec<f64>,
    m_r: Vec<f64>,
    c_r: Vec<Vec<f64>>,
    b_r: Vec<Vec<f64>>,
    area: Vec<f64>,
    boundary_layer: Option<boundary_layer::ReducedBoundaryLayer>,
}

/// The solution of the reduced driven problem at one frequency.
#[derive(Debug, Clone)]
pub struct ReducedResponse {
    /// Complex reduced coordinates `(re, im)`.
    pub amplitudes: Vec<[f64; 2]>,
    /// Outward volume flow of each port [m^3/s], complex `(re, im)`; the source port carries `-v_in * area`.
    pub outward_flow: Vec<[f64; 2]>,
    /// Mean pressure over each port [Pa], complex.
    pub mean_pressure: Vec<[f64; 2]>,
    /// Inflow `v_in * area(source)` [m^3/s].
    pub inflow: f64,
}

fn cmul(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

fn cdiv(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    let d = b[0] * b[0] + b[1] * b[1];
    [
        (a[0] * b[0] + a[1] * b[1]) / d,
        (a[1] * b[0] - a[0] * b[1]) / d,
    ]
}

/// Dense complex Gaussian elimination with partial pivoting; `None` when a pivot vanishes.
fn solve_complex(mut a: Vec<[f64; 2]>, mut rhs: Vec<[f64; 2]>, n: usize) -> Option<Vec<[f64; 2]>> {
    for col in 0..n {
        let mut piv = col;
        let mut best = 0.0;
        for r in col..n {
            let m = a[r * n + col][0].hypot(a[r * n + col][1]);
            if m > best {
                best = m;
                piv = r;
            }
        }
        if !(best.is_finite() && best > 0.0) {
            return None;
        }
        if piv != col {
            for c in 0..n {
                a.swap(piv * n + c, col * n + c);
            }
            rhs.swap(piv, col);
        }
        let d = a[col * n + col];
        for r in (col + 1)..n {
            let f = cdiv(a[r * n + col], d);
            if f == [0.0, 0.0] {
                continue;
            }
            for c in col..n {
                let t = cmul(f, a[col * n + c]);
                a[r * n + c][0] -= t[0];
                a[r * n + c][1] -= t[1];
            }
            let t = cmul(f, rhs[col]);
            rhs[r][0] -= t[0];
            rhs[r][1] -= t[1];
        }
    }
    let mut x = vec![[0.0, 0.0]; n];
    for r in (0..n).rev() {
        let mut s = rhs[r];
        for c in (r + 1)..n {
            let t = cmul(a[r * n + c], x[c]);
            s[0] -= t[0];
            s[1] -= t[1];
        }
        x[r] = cdiv(s, a[r * n + r]);
    }
    Some(x)
}

impl ReducedHelmholtz {
    /// Build the reduction on `modes` (and, when `static_vectors`, one `(K + M)^-1 b_port` residual-flexibility vector per
    /// port, M-orthonormalized against everything before it) for the given boundary-face sets (ports).
    ///
    /// # Errors
    /// A bad face, a basis above [`MAX_REDUCED_DIM`], or a sparse factorization refusal.
    pub fn new(
        fe: &HelmholtzFe<'_>,
        modes: &HelmholtzModes,
        ports: &[&[[u32; 3]]],
        static_vectors: bool,
    ) -> Result<Self, AcousticFeError> {
        Self::build(fe, modes, ports, static_vectors, None)
    }

    fn build(
        fe: &HelmholtzFe<'_>,
        modes: &HelmholtzModes,
        ports: &[&[[u32; 3]]],
        static_vectors: bool,
        wall: Option<(
            &[[u32; 3]],
            boundary_layer::BoundaryLayerTransport,
            boundary_layer::BoundaryLayerGeometry,
        )>,
    ) -> Result<Self, AcousticFeError> {
        Self::build_with_port_loads(fe, modes, ports, static_vectors, wall, None)
    }

    fn build_with_port_loads(
        fe: &HelmholtzFe<'_>,
        modes: &HelmholtzModes,
        ports: &[&[[u32; 3]]],
        static_vectors: bool,
        wall: Option<(
            &[[u32; 3]],
            boundary_layer::BoundaryLayerTransport,
            boundary_layer::BoundaryLayerGeometry,
        )>,
        supplied_loads: Option<&[Vec<f64>]>,
    ) -> Result<Self, AcousticFeError> {
        let mut loads = Vec::with_capacity(ports.len());
        let mut cs = Vec::with_capacity(ports.len());
        let mut area = Vec::with_capacity(ports.len());
        for faces in ports {
            let (b, a) = fe.boundary_load(faces)?;
            loads.push(b);
            area.push(a);
            cs.push(fe.boundary_mass(faces)?);
        }
        if let Some(supplied) = supplied_loads {
            if supplied.len() != ports.len()
                || supplied.iter().any(|b| {
                    b.len() != fe.ndof() || b.iter().any(|value| !value.is_finite())
                })
            {
                return Err(AcousticFeError::InvalidInput {
                    what: "spatial port load shape or nonfinite value",
                });
            }
            loads.clone_from_slice(supplied);
        }
        let mut basis: Vec<Vec<f64>> = modes.shapes.clone();
        if static_vectors {
            let ks = csr_sum(&fe.stiffness, &fe.mass);
            let sym = SymbolicLdlt::analyze(&ks, DirectOrdering::Amd)?;
            let factor = sym.factor(&ks, &LdltOptions::default())?;
            for b in &loads {
                let mut x = factor.solve(b);
                let norm0 = dot(&x, &spmv_new(&fe.mass, &x)).max(0.0).sqrt();
                for _ in 0..2 {
                    let mx = spmv_new(&fe.mass, &x);
                    for phi in &basis {
                        let c = dot(phi, &mx);
                        for (xi, pi) in x.iter_mut().zip(phi) {
                            *xi -= c * pi;
                        }
                    }
                }
                let norm = dot(&x, &spmv_new(&fe.mass, &x)).max(0.0).sqrt();
                if norm > 1e-10 * norm0 && norm0 > 0.0 {
                    for xi in &mut x {
                        *xi /= norm;
                    }
                    basis.push(x);
                }
            }
        }
        let n = basis.len();
        if n > MAX_REDUCED_DIM {
            return Err(AcousticFeError::ReducedBasisTooLarge {
                size: n,
                cap: MAX_REDUCED_DIM,
            });
        }
        let project = |a: &Csr| -> Vec<f64> {
            let mut out = vec![0.0; n * n];
            for j in 0..n {
                let y = spmv_new(a, &basis[j]);
                for i in 0..n {
                    out[i * n + j] = dot(&basis[i], &y);
                }
            }
            out
        };
        let k_r = project(&fe.stiffness);
        let m_r = project(&fe.mass);
        let c_r: Vec<Vec<f64>> = cs.iter().map(&project).collect();
        let b_r: Vec<Vec<f64>> = loads
            .iter()
            .map(|b| basis.iter().map(|phi| dot(phi, b)).collect::<Vec<f64>>())
            .collect();
        let boundary_layer = wall
            .map(|(faces, transport, geometry)| {
                boundary_layer::ReducedBoundaryLayer::new(
                    fe, faces, ports, transport, geometry, &project,
                )
            })
            .transpose()?;
        Ok(Self {
            medium: fe.medium,
            n,
            k_r,
            m_r,
            c_r,
            b_r,
            area,
            boundary_layer,
        })
    }

    /// Dimension of the reduced basis.
    #[must_use]
    pub const fn dim(&self) -> usize {
        self.n
    }

    /// Solve at `frequency_hz` with the port `source` driven at inward particle velocity `v_in` [m/s] and every other port
    /// loaded with the complex specific admittance `admittances[port] = rho c / Z` (the source entry is ignored: its face is
    /// otherwise rigid).
    ///
    /// # Errors
    /// Inconsistent inputs, or a singular reduced system at this frequency.
    pub fn response(
        &self,
        frequency_hz: f64,
        source: usize,
        v_in: f64,
        admittances: &[[f64; 2]],
    ) -> Result<ReducedResponse, AcousticFeError> {
        let ports = self.b_r.len();
        if source >= ports || admittances.len() != ports {
            return Err(AcousticFeError::InvalidInput {
                what: "source index or admittance count",
            });
        }
        if !(frequency_hz.is_finite() && frequency_hz > 0.0 && v_in.is_finite()) {
            return Err(AcousticFeError::InvalidInput {
                what: "frequency or inflow velocity",
            });
        }
        let n = self.n;
        let rho_c = self.medium.density_kg_m3 * self.medium.sound_speed_m_s;
        let k = 2.0 * core::f64::consts::PI * frequency_hz / self.medium.sound_speed_m_s;
        let mut a = vec![[0.0, 0.0]; n * n];
        for i in 0..n * n {
            a[i][0] = self.k_r[i] - k * k * self.m_r[i];
        }
        if let Some(wall) = &self.boundary_layer {
            wall.add_operator(self.medium, frequency_hz, &mut a)?;
        }
        for (p, y) in admittances.iter().enumerate() {
            if p == source {
                continue;
            }
            // j k y C: (re, im) of j k (yr + j yi) = (-k yi, k yr)
            let jy = [-k * y[1], k * y[0]];
            for i in 0..n * n {
                a[i][0] += jy[0] * self.c_r[p][i];
                a[i][1] += jy[1] * self.c_r[p][i];
            }
        }
        let scale = k * rho_c * v_in;
        let rhs: Vec<[f64; 2]> = self.b_r[source].iter().map(|b| [0.0, scale * b]).collect();
        let x = solve_complex(a, rhs, n)
            .ok_or(AcousticFeError::SingularAtFrequency { frequency_hz })?;
        let inflow = v_in * self.area[source];
        let mut outward_flow = Vec::with_capacity(ports);
        let mut mean_pressure = Vec::with_capacity(ports);
        for p in 0..ports {
            let mut s = [0.0, 0.0];
            for (b, xi) in self.b_r[p].iter().zip(&x) {
                s[0] += b * xi[0];
                s[1] += b * xi[1];
            }
            mean_pressure.push([s[0] / self.area[p], s[1] / self.area[p]]);
            if p == source {
                outward_flow.push([-inflow, 0.0]);
            } else {
                let y = admittances[p];
                let f = cmul([y[0] / rho_c, y[1] / rho_c], s);
                outward_flow.push(f);
            }
        }
        Ok(ReducedResponse {
            amplitudes: x,
            outward_flow,
            mean_pressure,
            inflow,
        })
    }
}

impl ReducedResponse {
    /// The volume-velocity transfer function from the driven port to port `port`: outward flow over inflow.
    #[must_use]
    pub fn transfer_to(&self, port: usize) -> [f64; 2] {
        [
            self.outward_flow[port][0] / self.inflow,
            self.outward_flow[port][1] / self.inflow,
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fixtures::kuhn_cube;

    const AIR: Medium = Medium {
        density_kg_m3: 1.2,
        sound_speed_m_s: 343.0,
    };

    fn scaled_box(n: usize, l: [f64; 3]) -> (TetComplex, Vec<[f64; 3]>) {
        let (c, p) = kuhn_cube(n);
        let p = p
            .iter()
            .map(|q| [q[0] * l[0], q[1] * l[1], q[2] * l[2]])
            .collect();
        (c, p)
    }

    fn plane_faces(
        fe: &HelmholtzFe<'_>,
        pos: &[[f64; 3]],
        axis: usize,
        value: f64,
    ) -> Vec<[u32; 3]> {
        fe.boundary_faces()
            .into_iter()
            .filter(|f| {
                f.iter()
                    .all(|&v| (pos[v as usize][axis] - value).abs() < 1e-12)
            })
            .collect()
    }

    /// Analytic rigid-box frequencies (l, m, n) != 0, ascending.
    fn box_frequencies(l: [f64; 3], count: usize) -> Vec<f64> {
        let mut f = Vec::new();
        for a in 0..8usize {
            for b in 0..8usize {
                for c in 0..8usize {
                    if a + b + c == 0 {
                        continue;
                    }
                    let q = (a as f64 / l[0]).powi(2)
                        + (b as f64 / l[1]).powi(2)
                        + (c as f64 / l[2]).powi(2);
                    f.push(0.5 * AIR.sound_speed_m_s * q.sqrt());
                }
            }
        }
        f.sort_by(|x, y| x.partial_cmp(y).expect("finite"));
        f.truncate(count);
        f
    }

    #[test]
    fn helmholtz_rigid_box_modes() {
        let dims = [1.0, 0.7, 0.4];
        let exact = box_frequencies(dims, 24);
        let window = (-1.0, 0.5 * (exact[19] + exact[20]));
        // Pre-registered gates, fixed before the order slopes were measured: the finest pair of rungs of every
        // order must show a slope of the lowest nonzero frequency error of at least 2r - 0.2 (asymptotic order
        // 2r), and the 20 modes at the finest rung of r = 4 must be within 1e-4. A shortfall fails the test
        // after every order has printed: report it as pre-asymptotic (add a finer rung) or as a defect, never
        // lower the target here.
        let mut shortfalls = Vec::new();
        for (order, rungs) in [
            (2usize, vec![3usize, 4, 6, 8]),
            (3, vec![2, 3, 4, 6]),
            (4, vec![2, 3, 4, 6]),
        ] {
            let mut first_mode_err = Vec::new();
            for &n in &rungs {
                let t0 = std::time::Instant::now();
                let (complex, pos) = scaled_box(n, dims);
                let fe = HelmholtzFe::new(&complex, &pos, order, AIR).expect("fe");
                let t_assemble = t0.elapsed().as_secs_f64();
                let modes = fe.modes(window, &SliceOptions::default()).expect("modes");
                let t_modes = t0.elapsed().as_secs_f64() - t_assemble;
                assert_eq!(modes.below_high - modes.below_low, modes.lambda.len());
                let found = modes.frequency_hz.len();
                println!(
                    "{{\"suite\":\"fs-feec/acoustics\",\"case\":\"rigid_box\",\"r\":{order},\"mesh_n\":{n},\"dofs\":{},\"modes\":{found},\"t_assemble_s\":{t_assemble:.3},\"t_modes_s\":{t_modes:.3}}}",
                    fe.ndof()
                );
                assert!(modes.frequency_hz[0] < 1e-3 * exact[0], "constant mode");
                if n == *rungs.last().expect("rungs") {
                    // certified count: the window holds exactly the analytic modes, multiplicities included
                    assert_eq!(
                        found, 21,
                        "r={order} n={n}: the window must hold the constant mode and 20 modes"
                    );
                }
                if found == 21 {
                    let worst = (1..21)
                        .map(|i| (modes.frequency_hz[i] / exact[i - 1] - 1.0).abs())
                        .fold(0.0, f64::max);
                    println!(
                        "{{\"suite\":\"fs-feec/acoustics\",\"case\":\"rigid_box\",\"r\":{order},\"mesh_n\":{n},\"worst_rel_err_20\":{worst:e}}}"
                    );
                    if order == 4 && n == *rungs.last().expect("rungs") && worst >= 1e-4 {
                        shortfalls.push(format!(
                            "r={order} n={n}: worst relative error {worst:e} >= 1e-4"
                        ));
                    }
                }
                first_mode_err.push((
                    1.0 / n as f64,
                    (modes.frequency_hz[1] / exact[0] - 1.0).abs(),
                ));
            }
            // observed order of the lowest nonzero frequency over every pair of rungs
            for pair in first_mode_err.windows(2) {
                let (h1, e1) = pair[0];
                let (h2, e2) = pair[1];
                if e1 > 1e-13 && e2 > 1e-13 {
                    let slope = (e1 / e2).ln() / (h1 / h2).ln();
                    println!(
                        "{{\"suite\":\"fs-feec/acoustics\",\"case\":\"rigid_box_order\",\"r\":{order},\"h_coarse\":{h1},\"h_fine\":{h2},\"observed_order\":{slope}}}"
                    );
                }
            }
            let target = 2.0 * order as f64 - 0.2;
            let (h1, e1) = first_mode_err[first_mode_err.len() - 2];
            let (h2, e2) = first_mode_err[first_mode_err.len() - 1];
            if e1 > 1e-13 && e2 > 1e-13 {
                let slope = (e1 / e2).ln() / (h1 / h2).ln();
                if slope < target {
                    shortfalls.push(format!("r={order}: observed order {slope} < {target}"));
                }
            } else {
                shortfalls.push(format!(
                    "r={order}: the finest pair is at the rounding floor ({e1:e}, {e2:e})"
                ));
            }
        }
        // every order has printed its slopes before any shortfall fails the test
        assert!(shortfalls.is_empty(), "{shortfalls:?}");
    }

    /// Exact transfer of a uniform duct driven by inflow velocity at x = 0 and closed by the specific admittance y at x = L.
    fn duct_exact(k: f64, len: f64, y: [f64; 2]) -> [f64; 2] {
        let gamma = cdiv([1.0 - y[0], -y[1]], [1.0 + y[0], y[1]]);
        let ejkl = [(k * len).cos(), (k * len).sin()];
        let emjkl = [(k * len).cos(), -(k * len).sin()];
        let g_e = cmul(gamma, emjkl);
        let den = [ejkl[0] - g_e[0], ejkl[1] - g_e[1]];
        cdiv([1.0 - gamma[0], -gamma[1]], den)
    }

    #[test]
    fn helmholtz_duct_driven_transfer_matches_plane_wave_theory() {
        let len = 0.5;
        let dims = [len, 0.05, 0.05];
        let y = [0.4, 0.3];
        let freqs = [150.0, 300.0, 450.0, 600.0];
        let mut last = None;
        let mut ladder = Vec::new();
        for &n in &[4usize, 6, 8] {
            let (complex, pos) = scaled_box(n, dims);
            let fe = HelmholtzFe::new(&complex, &pos, 2, AIR).expect("fe");
            let src = plane_faces(&fe, &pos, 0, 0.0);
            let end = plane_faces(&fe, &pos, 0, len);
            let t0 = std::time::Instant::now();
            // below the first cross mode c/(2*0.05) = 3430 Hz: the response is purely longitudinal
            let modes = fe
                .modes((-1.0, 3000.0), &SliceOptions::default())
                .expect("modes");
            let t_modes = t0.elapsed().as_secs_f64();
            for &with_static in &[false, true] {
                let rom =
                    ReducedHelmholtz::new(&fe, &modes, &[&src, &end], with_static).expect("rom");
                let mut worst = 0.0f64;
                for &f in &freqs {
                    let r = rom.response(f, 0, 1.0, &[[0.0, 0.0], y]).expect("response");
                    let h = r.transfer_to(1);
                    let k = 2.0 * core::f64::consts::PI * f / AIR.sound_speed_m_s;
                    let e = duct_exact(k, len, y);
                    let err = ((h[0] - e[0]).hypot(h[1] - e[1])) / e[0].hypot(e[1]);
                    worst = worst.max(err);
                }
                println!(
                    "{{\"suite\":\"fs-feec/acoustics\",\"case\":\"duct_driven\",\"mesh_n\":{n},\"static_vectors\":{with_static},\"modes\":{},\"basis\":{},\"worst_rel_err\":{worst:e},\"t_modes_s\":{t_modes:.2}}}",
                    modes.lambda.len(),
                    rom.dim()
                );
                if with_static {
                    last = Some(worst);
                    ladder.push(worst);
                }
            }
        }
        assert!(
            ladder.windows(2).all(|w| w[1] < w[0]),
            "the error with static vectors must fall with the mesh: {ladder:?}"
        );
        let worst = last.expect("ran");
        assert!(
            worst < 2e-3,
            "finest duct rung with static vectors: {worst:e}"
        );
    }

    #[test]
    fn helmholtz_refusals() {
        let (complex, pos) = scaled_box(2, [1.0, 1.0, 1.0]);
        assert!(HelmholtzFe::new(&complex, &pos, MAX_ORDER, AIR).is_ok());
        assert_eq!(
            HelmholtzFe::new(&complex, &pos, MAX_ORDER + 1, AIR).err(),
            Some(AcousticFeError::OrderOutOfRange {
                order: MAX_ORDER + 1,
                cap: MAX_ORDER
            })
        );
        assert!(matches!(
            HelmholtzFe::new(&complex, &pos, 0, AIR).err(),
            Some(AcousticFeError::OrderOutOfRange { .. })
        ));
        let bad = Medium {
            density_kg_m3: 0.0,
            sound_speed_m_s: 343.0,
        };
        assert!(matches!(
            HelmholtzFe::new(&complex, &pos, 2, bad).err(),
            Some(AcousticFeError::InvalidMedium { what: "density" })
        ));
        let mut flat = pos.clone();
        for p in &mut flat {
            p[2] = 0.0;
        }
        assert!(matches!(
            HelmholtzFe::new(&complex, &flat, 2, AIR).err(),
            Some(AcousticFeError::DegenerateTet { .. })
        ));
        assert!(matches!(
            HelmholtzFe::new(&complex, &pos[1..], 2, AIR).err(),
            Some(AcousticFeError::InvalidPositions { .. })
        ));
        let fe = HelmholtzFe::new(&complex, &pos, 2, AIR).expect("fe");
        // an interior face and a nonexistent face refuse as boundaries
        let boundary = fe.boundary_faces();
        let interior = complex
            .faces
            .iter()
            .find(|f| !boundary.contains(f))
            .copied()
            .expect("interior face");
        assert!(matches!(
            fe.boundary_load(&[interior]).err(),
            Some(AcousticFeError::FaceNotOnBoundary { .. })
        ));
        assert!(matches!(
            fe.boundary_load(&[[0, 1, 2_000_000]]).err(),
            Some(AcousticFeError::FaceNotInMesh { .. })
        ));
        assert!(matches!(
            fe.modes((5.0, 5.0), &SliceOptions::default()).err(),
            Some(AcousticFeError::Modal(ModalError::InvalidWindow { .. }))
        ));
        let modes = fe
            .modes((-1.0, 400.0), &SliceOptions::default())
            .expect("modes");
        let rom = ReducedHelmholtz::new(&fe, &modes, &[&boundary], false).expect("rom");
        assert!(matches!(
            rom.response(100.0, 3, 1.0, &[[0.0, 0.0]]).err(),
            Some(AcousticFeError::InvalidInput { .. })
        ));
        assert!(matches!(
            rom.response(-1.0, 0, 1.0, &[[0.0, 0.0]]).err(),
            Some(AcousticFeError::InvalidInput { .. })
        ));
    }

    #[test]
    fn linear_boundary_trace_integrates_the_product_and_is_vertex_order_invariant() {
        let (complex, positions) = scaled_box(1, [1.0, 1.0, 1.0]);
        for order in 1..=MAX_ORDER {
            let fe = HelmholtzFe::new(&complex, &positions, order, AIR).unwrap();
            let face = fe.boundary_faces()[0];
            let (_, _, area) = fe.boundary_geometry(face).unwrap();
            let values = [1.0, -2.0, 4.0];
            let load = fe.boundary_linear_load(&[(face, values)]).unwrap();
            // The three pressure vertex functions are barycentric coordinates:
            // int N_i N_j = A/6 on the diagonal and A/12 off diagonal.
            for i in 0..3 {
                let expected = area * (values.iter().sum::<f64>() + values[i]) / 12.0;
                assert!((load[face[i] as usize] - expected).abs() < 2e-13 * area);
            }
            let permuted = fe.boundary_linear_load(&[
                ([face[2], face[0], face[1]], [values[2], values[0], values[1]])
            ]).unwrap();
            assert_eq!(load, permuted);
            let uniform = fe.boundary_linear_load(&[(face, [1.0; 3])]).unwrap();
            let (original, _) = fe.boundary_load(&[face]).unwrap();
            for (a,b) in uniform.iter().zip(original) { assert!((a-b).abs() < 2e-13 * area); }
            assert!(fe.boundary_linear_load(&[(face, values), (face, values)]).is_err());
            assert!(fe.boundary_linear_load(&[(face, [f64::NAN, 0.0, 0.0])]).is_err());
        }
    }

    #[test]
    fn boundary_geometry_uses_the_owning_volume_for_its_outward_normal() {
        let (complex, positions) = scaled_box(1, [0.4, 0.3, 0.2]);
        let fe = HelmholtzFe::new(&complex, &positions, 1, AIR).unwrap();
        for face in fe.boundary_faces() {
            let (vertices, normal, area) = fe.boundary_geometry(face).unwrap();
            let center = [0.2, 0.15, 0.1];
            assert!(dot3(normal, sub3(vertices[0], center)) > 0.0);
            let (_, reversed, a) = fe.boundary_geometry([face[2], face[1], face[0]]).unwrap();
            assert!((area-a).abs() < 1e-14);
            for (n,r) in normal.iter().zip(reversed) { assert!((n-r).abs() < 1e-14); }
        }
    }

    #[test]
    fn helmholtz_modes_repeat_bitwise() {
        let (complex, pos) = scaled_box(3, [1.0, 0.7, 0.4]);
        let fe = HelmholtzFe::new(&complex, &pos, 2, AIR).expect("fe");
        let a = fe
            .modes((-1.0, 400.0), &SliceOptions::default())
            .expect("modes");
        let b = fe
            .modes((-1.0, 400.0), &SliceOptions::default())
            .expect("modes");
        assert_eq!(a.lambda.len(), b.lambda.len());
        for (x, y) in a.lambda.iter().zip(&b.lambda) {
            assert_eq!(x.to_bits(), y.to_bits());
        }
    }
}



