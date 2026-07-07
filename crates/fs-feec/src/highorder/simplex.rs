//! Hierarchical high-order H¹ on tet complexes (tfz.6 slice 3):
//! Szabó–Babuška-style entity hierarchy — vertex λ_a, edge
//! λ_aλ_b·P_k(λ_b−λ_a), face λ_aλ_bλ_c·P_i(λ_b−λ_a)·P_j(λ_c−λ_b),
//! interior λ_0λ_1λ_2λ_3·(products of Legendre kernels) — with the
//! DETERMINISTIC GLOBAL ORIENTATION convention: every entity's kernel
//! arguments use the entity's vertices in SORTED GLOBAL INDEX order
//! (fs-rep-mesh's canonical tables), so two elements sharing an edge
//! or face evaluate IDENTICAL trace polynomials without any sign or
//! permutation bookkeeping — the classic high-order FEEC bug farm
//! removed by construction, and the permutation battery (G3) checks
//! it rather than trusting it.
//!
//! Quadrature: collapsed Duffy tensor Gauss–Legendre (exactness
//! chosen per order). Assembly is straightforward element loops into
//! deterministic COO — matrix-free/perf paths are slice 4's concern.

use crate::highorder::quad1d::{gauss_legendre, legendre};
use fs_rep_mesh::TetComplex;
use fs_sparse::{Coo, Csr};

/// Per-entity dof counts at order r.
#[must_use]
pub fn entity_dofs(r: usize) -> (usize, usize, usize) {
    let e = r.saturating_sub(1);
    let f = if r >= 3 { (r - 1) * (r - 2) / 2 } else { 0 };
    let i = if r >= 4 {
        (r - 1) * (r - 2) * (r - 3) / 6
    } else {
        0
    };
    (e, f, i)
}

/// The global hierarchical H¹ space on a tet complex.
pub struct SimplexSpace<'c> {
    /// The complex (canonical sorted entity tables).
    pub complex: &'c TetComplex,
    /// Order r ≥ 1.
    pub r: usize,
    /// Dofs per edge / face / cell.
    pub per_edge: usize,
    /// Dofs per face.
    pub per_face: usize,
    /// Dofs per cell interior.
    pub per_cell: usize,
    /// Total dof count.
    pub ndof: usize,
    /// Offsets: edges start after vertices, then faces, then cells.
    pub edge_off: usize,
    /// Face block offset.
    pub face_off: usize,
    /// Cell block offset.
    pub cell_off: usize,
}

impl<'c> SimplexSpace<'c> {
    /// Build the space.
    ///
    /// # Panics
    /// If `r == 0`.
    #[must_use]
    pub fn new(complex: &'c TetComplex, r: usize) -> SimplexSpace<'c> {
        assert!(r >= 1, "H1 needs r >= 1");
        let (pe, pf, pi) = entity_dofs(r);
        let edge_off = complex.vertex_count;
        let face_off = edge_off + complex.edges.len() * pe;
        let cell_off = face_off + complex.faces.len() * pf;
        let ndof = cell_off + complex.tets.len() * pi;
        SimplexSpace {
            complex,
            r,
            per_edge: pe,
            per_face: pf,
            per_cell: pi,
            ndof,
            edge_off,
            face_off,
            cell_off,
        }
    }

    /// Local basis enumeration for one tet: returns (global dof, local
    /// descriptor) pairs. Local descriptors reference the tet's
    /// vertices by LOCAL index but with entity vertices already in
    /// sorted-GLOBAL order — the orientation convention.
    #[must_use]
    pub fn element_dofs(&self, t: usize) -> Vec<(usize, LocalFn)> {
        let tet = self.complex.tets[t];
        let mut out = Vec::new();
        // Vertices.
        for (l, &v) in tet.iter().enumerate() {
            out.push((v as usize, LocalFn::Vertex(l)));
        }
        // Edges: all 6 local pairs, sorted by global index.
        for p in 0..4 {
            for q in (p + 1)..4 {
                let (la, lb) = if tet[p] < tet[q] { (p, q) } else { (q, p) };
                let key = [tet[la], tet[lb]];
                let e = self
                    .complex
                    .edges
                    .binary_search(&key)
                    .expect("edge in table");
                for k in 0..self.per_edge {
                    out.push((
                        self.edge_off + e * self.per_edge + k,
                        LocalFn::Edge(la, lb, k),
                    ));
                }
            }
        }
        // Faces: 4 local triples, sorted by global index.
        if self.per_face > 0 {
            for omit in 0..4 {
                let mut locs: Vec<usize> = (0..4).filter(|&l| l != omit).collect();
                locs.sort_by_key(|&l| tet[l]);
                let key = {
                    let mut kk = [tet[locs[0]], tet[locs[1]], tet[locs[2]]];
                    kk.sort_unstable();
                    kk
                };
                let f = self
                    .complex
                    .faces
                    .binary_search(&key)
                    .expect("face in table");
                let mut k = 0usize;
                for i in 0..self.r.saturating_sub(2) {
                    for j in 0..self.r.saturating_sub(2) {
                        if i + j <= self.r.saturating_sub(3) {
                            out.push((
                                self.face_off + f * self.per_face + k,
                                LocalFn::Face(locs[0], locs[1], locs[2], i, j),
                            ));
                            k += 1;
                        }
                    }
                }
            }
        }
        // Interior.
        if self.per_cell > 0 {
            let mut k = 0usize;
            for i in 0..self.r.saturating_sub(3) {
                for j in 0..self.r.saturating_sub(3) {
                    for l in 0..self.r.saturating_sub(3) {
                        if i + j + l <= self.r.saturating_sub(4) {
                            out.push((
                                self.cell_off + t * self.per_cell + k,
                                LocalFn::Interior(i, j, l),
                            ));
                            k += 1;
                        }
                    }
                }
            }
        }
        out
    }

    /// Assemble the global stiffness ∫ ∇φ_i·∇φ_j (deterministic COO).
    #[must_use]
    pub fn stiffness(&self, positions: &[[f64; 3]]) -> Csr {
        let geo = crate::whitney::element_geometry(self.complex, positions);
        let quad = duffy_quadrature(self.r + 2);
        let mut coo = Coo::new(self.ndof, self.ndof);
        for t in 0..self.complex.tets.len() {
            let dofs = self.element_dofs(t);
            let nl = dofs.len();
            let vol_jac = 6.0 * geo.vol_signed[t].abs(); // reference volume 1/6
            // Gradients at each quadrature point: ∇φ = Σ_a (∂φ/∂λ_a)·∇λ_a.
            let mut ke = vec![0.0f64; nl * nl];
            for &(lam, w) in &quad {
                let mut grads = Vec::with_capacity(nl);
                for (_, f) in &dofs {
                    let dl = f.d_lambda(lam, self.r);
                    let mut g = [0.0f64; 3];
                    for (a, dla) in dl.iter().enumerate() {
                        for (c, gc) in g.iter_mut().enumerate() {
                            *gc = dla.mul_add(geo.grads[t][a][c], *gc);
                        }
                    }
                    grads.push(g);
                }
                for i in 0..nl {
                    for j in 0..nl {
                        let dot = grads[i][0].mul_add(
                            grads[j][0],
                            grads[i][1].mul_add(grads[j][1], grads[i][2] * grads[j][2]),
                        );
                        ke[i * nl + j] = (w * vol_jac).mul_add(dot, ke[i * nl + j]);
                    }
                }
            }
            for (i, &(gi, _)) in dofs.iter().enumerate() {
                for (j, &(gj, _)) in dofs.iter().enumerate() {
                    coo.push(gi, gj, ke[i * nl + j]);
                }
            }
        }
        coo.assemble()
    }

    /// Assemble the global mass ∫ φ_i·φ_j (deterministic COO).
    #[must_use]
    pub fn mass(&self, positions: &[[f64; 3]]) -> Csr {
        let geo = crate::whitney::element_geometry(self.complex, positions);
        let quad = duffy_quadrature(self.r + 2);
        let mut coo = Coo::new(self.ndof, self.ndof);
        for t in 0..self.complex.tets.len() {
            let dofs = self.element_dofs(t);
            let nl = dofs.len();
            let vol_jac = 6.0 * geo.vol_signed[t].abs();
            let mut me = vec![0.0f64; nl * nl];
            for &(lam, w) in &quad {
                let vals: Vec<f64> = dofs.iter().map(|(_, f)| f.eval(lam, self.r)).collect();
                for i in 0..nl {
                    for j in 0..nl {
                        me[i * nl + j] = (w * vol_jac * vals[i]).mul_add(vals[j], me[i * nl + j]);
                    }
                }
            }
            for (i, &(gi, _)) in dofs.iter().enumerate() {
                for (j, &(gj, _)) in dofs.iter().enumerate() {
                    coo.push(gi, gj, me[i * nl + j]);
                }
            }
        }
        coo.assemble()
    }

    /// Load vector ∫ f·φ_i by Duffy quadrature (two extra points for
    /// smooth forcing).
    #[must_use]
    pub fn load<F: Fn([f64; 3]) -> f64>(&self, positions: &[[f64; 3]], f: &F) -> Vec<f64> {
        let geo = crate::whitney::element_geometry(self.complex, positions);
        let quad = duffy_quadrature(self.r + 4);
        let mut b = vec![0.0f64; self.ndof];
        for (t, tet) in self.complex.tets.iter().enumerate() {
            let dofs = self.element_dofs(t);
            let vol_jac = 6.0 * geo.vol_signed[t].abs();
            let corners: Vec<[f64; 3]> = tet.iter().map(|&v| positions[v as usize]).collect();
            for &(lam, w) in &quad {
                let mut p = [0.0f64; 3];
                for a in 0..4 {
                    for c in 0..3 {
                        p[c] = lam[a].mul_add(corners[a][c], p[c]);
                    }
                }
                let fw = w * vol_jac * f(p);
                for (gi, lf) in &dofs {
                    b[*gi] = fw.mul_add(lf.eval(lam, self.r), b[*gi]);
                }
            }
        }
        b
    }

    /// L2 error of a dof vector against an analytic field.
    #[must_use]
    pub fn l2_error<F: Fn([f64; 3]) -> f64>(
        &self,
        positions: &[[f64; 3]],
        u: &[f64],
        u_exact: &F,
    ) -> f64 {
        let geo = crate::whitney::element_geometry(self.complex, positions);
        let quad = duffy_quadrature(self.r + 4);
        let mut total = 0.0f64;
        for (t, tet) in self.complex.tets.iter().enumerate() {
            let dofs = self.element_dofs(t);
            let vol_jac = 6.0 * geo.vol_signed[t].abs();
            let corners: Vec<[f64; 3]> = tet.iter().map(|&v| positions[v as usize]).collect();
            for &(lam, w) in &quad {
                let mut p = [0.0f64; 3];
                for a in 0..4 {
                    for c in 0..3 {
                        p[c] = lam[a].mul_add(corners[a][c], p[c]);
                    }
                }
                let mut uh = 0.0f64;
                for (gi, lf) in &dofs {
                    uh = u[*gi].mul_add(lf.eval(lam, self.r), uh);
                }
                let e = uh - u_exact(p);
                total += w * vol_jac * e * e;
            }
        }
        fs_math::det::sqrt(total)
    }

    /// Boundary dof mask for the unit-cube fixtures: a dof is boundary
    /// iff its entity lies in a boundary FACE (face incident to one
    /// tet). Interior bubbles are never boundary.
    #[must_use]
    pub fn boundary_mask(&self) -> Vec<bool> {
        let d2 = self.complex.d2();
        let mut face_use = vec![0usize; self.complex.faces.len()];
        for row in &d2.rows {
            for &(f, _) in row {
                face_use[f] += 1;
            }
        }
        let mut mask = vec![false; self.ndof];
        for (f, &uses) in face_use.iter().enumerate() {
            if uses == 1 {
                let tri = self.complex.faces[f];
                for &v in &tri {
                    mask[v as usize] = true;
                }
                for p in 0..3 {
                    for q in (p + 1)..3 {
                        let key = if tri[p] < tri[q] {
                            [tri[p], tri[q]]
                        } else {
                            [tri[q], tri[p]]
                        };
                        let e = self.complex.edges.binary_search(&key).expect("edge");
                        for k in 0..self.per_edge {
                            mask[self.edge_off + e * self.per_edge + k] = true;
                        }
                    }
                }
                for k in 0..self.per_face {
                    mask[self.face_off + f * self.per_face + k] = true;
                }
            }
        }
        mask
    }
}

/// A local basis function in barycentric form.
#[derive(Debug, Clone, Copy)]
pub enum LocalFn {
    /// λ_l.
    Vertex(usize),
    /// λ_a λ_b P_k(λ_b − λ_a) with (a, b) sorted globally.
    Edge(usize, usize, usize),
    /// λ_a λ_b λ_c P_i(λ_b − λ_a) P_j(λ_c − λ_b), sorted globally.
    Face(usize, usize, usize, usize, usize),
    /// λ_0λ_1λ_2λ_3 P_i(λ_1−λ_0) P_j(λ_2−λ_1) P_l(λ_3−λ_2).
    Interior(usize, usize, usize),
}

impl LocalFn {
    /// Value at barycentric point λ.
    #[must_use]
    pub fn eval(&self, lam: [f64; 4], _r: usize) -> f64 {
        match *self {
            LocalFn::Vertex(l) => lam[l],
            LocalFn::Edge(a, b, k) => {
                let (pk, _) = legendre(k, lam[b] - lam[a]);
                lam[a] * lam[b] * pk
            }
            LocalFn::Face(a, b, c, i, j) => {
                let (pi, _) = legendre(i, lam[b] - lam[a]);
                let (pj, _) = legendre(j, lam[c] - lam[b]);
                lam[a] * lam[b] * lam[c] * pi * pj
            }
            LocalFn::Interior(i, j, l) => {
                let (pi, _) = legendre(i, lam[1] - lam[0]);
                let (pj, _) = legendre(j, lam[2] - lam[1]);
                let (pl, _) = legendre(l, lam[3] - lam[2]);
                lam[0] * lam[1] * lam[2] * lam[3] * pi * pj * pl
            }
        }
    }

    /// Partial derivatives ∂φ/∂λ_a, a = 0..4 (barycentrics treated as
    /// independent — the assembly contracts with ∇λ_a, which encodes
    /// the constraint).
    #[must_use]
    pub fn d_lambda(&self, lam: [f64; 4], _r: usize) -> [f64; 4] {
        let mut d = [0.0f64; 4];
        match *self {
            LocalFn::Vertex(l) => d[l] = 1.0,
            LocalFn::Edge(a, b, k) => {
                let t = lam[b] - lam[a];
                let (pk, dpk) = legendre(k, t);
                d[a] = lam[b].mul_add(pk, -(lam[a] * lam[b] * dpk));
                d[b] = lam[a].mul_add(pk, lam[a] * lam[b] * dpk);
            }
            LocalFn::Face(a, b, c, i, j) => {
                let t1 = lam[b] - lam[a];
                let t2 = lam[c] - lam[b];
                let (pi, dpi) = legendre(i, t1);
                let (pj, dpj) = legendre(j, t2);
                let base = lam[a] * lam[b] * lam[c];
                d[a] = (lam[b] * lam[c]).mul_add(pi * pj, -(base * dpi * pj));
                d[b] = (lam[a] * lam[c]).mul_add(pi * pj, base * (dpi * pj - pi * dpj));
                d[c] = (lam[a] * lam[b]).mul_add(pi * pj, base * pi * dpj);
            }
            LocalFn::Interior(i, j, l) => {
                let t1 = lam[1] - lam[0];
                let t2 = lam[2] - lam[1];
                let t3 = lam[3] - lam[2];
                let (pi, dpi) = legendre(i, t1);
                let (pj, dpj) = legendre(j, t2);
                let (pl, dpl) = legendre(l, t3);
                let base = lam[0] * lam[1] * lam[2] * lam[3];
                let bub = |skip: usize| -> f64 {
                    (0..4).filter(|&a| a != skip).map(|a| lam[a]).product()
                };
                d[0] = bub(0).mul_add(pi * pj * pl, -(base * dpi * pj * pl));
                d[1] = bub(1).mul_add(pi * pj * pl, base * (dpi * pj * pl - pi * dpj * pl));
                d[2] = bub(2).mul_add(pi * pj * pl, base * (pi * dpj * pl - pi * pj * dpl));
                d[3] = bub(3).mul_add(pi * pj * pl, base * pi * pj * dpl);
            }
        }
        d
    }
}

/// Collapsed Duffy tensor quadrature on the reference tet:
/// barycentric points + weights summing to the REFERENCE volume 1/6,
/// applied against 6·|V| in assembly (so a constant integrates to
/// |V| exactly). `n` Gauss–Legendre points per direction.
#[must_use]
pub fn duffy_quadrature(n: usize) -> Vec<([f64; 4], f64)> {
    let (qx, qw) = gauss_legendre(n);
    // Map [-1,1] → [0,1].
    let map = |x: f64| f64::midpoint(1.0, x);
    let mut out = Vec::with_capacity(n * n * n);
    for (&xu, &wu) in qx.iter().zip(&qw) {
        let u = map(xu);
        for (&xv, &wv) in qx.iter().zip(&qw) {
            let v = map(xv);
            for (&xw, &ww) in qx.iter().zip(&qw) {
                let w = map(xw);
                let l1 = u;
                let l2 = v * (1.0 - u);
                let l3 = w * (1.0 - u) * (1.0 - v);
                let l0 = 1.0 - l1 - l2 - l3;
                // Jacobian of the collapse: (1−u)²(1−v); the [-1,1]→[0,1]
                // maps contribute (1/2)³. Reference-tet volume is 1/6;
                // normalize so Σw = 1/6 · 6 = 1 against 6|V|.
                let jac = (1.0 - u) * (1.0 - u) * (1.0 - v) / 8.0;
                out.push(([l0, l1, l2, l3], wu * wv * ww * jac));
            }
        }
    }
    out
}
