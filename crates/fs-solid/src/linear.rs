//! Small-strain elasticity on body-fitted meshes: plane strain or
//! plane stress, standard displacement elements, and the B-BAR
//! dilatation projection — the measured locking-free path. B-bar
//! replaces the dilatational part of the strain-displacement operator
//! with its element average (Hughes' projection; equivalent to a
//! condensed element-constant-pressure mixed method), which is what
//! removes the spurious volumetric constraints that lock standard
//! elements as ν → 1/2 and in thin bending. TDNNS-proper is a
//! recorded no-claim awaiting the simplicial H(div) families.

use crate::SolidError;
use crate::mesh2::{Mesh2, Patch, quad_points, shapes_at};
use fs_solver::krylov::CgState;
use fs_solver::op::CsrOp;
use fs_sparse::precond::Precond;
use fs_sparse::{Coo, Csr};
use std::collections::BTreeMap;

/// Plane reduction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlaneKind {
    /// Plane strain (ε₃₃ = 0).
    Strain,
    /// Plane stress (σ₃₃ = 0).
    Stress,
}

/// Displacement formulation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Formulation {
    /// Standard displacement elements (locks as ν → 1/2).
    Standard,
    /// B-bar dilatation projection (locking-free; the battery
    /// measures both).
    BBar,
}

/// A displacement/traction field on a patch.
pub type VecField<'m> = &'m dyn Fn(f64, f64) -> [f64; 2];

/// A linear elasticity problem on a body-fitted mesh.
pub struct LinearProblem<'m> {
    /// The mesh.
    pub mesh: &'m Mesh2,
    /// Young's modulus.
    pub youngs: f64,
    /// Poisson ratio.
    pub poisson: f64,
    /// Plane reduction.
    pub plane: PlaneKind,
    /// Formulation.
    pub formulation: Formulation,
    /// Body force.
    pub body_force: Option<&'m dyn Fn(f64, f64) -> [f64; 2]>,
    /// Dirichlet patches: (patch, prescribed displacement field).
    pub dirichlet: Vec<(Patch, VecField<'m>)>,
    /// Traction patches: (patch, traction field).
    pub traction: Vec<(Patch, VecField<'m>)>,
}

/// The Lamé pair for the plane reduction (plane stress uses the
/// standard effective λ).
#[must_use]
pub fn lame(youngs: f64, poisson: f64, plane: PlaneKind) -> (f64, f64) {
    let mu = 0.5 * youngs / (1.0 + poisson);
    let lambda = youngs * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    match plane {
        PlaneKind::Strain => (lambda, mu),
        PlaneKind::Stress => (2.0 * lambda * mu / (lambda + 2.0 * mu), mu),
    }
}

/// SPD diagonal preconditioner (shared by the frontends).
pub(crate) struct Jacobi {
    inv: Vec<f64>,
}

impl Jacobi {
    pub(crate) fn new(a: &Csr) -> Jacobi {
        Jacobi {
            inv: (0..a.nrows())
                .map(|i| {
                    let d = a.get(i, i);
                    if d > 0.0 { 1.0 / d } else { 1.0 }
                })
                .collect(),
        }
    }
}

impl Precond for Jacobi {
    fn apply(&self, r: &[f64], z: &mut [f64]) {
        for (zi, (ri, di)) in z.iter_mut().zip(r.iter().zip(&self.inv)) {
            *zi = ri * di;
        }
    }
}

impl LinearProblem<'_> {
    /// Assemble and solve; returns nodal displacements (2 per node).
    ///
    /// # Errors
    /// [`SolidError::SolveFailed`], [`SolidError::UnknownPatch`].
    pub fn solve(&self) -> Result<Vec<[f64; 2]>, SolidError> {
        let (lambda, mu) = lame(self.youngs, self.poisson, self.plane);
        let n = self.mesh.node_count();
        let ndof = 2 * n;
        // Strong Dirichlet values by node/component.
        let mut fixed: BTreeMap<usize, f64> = BTreeMap::new();
        for (patch, g) in &self.dirichlet {
            if self.mesh.patch_edges(*patch).is_none() {
                return Err(SolidError::UnknownPatch { patch: *patch });
            }
            for node in self.mesh.patch_nodes(*patch) {
                let p = self.mesh.nodes[node];
                let val = g(p[0], p[1]);
                fixed.insert(2 * node, val[0]);
                fixed.insert(2 * node + 1, val[1]);
            }
        }
        let mut coo = Coo::new(ndof, ndof);
        let mut rhs = vec![0.0f64; ndof];
        for conn in &self.mesh.elems {
            let (k, fl) = self.element(conn, lambda, mu);
            scatter_constrained(&mut coo, &mut rhs, conn, &k, &fl, &fixed);
        }
        // Traction terms (dead loads along patch edges, 2-pt Gauss).
        for (patch, t) in &self.traction {
            let edges = self
                .mesh
                .patch_edges(*patch)
                .ok_or(SolidError::UnknownPatch { patch: *patch })?;
            for &(a, b) in edges {
                let (pa, pb) = (self.mesh.nodes[a], self.mesh.nodes[b]);
                let len = ((pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2)).sqrt();
                let g = 0.5 / 3.0f64.sqrt();
                for s in [0.5 - g, 0.5 + g] {
                    let p = [pa[0] + s * (pb[0] - pa[0]), pa[1] + s * (pb[1] - pa[1])];
                    let tv = t(p[0], p[1]);
                    let w = 0.5 * len;
                    for (node, shape) in [(a, 1.0 - s), (b, s)] {
                        for (c, tc) in tv.iter().enumerate() {
                            let dof = 2 * node + c;
                            if !fixed.contains_key(&dof) {
                                rhs[dof] += w * shape * tc;
                            }
                        }
                    }
                }
            }
        }
        // Identity rows for fixed DOFs.
        for (&dof, &val) in &fixed {
            coo.push(dof, dof, 1.0);
            rhs[dof] = val;
        }
        let a = coo.assemble();
        let m = Jacobi::new(&a);
        let op = CsrOp::symmetric(a);
        let mut st = CgState::new(&op, &m, &rhs);
        let _ = st.run(&op, &m, 1e-12, 40_000);
        let rr = st.rel_residual();
        if !rr.is_finite() || rr > 1e-8 {
            return Err(SolidError::SolveFailed {
                iters: st.iters,
                rel_residual: rr,
            });
        }
        Ok((0..n).map(|i| [st.x[2 * i], st.x[2 * i + 1]]).collect())
    }

    /// One element's stiffness and body-force load (B-matrix form,
    /// Voigt [εxx, εyy, γxy]).
    fn element(&self, conn: &[usize], lambda: f64, mu: f64) -> (Vec<Vec<f64>>, Vec<f64>) {
        let nn = conn.len();
        let d = [
            [lambda + 2.0 * mu, lambda, 0.0],
            [lambda, lambda + 2.0 * mu, 0.0],
            [0.0, 0.0, mu],
        ];
        let pts = quad_points(nn);
        // B-bar: element-average dilatational gradient.
        let mut bbar = vec![[0.0f64; 2]; nn];
        let mut vol = 0.0f64;
        if self.formulation == Formulation::BBar {
            for &(xi, eta, w) in &pts {
                let (_, grads, det) = shapes_at(&self.mesh.nodes, conn, xi, eta);
                for (a, g) in grads.iter().enumerate() {
                    bbar[a][0] += w * det * g[0];
                    bbar[a][1] += w * det * g[1];
                }
                vol += w * det;
            }
            for g in &mut bbar {
                g[0] /= vol;
                g[1] /= vol;
            }
        }
        let mut k = vec![vec![0.0f64; 2 * nn]; 2 * nn];
        let mut fl = vec![0.0f64; 2 * nn];
        for &(xi, eta, w) in &pts {
            let (shape, grads, det) = shapes_at(&self.mesh.nodes, conn, xi, eta);
            let wq = w * det;
            // Rows of B per node/component, with the B-bar dilatation
            // swap when enabled.
            let b_of = |a: usize, c: usize| -> [f64; 3] {
                let g = grads[a];
                let mut b = if c == 0 {
                    [g[0], 0.0, g[1]]
                } else {
                    [0.0, g[1], g[0]]
                };
                if self.formulation == Formulation::BBar {
                    let corr = 0.5 * (bbar[a][c] - g[c]);
                    b[0] += corr;
                    b[1] += corr;
                }
                b
            };
            for a in 0..nn {
                for ca in 0..2 {
                    let ba = b_of(a, ca);
                    let db = [
                        d[0][0] * ba[0] + d[0][1] * ba[1],
                        d[1][0] * ba[0] + d[1][1] * ba[1],
                        d[2][2] * ba[2],
                    ];
                    for b in 0..nn {
                        for cb in 0..2 {
                            let bb = b_of(b, cb);
                            k[2 * a + ca][2 * b + cb] +=
                                wq * (db[0] * bb[0] + db[1] * bb[1] + db[2] * bb[2]);
                        }
                    }
                }
            }
            if let Some(bf) = self.body_force {
                let mut p = [0.0f64; 2];
                for (a, &node) in conn.iter().enumerate() {
                    p[0] += shape[a] * self.mesh.nodes[node][0];
                    p[1] += shape[a] * self.mesh.nodes[node][1];
                }
                let f = bf(p[0], p[1]);
                for a in 0..nn {
                    fl[2 * a] += wq * shape[a] * f[0];
                    fl[2 * a + 1] += wq * shape[a] * f[1];
                }
            }
        }
        (k, fl)
    }
}

/// Scatter an element matrix with strong-Dirichlet elimination
/// (symmetric: fixed columns move to the RHS, fixed rows drop).
pub(crate) fn scatter_constrained(
    coo: &mut Coo,
    rhs: &mut [f64],
    conn: &[usize],
    k: &[Vec<f64>],
    fl: &[f64],
    fixed: &BTreeMap<usize, f64>,
) {
    let nn = conn.len();
    let dof = |a: usize, c: usize| 2 * conn[a] + c;
    for a in 0..nn {
        for ca in 0..2 {
            let ia = dof(a, ca);
            if fixed.contains_key(&ia) {
                continue;
            }
            rhs[ia] += fl[2 * a + ca];
            for b in 0..nn {
                for cb in 0..2 {
                    let ib = dof(b, cb);
                    let v = k[2 * a + ca][2 * b + cb];
                    if v == 0.0 {
                        continue;
                    }
                    if let Some(&gb) = fixed.get(&ib) {
                        rhs[ia] -= v * gb;
                    } else {
                        coo.push(ia, ib, v);
                    }
                }
            }
        }
    }
}

/// L2 and H1-seminorm displacement errors against an exact field.
#[must_use]
pub fn l2_h1_error(
    mesh: &Mesh2,
    u: &[[f64; 2]],
    exact: &dyn Fn(f64, f64) -> [f64; 2],
    grad_exact: &dyn Fn(f64, f64) -> [[f64; 2]; 2],
) -> (f64, f64) {
    let mut l2 = 0.0f64;
    let mut h1 = 0.0f64;
    for conn in &mesh.elems {
        for &(xi, eta, w) in &quad_points(conn.len()) {
            let (shape, grads, det) = shapes_at(&mesh.nodes, conn, xi, eta);
            let wq = w * det;
            let mut p = [0.0f64; 2];
            let mut uh = [0.0f64; 2];
            let mut guh = [[0.0f64; 2]; 2];
            for (a, &node) in conn.iter().enumerate() {
                p[0] += shape[a] * mesh.nodes[node][0];
                p[1] += shape[a] * mesh.nodes[node][1];
                for c in 0..2 {
                    uh[c] += shape[a] * u[node][c];
                    guh[c][0] += grads[a][0] * u[node][c];
                    guh[c][1] += grads[a][1] * u[node][c];
                }
            }
            let ue = exact(p[0], p[1]);
            let ge = grad_exact(p[0], p[1]);
            for c in 0..2 {
                let e = ue[c] - uh[c];
                l2 += wq * e * e;
                for r in 0..2 {
                    let d = ge[c][r] - guh[c][r];
                    h1 += wq * d * d;
                }
            }
        }
    }
    (l2.max(0.0).sqrt(), h1.max(0.0).sqrt())
}
