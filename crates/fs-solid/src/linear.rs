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
use fs_solver::checked_cg::{CheckedCgConfig, CheckedCgError, checked_cg};
use fs_solver::op::{CsrOp, LinearOp};
use fs_sparse::precond::Precond;
use fs_sparse::{Coo, Csr};
use std::collections::BTreeMap;

/// Original quadrature stress, elastic storage and physical support forces.
pub mod response;

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
    /// SYMMETRY constraints: (patch, component) pins ONE displacement
    /// component to zero along a patch, leaving the other free — the
    /// half/quarter-model boundary condition every NAFEMS LE-class
    /// benchmark needs (bead g42o). An explicit `dirichlet` entry on
    /// the same node takes precedence (symmetry never overwrites it).
    pub symmetry: Vec<(Patch, usize)>,
}

/// The Lamé pair for the plane reduction (plane stress uses the
/// standard effective λ).
#[must_use]
pub fn lame(youngs: f64, poisson: f64, plane: PlaneKind) -> (f64, f64) {
    let mu = 0.5 * youngs / (1.0 + poisson);
    match plane {
        PlaneKind::Strain => (
            youngs * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson)),
            mu,
        ),
        // Reduce analytically before evaluating: forming lambda*mu needlessly
        // squares the stress scale and overflows/underflows for finite moduli.
        PlaneKind::Stress => (youngs * poisson / ((1.0 + poisson) * (1.0 - poisson)), mu),
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
    /// Invalid material, geometry, supports or loads return `InvalidInput`.
    /// Unknown patches and missed equilibrium gates propagate separately.
    pub fn solve(&self) -> Result<Vec<[f64; 2]>, SolidError> {
        self.solve_scaled(None, &|| false)
    }

    /// Solve with a positive, dimensionless Young's-modulus multiplier per
    /// original element. Poisson ratio, body forces, tractions and prescribed
    /// displacements retain their original values. Fixed-column lifting uses
    /// the scaled stiffness, including nonzero prescribed displacements.
    ///
    /// # Errors
    /// Invalid scale count, nonpositive/nonfinite or unrepresentable elastic
    /// coefficients return `InvalidInput`; equilibrium failures propagate.
    pub fn solve_with_element_scales(&self, scales: &[f64]) -> Result<Vec<[f64; 2]>, SolidError> {
        self.solve_scaled(Some(scales), &|| false)
    }

    pub(crate) fn solve_scaled(
        &self,
        scales: Option<&[f64]>,
        cancelled: &impl Fn() -> bool,
    ) -> Result<Vec<[f64; 2]>, SolidError> {
        let checkpoint = || {
            if cancelled() {
                Err(SolidError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        let (lambda, mu) = self.validate(scales)?;
        let n = self.mesh.node_count();
        let ndof = 2 * n;
        let fixed = self.fixed_values()?;
        let mut coo = Coo::new(ndof, ndof);
        let mut rhs = vec![0.0f64; ndof];
        for (element, conn) in self.mesh.elems.iter().enumerate() {
            checkpoint()?;
            let scale = scales.map_or(1.0, |values| values[element]);
            let (k, fl) = self.element(conn, lambda * scale, mu * scale);
            if k.iter().flatten().chain(&fl).any(|v| !v.is_finite()) {
                return Err(SolidError::InvalidInput {
                    what: format!("element {element} has nonfinite stiffness or body load"),
                });
            }
            scatter_constrained(&mut coo, &mut rhs, conn, &k, &fl, &fixed);
        }
        self.traction_forces(&checkpoint, |dof, force| {
            if !fixed.contains_key(&dof) {
                rhs[dof] += force;
            }
        })?;
        // Fixed columns have already entered the physical free-force RHS.
        // Keep identity rows homogeneous during the solve: known displacement
        // values in metres must not enter a force residual in newtons. Restore
        // these prescribed values exactly after the free equilibrium is checked.
        for &dof in fixed.keys() {
            coo.push(dof, dof, 1.0);
            rhs[dof] = 0.0;
        }
        let a = coo.assemble();
        let m = Jacobi::new(&a);
        let op = CsrOp::symmetric(a);
        let mut solved = checked_cg(
            &op,
            &m,
            &rhs,
            CheckedCgConfig {
                tolerance: 1e-12,
                max_iterations: 40_000,
                max_corrections: 2,
            },
            |_| checkpoint(),
        )
        .map_err(|error| match error {
            CheckedCgError::InvalidInput(what) => SolidError::InvalidInput { what: what.into() },
            CheckedCgError::Interrupted(error) => error,
        })?;
        let rr = solved.report.rel_residual;
        if !rr.is_finite() || (rr > 1e-8 && !within_rounding_floor(&op, &rhs, &solved.x, 1e-8)) {
            return Err(SolidError::SolveFailed {
                iters: solved.report.iters,
                rel_residual: rr,
            });
        }
        for (&dof, &value) in &fixed {
            solved.x[dof] = value;
        }
        checkpoint()?;
        Ok((0..n)
            .map(|i| [solved.x[2 * i], solved.x[2 * i + 1]])
            .collect())
    }

    // Solving and observing must admit the same physical material and mesh.
    // In particular the ordinary path cannot bypass checks by omitting scales.
    fn validate(&self, scales: Option<&[f64]>) -> Result<(f64, f64), SolidError> {
        let (lambda, mu) = lame(self.youngs, self.poisson, self.plane);
        if scales.is_some_and(|s| s.len() != self.mesh.elems.len())
            || self.symmetry.iter().any(|(_, component)| *component >= 2)
            || self.mesh.nodes.iter().flatten().any(|v| !v.is_finite())
            || self.mesh.elems.iter().any(|conn| {
                crate::mesh2::checked_reference_element(&self.mesh.nodes, conn).is_none()
            })
            || self.mesh.patches.iter().any(|(_, edges)| {
                edges
                    .iter()
                    .any(|&(a, b)| a >= self.mesh.node_count() || b >= self.mesh.node_count())
            })
            || !(self.youngs.is_finite() && self.youngs > 0.0)
            || !(self.poisson.is_finite() && self.poisson > -1.0 && self.poisson < 0.5)
            || !(mu.is_finite() && mu > 0.0 && lambda.is_finite())
            || !(lambda + 2.0 * mu).is_finite()
            || (0..self.mesh.elems.len()).any(|element| {
                let s = scales.map_or(1.0, |values| values[element]);
                !s.is_finite()
                    || s <= 0.0
                    || !(mu * s > 0.0)
                    || !(lambda * s).is_finite()
                    || !((lambda + 2.0 * mu) * s).is_finite()
            })
        {
            return Err(SolidError::InvalidInput {
                what: "linear elasticity requires valid P1/Q1 geometry and supports, finite positive element multipliers and representable stable isotropic coefficients".into(),
            });
        }
        Ok((lambda, mu))
    }

    fn fixed_values(&self) -> Result<BTreeMap<usize, f64>, SolidError> {
        let mut fixed = BTreeMap::new();
        for (patch, g) in &self.dirichlet {
            if self.mesh.patch_edges(*patch).is_none() {
                return Err(SolidError::UnknownPatch { patch: *patch });
            }
            for node in self.mesh.patch_nodes(*patch) {
                let p = self.mesh.nodes[node];
                let val = g(p[0], p[1]);
                if val.iter().any(|value| !value.is_finite()) {
                    return Err(SolidError::InvalidInput {
                        what: format!("prescribed displacement at node {node} must be finite"),
                    });
                }
                fixed.insert(2 * node, val[0]);
                fixed.insert(2 * node + 1, val[1]);
            }
        }
        for (patch, comp) in &self.symmetry {
            if self.mesh.patch_edges(*patch).is_none() {
                return Err(SolidError::UnknownPatch { patch: *patch });
            }
            for node in self.mesh.patch_nodes(*patch) {
                // entry-or-insert: explicit Dirichlet wins on shared nodes.
                fixed.entry(2 * node + comp).or_insert(0.0);
            }
        }
        Ok(fixed)
    }

    // Keep full external forces available to the reaction observer, including
    // loads on prescribed rows. Elimination is the solve caller's decision.
    fn traction_forces(
        &self,
        checkpoint: &impl Fn() -> Result<(), SolidError>,
        mut visit: impl FnMut(usize, f64),
    ) -> Result<(), SolidError> {
        for (patch, t) in &self.traction {
            let edges = self
                .mesh
                .patch_edges(*patch)
                .ok_or(SolidError::UnknownPatch { patch: *patch })?;
            for &(a, b) in edges {
                checkpoint()?;
                let (pa, pb) = (self.mesh.nodes[a], self.mesh.nodes[b]);
                let len = (pb[0] - pa[0]).hypot(pb[1] - pa[1]);
                let g = 0.5 / 3.0f64.sqrt();
                for s in [0.5 - g, 0.5 + g] {
                    let p = [pa[0] + s * (pb[0] - pa[0]), pa[1] + s * (pb[1] - pa[1])];
                    let tv = t(p[0], p[1]);
                    if tv.iter().any(|v| !v.is_finite()) || !len.is_finite() {
                        return Err(SolidError::InvalidInput {
                            what: "boundary traction and its physical measure must be finite"
                                .into(),
                        });
                    }
                    let w = 0.5 * len;
                    for (node, shape) in [(a, 1.0 - s), (b, s)] {
                        for (c, tc) in tv.iter().enumerate() {
                            let force = w * shape * tc;
                            if !force.is_finite() {
                                return Err(SolidError::InvalidInput {
                                    what: "unrepresentable integrated boundary force".into(),
                                });
                            }
                            visit(2 * node + c, force);
                        }
                    }
                }
            }
        }
        Ok(())
    }

    // The same assumed strain enters stiffness, recovered stress and virtual
    // work. Keep the mean gradient independent of the displacement field.
    fn mean_gradients(&self, conn: &[usize]) -> Option<Vec<[f64; 2]>> {
        if self.formulation != Formulation::BBar {
            return None;
        }
        let mut mean = vec![[0.0; 2]; conn.len()];
        let mut area = 0.0;
        for (xi, eta, w) in quad_points(conn.len()) {
            let (_, gradients, determinant) = shapes_at(&self.mesh.nodes, conn, xi, eta);
            for (average, gradient) in mean.iter_mut().zip(gradients) {
                for i in 0..2 {
                    average[i] += w * determinant * gradient[i];
                }
            }
            area += w * determinant;
        }
        for gradient in &mut mean {
            for value in gradient {
                *value /= area;
            }
        }
        Some(mean)
    }

    fn point_body_force(&self, conn: &[usize], shape: &[f64]) -> Option<[f64; 2]> {
        self.body_force.map(|body| {
            let mut point = [0.0; 2];
            for (a, &node) in conn.iter().enumerate() {
                for (i, value) in point.iter_mut().enumerate() {
                    *value += shape[a] * self.mesh.nodes[node][i];
                }
            }
            body(point[0], point[1])
        })
    }

    /// One element's stiffness and body-force load (B-matrix form,
    /// Voigt [εxx, εyy, γxy]).
    pub(crate) fn element(
        &self,
        conn: &[usize],
        lambda: f64,
        mu: f64,
    ) -> (Vec<Vec<f64>>, Vec<f64>) {
        let nn = conn.len();
        let d = [
            [lambda + 2.0 * mu, lambda, 0.0],
            [lambda, lambda + 2.0 * mu, 0.0],
            [0.0, 0.0, mu],
        ];
        let pts = quad_points(nn);
        let bbar = self.mean_gradients(conn);
        let mut k = vec![vec![0.0f64; 2 * nn]; 2 * nn];
        let mut fl = vec![0.0f64; 2 * nn];
        for &(xi, eta, w) in &pts {
            let (shape, grads, det) = shapes_at(&self.mesh.nodes, conn, xi, eta);
            let wq = w * det;
            for a in 0..nn {
                for ca in 0..2 {
                    let ba = strain_column(&grads, bbar.as_deref(), a, ca);
                    let db = [
                        d[0][0] * ba[0] + d[0][1] * ba[1],
                        d[1][0] * ba[0] + d[1][1] * ba[1],
                        d[2][2] * ba[2],
                    ];
                    for b in 0..nn {
                        for cb in 0..2 {
                            let bb = strain_column(&grads, bbar.as_deref(), b, cb);
                            k[2 * a + ca][2 * b + cb] +=
                                wq * (db[0] * bb[0] + db[1] * bb[1] + db[2] * bb[2]);
                        }
                    }
                }
            }
            if let Some(f) = self.point_body_force(conn, &shape) {
                for a in 0..nn {
                    fl[2 * a] += wq * shape[a] * f[0];
                    fl[2 * a + 1] += wq * shape[a] * f[1];
                }
            }
        }
        (k, fl)
    }
}

/// The plain gate `||b - Ax|| <= tol ||b||` can ask a stiff, nearly incompressible
/// system for less than the rounding error of computing `b - Ax`. On its refusal
/// path a residual inside `tol ||b|| + gamma_{m+1} (||b|| + || |A||x| ||)`,
/// `gamma_k = k u / (1 - k u)` (Higham, 2nd ed., section 3.5; the gate of
/// fs-topopt's sdf3 solves since f834daca9d), cannot be told from the residual of
/// an exact solution and is admitted. `m` and `|A||x|` are the operator's own
/// (`LinearOp::apply_abs`); an operator that states neither keeps the plain gate.
/// This is not a bound on the error of the displacement.
fn within_rounding_floor(op: &impl LinearOp, rhs: &[f64], x: &[f64], tol: f64) -> bool {
    let mut applied = vec![0.0; rhs.len()];
    op.apply(x, &mut applied);
    let mut absolute = vec![0.0; rhs.len()];
    let Some(chain) = op.apply_abs(x, &mut absolute) else {
        return false;
    };
    #[allow(clippy::cast_precision_loss)]
    let unit = (chain as f64 + 1.0) * (f64::EPSILON / 2.0);
    let scale = rhs
        .iter()
        .chain(&absolute)
        .fold(0.0_f64, |s, v| s.max(v.abs()));
    if unit >= 1.0 || !(scale > 0.0 && scale.is_finite()) {
        return false;
    }
    let norm = |values: &mut dyn Iterator<Item = f64>| {
        values
            .map(|v| (v / scale) * (v / scale))
            .sum::<f64>()
            .sqrt()
    };
    let residual = norm(&mut rhs.iter().zip(&applied).map(|(b, a)| b - a));
    let load = norm(&mut rhs.iter().copied());
    let floor = unit / (1.0 - unit) * (load + norm(&mut absolute.iter().copied()));
    residual.is_finite() && floor.is_finite() && residual <= tol * load + floor
}

// Engineering-shear column [epsilon_xx, epsilon_yy, gamma_xy].
fn strain_column(
    gradients: &[[f64; 2]],
    mean: Option<&[[f64; 2]]>,
    node: usize,
    component: usize,
) -> [f64; 3] {
    let g = gradients[node];
    let mut b = if component == 0 {
        [g[0], 0.0, g[1]]
    } else {
        [0.0, g[1], g[0]]
    };
    if let Some(mean) = mean {
        let correction = 0.5 * (mean[node][component] - g[component]);
        b[0] += correction;
        b[1] += correction;
    }
    b
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

#[cfg(test)]
mod tests {
    use super::*;

    fn relative_free_force_residual(
        problem: &LinearProblem<'_>,
        u: &[[f64; 2]],
        traction: f64,
    ) -> f64 {
        let (lambda, mu) = lame(problem.youngs, problem.poisson, problem.plane);
        let mut force = vec![0.0; 2 * u.len()];
        for conn in &problem.mesh.elems {
            let (k, _) = problem.element(conn, lambda, mu);
            for (a, &node_a) in conn.iter().enumerate() {
                for i in 0..2 {
                    for (b, &node_b) in conn.iter().enumerate() {
                        for j in 0..2 {
                            force[2 * node_a + i] += k[2 * a + i][2 * b + j] * u[node_b][j];
                        }
                    }
                }
            }
        }
        for &(a, b) in problem.mesh.patch_edges(Patch::Right).unwrap() {
            let [ax, ay] = problem.mesh.nodes[a];
            let [bx, by] = problem.mesh.nodes[b];
            let half_edge_force = 0.5 * (bx - ax).hypot(by - ay) * traction;
            force[2 * a] -= half_edge_force;
            force[2 * b] -= half_edge_force;
        }
        for (patch, _) in &problem.dirichlet {
            for node in problem.mesh.patch_nodes(*patch) {
                force[2 * node] = 0.0;
                force[2 * node + 1] = 0.0;
            }
        }
        for &(patch, component) in &problem.symmetry {
            for node in problem.mesh.patch_nodes(patch) {
                force[2 * node + component] = 0.0;
            }
        }
        force
            .into_iter()
            .fold(0.0_f64, |norm, value| norm.hypot(value / traction))
    }

    // A traction-loaded rectangle has an exact affine solution. Scaling every
    // stress-like quantity leaves the strain unchanged, so even the large-force
    // checks stay inside small-strain elasticity. Test both signs and the native
    // element/formulation choices without comparing two solver implementations.
    #[test]
    fn g1_static_elasticity_preserves_strain_across_signed_stress_scales() {
        for mesh in [
            Mesh2::quads(2.0, 1.0, 3, 2),
            Mesh2::triangles(2.0, 1.0, 3, 2),
        ] {
            for plane in [PlaneKind::Stress, PlaneKind::Strain] {
                for formulation in [Formulation::Standard, Formulation::BBar] {
                    for stress_scale in [1.0, 1e-180, 1e180] {
                        for sign in [-1.0, 1.0] {
                            let youngs = 2.4e6 * stress_scale;
                            let poisson = 0.25;
                            let traction = sign * 1e-3 * youngs;
                            let load = |_: f64, _: f64| [traction, 0.0];
                            let problem = LinearProblem {
                                mesh: &mesh,
                                youngs,
                                poisson,
                                plane,
                                formulation,
                                body_force: None,
                                dirichlet: vec![],
                                traction: vec![(Patch::Right, &load)],
                                symmetry: vec![(Patch::Left, 0), (Patch::Bottom, 1)],
                            };
                            let u = problem.solve().expect("finite small-strain traction solve");
                            assert!(relative_free_force_residual(&problem, &u, traction) < 1e-9);
                            let strain = match plane {
                                PlaneKind::Stress => [sign * 1e-3, -poisson * sign * 1e-3],
                                PlaneKind::Strain => [
                                    (1.0 - poisson * poisson) * sign * 1e-3,
                                    -poisson * (1.0 + poisson) * sign * 1e-3,
                                ],
                            };
                            for (point, actual) in mesh.nodes.iter().zip(&u) {
                                for c in 0..2 {
                                    let expected = strain[c] * point[c];
                                    assert!(
                                        (actual[c] - expected).abs() < 2e-12,
                                        "{plane:?}/{formulation:?}, scale {stress_scale:e}, \
                                         sign {sign}, component {c}: {} != {expected}",
                                        actual[c]
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn g1_nonzero_supports_do_not_mask_free_force_equilibrium() {
        let mesh = Mesh2::quads(2.0, 1.0, 3, 2);
        for stress_scale in [1.0, 1e-180, 1e180] {
            let youngs = 2.4e6 * stress_scale;
            let extension = -1e-3;
            let poisson = 0.25;
            let translation = 0.002;
            let support = |_: f64, y: f64| [translation, -poisson * extension * y];
            let load = |_: f64, _: f64| [youngs * extension, 0.0];
            let problem = LinearProblem {
                mesh: &mesh,
                youngs,
                poisson,
                plane: PlaneKind::Stress,
                formulation: Formulation::Standard,
                body_force: None,
                dirichlet: vec![(Patch::Left, &support)],
                traction: vec![(Patch::Right, &load)],
                symmetry: vec![],
            };
            let u = problem
                .solve()
                .expect("finite force and prescribed-displacement solve");
            assert!(relative_free_force_residual(&problem, &u, youngs * extension) < 1e-9);
            for (point, actual) in mesh.nodes.iter().zip(&u) {
                let expected = [
                    translation + extension * point[0],
                    -poisson * extension * point[1],
                ];
                for c in 0..2 {
                    assert!(
                        (actual[c] - expected[c]).abs() < 2e-12,
                        "scale {stress_scale:e}, component {c}: {} != {}",
                        actual[c],
                        expected[c]
                    );
                }
                if point[0] == 0.0 {
                    assert_eq!(*actual, support(point[0], point[1]));
                }
            }
        }
    }

    #[test]
    fn g3_tiny_traction_remains_linear_and_zero_load_remains_exact() {
        let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
        for amplitude in [0.0, 1e-180, -1e-180, 1e-80, -1.0, 1.0] {
            let load = |_: f64, _: f64| [amplitude, 0.0];
            let u = LinearProblem {
                mesh: &mesh,
                youngs: 1000.0,
                poisson: 0.0,
                plane: PlaneKind::Stress,
                formulation: Formulation::Standard,
                body_force: None,
                dirichlet: vec![],
                traction: vec![(Patch::Right, &load)],
                symmetry: vec![(Patch::Left, 0), (Patch::Bottom, 1)],
            }
            .solve()
            .expect("nonzero finite traction cannot disappear");
            for (point, actual) in mesh.nodes.iter().zip(&u) {
                if amplitude == 0.0 {
                    assert_eq!(*actual, [0.0; 2]);
                } else {
                    assert!(
                        (actual[0] / amplitude - point[0] / 1000.0).abs() < 2e-12,
                        "amplitude {amplitude:e}: scaled displacement {}",
                        actual[0] / amplitude
                    );
                    assert!((actual[1] / amplitude).abs() < 2e-12);
                }
            }
        }
    }

    #[test]
    fn g0_nonfinite_prescriptions_refuse_even_when_every_node_is_fixed() {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let invalid = |_: f64, _: f64| [f64::NAN, 0.0];
        let problem = LinearProblem {
            mesh: &mesh,
            youngs: 1000.0,
            poisson: 0.25,
            plane: PlaneKind::Stress,
            formulation: Formulation::Standard,
            body_force: None,
            dirichlet: vec![(Patch::Left, &invalid), (Patch::Right, &invalid)],
            traction: vec![],
            symmetry: vec![],
        };
        assert!(matches!(
            problem.solve(),
            Err(SolidError::InvalidInput { .. })
        ));
    }
}
