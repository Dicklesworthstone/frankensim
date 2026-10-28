//! Plane-strain finite-strain hyperelasticity: kinematics, weak forms,
//! and the globalized Newton loop. Constitutive behavior comes from
//! fs-material cards (exact AD Piola stress and consistent 9×9
//! tangent); this module embeds the 2D displacement gradient into the
//! 3D deformation gradient (F₃₃ = 1), integrates the residual
//! `∫ P : ∇δu`, and assembles the consistent tangent
//! `∫ ∇δu : A : ∇Δu`.
//!
//! Globalization: backtracking (Armijo) line search on the total
//! potential energy with dead loads, with material refusals
//! (det F ≤ 0 on a trial state) absorbed by step halving, plus load
//! stepping — the battery drives a buckling-adjacent fixture through
//! exactly this recipe.

use crate::SolidError;
use crate::mesh2::{Mesh2, Patch, quad_points, shapes_at};
use fs_material::hyper::Hyperelastic;
use fs_solver::krylov::MinresState;
use fs_solver::op::CsrOp;
use fs_sparse::Coo;
use std::collections::BTreeMap;

/// Newton controls.
#[derive(Debug, Clone, Copy)]
pub struct NewtonSettings {
    /// Residual-norm gate (absolute, on the free DOFs).
    pub tol: f64,
    /// Iteration cap per load step.
    pub max_iters: usize,
    /// Load steps (the Dirichlet data and tractions ramp linearly).
    pub load_steps: usize,
    /// Backtracking halvings before declaring a stall.
    pub max_backtracks: usize,
}

impl Default for NewtonSettings {
    fn default() -> Self {
        NewtonSettings {
            tol: 1e-10,
            max_iters: 30,
            load_steps: 1,
            max_backtracks: 25,
        }
    }
}

/// Convergence evidence for one solve (all load steps).
#[derive(Debug, Clone, Default)]
pub struct NewtonReport {
    /// Residual norms per iteration, per load step.
    pub histories: Vec<Vec<f64>>,
    /// Total line-search backtracks taken.
    pub backtracks: usize,
}

/// A displacement/traction field on a patch.
pub type VecField<'m> = &'m dyn Fn(f64, f64) -> [f64; 2];

/// A finite-strain problem on a body-fitted mesh (plane strain).
pub struct HyperProblem<'m> {
    /// The mesh.
    pub mesh: &'m Mesh2,
    /// The material card.
    pub material: &'m Hyperelastic,
    /// Dirichlet patches (ramped by the load factor).
    pub dirichlet: Vec<(Patch, VecField<'m>)>,
    /// Dead-load traction patches (ramped by the load factor).
    pub traction: Vec<(Patch, VecField<'m>)>,
    /// Newton controls.
    pub settings: NewtonSettings,
}

impl HyperProblem<'_> {
    /// Solve by load stepping + globalized Newton; returns nodal
    /// displacements and the convergence evidence.
    ///
    /// # Errors
    /// [`SolidError::NewtonStalled`], [`SolidError::SolveFailed`],
    /// [`SolidError::MaterialRefused`] (only if globalization cannot
    /// absorb a refusal), [`SolidError::UnknownPatch`].
    pub fn solve(&self) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        let n = self.mesh.node_count();
        let mut u = vec![0.0f64; 2 * n];
        let mut report = NewtonReport::default();
        for step in 1..=self.settings.load_steps {
            #[allow(clippy::cast_precision_loss)]
            let load = step as f64 / self.settings.load_steps as f64;
            let fixed = self.fixed_dofs(load)?;
            for (&dof, &val) in &fixed {
                u[dof] = val;
            }
            let mut history = Vec::new();
            let mut converged = false;
            for _ in 0..self.settings.max_iters {
                let (resid, tangent) = self.assemble(&u, load, &fixed)?;
                let rnorm = norm(&resid);
                history.push(rnorm);
                if rnorm < self.settings.tol {
                    converged = true;
                    break;
                }
                // Solve K δ = −R (MINRES: symmetric, indefinite-safe
                // near buckling).
                let b: Vec<f64> = resid.iter().map(|r| -r).collect();
                let op = CsrOp::symmetric(tangent);
                let mut st = MinresState::new(&op, &b);
                let _ = st.run(&op, 1e-10, 20_000);
                if !st.rel_residual().is_finite() || st.rel_residual() > 1e-6 {
                    return Err(SolidError::SolveFailed {
                        iters: st.iters,
                        rel_residual: st.rel_residual(),
                    });
                }
                // Armijo backtracking on the potential energy; material
                // refusals (det F ≤ 0) halve the step.
                let slope: f64 = resid.iter().zip(&st.x).map(|(r, d)| r * d).sum();
                let e0 = self.energy(&u, load);
                let mut alpha = 1.0f64;
                let mut taken = false;
                for _ in 0..self.settings.max_backtracks {
                    let trial: Vec<f64> = u
                        .iter()
                        .zip(&st.x)
                        .map(|(ui, di)| ui + alpha * di)
                        .collect();
                    if let Some(e1) = self.try_energy(&trial, load)
                        && (e1 <= e0 + 1e-4 * alpha * slope || rnorm < 1e-6)
                    {
                        u = trial;
                        taken = true;
                        break;
                    }
                    alpha *= 0.5;
                    report.backtracks += 1;
                }
                if !taken {
                    return Err(SolidError::NewtonStalled { history });
                }
            }
            if !converged {
                return Err(SolidError::NewtonStalled { history });
            }
            report.histories.push(history);
        }
        Ok(((0..n).map(|i| [u[2 * i], u[2 * i + 1]]).collect(), report))
    }

    /// Public probe for the consistency merge gate: residual and
    /// consistent tangent at an unconstrained state (no pinned DOFs).
    ///
    /// # Errors
    /// [`SolidError::MaterialRefused`] when det F ≤ 0 at `u`.
    pub fn residual_and_tangent(
        &self,
        u: &[f64],
        load: f64,
    ) -> Result<(Vec<f64>, fs_sparse::Csr), SolidError> {
        self.assemble(u, load, &BTreeMap::new())
    }

    /// Unconstrained force residual, without material tangents or sparse assembly.
    ///
    /// Contact line searches and refreshed friction lags need only forces.
    /// This uses the same Piola stresses, quadrature order, and dead tractions
    /// as [`Self::residual_and_tangent`], avoiding nine nested-dual sweeps per
    /// quadrature point when no Newton matrix will be consumed.
    ///
    /// # Errors
    /// [`SolidError::MaterialRefused`] when the material rejects `u`, or
    /// [`SolidError::UnknownPatch`] for an unknown traction patch.
    pub fn unconstrained_residual(&self, u: &[f64], load: f64) -> Result<Vec<f64>, SolidError> {
        let mut resid = vec![0.0; u.len()];
        for conn in &self.mesh.elems {
            for &(xi, eta, w) in &quad_points(conn.len()) {
                let (_, grads, det) = shapes_at(&self.mesh.nodes, conn, xi, eta);
                let wq = w * det;
                let f = Self::def_grad(&grads, conn, u);
                let p = self
                    .material
                    .piola(&f)
                    .map_err(|e| SolidError::MaterialRefused {
                        what: format!("{e:?}"),
                    })?;
                for (a, &na) in conn.iter().enumerate() {
                    for i in 0..2 {
                        resid[2 * na + i] +=
                            wq * (p[3 * i] * grads[a][0] + p[3 * i + 1] * grads[a][1]);
                    }
                }
            }
        }
        self.for_traction(load, |dof, v| resid[dof] -= v)?;
        Ok(resid)
    }

    /// Public probe: the total potential energy at a state (`None`
    /// when the material refuses it).
    #[must_use]
    pub fn potential_energy(&self, u: &[f64], load: f64) -> Option<f64> {
        self.try_energy(u, load)
    }

    fn fixed_dofs(&self, load: f64) -> Result<BTreeMap<usize, f64>, SolidError> {
        let mut fixed = BTreeMap::new();
        for (patch, g) in &self.dirichlet {
            if self.mesh.patch_edges(*patch).is_none() {
                return Err(SolidError::UnknownPatch { patch: *patch });
            }
            for node in self.mesh.patch_nodes(*patch) {
                let p = self.mesh.nodes[node];
                let val = g(p[0], p[1]);
                fixed.insert(2 * node, load * val[0]);
                fixed.insert(2 * node + 1, load * val[1]);
            }
        }
        Ok(fixed)
    }

    /// The 3D deformation gradient at a quadrature point (row-major,
    /// F₃₃ = 1 for plane strain).
    fn def_grad(grads: &[[f64; 2]], conn: &[usize], u: &[f64]) -> [f64; 9] {
        let mut h = [[0.0f64; 2]; 2];
        for (a, &node) in conn.iter().enumerate() {
            for i in 0..2 {
                h[i][0] += u[2 * node + i] * grads[a][0];
                h[i][1] += u[2 * node + i] * grads[a][1];
            }
        }
        [
            1.0 + h[0][0],
            h[0][1],
            0.0,
            h[1][0],
            1.0 + h[1][1],
            0.0,
            0.0,
            0.0,
            1.0,
        ]
    }

    /// Residual and consistent tangent on the free DOFs (fixed DOFs
    /// pinned with identity rows and zero residual).
    #[allow(clippy::type_complexity)]
    fn assemble(
        &self,
        u: &[f64],
        load: f64,
        fixed: &BTreeMap<usize, f64>,
    ) -> Result<(Vec<f64>, fs_sparse::Csr), SolidError> {
        let ndof = u.len();
        let mut resid = vec![0.0f64; ndof];
        let mut coo = Coo::new(ndof, ndof);
        for conn in &self.mesh.elems {
            let nn = conn.len();
            for &(xi, eta, w) in &quad_points(nn) {
                let (_, grads, det) = shapes_at(&self.mesh.nodes, conn, xi, eta);
                let wq = w * det;
                let f = Self::def_grad(&grads, conn, u);
                let p = self
                    .material
                    .piola(&f)
                    .map_err(|e| SolidError::MaterialRefused {
                        what: format!("{e:?}"),
                    })?;
                let a_mat = self
                    .material
                    .tangent(&f)
                    .map_err(|e| SolidError::MaterialRefused {
                        what: format!("{e:?}"),
                    })?;
                for (a, &na) in conn.iter().enumerate() {
                    for i in 0..2 {
                        let ia = 2 * na + i;
                        let ra = wq * (p[3 * i] * grads[a][0] + p[3 * i + 1] * grads[a][1]);
                        if !fixed.contains_key(&ia) {
                            resid[ia] += ra;
                        }
                        for (b, &nb) in conn.iter().enumerate() {
                            for k in 0..2 {
                                let ib = 2 * nb + k;
                                if fixed.contains_key(&ia) || fixed.contains_key(&ib) {
                                    continue;
                                }
                                let mut v = 0.0;
                                for jj in 0..2 {
                                    for ll in 0..2 {
                                        v += a_mat[3 * i + jj][3 * k + ll]
                                            * grads[a][jj]
                                            * grads[b][ll];
                                    }
                                }
                                coo.push(ia, ib, wq * v);
                            }
                        }
                    }
                }
            }
        }
        // Dead-load tractions (ramped).
        self.for_traction(load, |dof, v| {
            if !fixed.contains_key(&dof) {
                resid[dof] -= v;
            }
        })?;
        for &dof in fixed.keys() {
            coo.push(dof, dof, 1.0);
        }
        Ok((resid, coo.assemble()))
    }

    fn for_traction(&self, load: f64, mut sink: impl FnMut(usize, f64)) -> Result<(), SolidError> {
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
                    let w = 0.5 * len * load;
                    for (node, shape) in [(a, 1.0 - s), (b, s)] {
                        for (c, tc) in tv.iter().enumerate() {
                            sink(2 * node + c, w * shape * tc);
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Total potential energy (None when the material refuses the
    /// state — the line search treats that as "worse").
    fn try_energy(&self, u: &[f64], load: f64) -> Option<f64> {
        let mut internal = 0.0f64;
        for conn in &self.mesh.elems {
            for &(xi, eta, w) in &quad_points(conn.len()) {
                let (_, grads, det) = shapes_at(&self.mesh.nodes, conn, xi, eta);
                let f = Self::def_grad(&grads, conn, u);
                let j = f[0] * f[4] - f[1] * f[3];
                if !(j.is_finite() && j > 0.0) {
                    return None;
                }
                internal += w * det * self.material.energy(&f);
            }
        }
        let mut external = 0.0f64;
        let mut ok = true;
        let _ = self.for_traction(load, |dof, v| {
            if u.get(dof).is_some() {
                external += v * u[dof];
            } else {
                ok = false;
            }
        });
        ok.then_some(internal - external)
    }

    fn energy(&self, u: &[f64], load: f64) -> f64 {
        self.try_energy(u, load).unwrap_or(f64::INFINITY)
    }
}

/// Continuation plumbing: a [`HyperProblem`] with HOMOGENEOUS
/// Dirichlet data (clamps) is a [`crate::continuation::PathResidual`]
/// — the load parameter scales the dead tractions.
impl crate::continuation::PathResidual for HyperProblem<'_> {
    fn ndof(&self) -> usize {
        2 * self.mesh.node_count()
    }

    fn residual(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        let fixed = self.fixed_dofs(lambda)?;
        Ok(self.assemble(u, lambda, &fixed)?.0)
    }

    fn tangent(&self, u: &[f64], lambda: f64) -> Result<fs_sparse::Csr, SolidError> {
        let fixed = self.fixed_dofs(lambda)?;
        Ok(self.assemble(u, lambda, &fixed)?.1)
    }

    fn load_vector(&self) -> Vec<f64> {
        let mut f = vec![0.0f64; 2 * self.mesh.node_count()];
        let fixed = self.fixed_dofs(1.0).unwrap_or_default();
        let _ = self.for_traction(1.0, |dof, v| {
            if !fixed.contains_key(&dof) {
                f[dof] += v;
            }
        });
        f
    }
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::hyper::{FibreInvariants, HyperelasticModel, OgdenTerm, Volumetric};

    /// The cheaper force probe must preserve the fully assembled residual bits.
    #[test]
    fn residual_only_matches_all_materials_on_triangles_and_quads() {
        let models = [
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 5.0,
            },
            HyperelasticModel::MooneyRivlin {
                c10: 0.7,
                c01: 0.3,
                kappa: 8.0,
            },
            HyperelasticModel::Ogden {
                terms: vec![OgdenTerm {
                    mu_pa: 2.0,
                    alpha: 2.0,
                }],
                volumetric: Volumetric::LogSquared { kappa: 8.0 },
            },
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 0.7,
                k1: 1.2,
                k2: 2.0,
                dispersion: 0.1,
                fibres: vec![[1.0, 0.0, 0.0]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 8.0 },
            },
        ];
        let traction = |x: f64, y: f64| [0.2 + 0.1 * y, -0.15 * x];
        let prescribed = |_: f64, _: f64| [0.03, -0.01];
        for mesh in [
            Mesh2::quads(1.0, 0.7, 2, 2),
            Mesh2::triangles(1.0, 0.7, 2, 2),
        ] {
            for model in &models {
                let material = Hyperelastic::new(model.clone(), 3.0).unwrap();
                let problem = HyperProblem {
                    mesh: &mesh,
                    material: &material,
                    // Both public probes intentionally leave Dirichlet DOFs free.
                    dirichlet: vec![(Patch::Left, &prescribed)],
                    traction: vec![(Patch::Top, &traction), (Patch::Right, &traction)],
                    settings: NewtonSettings::default(),
                };
                for scale in [0.0, 0.1, -0.15] {
                    let u: Vec<f64> = mesh
                        .nodes
                        .iter()
                        .flat_map(|&[x, y]| {
                            [
                                scale * (0.2 * x + 0.1 * y + 0.15 * x * y),
                                scale * (-0.1 * x + 0.3 * y - 0.1 * x * y),
                            ]
                        })
                        .collect();
                    for load in [0.0, 0.6] {
                        let actual = problem.unconstrained_residual(&u, load).unwrap();
                        let expected = problem.residual_and_tangent(&u, load).unwrap().0;
                        assert_eq!(
                            actual.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                            expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                        );
                    }
                }
                let inverted: Vec<f64> = mesh
                    .nodes
                    .iter()
                    .flat_map(|&[x, _]| [-2.0 * x, 0.0])
                    .collect();
                assert!(matches!(
                    problem.unconstrained_residual(&inverted, 0.6),
                    Err(SolidError::MaterialRefused { .. })
                ));
                assert!(matches!(
                    problem.residual_and_tangent(&inverted, 0.6),
                    Err(SolidError::MaterialRefused { .. })
                ));
            }
        }
    }
}
