//! Finite-strain P1/Q1 plane stress using the original three-dimensional law.
//!
//! Every quadrature point relaxes the transverse deformation through
//! `fs-material`; the residual and Schur-complement tangent differentiate
//! that SAME reduced energy. Reference geometry is prepared once and the
//! existing load-stepped Newton solver owns global equilibrium. Loads,
//! forces and energies use unit original thickness, not current thickness.
//! This is an in-plane continuum, with no shell bending or wrinkling model.

use super::{HyperProblem, NewtonReport, Reduction};
use crate::SolidError;
use crate::mesh2::{quad_points, shapes_at};
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::plane_stress::{PlaneStressSettings, PlaneStressState};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeMap;

fn invalid(what: impl Into<String>) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}

fn material_error(element: usize, error: fs_material::MaterialError) -> SolidError {
    SolidError::MaterialRefused {
        what: format!("plane-stress element {element}: {error}"),
    }
}

struct Point {
    gradients: [[f64; 2]; 4],
    weight: f64,
}

struct Element {
    nodes: [usize; 4],
    arity: usize,
    reference: [[f64; 2]; 4],
    points: Vec<Point>,
}

impl Element {
    // A bilinear quad's in-plane Jacobian is affine in reference coordinates:
    // positivity at every corner excludes a folded corner hidden from Gauss
    // points. A triangle has a constant Jacobian. Use relative motion so a
    // common displacement does not enter any deformation calculation.
    fn positive_orientation(&self, u: Option<&[f64]>) -> bool {
        (0..self.arity).all(|a| {
            let b = (a + 1) % self.arity;
            let c = (a + self.arity - 1) % self.arity;
            let edge = |v: usize, i: usize| {
                self.reference[v][i] - self.reference[a][i]
                    + u.map_or(0.0, |u| u[2 * self.nodes[v] + i] - u[2 * self.nodes[a] + i])
            };
            let cross = edge(b, 0) * edge(c, 1) - edge(b, 1) * edge(c, 0);
            cross.is_finite() && cross > 0.0
        })
    }
}

pub(super) struct Prepared {
    ndof: usize,
    elements: Vec<Element>,
    local: PlaneStressSettings,
}

impl Prepared {
    fn new(
        problem: &HyperProblem<'_>,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Self, SolidError> {
        checkpoint(cx)?;
        let mesh = problem.mesh;
        if mesh.nodes.is_empty()
            || mesh.elems.is_empty()
            || mesh.nodes.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid("plane stress needs a nonempty finite P1/Q1 mesh"));
        }
        let ndof = mesh
            .nodes
            .len()
            .checked_mul(2)
            .ok_or_else(|| invalid("plane-stress displacement count overflows"))?;
        // Validate public material fields and local controls before allocating
        // reference geometry. This also checks stability at the original state.
        Hyperelastic::new(
            problem.material.model.clone(),
            problem.material.stretch_limit,
        )
        .and_then(|law| law.plane_stress_state(&[1.0, 0.0, 0.0, 1.0], local))
        .map_err(|error| material_error(0, error))?;
        for (_, edges) in &mesh.patches {
            if edges
                .iter()
                .any(|&(a, b)| a >= mesh.nodes.len() || b >= mesh.nodes.len() || a == b)
            {
                return Err(invalid("plane-stress boundary edge has invalid nodes"));
            }
        }
        let mut elements = Vec::with_capacity(mesh.elems.len());
        for (index, conn) in mesh.elems.iter().enumerate() {
            checkpoint(cx)?;
            let arity = conn.len();
            if !(arity == 3 || arity == 4)
                || conn.iter().any(|&a| a >= mesh.nodes.len())
                || (0..arity).any(|a| (0..a).any(|b| conn[a] == conn[b]))
            {
                return Err(invalid(format!(
                    "plane-stress element {index} must be a valid P1 triangle or Q1 quad"
                )));
            }
            let mut element = Element {
                nodes: [0; 4],
                arity,
                reference: [[0.0; 2]; 4],
                points: Vec::new(),
            };
            element.nodes[..arity].copy_from_slice(conn);
            for (a, &node) in conn.iter().enumerate() {
                for i in 0..2 {
                    element.reference[a][i] = mesh.nodes[node][i] - mesh.nodes[conn[0]][i];
                }
            }
            if !element.positive_orientation(None) {
                return Err(invalid(format!(
                    "plane-stress element {index} is folded or reversed"
                )));
            }
            // P1 has constant F and one homogeneous law per element, so one
            // centroid evaluation integrates its reduced energy and tangent
            // exactly. Q1 retains all four Gauss points.
            let rule = if arity == 3 {
                vec![(1.0 / 3.0, 1.0 / 3.0, 0.5)]
            } else {
                quad_points(4)
            };
            for (xi, eta, w) in rule {
                let (_, gradients, det) =
                    shapes_at(&element.reference, &[0, 1, 2, 3][..arity], xi, eta);
                let weight = w * det;
                if !weight.is_finite()
                    || weight <= 0.0
                    || gradients.iter().flatten().any(|v| !v.is_finite())
                {
                    return Err(invalid(format!(
                        "plane-stress element {index} has unrepresentable geometry"
                    )));
                }
                let mut point = Point {
                    gradients: [[0.0; 2]; 4],
                    weight,
                };
                point.gradients[..arity].copy_from_slice(&gradients);
                element.points.push(point);
            }
            elements.push(element);
        }
        checkpoint(cx)?;
        Ok(Self {
            ndof,
            elements,
            local,
        })
    }

    fn for_points(
        &self,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
        mut visit: impl FnMut(usize, &Element, &Point, [f64; 4]) -> Result<(), SolidError>,
    ) -> Result<(), SolidError> {
        checkpoint(cx)?;
        if u.len() != self.ndof || !load.is_finite() || u.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "plane-stress displacement coverage and load must be finite and exact",
            ));
        }
        for (index, element) in self.elements.iter().enumerate() {
            if !element.positive_orientation(Some(u)) {
                return Err(SolidError::MaterialRefused {
                    what: format!("plane-stress element {index} is inverted or folded"),
                });
            }
            for point in &element.points {
                checkpoint(cx)?;
                let f = HyperProblem::def_grad(
                    &point.gradients[..element.arity],
                    &element.nodes[..element.arity],
                    u,
                );
                visit(index, element, point, [f[0], f[1], f[3], f[4]])?;
            }
        }
        checkpoint(cx)
    }

    pub(super) fn assemble(
        &self,
        problem: &HyperProblem<'_>,
        u: &[f64],
        load: f64,
        fixed: &BTreeMap<usize, f64>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Csr), SolidError> {
        let mut residual = vec![0.0; self.ndof];
        let mut coo = Coo::new(self.ndof, self.ndof);
        self.for_points(u, load, cx, |index, element, point, f| {
            let response = problem
                .material
                .plane_stress(&f, self.local)
                .map_err(|error| material_error(index, error))?;
            for (a, &na) in element.nodes[..element.arity].iter().enumerate() {
                for i in 0..2 {
                    let ia = 2 * na + i;
                    if fixed.contains_key(&ia) {
                        continue;
                    }
                    residual[ia] += point.weight
                        * (response.piola[2 * i] * point.gradients[a][0]
                            + response.piola[2 * i + 1] * point.gradients[a][1]);
                    for (b, &nb) in element.nodes[..element.arity].iter().enumerate() {
                        for k in 0..2 {
                            let ib = 2 * nb + k;
                            if fixed.contains_key(&ib) {
                                continue;
                            }
                            let mut value = 0.0;
                            for j in 0..2 {
                                for l in 0..2 {
                                    value += response.tangent[2 * i + j][2 * k + l]
                                        * point.gradients[a][j]
                                        * point.gradients[b][l];
                                }
                            }
                            coo.push(ia, ib, point.weight * value);
                        }
                    }
                }
            }
            Ok(())
        })?;
        problem.for_traction(load, |dof, value| {
            if !fixed.contains_key(&dof) {
                residual[dof] -= value;
            }
        })?;
        for &dof in fixed.keys() {
            coo.push(dof, dof, 1.0);
        }
        let tangent = coo.assemble();
        if residual.iter().any(|v| !v.is_finite())
            || (0..self.ndof).any(|row| tangent.row(row).1.iter().any(|v| !v.is_finite()))
        {
            return Err(invalid("unrepresentable plane-stress force or tangent"));
        }
        checkpoint(cx)?;
        Ok((residual, tangent))
    }

    pub(super) fn potential_energy(
        &self,
        problem: &HyperProblem<'_>,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, SolidError> {
        let mut potential = 0.0;
        self.for_points(u, load, cx, |index, _, point, f| {
            let state = problem
                .material
                .plane_stress_state(&f, self.local)
                .map_err(|error| material_error(index, error))?;
            potential += point.weight * state.energy_density_pa;
            Ok(())
        })?;
        problem.for_traction(load, |dof, value| potential -= value * u[dof])?;
        if !potential.is_finite() {
            return Err(invalid("unrepresentable plane-stress potential energy"));
        }
        checkpoint(cx)?;
        Ok(potential)
    }
}

impl HyperProblem<'_> {
    /// Solve a thin P1/Q1 hyperelastic continuum with traction-free thickness.
    ///
    /// Reuses the original supports, dead tractions, load steps, Newton and
    /// backtracking. Geometry is cached once; local constitutive Newton is
    /// bounded by `local`. Forces and energy are per unit reference thickness.
    /// Cancellation, invalid geometry, material-domain loss and either failed
    /// equilibrium solve return no partial displacement field.
    pub fn solve_plane_stress(
        &self,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        if !self.settings.tol.is_finite()
            || self.settings.tol <= 0.0
            || self.settings.max_iters == 0
            || self.settings.load_steps == 0
            || self.settings.max_backtracks == 0
        {
            return Err(invalid(
                "plane-stress Newton controls must be finite and positive",
            ));
        }
        let prepared = Prepared::new(self, local, cx)?;
        self.solve_with_projection(Some(Reduction::PlaneStress(&prepared, cx)))
    }

    /// Unconstrained force residual and consistent reduced tangent. At a solved
    /// state, support entries of the residual contain the reaction forces.
    pub fn plane_stress_residual_and_tangent(
        &self,
        u: &[f64],
        load: f64,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Csr), SolidError> {
        Prepared::new(self, local, cx)?.assemble(self, u, load, &BTreeMap::new(), cx)
    }

    /// Reduced total potential, including work of the original dead tractions.
    /// Uses the stress/energy-only local path without full material Hessians.
    pub fn plane_stress_potential_energy(
        &self,
        u: &[f64],
        load: f64,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<f64, SolidError> {
        Prepared::new(self, local, cx)?.potential_energy(self, u, load, cx)
    }

    /// Recover the actual relaxed three-dimensional states, in element order
    /// and then quadrature order (one exact centroid for P1, four for Q1).
    pub fn plane_stress_states(
        &self,
        u: &[f64],
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Vec<PlaneStressState>, SolidError> {
        let prepared = Prepared::new(self, local, cx)?;
        let mut states = Vec::new();
        prepared.for_points(u, 0.0, cx, |index, _, _, f| {
            states.push(
                self.material
                    .plane_stress_state(&f, local)
                    .map_err(|error| material_error(index, error))?,
            );
            Ok(())
        })?;
        Ok(states)
    }
}
