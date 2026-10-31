//! Finite-strain Q1 mean dilatation, with element pressure and dilatation
//! condensed out of the three-field energy:
//!
//! `E_e(u) = sum_q w_q W_iso(F_q(u)) + V_e U(Jbar(u))`,
//! `Jbar = sum_q w_q det(F_q) / V_e`.
//!
//! The original material owns BOTH energy terms. In particular, its
//! Neo-Hookean remainder and full-invariant fibre terms are retained as
//! authored; they are not silently replaced with a different isochoric law.
//! Nodal automatic differentiation includes the complete first and second
//! derivatives of Jbar, including geometric stiffness and cross-quadrature
//! coupling. Geometry is prepared once per solve. The standard formulation
//! remains the default.
//!
//! This is a quasi-incompressible plane-strain Q1 path, not a TDNNS element,
//! a fully incompressible solver, a P1 locking remedy, or an hourglass/
//! arbitrary-mesh stability guarantee. The deviatoric/remainder energy uses
//! all four Gauss points; local material admission is never averaged away.
//! The three-field condensation is described in the deal.II step-44 tutorial
//! and Simo, Taylor & Pister (1985), doi:10.1016/0045-7825(85)90033-7.

use super::{HyperProblem, NewtonReport};
use crate::SolidError;
use crate::mesh2::{quad_points, shapes_at};
use fs_ad::{Dual, Dual64, Real};
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_sparse::{Coo, Csr};
use std::collections::BTreeMap;

/// Element-level variables of the condensed volumetric field.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeanDilatationState {
    /// Reference area, in squared mesh-coordinate units (per unit thickness).
    pub reference_area: f64,
    /// Reference-volume-weighted mean of the physical local det F.
    pub mean_volume_ratio: f64,
    /// Smallest physical det F at the four Gauss points.
    pub minimum_quadrature_volume_ratio: f64,
    /// Largest physical det F at the four Gauss points.
    pub maximum_quadrature_volume_ratio: f64,
    /// `U'(Jbar)`, positive in tension, in the material's stress units.
    ///
    /// This contributes `U'(Jbar) cof(F_q)` to first Piola stress. It is
    /// ONLY the condensed bulk term; the original material remainder can
    /// also carry spherical stress. Compression-positive pressure is the
    /// negative of this value. Evaluating the original full material at
    /// local F would report a different volumetric field.
    pub volumetric_cauchy_stress: f64,
}

fn invalid(what: impl Into<String>) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

fn refused(what: impl Into<String>) -> SolidError {
    SolidError::MaterialRefused { what: what.into() }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}

#[derive(Clone, Copy)]
struct Point {
    weight: f64,
    gradients: [[f64; 2]; 4],
}

struct Element {
    nodes: [usize; 4],
    reference: [[f64; 2]; 4],
    points: [Point; 4],
    area: f64,
}

pub(super) struct Prepared {
    ndof: usize,
    elements: Vec<Element>,
}

impl Prepared {
    fn new(problem: &HyperProblem<'_>, cx: &Cx<'_>) -> Result<Self, SolidError> {
        checkpoint(cx)?;
        let mesh = problem.mesh;
        if mesh.nodes.is_empty() || mesh.elems.is_empty() {
            return Err(invalid("mean dilatation needs a nonempty Q1 mesh"));
        }
        let ndof = mesh
            .nodes
            .len()
            .checked_mul(2)
            .ok_or_else(|| invalid("mean-dilatation displacement count overflows"))?;
        if mesh.nodes.iter().flatten().any(|x| !x.is_finite()) {
            return Err(invalid("mean-dilatation coordinates must be finite"));
        }
        // Public material fields can have changed after construction.
        Hyperelastic::new(
            problem.material.model.clone(),
            problem.material.stretch_limit,
        )
        .map_err(|e| refused(format!("{e:?}")))?;
        for (_, edges) in &mesh.patches {
            if edges
                .iter()
                .any(|&(a, b)| a >= mesh.nodes.len() || b >= mesh.nodes.len() || a == b)
            {
                return Err(invalid("mean-dilatation boundary edge has invalid nodes"));
            }
        }
        let rule = quad_points(4);
        let mut elements = Vec::with_capacity(mesh.elems.len());
        for (index, conn) in mesh.elems.iter().enumerate() {
            checkpoint(cx)?;
            let nodes: [usize; 4] = conn.as_slice().try_into().map_err(|_| {
                invalid("mean dilatation supports Q1 quadrilaterals only; P1 remains unmodified")
            })?;
            if nodes.iter().any(|&a| a >= mesh.nodes.len())
                || (0..4).any(|a| (0..a).any(|b| nodes[a] == nodes[b]))
            {
                return Err(invalid(format!(
                    "mean-dilatation element {index} has invalid nodes"
                )));
            }
            let reference = nodes.map(|a| mesh.nodes[a]);
            if !positive_corners(&reference, &[0.0; 8]) {
                return Err(invalid(format!(
                    "mean-dilatation element {index} is folded or reversed"
                )));
            }
            let mut points = [Point {
                weight: 0.0,
                gradients: [[0.0; 2]; 4],
            }; 4];
            let mut area = 0.0;
            for (point, &(xi, eta, weight)) in points.iter_mut().zip(&rule) {
                let (_, gradients, det) = shapes_at(&mesh.nodes, conn, xi, eta);
                point.weight = weight * det;
                point.gradients.copy_from_slice(&gradients);
                if !point.weight.is_finite()
                    || point.weight <= 0.0
                    || point.gradients.iter().flatten().any(|x| !x.is_finite())
                {
                    return Err(invalid(format!(
                        "mean-dilatation element {index} has invalid geometry"
                    )));
                }
                area += point.weight;
            }
            if !area.is_finite() || area <= 0.0 {
                return Err(invalid(
                    "mean-dilatation reference area is not finite and positive",
                ));
            }
            elements.push(Element {
                nodes,
                reference,
                points,
                area,
            });
        }
        checkpoint(cx)?;
        Ok(Self { ndof, elements })
    }

    fn validate_state(&self, u: &[f64], load: f64) -> Result<(), SolidError> {
        if u.len() != self.ndof || !load.is_finite() || u.iter().any(|x| !x.is_finite()) {
            return Err(invalid(
                "mean-dilatation displacement coverage and load must be finite and exact",
            ));
        }
        Ok(())
    }

    pub(super) fn assemble(
        &self,
        problem: &HyperProblem<'_>,
        u: &[f64],
        load: f64,
        fixed: &BTreeMap<usize, f64>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Csr), SolidError> {
        checkpoint(cx)?;
        self.validate_state(u, load)?;
        let mut residual = vec![0.0; self.ndof];
        let mut coo = Coo::new(self.ndof, self.ndof);
        for element in &self.elements {
            checkpoint(cx)?;
            let local = element.admitted_motion(problem.material, u)?;
            let mut force = [0.0; 8];
            let mut stiffness = [[0.0; 8]; 8];
            for column in 0..8 {
                checkpoint(cx)?;
                let variables = core::array::from_fn(|i| Dual {
                    re: Dual64::<8>::variable(local[i], i),
                    eps: [Dual64::<8>::constant(if i == column { 1.0 } else { 0.0 })],
                });
                let energy = element.energy(problem.material, &variables);
                if column == 0 {
                    force = energy.re.eps;
                }
                for (row, values) in stiffness.iter_mut().enumerate() {
                    values[column] = energy.eps[0].eps[row];
                }
            }
            if force
                .iter()
                .chain(stiffness.iter().flatten())
                .any(|x| !x.is_finite())
            {
                return Err(refused("mean-dilatation element derivative is nonfinite"));
            }
            for a in 0..8 {
                let ia = 2 * element.nodes[a / 2] + a % 2;
                if fixed.contains_key(&ia) {
                    continue;
                }
                residual[ia] += force[a];
                for b in 0..8 {
                    let ib = 2 * element.nodes[b / 2] + b % 2;
                    if !fixed.contains_key(&ib) {
                        // Store one Hessian triangle identically in both
                        // orientations; nested-dual roundoff is not physics.
                        coo.push(ia, ib, stiffness[a.min(b)][a.max(b)]);
                    }
                }
            }
        }
        problem.for_traction(load, |dof, value| {
            if !fixed.contains_key(&dof) {
                residual[dof] -= value;
            }
        })?;
        for &dof in fixed.keys() {
            coo.push(dof, dof, 1.0);
        }
        if residual.iter().any(|x| !x.is_finite()) {
            return Err(invalid("mean-dilatation assembled force is nonfinite"));
        }
        let tangent = coo.assemble();
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
        checkpoint(cx)?;
        self.validate_state(u, load)?;
        let mut value = 0.0;
        for element in &self.elements {
            checkpoint(cx)?;
            let local = element.admitted_motion(problem.material, u)?;
            value += element.energy(problem.material, &local);
        }
        problem.for_traction(load, |dof, force| value -= force * u[dof])?;
        if !value.is_finite() {
            return Err(refused("mean-dilatation energy is nonfinite"));
        }
        checkpoint(cx)?;
        Ok(value)
    }
}

impl Element {
    fn admitted_motion(&self, material: &Hyperelastic, u: &[f64]) -> Result<[f64; 8], SolidError> {
        let local = core::array::from_fn(|i| u[2 * self.nodes[i / 2] + i % 2]);
        if !positive_corners(&self.reference, &local) {
            return Err(refused(
                "mean-dilatation current Q1 element folds at a corner",
            ));
        }
        for point in &self.points {
            // Admission belongs to the original material. In particular, a
            // good mean determinant cannot hide one inadmissible local F.
            material
                .piola(&deformation(&local, &point.gradients))
                .map_err(|e| refused(format!("{e:?}")))?;
        }
        Ok(local)
    }

    fn energy<T: Real>(&self, material: &Hyperelastic, u: &[T; 8]) -> T {
        let mut remainder = T::zero();
        let mut volume = T::zero();
        for point in &self.points {
            let f = deformation(u, &point.gradients);
            let weight = T::from_f64(point.weight);
            volume = volume + weight * (f[0] * f[4] - f[1] * f[3]);
            remainder = remainder + weight * material.isochoric_energy(&f);
        }
        let area = T::from_f64(self.area);
        remainder + area * material.volumetric_energy(volume / area)
    }

    fn state(
        &self,
        material: &Hyperelastic,
        u: &[f64; 8],
    ) -> Result<MeanDilatationState, SolidError> {
        let mut volume = 0.0;
        let mut minimum = f64::INFINITY;
        let mut maximum = f64::NEG_INFINITY;
        for point in &self.points {
            let f = deformation(u, &point.gradients);
            let j = f[0] * f[4] - f[1] * f[3];
            volume += point.weight * j;
            minimum = minimum.min(j);
            maximum = maximum.max(j);
        }
        let mean = volume / self.area;
        let stress = material
            .volumetric_energy(Dual64::<1>::variable(mean, 0))
            .eps[0];
        if !mean.is_finite() || !stress.is_finite() {
            return Err(refused("mean-dilatation volumetric state is nonfinite"));
        }
        Ok(MeanDilatationState {
            reference_area: self.area,
            mean_volume_ratio: mean,
            minimum_quadrature_volume_ratio: minimum,
            maximum_quadrature_volume_ratio: maximum,
            volumetric_cauchy_stress: stress,
        })
    }
}

fn deformation<T: Real>(u: &[T; 8], gradients: &[[f64; 2]; 4]) -> [T; 9] {
    let mut f = [T::zero(); 9];
    f[0] = T::one();
    f[4] = T::one();
    f[8] = T::one();
    for a in 1..4 {
        for i in 0..2 {
            let relative = u[2 * a + i] - u[i];
            for j in 0..2 {
                f[3 * i + j] = f[3 * i + j] + relative * T::from_f64(gradients[a][j]);
            }
        }
    }
    f
}

fn positive_corners(reference: &[[f64; 2]; 4], u: &[f64; 8]) -> bool {
    (0..4).all(|a| {
        let b = (a + 1) % 4;
        let c = (a + 3) % 4;
        let edge = |v: usize, i: usize| {
            (reference[v][i] - reference[a][i]) + (u[2 * v + i] - u[2 * a + i])
        };
        let cross = edge(b, 0) * edge(c, 1) - edge(b, 1) * edge(c, 0);
        cross.is_finite() && cross > 0.0
    })
}

impl HyperProblem<'_> {
    /// Solve finite-strain plane compression using condensed Q1 mean dilatation.
    ///
    /// Uses the existing loads, material and Newton controls. Reference
    /// quadrature is prepared once; cancellation is polled between elements,
    /// derivative sweeps, line-search trials and bounded Krylov chunks. No
    /// partial displacement is returned after refusal or cancellation.
    ///
    /// # Errors
    /// Refuses non-Q1/invalid meshes, invalid controls, inadmissible local
    /// material states, nonconvergence and [`SolidError::Cancelled`].
    pub fn solve_mean_dilatation(
        &self,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        let prepared = Prepared::new(self, cx)?;
        if !self.settings.tol.is_finite()
            || self.settings.tol <= 0.0
            || self.settings.max_iters == 0
            || self.settings.load_steps == 0
            || self.settings.max_backtracks == 0
        {
            return Err(invalid(
                "mean-dilatation Newton controls must be finite and positive",
            ));
        }
        self.solve_with_projection(Some((&prepared, cx)))
    }

    /// Unconstrained forces and exact condensed-energy Hessian, per unit thickness.
    /// Support entries of the returned residual give reactions at a solved state.
    pub fn mean_dilatation_residual_and_tangent(
        &self,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Csr), SolidError> {
        Prepared::new(self, cx)?.assemble(self, u, load, &BTreeMap::new(), cx)
    }

    /// Condensed total potential energy per unit thickness, with original dead loads.
    pub fn mean_dilatation_potential_energy(
        &self,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, SolidError> {
        Prepared::new(self, cx)?.potential_energy(self, u, load, cx)
    }

    /// Recover the actual condensed bulk field in element order.
    /// Local determinants and the mean are returned separately.
    pub fn mean_dilatation_states(
        &self,
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<MeanDilatationState>, SolidError> {
        let prepared = Prepared::new(self, cx)?;
        prepared.validate_state(u, 0.0)?;
        let mut states = Vec::with_capacity(prepared.elements.len());
        for element in &prepared.elements {
            checkpoint(cx)?;
            let local = element.admitted_motion(self.material, u)?;
            states.push(element.state(self.material, &local)?);
        }
        checkpoint(cx)?;
        Ok(states)
    }
}
