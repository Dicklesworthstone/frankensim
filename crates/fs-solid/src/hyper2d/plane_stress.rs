//! Finite-strain P1/Q1 plane stress using the original three-dimensional law.
//!
//! Every quadrature point relaxes the transverse deformation through
//! `fs-material`; the residual and Schur-complement tangent differentiate
//! that SAME reduced energy. Reference geometry is prepared once and the
//! existing load-stepped Newton solver owns global equilibrium. Loads,
//! forces and energies use unit original thickness, not current thickness.
//! The crate-private motion entry scales this original continuum by its
//! declared thickness and adds consistent physical inertia and nodal forces.
//! This is an in-plane continuum, with no shell bending or wrinkling model.

use super::{
    HyperProblem, NewtonReport, PlaneMaterial, PlaneMaterialPoint, Reduction, material_quadrature,
};
use crate::SolidError;
use crate::inertia::{BackwardEulerInertia, InertiaError};
use crate::mesh2::{Mesh2, checked_reference_element, shapes_at};
use fs_exec::Cx;
use fs_material::hyper::incompressible_plane_stress::{
    self, IncompressiblePlaneStressResponse, IncompressiblePlaneStressState,
};
use fs_material::hyper::plane_stress::{
    self, PlaneStressResponse, PlaneStressSettings, PlaneStressState,
};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeMap;

#[cfg(test)]
mod dynamics_tests;

mod motion;
pub use motion::PlaneMaterialPointMotion;

fn invalid(what: impl Into<String>) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}

fn inertia_error(error: InertiaError) -> SolidError {
    match error {
        InertiaError::InvalidInput { what } => SolidError::InvalidInput { what },
        InertiaError::Cancelled => SolidError::Cancelled,
    }
}

fn material_error(element: usize, error: fs_material::MaterialError) -> SolidError {
    SolidError::MaterialRefused {
        what: format!("plane-stress element {element}: {error}"),
    }
}

fn reference_seed(material: &dyn PlaneMaterial, f: &[f64; 4]) -> [f64; 3] {
    // This is an initial guess only: the finite-bulk solve still relaxes
    // volume. Keeping the owner's natural volume avoids rejecting an expanded
    // material merely because unit total thickness lies outside its elastic
    // domain. The reducer retains its unit-thickness fallback for soft bulk.
    [
        0.0,
        0.0,
        material.incompressible_reference_volume_ratio() / (f[0] * f[3] - f[1] * f[2]),
    ]
}

fn finite_state(
    material: &dyn PlaneMaterial,
    f: &[f64; 4],
    settings: PlaneStressSettings,
    seed: Option<[f64; 3]>,
) -> Result<PlaneStressState, fs_material::MaterialError> {
    if let Some(seed) = seed {
        return plane_stress::plane_stress_state_with_seed(material, f, seed, settings);
    }
    if let Some(native) = material.native_hyperelastic() {
        native.plane_stress_state(f, settings)
    } else {
        plane_stress::plane_stress_state_with_seed(
            material,
            f,
            reference_seed(material, f),
            settings,
        )
    }
}

fn finite_response(
    material: &dyn PlaneMaterial,
    f: &[f64; 4],
    settings: PlaneStressSettings,
    seed: Option<[f64; 3]>,
) -> Result<PlaneStressResponse, fs_material::MaterialError> {
    if let Some(seed) = seed {
        return plane_stress::plane_stress_with_seed(material, f, seed, settings);
    }
    if let Some(native) = material.native_hyperelastic() {
        native.plane_stress(f, settings)
    } else {
        plane_stress::plane_stress_with_seed(material, f, reference_seed(material, f), settings)
    }
}

fn exact_state(
    material: &dyn PlaneMaterial,
    f: &[f64; 4],
    settings: PlaneStressSettings,
    seed: Option<[f64; 3]>,
) -> Result<IncompressiblePlaneStressState, fs_material::MaterialError> {
    if let Some(seed) = seed {
        return incompressible_plane_stress::incompressible_plane_stress_state_with_seed(
            material,
            f,
            [seed[0], seed[1]],
            settings,
        );
    }
    if let Some(native) = material.native_hyperelastic() {
        native.incompressible_plane_stress_state(f, settings)
    } else {
        incompressible_plane_stress::incompressible_plane_stress_state(material, f, settings)
    }
}

fn exact_response(
    material: &dyn PlaneMaterial,
    f: &[f64; 4],
    settings: PlaneStressSettings,
    seed: Option<[f64; 3]>,
) -> Result<IncompressiblePlaneStressResponse, fs_material::MaterialError> {
    if let Some(seed) = seed {
        return incompressible_plane_stress::incompressible_plane_stress_with_seed(
            material,
            f,
            [seed[0], seed[1]],
            settings,
        );
    }
    if let Some(native) = material.native_hyperelastic() {
        native.incompressible_plane_stress(f, settings)
    } else {
        incompressible_plane_stress::incompressible_plane_stress(material, f, settings)
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

#[derive(Clone, Copy)]
enum Constraint {
    FiniteBulk,
    Incompressible,
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

struct Geometry {
    ndof: usize,
    elements: Vec<Element>,
}

pub(super) struct Prepared<'a> {
    geometry: Geometry,
    local: PlaneStressSettings,
    constraint: Constraint,
    reference_thickness_m: f64,
    inertia: Option<&'a BackwardEulerInertia<'a>>,
    nodal_forces_n: &'a [(usize, f64)],
}

impl Geometry {
    fn new(mesh: &Mesh2, cx: &Cx<'_>) -> Result<Self, SolidError> {
        checkpoint(cx)?;
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
            let reference = checked_reference_element(&mesh.nodes, conn).ok_or_else(|| {
                invalid(format!(
                    "plane-stress element {index} must be an oriented P1 triangle or Q1 quad"
                ))
            })?;
            let mut element = Element {
                nodes: [0; 4],
                arity,
                reference,
                points: Vec::new(),
            };
            element.nodes[..arity].copy_from_slice(conn);
            for (xi, eta, w) in material_quadrature(arity) {
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
        Ok(Self { ndof, elements })
    }

    fn for_points(
        &self,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
        mut visit: impl FnMut(usize, usize, &Element, &Point, [f64; 4]) -> Result<(), SolidError>,
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
            for (local_point, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let f = HyperProblem::def_grad(
                    &point.gradients[..element.arity],
                    &element.nodes[..element.arity],
                    u,
                );
                let in_plane = [f[0], f[1], f[3], f[4]];
                if in_plane.iter().any(|value| !value.is_finite()) {
                    return Err(SolidError::MaterialRefused {
                        what: format!(
                            "element {index}, point {local_point} deformation is nonfinite"
                        ),
                    });
                }
                visit(index, local_point, element, point, in_plane)?;
            }
        }
        checkpoint(cx)
    }
}

pub(super) fn material_point_kinematics(
    mesh: &Mesh2,
    displacement: &[[f64; 2]],
    cx: &Cx<'_>,
) -> Result<Vec<PlaneMaterialPoint>, SolidError> {
    checkpoint(cx)?;
    if displacement.len() != mesh.nodes.len() {
        return Err(invalid(
            "material-point displacement must cover every mesh node",
        ));
    }
    let mut flat = Vec::with_capacity(2 * displacement.len());
    for (index, value) in displacement.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if value.iter().any(|value| !value.is_finite()) {
            return Err(invalid("material-point displacement must be finite"));
        }
        flat.extend_from_slice(value);
    }
    let geometry = Geometry::new(mesh, cx)?;
    let mut points = Vec::new();
    geometry.for_points(
        &flat,
        0.0,
        cx,
        |element, local_point, _, point, in_plane| {
            points.push(PlaneMaterialPoint {
                element,
                local_point,
                reference_area_m2: point.weight,
                in_plane,
            });
            Ok(())
        },
    )?;
    Ok(points)
}

impl<'a> Prepared<'a> {
    fn new(
        problem: &HyperProblem<'_>,
        local: PlaneStressSettings,
        constraint: Constraint,
        cx: &Cx<'_>,
    ) -> Result<Self, SolidError> {
        checkpoint(cx)?;
        // Identity may be outside a prescribed-reference/history domain.
        // Admit parameters here and actual states at their own quadrature points.
        problem.validate_materials_with(|| checkpoint(cx))?;
        Ok(Self {
            geometry: Geometry::new(problem.mesh, cx)?,
            local,
            constraint,
            reference_thickness_m: 1.0,
            inertia: None,
            nodal_forces_n: &[],
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn dynamic(
        problem: &HyperProblem<'_>,
        local: PlaneStressSettings,
        exact: bool,
        reference_thickness_m: f64,
        inertia: &'a BackwardEulerInertia<'a>,
        nodal_forces_n: &'a [(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<Self, SolidError> {
        if !reference_thickness_m.is_finite()
            || reference_thickness_m <= 0.0
            || problem.settings.load_steps != 1
        {
            return Err(invalid(
                "plane motion needs positive original thickness and one physical endpoint",
            ));
        }
        let constraint = if exact {
            Constraint::Incompressible
        } else {
            Constraint::FiniteBulk
        };
        let mut prepared = Self::new(problem, local, constraint, cx)?;
        for (index, &(dof, force)) in nodal_forces_n.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if dof >= prepared.geometry.ndof || !force.is_finite() {
                return Err(invalid(
                    "physical nodal force needs an admitted DOF and finite value",
                ));
            }
        }
        prepared.reference_thickness_m = reference_thickness_m;
        prepared.inertia = Some(inertia);
        prepared.nodal_forces_n = nodal_forces_n;
        Ok(prepared)
    }

    /// Complete the caller's free dead-traction vector in physical force units.
    /// Point/body forces are endpoint values and are never multiplied by the
    /// load ramp or by thickness a second time. Duplicate DOFs accumulate.
    pub(super) fn complete_applied_force(
        &self,
        applied: &mut [f64],
        fixed: &BTreeMap<usize, f64>,
        cx: &Cx<'_>,
    ) -> Result<(), SolidError> {
        for (dof, value) in applied.iter_mut().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            *value *= self.reference_thickness_m;
        }
        for (index, &(dof, value)) in self.nodal_forces_n.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            if !fixed.contains_key(&dof) {
                applied[dof] += value;
            }
        }
        Ok(())
    }

    pub(super) fn is_dynamic(&self) -> bool {
        self.inertia.is_some()
    }

    /// Use momentum and actual free forces, never stiffness times rigid travel.
    /// The same normalization is used for current and globalization trial states.
    pub(super) fn dynamic_force_scale(
        &self,
        u: &[f64],
        residual: &[f64],
        external: &[f64],
        free: &[usize],
        cx: &Cx<'_>,
    ) -> Result<Option<f64>, SolidError> {
        let Some(inertia) = self.inertia else {
            return Ok(None);
        };
        let mut inertial = vec![0.0; self.geometry.ndof];
        inertia
            .add_force(u, &mut inertial, cx)
            .map_err(inertia_error)?;
        let mut elastic_norm = 0.0_f64;
        let mut external_norm = 0.0_f64;
        for (index, &dof) in free.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            elastic_norm = elastic_norm.hypot(residual[dof] - inertial[dof] + external[dof]);
            external_norm = external_norm.hypot(external[dof]);
        }
        inertia
            .force_scale(u, free, elastic_norm, external_norm, cx)
            .map(Some)
            .map_err(inertia_error)
    }

    fn for_points(
        &self,
        problem: &HyperProblem<'_>,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
        mut visit: impl FnMut(
            usize,
            &Element,
            &Point,
            &dyn PlaneMaterial,
            [f64; 4],
            Option<[f64; 3]>,
        ) -> Result<(), SolidError>,
    ) -> Result<(), SolidError> {
        self.geometry.for_points(u, load, cx, |index, local_point, element, point, f| {
            let material = problem.point_material(index, local_point)?;
            let seed = problem.material.transverse_seed_at(index, local_point);
            if let Some(seed) = seed
                && (seed.iter().any(|value| !value.is_finite()) || seed[2] <= 0.0)
            {
                return Err(invalid(format!(
                    "element {index}, point {local_point} needs a finite positive-thickness transverse seed"
                )));
            }
            visit(index, element, point, material, f, seed)
        })
    }

    pub(super) fn assemble(
        &self,
        problem: &HyperProblem<'_>,
        u: &[f64],
        load: f64,
        fixed: &BTreeMap<usize, f64>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Csr), SolidError> {
        let mut residual = vec![0.0; self.geometry.ndof];
        let mut coo = Coo::new(self.geometry.ndof, self.geometry.ndof);
        self.for_points(
            problem,
            u,
            load,
            cx,
            |index, element, point, material, f, seed| {
                let (piola, tangent) = match self.constraint {
                    Constraint::FiniteBulk => finite_response(material, &f, self.local, seed)
                        .map(|response| (response.piola, response.tangent)),
                    Constraint::Incompressible => exact_response(material, &f, self.local, seed)
                        .map(|response| (response.piola, response.tangent)),
                }
                .map_err(|error| material_error(index, error))?;
                let weight = point.weight * self.reference_thickness_m;
                for (a, &na) in element.nodes[..element.arity].iter().enumerate() {
                    for i in 0..2 {
                        let ia = 2 * na + i;
                        if fixed.contains_key(&ia) {
                            continue;
                        }
                        residual[ia] += weight
                            * (piola[2 * i] * point.gradients[a][0]
                                + piola[2 * i + 1] * point.gradients[a][1]);
                        for (b, &nb) in element.nodes[..element.arity].iter().enumerate() {
                            for k in 0..2 {
                                let ib = 2 * nb + k;
                                if fixed.contains_key(&ib) {
                                    continue;
                                }
                                let mut value = 0.0;
                                for j in 0..2 {
                                    for l in 0..2 {
                                        value += tangent[2 * i + j][2 * k + l]
                                            * point.gradients[a][j]
                                            * point.gradients[b][l];
                                    }
                                }
                                coo.push(ia, ib, weight * value);
                            }
                        }
                    }
                }
                Ok(())
            },
        )?;
        problem.for_traction(load, |dof, value| {
            if !fixed.contains_key(&dof) {
                residual[dof] -= self.reference_thickness_m * value;
            }
        })?;
        for (index, &(dof, force)) in self.nodal_forces_n.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            residual[dof] -= force;
        }
        if let Some(inertia) = self.inertia {
            inertia
                .add_force(u, &mut residual, cx)
                .map_err(inertia_error)?;
            inertia
                .add_tangent_where(&mut coo, |dof| !fixed.contains_key(&dof), cx)
                .map_err(inertia_error)?;
        }
        for &dof in fixed.keys() {
            residual[dof] = 0.0;
            coo.push(dof, dof, 1.0);
        }
        let tangent = coo.assemble();
        if residual.iter().any(|v| !v.is_finite())
            || (0..self.geometry.ndof).any(|row| tangent.row(row).1.iter().any(|v| !v.is_finite()))
        {
            return Err(invalid("unrepresentable plane-stress force or tangent"));
        }
        checkpoint(cx)?;
        Ok((residual, tangent))
    }

    pub(super) fn force_residual(
        &self,
        problem: &HyperProblem<'_>,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SolidError> {
        let mut residual = vec![0.0; self.geometry.ndof];
        self.for_points(
            problem,
            u,
            load,
            cx,
            |index, element, point, material, f, seed| {
                let piola = match self.constraint {
                    Constraint::FiniteBulk => {
                        finite_state(material, &f, self.local, seed).map(|state| state.piola)
                    }
                    Constraint::Incompressible => {
                        exact_state(material, &f, self.local, seed).map(|state| state.piola)
                    }
                }
                .map_err(|error| material_error(index, error))?;
                let weight = point.weight * self.reference_thickness_m;
                for (a, &node) in element.nodes[..element.arity].iter().enumerate() {
                    for i in 0..2 {
                        residual[2 * node + i] += weight
                            * (piola[2 * i] * point.gradients[a][0]
                                + piola[2 * i + 1] * point.gradients[a][1]);
                    }
                }
                Ok(())
            },
        )?;
        problem.for_traction(load, |dof, force| {
            residual[dof] -= self.reference_thickness_m * force;
        })?;
        for (index, &(dof, force)) in self.nodal_forces_n.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            residual[dof] -= force;
        }
        if let Some(inertia) = self.inertia {
            inertia
                .add_force(u, &mut residual, cx)
                .map_err(inertia_error)?;
        }
        if residual.iter().any(|value| !value.is_finite()) {
            return Err(invalid("plane-stress assembled force is nonfinite"));
        }
        checkpoint(cx)?;
        Ok(residual)
    }

    pub(super) fn potential_energy(
        &self,
        problem: &HyperProblem<'_>,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, SolidError> {
        let mut potential = 0.0;
        self.for_points(
            problem,
            u,
            load,
            cx,
            |index, _, point, material, f, seed| {
                let energy_density = match self.constraint {
                    Constraint::FiniteBulk => finite_state(material, &f, self.local, seed)
                        .map(|state| state.energy_density_pa),
                    Constraint::Incompressible => exact_state(material, &f, self.local, seed)
                        .map(|state| state.energy_density_pa),
                }
                .map_err(|error| material_error(index, error))?;
                potential += (point.weight * self.reference_thickness_m) * energy_density;
                Ok(())
            },
        )?;
        problem.for_traction(load, |dof, value| {
            potential -= (self.reference_thickness_m * value) * u[dof];
        })?;
        for (index, &(dof, force)) in self.nodal_forces_n.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            potential -= force * u[dof];
        }
        if let Some(inertia) = self.inertia {
            potential += inertia.potential(u, cx).map_err(inertia_error)?;
        }
        if !potential.is_finite() {
            return Err(invalid("unrepresentable plane-stress potential energy"));
        }
        checkpoint(cx)?;
        Ok(potential)
    }
}

impl HyperProblem<'_> {
    /// One physical endpoint through the existing native Newton owner. Material
    /// and dead tractions use original thickness; M and nodal loads are already
    /// physical. No material history is accepted by this solve.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_plane_stress_dynamic_from(
        &self,
        initial_displacement: &[[f64; 2]],
        local: PlaneStressSettings,
        exact: bool,
        reference_thickness_m: f64,
        inertia: &BackwardEulerInertia<'_>,
        nodal_forces_n: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        let prepared = Prepared::dynamic(
            self,
            local,
            exact,
            reference_thickness_m,
            inertia,
            nodal_forces_n,
            cx,
        )?;
        self.solve_with_projection_from(
            Some(Reduction::PlaneStress(&prepared, cx)),
            Some(initial_displacement),
        )
    }

    /// Full physical residual, retaining constrained-row reactions. Explicit
    /// nodal forces are endpoint values; only dead tractions follow load.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn plane_stress_dynamic_unconstrained_residual(
        &self,
        u: &[f64],
        load: f64,
        local: PlaneStressSettings,
        exact: bool,
        reference_thickness_m: f64,
        inertia: &BackwardEulerInertia<'_>,
        nodal_forces_n: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SolidError> {
        Prepared::dynamic(
            self,
            local,
            exact,
            reference_thickness_m,
            inertia,
            nodal_forces_n,
            cx,
        )?
        .force_residual(self, u, load, cx)
    }

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
        self.solve_plane_stress_constrained(local, Constraint::FiniteBulk, None, cx)
    }

    /// Continue plane-stress equilibrium from an admitted displacement field.
    ///
    /// The same transverse reduction and global Newton owner admit and evolve
    /// the starting geometry. Load ramps and absolute boundary lifting retain
    /// the semantics of [`Self::solve_from`]; no reference or history is reset.
    pub fn solve_plane_stress_from(
        &self,
        initial_displacement: &[[f64; 2]],
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        self.solve_plane_stress_constrained(
            local,
            Constraint::FiniteBulk,
            Some(initial_displacement),
            cx,
        )
    }

    /// Solve an exactly volume-preserving thin hyperelastic continuum.
    ///
    /// The original material remainder supplies energy and transverse shear
    /// relaxation; `det F = J_ref` fixes thickness and traction equilibrium fixes
    /// the pressure multiplier. No bulk penalty or pressure mesh is needed for
    /// this plane-stress reduction. Prescribed natural volume is retained.
    /// This is not a general 3-D mixed element.
    /// Material stretch bounds, global equilibrium and cancellation retain the
    /// same gates as [`Self::solve_plane_stress`].
    pub fn solve_incompressible_plane_stress(
        &self,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        self.solve_plane_stress_constrained(local, Constraint::Incompressible, None, cx)
    }

    /// Continue exact elastic-volume plane stress from an admitted displacement.
    ///
    /// The original reference-volume constraint and material history remain
    /// unchanged. The borrowed starting field is checked before boundary lifting;
    /// load ramps retain the semantics of [`Self::solve_from`].
    pub fn solve_incompressible_plane_stress_from(
        &self,
        initial_displacement: &[[f64; 2]],
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        self.solve_plane_stress_constrained(
            local,
            Constraint::Incompressible,
            Some(initial_displacement),
            cx,
        )
    }

    fn solve_plane_stress_constrained(
        &self,
        local: PlaneStressSettings,
        constraint: Constraint,
        initial_displacement: Option<&[[f64; 2]]>,
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
        let prepared = Prepared::new(self, local, constraint, cx)?;
        self.solve_with_projection_from(
            Some(Reduction::PlaneStress(&prepared, cx)),
            initial_displacement,
        )
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
        Prepared::new(self, local, Constraint::FiniteBulk, cx)?.assemble(
            self,
            u,
            load,
            &BTreeMap::new(),
            cx,
        )
    }

    /// Unconstrained plane-stress forces, including support reactions and dead
    /// tractions. No global tangent is formed; local transverse relaxation
    /// retains its original constitutive derivatives and acceptance gates.
    pub fn plane_stress_unconstrained_residual(
        &self,
        u: &[f64],
        load: f64,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SolidError> {
        Prepared::new(self, local, Constraint::FiniteBulk, cx)?.force_residual(self, u, load, cx)
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
        Prepared::new(self, local, Constraint::FiniteBulk, cx)?.potential_energy(self, u, load, cx)
    }

    /// Recover the actual relaxed three-dimensional states, in element order
    /// and then quadrature order (one exact centroid for P1, four for Q1).
    pub fn plane_stress_states(
        &self,
        u: &[f64],
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Vec<PlaneStressState>, SolidError> {
        let prepared = Prepared::new(self, local, Constraint::FiniteBulk, cx)?;
        let mut states = Vec::new();
        prepared.for_points(self, u, 0.0, cx, |index, _, _, material, f, seed| {
            states.push(
                finite_state(material, &f, local, seed)
                    .map_err(|error| material_error(index, error))?,
            );
            Ok(())
        })?;
        Ok(states)
    }

    /// Unconstrained forces and the Hessian of the exactly incompressible
    /// reduced potential, per unit reference thickness. Support residuals at
    /// equilibrium are reactions, including the recovered pressure contribution.
    pub fn incompressible_plane_stress_residual_and_tangent(
        &self,
        u: &[f64],
        load: f64,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Csr), SolidError> {
        Prepared::new(self, local, Constraint::Incompressible, cx)?.assemble(
            self,
            u,
            load,
            &BTreeMap::new(),
            cx,
        )
    }

    /// Unconstrained exact elastic-volume forces without a global tangent.
    /// Each point retains its original remainder, volume and transverse seed.
    pub fn incompressible_plane_stress_unconstrained_residual(
        &self,
        u: &[f64],
        load: f64,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SolidError> {
        Prepared::new(self, local, Constraint::Incompressible, cx)?
            .force_residual(self, u, load, cx)
    }

    /// Exactly incompressible reduced energy minus original dead-load work.
    /// Neither the energy nor its derivatives use a finite bulk modulus.
    pub fn incompressible_plane_stress_potential_energy(
        &self,
        u: &[f64],
        load: f64,
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<f64, SolidError> {
        Prepared::new(self, local, Constraint::Incompressible, cx)?
            .potential_energy(self, u, load, cx)
    }

    /// Recover volume-preserving 3-D states, physical first Piola stress and
    /// the material remainder's pressure multiplier, in element/quadrature
    /// order. The multiplier's gauge is not the mean Cauchy pressure.
    pub fn incompressible_plane_stress_states(
        &self,
        u: &[f64],
        local: PlaneStressSettings,
        cx: &Cx<'_>,
    ) -> Result<Vec<IncompressiblePlaneStressState>, SolidError> {
        let prepared = Prepared::new(self, local, Constraint::Incompressible, cx)?;
        let mut states = Vec::new();
        prepared.for_points(self, u, 0.0, cx, |index, _, _, material, f, seed| {
            states.push(
                exact_state(material, &f, local, seed)
                    .map_err(|error| material_error(index, error))?,
            );
            Ok(())
        })?;
        Ok(states)
    }
}
