//! Finite-plastic P1/Q1 sheets with traction-free thickness and retained history.
//!
//! The native three-dimensional J2 return map relaxes the transverse column at
//! each material point. Its condensed tangent drives the existing shared
//! increment/Newton solver. Every trial uses the same accepted history and the
//! same original law; a rejected endpoint cannot deform or harden the sheet.
//! Reference area times original thickness supplies every integration weight.
//!
//! This is an in-plane, isothermal continuum with native rate-independent,
//! overstress and creep evolution. Kinetic laws retain their physical interval
//! through every trial; only an accepted endpoint advances the clock. The recovered
//! thickness is a solved stretch, not an imposed total-volume constraint.
//! There is no out-of-plane bending, wrinkling, necking, contact, damage or
//! forming-limit prediction. Work quadrature defects are reported separately
//! from irreversible plastic/creep work and are never added to material heat.

use std::collections::BTreeSet;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::Voigt;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::plastic::finite::plane_stress::{
    FinitePlaneStressMaterial, FinitePlaneStressResponse,
};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Error, FiniteJ2State};
use fs_sparse::{Coo, Csr};

use crate::hyper3::plastic::{FinitePlasticLawField, FinitePlasticStepControls};
use crate::mesh2::{Mesh2, quad_points, shapes_at};
use crate::plastic3::increment::{self, IncrementControls, IncrementOperator};
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

/// Source-bound, frozen-temperature forming with native thermal work accounts.
pub mod thermal;

/// Explicit size and reference-element quality envelope.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SheetAssemblyBudget {
    pub maximum_nodes: usize,
    pub maximum_elements: usize,
    pub maximum_free_dofs: usize,
    /// Minimum positive corner cross product divided by the largest edge².
    pub minimum_scaled_jacobian: f64,
}

/// Native solver, material-point or independently checked energy refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum FinitePlasticSheetError {
    Solid(SmallStrainTetError),
    Material {
        element: usize,
        point: usize,
        source: FiniteJ2Error,
    },
    Energy {
        quantity: &'static str,
        error_j: f64,
        maximum_j: f64,
    },
}

impl From<SmallStrainTetError> for FinitePlasticSheetError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Solid(error)
    }
}
impl core::fmt::Display for FinitePlasticSheetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite plastic sheet refused: {self:?}")
    }
}
impl std::error::Error for FinitePlasticSheetError {}

/// Whole-sheet physical accounts, all in joules.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FinitePlasticSheetEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub plastic_dissipation_j: f64,
    pub creep_dissipation_j: f64,
    /// Reduced incremental potential, distinct from endpoint stored energy.
    pub incremental_potential_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    /// Sum of original reference volume times P_after : (F_after - F_before).
    pub internal_endpoint_work_j: f64,
    /// External minus internal endpoint work, signed and independently gated.
    pub equilibrium_work_residual_j: f64,
    /// Internal endpoint work minus storage change and physical dissipation.
    pub endpoint_work_defect_j: f64,
}

impl FinitePlasticSheetEnergy {
    /// Physical irreversible work; the kinetic potential is not material heat.
    pub fn inelastic_dissipation_j(&self) -> f64 {
        self.plastic_dissipation_j + self.creep_dissipation_j
    }
}

/// One frozen-history force evaluation, with an optional consistent tangent.
#[derive(Debug)]
pub struct FinitePlasticSheetEvaluation {
    pub internal_force_n: Vec<[f64; 2]>,
    pub tangent_n_m: Option<Csr>,
    pub incremental_potential_j: f64,
}

#[derive(Debug)]
struct Point {
    gradients: [[f64; 2]; 4],
    volume_m3: f64,
}

#[derive(Debug)]
struct Element {
    nodes: [usize; 4],
    arity: usize,
    reference: [[f64; 2]; 4],
    points: Vec<Point>,
    offset: usize,
}

impl Element {
    // A Q1 map has an affine in-plane Jacobian. All positive corner values
    // exclude an inverted corner hidden from the four Gauss points.
    fn minimum_corner_cross(&self, u: Option<&[f64]>) -> f64 {
        let mut minimum = f64::INFINITY;
        for a in 0..self.arity {
            let b = (a + 1) % self.arity;
            let c = (a + self.arity - 1) % self.arity;
            let edge = |v: usize, i: usize| {
                self.reference[v][i] - self.reference[a][i]
                    + u.map_or(0.0, |u| u[2 * self.nodes[v] + i] - u[2 * self.nodes[a] + i])
            };
            let cross = edge(b, 0) * edge(c, 1) - edge(b, 1) * edge(c, 0);
            if !cross.is_finite() {
                return f64::NAN;
            }
            minimum = minimum.min(cross);
        }
        minimum
    }

    fn in_plane(&self, point: &Point, u: &[f64]) -> [f64; 4] {
        let mut f = [1.0, 0.0, 0.0, 1.0];
        // Remove translation before contraction with reference gradients.
        for a in 1..self.arity {
            for i in 0..2 {
                let du = u[2 * self.nodes[a] + i] - u[2 * self.nodes[0] + i];
                for j in 0..2 {
                    f[2 * i + j] = du.mul_add(point.gradients[a][j], f[2 * i + j]);
                }
            }
        }
        f
    }
}

/// Immutable accepted sheet. Geometry and original native/source-resolved laws
/// are borrowed; prepared geometry is shared across accepted endpoints.
#[derive(Debug)]
pub struct FinitePlasticSheet<'a, L: FinitePlaneStressMaterial = FiniteJ2> {
    mesh: &'a Mesh2,
    laws: FinitePlasticLawField<'a, L>,
    thickness_m: &'a [f64],
    elements: Arc<Vec<Element>>,
    local: PlaneStressSettings,
    budget: SheetAssemblyBudget,
    displacement_m: Vec<[f64; 2]>,
    states: Vec<FiniteJ2State>,
    deformation: Vec<[f64; 9]>,
    piola_pa: Vec<[f64; 9]>,
    stored_energy_j: f64,
    physical_time_s: Option<f64>,
}

/// A complete accepted endpoint; the previous borrowed sheet is unchanged.
#[derive(Debug)]
pub struct FinitePlasticSheetStep<'a, L: FinitePlaneStressMaterial = FiniteJ2> {
    pub sheet: FinitePlasticSheet<'a, L>,
    pub equilibrium: IncrementRecord,
    pub energy: FinitePlasticSheetEnergy,
}

fn invalid(what: impl Into<String>) -> SmallStrainTetError {
    SmallStrainTetError::InvalidInput { what: what.into() }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), SmallStrainTetError> {
    cx.checkpoint().map_err(|_| SmallStrainTetError::Cancelled)
}
fn energy_gate(
    quantity: &'static str,
    error_j: f64,
    maximum_j: f64,
) -> Result<(), FinitePlasticSheetError> {
    if !error_j.is_finite() || error_j.abs() > maximum_j {
        return Err(FinitePlasticSheetError::Energy {
            quantity,
            error_j,
            maximum_j,
        });
    }
    Ok(())
}

fn shared_duration_s<L: FinitePlaneStressMaterial>(
    laws: FinitePlasticLawField<'_, L>,
    cx: &Cx<'_>,
) -> Result<Option<f64>, SmallStrainTetError> {
    let mut duration = None;
    let mut inspect = |law: &L| {
        checkpoint(cx)?;
        if let Some(dt) = law.duration_s() {
            if !dt.is_finite() || dt <= 0.0 {
                return Err(invalid(
                    "sheet material duration must be finite and positive",
                ));
            }
            if duration.is_some_and(|accepted| accepted != dt) {
                return Err(invalid(
                    "sheet kinetic laws must share one physical duration",
                ));
            }
            duration = Some(dt);
        }
        Ok(())
    };
    match laws {
        FinitePlasticLawField::Uniform(law) => inspect(law)?,
        FinitePlasticLawField::PerElement(laws) => {
            for law in laws {
                inspect(law)?;
            }
        }
    }
    Ok(duration)
}

fn prepare_elements(
    mesh: &Mesh2,
    thickness_m: &[f64],
    budget: SheetAssemblyBudget,
    cx: &Cx<'_>,
) -> Result<Arc<Vec<Element>>, SmallStrainTetError> {
    checkpoint(cx)?;
    if mesh.nodes.is_empty()
        || mesh.elems.is_empty()
        || mesh.nodes.len() > budget.maximum_nodes
        || mesh.elems.len() > budget.maximum_elements
        || mesh.nodes.len().checked_mul(2).is_none()
        || mesh.elems.len().checked_mul(4).is_none()
        || thickness_m.len() != mesh.elems.len()
        || thickness_m.iter().any(|x| !x.is_finite() || *x <= 0.0)
        || mesh.nodes.iter().flatten().any(|x| !x.is_finite())
        || !budget.minimum_scaled_jacobian.is_finite()
        || budget.minimum_scaled_jacobian <= 0.0
    {
        return Err(invalid(
            "finite sheet mesh, exact positive thicknesses and explicit size/quality bounds required",
        ));
    }
    let mut elements = Vec::with_capacity(mesh.elems.len());
    let mut reference_volume_m3 = 0.0;
    let mut point_count = 0;
    let mut used = vec![false; mesh.nodes.len()];
    let mut seen = BTreeSet::new();
    for (e, conn) in mesh.elems.iter().enumerate() {
        checkpoint(cx)?;
        let arity = conn.len();
        if !(arity == 3 || arity == 4)
            || conn.iter().any(|&a| a >= mesh.nodes.len())
            || (0..arity).any(|a| (0..a).any(|b| conn[a] == conn[b]))
        {
            return Err(invalid(format!(
                "sheet element {e} must be a P1 triangle or Q1 quad with distinct in-range nodes"
            )));
        }
        let mut key = conn.clone();
        key.sort_unstable();
        if !seen.insert(key) {
            return Err(invalid("duplicate sheet element"));
        }
        let mut element = Element {
            nodes: [0; 4],
            arity,
            reference: [[0.0; 2]; 4],
            points: Vec::new(),
            offset: point_count,
        };
        element.nodes[..arity].copy_from_slice(conn);
        for (a, &node) in conn.iter().enumerate() {
            used[node] = true;
            element.reference[a] = [0, 1].map(|i| mesh.nodes[node][i] - mesh.nodes[conn[0]][i]);
        }
        let mut edge_squared = 0.0_f64;
        for a in 0..arity {
            let b = (a + 1) % arity;
            let d = [0, 1].map(|i| element.reference[b][i] - element.reference[a][i]);
            edge_squared = edge_squared.max(d[0] * d[0] + d[1] * d[1]);
        }
        let scaled = element.minimum_corner_cross(None) / edge_squared;
        if !scaled.is_finite() || scaled < budget.minimum_scaled_jacobian {
            return Err(invalid(format!(
                "sheet element {e} has inadmissible corner scaled Jacobian {scaled:e}"
            )));
        }
        // Constant P1 F is exactly integrated at one centroid. Q1 retains
        // independent histories at all four native Gauss points.
        let rule = if arity == 3 {
            vec![(1.0 / 3.0, 1.0 / 3.0, 0.5)]
        } else {
            quad_points(4)
        };
        for (xi, eta, weight) in rule {
            let (_, gradients, determinant) =
                shapes_at(&element.reference, &[0, 1, 2, 3][..arity], xi, eta);
            let volume_m3 = weight * determinant * thickness_m[e];
            if !volume_m3.is_finite()
                || volume_m3 <= 0.0
                || gradients.iter().flatten().any(|x| !x.is_finite())
            {
                return Err(invalid(format!(
                    "sheet element {e} has unrepresentable reference geometry"
                )));
            }
            let mut point = Point {
                gradients: [[0.0; 2]; 4],
                volume_m3,
            };
            point.gradients[..arity].copy_from_slice(&gradients);
            element.points.push(point);
            reference_volume_m3 += volume_m3;
            point_count += 1;
        }
        elements.push(element);
    }
    if used.iter().any(|used| !used) || !reference_volume_m3.is_finite() {
        return Err(invalid(
            "sheet contains an unused node or unrepresentable volume",
        ));
    }
    checkpoint(cx)?;
    Ok(Arc::new(elements))
}

fn validate_step_inputs(
    dofs: usize,
    budget: SheetAssemblyBudget,
    prescribed_m: &[(usize, f64)],
    nodal_forces_n: &[(usize, f64)],
    settings: IncrementSettings,
    controls: FinitePlasticStepControls,
) -> Result<Vec<f64>, SmallStrainTetError> {
    if settings.max_newton_iterations == 0
        || settings.max_linear_iterations == 0
        || !settings.relative_tolerance.is_finite()
        || settings.relative_tolerance <= 0.0
        || settings.relative_tolerance >= 1.0
        || !settings.linear_tolerance.is_finite()
        || settings.linear_tolerance <= 0.0
        || settings.linear_tolerance >= 1.0
        || [
            controls.equilibrium_work_tolerance_j,
            controls.maximum_endpoint_work_defect_j,
            controls.yield_tolerance_pa,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(invalid(
            "sheet needs finite positive solver controls and nonnegative energy/yield gates",
        ));
    }
    let mut fixed = vec![false; dofs];
    for &(dof, value) in prescribed_m {
        if dof >= dofs || !value.is_finite() || fixed[dof] {
            return Err(invalid("invalid or repeated sheet prescribed DOF"));
        }
        fixed[dof] = true;
    }
    if dofs - prescribed_m.len() > budget.maximum_free_dofs {
        return Err(invalid("sheet free displacement count exceeds budget"));
    }
    let mut external = vec![0.0; dofs];
    for &(dof, value) in nodal_forces_n {
        if dof >= dofs || !value.is_finite() {
            return Err(invalid("invalid sheet dead nodal force"));
        }
        external[dof] += value;
    }
    if external.iter().any(|x| !x.is_finite()) {
        return Err(invalid("sheet summed force is not representable"));
    }
    Ok(external)
}

impl<'a, L: FinitePlaneStressMaterial> FinitePlasticSheet<'a, L> {
    /// Admit a virgin P1/Q1 reference sheet. Thickness has one positive value
    /// per element; a graded sheet retains those original volume weights.
    /// Kinetic owners must share one duration. Construction cannot consume a
    /// physical interval or hide nonzero inelastic initialization flow.
    pub fn try_new(
        mesh: &'a Mesh2,
        laws: FinitePlasticLawField<'a, L>,
        thickness_m: &'a [f64],
        local: PlaneStressSettings,
        budget: SheetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FinitePlasticSheetError> {
        let elements = prepare_elements(mesh, thickness_m, budget, cx)?;
        if let FinitePlasticLawField::PerElement(materials) = laws
            && materials.len() != mesh.elems.len()
        {
            return Err(SmallStrainTetError::LawCountMismatch {
                expected: mesh.elems.len(),
                got: materials.len(),
            }
            .into());
        }
        let physical_time_s = shared_duration_s(laws, cx)?.map(|_| 0.0);
        let mut states = Vec::new();
        let mut deformation = Vec::new();
        let mut piola_pa = Vec::new();
        let mut stored_energy_j = 0.0;
        for (e, element) in elements.iter().enumerate() {
            checkpoint(cx)?;
            let material = match laws {
                FinitePlasticLawField::Uniform(law) => law,
                FinitePlasticLawField::PerElement(laws) => &laws[e],
            };
            let virgin = material.initial_state();
            let initial = material
                .plane_stress_response(&[1.0, 0.0, 0.0, 1.0], &virgin, [0.0, 0.0, 1.0], local)
                .map_err(|source| FinitePlasticSheetError::Material {
                    element: e,
                    point: 0,
                    source,
                })?;
            if initial.material.state != virgin
                || initial.material.plastic_increment != 0.0
                || initial.material.creep_increment != 0.0
                || initial.material.plastic_dissipation_j_m3 != 0.0
                || initial.material.creep_dissipation_j_m3 != 0.0
            {
                return Err(invalid(format!(
                    "sheet element {e} initialization requires inelastic flow"
                ))
                .into());
            }
            for point in &element.points {
                states.push(initial.material.state.clone());
                deformation.push(initial.deformation_gradient);
                piola_pa.push(initial.material.first_piola_pa);
                stored_energy_j += point.volume_m3 * initial.material.stored_energy_j_m3;
            }
        }
        if !stored_energy_j.is_finite() {
            return Err(invalid("sheet initial energy is not representable").into());
        }
        checkpoint(cx)?;
        Ok(Self {
            mesh,
            laws,
            thickness_m,
            elements,
            local,
            budget,
            displacement_m: vec![[0.0; 2]; mesh.nodes.len()],
            states,
            deformation,
            piola_pa,
            stored_energy_j,
            physical_time_s,
        })
    }

    /// Exact original law binding, retained after every accepted endpoint.
    pub fn laws(&self) -> FinitePlasticLawField<'a, L> {
        self.laws
    }
    pub fn displacement_m(&self) -> &[[f64; 2]] {
        &self.displacement_m
    }
    pub fn quadrature_states(&self) -> &[FiniteJ2State] {
        &self.states
    }
    pub fn quadrature_deformation(&self) -> &[[f64; 9]] {
        &self.deformation
    }
    pub fn quadrature_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.piola_pa
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.stored_energy_j
    }
    /// Accepted physical time for kinetic sheets; static endpoints have no clock.
    pub fn physical_time_s(&self) -> Option<f64> {
        self.physical_time_s
    }
    pub fn reference_volume_m3(&self, element: usize) -> Option<f64> {
        self.elements
            .get(element)
            .map(|element| element.points.iter().map(|p| p.volume_m3).sum())
    }
    /// Positive current normal thickness in original element/quadrature order.
    pub fn quadrature_thickness_m(&self) -> Vec<f64> {
        self.elements
            .iter()
            .enumerate()
            .flat_map(|(e, element)| {
                (0..element.points.len())
                    .map(move |q| self.thickness_m[e] * self.deformation[element.offset + q][8])
            })
            .collect()
    }

    fn law(&self, element: usize) -> &L {
        match self.laws {
            FinitePlasticLawField::Uniform(law) => law,
            FinitePlasticLawField::PerElement(laws) => &laws[element],
        }
    }

    fn response(
        &self,
        e: usize,
        q: usize,
        u: &[f64],
        states: &[FiniteJ2State],
    ) -> Result<FinitePlaneStressResponse, FiniteJ2Error> {
        let element = &self.elements[e];
        let index = element.offset + q;
        self.law(e).plane_stress_response(
            &element.in_plane(&element.points[q], u),
            &states[index],
            [2, 5, 8].map(|i| self.deformation[index][i]),
            self.local,
        )
    }

    /// Evaluate an arbitrary trial from this accepted history without advancing
    /// it. The tangent includes the full local traction-free condensation.
    pub fn evaluate_trial(
        &self,
        displacement_m: &[[f64; 2]],
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<FinitePlasticSheetEvaluation, FinitePlasticSheetError> {
        let result = SheetOperator { sheet: self }
            .assemble(
                &self.states,
                displacement_m.as_flattened(),
                None,
                need_tangent,
                cx,
            )?
            .0;
        Ok(FinitePlasticSheetEvaluation {
            internal_force_n: result.force.chunks_exact(2).map(|u| [u[0], u[1]]).collect(),
            tangent_n_m: result.tangent,
            incremental_potential_j: result.potential,
        })
    }

    /// Solve one physical endpoint from the accepted displacement and complete
    /// point histories. DOFs are `2 * node + component`; forces are dead nodal
    /// forces in N and prescribed values are absolute displacements in m.
    /// Native kinetic laws supply one shared duration for this whole endpoint.
    /// Every trial retains that same duration and the original accepted history.
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePlasticSheetStep<'a, L>, FinitePlasticSheetError> {
        checkpoint(cx)?;
        let dofs = 2 * self.mesh.nodes.len();
        let external = validate_step_inputs(
            dofs,
            self.budget,
            prescribed_m,
            nodal_forces_n,
            settings,
            controls,
        )?;
        let physical_time_s = match (self.physical_time_s, shared_duration_s(self.laws, cx)?) {
            (None, None) => None,
            (Some(time), Some(dt)) => {
                let next = time + dt;
                if !next.is_finite() || next <= time {
                    return Err(invalid("sheet physical time cannot advance representably").into());
                }
                Some(next)
            }
            _ => return Err(invalid("sheet material clock changed after admission").into()),
        };
        let mut old_storage = 0.0;
        for (e, element) in self.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = element.offset + q;
                let material = self.law(e);
                let map = |source| FinitePlasticSheetError::Material {
                    element: e,
                    point: q,
                    source,
                };
                material
                    .validate_committed(
                        &self.deformation[index],
                        &self.states[index],
                        controls.yield_tolerance_pa,
                    )
                    .map_err(map)?;
                old_storage += point.volume_m3
                    * material
                        .stored_energy(&self.deformation[index], &self.states[index])
                        .map_err(map)?;
            }
        }
        energy_gate(
            "accepted storage",
            old_storage - self.stored_energy_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        let operator = SheetOperator { sheet: self };
        let solved = increment::solve(
            &operator,
            IncrementControls {
                dofs,
                prescribed_m,
                nodal_forces_n,
                load_path: &[1.0],
                settings,
            },
            Some((self.displacement_m.as_flattened(), &self.states)),
            None,
            cx,
        )?;
        let equilibrium = solved
            .increments
            .into_iter()
            .next()
            .expect("one physical endpoint");
        let assembled = operator
            .assemble(&self.states, &solved.displacement_m, None, false, cx)?
            .0;
        let mut energy = FinitePlasticSheetEnergy {
            stored_before_j: old_storage,
            ..Default::default()
        };
        let mut deformation = Vec::with_capacity(self.states.len());
        let mut piola_pa = Vec::with_capacity(self.states.len());
        for (e, element) in self.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = element.offset + q;
                let map = |source| FinitePlasticSheetError::Material {
                    element: e,
                    point: q,
                    source,
                };
                // Audit from the ORIGINAL history, never perform a second return
                // from the newly accepted state when computing plastic work.
                let response = self
                    .response(e, q, &solved.displacement_m, &self.states)
                    .map_err(map)?;
                if response.material.state != solved.point_states[index] {
                    return Err(invalid(
                        "sheet endpoint audit differs from accepted frozen-history return",
                    )
                    .into());
                }
                self.law(e)
                    .validate_committed(
                        &response.deformation_gradient,
                        &response.material.state,
                        controls.yield_tolerance_pa,
                    )
                    .map_err(map)?;
                let stored = self
                    .law(e)
                    .stored_energy(&response.deformation_gradient, &response.material.state)
                    .map_err(map)?;
                energy_gate(
                    "point storage",
                    point.volume_m3 * (stored - response.material.stored_energy_j_m3),
                    controls.equilibrium_work_tolerance_j,
                )?;
                if [
                    response.material.plastic_dissipation_j_m3,
                    response.material.creep_dissipation_j_m3,
                ]
                .iter()
                .any(|x| !x.is_finite() || *x < 0.0)
                {
                    return Err(invalid(format!(
                        "sheet element {e} point {q} has invalid physical dissipation"
                    ))
                    .into());
                }
                energy.stored_after_j += point.volume_m3 * stored;
                energy.plastic_dissipation_j +=
                    point.volume_m3 * response.material.plastic_dissipation_j_m3;
                energy.creep_dissipation_j +=
                    point.volume_m3 * response.material.creep_dissipation_j_m3;
                energy.incremental_potential_j +=
                    point.volume_m3 * response.material.incremental_potential_j_m3;
                energy.internal_endpoint_work_j += point.volume_m3
                    * (0..9)
                        .map(|i| {
                            response.material.first_piola_pa[i]
                                * (response.deformation_gradient[i] - self.deformation[index][i])
                        })
                        .sum::<f64>();
                deformation.push(response.deformation_gradient);
                piola_pa.push(response.material.first_piola_pa);
            }
        }
        let old_u = self.displacement_m.as_flattened();
        energy.dead_load_endpoint_work_j = external
            .iter()
            .enumerate()
            .map(|(dof, force)| force * (solved.displacement_m[dof] - old_u[dof]))
            .sum();
        energy.support_endpoint_work_j = prescribed_m
            .iter()
            .zip(&equilibrium.reactions_n)
            .map(|(&(dof, _), reaction)| reaction * (solved.displacement_m[dof] - old_u[dof]))
            .sum();
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j
            - (energy.stored_after_j - energy.stored_before_j)
            - energy.inelastic_dissipation_j();
        energy_gate(
            "equilibrium endpoint work",
            energy.equilibrium_work_residual_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "endpoint work quadrature defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        energy_gate(
            "assembled incremental potential",
            assembled.potential - energy.incremental_potential_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        checkpoint(cx)?;
        Ok(FinitePlasticSheetStep {
            sheet: Self {
                mesh: self.mesh,
                laws: self.laws,
                thickness_m: self.thickness_m,
                elements: Arc::clone(&self.elements),
                local: self.local,
                budget: self.budget,
                displacement_m: solved
                    .displacement_m
                    .chunks_exact(2)
                    .map(|u| [u[0], u[1]])
                    .collect(),
                states: solved.point_states,
                deformation,
                piola_pa,
                stored_energy_j: energy.stored_after_j,
                physical_time_s,
            },
            equilibrium,
            energy,
        })
    }
}

// Both fixed-temperature and thermal owners use this same geometry, assembly
// and Newton adapter. The thermal owner retains its native error and work audit.
trait SheetPointResponse {
    type Error: core::fmt::Display;
    fn node_count(&self) -> usize;
    fn elements(&self) -> &[Element];
    fn states(&self) -> &[FiniteJ2State];
    fn response(
        &self,
        element: usize,
        point: usize,
        displacement: &[f64],
        states: &[FiniteJ2State],
    ) -> Result<FinitePlaneStressResponse, Self::Error>;
}

impl<L: FinitePlaneStressMaterial> SheetPointResponse for FinitePlasticSheet<'_, L> {
    type Error = FiniteJ2Error;
    fn node_count(&self) -> usize {
        self.mesh.nodes.len()
    }
    fn elements(&self) -> &[Element] {
        &self.elements
    }
    fn states(&self) -> &[FiniteJ2State] {
        &self.states
    }
    fn response(
        &self,
        element: usize,
        point: usize,
        displacement: &[f64],
        states: &[FiniteJ2State],
    ) -> Result<FinitePlaneStressResponse, Self::Error> {
        FinitePlasticSheet::response(self, element, point, displacement, states)
    }
}

struct SheetOperator<'s, S: SheetPointResponse> {
    sheet: &'s S,
}
struct Assembled {
    force: Vec<f64>,
    tangent: Option<Csr>,
    potential: f64,
}

impl<S: SheetPointResponse> SheetOperator<'_, S> {
    fn assemble(
        &self,
        states: &[FiniteJ2State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        need_tangent: bool,
        cx: &Cx<'_>,
    ) -> Result<(Assembled, f64), SmallStrainTetError> {
        checkpoint(cx)?;
        let sheet = self.sheet;
        let dofs = 2 * sheet.node_count();
        if u.len() != dofs
            || states.len() != sheet.states().len()
            || u.iter().any(|x| !x.is_finite())
        {
            return Err(invalid(
                "sheet trial needs complete finite displacements and point history",
            ));
        }
        let mut force = vec![0.0; dofs];
        let mut matrix = need_tangent.then(|| Coo::new(dofs, dofs));
        let mut potential = 0.0;
        let mut force_cap = 0.0;
        for (e, element) in sheet.elements().iter().enumerate() {
            let orientation = element.minimum_corner_cross(Some(u));
            if !orientation.is_finite() || orientation <= 0.0 {
                return Err(invalid(format!(
                    "sheet trial element {e} is inverted or folded"
                )));
            }
            let mut local_force = [0.0; 8];
            let mut local_tangent = [[0.0; 8]; 8];
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let response = sheet
                    .response(e, q, u, states)
                    .map_err(|error| invalid(format!("sheet element {e} point {q}: {error}")))?;
                potential += point.volume_m3 * response.material.incremental_potential_j_m3;
                for a in 0..element.arity {
                    for i in 0..2 {
                        let row = 2 * a + i;
                        local_force[row] += point.volume_m3
                            * (0..2)
                                .map(|j| response.piola_pa[2 * i + j] * point.gradients[a][j])
                                .sum::<f64>();
                        if need_tangent {
                            for b in 0..element.arity {
                                for k in 0..2 {
                                    let mut value = 0.0;
                                    for j in 0..2 {
                                        for l in 0..2 {
                                            value += point.gradients[a][j]
                                                * response.tangent_pa[2 * i + j][2 * k + l]
                                                * point.gradients[b][l];
                                        }
                                    }
                                    local_tangent[row][2 * b + k] += point.volume_m3 * value;
                                }
                            }
                        }
                    }
                }
            }
            for i in 0..2 * element.arity {
                let row = 2 * element.nodes[i / 2] + i % 2;
                force[row] += local_force[i];
                if prescribed.is_none_or(|fixed| !fixed[row]) {
                    force_cap += local_force[i].abs();
                }
                if let Some(matrix) = &mut matrix {
                    for j in 0..2 * element.arity {
                        let value = local_tangent[i][j];
                        if value != 0.0 {
                            matrix.push(row, 2 * element.nodes[j / 2] + j % 2, value);
                        }
                    }
                }
            }
        }
        let tangent = matrix.map(|m| m.assemble());
        if !potential.is_finite()
            || !force_cap.is_finite()
            || force.iter().any(|x| !x.is_finite())
            || tangent
                .as_ref()
                .is_some_and(|m| (0..m.nrows()).any(|r| m.row(r).1.iter().any(|x| !x.is_finite())))
        {
            return Err(invalid(
                "sheet assembled force, tangent or potential is not representable",
            ));
        }
        checkpoint(cx)?;
        Ok((
            Assembled {
                force,
                tangent,
                potential,
            },
            force_cap,
        ))
    }
}

impl<S: SheetPointResponse> IncrementOperator for SheetOperator<'_, S> {
    type State = FiniteJ2State;
    fn initial_states(&self) -> Vec<Self::State> {
        self.sheet.states().to_vec()
    }
    fn components_per_node(&self) -> usize {
        2
    }
    fn internal(
        &self,
        states: &[Self::State],
        u: &[f64],
        prescribed: Option<&[bool]>,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), SmallStrainTetError> {
        let (result, cap) = self.assemble(states, u, prescribed, false, cx)?;
        Ok((result.force, cap))
    }
    fn tangent(
        &self,
        states: &[Self::State],
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Csr, SmallStrainTetError> {
        Ok(self
            .assemble(states, u, None, true, cx)?
            .0
            .tangent
            .expect("tangent requested"))
    }
    fn tangent_symmetric(&self) -> bool {
        true
    }
    fn commit(
        &self,
        states: &mut [Self::State],
        u: &[f64],
        _final_stress: Option<&mut Vec<Voigt>>,
        cx: &Cx<'_>,
    ) -> Result<usize, SmallStrainTetError> {
        let mut changed = 0;
        for (e, element) in self.sheet.elements().iter().enumerate() {
            let mut evolving = false;
            for q in 0..element.points.len() {
                checkpoint(cx)?;
                let index = element.offset + q;
                let response = self
                    .sheet
                    .response(e, q, u, states)
                    .map_err(|error| invalid(format!("sheet element {e} point {q}: {error}")))?;
                evolving |= response.material.state != states[index];
                states[index] = response.material.state;
            }
            changed += usize::from(evolving);
        }
        Ok(changed)
    }
    fn recoverable_trial_error(&self, error: &SmallStrainTetError) -> bool {
        matches!(error, SmallStrainTetError::InvalidInput { .. })
    }
}
