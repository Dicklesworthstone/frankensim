//! Reference-mass enthalpy transport with explicit material assignments.
//!
//! Each tetrahedron names a declared phase chart and frozen reference density.
//! A shared vertex must have exactly one material identity. Different materials
//! therefore meet through separate traces and an explicitly declared thermal
//! contact, not through an inferred mixture of enthalpy charts. Conductivity
//! assignments remain owned by `ConductionProblem::element_materials`.
//!
//! This is a stationary reference-mass model: equilibrium density does not
//! change storage, geometry does not move, and no phase-chart blending,
//! material-selection derivative or Dirichlet enthalpy is inferred.

use std::fmt;

use fs_exec::Cx;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;

use super::adjoint::{EnthalpyAdjointError, EnthalpyStepGradient, EnthalpyStepLinearization};
use super::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
    finite, poll,
};
use crate::{
    ConductionMesh, ConductionProblem, LinearConfig, ThermalInterfaces, assemble::ASSEMBLY_TILE,
};

/// One immutable constitutive chart and its separately declared storage density.
#[derive(Debug, Clone, Copy)]
pub struct ReferenceEnthalpyMaterial<'c> {
    /// Validated equilibrium enthalpy/temperature/phase chart.
    pub curve: &'c EquilibriumEnthalpyPhaseCurve,
    /// Frozen reference mass per mesh volume, kg/m3. Not equilibrium density.
    pub reference_density_kg_m3: f64,
}

/// A refused heterogeneous assignment publishes no prepared storage model.
#[derive(Debug, Clone, PartialEq)]
pub enum HeterogeneousEnthalpyError {
    /// Mesh, density, work, arithmetic or cancellation admission failed.
    Enthalpy(EnthalpyError),
    /// A tetrahedron refers outside the supplied material table.
    UnknownMaterial {
        /// Offending tetrahedron.
        element: usize,
        /// Supplied material-table index.
        material: usize,
        /// Number of available material records.
        material_count: usize,
    },
    /// Two material identities meet at one enthalpy degree of freedom.
    SharedVertex {
        /// Conflicting mesh vertex.
        vertex: usize,
        /// Material already assigned by an earlier tetrahedron.
        first_material: usize,
        /// Conflicting material assigned by the current tetrahedron.
        next_material: usize,
    },
}

impl fmt::Display for HeterogeneousEnthalpyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Enthalpy(error) => error.fmt(f),
            Self::UnknownMaterial {
                element,
                material,
                material_count,
            } => write!(
                f,
                "enthalpy element {element} selects material {material}, but only \
                 {material_count} materials were supplied"
            ),
            Self::SharedVertex {
                vertex,
                first_material,
                next_material,
            } => write!(
                f,
                "enthalpy vertex {vertex} belongs to materials {first_material} and \
                 {next_material}; use distinct interface vertices and an explicit \
                 thermal contact instead of an undeclared phase-chart mixture"
            ),
        }
    }
}

impl std::error::Error for HeterogeneousEnthalpyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Enthalpy(error) => Some(error),
            Self::UnknownMaterial { .. } | Self::SharedVertex { .. } => None,
        }
    }
}

impl From<EnthalpyError> for HeterogeneousEnthalpyError {
    fn from(error: EnthalpyError) -> Self {
        Self::Enthalpy(error)
    }
}

/// Prepared heterogeneous masses and nodal phase-chart ownership.
///
/// The uniform stepper is private, so its single-chart accessor cannot be used
/// to misrepresent this model. All public chart inspection is vertex-specific.
#[derive(Debug)]
pub struct HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    inner: EnthalpyBackwardEuler<'m, 'c>,
    element_material_ids: Vec<usize>,
    vertex_material_ids: Vec<usize>,
    material_count: usize,
}

impl<'m, 'c> HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    /// Prepare `m_i = sum_(e containing i) rho_e * V_e / 4` and chart ownership.
    ///
    /// `element_material_ids` is in mesh tetrahedron order. Material indices
    /// are identities, even when two records happen to reference equal curves.
    /// The existing element budget also caps the material table, including
    /// unused entries, so validation and later material pullbacks are bounded.
    /// Density validation, shape checks and checked size products precede the
    /// spatial allocations. Cancellation is polled throughout preparation.
    ///
    /// # Errors
    /// Refuses invalid densities, exhausted budgets, missing assignments,
    /// unknown material IDs, orphan vertices and conflicting shared vertices.
    #[allow(clippy::too_many_lines)]
    pub fn new(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        materials: &[ReferenceEnthalpyMaterial<'c>],
        element_material_ids: &[usize],
        budget: EnthalpyBudget,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0)?;
        let n = mesh.vertex_count();
        let elements = mesh.element_count();
        for (resource, required, limit) in [
            ("vertices", n, budget.max_vertices),
            ("elements", elements, budget.max_elements),
            ("materials", materials.len(), budget.max_elements),
        ] {
            if required > limit {
                return Err(EnthalpyError::Budget {
                    resource,
                    required,
                    limit,
                }
                .into());
            }
        }
        if n == 0 || elements == 0 || materials.is_empty() {
            return Err(EnthalpyError::InvalidInput(
                "nonempty mesh and reference-material table required",
            )
            .into());
        }
        if element_material_ids.len() != elements {
            return Err(EnthalpyError::InvalidInput(
                "one reference-material assignment per tetrahedron required",
            )
            .into());
        }
        n.checked_mul(n)
            .and_then(|v| v.checked_mul(24))
            .and_then(|_| elements.checked_mul(32 * 24))
            .ok_or(EnthalpyError::InvalidInput(
                "spatial allocation size overflow",
            ))?;
        for (index, material) in materials.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            if !material.reference_density_kg_m3.is_finite()
                || material.reference_density_kg_m3 <= 0.0
            {
                return Err(EnthalpyError::InvalidInput(
                    "every material requires positive finite reference density",
                )
                .into());
            }
        }
        let mut masses = vec![0.0; n];
        let mut nodal_curves = vec![materials[0].curve; n];
        let mut vertex_material_ids = vec![usize::MAX; n];
        for (element, vertices) in mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let id = element_material_ids[element];
            let material =
                materials
                    .get(id)
                    .ok_or(HeterogeneousEnthalpyError::UnknownMaterial {
                        element,
                        material: id,
                        material_count: materials.len(),
                    })?;
            let contribution =
                finite(material.reference_density_kg_m3 * mesh.element_volume(element) / 4.0)?;
            if contribution <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "positive reference mass is not representable",
                )
                .into());
            }
            for &vertex in vertices {
                let vertex = vertex as usize;
                let previous = vertex_material_ids[vertex];
                if previous != usize::MAX && previous != id {
                    return Err(HeterogeneousEnthalpyError::SharedVertex {
                        vertex,
                        first_material: previous,
                        next_material: id,
                    });
                }
                vertex_material_ids[vertex] = id;
                nodal_curves[vertex] = material.curve;
                masses[vertex] = finite(masses[vertex] + contribution)?;
            }
        }
        for (vertex, &mass) in masses.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            if mass <= 0.0 {
                return Err(EnthalpyError::InvalidInput(
                    "every vertex requires positive reference mass",
                )
                .into());
            }
        }
        poll(cx, n)?;
        Ok(Self {
            inner: EnthalpyBackwardEuler {
                mesh,
                curve: materials[0].curve,
                nodal_curves: Some(nodal_curves),
                masses,
                phase_conductivity: None,
            },
            element_material_ids: element_material_ids.to_vec(),
            vertex_material_ids,
            material_count: materials.len(),
        })
    }

    /// Bind one explicit liquid-mass-fraction conductivity multiplier per
    /// tetrahedron. Chart and conductivity material assignments retain their
    /// separate identities; this scales each element's existing base tensor.
    pub fn with_phase_conductivity(
        mut self,
        cx: &Cx<'_>,
        element_laws: Vec<crate::material::LiquidMassFractionConductivity>,
    ) -> Result<Self, EnthalpyError> {
        self.inner = self.inner.with_phase_conductivity(cx, element_laws)?;
        Ok(self)
    }

    /// Invariant reference masses in nodal mesh order, kg.
    #[must_use]
    pub fn reference_nodal_masses_kg(&self) -> &[f64] {
        self.inner.reference_nodal_masses_kg()
    }

    /// Exact chart assigned to a mesh vertex; out-of-range vertices return None.
    #[must_use]
    pub fn phase_curve_at(&self, vertex: usize) -> Option<&EquilibriumEnthalpyPhaseCurve> {
        self.vertex_material_ids
            .get(vertex)
            .map(|_| self.inner.curve_for_vertex(vertex))
    }

    /// Reference-material identities in mesh tetrahedron order.
    #[must_use]
    pub fn element_material_ids(&self) -> &[usize] {
        &self.element_material_ids
    }

    /// Reference-material identities in nodal mesh order.
    #[must_use]
    pub fn vertex_material_ids(&self) -> &[usize] {
        &self.vertex_material_ids
    }

    /// Number of admitted material records, including unused records.
    #[must_use]
    pub const fn material_count(&self) -> usize {
        self.material_count
    }

    /// Advance using the shared nonlinear transport, contact and energy gates.
    /// History is immutable; phase fields use each vertex's assigned chart.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        self.inner
            .advance(cx, problem, interfaces, old_specific_h, dt_s, config)
    }

    /// Add implicit ambient radiation while retaining each vertex's declared
    /// phase chart and immutable reference mass. The shared driver checks the
    /// complete endpoint residual and nonlinear radiative energy balance.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        patches: &[crate::AmbientRadiationPatch],
        radiation_config: crate::AmbientRadiationConfig,
    ) -> Result<super::radiation::EnthalpyRadiationStepSolution, EnthalpyError> {
        self.inner.advance_with_ambient_radiation(
            cx,
            problem,
            interfaces,
            old_specific_h,
            dt_s,
            config,
            patches,
            radiation_config,
        )
    }

    /// Bind total ambient-radiation derivatives while retaining each vertex's
    /// phase chart and reference mass. Reference densities remain fixed here;
    /// the returned controls are history, source, dt, reservoir and emissivity.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_step_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        patches: &[crate::AmbientRadiationPatch],
        radiation_config: crate::AmbientRadiationConfig,
        max_feedback_entries: usize,
    ) -> Result<super::EnthalpyRadiationStepLinearization<'m>, EnthalpyAdjointError> {
        self.inner.linearize_step_with_ambient_radiation(
            cx,
            problem,
            interfaces,
            old_specific_h,
            dt_s,
            config,
            patches,
            radiation_config,
            max_feedback_entries,
        )
    }

    /// Recheck the physical radiative endpoint before binding its derivatives.
    /// Original convection data and heterogeneous reference masses are fixed.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        patches: &[crate::AmbientRadiationPatch],
        accepted: EnthalpyStepSolution,
        max_feedback_entries: usize,
    ) -> Result<super::EnthalpyRadiationStepLinearization<'m>, EnthalpyAdjointError> {
        self.inner.linearize_accepted_with_ambient_radiation(
            cx,
            problem,
            interfaces,
            old_specific_h,
            dt_s,
            config,
            patches,
            accepted,
            max_feedback_entries,
        )
    }
}

/// Physical history/source and frozen-reference-density pullbacks of one step.
#[derive(Debug, Clone, PartialEq)]
pub struct HeterogeneousEnthalpyStepGradient {
    /// Existing history/source gradients and independently checked transpose solve.
    pub transport: EnthalpyStepGradient,
    /// Derivative with respect to each tetrahedron's independent reference
    /// density [kg/m3], holding specific enthalpy history and the chart fixed.
    pub element_reference_density: Vec<f64>,
    /// Sum of element derivatives for each declared material-table record.
    /// Unused material records receive exactly zero, not an inferred effect.
    pub material_reference_density: Vec<f64>,
}

/// Rechecked heterogeneous endpoint with a bound reference-storage derivative.
///
/// The immutable history difference is captured when the primal is accepted;
/// density pullbacks cannot accidentally combine another history or assignment
/// with this tangent. Geometry, phase charts, conductivity and contact data are
/// frozen. This differentiates reference storage, not density-dependent material
/// properties, phase-chart knots, geometry motion or material selection.
#[derive(Debug)]
pub struct HeterogeneousEnthalpyStepLinearization<'m> {
    inner: EnthalpyStepLinearization<'m>,
    mesh: &'m ConductionMesh,
    specific_enthalpy_change: Vec<f64>,
    element_material_ids: Vec<usize>,
    material_count: usize,
}

impl<'m> HeterogeneousEnthalpyBackwardEuler<'m, '_> {
    /// Solve and revalidate a heterogeneous endpoint before attaching gradients.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_step(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        let inner = self
            .inner
            .linearize_step(cx, problem, interfaces, old_h, dt_s, config)?;
        self.bind_linearization(cx, old_h, inner)
    }

    /// Recompute the supplied endpoint's residual, energy and vertex-specific
    /// chart/branch guards. Public endpoint fields are not acceptance evidence.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        accepted: EnthalpyStepSolution,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        let inner = self
            .inner
            .linearize_accepted(cx, problem, interfaces, old_h, dt_s, config, accepted)?;
        self.bind_linearization(cx, old_h, inner)
    }

    fn bind_linearization(
        &self,
        cx: &Cx<'_>,
        old_h: &[f64],
        inner: EnthalpyStepLinearization<'m>,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        let mut specific_enthalpy_change = Vec::with_capacity(old_h.len());
        for (vertex, (&new, &old)) in inner
            .primal()
            .specific_enthalpy_j_kg
            .iter()
            .zip(old_h)
            .enumerate()
        {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            specific_enthalpy_change.push(finite(new - old)?);
        }
        poll(cx, old_h.len())?;
        Ok(HeterogeneousEnthalpyStepLinearization {
            inner,
            mesh: self.inner.mesh,
            specific_enthalpy_change,
            element_material_ids: self.element_material_ids.clone(),
            material_count: self.material_count,
        })
    }
}

impl HeterogeneousEnthalpyStepLinearization<'_> {
    /// Shared exact endpoint Jacobian, transpose, temperature pullback and
    /// accepted primal. Retains the existing checked sparse-action interface.
    #[must_use]
    pub const fn transport(&self) -> &EnthalpyStepLinearization<'_> {
        &self.inner
    }

    /// Rechecked endpoint with the heterogeneous chart-resolved phase fields.
    #[must_use]
    pub const fn primal(&self) -> &EnthalpyStepSolution {
        self.inner.primal()
    }

    /// Apply `R_rho * direction` at fixed current and previous specific enthalpy.
    /// This is a RESIDUAL Jacobian action, not an implicit endpoint state JVP.
    /// Directions are independent reference-density perturbations per element.
    pub fn apply_reference_density_jacobian(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        poll(cx, 0)?;
        if direction.len() != self.mesh.element_count() {
            return Err(EnthalpyAdjointError::InvalidInput(
                "one reference-density direction per tetrahedron required",
            ));
        }
        let mut output = vec![0.0; self.specific_enthalpy_change.len()];
        for (element, &density_direction) in direction.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let mass_direction =
                finite(self.mesh.element_volume(element) / 4.0 * density_direction)?;
            for &vertex in &self.mesh.complex().tets[element] {
                let vertex = vertex as usize;
                let contribution = finite(mass_direction * self.specific_enthalpy_change[vertex])?;
                output[vertex] = finite(output[vertex] + contribution)?;
            }
        }
        poll(cx, direction.len())?;
        Ok(output)
    }

    /// Solve the physical transpose once, then return source/history gradients
    /// and `-lambda^T R_rho` for both element and shared material densities.
    ///
    /// For each tetrahedron, `dJ/drho_e = -(V_e/4) sum_i lambda_i*(h_i-h_old_i)`.
    /// Every term uses the bound primal history. No residual is divided by a
    /// density or a temperature slope, so latent-plateau storage remains valid.
    /// No partial result is published after cancellation or nonfinite arithmetic.
    pub fn pullback(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        config: LinearConfig,
    ) -> Result<HeterogeneousEnthalpyStepGradient, EnthalpyAdjointError> {
        let transport = self.inner.pullback(cx, seed, config)?;
        let mut element_reference_density = Vec::with_capacity(self.mesh.element_count());
        let mut material_reference_density = vec![0.0; self.material_count];
        for (element, &material) in self.element_material_ids.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let mut sum = 0.0;
            for &vertex in &self.mesh.complex().tets[element] {
                let vertex = vertex as usize;
                let term =
                    finite(transport.adjoint[vertex] * self.specific_enthalpy_change[vertex])?;
                sum = finite(sum + term)?;
            }
            let derivative = finite(-self.mesh.element_volume(element) / 4.0 * sum)?;
            element_reference_density.push(derivative);
            material_reference_density[material] =
                finite(material_reference_density[material] + derivative)?;
        }
        poll(cx, self.mesh.element_count())?;
        Ok(HeterogeneousEnthalpyStepGradient {
            transport,
            element_reference_density,
            material_reference_density,
        })
    }
}
