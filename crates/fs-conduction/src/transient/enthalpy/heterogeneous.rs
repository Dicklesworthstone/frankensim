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
//! material-selection derivative is inferred. Prescribed temperatures use each
//! boundary node's own chart and must select unique enthalpies.

use std::fmt;

use fs_exec::Cx;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;

use super::adjoint::{EnthalpyAdjointError, EnthalpyStepGradient, EnthalpyStepLinearization};
use super::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyError, EnthalpyStepConfig, EnthalpyStepSolution,
    finite, poll,
};
use super::{EnthalpyRadiationStepGradient, EnthalpyRadiationStepLinearization};
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

    /// Check scheduled bath values without changing the accepted caloric state.
    /// The original transport owner supplies chart inversion and ambiguity rules.
    pub(crate) fn admit_prescribed_boundary(
        &self,
        cx: &Cx<'_>,
        boundary: &crate::ThermalBoundary,
    ) -> Result<(), HeterogeneousEnthalpyError> {
        poll(cx, 0)?;
        for (index, &(vertex, temperature)) in boundary.dirichlet().iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 { poll(cx, index)?; }
            if vertex >= self.vertex_material_ids.len() {
                return Err(EnthalpyError::InvalidInput(
                    "prescribed vertex is outside the reference mesh",
                ).into());
            }
            self.inner.prescribed_enthalpy(vertex, temperature)?;
        }
        Ok(())
    }

    /// Build the same boundary-compatible candidate used to initialize a
    /// physical endpoint solve, without advancing or replacing history.
    /// Coupled boundary owners can evaluate their initial residual at this
    /// candidate while retaining the original `old_h` for storage accounting.
    /// A latent bath temperature never selects a liquid fraction implicitly.
    pub fn prescribed_endpoint_initial_guess(
        &self,
        cx: &Cx<'_>,
        boundary: &crate::ThermalBoundary,
        old_h: &[f64],
    ) -> Result<Vec<f64>, HeterogeneousEnthalpyError> {
        self.admit_prescribed_boundary(cx, boundary)?;
        if old_h.len() != self.vertex_material_ids.len() {
            return Err(EnthalpyError::InvalidInput("initial enthalpy must match the reference mesh").into());
        }
        for (vertex, &h) in old_h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 { poll(cx, vertex)?; }
            self.inner.curve_for_vertex(vertex).state_at_specific_enthalpy(h)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
        }
        let mut candidate = old_h.to_vec();
        for (index, &(vertex, temperature)) in boundary.dirichlet().iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 { poll(cx, index)?; }
            candidate[vertex] = self.inner.prescribed_enthalpy(vertex, temperature)?;
        }
        poll(cx, old_h.len())?;
        Ok(candidate)
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

    /// Evaluate one supplied physical endpoint using its assigned phase charts
    /// and reference masses, without advancing history or taking a Newton step.
    /// The caller owns residual and energy acceptance for the returned values.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        patches: Option<&[crate::AmbientRadiationPatch]>,
    ) -> Result<super::EnthalpyEndpointEvaluation, EnthalpyError> {
        self.inner.evaluate_endpoint(
            cx,
            problem,
            interfaces,
            old_specific_h,
            candidate_specific_h,
            dt_s,
            config,
            patches,
        )
    }

    /// Audit a supplied heterogeneous endpoint with fresh reflected enclosure
    /// exchange, retaining its original external heat and immutable history.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_enclosure(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        candidate_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        enclosure: &crate::radiation::GrayDiffuseEnclosure,
        radiation_config: crate::AmbientRadiationConfig,
    ) -> Result<super::EnthalpyEndpointEvaluation, EnthalpyError> {
        self.inner.evaluate_endpoint_with_enclosure(
            cx, problem, interfaces, old_h, candidate_h, dt_s, config,
            enclosure, radiation_config,
        )
    }

    /// Measure local ambient radiation and the complete physical balance of a
    /// supplied heterogeneous endpoint without advancing its immutable history.
    /// Returns original convection separately from radiation; the caller owns
    /// acceptance of both the local residual and global energy defect.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_pointwise_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        candidate_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &crate::PointwiseAmbientRadiation<'_>,
    ) -> Result<super::EnthalpyEndpointEvaluation, EnthalpyError> {
        self.inner.evaluate_endpoint_with_pointwise_ambient_radiation(
            cx, problem, interfaces, old_h, candidate_h, dt_s, config, radiation,
        )
    }

    /// Measure local emission and reflected enclosure exchange with each
    /// vertex's own chart and reference mass. No tangent or Newton solve is
    /// required; a zero total enclosure heat does not imply nodal equilibrium.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_pointwise_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        candidate_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &crate::PointwiseEnclosureRadiation<'_>,
    ) -> Result<super::EnthalpyEndpointEvaluation, EnthalpyError> {
        self.inner.evaluate_endpoint_with_pointwise_enclosure_radiation(
            cx, problem, interfaces, old_h, candidate_h, dt_s, config, radiation,
        )
    }

    /// Integrate ambient radiation over the local P1 temperature field within
    /// the enthalpy Newton solve, retaining heterogeneous phase charts and masses.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pointwise_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &crate::PointwiseAmbientRadiation<'_>,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        self.inner.advance_with_pointwise_ambient_radiation(
            cx, problem, interfaces, old_specific_h, dt_s, config, radiation,
        )
    }

    /// Resolve local emission and reflected patch irradiation in the original
    /// nonsymmetric enthalpy Newton solve, with each vertex's own phase chart.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pointwise_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &crate::PointwiseEnclosureRadiation<'_>,
    ) -> Result<EnthalpyStepSolution, EnthalpyError> {
        self.inner.advance_with_pointwise_enclosure_radiation(
            cx, problem, interfaces, old_specific_h, dt_s, config, radiation,
        )
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

    /// Recheck an ambient-radiation endpoint and bind its reference-storage
    /// derivatives to this model's exact material assignments and history.
    ///
    /// The full physical radiation tangent, including its mean-temperature
    /// feedback, is retained in the shared transport linearization. History,
    /// source and element/material reference-density pullbacks use that same
    /// checked transpose. Radiation parameters and dt are held fixed; callers
    /// needing their controls can use the separate radiative linearization.
    /// No additional nonlinear endpoint solve is performed.
    ///
    /// `problem.boundary` must contain original convection, without radiation
    /// already folded into it. The physical residual, energy, chart and
    /// emissivity-support gates and two-vector feedback cap are unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_storage_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        patches: &[crate::AmbientRadiationPatch],
        accepted: EnthalpyStepSolution,
        max_feedback_entries: usize,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        Ok(self
            .linearize_accepted_parameters_with_ambient_radiation(
                cx,
                problem,
                interfaces,
                old_h,
                dt_s,
                config,
                patches,
                accepted,
                max_feedback_entries,
            )?
            .into_storage())
    }

    /// Bind reference-storage and radiation-parameter derivatives together.
    ///
    /// The returned owner reuses one physical transpose for history, source,
    /// reference density, reservoir temperature, emissivity and dt controls.
    /// All controls refer to the same revalidated endpoint and immutable
    /// history. Original convection, phase charts and geometry remain fixed.
    /// The same physical residual, energy, smoothness and feedback-entry gates
    /// as the ambient-radiation binder apply; no new primal solve is performed.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_parameters_with_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        patches: &[crate::AmbientRadiationPatch],
        accepted: EnthalpyStepSolution,
        max_feedback_entries: usize,
    ) -> Result<HeterogeneousEnthalpyRadiationStepLinearization<'m>, EnthalpyAdjointError> {
        let inner = self.inner.linearize_accepted_with_ambient_radiation(
            cx,
            problem,
            interfaces,
            old_h,
            dt_s,
            config,
            patches,
            accepted,
            max_feedback_entries,
        )?;
        let density = self.bind_reference_density(cx, old_h, inner.transport())?;
        Ok(HeterogeneousEnthalpyRadiationStepLinearization { inner, density })
    }

    pub(crate) fn linearize_accepted_parameters_with_ambient_radiation_target(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        patches: &[crate::AmbientRadiationPatch],
        accepted: EnthalpyStepSolution,
        max_feedback_entries: usize,
        physical_target_j: f64,
    ) -> Result<HeterogeneousEnthalpyRadiationStepLinearization<'m>, EnthalpyAdjointError> {
        if !physical_target_j.is_finite() || physical_target_j <= 0.0 {
            return Err(EnthalpyAdjointError::InvalidInput(
                "positive finite outer radiative/contact residual target required",
            ));
        }
        let inner = self.inner.linearize_accepted_with_ambient_radiation_target(
            cx, problem, interfaces, old_h, dt_s, config, patches, accepted,
            max_feedback_entries, Some(physical_target_j),
        )?;
        let density = self.bind_reference_density(cx, old_h, inner.transport())?;
        Ok(HeterogeneousEnthalpyRadiationStepLinearization { inner, density })
    }

    /// Advance closed enclosure exchange with heterogeneous charts and the
    /// same fixed reference masses, contacts and phase conductivity laws.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        enclosure: &crate::radiation::GrayDiffuseEnclosure,
        radiation_config: crate::AmbientRadiationConfig,
    ) -> Result<super::enclosure::EnthalpyEnclosureStepSolution, EnthalpyError> {
        self.inner.advance_with_enclosure_radiation(
            cx,
            problem,
            interfaces,
            old_h,
            dt_s,
            config,
            enclosure,
            radiation_config,
        )
    }

    /// Bind complete reflected exchange and the same immutable density history.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_parameters_with_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        enclosure: &crate::radiation::GrayDiffuseEnclosure,
        radiation_config: crate::AmbientRadiationConfig,
        accepted: EnthalpyStepSolution,
        max_feedback_entries: usize,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        let inner = self.inner.linearize_accepted_with_enclosure_radiation(
            cx, problem, interfaces, old_h, dt_s, config, enclosure,
            radiation_config, accepted, max_feedback_entries,
        )?;
        self.bind_linearization(cx, old_h, inner)
    }

    /// Recheck and retain the complete enclosure tangent in h coordinates;
    /// original convection ports remain distinct from radiation exchange.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_with_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        enclosure: &crate::radiation::GrayDiffuseEnclosure,
        radiation_config: crate::AmbientRadiationConfig,
        accepted: EnthalpyStepSolution,
        max_feedback_entries: usize,
    ) -> Result<EnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        self.inner.linearize_accepted_with_enclosure_radiation(
            cx,
            problem,
            interfaces,
            old_h,
            dt_s,
            config,
            enclosure,
            radiation_config,
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
    density: ReferenceDensityLinearization<'m>,
}

/// Physical storage and radiative-parameter derivatives from one transpose.
/// The spatial source/history multiplier is stored only once, in `storage`.
#[derive(Debug, Clone, PartialEq)]
pub struct HeterogeneousEnthalpyRadiationStepGradient {
    /// History, P1 source and bound element/material reference-density controls.
    pub storage: HeterogeneousEnthalpyStepGradient,
    /// Derivative per second of dt, holding endpoint source/boundary data fixed.
    /// A schedule or horizon derivative requires its separate data chain rule.
    pub time_step_s: f64,
    /// Radiating regions in the original patch declaration order.
    pub radiation_regions: Vec<String>,
    /// Derivative per kelvin of each black-reservoir temperature, including
    /// coefficient and explicit reference-temperature changes.
    pub reservoir_temperatures: Vec<f64>,
    /// Derivative per unit absolute hemispherical emissivity; the source claim
    /// and query support are fixed. At epsilon=1 only feasible directions apply.
    pub emissivities: Vec<f64>,
}

/// Rechecked radiative endpoint with bound storage and parameter derivatives.
/// Both radiation metadata and reference-density history belong to this exact
/// endpoint; callers cannot combine them with an unrelated tangent or state.
#[derive(Debug)]
pub struct HeterogeneousEnthalpyRadiationStepLinearization<'m> {
    inner: EnthalpyRadiationStepLinearization<'m>,
    density: ReferenceDensityLinearization<'m>,
}

#[derive(Debug)]
struct ReferenceDensityLinearization<'m> {
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

    pub(crate) fn linearize_accepted_with_target(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        accepted: EnthalpyStepSolution,
        physical_target_j: f64,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        if !physical_target_j.is_finite() || physical_target_j <= 0.0 {
            return Err(EnthalpyAdjointError::InvalidInput(
                "positive finite outer enthalpy residual target required",
            ));
        }
        let inner = self.inner.linearize_accepted_with_target(
            cx, problem, interfaces, old_h, dt_s, config, accepted,
            Some(physical_target_j),
        )?;
        self.bind_linearization(cx, old_h, inner)
    }

    fn bind_linearization(
        &self,
        cx: &Cx<'_>,
        old_h: &[f64],
        inner: EnthalpyStepLinearization<'m>,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        let density = self.bind_reference_density(cx, old_h, &inner)?;
        Ok(HeterogeneousEnthalpyStepLinearization { inner, density })
    }

    fn bind_reference_density(
        &self,
        cx: &Cx<'_>,
        old_h: &[f64],
        step: &EnthalpyStepLinearization<'_>,
    ) -> Result<ReferenceDensityLinearization<'m>, EnthalpyAdjointError> {
        let mut specific_enthalpy_change = Vec::with_capacity(old_h.len());
        for (vertex, (&new, &old)) in step.primal().specific_enthalpy_j_kg.iter().zip(old_h).enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            // At a fixed row dR/drho scales h-h_bath, which is exactly zero;
            // its physical storage jump belongs to the reservoir reaction.
            specific_enthalpy_change.push(if step.is_prescribed(vertex) { 0.0 } else { finite(new - old)? });
        }
        poll(cx, old_h.len())?;
        Ok(ReferenceDensityLinearization {
            mesh: self.inner.mesh,
            specific_enthalpy_change,
            element_material_ids: self.element_material_ids.clone(),
            material_count: self.material_count,
        })
    }
}

impl HeterogeneousEnthalpyStepLinearization<'_> {
    pub(crate) fn with_contact_temperature_feedback(
        mut self,
        cx: &Cx<'_>,
        contact: &crate::interface::temperature::tangent::ContactTemperatureLinearization,
    ) -> Result<Self, EnthalpyAdjointError> {
        self.inner = self.inner.with_contact_temperature_feedback(cx, contact)?;
        Ok(self)
    }

    /// Attach original convective wall feedback without detaching the bound
    /// reference-density history. Ordering, units and the total nodal-factor
    /// cap are those of [`EnthalpyStepLinearization::with_robin_feedback`].
    /// All source, history and density controls retain one complete multiplier.
    pub fn with_robin_feedback(
        mut self,
        cx: &Cx<'_>,
        regions: &[&str],
        reference_derivatives: &[Vec<f64>],
        htc_wall_derivatives: &[f64],
        max_feedback_entries: usize,
    ) -> Result<Self, EnthalpyAdjointError> {
        self.inner = self.inner.with_robin_feedback(
            cx,
            regions,
            reference_derivatives,
            htc_wall_derivatives,
            max_feedback_entries,
        )?;
        Ok(self)
    }

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
        self.density.apply_jacobian(cx, direction)
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
        self.density.pullback(cx, transport)
    }
}

impl<'m> HeterogeneousEnthalpyRadiationStepLinearization<'m> {
    pub(crate) fn with_contact_temperature_feedback(
        mut self,
        cx: &Cx<'_>,
        contact: &crate::interface::temperature::tangent::ContactTemperatureLinearization,
    ) -> Result<Self, EnthalpyAdjointError> {
        self.inner = self.inner.with_contact_temperature_feedback(cx, contact)?;
        Ok(self)
    }

    /// Attach original convection laws to the full radiation/storage tangent.
    /// [`EnthalpyStepLinearization::with_robin_feedback`] defines the matrix
    /// ordering, units and total nodal-factor cap. Radiation parameters and
    /// reference-density controls remain bound to this endpoint and one dual.
    pub fn with_robin_feedback(
        mut self,
        cx: &Cx<'_>,
        regions: &[&str],
        reference_derivatives: &[Vec<f64>],
        htc_wall_derivatives: &[f64],
        max_feedback_entries: usize,
    ) -> Result<Self, EnthalpyAdjointError> {
        self.inner = self.inner.with_robin_feedback(
            cx,
            regions,
            reference_derivatives,
            htc_wall_derivatives,
            max_feedback_entries,
        )?;
        Ok(self)
    }

    /// Rechecked physical endpoint with chart-resolved phase fields.
    #[must_use]
    pub fn primal(&self) -> &EnthalpyStepSolution {
        self.inner.primal()
    }

    /// Shared endpoint Jacobian and temperature pullback, including complete
    /// nonlinear radiation feedback in specific-enthalpy coordinates.
    #[must_use]
    pub const fn transport(&self) -> &EnthalpyStepLinearization<'m> {
        self.inner.transport()
    }

    /// Retain h/source/reference-density derivatives while discarding only
    /// the additional radiation-parameter metadata. No state feedback is lost.
    #[must_use]
    pub fn into_storage(self) -> HeterogeneousEnthalpyStepLinearization<'m> {
        HeterogeneousEnthalpyStepLinearization {
            inner: self.inner.into_transport(),
            density: self.density,
        }
    }

    /// Apply the bound reference-density residual Jacobian at fixed h history.
    pub fn apply_reference_density_jacobian(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.density.apply_jacobian(cx, direction)
    }

    /// Solve the complete radiative transpose once. Contract radiation controls
    /// and the bound reference-density controls with that same multiplier.
    pub fn pullback(
        &self,
        cx: &Cx<'_>,
        seed: &[f64],
        config: LinearConfig,
    ) -> Result<HeterogeneousEnthalpyRadiationStepGradient, EnthalpyAdjointError> {
        let EnthalpyRadiationStepGradient {
            transport,
            time_step_s,
            radiation_regions,
            reservoir_temperatures,
            emissivities,
        } = self.inner.pullback(cx, seed, config)?;
        let storage = self.density.pullback(cx, transport)?;
        Ok(HeterogeneousEnthalpyRadiationStepGradient {
            storage,
            time_step_s,
            radiation_regions,
            reservoir_temperatures,
            emissivities,
        })
    }
}

impl ReferenceDensityLinearization<'_> {
    fn apply_jacobian(
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

    // This multiplier comes only from the paired private transport owner.
    // Move it into the result so all controls report one transpose solve.
    fn pullback(
        &self,
        cx: &Cx<'_>,
        transport: EnthalpyStepGradient,
    ) -> Result<HeterogeneousEnthalpyStepGradient, EnthalpyAdjointError> {
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
