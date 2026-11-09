//! Reference-mass enthalpy transport with explicit material assignments.
//!
//! Each tetrahedron names a declared phase chart and frozen reference density.
//! An explicitly bound positive P1 nodal profile can replace the initial
//! material-constant density without replacing any vertex's caloric chart.
//! A shared vertex must have exactly one material identity. Different materials
//! therefore meet through separate traces and an explicitly declared thermal
//! contact, not through an inferred mixture of enthalpy charts. Conductivity
//! assignments remain owned by `ConductionProblem::element_materials`, unless
//! an original phase-source assignment explicitly binds them to these same IDs.
//!
//! This is a stationary reference-mass model: equilibrium density does not
//! change storage, geometry does not move, and no phase-chart blending,
//! material-selection derivative is inferred. Prescribed temperatures use each
//! boundary node's own chart and must select unique enthalpies.

use std::{fmt, sync::Arc};

mod radiation;
mod material_ambient;
pub use material_ambient::HeterogeneousMaterialAmbientEnthalpyGradient;
mod material_enclosure;
pub use material_enclosure::HeterogeneousMaterialEnclosureEnthalpyGradient;
pub use radiation::{EnthalpyHeatStep, EnthalpyRadiation};

use fs_exec::Cx;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;

use super::adjoint::{EnthalpyAdjointError, EnthalpyStepGradient, EnthalpyStepLinearization};
use super::material::phase::{MaterialPhaseTransport, PhaseFieldTransport};
use super::material::melting::{MaterialMeltingTransport, MeltingFieldTransport};
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

fn admit_nodal_reference_density(
    cx: &Cx<'_>,
    vertices: usize,
    density_kg_m3: &[f64],
) -> Result<(), EnthalpyError> {
    if density_kg_m3.len() != vertices {
        return Err(EnthalpyError::InvalidInput(
            "one nodal reference density per mesh vertex required",
        ));
    }
    for (vertex, &density) in density_kg_m3.iter().enumerate() {
        if vertex % ASSEMBLY_TILE == 0 {
            poll(cx, vertex)?;
        }
        if !density.is_finite() || density <= 0.0 {
            return Err(EnthalpyError::InvalidInput(
                "positive finite nodal reference densities required",
            ));
        }
    }
    Ok(())
}

/// Prepared heterogeneous masses and nodal phase-chart ownership.
///
/// The uniform stepper is private, so its single-chart accessor cannot be used
/// to misrepresent this model. All public chart inspection is vertex-specific.
#[derive(Debug)]
pub struct HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    inner: EnthalpyBackwardEuler<'m, 'c>,
    element_material_ids: Arc<[usize]>,
    vertex_material_ids: Arc<[usize]>,
    material_count: usize,
    reference_materials: Vec<ReferenceEnthalpyMaterial<'c>>,
    nodal_reference_density_kg_m3: Option<Vec<f64>>,
}

impl<'m, 'c> HeterogeneousEnthalpyBackwardEuler<'m, 'c> {
    /// Add sensible, perfectly bonded film storage while retaining every
    /// original host material chart and reference mass.
    pub fn with_surface_heat_capacity(
        mut self,
        cx: &Cx<'_>,
        surfaces: &[crate::transient::backward_euler::SurfaceHeatCapacity<'_>],
    ) -> Result<Self, EnthalpyError> {
        self.inner = self.inner.with_surface_heat_capacity(cx, surfaces)?;
        Ok(self)
    }

    /// Add original nonlinear face inventories to the same heterogeneous host
    /// storage, without averaging charts or changing reference masses.
    pub fn with_surface_heat_storage(
        mut self, cx: &Cx<'_>, surfaces: &[super::SurfaceHeatStorage<'c>],
    ) -> Result<Self, EnthalpyError> {
        self.inner = self.inner.with_surface_heat_storage(cx, surfaces)?;
        Ok(self)
    }

    /// Original nonlinear face inventories, preserving source/table identity.
    #[must_use]
    pub fn surface_heat_storage(&self) -> &[super::SurfaceHeatStorage<'c>] {
        self.inner.surface_heat_storage()
    }

    /// Constant film row capacities in original mesh vertex order, J/K.
    /// Nonlinear `surface_heat_storage` inventories contribute additionally.
    #[must_use]
    pub fn surface_nodal_heat_capacity_j_per_k(&self) -> Option<&[f64]> {
        self.inner.surface_nodal_heat_capacity_j_per_k()
    }

    /// Original constant density/Cp film sources.
    #[must_use]
    pub fn surface_material_capacities(&self) -> &[crate::transient::MaterialHeatCapacity] {
        self.inner.surface_material_capacities()
    }

    /// Second-order SDIRK2 with the same original heterogeneous caloric charts,
    /// frozen reference masses and optional complete source conductivity.
    pub fn advance_sdirk2(&self, cx: &Cx<'_>, inputs: [super::sdirk::EnthalpySdirk2Input<'_>; 2],
        old_h: &[f64], dt_s: f64, config: EnthalpyStepConfig)
        -> Result<super::sdirk::EnthalpySdirk2Step, EnthalpyError>
    {
        self.inner.advance_sdirk2(cx, inputs, old_h, dt_s, config)
    }

    /// Advance both original heterogeneous SDIRK2 stages with additional
    /// physical entropy admission on their native time quadrature.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_sdirk2_with_entropy(
        &self, cx: &Cx<'_>, inputs: [super::sdirk::EnthalpySdirk2Input<'_>; 2],
        old_h: &[f64], dt_s: f64, config: EnthalpyStepConfig, entropy_tolerance_j_k: f64,
    ) -> Result<(super::sdirk::EnthalpySdirk2Step, super::entropy::EnthalpyEntropyBalance),
        super::entropy::EnthalpyEntropyError> {
        self.inner.advance_sdirk2_with_entropy(cx, inputs, old_h, dt_s, config, entropy_tolerance_j_k)
    }

    /// Solve and bind both SDIRK2 stages on the original heterogeneous charts.
    /// Reference inventory, material parameters and the time grid stay fixed.
    /// The radiation entry allowance covers both retained stage tangents.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_sdirk2(&self, cx: &Cx<'_>, inputs: [super::sdirk::EnthalpySdirk2Input<'_>; 2],
        old_h: &[f64], dt_s: f64, config: EnthalpyStepConfig, max_radiation_entries: usize)
        -> Result<super::adjoint::sdirk::EnthalpySdirk2Linearization<'m>, super::adjoint::EnthalpyAdjointError>
    {
        self.inner.linearize_sdirk2(cx, inputs, old_h, dt_s, config, max_radiation_entries)
    }

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
    pub fn new(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        materials: &[ReferenceEnthalpyMaterial<'c>],
        element_material_ids: &[usize],
        budget: EnthalpyBudget,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        Self::prepare_reference_inventory(cx, mesh, materials, element_material_ids, None, budget)
    }

    #[allow(clippy::too_many_lines)]
    fn prepare_reference_inventory(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        materials: &[ReferenceEnthalpyMaterial<'c>],
        element_material_ids: &[usize],
        nodal_density: Option<&[f64]>,
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
        if let Some(density) = nodal_density {
            admit_nodal_reference_density(cx, n, density)?;
        }
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
            let local = if let Some(density) = nodal_density {
                super::reference_mass::p1_nodal_masses(
                    mesh.element_volume(element),
                    vertices.map(|vertex| density[vertex as usize]),
                )
                .ok_or(EnthalpyError::InvalidInput(
                    "positive nodal reference mass is not representable",
                ))?
            } else {
                let contribution =
                    finite(material.reference_density_kg_m3 * mesh.element_volume(element) / 4.0)?;
                if contribution <= 0.0 {
                    return Err(EnthalpyError::InvalidInput(
                        "positive reference mass is not representable",
                    )
                    .into());
                }
                [contribution; 4]
            };
            for (&vertex, contribution) in vertices.iter().zip(local) {
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
                surface_capacity: None,
                phase_conductivity: None,
                sourced_melting: None,
                sourced_phase: None,
            },
            element_material_ids: Arc::from(element_material_ids),
            vertex_material_ids: vertex_material_ids.into(),
            material_count: materials.len(),
            reference_materials: materials.to_vec(),
            nodal_reference_density_kg_m3: nodal_density.map(<[f64]>::to_vec),
        })
    }

    /// Replace prepared storage by a frozen positive P1 reference-density field.
    ///
    /// Exact row masses are `sum_e V_e*(sum_b rho_b + rho_a)/20`, retaining
    /// density/enthalpy correlation. Original material IDs and borrowed charts
    /// remain unchanged. This is a reference-inventory declaration before
    /// evolution, not a current-density update of an accepted heated specimen.
    /// Scalar material-source bindings deliberately become unavailable through
    /// [`Self::reference_material`]; source owners must explicitly support the
    /// declared profile instead of treating the old scalar as its inventory.
    /// Existing preparation budgets bound this mesh-sized work and allocation.
    ///
    /// # Errors
    /// Refuses an already bound scalar phase-source inventory, wrong length,
    /// nonpositive/nonfinite density, unrepresentable
    /// positive mass, arithmetic overflow or cancellation before publishing.
    pub fn with_nodal_reference_density(
        mut self,
        cx: &Cx<'_>,
        density_kg_m3: &[f64],
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0)?;
        if self.inner.sourced_phase.is_some() || self.inner.sourced_melting.is_some() {
            return Err(EnthalpyError::InvalidInput(
                "nodal reference density requires an explicitly profile-aware original source",
            ).into());
        }
        let n = self.inner.mesh.vertex_count();
        admit_nodal_reference_density(cx, n, density_kg_m3)?;
        let mut masses = vec![0.0; n];
        for (element, vertices) in self.inner.mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 { poll(cx, element)?; }
            let density = vertices.map(|vertex| density_kg_m3[vertex as usize]);
            let local = super::reference_mass::p1_nodal_masses(
                self.inner.mesh.element_volume(element), density,
            ).ok_or(EnthalpyError::InvalidInput(
                "positive nodal reference mass is not representable",
            ))?;
            for (&vertex, contribution) in vertices.iter().zip(local) {
                let vertex = vertex as usize;
                masses[vertex] = finite(masses[vertex] + contribution)?;
            }
        }
        poll(cx, n)?;
        self.inner.masses = masses;
        self.nodal_reference_density_kg_m3 = Some(density_kg_m3.to_vec());
        Ok(self)
    }

    /// Original explicitly declared P1 density in mesh vertex order, kg/m3.
    /// Material-constant storage returns None. Equilibrium densities never
    /// replace this immutable reference profile.
    #[must_use]
    pub fn nodal_reference_density_kg_m3(&self) -> Option<&[f64]> {
        self.nodal_reference_density_kg_m3.as_deref()
    }

    /// Bind each region's original sharp-melting calorics, phase-specific k(T),
    /// volume-fraction coexistence closure and single-phase reference density.
    /// Distinct materials retain separate nodal traces and declared contact;
    /// their enthalpy datums and phase charts are never averaged.
    pub fn from_melting_materials(
        cx: &Cx<'_>, mesh: &'m ConductionMesh, materials: &'c [MaterialMeltingTransport],
        element_material_ids: &[usize], budget: EnthalpyBudget,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        Self::prepare_melting_materials(cx, mesh, materials, element_material_ids, None, budget)
    }

    /// Bind the same original melting sources with an explicit positive P1
    /// reference-density field. This initial inventory supersedes the source's
    /// scalar reference-density answer; equilibrium density still affects only
    /// mass-to-volume phase conversion, never the closed specimen's mass.
    pub fn from_melting_materials_with_nodal_reference_density(
        cx: &Cx<'_>, mesh: &'m ConductionMesh, materials: &'c [MaterialMeltingTransport],
        element_material_ids: &[usize], density_kg_m3: &[f64], budget: EnthalpyBudget,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        Self::prepare_melting_materials(cx, mesh, materials, element_material_ids, Some(density_kg_m3), budget)
    }

    fn prepare_melting_materials(
        cx: &Cx<'_>, mesh: &'m ConductionMesh, materials: &'c [MaterialMeltingTransport],
        element_material_ids: &[usize], nodal_density: Option<&[f64]>, budget: EnthalpyBudget,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0)?;
        if materials.len() > budget.max_elements {
            return Err(EnthalpyError::Budget { resource: "materials", required: materials.len(),
                limit: budget.max_elements }.into());
        }
        let mut reference = Vec::with_capacity(materials.len());
        for (index, material) in materials.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 { poll(cx, index)?; }
            reference.push(ReferenceEnthalpyMaterial { curve: material.caloric().curve(),
                reference_density_kg_m3: material.reference_density_kg_m3() });
        }
        let mut storage = Self::prepare_reference_inventory(cx, mesh, &reference,
            element_material_ids, nodal_density, budget)?;
        storage.inner.sourced_melting = Some(MeltingFieldTransport::Assigned {
            materials, vertex_material_ids: Arc::clone(&storage.vertex_material_ids),
        });
        Ok(storage)
    }

    /// Original caloric, conductivity, coexistence and density-source answers.
    /// With an explicit nodal density profile, the source's scalar reference
    /// answer is provenance only; `reference_nodal_masses_kg` owns the inventory.
    #[must_use]
    pub fn melting_transport_source(&self, material: usize) -> Option<&'c MaterialMeltingTransport> {
        self.inner.sourced_melting.as_ref()?.material(material)
    }

    /// Physical problem using the exact original melting-source assignment.
    /// The shared unit tensor is only the source assembly's numerical basis.
    /// Extra material tensors or phase multipliers remain explicitly refused.
    pub fn melting_transport_problem<'s>(
        &'s self, boundary: &'s crate::ThermalBoundary, source: &'s crate::ScalarField,
    ) -> Result<ConductionProblem<'s>, EnthalpyError> {
        let melting = self.inner.sourced_melting.as_ref().ok_or(EnthalpyError::InvalidInput(
            "heterogeneous transport has no bound original melting sources"))?;
        Ok(ConductionProblem { mesh: self.inner.mesh, boundary,
            material: melting.numerical_basis(), element_materials: None, source })
    }

    /// Actual element-average heat flux using the same source-specific nodal
    /// k(h), phase volume mixture and original T(h) as every native step.
    /// Contact flux retains its existing interface report and sign convention.
    pub fn melting_transport_heat_flux(&self, cx: &Cx<'_>, h: &[f64])
        -> Result<Vec<[f64; 3]>, EnthalpyError>
    {
        self.inner.sourced_melting.as_ref().ok_or(EnthalpyError::InvalidInput(
            "heterogeneous transport has no bound original melting sources"))?
            .element_heat_flux(cx, self.inner.mesh, h)
    }

    /// Bind one original isotropic phase source and fixed reference density
    /// per material to the existing heterogeneous enthalpy owner.
    ///
    /// Storage and conductivity share the same admitted vertex/element IDs.
    /// Different material IDs require disjoint traces and declared contacts,
    /// even when their source values or caloric charts happen to coincide.
    /// Sources retain their own enthalpy datum, support and original chart.
    pub fn from_phase_materials(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        materials: &'c [MaterialPhaseTransport<'c>],
        element_material_ids: &[usize],
        budget: EnthalpyBudget,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        Self::prepare_phase_materials(cx, mesh, materials, element_material_ids, None, budget)
    }

    /// Bind original isotropic k(T,f) sources to an explicitly declared initial
    /// P1 reference-density inventory on the same heterogeneous mesh.
    ///
    /// The positive nodal profile supplies the actual exact reference masses
    /// `m_i = sum_e V_e*(sum_a rho_a + rho_i)/20`. Material scalar densities
    /// are superseded at preparation; they never enter mass assembly or any
    /// subsequent solve. The original caloric charts, source claims, support,
    /// material assignments and conductivity moments remain unchanged.
    /// Density does not infer a porosity, conductivity or mixture correction.
    ///
    /// This declares inventory BEFORE evolution. Neither equilibrium rho(T)
    /// nor later source rebinding changes the closed specimen's mass. Scalar
    /// inventory accessors return None; [`Self::phase_conductivity_source`]
    /// exposes each retained source independently of the spatial mass field.
    /// Existing preparation budgets and cancellation cover all added work.
    /// Invalid profiles, material assignments or unrepresentable mass refuse
    /// before returning a prepared owner.
    pub fn from_phase_materials_with_nodal_reference_density(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        materials: &'c [MaterialPhaseTransport<'c>],
        element_material_ids: &[usize],
        density_kg_m3: &[f64],
        budget: EnthalpyBudget,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        Self::prepare_phase_materials(
            cx,
            mesh,
            materials,
            element_material_ids,
            Some(density_kg_m3),
            budget,
        )
    }

    fn prepare_phase_materials(
        cx: &Cx<'_>,
        mesh: &'m ConductionMesh,
        materials: &'c [MaterialPhaseTransport<'c>],
        element_material_ids: &[usize],
        nodal_density: Option<&[f64]>,
        budget: EnthalpyBudget,
    ) -> Result<Self, HeterogeneousEnthalpyError> {
        poll(cx, 0)?;
        // Bound the temporary reference inventory before allocating it; the
        // original constructor admits the complete spatial assignment below.
        if materials.len() > budget.max_elements {
            return Err(EnthalpyError::Budget {
                resource: "materials",
                required: materials.len(),
                limit: budget.max_elements,
            }.into());
        }
        let mut reference = Vec::with_capacity(materials.len());
        for (index, material) in materials.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            reference.push(ReferenceEnthalpyMaterial {
                curve: material.conductivity_source().caloric_curve(),
                reference_density_kg_m3: material.reference_density_kg_m3(),
            });
        }
        let mut storage = Self::prepare_reference_inventory(
            cx,
            mesh,
            &reference,
            element_material_ids,
            nodal_density,
            budget,
        )?;
        storage.inner.sourced_phase = Some(PhaseFieldTransport::Assigned {
            materials,
            element_material_ids: Arc::clone(&storage.element_material_ids),
            vertex_material_ids: Arc::clone(&storage.vertex_material_ids),
        });
        Ok(storage)
    }

    /// Original source and reference density at a retained material-table ID.
    /// Returns None for an unknown ID, absent phase sources or an explicit
    /// nodal profile: the source's scalar density is not that field's inventory.
    #[must_use]
    pub fn phase_transport_material(
        &self,
        material: usize,
    ) -> Option<&'c MaterialPhaseTransport<'c>> {
        if self.nodal_reference_density_kg_m3.is_some() {
            None
        } else {
            self.inner.sourced_phase.as_ref()?.material(material)
        }
    }

    /// Original isotropic or first principal conductivity/caloric source,
    /// without a scalar mass assertion. For the complete directional source
    /// set and its shared support use [`Self::phase_conductivity_sources`] and
    /// [`Self::phase_conductivity_enthalpy_bounds_j_kg`].
    /// Available for both scalar and P1 reference inventories; returns None
    /// for an unknown material ID or storage without bound phase sources.
    #[must_use]
    pub fn phase_conductivity_source(
        &self,
        material: usize,
    ) -> Option<&'c PreparedPhaseConductivity<'c>> {
        Some(
            self.inner
                .sourced_phase
                .as_ref()?
                .material(material)?
                .conductivity_source(),
        )
    }

    /// Every original directional conductivity source, in its declared order.
    /// One entry means explicit isotropy; three retain the material's fixed
    /// principal frame. Scalar and P1 reference inventories share these sources.
    #[must_use]
    pub fn phase_conductivity_sources(
        &self,
        material: usize,
    ) -> Option<&'c [&'c PreparedPhaseConductivity<'c>]> {
        Some(
            self.inner
                .sourced_phase
                .as_ref()?
                .material(material)?
                .conductivity_sources(),
        )
    }

    /// Complete conductivity support at this material, intersecting all
    /// principal sources. It is independent of the reference-density profile.
    #[must_use]
    pub fn phase_conductivity_enthalpy_bounds_j_kg(&self, material: usize) -> Option<[f64; 2]> {
        Some(
            self.inner
                .sourced_phase
                .as_ref()?
                .material(material)?
                .enthalpy_bounds_j_kg(),
        )
    }

    /// Physical problem for the original phase-source assignment and exact mesh.
    /// The shared unit numerical basis does not replace any physical source.
    /// Separate tensor assignments or phase multipliers are refused at admission.
    pub fn phase_transport_problem<'s>(
        &'s self,
        boundary: &'s crate::ThermalBoundary,
        source: &'s crate::ScalarField,
    ) -> Result<ConductionProblem<'s>, EnthalpyError> {
        let phase = self.inner.sourced_phase.as_ref().ok_or(EnthalpyError::InvalidInput(
            "heterogeneous transport has no bound original phase sources",
        ))?;
        Ok(ConductionProblem {
            mesh: self.inner.mesh,
            boundary,
            material: phase.numerical_basis(),
            element_materials: None,
            source,
        })
    }

    /// Element-average physical heat flux from each material's original k(h_h)
    /// and its own nodal chart temperatures. Contact exchange remains in the
    /// existing contact reports; this is not a face-conservative recovery.
    pub fn phase_transport_heat_flux(
        &self,
        cx: &Cx<'_>,
        h: &[f64],
    ) -> Result<Vec<[f64; 3]>, EnthalpyError> {
        let phase = self.inner.sourced_phase.as_ref().ok_or(EnthalpyError::InvalidInput(
            "heterogeneous transport has no bound original phase sources",
        ))?;
        phase.element_heat_flux(cx, self.inner.mesh, h)
    }

    /// Exact stationary reference mesh used by storage and transport assembly.
    #[must_use]
    pub const fn mesh(&self) -> &'m ConductionMesh {
        self.inner.mesh
    }

    /// Original borrowed caloric chart at a retained material-table ID.
    ///
    /// Available for scalar and P1 reference inventories, including unused
    /// table entries. This does not assert a scalar density or mass inventory;
    /// those remain owned by the original reference-material or nodal profile.
    #[must_use]
    pub fn reference_curve(&self, material: usize) -> Option<&'c EquilibriumEnthalpyPhaseCurve> {
        self.reference_materials.get(material).map(|entry| entry.curve)
    }

    /// Original borrowed caloric chart and frozen reference density. Coupled
    /// Lagrangian owners must use this inventory, not a fresh density at T.
    /// Returns None for an explicit nodal profile: no scalar material density
    /// represents that declared inventory, even when its values happen to match.
    #[must_use]
    pub fn reference_material(&self, material: usize) -> Option<ReferenceEnthalpyMaterial<'c>> {
        if self.nodal_reference_density_kg_m3.is_some() {
            None
        } else {
            self.reference_materials.get(material).copied()
        }
    }

    /// Whether a separately declared phase multiplier modifies conductivity.
    /// Source-bound consumers can refuse an additional, unbound coefficient.
    #[must_use]
    pub const fn has_phase_conductivity(&self) -> bool {
        self.inner.phase_conductivity.is_some()
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

    /// Prepare adaptive time integration with the same per-node phase charts,
    /// reference masses and phase-dependent conductivity assignments.
    pub fn adaptive_state<'s>(
        &'s self,
        cx: &Cx<'_>,
        initial_specific_h: &[f64],
        initial_time_s: f64,
        initial_step_s: f64,
        config: super::adaptive::AdaptiveEnthalpyConfig,
    ) -> Result<super::adaptive::AdaptiveEnthalpyState<'s, 'm, 'c>, EnthalpyError> {
        self.inner.adaptive_state(
            cx,
            initial_specific_h,
            initial_time_s,
            initial_step_s,
            config,
        )
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
            self.inner.curve_for_vertex(vertex).values_at_specific_enthalpy(h)
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

    /// Stationary isolated equilibrium with each original chart, reference mass
    /// and sensible film retained. Connected inventories come from the native
    /// physical heat operator; latent phase allocations remain explicit.
    pub fn isolated_equilibrium<'s>(
        &'s self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'s>,
        interfaces: Option<&'s ThermalInterfaces>,
        old_h: &'s [f64],
        config: super::equilibrium::IsolatedEquilibriumConfig,
    ) -> Result<super::equilibrium::IsolatedEnthalpyEquilibrium<'s, 'm, 'c>,
        super::equilibrium::IsolatedEquilibriumError> {
        self.inner.isolated_equilibrium(cx, problem, interfaces, old_h, config)
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

    /// Add signed assembled interval heat to the original heterogeneous nodal
    /// residual, keeping every vertex's original caloric chart and inventory.
    /// Existing volumetric sources, contact and boundary heat remain additive.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_nodal_heat(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[fs_math::dd::Dd],
    ) -> Result<super::nodal_heat::NodalHeatStep, EnthalpyError> {
        self.inner.advance_with_nodal_heat(
            cx, problem, interfaces, old_specific_h, dt_s, config, nodal_heat_j,
        )
    }

    /// Measure an endpoint with freshly evaluated signed nodal interval heat.
    /// This performs no solve and applies no residual/energy acceptance gate;
    /// an outer constitutive owner must check both local and total balances.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_nodal_heat(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[fs_math::dd::Dd],
    ) -> Result<super::nodal_heat::NodalHeatEndpoint, EnthalpyError> {
        self.inner.evaluate_endpoint_with_nodal_heat(
            cx, problem, interfaces, old_specific_h, candidate_specific_h,
            dt_s, config, nodal_heat_j,
        )
    }

    /// Combine original heterogeneous storage and signed nodal interval heat
    /// with native pointwise radiation. Every material-emissivity owner must
    /// retain the same original nodal caloric objects; no averaged chart is used.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_nodal_heat_and_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[fs_math::dd::Dd],
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<super::nodal_heat::NodalHeatStep, EnthalpyError> {
        self.inner.advance_with_nodal_heat_and_radiation(
            cx, problem, interfaces, old_specific_h, dt_s, config, nodal_heat_j, radiation,
        )
    }

    /// Independently measure the actual radiative endpoint with a freshly
    /// supplied nodal heat vector, retaining original heterogeneous inventories.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_nodal_heat_and_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        nodal_heat_j: &[fs_math::dd::Dd],
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<super::nodal_heat::NodalHeatEndpoint, EnthalpyError> {
        self.inner.evaluate_endpoint_with_nodal_heat_and_radiation(
            cx, problem, interfaces, old_specific_h, candidate_specific_h,
            dt_s, config, nodal_heat_j, radiation,
        )
    }

    /// Combine original heterogeneous storage and signed nodal heat with live
    /// source-backed temperature/phase contacts and optional pointwise radiation.
    /// All original nodal caloric objects and held contact coordinates survive.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_nodal_heat_and_temperature_contacts(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &crate::interface::temperature::TemperatureContacts<'_, '_>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig,
        nodal_heat_j: &[fs_math::dd::Dd],
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<super::nodal_heat::NodalHeatStep, EnthalpyError> {
        self.inner.advance_with_nodal_heat_and_temperature_contacts(
            cx, problem, contacts, old_specific_h, dt_s, config, nodal_heat_j, radiation,
        )
    }

    /// Independently re-query live contacts and radiation at the actual endpoint.
    /// This measures the fresh supplied heat without solving or admitting it.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_nodal_heat_and_temperature_contacts(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &crate::interface::temperature::TemperatureContacts<'_, '_>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig,
        nodal_heat_j: &[fs_math::dd::Dd],
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<super::nodal_heat::NodalHeatEndpoint, EnthalpyError> {
        self.inner.evaluate_endpoint_with_nodal_heat_and_temperature_contacts(
            cx, problem, contacts, old_specific_h, candidate_specific_h,
            dt_s, config, nodal_heat_j, radiation,
        )
    }

    /// Advance with the original per-vertex caloric charts and an additional
    /// discrete entropy admission. This shares the uniform owner's operator,
    /// storage and transaction; no homogenized material or entropy is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_entropy(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_specific_h: &[f64],
        dt_s: f64, config: EnthalpyStepConfig, entropy_tolerance_j_k: f64,
    ) -> Result<(EnthalpyStepSolution, super::entropy::EnthalpyEntropyBalance),
        super::entropy::EnthalpyEntropyError> {
        self.inner.advance_with_entropy(cx, problem, interfaces, old_specific_h,
            dt_s, config, entropy_tolerance_j_k)
    }

    /// Exact original-chart entropy change, including any retained surface films.
    /// This evaluates storage only; it does not assemble or solve a heat equation.
    pub fn thermal_entropy_change_j_k(
        &self,
        cx: &Cx<'_>,
        old_specific_h: &[f64],
        new_specific_h: &[f64],
    ) -> Result<f64, EnthalpyError> {
        self.inner
            .thermal_entropy_change_j_k(cx, old_specific_h, new_specific_h)
    }

    /// Advance with the original source-qualified radiation and caloric charts,
    /// admitting the native endpoint's physical entropy account.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_radiation_and_entropy(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: EnthalpyRadiation<'_>,
        entropy_tolerance_j_k: f64,
    ) -> Result<
        (EnthalpyStepSolution, super::entropy::EnthalpyEntropyBalance),
        super::entropy::EnthalpyEntropyError,
    > {
        self.inner.advance_with_radiation_and_entropy(
            cx, problem, interfaces, old_specific_h, dt_s, config,
            radiation, entropy_tolerance_j_k,
        )
    }

    /// Measure a supplied endpoint with the original physical radiation owner.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_radiation_and_entropy(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: EnthalpyRadiation<'_>,
    ) -> Result<
        (super::EnthalpyEndpointEvaluation, super::entropy::EnthalpyEntropyBalance),
        EnthalpyError,
    > {
        self.inner.evaluate_endpoint_with_radiation_and_entropy(
            cx, problem, interfaces, old_specific_h, candidate_specific_h,
            dt_s, config, radiation,
        )
    }

    /// Measure a supplied nonradiative endpoint using each original chart.
    /// Signed discrete generation and residuals remain available without a gate.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_entropy(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old_specific_h: &[f64],
        candidate_specific_h: &[f64], dt_s: f64, config: EnthalpyStepConfig,
    ) -> Result<(super::EnthalpyEndpointEvaluation, super::entropy::EnthalpyEntropyBalance),
        EnthalpyError> {
        self.inner.evaluate_endpoint_with_entropy(cx, problem, interfaces,
            old_specific_h, candidate_specific_h, dt_s, config)
    }

    pub(crate) fn advance_with_temperature_contacts(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &crate::interface::temperature::TemperatureContacts<'_, '_>,
        old_h: &[f64],
        dt_s: f64,
        config: crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig,
    ) -> Result<(EnthalpyStepSolution, crate::interface::temperature::TemperatureContactEvaluation), EnthalpyError> {
        self.inner.advance_with_temperature_contacts(cx, problem, contacts, old_h, dt_s, config, None)
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn advance_with_temperature_contacts_and_pointwise_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &crate::interface::temperature::TemperatureContacts<'_, '_>,
        old_h: &[f64],
        dt_s: f64,
        config: crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig,
        radiation: &crate::PointwiseAmbientRadiation<'_>,
    ) -> Result<(EnthalpyStepSolution, crate::interface::temperature::TemperatureContactEvaluation), EnthalpyError> {
        self.inner.advance_with_temperature_contacts(cx, problem, contacts, old_h, dt_s, config,
            Some(super::PointwiseBoundary::Ambient(radiation)))
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn advance_with_temperature_contacts_and_pointwise_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &crate::interface::temperature::TemperatureContacts<'_, '_>,
        old_h: &[f64],
        dt_s: f64,
        config: crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig,
        radiation: &crate::PointwiseEnclosureRadiation<'_>,
    ) -> Result<(EnthalpyStepSolution, crate::interface::temperature::TemperatureContactEvaluation), EnthalpyError> {
        self.inner.advance_with_temperature_contacts(cx, problem, contacts, old_h, dt_s, config,
            Some(super::PointwiseBoundary::Enclosure(radiation)))
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
    /// Derivative for an additive constant density [kg/m3] in each tetrahedron,
    /// holding its reference profile, specific enthalpy history and chart fixed.
    /// These element controls remain independent for an explicit nodal profile.
    pub element_reference_density: Vec<f64>,
    /// Sum of element derivatives for each declared material-table record.
    /// Unused material records receive exactly zero, not an inferred effect.
    pub material_reference_density: Vec<f64>,
    /// Derivative per nodal P1 reference density, in mesh vertex order.
    /// Present only when that independent profile was explicitly declared.
    pub nodal_reference_density: Option<Vec<f64>>,
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
    element_material_ids: Arc<[usize]>,
    material_count: usize,
    nodal_profile: bool,
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

    /// Revalidate an unchanged endpoint against its original live contact laws
    /// and attach their complete specific-enthalpy tangent, without a new solve.
    ///
    /// The supplied absolute residual target is the original forward physical
    /// gate; it is never reconstructed from endpoint-frozen resistance or Robin
    /// coefficients. The existing energy gate, strict chart/material/source
    /// smoothness checks and exact mesh ownership remain active. Only the mesh
    /// is borrowed by the result, which retains the accepted enthalpy field and
    /// Newton history while recomputing its physical residual and energy.
    ///
    /// The tangent includes frozen contact transport once. The additional source
    /// feedback includes both resistance temperature partials times `dT/dh`
    /// and any independent phase-source derivative. Free rows retain their
    /// prescribed columns, so the existing bath, history, source, reference
    /// density and conductivity pullbacks all use the same complete transpose.
    /// Physical Robin state feedback may subsequently be attached through
    /// [`HeterogeneousEnthalpyStepLinearization::with_robin_feedback`].
    ///
    /// This entry point owns no radiation. `max_face_pairs` includes fixed
    /// contacts and bounds their tangent work, not the total solver workspace.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_temperature_contact_endpoint(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &crate::interface::temperature::TemperatureContacts<'_, '_>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        accepted: EnthalpyStepSolution,
        residual_tolerance_j: f64,
        max_face_pairs: usize,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        self.inner.admit_step(cx, problem, old_h, dt_s, config)?;
        contacts.check_mesh(problem.mesh)?;
        if !residual_tolerance_j.is_finite() || residual_tolerance_j <= 0.0 {
            return Err(EnthalpyAdjointError::InvalidInput(
                "positive finite original contact residual target required",
            ));
        }
        let contact = contacts.linearize_enthalpy(
            cx, problem.boundary, self, &accepted.specific_enthalpy_j_kg, max_face_pairs,
        )?;
        self.linearize_accepted_with_target(
            cx, problem, Some(contact.physical_contacts().frozen_interfaces()),
            old_h, dt_s, config, accepted, residual_tolerance_j,
        )?.with_contact_temperature_feedback(cx, &contact)
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

    /// Bind the actual pointwise ambient endpoint with heterogeneous charts
    /// and frozen reference densities. An explicit coupled-owner residual target
    /// is retained; `None` uses the native solve's original absolute/relative
    /// target, including its prescribed-temperature initial state.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_with_pointwise_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &crate::PointwiseAmbientRadiation<'_>,
        accepted: EnthalpyStepSolution,
        physical_target_j: impl Into<Option<f64>>,
        max_radiation_entries: usize,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        let inner = self.inner.linearize_accepted_with_pointwise_ambient_radiation(
            cx, problem, interfaces, old_h, dt_s, config, radiation, accepted,
            physical_target_j, max_radiation_entries,
        )?;
        self.bind_linearization(cx, old_h, inner)
    }

    /// Bind local P1 emission and full reflected irradiation while preserving
    /// the existing storage-density, source, bath and transport reverse owners.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_accepted_with_pointwise_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &crate::PointwiseEnclosureRadiation<'_>,
        accepted: EnthalpyStepSolution,
        physical_target_j: f64,
        max_radiation_entries: usize,
    ) -> Result<HeterogeneousEnthalpyStepLinearization<'m>, EnthalpyAdjointError> {
        let inner = self.inner.linearize_accepted_with_pointwise_enclosure_radiation(
            cx, problem, interfaces, old_h, dt_s, config, radiation, accepted,
            physical_target_j, max_radiation_entries,
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
            nodal_profile: self.nodal_reference_density_kg_m3.is_some(),
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
    /// Directions are additive constant reference-density perturbations per
    /// element about the frozen profile, including for a nodal profile.
    pub fn apply_reference_density_jacobian(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.density.apply_jacobian(cx, direction)
    }

    /// Residual action for the explicitly declared nodal P1 density controls.
    /// History, endpoint h, chart ownership and geometry remain fixed.
    pub fn apply_nodal_reference_density_jacobian(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.density.apply_nodal_jacobian(cx, direction)
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

    /// Residual action for the bound nodal density profile, retaining the same
    /// radiative endpoint and prescribed-row mask as all other controls.
    pub fn apply_nodal_reference_density_jacobian(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        self.density.apply_nodal_jacobian(cx, direction)
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
    fn apply_nodal_jacobian(
        &self,
        cx: &Cx<'_>,
        direction: &[f64],
    ) -> Result<Vec<f64>, EnthalpyAdjointError> {
        poll(cx, 0)?;
        if !self.nodal_profile || direction.len() != self.specific_enthalpy_change.len() {
            return Err(EnthalpyAdjointError::InvalidInput(
                "one direction per vertex of an explicit nodal density profile required",
            ));
        }
        let mut output = vec![0.0; direction.len()];
        for (element, vertices) in self.mesh.complex().tets.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 { poll(cx, element)?; }
            let local = super::reference_mass::p1_nodal_mass_direction(
                self.mesh.element_volume(element),
                vertices.map(|vertex| direction[vertex as usize]),
            ).ok_or(EnthalpyAdjointError::InvalidInput(
                "finite representable nodal density direction required",
            ))?;
            for (&vertex, mass_direction) in vertices.iter().zip(local) {
                let vertex = vertex as usize;
                let contribution = finite(mass_direction * self.specific_enthalpy_change[vertex])?;
                output[vertex] = finite(output[vertex] + contribution)?;
            }
        }
        poll(cx, direction.len())?;
        Ok(output)
    }

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
        let mut nodal_reference_density = self.nodal_profile
            .then(|| vec![0.0; self.specific_enthalpy_change.len()]);
        for (element, &material) in self.element_material_ids.iter().enumerate() {
            if element % ASSEMBLY_TILE == 0 {
                poll(cx, element)?;
            }
            let mut sum = 0.0;
            let vertices = &self.mesh.complex().tets[element];
            let mut weighted_change = [0.0; 4];
            for (local, &vertex) in vertices.iter().enumerate() {
                let vertex = vertex as usize;
                let term =
                    finite(transport.adjoint[vertex] * self.specific_enthalpy_change[vertex])?;
                weighted_change[local] = term;
                sum = finite(sum + term)?;
            }
            let derivative = finite(-self.mesh.element_volume(element) / 4.0 * sum)?;
            element_reference_density.push(derivative);
            material_reference_density[material] =
                finite(material_reference_density[material] + derivative)?;
            if let Some(nodal) = nodal_reference_density.as_mut() {
                let local = super::reference_mass::p1_nodal_mass_direction(
                    self.mesh.element_volume(element), weighted_change,
                ).ok_or(EnthalpyAdjointError::InvalidInput(
                    "nodal reference-density pullback is not representable",
                ))?;
                for (&vertex, contribution) in vertices.iter().zip(local) {
                    let vertex = vertex as usize;
                    nodal[vertex] = finite(nodal[vertex] - contribution)?;
                }
            }
        }
        poll(cx, self.mesh.element_count())?;
        Ok(HeterogeneousEnthalpyStepGradient {
            transport,
            element_reference_density,
            material_reference_density,
            nodal_reference_density,
        })
    }
}

