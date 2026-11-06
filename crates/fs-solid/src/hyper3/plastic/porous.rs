//! Retained finite pressure-sensitive compaction on the native P1 tetrahedra.
//!
//! The original porous law owns its non-isochoric plastic metric, associated
//! return and analytic tangent. The existing tetrahedral geometry, Newton and
//! globalization own equilibrium. Every trial borrows the same accepted
//! history; only a fully admitted endpoint publishes a complete next solid.
//! Plastic volume is never normalized or reconstructed from total volume.
//!
//! This is the fixed-condition Deshpande--Fleck/Hencky model declared by
//! `fs_material::plastic::finite::porous`. Its pressure sensitivity does not
//! imply evolving porosity, cell-wall densification, asymmetric tension and
//! compression, damage, pore-fluid pressure or measured foam calibration.
//! Reference-volume physical yield loss is separate from the signed endpoint
//! work defect. [`dynamics`] adds the original consistent reference mass and
//! backward-Euler momentum owner. No contact, thermal coupling or mixed-pressure
//! element is added here.

pub mod dynamics;

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::plastic::finite::porous::{
    FinitePorous, FinitePorousError, FinitePorousMaterial, FinitePorousState, FrozenPorousTrial,
};

use super::{
    AcceptedMotion, FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetEnergy,
    FinitePlasticTetError, FinitePlasticTetState, IDENTITY, checkpoint, energy_gate, invalid,
    reference_size, validate_controls,
};
use crate::hyper3::{
    HyperTetError, HyperTetProblem, HyperTetSettings, HyperTetSolution, HyperelasticLawField,
};
use crate::linear3::TetAssemblyBudget;

/// Original geometry, laws and complete accepted porous material histories.
///
/// The virgin constructor and immutable `advance` are the only state-changing
/// paths. Neither a partial return nor a failed global solve can change an
/// accepted plastic metric. Prepared reference elements are shared across
/// candidates, while each original tetrahedron retains its own history.
#[derive(Debug)]
pub struct FinitePorousTetSolid<'a, L: FinitePorousMaterial = FinitePorous> {
    material: FinitePlasticLawField<'a, L>,
    state: FinitePlasticTetState<'a, FinitePorousState>,
    jacobians: Vec<f64>,
    plastic_jacobians: Vec<f64>,
}

/// Complete accepted endpoint plus its native equilibrium and physical work.
#[derive(Debug)]
pub struct FinitePorousTetStep<'a, L: FinitePorousMaterial = FinitePorous> {
    pub solid: FinitePorousTetSolid<'a, L>,
    /// `strain_energy_j` is the reduced incremental potential, not storage.
    pub equilibrium: HyperTetSolution,
    /// Physical yield heat, retained storage and numerical work remain distinct.
    pub energy: FinitePlasticTetEnergy,
}

/// Refusal from the original mechanical owner or native porous law.
#[derive(Clone, Debug, PartialEq)]
pub enum FinitePorousTetError {
    /// Shared geometry, equilibrium, budget, cancellation or work-gate refusal.
    Solid(FinitePlasticTetError),
    /// The original law refused an accepted or proposed element history.
    Material {
        element: usize,
        source: FinitePorousError,
    },
}

impl core::fmt::Display for FinitePorousTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(error) => error.fmt(f),
            Self::Material { element, source } => {
                write!(f, "finite porous tetrahedron {element}: {source}")
            }
        }
    }
}

impl std::error::Error for FinitePorousTetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            Self::Material { source, .. } => Some(source),
        }
    }
}

impl From<FinitePlasticTetError> for FinitePorousTetError {
    fn from(error: FinitePlasticTetError) -> Self {
        Self::Solid(error)
    }
}

impl From<HyperTetError> for FinitePorousTetError {
    fn from(error: HyperTetError) -> Self {
        Self::Solid(error.into())
    }
}

impl<'a, L: FinitePorousMaterial> FinitePorousTetSolid<'a, L> {
    /// Prepare an undeformed, stress-free solid with original per-element laws.
    /// Native and source-resolved owners remain borrowed; complete source cards
    /// and fixed physical contexts are retained without repeated source queries.
    ///
    /// # Errors
    /// Refuses invalid geometry, law cardinality, budget, or cancellation before
    /// exposing a state. Material coefficients and bounds are never inferred.
    pub fn try_new(
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        material: FinitePlasticLawField<'a, L>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FinitePorousTetError> {
        checkpoint(cx)?;
        material.validate_count(tetrahedra.len())?;
        reference_size(nodes_m.len(), tetrahedra.len(), budget)?;
        let mut element_states = Vec::with_capacity(tetrahedra.len());
        for index in 0..tetrahedra.len() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            element_states.push(material.at(index).finite_porous().initial_state());
        }
        let frozen = Self::freeze(&material, &element_states, cx)?;
        let elements = FinitePlasticTetState::reference_elements(
            nodes_m,
            tetrahedra,
            HyperelasticLawField::PerElement(&frozen),
            budget,
            cx,
        )?;
        Ok(Self {
            material,
            state: FinitePlasticTetState {
                nodes_m,
                tetrahedra,
                elements,
                budget,
                displacement_m: vec![[0.0; 3]; nodes_m.len()],
                element_deformation: vec![IDENTITY; tetrahedra.len()],
                element_states,
                element_first_piola_pa: vec![[0.0; 9]; tetrahedra.len()],
                stored_energy_j: 0.0,
            },
            jacobians: vec![1.0; tetrahedra.len()],
            plastic_jacobians: vec![1.0; tetrahedra.len()],
        })
    }

    fn freeze<'s>(
        material: &'s FinitePlasticLawField<'_, L>,
        states: &'s [FinitePorousState],
        cx: &Cx<'_>,
    ) -> Result<Vec<FrozenPorousTrial<'s>>, FinitePorousTetError> {
        let mut frozen = Vec::with_capacity(states.len());
        for (element, state) in states.iter().enumerate() {
            if element % 1024 == 0 {
                checkpoint(cx)?;
            }
            frozen.push(
                material
                    .at(element)
                    .finite_porous()
                    .frozen(state)
                    .map_err(|source| FinitePorousTetError::Material { element, source })?,
            );
        }
        Ok(frozen)
    }

    /// Original immutable reference coordinates, in metres.
    pub fn reference_nodes_m(&self) -> &'a [[f64; 3]] {
        self.state.nodes_m
    }

    pub fn reference_tetrahedra(&self) -> &'a [[usize; 4]] {
        self.state.tetrahedra
    }

    pub fn assembly_budget(&self) -> TetAssemblyBudget {
        self.state.budget
    }

    /// Accepted total displacement in original node order.
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.state.displacement_m
    }

    /// Accepted total deformation gradients, row-major, in element order.
    pub fn element_deformation(&self) -> &[[f64; 9]] {
        &self.state.element_deformation
    }

    /// Complete non-isochoric plastic metrics and work-conjugate accumulation.
    pub fn element_states(&self) -> &[FinitePorousState] {
        &self.state.element_states
    }

    pub fn element_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.state.element_first_piola_pa
    }

    /// Total current/reference volume ratios, directly from the native response.
    pub fn element_jacobians(&self) -> &[f64] {
        &self.jacobians
    }

    /// Permanent plastic volume ratios det(Fp), independent of elastic volume.
    pub fn element_plastic_jacobians(&self) -> &[f64] {
        &self.plastic_jacobians
    }

    /// Transform an explicitly supplied original reference density field to
    /// the accepted current configuration using mass conservation rho=rho_ref/J.
    ///
    /// One positive finite value is required per original tetrahedron. This is
    /// a query, not a mass update or material-source calibration. The caller
    /// must supply the same original inventory throughout its mechanical path.
    /// Neither the relaxed plastic volume Jp nor current density replaces that
    /// original reference inventory. No partial field is returned on refusal.
    pub fn current_density_kg_m3(
        &self,
        reference_density_kg_m3: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FinitePorousTetError> {
        checkpoint(cx)?;
        if reference_density_kg_m3.len() != self.jacobians.len() {
            return Err(invalid("one original reference density per tetrahedron required").into());
        }
        let mut densities = Vec::with_capacity(self.jacobians.len());
        for (index, (&reference, &j)) in reference_density_kg_m3
            .iter()
            .zip(&self.jacobians)
            .enumerate()
        {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            let density = reference / j;
            if !reference.is_finite() || reference <= 0.0 || !density.is_finite() || density <= 0.0
            {
                return Err(
                    invalid("reference and current densities must be finite and positive").into(),
                );
            }
            densities.push(density);
        }
        checkpoint(cx)?;
        Ok(densities)
    }

    /// Elastic and energetic-hardening storage over the reference mesh, J.
    pub fn stored_energy_j(&self) -> f64 {
        self.state.stored_energy_j
    }

    pub fn material(&self) -> FinitePlasticLawField<'a, L> {
        self.material
    }

    /// Solve one absolute support/dead-load endpoint from accepted history.
    ///
    /// Loading, unloading and reversal use the same original reference and laws.
    /// All material candidates, work and yield gates must pass before a complete
    /// next solid is returned. `self` remains valid and unchanged on success as
    /// well as every refusal. Body force is per original reference volume.
    ///
    /// # Errors
    /// Requires one native load step and finite nonnegative work/yield bounds.
    /// Native material, geometry, equilibrium, cancellation and work refusals
    /// preserve the previously accepted solid.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<FinitePorousTetStep<'a, L>, FinitePorousTetError> {
        self.advance_loaded(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            controls,
            None,
            cx,
        )
        .map(|(step, _)| step)
    }

    #[allow(
        clippy::too_many_arguments,
        clippy::too_many_lines,
        clippy::type_complexity
    )]
    fn advance_loaded(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        settings: HyperTetSettings,
        controls: FinitePlasticStepControls,
        inertia: Option<&crate::inertia::BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<(FinitePorousTetStep<'a, L>, Option<AcceptedMotion>), FinitePorousTetError> {
        checkpoint(cx)?;
        validate_controls(settings, controls)?;
        for (element, state) in self.state.element_states.iter().enumerate() {
            if element % 1024 == 0 {
                checkpoint(cx)?;
            }
            self.material
                .at(element)
                .finite_porous()
                .validate_committed(
                    &self.state.element_deformation[element],
                    state,
                    controls.yield_tolerance_pa,
                )
                .map_err(|source| FinitePorousTetError::Material { element, source })?;
        }
        let frozen = Self::freeze(&self.material, &self.state.element_states, cx)?;
        let (equilibrium, force) = self.state.solve_endpoint(
            HyperelasticLawField::PerElement(&frozen),
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            settings,
            inertia,
            cx,
        )?;
        let new_u = equilibrium.displacement_m.as_flattened();
        let motion = inertia
            .map(|inertia| inertia.motion(new_u, cx))
            .transpose()
            .map_err(crate::hyper3::inertia_error)?;
        let mut element_states = Vec::with_capacity(self.state.elements.len());
        let mut element_deformation = Vec::with_capacity(self.state.elements.len());
        let mut element_first_piola_pa = Vec::with_capacity(self.state.elements.len());
        let mut jacobians = Vec::with_capacity(self.state.elements.len());
        let mut plastic_jacobians = Vec::with_capacity(self.state.elements.len());
        let mut energy = FinitePlasticTetEnergy {
            stored_before_j: self.state.stored_energy_j,
            stored_after_j: 0.0,
            plastic_dissipation_j: 0.0,
            incremental_potential_j: 0.0,
            dead_load_endpoint_work_j: 0.0,
            support_endpoint_work_j: 0.0,
            internal_endpoint_work_j: 0.0,
            equilibrium_work_residual_j: 0.0,
            endpoint_work_defect_j: 0.0,
        };
        for (index, element) in self.state.elements.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            let law = self.material.at(index).finite_porous();
            let f = HyperTetProblem::<FrozenPorousTrial<'_>>::deformation_gradient(element, new_u);
            let refused = |source| FinitePorousTetError::Material {
                element: index,
                source,
            };
            let response = law
                .response(&f, &self.state.element_states[index])
                .map_err(refused)?;
            law.validate_committed(&f, &response.state, controls.yield_tolerance_pa)
                .map_err(refused)?;
            // Evaluate storage at the candidate without another return. In
            // particular, its plastic determinant must not be projected to one.
            let stored = law.stored_energy(&f, &response.state).map_err(refused)?;
            energy_gate(
                "porous constitutive stored-energy consistency",
                element.volume * (stored - response.stored_energy_j_m3),
                controls.equilibrium_work_tolerance_j,
            )?;
            energy.stored_after_j = element.volume.mul_add(stored, energy.stored_after_j);
            energy.plastic_dissipation_j = element.volume.mul_add(
                response.plastic_dissipation_j_m3,
                energy.plastic_dissipation_j,
            );
            energy.incremental_potential_j = element.volume.mul_add(
                response.incremental_potential_j_m3,
                energy.incremental_potential_j,
            );
            element_deformation.push(f);
            element_first_piola_pa.push(response.first_piola_pa);
            element_states.push(response.state);
            jacobians.push(response.jacobian);
            plastic_jacobians.push(response.plastic_jacobian);
        }
        let work = self.state.endpoint_work(
            &equilibrium,
            &force,
            &element_deformation,
            &element_first_piola_pa,
            motion.as_ref().map(|(_, account)| account.inertial_work_j),
            cx,
        )?;
        energy.internal_endpoint_work_j = work.internal_j;
        energy.dead_load_endpoint_work_j = work.dead_load_j;
        energy.support_endpoint_work_j = work.support_j;
        energy.equilibrium_work_residual_j = work.residual_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j
            - (energy.stored_after_j - energy.stored_before_j)
            - energy.plastic_dissipation_j;
        energy_gate(
            "porous assembled incremental-potential consistency",
            equilibrium.strain_energy_j - energy.incremental_potential_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "porous storage/dissipation potential consistency",
            energy.incremental_potential_j - energy.stored_after_j - energy.plastic_dissipation_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "porous support/dead-load equilibrium work",
            energy.equilibrium_work_residual_j,
            controls.equilibrium_work_tolerance_j,
        )?;
        energy_gate(
            "porous signed endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if !energy.stored_after_j.is_finite()
            || energy.stored_after_j < 0.0
            || !energy.plastic_dissipation_j.is_finite()
            || energy.plastic_dissipation_j < 0.0
        {
            return Err(invalid(
                "finite nonnegative porous storage and plastic dissipation required",
            )
            .into());
        }
        checkpoint(cx)?;
        Ok((
            FinitePorousTetStep {
                solid: Self {
                    material: self.material,
                    state: FinitePlasticTetState {
                        nodes_m: self.state.nodes_m,
                        tetrahedra: self.state.tetrahedra,
                        elements: Arc::clone(&self.state.elements),
                        budget: self.state.budget,
                        displacement_m: equilibrium.displacement_m.clone(),
                        element_deformation,
                        element_states,
                        element_first_piola_pa,
                        stored_energy_j: energy.stored_after_j,
                    },
                    jacobians,
                    plastic_jacobians,
                },
                equilibrium,
                energy,
            },
            motion,
        ))
    }
}
