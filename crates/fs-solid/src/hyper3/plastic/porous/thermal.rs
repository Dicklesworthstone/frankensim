//! Thermal porous compaction on the original finite tetrahedral owner.
//!
//! The native porous thermal family maps total deformation through its thermal
//! stretch while retaining the non-isochoric plastic metric. Each Newton solve
//! freezes only endpoint temperature; old history is admitted at its OLD
//! temperature. Heating can therefore drive additional collapse at fixed total
//! deformation. Reference geometry, optional inertia and the original source
//! families survive every accepted endpoint.
//!
//! Shared thermoplastic energy accounts separate mechanical free/internal
//! storage, physical yield heat, reciprocal thermal exchange, signed endpoint
//! quadrature error and numerical kinetic loss. Only physical yield loss plus
//! reversible heat enters a caloric owner, once. This mechanical endpoint does
//! not itself solve heat transport. P1 locking, dry associated compaction and
//! the declared affine-temperature material limits remain explicit.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::plastic::finite::porous::FinitePorousState;
use fs_material::plastic::finite::porous::thermal::{
    AffineThermalFinitePorous, FinitePorousThermalError, FinitePorousThermalIncrement,
    FiniteThermalPorousMaterial, FrozenThermalFinitePorous,
};
use fs_math::dd::Dd;

use super::super::{
    FinitePlasticTetError, FinitePlasticTetState, IDENTITY, checkpoint, energy_gate, invalid,
    reference_size, validate_controls,
};
use super::FinitePorousPressureEnergy;
use crate::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use crate::hyper3::{
    HyperTetError, HyperTetProblem, HyperTetSettings, HyperTetSolution, HyperelasticLawField,
    inertia_error,
};
use crate::inertia::{BackwardEulerInertia, P1ReferenceMass};
use crate::linear3::TetAssemblyBudget;

/// The existing material and numerical-inertia budgets retain their meanings.
pub use super::super::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticElementEnergy, FiniteThermoplasticTetEnergy,
    InertialEnergy,
};

#[derive(Debug)]
struct Motion {
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_nodal_density_kg_m3: Option<Arc<[f64]>>,
    mass: Arc<P1ReferenceMass>,
}

/// Accepted total deformation, temperature, plastic volume and original sources.
#[derive(Debug)]
pub struct FiniteThermalPorousTetSolid<
    'a,
    L: FiniteThermalPorousMaterial = AffineThermalFinitePorous,
> {
    state: FinitePlasticTetState<'a, FinitePorousState>,
    families: Arc<[L]>,
    temperatures_k: Vec<f64>,
    internal_energy_j: f64,
    jacobians: Vec<f64>,
    plastic_jacobians: Vec<f64>,
    motion: Option<Motion>,
}

/// Complete proposal; failure leaves the borrowed original solid unchanged.
#[derive(Debug)]
pub struct FiniteThermalPorousTetStep<
    'a,
    L: FiniteThermalPorousMaterial = AffineThermalFinitePorous,
> {
    pub solid: FiniteThermalPorousTetSolid<'a, L>,
    /// Native forces and reactions. `strain_energy_j` is the return potential.
    pub equilibrium: HyperTetSolution,
    /// Shared reference-volume accounts; the creep subaccount is zero here.
    pub element_energy: Vec<FiniteThermoplasticElementEnergy>,
    pub energy: FiniteThermoplasticTetEnergy,
    pub inertia: Option<InertialEnergy>,
    /// Original current-area pressure work and its separate geometric endpoint
    /// correction. Material endpoint defect and physical caloric heat are
    /// unchanged; the complete energy balance includes this pressure work.
    pub pressure: Option<FinitePorousPressureEnergy>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteThermalPorousTetError {
    Solid(FinitePlasticTetError),
    Material {
        element: usize,
        source: FinitePorousThermalError,
    },
}
impl From<FinitePlasticTetError> for FiniteThermalPorousTetError {
    fn from(error: FinitePlasticTetError) -> Self {
        Self::Solid(error)
    }
}
impl From<HyperTetError> for FiniteThermalPorousTetError {
    fn from(error: HyperTetError) -> Self {
        Self::Solid(error.into())
    }
}
impl core::fmt::Display for FiniteThermalPorousTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(error) => error.fmt(f),
            Self::Material { element, source } => {
                write!(f, "finite thermal porous tetrahedron {element}: {source}")
            }
        }
    }
}
impl std::error::Error for FiniteThermalPorousTetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            Self::Material { source, .. } => Some(source),
        }
    }
}

struct FrozenPoint<'a> {
    law: &'a FrozenThermalFinitePorous,
    accepted: &'a FinitePorousState,
}
fn material_error(error: FinitePorousThermalError) -> MaterialError {
    MaterialError::State {
        what: error.to_string(),
    }
}
impl FiniteStrainPotential for FrozenPoint<'_> {
    fn energy_density(&self, f: &[f64; 9]) -> f64 {
        self.law
            .stress_and_potential(f, self.accepted)
            .map_or(f64::NAN, |(_, energy)| energy)
    }
    fn first_piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        self.law
            .stress_and_potential(f, self.accepted)
            .map(|(stress, _)| stress)
            .map_err(material_error)
    }
    fn consistent_tangent(&self, f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        self.law
            .response(f, self.accepted)
            .map(|response| response.tangent_pa)
            .map_err(material_error)
    }
    fn stress_and_potential(&self, f: &[f64; 9]) -> Result<([f64; 9], f64), MaterialError> {
        self.law
            .stress_and_potential(f, self.accepted)
            .map_err(material_error)
    }
}

fn account(value: &FinitePorousThermalIncrement) -> FiniteThermoplasticElementEnergy {
    FiniteThermoplasticElementEnergy {
        free_before_j_m3: value.old.free_energy_j_m3,
        free_after_j_m3: value.new.free_energy_j_m3,
        internal_before_j_m3: value.old.internal_energy_j_m3,
        internal_after_j_m3: value.new.internal_energy_j_m3,
        plastic_dissipation_j_m3: value.mechanical.plastic_dissipation_j_m3,
        creep_dissipation_j_m3: 0.0,
        reversible_heat_j_m3: value.reversible_heat_j_m3,
        coefficient_storage_change_j_m3: value.coefficient_storage_change_j_m3,
        thermoelastic_storage_change_j_m3: value.thermoelastic_storage_change_j_m3,
        endpoint_work_defect_j_m3: value.endpoint_work_defect_j_m3,
        incremental_potential_j_m3: value.mechanical.incremental_potential_j_m3,
        internal_endpoint_work_j_m3: value.mechanical_work_j_m3,
        mechanical_coordinate_work_j_m3: value.mechanical_coordinate_work_j_m3,
        thermal_deformation_work_j_m3: value.thermal_deformation_work_j_m3,
        balance_error_j_m3: value.balance_error_j_m3,
    }
}

impl<'a, L: FiniteThermalPorousMaterial> FiniteThermalPorousTetSolid<'a, L> {
    /// Admit a virgin reference body with one complete thermal family per cell.
    /// Initial temperatures must equal each family's stress-free reference.
    pub fn try_new(
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        families: Vec<L>,
        initial_temperatures_k: &[f64],
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermalPorousTetError> {
        checkpoint(cx)?;
        reference_size(nodes_m.len(), tetrahedra.len(), budget)?;
        if families.len() != tetrahedra.len() || initial_temperatures_k.len() != tetrahedra.len() {
            return Err(
                invalid("one original family and temperature per porous cell required").into(),
            );
        }
        let mut laws = Vec::with_capacity(families.len());
        let mut states = Vec::with_capacity(families.len());
        for (element, family) in families.iter().enumerate() {
            if element % 256 == 0 {
                checkpoint(cx)?;
            }
            let family = family.thermal_finite_porous();
            let temperature = initial_temperatures_k[element];
            if temperature != family.parameters().reference_temperature_k {
                return Err(
                    invalid("virgin porous temperature must equal its source reference").into(),
                );
            }
            laws.push(
                family
                    .at_temperature(temperature)
                    .map_err(|source| FiniteThermalPorousTetError::Material { element, source })?,
            );
            states.push(family.initial_state());
        }
        let frozen: Vec<_> = laws
            .iter()
            .zip(&states)
            .map(|(law, accepted)| FrozenPoint { law, accepted })
            .collect();
        let elements = FinitePlasticTetState::reference_elements(
            nodes_m,
            tetrahedra,
            HyperelasticLawField::PerElement(&frozen),
            budget,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(Self {
            state: FinitePlasticTetState {
                nodes_m,
                tetrahedra,
                elements,
                budget,
                displacement_m: vec![[0.0; 3]; nodes_m.len()],
                element_deformation: vec![IDENTITY; tetrahedra.len()],
                element_states: states,
                element_first_piola_pa: vec![[0.0; 9]; tetrahedra.len()],
                stored_energy_j: 0.0,
            },
            families: Arc::from(families),
            temperatures_k: initial_temperatures_k.to_vec(),
            internal_energy_j: 0.0,
            jacobians: vec![1.0; tetrahedra.len()],
            plastic_jacobians: vec![1.0; tetrahedra.len()],
            motion: None,
        })
    }

    fn admit_motion(
        &self,
        velocity: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<(), FiniteThermalPorousTetError> {
        checkpoint(cx)?;
        if self.motion.is_some() || velocity.len() != self.state.nodes_m.len() {
            return Err(
                invalid("unbound inertia and one velocity per original node required").into(),
            );
        }
        for (i, value) in velocity.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            if value.iter().any(|v| !v.is_finite()) {
                return Err(invalid("initial porous velocity must be finite").into());
            }
        }
        Ok(())
    }

    /// Bind original cellwise density once; compaction and heating conserve mass.
    pub fn with_inertia(
        mut self,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermalPorousTetError> {
        self.admit_motion(velocity_mps, cx)?;
        let mass = P1ReferenceMass::new(
            self.state.elements.iter().map(|e| (e.nodes, e.volume)),
            reference_density_kg_m3,
            self.state.displacement_m.as_flattened().len(),
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        self.motion = Some(Motion {
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(reference_density_kg_m3),
            reference_nodal_density_kg_m3: None,
            mass: Arc::new(mass),
        });
        Ok(self)
    }

    /// Bind exact consistent P1 nodal-density inertia on the original mesh.
    pub fn with_nodal_inertia(
        mut self,
        reference_nodal_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermalPorousTetError> {
        self.admit_motion(velocity_mps, cx)?;
        let (mass, means) = P1ReferenceMass::from_nodal_density(
            self.state.elements.iter().map(|e| (e.nodes, e.volume)),
            reference_nodal_density_kg_m3,
            self.state.displacement_m.as_flattened().len(),
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        self.motion = Some(Motion {
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(means),
            reference_nodal_density_kg_m3: Some(Arc::from(reference_nodal_density_kg_m3)),
            mass: Arc::new(mass),
        });
        Ok(self)
    }

    pub fn reference_nodes_m(&self) -> &'a [[f64; 3]] {
        self.state.nodes_m
    }
    pub fn reference_tetrahedra(&self) -> &'a [[usize; 4]] {
        self.state.tetrahedra
    }
    pub fn assembly_budget(&self) -> TetAssemblyBudget {
        self.state.budget
    }
    pub fn families(&self) -> &[L] {
        &self.families
    }
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.state.displacement_m
    }
    pub fn element_deformation(&self) -> &[[f64; 9]] {
        &self.state.element_deformation
    }
    pub fn element_states(&self) -> &[FinitePorousState] {
        &self.state.element_states
    }
    pub fn element_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.state.element_first_piola_pa
    }
    pub fn element_jacobians(&self) -> &[f64] {
        &self.jacobians
    }
    pub fn element_plastic_jacobians(&self) -> &[f64] {
        &self.plastic_jacobians
    }
    pub fn free_energy_j(&self) -> f64 {
        self.state.stored_energy_j
    }
    pub fn internal_energy_j(&self) -> f64 {
        self.internal_energy_j
    }
    pub fn reference_volume_m3(&self, element: usize) -> Option<f64> {
        self.state.elements.get(element).map(|value| value.volume)
    }
    pub fn velocity_mps(&self) -> Option<&[[f64; 3]]> {
        self.motion
            .as_ref()
            .map(|motion| motion.velocity_mps.as_slice())
    }
    pub fn reference_density_kg_m3(&self) -> Option<&[f64]> {
        self.motion
            .as_ref()
            .map(|motion| motion.reference_density_kg_m3.as_ref())
    }
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.motion
            .as_ref()
            .and_then(|motion| motion.reference_nodal_density_kg_m3.as_deref())
    }
    pub fn reference_mass_kg(&self) -> Option<f64> {
        self.motion
            .as_ref()
            .map(|motion| motion.mass.reference_mass_kg())
    }

    /// Apply mass conservation to the actual total volume, including expansion.
    /// Reference density is an explicit inventory; the thermal law adds no mass.
    pub fn current_density_kg_m3(
        &self,
        reference_density_kg_m3: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteThermalPorousTetError> {
        checkpoint(cx)?;
        if reference_density_kg_m3.len() != self.jacobians.len() {
            return Err(invalid("one original density per porous cell required").into());
        }
        let mut result = Vec::with_capacity(self.jacobians.len());
        for (i, (&rho, &j)) in reference_density_kg_m3
            .iter()
            .zip(&self.jacobians)
            .enumerate()
        {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            let current = rho / j;
            if !rho.is_finite() || rho <= 0.0 || !current.is_finite() || current <= 0.0 {
                return Err(
                    invalid("positive finite reference and current densities required").into(),
                );
            }
            result.push(current);
        }
        checkpoint(cx)?;
        Ok(result)
    }

    /// Project a prescribed original-node acceleration with the retained M*a.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FiniteThermalPorousTetError> {
        checkpoint(cx)?;
        let motion = self
            .motion
            .as_ref()
            .ok_or_else(|| invalid("porous acceleration requires bound reference inertia"))?;
        if acceleration_m_s2.len() != self.state.nodes_m.len()
            || acceleration_m_s2.iter().flatten().any(|a| !a.is_finite())
        {
            return Err(invalid("one finite acceleration per original node required").into());
        }
        let forces = motion
            .mass
            .action(acceleration_m_s2.as_flattened(), cx)
            .map_err(inertia_error)?;
        let mut loads = Vec::with_capacity(forces.len());
        for (dof, force) in forces.into_iter().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            if force != 0.0 {
                loads.push((dof, force));
            }
        }
        checkpoint(cx)?;
        Ok(loads)
    }

    /// Propose a mechanical endpoint at supplied element temperatures.
    /// The caloric owner closes the returned physical heat against its own
    /// storage and transport; numerical losses are never deposited into it.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        settings: HyperTetSettings,
        controls: FiniteThermoplasticControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermalPorousTetStep<'a, L>, FiniteThermalPorousTetError> {
        self.advance_loaded(
            new_temperatures_k,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            dt_s,
            settings,
            controls,
            None,
            cx,
        )
    }

    /// Retain thermal compaction and optional original inertia under pressure
    /// applied to the CURRENT surface, using the native force and follower tangent.
    ///
    /// The prescribed pressure field is held over the interval; its original
    /// owner supplies uniform -p delta V or declared-path work for other fields.
    /// Original mesh binding, temperature/history admission and thermal accounts
    /// remain unchanged. Physical caloric heat is still yield heat plus reciprocal
    /// thermal exchange: pressure work, geometric work correction and numerical
    /// kinetic loss are not additional heat sources.
    ///
    /// The material work budget bounds the sum of absolute material endpoint
    /// defect and pressure correction, without cancellation. Every refused
    /// temperature/source, solve, final work or cancellation gate leaves the
    /// accepted solid and its original inertia unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_pressure<'p>(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        surface: &ClosedTetPressureSurface,
        pressure: impl Into<TetPressure<'p>>,
        settings: HyperTetSettings,
        controls: FiniteThermoplasticControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermalPorousTetStep<'a, L>, FiniteThermalPorousTetError> {
        self.advance_loaded(
            new_temperatures_k,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            dt_s,
            settings,
            controls,
            Some((surface, pressure.into())),
            cx,
        )
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_loaded(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        dt_s: f64,
        settings: HyperTetSettings,
        controls: FiniteThermoplasticControls,
        pressure: Option<(&ClosedTetPressureSurface, TetPressure<'_>)>,
        cx: &Cx<'_>,
    ) -> Result<FiniteThermalPorousTetStep<'a, L>, FiniteThermalPorousTetError> {
        checkpoint(cx)?;
        validate_controls(settings, controls.material)?;
        if let Some((surface, field)) = pressure {
            if !surface.matches_reference(self.state.nodes_m, self.state.tetrahedra) {
                return Err(invalid(
                    "thermal porous pressure surface must bind the exact original mesh",
                )
                .into());
            }
            surface.validate_pressure(field, cx)?;
        }
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || !controls.maximum_numerical_kinetic_loss_j.is_finite()
            || controls.maximum_numerical_kinetic_loss_j < 0.0
            || new_temperatures_k.len() != self.families.len()
        {
            return Err(invalid("positive interval, nonnegative kinetic budget and one temperature per porous cell required").into());
        }
        let mut laws = Vec::with_capacity(self.families.len());
        for (element, family) in self.families.iter().enumerate() {
            if element % 256 == 0 {
                checkpoint(cx)?;
            }
            let family = family.thermal_finite_porous();
            let refused = |source| FiniteThermalPorousTetError::Material { element, source };
            family
                .validate_committed(
                    &self.state.element_deformation[element],
                    &self.state.element_states[element],
                    self.temperatures_k[element],
                    controls.material.yield_tolerance_pa,
                )
                .map_err(refused)?;
            laws.push(
                family
                    .at_temperature(new_temperatures_k[element])
                    .map_err(refused)?,
            );
        }
        let frozen: Vec<_> = laws
            .iter()
            .zip(&self.state.element_states)
            .map(|(law, accepted)| FrozenPoint { law, accepted })
            .collect();
        let inertia = self
            .motion
            .as_ref()
            .map(|motion| {
                BackwardEulerInertia::new(
                    &motion.mass,
                    self.state.displacement_m.as_flattened(),
                    motion.velocity_mps.as_flattened(),
                    dt_s,
                    cx,
                )
            })
            .transpose()
            .map_err(inertia_error)?;
        let (equilibrium, force) = if pressure.is_some() {
            let problem = HyperTetProblem {
                nodes_m: self.state.nodes_m,
                tetrahedra: self.state.tetrahedra,
                material: HyperelasticLawField::PerElement(&frozen),
                prescribed_m,
                nodal_forces_n,
                body_force_n_m3,
                budget: self.state.budget,
                settings,
            };
            let equilibrium = problem.solve_loaded_from(
                Some(self.state.displacement_m.as_flattened()),
                Some(&self.state.elements),
                pressure,
                None,
                inertia.as_ref(),
                cx,
            )?;
            (equilibrium, problem.external(&self.state.elements, 1.0))
        } else {
            self.state.solve_endpoint(
                HyperelasticLawField::PerElement(&frozen),
                prescribed_m,
                nodal_forces_n,
                body_force_n_m3,
                settings,
                inertia.as_ref(),
                cx,
            )?
        };
        let new_u = equilibrium.displacement_m.as_flattened();
        let motion = inertia
            .as_ref()
            .map(|value| value.motion(new_u, cx))
            .transpose()
            .map_err(inertia_error)?;
        let count = self.state.elements.len();
        let mut states = Vec::with_capacity(count);
        let mut deformation = Vec::with_capacity(count);
        let mut stress = Vec::with_capacity(count);
        let mut jacobians = Vec::with_capacity(count);
        let mut plastic_jacobians = Vec::with_capacity(count);
        let mut element_energy = Vec::with_capacity(count);
        let mut energy = FiniteThermoplasticTetEnergy::default();
        let (mut old_free, mut new_free, mut old_internal, mut new_internal) = (0.0, 0.0, 0.0, 0.0);
        let mut point_work = 0.0;
        let tolerance = controls.material.equilibrium_work_tolerance_j;
        for (element, geometry) in self.state.elements.iter().enumerate() {
            if element % 256 == 0 {
                checkpoint(cx)?;
            }
            let f = HyperTetProblem::<FrozenPoint<'_>>::deformation_gradient(geometry, new_u);
            let value = self.families[element]
                .thermal_finite_porous()
                .checked_increment(
                    &self.state.element_deformation[element],
                    &self.state.element_states[element],
                    self.temperatures_k[element],
                    &f,
                    new_temperatures_k[element],
                    controls.material.yield_tolerance_pa,
                )
                .map_err(|source| FiniteThermalPorousTetError::Material { element, source })?;
            let a = account(&value);
            let v = geometry.volume;
            energy_gate(
                "porous thermal point identity",
                v * a.balance_error_j_m3,
                tolerance,
            )?;
            energy_gate(
                "porous thermal free energy",
                v * (a.free_after_j_m3 - value.mechanical.stored_energy_j_m3),
                tolerance,
            )?;
            old_free = v.mul_add(a.free_before_j_m3, old_free);
            new_free = v.mul_add(a.free_after_j_m3, new_free);
            old_internal = v.mul_add(a.internal_before_j_m3, old_internal);
            new_internal = v.mul_add(a.internal_after_j_m3, new_internal);
            energy.free_energy_change_j = v.mul_add(
                a.free_after_j_m3 - a.free_before_j_m3,
                energy.free_energy_change_j,
            );
            energy.internal_energy_change_j = v.mul_add(
                a.internal_energy_change_j_m3(),
                energy.internal_energy_change_j,
            );
            energy.plastic_dissipation_j =
                v.mul_add(a.plastic_dissipation_j_m3, energy.plastic_dissipation_j);
            energy.reversible_heat_j = v.mul_add(a.reversible_heat_j_m3, energy.reversible_heat_j);
            energy.caloric_heat_j = v.mul_add(a.caloric_heat_j_m3(), energy.caloric_heat_j);
            energy.coefficient_storage_change_j = v.mul_add(
                a.coefficient_storage_change_j_m3,
                energy.coefficient_storage_change_j,
            );
            energy.thermoelastic_storage_change_j = v.mul_add(
                a.thermoelastic_storage_change_j_m3,
                energy.thermoelastic_storage_change_j,
            );
            energy.endpoint_work_defect_j =
                v.mul_add(a.endpoint_work_defect_j_m3, energy.endpoint_work_defect_j);
            energy.incremental_potential_j =
                v.mul_add(a.incremental_potential_j_m3, energy.incremental_potential_j);
            point_work = v.mul_add(a.internal_endpoint_work_j_m3, point_work);
            jacobians.push(value.new.jacobian);
            plastic_jacobians.push(value.new.plastic_jacobian);
            deformation.push(f);
            stress.push(value.mechanical.first_piola_pa);
            states.push(value.mechanical.state);
            element_energy.push(a);
        }
        let work = self.state.endpoint_work(
            &equilibrium,
            &force,
            &deformation,
            &stress,
            motion.as_ref().map(|(_, value)| value.inertial_work_j),
            cx,
        )?;
        energy.internal_endpoint_work_j = work.internal_j;
        energy.dead_load_endpoint_work_j = work.dead_load_j;
        energy.support_endpoint_work_j = work.support_j;
        energy.equilibrium_work_residual_j = work.residual_j;
        let pressure = pressure
            .map(
                |(surface, field)| -> Result<_, FiniteThermalPorousTetError> {
                    let old_u = self.state.displacement_m.as_flattened();
                    let pressure_force = surface.force_n(new_u, field, cx)?;
                    let mut endpoint_work = Dd::ZERO;
                    for (dof, (&force, (&after, &before))) in pressure_force
                        .iter()
                        .zip(new_u.iter().zip(old_u))
                        .enumerate()
                    {
                        if dof % 256 == 0 {
                            checkpoint(cx)?;
                        }
                        endpoint_work = endpoint_work
                            + Dd::from_f64(force) * (Dd::from_f64(after) - Dd::from_f64(before));
                    }
                    let endpoint_work_j = endpoint_work.to_f64();
                    let work_j = surface.work_j(old_u, new_u, field, cx)?;
                    let correction = work_j - endpoint_work_j;
                    energy.equilibrium_work_residual_j += endpoint_work_j;
                    energy_gate(
                        "thermal porous pressure/constitutive endpoint integration defect",
                        energy.endpoint_work_defect_j.abs() + correction.abs(),
                        controls.material.maximum_endpoint_work_defect_j,
                    )?;
                    Ok(FinitePorousPressureEnergy {
                        endpoint_work_j,
                        work_j,
                        endpoint_work_correction_j: correction,
                        external_work_defect_j: energy.endpoint_work_defect_j + correction,
                    })
                },
            )
            .transpose()?;
        energy.balance_error_j = work.dead_load_j + work.support_j
            - energy.internal_energy_change_j
            - energy.caloric_heat_j
            - energy.endpoint_work_defect_j;
        if let Some(pressure) = pressure {
            // Exact/path pressure work and its geometric correction remain a
            // separate mechanical account. Neither changes caloric heat.
            energy.balance_error_j = (Dd::from_f64(work.dead_load_j)
                + Dd::from_f64(work.support_j)
                + Dd::from_f64(pressure.work_j)
                - Dd::from_f64(energy.internal_energy_change_j)
                - Dd::from_f64(energy.caloric_heat_j)
                - Dd::from_f64(energy.endpoint_work_defect_j)
                - Dd::from_f64(pressure.endpoint_work_correction_j))
            .to_f64();
        }
        if let Some((_, value)) = &motion {
            energy_gate(
                "porous backward-Euler kinetic loss",
                value.numerical_loss_j,
                controls.maximum_numerical_kinetic_loss_j,
            )?;
            energy_gate("porous kinetic identity", value.balance_error_j, tolerance)?;
            energy.balance_error_j -= value.kinetic_change_j + value.numerical_loss_j;
        }
        for (quantity, value) in [
            (
                "accepted porous thermal free energy",
                old_free - self.state.stored_energy_j,
            ),
            (
                "accepted porous thermal internal energy",
                old_internal - self.internal_energy_j,
            ),
            (
                "porous assembled thermal potential",
                equilibrium.strain_energy_j - energy.incremental_potential_j,
            ),
            (
                "porous free energy/yield potential",
                energy.incremental_potential_j - new_free - energy.plastic_dissipation_j,
            ),
            (
                "porous thermal point/assembled work",
                point_work - work.internal_j,
            ),
            (
                "porous thermal free-energy increment",
                new_free - old_free - energy.free_energy_change_j,
            ),
            (
                "porous thermal internal-energy increment",
                new_internal - old_internal - energy.internal_energy_change_j,
            ),
            (
                "porous physical caloric decomposition",
                energy.caloric_heat_j - energy.plastic_dissipation_j - energy.reversible_heat_j,
            ),
            (
                "porous thermal equilibrium work",
                energy.equilibrium_work_residual_j,
            ),
            (
                "porous thermal complete work identity",
                energy.balance_error_j,
            ),
        ] {
            energy_gate(quantity, value, tolerance)?;
        }
        energy_gate(
            "signed porous thermal endpoint defect",
            energy.endpoint_work_defect_j,
            controls.material.maximum_endpoint_work_defect_j,
        )?;
        if !new_free.is_finite()
            || new_free < 0.0
            || !new_internal.is_finite()
            || !energy.plastic_dissipation_j.is_finite()
            || energy.plastic_dissipation_j < 0.0
            || [
                energy.coefficient_storage_change_j,
                energy.thermoelastic_storage_change_j,
                energy.reversible_heat_j,
                energy.caloric_heat_j,
            ]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite porous thermal storage/heat and nonnegative yield loss required",
            )
            .into());
        }
        checkpoint(cx)?;
        let (next_motion, inertia) = match (motion, self.motion.as_ref()) {
            (Some((velocity_mps, account)), Some(old)) => (
                Some(Motion {
                    velocity_mps,
                    reference_density_kg_m3: Arc::clone(&old.reference_density_kg_m3),
                    reference_nodal_density_kg_m3: old
                        .reference_nodal_density_kg_m3
                        .as_ref()
                        .map(Arc::clone),
                    mass: Arc::clone(&old.mass),
                }),
                Some(account),
            ),
            (None, None) => (None, None),
            _ => unreachable!("original inertia and motion share one owner"),
        };
        Ok(FiniteThermalPorousTetStep {
            solid: Self {
                state: FinitePlasticTetState {
                    nodes_m: self.state.nodes_m,
                    tetrahedra: self.state.tetrahedra,
                    elements: Arc::clone(&self.state.elements),
                    budget: self.state.budget,
                    displacement_m: equilibrium.displacement_m.clone(),
                    element_deformation: deformation,
                    element_states: states,
                    element_first_piola_pa: stress,
                    stored_energy_j: new_free,
                },
                families: Arc::clone(&self.families),
                temperatures_k: new_temperatures_k.to_vec(),
                internal_energy_j: new_internal,
                jacobians,
                plastic_jacobians,
                motion: next_motion,
            },
            equilibrium,
            element_energy,
            energy,
            inertia,
            pressure,
        })
    }
}
