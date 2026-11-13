//! Finite thermoplastic endpoints on the original tetrahedral Newton solve.
//!
//! Each source family remains owned in original element order. Temperature is
//! frozen only for an individual mechanical proposal; the accepted temperature,
//! total deformation, inelastic metric, accumulated flow and any recovering
//! kinetic resistance advance together. The complete native history type is
//! retained without reducing recovering creep to cumulative strain.
//! Old states are checked at their OLD temperature: new thermal softening may
//! legitimately produce plastic flow without any imposed displacement change.
//! Explicit Perzyna families use the actual physical interval in every frozen
//! trial; accepted positive overstress is retained for later relaxation. The
//! variational flow potential and physical yield/viscous heat remain distinct.
//!
//! Physical free/internal storage, yield dissipation, reversible caloric transfer
//! and signed mechanical endpoint defect come from the native thermal owner.
//! Only dissipation plus reversible transfer belongs in a separate caloric
//! solve. Neither the signed defect nor backward-Euler kinetic loss is heat.
//! The caller owns the thermal equation and source-specific caloric admission.
//!
//! Optional inertia uses the same cached original-reference consistent mass and
//! original momentum solve as fixed-temperature finite plasticity. This module
//! does not infer density from a source or update mass with thermal expansion.
//! It retains the parent P1 locking limitations, and supplies no evolving
//! conductance, phase change or new thermal solver. Its `contact` child
//! admits frictionless translating rigid planes through the original sustained
//! complementarity owner. Held follower
//! pressure reuses the original current-area load, tangent and path-work owner.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::MaterialError;
use fs_material::hyper::cards::FiniteStrainPotential;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2ThermalError, FiniteJ2ThermalIncrement,
    FiniteThermalInelasticMaterial, FiniteThermalInelasticPoint,
};
use fs_math::dd::Dd;

use super::super::{
    HyperTetError, HyperTetProblem, HyperTetSettings, HyperTetSolution, HyperelasticLawField,
    inertia_error,
};
use super::{
    FinitePlasticTetError, FinitePlasticTetState, IDENTITY, checkpoint, energy_gate, invalid,
    reference_size, validate_controls,
};
use crate::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use crate::inertia::{BackwardEulerInertia, P1ReferenceMass};
use crate::linear3::TetAssemblyBudget;

/// Material work/yield and numerical kinetic-loss budgets remain independent.
pub use super::dynamics::FinitePlasticMotionControls as FiniteThermoplasticControls;
/// Shared held-pressure work account; no constitutive or caloric state.
pub use super::porous::FinitePorousPressureEnergy as FiniteThermoplasticPressureWork;
/// Shared kinetic account; numerical loss is never caloric deposition.
pub use crate::inertia::InertialEnergy;

/// Finite-duration frictionless rigid tooling through the native thermal solid.
#[cfg(feature = "contact")]
pub mod contact;
/// Zero-duration momentum transfer with complete thermal and plastic memory.
pub mod impulse;
pub use crate::inertia::{
    SolidImpulseControls, SolidImpulseReport, SolidMomentum, SolidPointKinematics,
    SolidPointMobility, SolidPointMobilityMatrix, TetImpulsePoint,
};
pub use impulse::FiniteThermoplasticTetImpulse;

/// Compact original-reference-volume account; no tangent/history duplication.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteThermoplasticElementEnergy {
    pub free_before_j_m3: f64,
    pub free_after_j_m3: f64,
    pub internal_before_j_m3: f64,
    pub internal_after_j_m3: f64,
    /// Native primary loss: plastic work for J2/combined laws, or creep work
    /// for a pure Norton family (including its recovering variant).
    pub plastic_dissipation_j_m3: f64,
    /// Separate Norton loss for the simultaneous creep/plastic selection.
    pub creep_dissipation_j_m3: f64,
    pub reversible_heat_j_m3: f64,
    pub coefficient_storage_change_j_m3: f64,
    /// Delta(T*s), not an additional caloric source.
    pub thermoelastic_storage_change_j_m3: f64,
    /// Signed mechanical-coordinate work quadrature error, never heat.
    pub endpoint_work_defect_j_m3: f64,
    pub incremental_potential_j_m3: f64,
    pub internal_endpoint_work_j_m3: f64,
    pub mechanical_coordinate_work_j_m3: f64,
    pub thermal_deformation_work_j_m3: f64,
    pub balance_error_j_m3: f64,
}

impl FiniteThermoplasticElementEnergy {
    /// Physical inelastic loss, excluding reversible and numerical accounts.
    pub fn inelastic_dissipation_j_m3(&self) -> f64 {
        self.plastic_dissipation_j_m3 + self.creep_dissipation_j_m3
    }
    /// Deposit this signed quantity once into the original-volume caloric owner.
    pub fn caloric_heat_j_m3(&self) -> f64 {
        self.inelastic_dissipation_j_m3() + self.reversible_heat_j_m3
    }
    /// Stable mechanical internal-energy change from the native entropy account.
    pub fn internal_energy_change_j_m3(&self) -> f64 {
        self.free_after_j_m3 - self.free_before_j_m3 + self.thermoelastic_storage_change_j_m3
    }

    pub(crate) fn from_increment<S>(value: &FiniteJ2ThermalIncrement<S>) -> Self {
        Self {
            free_before_j_m3: value.old.free_energy_j_m3,
            free_after_j_m3: value.new.free_energy_j_m3,
            internal_before_j_m3: value.old.internal_energy_j_m3,
            internal_after_j_m3: value.new.internal_energy_j_m3,
            plastic_dissipation_j_m3: value.mechanical.plastic_dissipation_j_m3,
            creep_dissipation_j_m3: value.mechanical.creep_dissipation_j_m3,
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
}

/// Whole-original-mesh energies. Mechanical internal storage excludes caloric h.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FiniteThermoplasticTetEnergy {
    pub free_energy_change_j: f64,
    pub internal_energy_change_j: f64,
    /// Native primary loss; a pure Norton family's sole inelastic loss is
    /// creep work. Combined laws keep plastic and creep work separate.
    pub plastic_dissipation_j: f64,
    /// Separate Norton loss for the simultaneous creep/plastic selection.
    pub creep_dissipation_j: f64,
    pub reversible_heat_j: f64,
    pub caloric_heat_j: f64,
    pub coefficient_storage_change_j: f64,
    pub thermoelastic_storage_change_j: f64,
    /// Signed mechanical-coordinate quadrature defect; not physical dissipation.
    pub endpoint_work_defect_j: f64,
    pub incremental_potential_j: f64,
    pub internal_endpoint_work_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    /// External minus material and, when enabled, inertial endpoint work.
    pub equilibrium_work_residual_j: f64,
    /// External work minus ΔU, caloric transfer, signed defect and optional
    /// kinetic change/numerical loss. Independently checked before acceptance.
    pub balance_error_j: f64,
}

impl FiniteThermoplasticTetEnergy {
    /// Physical inelastic loss over the original mesh, with no numerical loss.
    pub fn inelastic_dissipation_j(&self) -> f64 {
        self.plastic_dissipation_j + self.creep_dissipation_j
    }
}

#[derive(Debug)]
struct ThermalMotion {
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    reference_nodal_density_kg_m3: Option<Arc<[f64]>>,
    mass: Arc<P1ReferenceMass>,
}

/// Original mesh/source families and an immutable complete thermal-plastic state.
#[derive(Debug)]
pub struct FiniteThermoplasticTetSolid<
    'a,
    L: FiniteThermalInelasticMaterial = AffineThermalFiniteJ2,
> {
    state: FinitePlasticTetState<'a, L::State>,
    families: Arc<[L]>,
    temperatures_k: Vec<f64>,
    internal_energy_j: f64,
    motion: Option<ThermalMotion>,
}

/// Fully checked proposal; no part of the borrowed old solid has changed.
#[derive(Debug)]
pub struct FiniteThermoplasticTetStep<'a, L: FiniteThermalInelasticMaterial = AffineThermalFiniteJ2>
{
    pub solid: FiniteThermoplasticTetSolid<'a, L>,
    /// Native forces/reactions and convergence; strain_energy_j is potential.
    pub equilibrium: HyperTetSolution,
    /// Compact accounts in original element order, per reference volume.
    pub element_energy: Vec<FiniteThermoplasticElementEnergy>,
    pub energy: FiniteThermoplasticTetEnergy,
    pub inertia: Option<InertialEnergy>,
    /// Current-area pressure path work and its signed endpoint correction.
    /// The material endpoint defect and physical caloric heat stay separate.
    pub pressure: Option<FiniteThermoplasticPressureWork>,
}

/// Preserve the original spatial and thermal material owners' diagnoses.
#[derive(Clone, Debug, PartialEq)]
pub enum FiniteThermoplasticTetError {
    Solid(FinitePlasticTetError),
    Material {
        element: usize,
        source: FiniteJ2ThermalError,
    },
}
impl From<FinitePlasticTetError> for FiniteThermoplasticTetError {
    fn from(error: FinitePlasticTetError) -> Self {
        Self::Solid(error)
    }
}
impl From<HyperTetError> for FiniteThermoplasticTetError {
    fn from(error: HyperTetError) -> Self {
        Self::Solid(error.into())
    }
}
impl core::fmt::Display for FiniteThermoplasticTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Solid(error) => error.fmt(f),
            Self::Material { element, source } => {
                write!(f, "finite thermal tetrahedron {element}: {source}")
            }
        }
    }
}
impl std::error::Error for FiniteThermoplasticTetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Solid(error) => Some(error),
            Self::Material { source, .. } => Some(source),
        }
    }
}

struct FrozenThermalPoint<'a, L: FiniteThermalInelasticPoint> {
    law: &'a L,
    accepted: &'a L::State,
}
fn material_error(source: FiniteJ2ThermalError) -> MaterialError {
    MaterialError::State {
        what: source.to_string(),
    }
}
impl<L: FiniteThermalInelasticPoint> FiniteStrainPotential for FrozenThermalPoint<'_, L> {
    fn energy_density(&self, f: &[f64; 9]) -> f64 {
        self.law
            .stress_and_potential(f, self.accepted)
            .map_or(f64::NAN, |(_, value)| value)
    }
    fn energy_density_increment(&self, before: &[f64; 9], after: &[f64; 9]) -> f64 {
        super::frozen_potential_increment(before, after, |f| {
            self.law.stress_and_potential(f, self.accepted).ok()
        })
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

impl<'a, L: FiniteThermalInelasticMaterial> FiniteThermoplasticTetSolid<'a, L> {
    /// Admit the virgin reference state and own one complete family per element.
    /// Initial temperatures MUST equal each family's stress-free reference.
    /// Other initial states require an actual accepted loading/thermal path.
    ///
    /// # Errors
    /// Count/reference-temperature, native mesh/budget or source-domain refusal.
    pub fn try_new(
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        families: Vec<L>,
        initial_temperatures_k: &[f64],
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticTetError> {
        checkpoint(cx)?;
        reference_size(nodes_m.len(), tetrahedra.len(), budget)?;
        if families.len() != tetrahedra.len() || initial_temperatures_k.len() != tetrahedra.len() {
            return Err(invalid(
                "one original family and initial temperature per tetrahedron required",
            )
            .into());
        }
        let mut laws = Vec::with_capacity(families.len());
        let mut states = Vec::with_capacity(families.len());
        for (element, family) in families.iter().enumerate() {
            if element % 1024 == 0 {
                checkpoint(cx)?;
            }
            let temperature = initial_temperatures_k[element];
            if temperature != family.reference_temperature_k() {
                return Err(invalid(
                    "virgin finite thermal initial temperature must equal its source reference",
                )
                .into());
            }
            laws.push(
                // Reference preparation only checks geometry and material
                // cardinality. This inspection view cannot advance kinetics
                // and introduces no invented physical interval.
                family
                    .at_temperature_for_inspection(temperature)
                    .map_err(|source| FiniteThermoplasticTetError::Material { element, source })?,
            );
            states.push(family.initial_state());
        }
        let frozen: Vec<_> = laws
            .iter()
            .zip(&states)
            .map(|(law, accepted)| FrozenThermalPoint { law, accepted })
            .collect();
        let elements = FinitePlasticTetState::reference_elements(
            nodes_m,
            tetrahedra,
            HyperelasticLawField::PerElement(&frozen),
            budget,
            cx,
        )?;
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
            motion: None,
        })
    }

    /// Bind explicit original-reference mass once, retaining the complete state.
    /// Density has one value per original element; velocity one per original node.
    /// This inventory alone does not imply a source or caloric density claim.
    ///
    /// # Errors
    /// Rebinding mass, invalid velocity/density, mass budget or cancellation.
    pub fn with_inertia(
        mut self,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticTetError> {
        self.admit_initial_motion(velocity_mps, cx)?;
        let mass = P1ReferenceMass::new(
            self.state.elements.iter().map(|e| (e.nodes, e.volume)),
            reference_density_kg_m3,
            self.state.displacement_m.as_flattened().len(),
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        self.motion = Some(ThermalMotion {
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(reference_density_kg_m3),
            reference_nodal_density_kg_m3: None,
            mass: Arc::new(mass),
        });
        Ok(self)
    }

    /// Bind positive P1 density on the original nodes to the native momentum solve.
    ///
    /// The shared consistent mass integrates `rho_0*N_a*N_b` exactly. Cell
    /// means describe inventory; they never replace the graded mass operator.
    /// All accepted deformation, temperature and plastic/creep history survives
    /// binding. Subsequent heating and expansion retain this original mass and
    /// nodal field. This authored inventory implies no new caloric/source claim.
    ///
    /// # Errors
    /// Rebinding, invalid density/velocity, unrepresentable mass or cancellation.
    pub fn with_nodal_inertia(
        mut self,
        reference_nodal_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteThermoplasticTetError> {
        self.admit_initial_motion(velocity_mps, cx)?;
        let (mass, means) = P1ReferenceMass::from_nodal_density(
            self.state.elements.iter().map(|e| (e.nodes, e.volume)),
            reference_nodal_density_kg_m3,
            self.state.displacement_m.as_flattened().len(),
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        self.motion = Some(ThermalMotion {
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(means),
            reference_nodal_density_kg_m3: Some(Arc::from(reference_nodal_density_kg_m3)),
            mass: Arc::new(mass),
        });
        Ok(self)
    }

    fn admit_initial_motion(
        &self,
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<(), FiniteThermoplasticTetError> {
        checkpoint(cx)?;
        if self.motion.is_some() {
            return Err(invalid("finite thermal reference mass is already bound").into());
        }
        if velocity_mps.len() != self.state.nodes_m.len()
            || velocity_mps.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid("one finite initial velocity per original node required").into());
        }
        Ok(())
    }

    /// Original immutable reference coordinates [m], in native node order.
    pub fn reference_nodes_m(&self) -> &'a [[f64; 3]] {
        self.state.nodes_m
    }
    /// Original ordered material connectivity; deformation never rewrites it.
    pub fn reference_tetrahedra(&self) -> &'a [[usize; 4]] {
        self.state.tetrahedra
    }
    /// Original size and geometric-quality admission bounds.
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
    pub fn element_states(&self) -> &[L::State] {
        &self.state.element_states
    }
    pub fn element_first_piola_pa(&self) -> &[[f64; 9]] {
        &self.state.element_first_piola_pa
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
        self.motion.as_ref().map(|m| m.velocity_mps.as_slice())
    }
    /// Original cell densities or cell means of the retained P1 nodal field.
    /// Means alone do not reproduce a graded body's consistent inertia.
    pub fn reference_density_kg_m3(&self) -> Option<&[f64]> {
        self.motion
            .as_ref()
            .map(|m| m.reference_density_kg_m3.as_ref())
    }
    /// Complete original nodal density, when inertia was bound with that field.
    /// None denotes cellwise inertia or an unbound quasistatic state.
    pub fn reference_nodal_density_kg_m3(&self) -> Option<&[f64]> {
        self.motion
            .as_ref()
            .and_then(|m| m.reference_nodal_density_kg_m3.as_deref())
    }
    pub fn reference_mass_kg(&self) -> Option<f64> {
        self.motion.as_ref().map(|m| m.mass.reference_mass_kg())
    }

    /// Project an original-node P1 acceleration through the retained mass, M*a.
    ///
    /// The returned full nodal loads [N] include supported nodes and preserve
    /// the density-weighted force and moment, including gravity on graded
    /// material. Supply these held loads to `advance` exactly once. This does
    /// not resample density at the new temperature or deformed configuration.
    /// Evolving acceleration, potential and tangent belong to the field owner.
    ///
    /// # Errors
    /// Unbound inertia, invalid shape, nonfinite acceleration/force or cancellation.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, FiniteThermoplasticTetError> {
        checkpoint(cx)?;
        let motion = self.motion.as_ref().ok_or_else(|| {
            invalid("finite thermal acceleration requires already bound reference inertia")
        })?;
        if acceleration_m_s2.len() != self.state.nodes_m.len()
            || acceleration_m_s2.iter().flatten().any(|a| !a.is_finite())
        {
            return Err(
                invalid("one finite acceleration vector per original node required").into(),
            );
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

    /// Propose one endpoint at supplied element temperatures and held dead loads.
    /// Every Newton trial sees identical old plastic history and endpoint T.
    /// Thermal, mechanical and optional momentum state are accepted atomically.
    /// `dt_s` is a finite positive physical interval even in quasistatic mode.
    ///
    /// # Errors
    /// Temperature/domain, native solve, work/yield/kinetic gate or cancellation.
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
    ) -> Result<FiniteThermoplasticTetStep<'a, L>, FiniteThermoplasticTetError> {
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

    /// Advance the original thermal plastic history under a held current-area
    /// pressure field. Native pressure forces and follower tangents participate
    /// in every Newton trial; original material and motion histories stay frozen.
    ///
    /// Uniform pressure uses exact -p delta(V) work. Other pressure variants
    /// retain their original conservative or declared-path work conventions.
    /// The sum of absolute constitutive and pressure integration defects must
    /// fit the material work budget. Neither correction becomes physical heat.
    ///
    /// # Errors
    /// Foreign reference surface, invalid pressure, or any original thermal,
    /// mechanical, inertia, energy-budget or cancellation refusal.
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
    ) -> Result<FiniteThermoplasticTetStep<'a, L>, FiniteThermoplasticTetError> {
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
    ) -> Result<FiniteThermoplasticTetStep<'a, L>, FiniteThermoplasticTetError> {
        checkpoint(cx)?;
        validate_controls(settings, controls.material)?;
        if let Some((surface, field)) = pressure {
            if !surface.matches_reference(self.state.nodes_m, self.state.tetrahedra) {
                return Err(invalid(
                    "thermal plastic pressure surface must bind the original mesh",
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
            return Err(invalid("finite positive interval, finite nonnegative kinetic budget and one temperature per element required").into());
        }
        let mut laws = Vec::with_capacity(self.families.len());
        for (element, family) in self.families.iter().enumerate() {
            if element % 1024 == 0 {
                checkpoint(cx)?;
            }
            let refused = |source| FiniteThermoplasticTetError::Material { element, source };
            family
                .at_temperature_for_inspection(self.temperatures_k[element])
                .map_err(refused)?
                .validate_committed(
                    &self.state.element_deformation[element],
                    &self.state.element_states[element],
                    controls.material.yield_tolerance_pa,
                )
                .map_err(refused)?;
            laws.push(
                family
                    .at_temperature_for_step(new_temperatures_k[element], dt_s)
                    .map_err(refused)?,
            );
        }
        let frozen: Vec<_> = laws
            .iter()
            .zip(&self.state.element_states)
            .map(|(law, accepted)| FrozenThermalPoint { law, accepted })
            .collect();
        let inertia = self
            .motion
            .as_ref()
            .map(|m| {
                BackwardEulerInertia::new(
                    &m.mass,
                    self.state.displacement_m.as_flattened(),
                    m.velocity_mps.as_flattened(),
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
            .map(|m| m.motion(new_u, cx))
            .transpose()
            .map_err(inertia_error)?;
        let count = self.state.elements.len();
        let mut states = Vec::with_capacity(count);
        let mut deformation = Vec::with_capacity(count);
        let mut stress = Vec::with_capacity(count);
        let mut element_energy = Vec::with_capacity(count);
        let mut energy = FiniteThermoplasticTetEnergy::default();
        let (mut old_free, mut new_free, mut old_internal, mut new_internal) = (0.0, 0.0, 0.0, 0.0);
        let mut point_work = 0.0;
        let mut flow_potential_j = 0.0;
        let tolerance = controls.material.equilibrium_work_tolerance_j;
        for (element, geometry) in self.state.elements.iter().enumerate() {
            if element % 1024 == 0 {
                checkpoint(cx)?;
            }
            let f = HyperTetProblem::<FrozenThermalPoint<'_, L::Frozen>>::deformation_gradient(
                geometry, new_u,
            );
            let value = self.families[element]
                .checked_increment_for_step(
                    &self.state.element_deformation[element],
                    &self.state.element_states[element],
                    self.temperatures_k[element],
                    &f,
                    new_temperatures_k[element],
                    controls.material.yield_tolerance_pa,
                    dt_s,
                )
                .map_err(|source| FiniteThermoplasticTetError::Material { element, source })?;
            let account = FiniteThermoplasticElementEnergy::from_increment(&value);
            let v = geometry.volume;
            // The native law integrates its inverse kinetics, including any
            // retained primary-creep resistance. This potential is distinct
            // from the full physical dissipation deposited as caloric heat.
            flow_potential_j = v.mul_add(value.mechanical.flow_potential_j_m3, flow_potential_j);
            energy_gate(
                "finite thermal point energy identity",
                v * account.balance_error_j_m3,
                tolerance,
            )?;
            energy_gate(
                "finite thermal constitutive free energy",
                v * (account.free_after_j_m3 - value.mechanical.stored_energy_j_m3),
                tolerance,
            )?;
            old_free = v.mul_add(account.free_before_j_m3, old_free);
            new_free = v.mul_add(account.free_after_j_m3, new_free);
            old_internal = v.mul_add(account.internal_before_j_m3, old_internal);
            new_internal = v.mul_add(account.internal_after_j_m3, new_internal);
            energy.free_energy_change_j = v.mul_add(
                account.free_after_j_m3 - account.free_before_j_m3,
                energy.free_energy_change_j,
            );
            energy.internal_energy_change_j = v.mul_add(
                account.internal_energy_change_j_m3(),
                energy.internal_energy_change_j,
            );
            energy.plastic_dissipation_j = v.mul_add(
                account.plastic_dissipation_j_m3,
                energy.plastic_dissipation_j,
            );
            energy.creep_dissipation_j =
                v.mul_add(account.creep_dissipation_j_m3, energy.creep_dissipation_j);
            energy.reversible_heat_j =
                v.mul_add(account.reversible_heat_j_m3, energy.reversible_heat_j);
            energy.caloric_heat_j = v.mul_add(account.caloric_heat_j_m3(), energy.caloric_heat_j);
            energy.coefficient_storage_change_j = v.mul_add(
                account.coefficient_storage_change_j_m3,
                energy.coefficient_storage_change_j,
            );
            energy.thermoelastic_storage_change_j = v.mul_add(
                account.thermoelastic_storage_change_j_m3,
                energy.thermoelastic_storage_change_j,
            );
            energy.endpoint_work_defect_j = v.mul_add(
                account.endpoint_work_defect_j_m3,
                energy.endpoint_work_defect_j,
            );
            energy.incremental_potential_j = v.mul_add(
                account.incremental_potential_j_m3,
                energy.incremental_potential_j,
            );
            point_work = v.mul_add(account.internal_endpoint_work_j_m3, point_work);
            deformation.push(f);
            stress.push(value.mechanical.first_piola_pa);
            states.push(value.mechanical.state);
            element_energy.push(account);
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
                |(surface, field)| -> Result<_, FiniteThermoplasticTetError> {
                    let old_u = self.state.displacement_m.as_flattened();
                    let force = surface.force_n(new_u, field, cx)?;
                    let mut endpoint_work = Dd::ZERO;
                    for (dof, (&value, (&after, &before))) in
                        force.iter().zip(new_u.iter().zip(old_u)).enumerate()
                    {
                        if dof % 256 == 0 {
                            checkpoint(cx)?;
                        }
                        endpoint_work = endpoint_work
                            + Dd::from_f64(value) * (Dd::from_f64(after) - Dd::from_f64(before));
                    }
                    let endpoint_work_j = endpoint_work.to_f64();
                    let work_j = surface.work_j(old_u, new_u, field, cx)?;
                    let correction = work_j - endpoint_work_j;
                    energy.equilibrium_work_residual_j += endpoint_work_j;
                    energy_gate(
                        "thermal plastic pressure/constitutive endpoint integration defect",
                        energy.endpoint_work_defect_j.abs() + correction.abs(),
                        controls.material.maximum_endpoint_work_defect_j,
                    )?;
                    Ok(FiniteThermoplasticPressureWork {
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
                "backward-Euler numerical kinetic loss",
                value.numerical_loss_j,
                controls.maximum_numerical_kinetic_loss_j,
            )?;
            energy_gate(
                "finite thermal kinetic identity",
                value.balance_error_j,
                tolerance,
            )?;
            energy.balance_error_j -= value.kinetic_change_j + value.numerical_loss_j;
        }
        for (quantity, value) in [
            (
                "accepted finite thermal free energy",
                old_free - self.state.stored_energy_j,
            ),
            (
                "accepted finite thermal internal energy",
                old_internal - self.internal_energy_j,
            ),
            (
                "assembled finite thermal potential",
                equilibrium.strain_energy_j - energy.incremental_potential_j,
            ),
            (
                "finite thermal stored energy/flow potential",
                energy.incremental_potential_j - new_free - flow_potential_j,
            ),
            (
                "finite thermal point/assembled work",
                point_work - work.internal_j,
            ),
            (
                "finite thermal free-energy increment",
                new_free - old_free - energy.free_energy_change_j,
            ),
            (
                "finite thermal internal-energy increment",
                new_internal - old_internal - energy.internal_energy_change_j,
            ),
            (
                "finite thermal caloric decomposition",
                energy.caloric_heat_j - energy.inelastic_dissipation_j() - energy.reversible_heat_j,
            ),
            (
                "finite thermal support/dead-load work",
                energy.equilibrium_work_residual_j,
            ),
            (
                "finite thermal total endpoint identity",
                energy.balance_error_j,
            ),
        ] {
            energy_gate(quantity, value, tolerance)?;
        }
        energy_gate(
            "signed thermal mechanical endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.material.maximum_endpoint_work_defect_j,
        )?;
        if !new_free.is_finite()
            || new_free < 0.0
            || !new_internal.is_finite()
            || !energy.plastic_dissipation_j.is_finite()
            || energy.plastic_dissipation_j < 0.0
            || !energy.creep_dissipation_j.is_finite()
            || energy.creep_dissipation_j < 0.0
            || [
                energy.coefficient_storage_change_j,
                energy.thermoelastic_storage_change_j,
                energy.reversible_heat_j,
                energy.caloric_heat_j,
            ]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(invalid("finite storage/internal energy and nonnegative physical plastic and creep dissipation required").into());
        }
        checkpoint(cx)?;
        let (next_motion, inertia) = match (motion, self.motion.as_ref()) {
            (Some((velocity_mps, account)), Some(old)) => (
                Some(ThermalMotion {
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
            _ => unreachable!("optional native inertia and accepted motion share one owner"),
        };
        let solid = Self {
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
            motion: next_motion,
        };
        Ok(FiniteThermoplasticTetStep {
            solid,
            equilibrium,
            element_energy,
            energy,
            inertia,
            pressure,
        })
    }
}
