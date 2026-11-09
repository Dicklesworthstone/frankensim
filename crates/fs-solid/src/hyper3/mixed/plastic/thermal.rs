//! Evolving temperature and complete finite inelastic history on P2/P1 solids.
//!
//! Each original element owns its authored thermal family and supplied endpoint
//! temperature. Old history is admitted at its OLD temperature and pressure;
//! Newton trials borrow the same complete history, new temperature and physical
//! interval. The original mixed solver retains pressure as an independent field.
//! Elastic incompressibility constrains the thermally expanded volume weakly,
//! rather than suppressing thermal expansion.
//!
//! Native original-volume free energy, pressure-dependent entropy, mechanical
//! internal energy, physical inelastic heat and reversible caloric exchange are
//! retained separately. Signed endpoint errors and backward-Euler kinetic loss
//! never become caloric heat. The per-element heat receipt can be consumed once
//! by an external thermal owner. Supplying temperatures does not solve heat
//! transport or establish a coupled heat-capacity/passivity condition.
//!
//! Optional inertia uses the same cached consistent P2 reference mass as the
//! isothermal mixed owner. Heating and deformation never resample density.
//! Complete successors publish only after material, equilibrium and work gates;
//! cancellation and refusal leave the accepted object available for retry.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::mixed::VolumetricTreatment;
use fs_material::plastic::finite::thermal::mixed::{
    MixedThermalInelasticPoint, checked_mixed_increment,
};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2ThermalError, FiniteThermalInelasticMaterial,
};

use super::super::element::PreparedMixedMesh;
use super::super::{
    MixedLawField, MixedUpError, MixedUpProblem, MixedUpSettings, MixedUpSolution, checkpoint,
    inertia_error, invalid, point_deformation,
};
use crate::inertia::{BackwardEulerInertia, ReferenceMass};
use crate::linear3::TetAssemblyBudget;

/// Independent material-work and numerical kinetic-loss allowances.
pub use crate::hyper3::plastic::dynamics::FinitePlasticMotionControls as MixedThermoplasticControls;
/// Native original-mass kinetic accounts, distinct from physical material heat.
pub use crate::inertia::InertialEnergy;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// Whole-original-mesh mechanical and caloric accounts, all in joules.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MixedThermoplasticTetEnergy {
    pub free_energy_change_j: f64,
    /// Mechanical internal storage includes T*s and excludes caloric storage.
    pub internal_energy_change_j: f64,
    /// Native primary loss; a standalone Norton law reports its creep here.
    pub plastic_dissipation_j: f64,
    /// Separate Norton heat for simultaneous creep and plasticity.
    pub creep_dissipation_j: f64,
    /// Signed reversible exchange with the caloric owner.
    pub reversible_heat_j: f64,
    /// Physical inelastic heat plus reversible exchange; deposit exactly once.
    pub caloric_heat_j: f64,
    pub coefficient_storage_change_j: f64,
    /// Delta(T*s), already contained in the internal-energy change.
    pub thermoelastic_storage_change_j: f64,
    /// Signed material/volume endpoint quadrature error; never physical heat.
    pub endpoint_work_defect_j: f64,
    /// Sum of absolute point defects; prevents cancellation in work admission.
    pub absolute_endpoint_work_defect_j: f64,
    /// Complete mixed saddle potential, including the pressure dual.
    pub incremental_potential_j: f64,
    /// Native variational flow potential, distinct from physical heat.
    pub flow_potential_j: f64,
    pub internal_endpoint_work_j: f64,
    pub thermal_deformation_work_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    /// Zero when inertia has not been bound.
    pub inertial_endpoint_work_j: f64,
    /// Integral of p times its weak volume residual; not physical heat.
    pub pressure_constraint_work_j: f64,
    /// External minus material and optional inertial endpoint work.
    pub equilibrium_work_residual_j: f64,
    /// External work minus internal change, caloric exchange, signed error,
    /// kinetic change and numerical kinetic loss.
    pub balance_error_j: f64,
}

impl MixedThermoplasticTetEnergy {
    /// Physical inelastic heat without reversible or numerical contributions.
    pub fn inelastic_dissipation_j(&self) -> f64 {
        self.plastic_dissipation_j + self.creep_dissipation_j
    }
}

#[derive(Debug)]
struct ThermalMotion {
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    mass: Arc<ReferenceMass>,
}

/// Immutable accepted mesh, source families, temperature, pressure and memory.
#[derive(Debug)]
pub struct MixedThermoplasticTetSolid<'a, L: FiniteThermalInelasticMaterial = AffineThermalFiniteJ2>
{
    nodes_m: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    families: Arc<[L]>,
    temperatures_k: Vec<f64>,
    volumetric: VolumetricTreatment,
    mesh: Arc<PreparedMixedMesh>,
    budget: TetAssemblyBudget,
    quadrature_points_per_dir: usize,
    maximum_quadrature_points: usize,
    displacement_m: Vec<[f64; 3]>,
    pressure_pa: Vec<f64>,
    states: Vec<L::State>,
    deformation: Vec<[f64; 9]>,
    first_piola_pa: Vec<[f64; 9]>,
    free_energy_j: f64,
    internal_energy_j: f64,
    entropy_j_k: f64,
    motion: Option<ThermalMotion>,
}

/// Complete candidate endpoint and its original-volume transfer receipts.
#[derive(Debug)]
pub struct MixedThermoplasticTetStep<'a, L: FiniteThermalInelasticMaterial = AffineThermalFiniteJ2>
{
    pub solid: MixedThermoplasticTetSolid<'a, L>,
    pub equilibrium: MixedUpSolution,
    pub energy: MixedThermoplasticTetEnergy,
    pub inertia: Option<InertialEnergy>,
    /// Integrated physical heat for each original element, in joules. Includes
    /// signed reversible exchange; excludes both numerical-loss channels.
    pub element_caloric_heat_j: Vec<f64>,
}

/// Preserve native geometry, constitutive, work-budget and cancellation refusals.
#[derive(Clone, Debug, PartialEq)]
pub enum MixedThermoplasticTetError {
    Native(MixedUpError),
    Material {
        element: usize,
        point: usize,
        source: FiniteJ2ThermalError,
    },
    EnergyGate {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}

impl From<MixedUpError> for MixedThermoplasticTetError {
    fn from(error: MixedUpError) -> Self {
        Self::Native(error)
    }
}

impl core::fmt::Display for MixedThermoplasticTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => error.fmt(f),
            Self::Material {
                element,
                point,
                source,
            } => {
                write!(
                    f,
                    "mixed thermal element {element}, point {point}: {source}"
                )
            }
            Self::EnergyGate {
                quantity,
                value_j,
                maximum_abs_j,
            } => write!(
                f,
                "mixed thermal {quantity}: {value_j:e} J exceeds absolute allowance {maximum_abs_j:e} J"
            ),
        }
    }
}

impl core::error::Error for MixedThermoplasticTetError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Native(error) => Some(error),
            Self::Material { source, .. } => Some(source),
            Self::EnergyGate { .. } => None,
        }
    }
}

fn gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), MixedThermoplasticTetError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(MixedThermoplasticTetError::EnergyGate {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

impl<'a, L> MixedThermoplasticTetSolid<'a, L>
where
    L: FiniteThermalInelasticMaterial,
    L::Frozen: MixedThermalInelasticPoint,
{
    /// Admit a virgin reference with one original family and temperature per
    /// element. Initial temperatures must equal each family's reference; an
    /// already heated or loaded body requires an actual accepted path.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        nodes_m: &'a [[f64; 3]],
        elements: &'a [[usize; 10]],
        families: Vec<L>,
        initial_temperatures_k: &[f64],
        volumetric: VolumetricTreatment,
        budget: TetAssemblyBudget,
        quadrature_points_per_dir: usize,
        maximum_quadrature_points: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, MixedThermoplasticTetError> {
        checkpoint(cx)?;
        if families.len() != elements.len() || initial_temperatures_k.len() != elements.len() {
            return Err(invalid(
                "one original thermal family and temperature per element required",
            )
            .into());
        }
        let mesh = Arc::new(
            PreparedMixedMesh::new(
                nodes_m,
                elements,
                budget,
                quadrature_points_per_dir,
                maximum_quadrature_points,
                cx,
            )
            .map_err(MixedUpError::Geometry)?,
        );
        let count = elements
            .len()
            .checked_mul(mesh.points_per_element)
            .ok_or_else(|| invalid("mixed thermal history count overflow"))?;
        let mut states = Vec::with_capacity(count);
        let mut first_piola_pa = Vec::with_capacity(count);
        let (mut free_energy_j, mut internal_energy_j, mut entropy_j_k) = (0.0, 0.0, 0.0);
        for (e, element) in mesh.elements.iter().enumerate() {
            checkpoint(cx)?;
            let family = &families[e];
            let temperature = initial_temperatures_k[e];
            if temperature != family.reference_temperature_k() {
                return Err(invalid(
                    "virgin mixed thermal temperature must equal the original family reference",
                )
                .into());
            }
            let law = family
                .at_temperature_for_inspection(temperature)
                .map_err(|source| MixedThermoplasticTetError::Material {
                    element: e,
                    point: 0,
                    source,
                })?;
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let refused = |source| MixedThermoplasticTetError::Material {
                    element: e,
                    point: q,
                    source,
                };
                let state = family.initial_state();
                law.validate_mixed_committed(&IDENTITY, &state, 0.0, volumetric, 0.0)
                    .map_err(refused)?;
                let account = law
                    .mixed_thermodynamics(&IDENTITY, &state, 0.0, volumetric)
                    .map_err(refused)?;
                free_energy_j += point.weight_m3 * account.free_energy_j_m3;
                internal_energy_j += point.weight_m3 * account.internal_energy_j_m3;
                entropy_j_k += point.weight_m3 * account.entropy_j_m3_k;
                first_piola_pa.push(account.first_piola_pa);
                states.push(state);
            }
        }
        if ![free_energy_j, internal_energy_j, entropy_j_k]
            .iter()
            .all(|x| x.is_finite())
        {
            return Err(invalid("mixed thermal initial storage must be representable").into());
        }
        checkpoint(cx)?;
        Ok(Self {
            nodes_m,
            elements,
            families: Arc::from(families),
            temperatures_k: initial_temperatures_k.to_vec(),
            volumetric,
            pressure_pa: vec![0.0; mesh.pressure_nodes.len()],
            mesh,
            budget,
            quadrature_points_per_dir,
            maximum_quadrature_points,
            displacement_m: vec![[0.0; 3]; nodes_m.len()],
            states,
            deformation: vec![IDENTITY; count],
            first_piola_pa,
            free_energy_j,
            internal_energy_j,
            entropy_j_k,
            motion: None,
        })
    }

    /// Bind immutable original element densities and nodal velocities once.
    /// Signed consistent P2 mass uses the original cached reference quadrature;
    /// neither temperature nor current volume changes this inventory.
    pub fn with_inertia(
        mut self,
        reference_density_kg_m3: &[f64],
        velocity_mps: &[[f64; 3]],
        maximum_mass_assembly_bytes: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, MixedThermoplasticTetError> {
        checkpoint(cx)?;
        if self.motion.is_some() {
            return Err(invalid("mixed thermal original mass is already bound").into());
        }
        if velocity_mps.len() != self.nodes_m.len()
            || velocity_mps.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(invalid("one finite velocity per original quadratic node required").into());
        }
        let mass = ReferenceMass::from_tet10_quadrature(
            self.mesh.elements.iter().map(|e| {
                (
                    e.connectivity,
                    e.points.iter().map(|q| (q.shapes, q.weight_m3)),
                )
            }),
            reference_density_kg_m3,
            self.displacement_m.as_flattened().len(),
            maximum_mass_assembly_bytes,
            cx,
        )
        .map_err(inertia_error)?;
        checkpoint(cx)?;
        self.motion = Some(ThermalMotion {
            velocity_mps: velocity_mps.to_vec(),
            reference_density_kg_m3: Arc::from(reference_density_kg_m3),
            mass: Arc::new(mass),
        });
        Ok(self)
    }

    pub fn reference_nodes_m(&self) -> &'a [[f64; 3]] {
        self.nodes_m
    }
    pub fn reference_elements(&self) -> &'a [[usize; 10]] {
        self.elements
    }
    pub fn families(&self) -> &[L] {
        &self.families
    }
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures_k
    }
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.displacement_m
    }
    pub fn pressure_nodes(&self) -> &[usize] {
        &self.mesh.pressure_nodes
    }
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
    pub fn quadrature_states(&self) -> &[L::State] {
        &self.states
    }
    pub fn deformation_gradients(&self) -> &[[f64; 9]] {
        &self.deformation
    }
    pub fn first_piola_pa(&self) -> &[[f64; 9]] {
        &self.first_piola_pa
    }
    pub fn quadrature_points_per_element(&self) -> usize {
        self.mesh.points_per_element
    }
    pub fn reference_volume_m3(&self) -> f64 {
        self.mesh.total_volume_m3
    }
    pub const fn volumetric_treatment(&self) -> VolumetricTreatment {
        self.volumetric
    }
    pub const fn free_energy_j(&self) -> f64 {
        self.free_energy_j
    }
    pub const fn internal_energy_j(&self) -> f64 {
        self.internal_energy_j
    }
    pub const fn entropy_j_k(&self) -> f64 {
        self.entropy_j_k
    }
    pub fn velocity_mps(&self) -> Option<&[[f64; 3]]> {
        self.motion.as_ref().map(|m| m.velocity_mps.as_slice())
    }
    pub fn reference_density_kg_m3(&self) -> Option<&[f64]> {
        self.motion
            .as_ref()
            .map(|m| m.reference_density_kg_m3.as_ref())
    }
    pub fn reference_mass_kg(&self) -> Option<f64> {
        self.motion.as_ref().map(|m| m.mass.reference_mass_kg())
    }

    /// Project original-node acceleration through the retained mass, M*a.
    /// Returned dead loads include supported components and are applied once.
    pub fn reference_acceleration_load_n(
        &self,
        acceleration_m_s2: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<(usize, f64)>, MixedThermoplasticTetError> {
        checkpoint(cx)?;
        let motion = self.motion.as_ref().ok_or_else(|| {
            invalid("mixed thermal acceleration requires bound reference inertia")
        })?;
        if acceleration_m_s2.len() != self.nodes_m.len()
            || acceleration_m_s2.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(
                invalid("one finite acceleration per original quadratic node required").into(),
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

    /// Propose one physical interval at supplied original-element temperatures.
    /// Frozen laws are prepared once per element and shared across quadrature;
    /// all nonlinear and thermal retries reuse the same accepted history.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        new_temperatures_k: &[f64],
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_gauge: Option<(usize, f64)>,
        dt_s: f64,
        settings: MixedUpSettings,
        controls: MixedThermoplasticControls,
        cx: &Cx<'_>,
    ) -> Result<MixedThermoplasticTetStep<'a, L>, MixedThermoplasticTetError> {
        checkpoint(cx)?;
        if settings.load_steps != 1
            || !dt_s.is_finite()
            || dt_s <= 0.0
            || new_temperatures_k.len() != self.families.len()
            || ![
                controls.material.equilibrium_work_tolerance_j,
                controls.material.maximum_endpoint_work_defect_j,
                controls.material.yield_tolerance_pa,
                controls.maximum_numerical_kinetic_loss_j,
            ]
            .iter()
            .all(|x| x.is_finite() && *x >= 0.0)
        {
            return Err(invalid("one physical endpoint, positive duration, one temperature per element and finite nonnegative allowances required").into());
        }
        let points = self.mesh.points_per_element;
        let mut old_laws = Vec::with_capacity(self.families.len());
        let mut new_laws = Vec::with_capacity(self.families.len());
        for (e, family) in self.families.iter().enumerate() {
            checkpoint(cx)?;
            let refused = |source| MixedThermoplasticTetError::Material {
                element: e,
                point: 0,
                source,
            };
            let old = family
                .at_temperature_for_inspection(self.temperatures_k[e])
                .map_err(refused)?;
            for (q, point) in self.mesh.elements[e].points.iter().enumerate() {
                checkpoint(cx)?;
                let index = e * points + q;
                let pressure = (0..4)
                    .map(|a| {
                        self.pressure_pa[self.mesh.elements[e].pressure_indices[a]]
                            * point.pressure_shapes[a]
                    })
                    .sum();
                old.validate_mixed_committed(
                    &self.deformation[index],
                    &self.states[index],
                    pressure,
                    self.volumetric,
                    controls.material.yield_tolerance_pa,
                )
                .map_err(|source| MixedThermoplasticTetError::Material {
                    element: e,
                    point: q,
                    source,
                })?;
            }
            old_laws.push(old);
            new_laws.push(
                family
                    .at_temperature_for_step(new_temperatures_k[e], dt_s)
                    .map_err(refused)?,
            );
        }
        let frozen = self
            .states
            .iter()
            .enumerate()
            .map(|(index, accepted)| {
                checkpoint(cx)?;
                let (e, q) = (index / points, index % points);
                new_laws[e].freeze_mixed(accepted).map_err(|source| {
                    MixedThermoplasticTetError::Material {
                        element: e,
                        point: q,
                        source,
                    }
                })
            })
            .collect::<Result<Vec<_>, MixedThermoplasticTetError>>()?;
        let inertia = self
            .motion
            .as_ref()
            .map(|m| {
                BackwardEulerInertia::new(
                    &m.mass,
                    self.displacement_m.as_flattened(),
                    m.velocity_mps.as_flattened(),
                    dt_s,
                    cx,
                )
            })
            .transpose()
            .map_err(inertia_error)?;
        let problem = MixedUpProblem {
            nodes_m: self.nodes_m,
            elements: self.elements,
            material: MixedLawField::PerQuadrature(&frozen),
            volumetric: self.volumetric,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            pressure_gauge,
            budget: self.budget,
            quadrature_points_per_dir: self.quadrature_points_per_dir,
            maximum_quadrature_points: self.maximum_quadrature_points,
            settings,
        };
        let prepared = problem.prepare_with_mesh(Arc::clone(&self.mesh), cx)?;
        let equilibrium = prepared.solve_loaded_from(
            &self.displacement_m,
            &self.pressure_pa,
            inertia.as_ref(),
            cx,
        )?;
        let new_u = equilibrium.displacement_m.as_flattened();
        let next_motion = inertia
            .as_ref()
            .map(|owner| owner.motion(new_u, cx))
            .transpose()
            .map_err(inertia_error)?;
        let mut energy = MixedThermoplasticTetEnergy::default();
        let tolerance = controls.material.equilibrium_work_tolerance_j;
        let mut states = Vec::with_capacity(self.states.len());
        let mut deformation = Vec::with_capacity(self.states.len());
        let mut first_piola_pa = Vec::with_capacity(self.states.len());
        let mut element_caloric_heat_j = vec![0.0; self.elements.len()];
        let (mut old_free, mut new_free, mut old_internal, mut new_internal, mut new_entropy) =
            (0.0, 0.0, 0.0, 0.0, 0.0);
        for (e, element) in self.mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = e * points + q;
                let f = point_deformation(element, point, &equilibrium.displacement_m);
                let pressure = |values: &[f64]| {
                    (0..4)
                        .map(|a| values[element.pressure_indices[a]] * point.pressure_shapes[a])
                        .sum()
                };
                let new_pressure = pressure(&equilibrium.pressure_pa);
                let value = checked_mixed_increment(
                    &old_laws[e],
                    &new_laws[e],
                    &self.deformation[index],
                    &self.states[index],
                    pressure(&self.pressure_pa),
                    &f,
                    new_pressure,
                    self.volumetric,
                    controls.material.yield_tolerance_pa,
                )
                .map_err(|source| MixedThermoplasticTetError::Material {
                    element: e,
                    point: q,
                    source,
                })?;
                let w = point.weight_m3;
                gate(
                    "constitutive thermal identity",
                    w * value.balance_error_j_m3,
                    tolerance,
                )?;
                gate(
                    "constitutive remainder storage",
                    w * (value.new.free_energy_j_m3
                        - value.new.bulk_free_energy_j_m3
                        - value.mechanical.stored_energy_j_m3),
                    tolerance,
                )?;
                gate(
                    "constitutive reduced potential",
                    w * (value.mixed.remainder_potential_pa
                        - value.mechanical.stored_energy_j_m3
                        - value.mechanical.flow_potential_j_m3),
                    tolerance,
                )?;
                old_free += w * value.old.free_energy_j_m3;
                new_free += w * value.new.free_energy_j_m3;
                old_internal += w * value.old.internal_energy_j_m3;
                new_internal += w * value.new.internal_energy_j_m3;
                new_entropy += w * value.new.entropy_j_m3_k;
                energy.free_energy_change_j +=
                    w * (value.new.free_energy_j_m3 - value.old.free_energy_j_m3);
                energy.internal_energy_change_j += w * value.internal_energy_change_j_m3();
                energy.plastic_dissipation_j += w * value.mechanical.plastic_dissipation_j_m3;
                energy.creep_dissipation_j += w * value.mechanical.creep_dissipation_j_m3;
                energy.reversible_heat_j += w * value.reversible_heat_j_m3;
                let heat = w * value.caloric_heat_j_m3();
                energy.caloric_heat_j += heat;
                element_caloric_heat_j[e] += heat;
                energy.coefficient_storage_change_j += w * value.coefficient_storage_change_j_m3;
                energy.thermoelastic_storage_change_j +=
                    w * value.thermoelastic_storage_change_j_m3;
                energy.endpoint_work_defect_j += w * value.endpoint_work_defect_j_m3;
                energy.absolute_endpoint_work_defect_j += w * value.endpoint_work_defect_j_m3.abs();
                energy.incremental_potential_j += w * value.mixed.potential_pa;
                energy.flow_potential_j += w * value.mechanical.flow_potential_j_m3;
                energy.internal_endpoint_work_j += w * value.mechanical_work_j_m3;
                energy.thermal_deformation_work_j += w * value.thermal_deformation_work_j_m3;
                energy.pressure_constraint_work_j +=
                    w * new_pressure * value.mixed.pressure_residual;
                first_piola_pa.push(value.mixed.first_piola_pa);
                deformation.push(f);
                states.push(value.mechanical.state);
            }
        }
        let old_u = self.displacement_m.as_flattened();
        energy.dead_load_endpoint_work_j = prepared
            .external_n
            .iter()
            .zip(new_u.iter().zip(old_u))
            .fold(0.0, |sum, (force, (after, before))| {
                force.mul_add(after - before, sum)
            });
        energy.support_endpoint_work_j =
            equilibrium.reactions_n.iter().fold(0.0, |sum, (d, force)| {
                force.mul_add(new_u[*d] - old_u[*d], sum)
            });
        energy.inertial_endpoint_work_j =
            next_motion.as_ref().map_or(0.0, |(_, a)| a.inertial_work_j);
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j
            - energy.inertial_endpoint_work_j;
        let kinetic = next_motion
            .as_ref()
            .map_or(0.0, |(_, a)| a.kinetic_change_j + a.numerical_loss_j);
        energy.balance_error_j = energy.dead_load_endpoint_work_j + energy.support_endpoint_work_j
            - energy.internal_energy_change_j
            - energy.caloric_heat_j
            - energy.endpoint_work_defect_j
            - kinetic;
        for (quantity, value) in [
            ("accepted free energy", old_free - self.free_energy_j),
            (
                "accepted internal energy",
                old_internal - self.internal_energy_j,
            ),
            (
                "assembled mixed potential",
                energy.incremental_potential_j - equilibrium.potential_j,
            ),
            ("weak pressure work", energy.pressure_constraint_work_j),
            (
                "support/dead-load equilibrium work",
                energy.equilibrium_work_residual_j,
            ),
            (
                "complete thermal mechanical balance",
                energy.balance_error_j,
            ),
        ] {
            gate(quantity, value, tolerance)?;
        }
        gate(
            "absolute thermal endpoint integration defect",
            energy.absolute_endpoint_work_defect_j,
            controls.material.maximum_endpoint_work_defect_j,
        )?;
        if let Some((_, account)) = &next_motion {
            gate(
                "backward-Euler numerical kinetic loss",
                account.numerical_loss_j,
                controls.maximum_numerical_kinetic_loss_j,
            )?;
            gate(
                "native inertial work identity",
                account.balance_error_j,
                tolerance,
            )?;
        }
        if ![
            new_free,
            energy.plastic_dissipation_j,
            energy.creep_dissipation_j,
            energy.flow_potential_j,
        ]
        .iter()
        .all(|x| x.is_finite() && *x >= 0.0)
            || ![
                new_internal,
                new_entropy,
                energy.caloric_heat_j,
                energy.reversible_heat_j,
                energy.coefficient_storage_change_j,
                energy.thermoelastic_storage_change_j,
                energy.thermal_deformation_work_j,
            ]
            .iter()
            .all(|x| x.is_finite())
            || element_caloric_heat_j.iter().any(|x| !x.is_finite())
        {
            return Err(
                invalid("representable mixed thermal storage, heat and work required").into(),
            );
        }
        checkpoint(cx)?;
        let (motion, inertia) = match (next_motion, &self.motion) {
            (Some((velocity_mps, account)), Some(old)) => (
                Some(ThermalMotion {
                    velocity_mps,
                    reference_density_kg_m3: Arc::clone(&old.reference_density_kg_m3),
                    mass: Arc::clone(&old.mass),
                }),
                Some(account),
            ),
            (None, None) => (None, None),
            _ => unreachable!("native inertia and retained motion are bound together"),
        };
        let solid = Self {
            nodes_m: self.nodes_m,
            elements: self.elements,
            families: Arc::clone(&self.families),
            temperatures_k: new_temperatures_k.to_vec(),
            volumetric: self.volumetric,
            mesh: Arc::clone(&self.mesh),
            budget: self.budget,
            quadrature_points_per_dir: self.quadrature_points_per_dir,
            maximum_quadrature_points: self.maximum_quadrature_points,
            displacement_m: equilibrium.displacement_m.clone(),
            pressure_pa: equilibrium.pressure_pa.clone(),
            states,
            deformation,
            first_piola_pa,
            free_energy_j: new_free,
            internal_energy_j: new_internal,
            entropy_j_k: new_entropy,
            motion,
        };
        Ok(MixedThermoplasticTetStep {
            solid,
            equilibrium,
            energy,
            inertia,
            element_caloric_heat_j,
        })
    }
}
