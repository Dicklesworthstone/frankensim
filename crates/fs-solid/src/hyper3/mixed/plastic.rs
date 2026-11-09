//! Retained finite inelasticity on the native P2/P1 mixed tetrahedra.
//!
//! The original material owns the isochoric return, kinetics, hardening and
//! consistent tangent. Each quadrature point keeps its complete native history,
//! including independent recovering resistance when selected. Newton trials
//! borrow that same history and physical interval; only a complete
//! converged, work-checked endpoint publishes new displacement, pressure and
//! history together. Original source-resolved laws remain borrowed intact.
//!
//! Physical storage is the returned deviatoric/hardening energy plus the
//! original bulk energy at the pressure-reconstructed volume, `U(J(p))`.
//! Continuous P1 pressure enforces volume weakly; substituting `U(det F)` would
//! restore a displacement penalty and change this mixed formulation. Explicit
//! exact incompressibility has no bulk storage. Flow potential, physical
//! inelastic heat, pressure-constraint work and endpoint quadrature error remain
//! separate. This is an isothermal quasistatic path for rate-independent J2,
//! Perzyna overstress, Norton creep and their native combined/recovering laws.
//! All kinetic points share one original physical interval. The original
//! finite-history limits and local convex bulk branch remain enforced. Mesh
//! stability, locking accuracy and injectivity require their own evidence.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::plastic::finite::mixed::MixedFiniteInelasticMaterial;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Error};

use super::element::PreparedMixedMesh;
use super::{
    MixedLawField, MixedUpError, MixedUpProblem, MixedUpSettings, MixedUpSolution, checkpoint,
    invalid, point_deformation,
};
use crate::hyper3::plastic::FinitePlasticStepControls;
use crate::linear3::TetAssemblyBudget;

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// Original laws without property averaging or replacement of their sources.
#[derive(Debug)]
pub enum MixedPlasticLawField<'a, L: MixedFiniteInelasticMaterial = FiniteJ2> {
    /// One original material for all points; histories are still distinct.
    Uniform(&'a L),
    /// One original material per tetrahedron.
    PerElement(&'a [L]),
    /// Element-major, then the native deterministic Duffy quadrature order.
    PerQuadrature(&'a [L]),
}

impl<L: MixedFiniteInelasticMaterial> Copy for MixedPlasticLawField<'_, L> {}
impl<L: MixedFiniteInelasticMaterial> Clone for MixedPlasticLawField<'_, L> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<L: MixedFiniteInelasticMaterial> MixedPlasticLawField<'_, L> {
    fn at(&self, element: usize, point: usize, points: usize) -> &L {
        match self {
            Self::Uniform(law) => law,
            Self::PerElement(laws) => &laws[element],
            Self::PerQuadrature(laws) => &laws[element * points + point],
        }
    }
}

/// Independent physical and incremental accounts, integrated over reference volume.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MixedPlasticTetEnergy {
    /// Accepted deviatoric elastic and energetic-hardening storage, J.
    pub remainder_before_j: f64,
    /// Returned deviatoric elastic and energetic-hardening storage, J.
    pub remainder_after_j: f64,
    /// Accepted physical bulk storage reconstructed from pressure, J.
    pub bulk_before_j: f64,
    /// Returned physical bulk storage reconstructed from pressure, J.
    pub bulk_after_j: f64,
    /// Total accepted physical storage, J.
    pub stored_before_j: f64,
    /// Total returned physical storage, J.
    pub stored_after_j: f64,
    /// Native plastic/viscous channel, excluding numerical defects. In the
    /// native sole-Norton convention this is creep heat; combined laws retain
    /// creep in its independent channel below.
    pub plastic_dissipation_j: f64,
    /// Physical creep heat from the native combined creep/plastic channel.
    pub creep_dissipation_j: f64,
    /// Integrated native flow potential, distinct from stored energy.
    pub flow_potential_j: f64,
    /// Bulk-free reduced potential used during Newton, J.
    pub incremental_remainder_potential_j: f64,
    /// Complete mixed saddle potential, including the pressure dual, J.
    pub mixed_potential_j: f64,
    /// Integral of p times the weak volume residual, J; not material heat.
    pub pressure_constraint_work_j: f64,
    /// Integral of endpoint P contracted with the actual change of F, J.
    pub internal_endpoint_work_j: f64,
    /// Original dead nodal/body load dotted with the accepted displacement change.
    pub dead_load_endpoint_work_j: f64,
    /// Endpoint support reaction dotted with its prescribed displacement change.
    pub support_endpoint_work_j: f64,
    /// External endpoint work minus internal endpoint work, J.
    pub equilibrium_work_residual_j: f64,
    /// Internal endpoint work minus storage change and total physical heat, J.
    /// Signed first-order quadrature defect; never added to material heat.
    pub endpoint_work_defect_j: f64,
}

impl MixedPlasticTetEnergy {
    /// Complete physical inelastic heat. Neither the incremental flow potential
    /// nor signed endpoint integration defects are caloric dissipation.
    pub fn inelastic_dissipation_j(&self) -> f64 {
        self.plastic_dissipation_j + self.creep_dissipation_j
    }
}

/// Native numerical/material refusal or an explicit physical-work gate.
#[derive(Clone, Debug, PartialEq)]
pub enum MixedPlasticTetError {
    /// Original mixed geometry, solve, budget or cancellation diagnosis.
    Native(MixedUpError),
    /// Original material refusal at one retained quadrature point.
    Material {
        element: usize,
        point: usize,
        source: FiniteJ2Error,
    },
    /// A nonfinite result or an exceeded caller-supplied absolute work allowance.
    EnergyGate {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}

impl From<MixedUpError> for MixedPlasticTetError {
    fn from(error: MixedUpError) -> Self {
        Self::Native(error)
    }
}

impl core::fmt::Display for MixedPlasticTetError {
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
                    "mixed plastic element {element}, point {point}: {source}"
                )
            }
            Self::EnergyGate {
                quantity,
                value_j,
                maximum_abs_j,
            } => write!(
                f,
                "mixed plastic {quantity}: {value_j:e} J exceeds absolute allowance {maximum_abs_j:e} J"
            ),
        }
    }
}

impl core::error::Error for MixedPlasticTetError {
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
) -> Result<(), MixedPlasticTetError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(MixedPlasticTetError::EnergyGate {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

/// Immutable accepted mixed configuration and complete original plastic history.
/// Reference geometry is prepared once and shared by accepted successors.
#[derive(Debug)]
pub struct MixedPlasticTetSolid<'a, L: MixedFiniteInelasticMaterial = FiniteJ2> {
    nodes_m: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: MixedPlasticLawField<'a, L>,
    duration_s: Option<f64>,
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
    remainder_energy_j: f64,
    bulk_energy_j: f64,
}

/// A complete proposed endpoint. Dropping it preserves the previous solid.
#[derive(Debug)]
pub struct MixedPlasticTetStep<'a, L: MixedFiniteInelasticMaterial = FiniteJ2> {
    /// Displacement, pressure, material binding and accepted histories together.
    pub solid: MixedPlasticTetSolid<'a, L>,
    /// Native pressure/displacement equilibrium and reactions.
    pub equilibrium: MixedUpSolution,
    /// Independently evaluated physical storage, inelastic heat and work.
    pub energy: MixedPlasticTetEnergy,
}

impl<'a, L: MixedFiniteInelasticMaterial> MixedPlasticTetSolid<'a, L> {
    /// Prepare the virgin, stress-free reference configuration.
    /// Finite compliance retains the original bulk modulus; exact
    /// incompressibility is an explicit alternative with no bulk penalty.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        nodes_m: &'a [[f64; 3]],
        elements: &'a [[usize; 10]],
        material: MixedPlasticLawField<'a, L>,
        volumetric: VolumetricTreatment,
        budget: TetAssemblyBudget,
        quadrature_points_per_dir: usize,
        maximum_quadrature_points: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, MixedPlasticTetError> {
        checkpoint(cx)?;
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
            .ok_or_else(|| invalid("mixed plastic material count overflow"))?;
        let counts = match material {
            MixedPlasticLawField::Uniform(_) => None,
            MixedPlasticLawField::PerElement(laws) => Some((laws.len(), elements.len())),
            MixedPlasticLawField::PerQuadrature(laws) => Some((laws.len(), count)),
        };
        if counts.is_some_and(|(got, expected)| got != expected) {
            return Err(
                invalid("mixed plastic laws must cover every admitted element or point").into(),
            );
        }
        let mut states = Vec::with_capacity(count);
        let mut first_piola_pa = Vec::with_capacity(count);
        let mut duration_s = None;
        let mut remainder_energy_j = 0.0;
        let mut bulk_energy_j = 0.0;
        for (e, element) in mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let law = material.at(e, q, mesh.points_per_element);
                if let Some(duration) = law.duration_s() {
                    if !duration.is_finite()
                        || duration <= 0.0
                        || duration_s.is_some_and(|old| old != duration)
                    {
                        return Err(invalid(
                            "every kinetic quadrature point must share one finite positive physical duration",
                        )
                        .into());
                    }
                    duration_s = Some(duration);
                }
                let state = law.initial_state();
                let source = |source| MixedPlasticTetError::Material {
                    element: e,
                    point: q,
                    source,
                };
                let frozen = law.freeze_remainder(&state).map_err(source)?;
                let native = |error| MixedUpError::Material {
                    element: e,
                    point: q,
                    error,
                };
                frozen.validate_mixed().map_err(native)?;
                let mixed = frozen
                    .mixed_state(&IDENTITY, 0.0, volumetric)
                    .map_err(native)?;
                remainder_energy_j += point.weight_m3
                    * law
                        .stored_remainder_energy(&IDENTITY, &state)
                        .map_err(source)?;
                bulk_energy_j += point.weight_m3 * mixed.bulk_stored_pa;
                first_piola_pa.push(mixed.first_piola_pa);
                states.push(state);
            }
        }
        if !(remainder_energy_j + bulk_energy_j).is_finite() {
            return Err(invalid("mixed plastic initial storage is not representable").into());
        }
        checkpoint(cx)?;
        Ok(Self {
            nodes_m,
            elements,
            material,
            duration_s,
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
            remainder_energy_j,
            bulk_energy_j,
        })
    }

    /// Accepted node-major displacement, in metres.
    pub fn displacement_m(&self) -> &[[f64; 3]] {
        &self.displacement_m
    }
    /// Sorted original corner-node IDs owning the continuous pressure field.
    pub fn pressure_nodes(&self) -> &[usize] {
        &self.mesh.pressure_nodes
    }
    /// Tensile-positive bulk pressure, in the order of [`Self::pressure_nodes`].
    pub fn pressure_pa(&self) -> &[f64] {
        &self.pressure_pa
    }
    /// Complete native histories in element-major Duffy order, including
    /// separately recovering resistance for the native recovery owners.
    pub fn quadrature_states(&self) -> &[L::State] {
        &self.states
    }
    /// Shared original physical interval for explicit kinetics, or `None` for
    /// a completely rate-independent field. Only an accepted advance consumes
    /// an interval; inspection and Newton trials never age the accepted state.
    pub const fn duration_s(&self) -> Option<f64> {
        self.duration_s
    }
    /// Accepted original deformation gradient at each quadrature point.
    pub fn deformation_gradients(&self) -> &[[f64; 9]] {
        &self.deformation
    }
    /// Accepted full mixed first Piola stress, including pressure, in Pa.
    pub fn first_piola_pa(&self) -> &[[f64; 9]] {
        &self.first_piola_pa
    }
    /// Native quadrature count for each element.
    pub fn quadrature_points_per_element(&self) -> usize {
        self.mesh.points_per_element
    }
    /// Original mesh volume, in cubic metres.
    pub fn reference_volume_m3(&self) -> f64 {
        self.mesh.total_volume_m3
    }
    /// Deviatoric elastic and energetic-hardening storage, in joules.
    pub fn remainder_energy_j(&self) -> f64 {
        self.remainder_energy_j
    }
    /// Physical original bulk storage, in joules; zero for exact incompressibility.
    pub fn bulk_energy_j(&self) -> f64 {
        self.bulk_energy_j
    }
    /// Complete physical stored energy, in joules.
    pub fn stored_energy_j(&self) -> f64 {
        self.remainder_energy_j + self.bulk_energy_j
    }

    /// Solve one physical loading endpoint from the immutable accepted history.
    /// `settings.load_steps` must be one: physical substeps belong to the caller,
    /// so numerical Newton iterations never silently advance the plastic history.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_gauge: Option<(usize, f64)>,
        settings: MixedUpSettings,
        controls: FinitePlasticStepControls,
        cx: &Cx<'_>,
    ) -> Result<MixedPlasticTetStep<'a, L>, MixedPlasticTetError> {
        checkpoint(cx)?;
        if settings.load_steps != 1
            || ![
                controls.equilibrium_work_tolerance_j,
                controls.maximum_endpoint_work_defect_j,
                controls.yield_tolerance_pa,
            ]
            .iter()
            .all(|value| value.is_finite() && *value >= 0.0)
        {
            return Err(invalid(
                "one physical endpoint and finite nonnegative work/yield allowances required",
            )
            .into());
        }
        let points = self.mesh.points_per_element;
        let frozen = self
            .states
            .iter()
            .enumerate()
            .map(|(index, state)| {
                checkpoint(cx)?;
                let (element, point) = (index / points, index % points);
                let law = self.material.at(element, point, points);
                let refused = |source| MixedPlasticTetError::Material {
                    element,
                    point,
                    source,
                };
                law.validate_remainder_committed(
                    &self.deformation[index],
                    state,
                    controls.yield_tolerance_pa,
                )
                .map_err(refused)?;
                law.freeze_remainder(state).map_err(refused)
            })
            .collect::<Result<Vec<_>, MixedPlasticTetError>>()?;
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
        let equilibrium = prepared.solve_from(&self.displacement_m, &self.pressure_pa, cx)?;
        let u = &equilibrium.displacement_m;
        let mut energy = MixedPlasticTetEnergy {
            remainder_before_j: self.remainder_energy_j,
            bulk_before_j: self.bulk_energy_j,
            stored_before_j: self.stored_energy_j(),
            ..Default::default()
        };
        let mut states = Vec::with_capacity(self.states.len());
        let mut deformation = Vec::with_capacity(self.states.len());
        let mut first_piola_pa = Vec::with_capacity(self.states.len());
        for (e, element) in self.mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = e * points + q;
                let f = point_deformation(element, point, u);
                let pressure = (0..4)
                    .map(|a| {
                        equilibrium.pressure_pa[element.pressure_indices[a]]
                            * point.pressure_shapes[a]
                    })
                    .sum();
                let mixed = frozen[index]
                    .mixed_state(&f, pressure, self.volumetric)
                    .map_err(|error| MixedUpError::Material {
                        element: e,
                        point: q,
                        error,
                    })?;
                let law = self.material.at(e, q, points);
                let refused = |source| MixedPlasticTetError::Material {
                    element: e,
                    point: q,
                    source,
                };
                let response = law
                    .remainder_response(&f, &self.states[index])
                    .map_err(refused)?;
                law.validate_remainder_committed(&f, &response.state, controls.yield_tolerance_pa)
                    .map_err(refused)?;
                let stored = law
                    .stored_remainder_energy(&f, &response.state)
                    .map_err(refused)?;
                let w = point.weight_m3;
                gate(
                    "constitutive stored-energy consistency",
                    w * (stored - response.stored_energy_j_m3),
                    controls.equilibrium_work_tolerance_j,
                )?;
                gate(
                    "constitutive reduced-potential consistency",
                    w * (mixed.remainder_potential_pa - stored - response.flow_potential_j_m3),
                    controls.equilibrium_work_tolerance_j,
                )?;
                energy.remainder_after_j += w * stored;
                energy.bulk_after_j += w * mixed.bulk_stored_pa;
                energy.plastic_dissipation_j += w * response.plastic_dissipation_j_m3;
                energy.creep_dissipation_j += w * response.creep_dissipation_j_m3;
                energy.flow_potential_j += w * response.flow_potential_j_m3;
                energy.incremental_remainder_potential_j += w * mixed.remainder_potential_pa;
                energy.mixed_potential_j += w * mixed.potential_pa;
                energy.pressure_constraint_work_j += w * pressure * mixed.pressure_residual;
                energy.internal_endpoint_work_j += w
                    * (0..9).fold(0.0, |sum, c| {
                        mixed.first_piola_pa[c].mul_add(f[c] - self.deformation[index][c], sum)
                    });
                first_piola_pa.push(mixed.first_piola_pa);
                deformation.push(f);
                states.push(response.state);
            }
        }
        energy.stored_after_j = energy.remainder_after_j + energy.bulk_after_j;
        let old_u = self.displacement_m.as_flattened();
        let new_u = u.as_flattened();
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
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j
            - (energy.stored_after_j - energy.stored_before_j)
            - energy.inelastic_dissipation_j();
        for (quantity, value) in [
            (
                "assembled mixed potential",
                energy.mixed_potential_j - equilibrium.potential_j,
            ),
            (
                "pressure constraint work",
                energy.pressure_constraint_work_j,
            ),
            (
                "support/dead-load equilibrium work",
                energy.equilibrium_work_residual_j,
            ),
        ] {
            gate(quantity, value, controls.equilibrium_work_tolerance_j)?;
        }
        gate(
            "signed endpoint integration defect",
            energy.endpoint_work_defect_j,
            controls.maximum_endpoint_work_defect_j,
        )?;
        if ![
            energy.stored_after_j,
            energy.plastic_dissipation_j,
            energy.creep_dissipation_j,
            energy.flow_potential_j,
        ]
        .iter()
        .all(|value| value.is_finite() && *value >= 0.0)
        {
            return Err(invalid(
                "finite nonnegative plastic storage, flow potential and physical loss required",
            )
            .into());
        }
        checkpoint(cx)?;
        let solid = Self {
            nodes_m: self.nodes_m,
            elements: self.elements,
            material: self.material,
            duration_s: self.duration_s,
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
            remainder_energy_j: energy.remainder_after_j,
            bulk_energy_j: energy.bulk_after_j,
        };
        Ok(MixedPlasticTetStep {
            solid,
            equilibrium,
            energy,
        })
    }
}
