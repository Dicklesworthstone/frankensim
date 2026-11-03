//! Retained finite-strain viscoelasticity on the native P2/P1 mixed solid.
//!
//! Every quadrature point retains its original law and remainder history. The
//! accepted history and interval are frozen throughout the native mixed Newton
//! solve. A successful candidate publishes displacement, pressure, histories,
//! stresses, physical storage and time together. The reference mesh is shared.
//!
//! Physical storage is the native remainder plus `U(J(p))`, where pressure
//! reconstructs the bulk volume. It is not `U(det F)`: the continuous pressure
//! equations impose volume weakly. Exact incompressibility never evaluates the
//! unused bulk penalty, including during accepted material updates. The mixed
//! saddle potential, native relaxation loss and signed endpoint-work defect
//! remain separate. This is an isothermal quasistatic path; no loss is heat.
//! [`dynamics`] adds original consistent P2 mass to the same endpoint solve.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::mixed::{MixedFiniteStrainPotential, VolumetricTreatment};
use fs_material::visco_hyper::{
    FrozenViscoRemainder, MAX_FINITE_STRAIN_BRANCHES, ViscoHyperelastic, ViscoRemainderState,
};

use super::element::PreparedMixedMesh;
use super::{
    MixedLawField, MixedUpError, MixedUpProblem, MixedUpSettings, MixedUpSolution, checkpoint,
    invalid,
};
use crate::hyper3::visco::FiniteViscoStepControls;
use crate::linear3::TetAssemblyBudget;

/// Native mixed motion with independent physical and numerical loss accounts.
pub mod dynamics;

type AcceptedMotion = (Vec<[f64; 3]>, crate::inertia::InertialEnergy);

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

/// Immutable original spectra; histories are always distinct at every point.
#[derive(Clone, Copy, Debug)]
pub enum MixedViscoLawField<'a> {
    Uniform(&'a ViscoHyperelastic),
    PerElement(&'a [ViscoHyperelastic]),
    /// Element-major, then the native deterministic Duffy quadrature order.
    PerQuadrature(&'a [ViscoHyperelastic]),
}

impl<'a> MixedViscoLawField<'a> {
    fn at(self, element: usize, point: usize, points: usize) -> &'a ViscoHyperelastic {
        match self {
            Self::Uniform(law) => law,
            Self::PerElement(laws) => &laws[element],
            Self::PerQuadrature(laws) => &laws[element * points + point],
        }
    }
}

/// Physical and incremental accounts for one accepted quasistatic interval.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MixedViscoTetEnergy {
    pub remainder_before_j: f64,
    pub remainder_after_j: f64,
    pub bulk_before_j: f64,
    pub bulk_after_j: f64,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Native nonnegative branch loss for this interval, not cumulative.
    pub dissipated_j: f64,
    /// Frozen-history loading work of the original bulk-free remainder.
    pub remainder_loading_work_j: f64,
    /// Native reduced remainder potential; not physical storage.
    pub incremental_remainder_potential_j: f64,
    /// Complete mixed saddle potential, including the pressure dual.
    pub mixed_potential_j: f64,
    /// Integral of `p*(det F-J(p))`; zero at exact pressure equilibrium.
    pub pressure_constraint_work_j: f64,
    /// Remainder loading work minus remainder storage change and native loss.
    pub material_balance_error_j: f64,
    pub internal_endpoint_work_j: f64,
    pub dead_load_endpoint_work_j: f64,
    pub support_endpoint_work_j: f64,
    /// Native inertial endpoint work; zero for quasistatic intervals.
    pub inertial_endpoint_work_j: f64,
    /// Dead/support work minus material and inertial endpoint work.
    pub equilibrium_work_residual_j: f64,
    /// Internal endpoint work minus total physical storage change and loss.
    /// Signed quadrature defect, explicitly gated and never counted as heat.
    pub endpoint_work_defect_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum MixedViscoTetError {
    Native(MixedUpError),
    /// A mass-only impulse cannot enforce the exact pressure constraint.
    ExactIncompressibleImpulse,
    EnergyGate {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}

impl From<MixedUpError> for MixedViscoTetError {
    fn from(error: MixedUpError) -> Self {
        Self::Native(error)
    }
}
impl core::fmt::Display for MixedViscoTetError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => error.fmt(f),
            Self::ExactIncompressibleImpulse => write!(
                f,
                "exact-incompressible impulses require a constrained pressure solve"
            ),
            Self::EnergyGate {
                quantity,
                value_j,
                maximum_abs_j,
            } => write!(
                f,
                "mixed viscoelastic {quantity}: {value_j:e} J exceeds {maximum_abs_j:e} J"
            ),
        }
    }
}
impl core::error::Error for MixedViscoTetError {}

fn gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), MixedViscoTetError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(MixedViscoTetError::EnergyGate {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}

/// Immutable accepted solid. A rejected candidate cannot advance any history.
#[derive(Debug)]
pub struct MixedViscoTetSolid<'a> {
    nodes_m: &'a [[f64; 3]],
    elements: &'a [[usize; 10]],
    material: MixedViscoLawField<'a>,
    volumetric: VolumetricTreatment,
    mesh: Arc<PreparedMixedMesh>,
    budget: TetAssemblyBudget,
    quadrature_points_per_dir: usize,
    maximum_quadrature_points: usize,
    displacement_m: Vec<[f64; 3]>,
    pressure_pa: Vec<f64>,
    states: Vec<ViscoRemainderState>,
    first_piola_pa: Vec<[f64; 9]>,
    remainder_energy_j: f64,
    bulk_energy_j: f64,
    time_s: f64,
}

/// Complete proposed endpoint. Dropping it leaves the previous solid intact.
#[derive(Debug)]
pub struct MixedViscoTetStep<'a> {
    pub solid: MixedViscoTetSolid<'a>,
    pub equilibrium: MixedUpSolution,
    pub energy: MixedViscoTetEnergy,
}

impl<'a> MixedViscoTetSolid<'a> {
    /// Prepare the undeformed relaxed solid and its immutable original geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        nodes_m: &'a [[f64; 3]],
        elements: &'a [[usize; 10]],
        material: MixedViscoLawField<'a>,
        volumetric: VolumetricTreatment,
        budget: TetAssemblyBudget,
        quadrature_points_per_dir: usize,
        maximum_quadrature_points: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, MixedViscoTetError> {
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
        let count = elements.len() * mesh.points_per_element;
        let expected = match material {
            MixedViscoLawField::Uniform(_) => None,
            MixedViscoLawField::PerElement(laws) => Some((laws.len(), elements.len())),
            MixedViscoLawField::PerQuadrature(laws) => Some((laws.len(), count)),
        };
        if expected.is_some_and(|(got, wanted)| got != wanted) {
            return Err(
                invalid("mixed viscoelastic law field must cover the admitted points").into(),
            );
        }
        let mut states = Vec::with_capacity(count);
        let mut first_piola_pa = Vec::with_capacity(count);
        let mut remainder_energy_j = 0.0;
        for (e, element) in mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let law = material.at(e, q, mesh.points_per_element);
                if law.branches.len() > MAX_FINITE_STRAIN_BRANCHES {
                    return Err(invalid("mixed viscous branch count exceeds the native cap").into());
                }
                let material_error = |error| MixedUpError::Material {
                    element: e,
                    point: q,
                    error,
                };
                let initial = law.initial_remainder_state();
                let frozen = law
                    .freeze_remainder(&initial, 0.0)
                    .map_err(material_error)?;
                frozen.validate_mixed().map_err(material_error)?;
                let mixed = frozen
                    .mixed_state(&IDENTITY, 0.0, volumetric)
                    .map_err(material_error)?;
                let response = law
                    .step_remainder_response(&initial, &IDENTITY, 0.0)
                    .map_err(material_error)?;
                remainder_energy_j += point.weight_m3 * response.state.stored;
                first_piola_pa.push(mixed.first_piola_pa);
                states.push(response.state);
            }
        }
        if !remainder_energy_j.is_finite() {
            return Err(invalid("mixed initial physical storage is unrepresentable").into());
        }
        let pressure_pa = vec![0.0; mesh.pressure_nodes.len()];
        checkpoint(cx)?;
        Ok(Self {
            nodes_m,
            elements,
            material,
            volumetric,
            mesh,
            budget,
            quadrature_points_per_dir,
            maximum_quadrature_points,
            displacement_m: vec![[0.0; 3]; nodes_m.len()],
            pressure_pa,
            states,
            first_piola_pa,
            remainder_energy_j,
            bulk_energy_j: 0.0,
            time_s: 0.0,
        })
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
    /// Element-major then Duffy order; ledgers contain ONLY physical remainder.
    pub fn quadrature_states(&self) -> &[ViscoRemainderState] {
        &self.states
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
    pub fn remainder_energy_j(&self) -> f64 {
        self.remainder_energy_j
    }
    pub fn bulk_energy_j(&self) -> f64 {
        self.bulk_energy_j
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.remainder_energy_j + self.bulk_energy_j
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }

    /// Solve one physical endpoint with original per-point histories frozen.
    /// `dt_s=0` is an instantaneous elastic jump. Use one solver load endpoint;
    /// accepted physical substeps belong to the caller, never to Newton trials.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_gauge: Option<(usize, f64)>,
        settings: MixedUpSettings,
        controls: FiniteViscoStepControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoTetStep<'a>, MixedViscoTetError> {
        self.advance_loaded(
            dt_s,
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            pressure_gauge,
            settings,
            controls,
            None,
            cx,
        )
        .map(|(step, _)| step)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn advance_loaded(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        pressure_gauge: Option<(usize, f64)>,
        settings: MixedUpSettings,
        controls: FiniteViscoStepControls,
        inertia: Option<&crate::inertia::BackwardEulerInertia<'_>>,
        cx: &Cx<'_>,
    ) -> Result<(MixedViscoTetStep<'a>, Option<AcceptedMotion>), MixedViscoTetError> {
        checkpoint(cx)?;
        let time_s = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s < 0.0
            || !time_s.is_finite()
            || settings.load_steps != 1
            || ![
                controls.equilibrium_work_tolerance_j,
                controls.maximum_endpoint_work_defect_j,
            ]
            .iter()
            .all(|v| v.is_finite() && *v >= 0.0)
        {
            return Err(invalid(
                "one endpoint, nonnegative finite interval and finite work allowances required",
            )
            .into());
        }
        let points = self.mesh.points_per_element;
        let frozen: Vec<FrozenViscoRemainder<'_>> = self
            .states
            .iter()
            .enumerate()
            .map(|(index, state)| {
                checkpoint(cx)?;
                let (element, point) = (index / points, index % points);
                self.material
                    .at(element, point, points)
                    .freeze_remainder(state, dt_s)
                    .map_err(|error| MixedUpError::Material {
                        element,
                        point,
                        error,
                    })
            })
            .collect::<Result<_, _>>()?;
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
        let equilibrium =
            prepared.solve_loaded_from(&self.displacement_m, &self.pressure_pa, inertia, cx)?;
        let u = &equilibrium.displacement_m;
        let motion = inertia
            .map(|owner| {
                owner
                    .motion(u.as_flattened(), cx)
                    .map_err(super::inertia_error)
            })
            .transpose()?;
        let mut energy = MixedViscoTetEnergy {
            remainder_before_j: self.remainder_energy_j,
            bulk_before_j: self.bulk_energy_j,
            stored_before_j: self.stored_energy_j(),
            inertial_endpoint_work_j: motion
                .as_ref()
                .map_or(0.0, |(_, account)| account.inertial_work_j),
            ..Default::default()
        };
        let mut states = Vec::with_capacity(self.states.len());
        let mut first_piola_pa = Vec::with_capacity(self.states.len());
        for (e, element) in self.mesh.elements.iter().enumerate() {
            for (q, point) in element.points.iter().enumerate() {
                checkpoint(cx)?;
                let index = e * points + q;
                let mut deformation = IDENTITY;
                for a in 1..10 {
                    for i in 0..3 {
                        for j in 0..3 {
                            deformation[3 * i + j] += (u[element.connectivity[a]][i]
                                - u[element.connectivity[0]][i])
                                * point.gradients[a][j];
                        }
                    }
                }
                let pressure = (0..4)
                    .map(|a| {
                        equilibrium.pressure_pa[element.pressure_indices[a]]
                            * point.pressure_shapes[a]
                    })
                    .sum();
                let material_error = |error| MixedUpError::Material {
                    element: e,
                    point: q,
                    error,
                };
                let mixed = frozen[index]
                    .mixed_state(&deformation, pressure, self.volumetric)
                    .map_err(material_error)?;
                let response = self
                    .material
                    .at(e, q, points)
                    .step_remainder_response(&self.states[index], &deformation, dt_s)
                    .map_err(material_error)?;
                let w = point.weight_m3;
                energy.remainder_after_j += w * response.state.stored;
                energy.bulk_after_j += w * mixed.bulk_stored_pa;
                energy.dissipated_j += w * response.dissipated_pa;
                energy.remainder_loading_work_j += w * response.loading_work_pa;
                energy.incremental_remainder_potential_j += w * mixed.remainder_potential_pa;
                energy.mixed_potential_j += w * mixed.potential_pa;
                energy.pressure_constraint_work_j += w * pressure * mixed.pressure_residual;
                energy.internal_endpoint_work_j += w
                    * (0..9).fold(0.0, |sum, c| {
                        mixed.first_piola_pa[c].mul_add(
                            deformation[c] - self.states[index].deformation_gradient[c],
                            sum,
                        )
                    });
                first_piola_pa.push(mixed.first_piola_pa);
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
        energy.material_balance_error_j = energy.remainder_loading_work_j
            - (energy.remainder_after_j - energy.remainder_before_j)
            - energy.dissipated_j;
        energy.equilibrium_work_residual_j = energy.dead_load_endpoint_work_j
            + energy.support_endpoint_work_j
            - energy.internal_endpoint_work_j
            - energy.inertial_endpoint_work_j;
        energy.endpoint_work_defect_j = energy.internal_endpoint_work_j
            - (energy.stored_after_j - energy.stored_before_j)
            - energy.dissipated_j;
        for (quantity, value) in [
            ("native remainder balance", energy.material_balance_error_j),
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
        if !energy.stored_after_j.is_finite()
            || !energy.dissipated_j.is_finite()
            || energy.dissipated_j < 0.0
        {
            return Err(invalid(
                "finite physical storage and nonnegative native relaxation loss required",
            )
            .into());
        }
        checkpoint(cx)?;
        let solid = Self {
            nodes_m: self.nodes_m,
            elements: self.elements,
            material: self.material,
            volumetric: self.volumetric,
            mesh: Arc::clone(&self.mesh),
            budget: self.budget,
            quadrature_points_per_dir: self.quadrature_points_per_dir,
            maximum_quadrature_points: self.maximum_quadrature_points,
            displacement_m: equilibrium.displacement_m.clone(),
            pressure_pa: equilibrium.pressure_pa.clone(),
            states,
            first_piola_pa,
            remainder_energy_j: energy.remainder_after_j,
            bulk_energy_j: energy.bulk_after_j,
            time_s,
        };
        Ok((
            MixedViscoTetStep {
                solid,
                equilibrium,
                energy,
            },
            motion,
        ))
    }
}
