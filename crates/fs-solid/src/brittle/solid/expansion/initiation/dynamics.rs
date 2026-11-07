//! Thermal elastic motion, source-detected tensile onset and midpoint replay.
//!
//! The native centered-inertia solver uses the original consistent reference
//! mass and the exact mean stress of the next elastic coefficients along the
//! supplied MECHANICAL strain chord. Total displacement and integrated thermal
//! expansion both contribute to that chord. Physical endpoint stress is restored
//! before the source detector inspects it. A coefficient-only activation is a
//! distinct result obtained before an elastic predictor is solved.
//!
//! Both triggers transfer the OLD accepted thermal state, velocity, mass, kinetic
//! storage, clock and consumed step budget to the existing thermal brittle motion
//! owner. Its midpoint solver replays the complete requested physical interval
//! with the original interval-MEAN forces and endpoint supports. No over-strength
//! elastic predictor is accepted and no strain fraction becomes a physical time.
//!
//! Parameter work, eigenstrain work and fracture loss retain separate accounts.
//! Centered inertia introduces no numerical kinetic damping; the coefficient-first
//! thermal split still makes no second-order or continuous-trajectory claim.
//! Sources, reversals, thermal changes and force means must be resolved by the
//! caller. This is a small-strain fixed-plane tensile model, not a heat equation.

use std::sync::Arc;

use super::{
    ThermalTensileElasticTetSolid, ThermalTensileEndpoint, ThermalTensileTetActivation,
    ThermalTensileTetBands, ThermalTensileTetInitiationError, ThermalTensileTetOnset,
    ThermalTensileTetRequest, mechanical,
};
use crate::brittle::solid::dynamics::BrittleTetDynamicError;
use crate::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use crate::brittle::solid::expansion::dynamics::ThermalBrittleTetMotion;
use crate::brittle::solid::expansion::dynamics::midpoint::{
    ThermalBrittleTetMidpointEnergy, ThermalBrittleTetMidpointStep,
};
use crate::brittle::solid::{BrittleTetError, checkpoint, energy_gate};
use crate::inertia::{InertialEnergy, P1ReferenceMass};
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    SmallStrainTetSolution,
};
use fs_evidence::ModelCard;
use fs_exec::Cx;
use fs_material::brittle::initiation::TensileInitiationSettings;
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;
use fs_material::brittle::source::initiation::{ResolvedTensileActivation, ResolvedTensileOnset};
use fs_material::tensor::contract;
use fs_material::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};

#[derive(Clone, Debug)]
struct RetainedThermalInertia {
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    mass: Arc<P1ReferenceMass>,
    kinetic_energy_j: f64,
    time_s: f64,
    maximum_steps: usize,
    accepted_steps: usize,
}

/// Accepted elastic thermal motion before physical crack planes are selected.
#[derive(Debug)]
pub struct ThermalTensileElasticTetMotion<'a> {
    solid: ThermalTensileElasticTetSolid<'a>,
    inertia: RetainedThermalInertia,
    initiation: TensileInitiationSettings,
}

#[derive(Debug)]
pub struct ThermalTensileElasticTetMidpointStep<'a> {
    pub motion: ThermalTensileElasticTetMotion<'a>,
    /// Work-conjugate interval response, separate from physical endpoint stress.
    pub mean_element_stress_pa: Vec<Voigt>,
    pub energy: ThermalBrittleTetMidpointEnergy,
}

#[derive(Debug)]
pub enum ThermalTensileTetMidpointProbe<'a> {
    Elastic(ThermalTensileElasticTetMidpointStep<'a>),
    Onset(ThermalTensileTetMidpointOnset<'a>),
    Activation(ThermalTensileTetMidpointActivation<'a>),
}

/// Disposable elastic prediction retaining a MEAN-force replay request.
#[derive(Debug)]
pub struct ThermalTensileTetMidpointOnset<'a> {
    // Private: a static endpoint-force replay must not consume interval means.
    onset: ThermalTensileTetOnset<'a>,
    inertia: RetainedThermalInertia,
    dt_s: f64,
}

/// Next coefficients activate tensile failure at held old mechanical strain.
/// There is no elastic endpoint predictor or physical activation time here.
#[derive(Debug)]
pub struct ThermalTensileTetMidpointActivation<'a> {
    activation: ThermalTensileTetActivation<'a>,
    inertia: RetainedThermalInertia,
    dt_s: f64,
}

#[derive(Debug)]
pub struct ThermalBrittleTetMotionHandoff<'a> {
    pub motion: ThermalBrittleTetMotion<'a>,
    pub elastic_storage_j: f64,
    pub storage_defect_j: f64,
    pub maximum_stress_defect_pa: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ThermalTensileTetDynamicError {
    Invalid(&'static str),
    Initiation(ThermalTensileTetInitiationError),
    Motion(BrittleTetDynamicError),
}

impl From<ThermalTensileTetInitiationError> for ThermalTensileTetDynamicError {
    fn from(error: ThermalTensileTetInitiationError) -> Self {
        Self::Initiation(error)
    }
}
impl From<BrittleTetDynamicError> for ThermalTensileTetDynamicError {
    fn from(error: BrittleTetDynamicError) -> Self {
        Self::Motion(error)
    }
}
impl From<BrittleTetError> for ThermalTensileTetDynamicError {
    fn from(error: BrittleTetError) -> Self {
        Self::Initiation(error.into())
    }
}
impl From<SmallStrainTetError> for ThermalTensileTetDynamicError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Initiation(error.into())
    }
}
impl core::fmt::Display for ThermalTensileTetDynamicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native thermal tensile motion refused: {self:?}")
    }
}
impl std::error::Error for ThermalTensileTetDynamicError {}

/// Algebraic mean of the existing source-bound elastic law, not another law.
struct ThermalElasticMean<'a> {
    owner: &'a ResolvedThermalBrittle,
    mechanical_start: Voigt,
}

impl SmallStrainLaw for ThermalElasticMean<'_> {
    type State = ();

    fn initial_state(&self) {}
    fn stress(&self, total_end: &Voigt, state: &()) -> Voigt {
        self.owner
            .mechanical_strain(total_end)
            .map(|end| {
                let mean = core::array::from_fn(|c| 0.5 * self.mechanical_start[c] + 0.5 * end[c]);
                self.owner.band().law().elastic().stress(&mean, state)
            })
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, total_end: &Voigt, state: &()) -> Tangent6 {
        self.owner
            .band()
            .law()
            .elastic()
            .tangent(total_end, state)
            .map(|row| row.map(|v| 0.5 * v))
    }
    fn update_state(&self, _end: &Voigt, _state: &()) {}
    fn admissibility(&self) -> MaterialAdmissibility {
        self.owner.band().law().elastic().admissibility()
    }
    fn card(&self) -> ModelCard {
        self.owner.band().law().elastic().card()
    }
}

impl<'a> ThermalTensileElasticTetMotion<'a> {
    /// Bind the original source density to the reference geometry, preserving
    /// the supplied elastic thermal preload and physical initial velocity.
    /// Later rho(T) changes cannot alter this fixed reference mass inventory.
    /// Starting observation performs no impulse or work and requires an elastic
    /// initial tensile state under the fixed supplied initiation settings.
    pub fn try_from_solid(
        solid: &ThermalTensileElasticTetSolid<'a>,
        initial_velocity_mps: &[[f64; 3]],
        initial_time_s: f64,
        maximum_steps: usize,
        initiation: TensileInitiationSettings,
        cx: &Cx<'_>,
    ) -> Result<Self, ThermalTensileTetDynamicError> {
        checkpoint(cx)?;
        if !initial_time_s.is_finite()
            || initial_time_s < 0.0
            || maximum_steps == 0
            || initial_velocity_mps.len() != solid.geometry.nodes_m().len()
            || initial_velocity_mps
                .iter()
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(ThermalTensileTetDynamicError::Invalid(
                "finite full initial velocity, nonnegative time and positive step budget required",
            ));
        }
        let mut densities = Vec::with_capacity(solid.materials.len());
        for (element, source) in solid.materials.iter().enumerate() {
            checkpoint(cx)?;
            let strain = mechanical(source, &solid.solution.element_strain[element], element)?;
            let detector = source
                .band()
                .tensile_initiator(initiation)
                .map_err(|source| ThermalTensileTetInitiationError::Element { element, source })?;
            if detector
                .between_controlled(&strain, &strain, None, &mut || cx.checkpoint().is_err())
                .map_err(|source| ThermalTensileTetInitiationError::Element { element, source })?
                .is_some()
            {
                return Err(ThermalTensileTetDynamicError::Invalid(
                    "initial thermal tensile state must be elastic",
                ));
            }
            densities.push(
                source
                    .band()
                    .reference_density_kg_m3()
                    .map_err(|source| BrittleTetDynamicError::Source { element, source })?,
            );
        }
        let mass = P1ReferenceMass::new(
            solid
                .geometry
                .tetrahedra()
                .iter()
                .zip(solid.geometry.element_volumes_m3())
                .map(|(&nodes, &volume)| (nodes, volume)),
            &densities,
            initial_velocity_mps.len() * 3,
            cx,
        )
        .map_err(SmallStrainTetError::from)?;
        let momentum = mass
            .action(initial_velocity_mps.as_flattened(), cx)
            .map_err(SmallStrainTetError::from)?;
        let mut kinetic_energy_j = 0.0;
        for (dof, (&velocity, &momentum)) in initial_velocity_mps
            .as_flattened()
            .iter()
            .zip(&momentum)
            .enumerate()
        {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            kinetic_energy_j = velocity.mul_add(momentum, kinetic_energy_j);
        }
        kinetic_energy_j *= 0.5;
        if !kinetic_energy_j.is_finite() || kinetic_energy_j < 0.0 {
            return Err(ThermalTensileTetDynamicError::Invalid(
                "representable nonnegative initial kinetic storage required",
            ));
        }
        checkpoint(cx)?;
        Ok(Self {
            solid: ThermalTensileElasticTetSolid {
                geometry: solid.geometry,
                materials: solid.materials,
                budget: solid.budget,
                solution: solid.solution.clone(),
                stored_energy_j: solid.stored_energy_j,
                supports: solid.supports.clone(),
            },
            inertia: RetainedThermalInertia {
                velocity_mps: initial_velocity_mps.to_vec(),
                reference_density_kg_m3: densities.into(),
                mass: Arc::new(mass),
                kinetic_energy_j,
                time_s: initial_time_s,
                maximum_steps,
                accepted_steps: 0,
            },
            initiation,
        })
    }

    pub fn solution(&self) -> &SmallStrainTetSolution<()> {
        &self.solid.solution
    }
    pub fn materials(&self) -> &[ResolvedThermalBrittle] {
        self.solid.materials
    }
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.inertia.velocity_mps
    }
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.inertia.reference_density_kg_m3
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.inertia.mass.reference_mass_kg()
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        self.inertia.kinetic_energy_j
    }
    pub fn stored_energy_j(&self) -> f64 {
        self.solid.stored_energy_j
    }
    pub fn time_s(&self) -> f64 {
        self.inertia.time_s
    }
    pub fn accepted_steps(&self) -> usize {
        self.inertia.accepted_steps
    }
    pub fn remaining_steps(&self) -> usize {
        self.inertia.maximum_steps - self.inertia.accepted_steps
    }

    /// Solve a centered-inertia elastic candidate and inspect its mechanical
    /// strain chord, or return an earlier coefficient-stage activation.
    /// Supports are absolute endpoints; forces are actual interval MEANS,
    /// including forces applied on moving supports. Plane hints are physical
    /// inputs for repeated principal tension, never inferred coordinate axes.
    ///
    /// Only an entirely elastic result advances time and consumes a step. Any
    /// trigger, source refusal, solve failure or late work refusal leaves this
    /// accepted motion usable. A refined interval requires a fresh native solve
    /// with its own explicitly supplied source, support and mean-force request.
    #[allow(clippy::too_many_arguments)]
    pub fn probe_midpoint<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetMidpointControls,
        direction_hints: Option<&[Option<[f64; 3]>]>,
        cx: &Cx<'_>,
    ) -> Result<ThermalTensileTetMidpointProbe<'next>, ThermalTensileTetDynamicError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        let old = &self.inertia;
        let next_time = old.admit_step(dt_s, controls)?;
        let (supports, activation) = self.solid.prepare_endpoint(
            next_materials,
            prescribed_m,
            mean_nodal_forces_n,
            settings,
            self.initiation,
            direction_hints,
            cx,
        )?;
        if let Some(activation) = activation {
            checkpoint(cx)?;
            return Ok(ThermalTensileTetMidpointProbe::Activation(
                ThermalTensileTetMidpointActivation {
                    activation,
                    inertia: old.clone(),
                    dt_s,
                },
            ));
        }
        let mut laws = Vec::with_capacity(next_materials.len());
        for (element, owner) in next_materials.iter().enumerate() {
            checkpoint(cx)?;
            laws.push(ThermalElasticMean {
                owner,
                mechanical_start: mechanical(
                    &self.solid.materials[element],
                    &self.solid.solution.element_strain[element],
                    element,
                )?,
            });
        }
        let problem = SmallStrainTetProblem {
            nodes_m: self.solid.geometry.nodes_m(),
            tetrahedra: self.solid.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m,
            nodal_forces_n: mean_nodal_forces_n,
            load_path: &[1.0],
            budget: self.solid.budget,
            settings,
        };
        // The elastic mean is linear with analytic C/2: no numerical fallback
        // is needed. The selected fallback remains available to brittle replay.
        let step = problem.advance_midpoint_with_mass(
            self.solid.solution.displacement_m.as_flattened(),
            &old.velocity_mps,
            &self.solid.solution.element_states,
            &old.mass,
            dt_s,
            None,
            cx,
        )?;
        let mut solution = step.solution;
        let mean_element_stress_pa = core::mem::take(&mut solution.element_stress);
        solution.element_stress = Vec::with_capacity(laws.len());
        for (element, law) in laws.iter().enumerate() {
            checkpoint(cx)?;
            let total = &solution.element_strain[element];
            if law.stress(total, &()) != mean_element_stress_pa[element] {
                return Err(ThermalTensileTetDynamicError::Invalid(
                    "native thermal elastic mean response mismatch",
                ));
            }
            let end = mechanical(law.owner, total, element)?;
            solution
                .element_stress
                .push(law.owner.band().law().elastic().stress(&end, &()));
        }
        let endpoint = self.solid.inspect_endpoint(
            next_materials,
            solution,
            prescribed_m,
            mean_nodal_forces_n,
            settings,
            self.initiation,
            direction_hints,
            cx,
        )?;
        let solution = match &endpoint {
            ThermalTensileEndpoint::Elastic(solution) => solution,
            ThermalTensileEndpoint::Onset(onset) => onset.elastic_predictor(),
        };
        let energy = self.admit_energy(
            next_materials,
            solution,
            &mean_element_stress_pa,
            step.energy,
            prescribed_m,
            mean_nodal_forces_n,
            controls,
            cx,
        )?;
        checkpoint(cx)?;
        match endpoint {
            ThermalTensileEndpoint::Onset(onset) => Ok(ThermalTensileTetMidpointProbe::Onset(
                ThermalTensileTetMidpointOnset {
                    onset,
                    inertia: old.clone(),
                    dt_s,
                },
            )),
            ThermalTensileEndpoint::Elastic(solution) => Ok(
                ThermalTensileTetMidpointProbe::Elastic(ThermalTensileElasticTetMidpointStep {
                    motion: ThermalTensileElasticTetMotion {
                        solid: ThermalTensileElasticTetSolid {
                            geometry: self.solid.geometry,
                            materials: next_materials,
                            budget: self.solid.budget,
                            solution,
                            stored_energy_j: energy.stored_after_j,
                            supports: Some(supports),
                        },
                        inertia: RetainedThermalInertia {
                            velocity_mps: step.velocity_mps,
                            reference_density_kg_m3: Arc::clone(&old.reference_density_kg_m3),
                            mass: Arc::clone(&old.mass),
                            kinetic_energy_j: step.energy.kinetic_after_j,
                            time_s: next_time,
                            maximum_steps: old.maximum_steps,
                            accepted_steps: old.accepted_steps + 1,
                        },
                        initiation: self.initiation,
                    },
                    mean_element_stress_pa,
                    energy,
                }),
            ),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_energy(
        &self,
        next_materials: &[ResolvedThermalBrittle],
        endpoint: &SmallStrainTetSolution<()>,
        mean_stress: &[Voigt],
        inertia: InertialEnergy,
        prescribed_m: &[(usize, f64)],
        mean_forces_n: &[(usize, f64)],
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMidpointEnergy, ThermalTensileTetDynamicError> {
        let material = self
            .solid
            .elastic_material_energy(next_materials, endpoint, cx)?;
        let mut mean_internal_work_j = 0.0;
        for (element, response) in mean_stress.iter().enumerate() {
            checkpoint(cx)?;
            let total_change = core::array::from_fn(|c| {
                endpoint.element_strain[element][c] - self.solid.solution.element_strain[element][c]
            });
            mean_internal_work_j = self.solid.geometry.element_volumes_m3()[element]
                .mul_add(contract(response, &total_change), mean_internal_work_j);
        }
        let old_u = self.solid.solution.displacement_m.as_flattened();
        let new_u = endpoint.displacement_m.as_flattened();
        let mut dead = 0.0;
        for &(dof, force) in mean_forces_n {
            checkpoint(cx)?;
            dead = force.mul_add(new_u[dof] - old_u[dof], dead);
        }
        let mut support = 0.0;
        for (&(dof, _), &reaction) in prescribed_m.iter().zip(&endpoint.increments[0].reactions_n) {
            checkpoint(cx)?;
            support = reaction.mul_add(new_u[dof] - old_u[dof], support);
        }
        let external = dead + support;
        let storage = material.stored_after_j - material.stored_before_j;
        let energy = ThermalBrittleTetMidpointEnergy {
            inertia,
            stored_before_j: material.stored_before_j,
            stored_after_j: material.stored_after_j,
            parameter_work_j: material.parameter_work_j,
            mechanical_work_j: material.mechanical_work_j,
            eigenstrain_work_j: material.eigenstrain_work_j,
            constitutive_strain_work_j: material.constitutive_strain_work_j,
            fracture_dissipation_j: 0.0,
            parameter_fracture_dissipation_j: 0.0,
            viscous_dissipation_j: 0.0,
            strain_split_error_j: material.strain_split_error_j,
            material_balance_error_j: material.material_balance_error_j,
            mean_internal_work_j,
            mean_work_error_j: mean_internal_work_j - material.mechanical_work_j,
            dead_load_mean_work_j: dead,
            support_mean_work_j: support,
            equilibrium_work_residual_j: (external - mean_internal_work_j)
                - inertia.inertial_work_j,
            total_balance_error_j: ((external
                + material.parameter_work_j
                + material.eigenstrain_work_j)
                - storage)
                - inertia.kinetic_change_j,
        };
        for (account, defect, tolerance) in [
            (
                "thermal elastic midpoint strain split",
                energy.strain_split_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "thermal elastic midpoint material work",
                energy.material_balance_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "thermal elastic midpoint mean response",
                energy.mean_work_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "thermal elastic midpoint equilibrium work",
                energy.equilibrium_work_residual_j,
                controls.equilibrium_work_tolerance_j,
            ),
            (
                "thermal elastic midpoint native kinetic identity",
                inertia.balance_error_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "thermal elastic midpoint retained kinetic storage",
                inertia.kinetic_before_j - self.inertia.kinetic_energy_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "thermal elastic midpoint complete work balance",
                energy.total_balance_error_j,
                controls.total_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        if inertia.numerical_loss_j != 0.0 {
            return Err(ThermalTensileTetDynamicError::Invalid(
                "thermal midpoint requires zero native numerical kinetic loss",
            ));
        }
        Ok(energy)
    }
}

impl RetainedThermalInertia {
    fn admit_step(
        &self,
        dt_s: f64,
        controls: BrittleTetMidpointControls,
    ) -> Result<f64, ThermalTensileTetDynamicError> {
        if self.accepted_steps >= self.maximum_steps {
            return Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: self.accepted_steps,
                maximum_steps: self.maximum_steps,
            }
            .into());
        }
        let next_time = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || !next_time.is_finite()
            || next_time <= self.time_s
            || [
                controls.material_energy_tolerance_j,
                controls.equilibrium_work_tolerance_j,
                controls.kinetic_work_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            || controls
                .difference_step_strain
                .is_some_and(|h| !h.is_finite() || h <= 0.0)
        {
            return Err(ThermalTensileTetDynamicError::Invalid(
                "representable positive midpoint time and finite nonnegative work gates required",
            ));
        }
        Ok(next_time)
    }

    fn handoff<'a, 'next>(
        &self,
        request: &ThermalTensileTetRequest<'a>,
        materials: &'next [ResolvedThermalBrittle],
        tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMotionHandoff<'next>, ThermalTensileTetDynamicError>
    where
        'a: 'next,
    {
        let handoff = request.handoff(materials, tolerance_j, cx)?;
        let motion = ThermalBrittleTetMotion::from_tensile_handoff(
            handoff.solid,
            self.velocity_mps.clone(),
            Arc::clone(&self.reference_density_kg_m3),
            Arc::clone(&self.mass),
            self.kinetic_energy_j,
            self.time_s,
            self.maximum_steps,
            self.accepted_steps,
            cx,
        )?;
        Ok(ThermalBrittleTetMotionHandoff {
            motion,
            elastic_storage_j: handoff.elastic_storage_j,
            storage_defect_j: handoff.storage_defect_j,
            maximum_stress_defect_pa: handoff.maximum_stress_defect_pa,
        })
    }

    fn replay<'a, 'next>(
        &self,
        request: &ThermalTensileTetRequest<'a>,
        dt_s: f64,
        bands: &'next ThermalTensileTetBands,
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMidpointStep<'next>, ThermalTensileTetDynamicError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        self.admit_step(dt_s, controls)?;
        request.admit_next_materials(&bands.next, cx)?;
        let handoff = self.handoff(
            request,
            &bands.initial,
            controls.material_energy_tolerance_j,
            cx,
        )?;
        Ok(handoff.motion.advance_midpoint(
            &bands.next,
            dt_s,
            &request.prescribed_end_m,
            &request.nodal_forces_end_n,
            request.solver_settings,
            controls,
            cx,
        )?)
    }
}

impl<'a> ThermalTensileTetMidpointOnset<'a> {
    pub fn events(&self) -> &[Option<ResolvedTensileOnset<'a>>] {
        self.onset.events()
    }
    pub fn initial_solution(&self) -> &SmallStrainTetSolution<()> {
        self.onset.initial_solution()
    }
    /// Disposable physical endpoint stress, not the interval-mean response.
    pub fn elastic_predictor(&self) -> &SmallStrainTetSolution<()> {
        self.onset.elastic_predictor()
    }
    /// Mechanical strain-chord coordinate only, never a physical event time.
    pub fn first_fraction_bracket(&self) -> [f64; 2] {
        self.onset.first_fraction_bracket()
    }
    pub fn first_candidate_element(&self) -> usize {
        self.onset.first_candidate_element()
    }
    pub fn requested_dt_s(&self) -> f64 {
        self.dt_s
    }
    pub fn physical_interval_s(&self) -> [f64; 2] {
        [self.inertia.time_s, self.inertia.time_s + self.dt_s]
    }
    pub fn prescribed_end_m(&self) -> &[(usize, f64)] {
        &self.onset.request.prescribed_end_m
    }
    pub fn mean_nodal_forces_n(&self) -> &[(usize, f64)] {
        &self.onset.request.nodal_forces_end_n
    }
    pub fn solver_settings(&self) -> IncrementSettings {
        self.onset.request.solver_settings
    }
    pub fn resolve_bands(
        &self,
        cx: &Cx<'_>,
    ) -> Result<ThermalTensileTetBands, ThermalTensileTetDynamicError> {
        Ok(self.onset.request.resolve_bands(cx)?)
    }
    /// Preserve the OLD accepted total/mechanical strain and original inertia.
    pub fn handoff<'next>(
        &self,
        materials: &'next [ResolvedThermalBrittle],
        tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMotionHandoff<'next>, ThermalTensileTetDynamicError>
    where
        'a: 'next,
    {
        self.inertia
            .handoff(&self.onset.request, materials, tolerance_j, cx)
    }
    /// Replay the original source/support/MEAN-force interval through native
    /// thermal brittle midpoint motion, consuming exactly one original step.
    pub fn replay<'next>(
        &self,
        bands: &'next ThermalTensileTetBands,
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMidpointStep<'next>, ThermalTensileTetDynamicError>
    where
        'a: 'next,
    {
        self.inertia
            .replay(&self.onset.request, self.dt_s, bands, controls, cx)
    }
}

impl<'a> ThermalTensileTetMidpointActivation<'a> {
    pub fn activations(&self) -> &[Option<ResolvedTensileActivation<'a>>] {
        self.activation.activations()
    }
    pub fn initial_solution(&self) -> &SmallStrainTetSolution<()> {
        self.activation.initial_solution()
    }
    pub fn requested_dt_s(&self) -> f64 {
        self.dt_s
    }
    pub fn physical_interval_s(&self) -> [f64; 2] {
        [self.inertia.time_s, self.inertia.time_s + self.dt_s]
    }
    pub fn prescribed_end_m(&self) -> &[(usize, f64)] {
        &self.activation.request.prescribed_end_m
    }
    pub fn mean_nodal_forces_n(&self) -> &[(usize, f64)] {
        &self.activation.request.nodal_forces_end_n
    }
    pub fn solver_settings(&self) -> IncrementSettings {
        self.activation.request.solver_settings
    }
    pub fn resolve_bands(
        &self,
        cx: &Cx<'_>,
    ) -> Result<ThermalTensileTetBands, ThermalTensileTetDynamicError> {
        Ok(self.activation.request.resolve_bands(cx)?)
    }
    /// Transfer the accepted OLD coefficients and motion; activation itself
    /// consumes no time, step, impulse, parameter work or fracture energy.
    pub fn handoff<'next>(
        &self,
        materials: &'next [ResolvedThermalBrittle],
        tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMotionHandoff<'next>, ThermalTensileTetDynamicError>
    where
        'a: 'next,
    {
        self.inertia
            .handoff(&self.activation.request, materials, tolerance_j, cx)
    }
    /// The native thermal owner performs coefficient/history transfer and the
    /// actual midpoint solve together with their separate physical work accounts.
    pub fn replay<'next>(
        &self,
        bands: &'next ThermalTensileTetBands,
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMidpointStep<'next>, ThermalTensileTetDynamicError>
    where
        'a: 'next,
    {
        self.inertia
            .replay(&self.activation.request, self.dt_s, bands, controls, cx)
    }
}
