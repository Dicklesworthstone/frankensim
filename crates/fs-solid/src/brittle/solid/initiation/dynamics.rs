//! Native elastic inertia, tensile prediction and conservative fracture replay.
//!
//! Backward Euler supplies the actual candidate displacement and velocity using
//! the original consistent reference mass. The source-bound Rankine detector
//! inspects the affine STRAIN CHORD between accepted and candidate endpoints.
//! An entirely elastic chord admits the native time step; a crossing returns a
//! disposable predictor. Its chord fractions are never interpreted as physical
//! event times. Time resolution must resolve reversals and inertial oscillations;
//! this first-order endpoint model does not certify the continuous trajectory.
//!
//! At a crossing, the existing geometry/source owner selects fixed crack bands
//! and seeds them at the OLD accepted strain. The same mass, velocity, kinetic
//! storage, clock and consumed step budget transfer without an impulse. Replay
//! solves the original physical interval with the original nonlinear brittle
//! motion owner; the over-strength elastic predictor is never accepted.

/// Exact elastic mean force and centered inertia, with native midpoint replay.
pub mod midpoint;

use std::sync::Arc;

use super::{TensileElasticTetSolid, TensileTetInitiationError, TensileTetOnset, TensileTetProbe};
use crate::brittle::solid::dynamics::{
    BrittleTetDynamicControls, BrittleTetDynamicEnergy, BrittleTetDynamicError, BrittleTetMotion,
    BrittleTetMotionStep,
};
use crate::brittle::solid::{BrittleTetError, checkpoint, energy_gate};
use crate::inertia::{InertialEnergy, P1ReferenceMass};
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    SmallStrainTetSolution,
};
use fs_exec::Cx;
use fs_material::brittle::initiation::TensileInitiationSettings;
use fs_material::brittle::source::ResolvedBrittleCrackBand;
use fs_material::brittle::source::initiation::ResolvedTensileOnset;
use fs_material::tensor::contract;

/// Private accepted inertia; cloning retains the SAME immutable mass operator.
#[derive(Clone, Debug)]
struct RetainedTensileInertia {
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    mass: Arc<P1ReferenceMass>,
    kinetic_energy_j: f64,
    time_s: f64,
    maximum_steps: usize,
    accepted_steps: usize,
}

/// Accepted tensile-only elastic motion before fixed crack planes are selected.
#[derive(Debug)]
pub struct TensileElasticTetMotion<'a> {
    solid: TensileElasticTetSolid<'a>,
    inertia: RetainedTensileInertia,
    initiation: TensileInitiationSettings,
}

/// Complete elastic step; physical fracture dissipation is exactly zero.
#[derive(Debug)]
pub struct TensileElasticTetMotionStep<'a> {
    pub motion: TensileElasticTetMotion<'a>,
    /// Endpoint elastic quadrature and backward-Euler kinetic loss remain
    /// separate numerical accounts, never fracture energy or thermal input.
    pub energy: BrittleTetDynamicEnergy,
}

#[derive(Debug)]
pub enum TensileTetMotionProbe<'a> {
    Elastic(TensileElasticTetMotionStep<'a>),
    Onset(TensileTetMotionOnset<'a>),
}

/// Disposable native elastic prediction retaining the original replay request.
#[derive(Debug)]
pub struct TensileTetMotionOnset<'a> {
    onset: TensileTetOnset<'a>,
    inertia: RetainedTensileInertia,
    dt_s: f64,
}

/// The original accepted motion under native brittle constitutive ownership.
#[derive(Debug)]
pub struct BrittleTetMotionHandoff<'a> {
    pub motion: BrittleTetMotion<'a>,
    pub elastic_storage_j: f64,
    pub storage_defect_j: f64,
    pub maximum_stress_defect_pa: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum TensileTetDynamicError {
    Invalid(&'static str),
    Initiation(TensileTetInitiationError),
    Motion(BrittleTetDynamicError),
}

impl From<TensileTetInitiationError> for TensileTetDynamicError {
    fn from(error: TensileTetInitiationError) -> Self {
        Self::Initiation(error)
    }
}
impl From<BrittleTetDynamicError> for TensileTetDynamicError {
    fn from(error: BrittleTetDynamicError) -> Self {
        Self::Motion(error)
    }
}
impl From<BrittleTetError> for TensileTetDynamicError {
    fn from(error: BrittleTetError) -> Self {
        Self::Initiation(error.into())
    }
}
impl From<SmallStrainTetError> for TensileTetDynamicError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Initiation(error.into())
    }
}
impl core::fmt::Display for TensileTetDynamicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native inertial tensile initiation refused: {self:?}")
    }
}
impl std::error::Error for TensileTetDynamicError {}

impl<'a> TensileElasticTetMotion<'a> {
    /// Bind source density to the original geometry and retain elastic preload.
    ///
    /// The accepted strain must remain strictly below the supplied tensile
    /// criterion. Numerical initiation settings then stay fixed for this motion.
    /// Initial velocity is actual physical state; starting observation neither
    /// applies an impulse nor charges its kinetic/elastic storage as new work.
    pub fn try_from_solid(
        solid: &TensileElasticTetSolid<'a>,
        initial_velocity_mps: &[[f64; 3]],
        initial_time_s: f64,
        maximum_steps: usize,
        initiation: TensileInitiationSettings,
        cx: &Cx<'_>,
    ) -> Result<Self, TensileTetDynamicError> {
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
            return Err(TensileTetDynamicError::Invalid(
                "finite full initial velocity, nonnegative time and positive step budget required",
            ));
        }
        let mut densities = Vec::with_capacity(solid.materials.len());
        for (element, source) in solid.materials.iter().enumerate() {
            checkpoint(cx)?;
            let detector = source
                .tensile_initiator(initiation)
                .map_err(|source| TensileTetInitiationError::Element { element, source })?;
            let strain = &solid.solution.element_strain[element];
            if detector
                .between_controlled(strain, strain, None, &mut || cx.checkpoint().is_err())
                .map_err(|source| TensileTetInitiationError::Element { element, source })?
                .is_some()
            {
                return Err(TensileTetDynamicError::Invalid(
                    "initial tensile state must be elastic",
                ));
            }
            densities.push(
                source
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
        for (dof, (&v, &p)) in initial_velocity_mps
            .as_flattened()
            .iter()
            .zip(&momentum)
            .enumerate()
        {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            kinetic_energy_j = v.mul_add(p, kinetic_energy_j);
        }
        kinetic_energy_j *= 0.5;
        if !kinetic_energy_j.is_finite() || kinetic_energy_j < 0.0 {
            return Err(TensileTetDynamicError::Invalid(
                "representable nonnegative initial kinetic storage required",
            ));
        }
        checkpoint(cx)?;
        Ok(Self {
            solid: TensileElasticTetSolid {
                geometry: solid.geometry,
                materials: solid.materials,
                budget: solid.budget,
                solution: solid.solution.clone(),
                stored_energy_j: solid.stored_energy_j,
                supports: solid.supports.clone(),
                nodal_forces_n: solid.nodal_forces_n.clone(),
            },
            inertia: RetainedTensileInertia {
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

    /// Solve the native elastic inertial endpoint, then inspect its strain chord.
    ///
    /// Prescribed DOFs remain the same as in the accepted elastic history;
    /// their values and all dead forces are actual absolute endpoint inputs.
    /// Free bodies need no artificial restraints. Repeated principal tension
    /// requires physical direction hints, exactly as in the native detector.
    ///
    /// A crossing retains the OLD motion and exact original request for replay.
    /// No linear interpolation of displacement or chord fraction is accepted
    /// as an inertial endpoint or physical event time. To refine a prediction,
    /// call probe again with a smaller explicitly supplied physical interval
    /// and its corresponding original loads/support values.
    #[allow(clippy::too_many_arguments)]
    pub fn probe(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetDynamicControls,
        direction_hints: Option<&[Option<[f64; 3]>]>,
        cx: &Cx<'_>,
    ) -> Result<TensileTetMotionProbe<'a>, TensileTetDynamicError> {
        checkpoint(cx)?;
        let old = &self.inertia;
        if old.accepted_steps >= old.maximum_steps {
            return Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: old.accepted_steps,
                maximum_steps: old.maximum_steps,
            }
            .into());
        }
        let next_time = old.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || !next_time.is_finite()
            || next_time <= old.time_s
            || [
                controls.work.material_energy_tolerance_j,
                controls.work.equilibrium_work_tolerance_j,
                controls.work.maximum_endpoint_work_defect_j,
                controls.kinetic_work_tolerance_j,
                controls.maximum_numerical_kinetic_loss_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(TensileTetDynamicError::Invalid(
                "representable positive time and finite nonnegative work gates required",
            ));
        }
        let supports = self
            .solid
            .admit_endpoint(prescribed_m, nodal_forces_n, direction_hints)?;
        let mut laws = Vec::with_capacity(self.solid.materials.len());
        for source in self.solid.materials {
            checkpoint(cx)?;
            laws.push(source.law().elastic());
        }
        let problem = SmallStrainTetProblem {
            nodes_m: self.solid.geometry.nodes_m(),
            tetrahedra: self.solid.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m,
            nodal_forces_n,
            load_path: &[1.0],
            budget: self.solid.budget,
            settings,
        };
        let step = problem.advance_backward_euler_with_mass(
            self.solid.solution.displacement_m.as_flattened(),
            &old.velocity_mps,
            &self.solid.solution.element_states,
            &old.mass,
            dt_s,
            cx,
        )?;
        let (probe, stored_after_j) = self.solid.inspect_endpoint(
            step.solution,
            supports,
            prescribed_m,
            nodal_forces_n,
            settings,
            self.initiation,
            direction_hints,
            cx,
        )?;
        let endpoint = match &probe {
            TensileTetProbe::Elastic(solid) => &solid.solution,
            TensileTetProbe::Onset(onset) => &onset.predictor,
        };
        let energy = self.admit_energy(
            endpoint,
            stored_after_j,
            step.energy,
            prescribed_m,
            nodal_forces_n,
            controls,
            cx,
        )?;
        checkpoint(cx)?;
        match probe {
            TensileTetProbe::Elastic(solid) => Ok(TensileTetMotionProbe::Elastic(
                TensileElasticTetMotionStep {
                    motion: Self {
                        solid,
                        inertia: RetainedTensileInertia {
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
                    energy,
                },
            )),
            TensileTetProbe::Onset(onset) => {
                Ok(TensileTetMotionProbe::Onset(TensileTetMotionOnset {
                    onset,
                    inertia: old.clone(),
                    dt_s,
                }))
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_energy(
        &self,
        endpoint: &SmallStrainTetSolution<()>,
        stored_after_j: f64,
        inertia: InertialEnergy,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        controls: BrittleTetDynamicControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetDynamicEnergy, TensileTetDynamicError> {
        let mut material_work_j = 0.0;
        let mut internal_endpoint_work_j = 0.0;
        for element in 0..self.solid.materials.len() {
            checkpoint(cx)?;
            let increment = core::array::from_fn(|c| {
                endpoint.element_strain[element][c] - self.solid.solution.element_strain[element][c]
            });
            let mean = core::array::from_fn(|c| {
                0.5 * self.solid.solution.element_stress[element][c]
                    + 0.5 * endpoint.element_stress[element][c]
            });
            let volume = self.solid.geometry.element_volumes_m3()[element];
            material_work_j = volume.mul_add(contract(&mean, &increment), material_work_j);
            internal_endpoint_work_j = volume.mul_add(
                contract(&endpoint.element_stress[element], &increment),
                internal_endpoint_work_j,
            );
        }
        let old_u = self.solid.solution.displacement_m.as_flattened();
        let new_u = endpoint.displacement_m.as_flattened();
        let mut dead_load_endpoint_work_j = 0.0;
        for &(dof, force) in nodal_forces_n {
            checkpoint(cx)?;
            dead_load_endpoint_work_j =
                force.mul_add(new_u[dof] - old_u[dof], dead_load_endpoint_work_j);
        }
        let mut support_endpoint_work_j = 0.0;
        for (&(dof, _), &reaction) in prescribed_m.iter().zip(&endpoint.increments[0].reactions_n) {
            checkpoint(cx)?;
            support_endpoint_work_j =
                reaction.mul_add(new_u[dof] - old_u[dof], support_endpoint_work_j);
        }
        let external = dead_load_endpoint_work_j + support_endpoint_work_j;
        let storage_change = stored_after_j - self.solid.stored_energy_j;
        let constitutive_endpoint_defect_j = internal_endpoint_work_j - material_work_j;
        let energy = BrittleTetDynamicEnergy {
            inertia,
            stored_before_j: self.solid.stored_energy_j,
            stored_after_j,
            material_work_j,
            fracture_dissipation_j: 0.0,
            material_balance_error_j: material_work_j - storage_change,
            internal_endpoint_work_j,
            dead_load_endpoint_work_j,
            support_endpoint_work_j,
            equilibrium_work_residual_j: (external - internal_endpoint_work_j)
                - inertia.inertial_work_j,
            constitutive_endpoint_defect_j,
            total_balance_error_j: (((external - storage_change) - inertia.kinetic_change_j)
                - inertia.numerical_loss_j)
                - constitutive_endpoint_defect_j,
        };
        for (account, defect, tolerance) in [
            (
                "elastic inertial work/storage",
                energy.material_balance_error_j,
                controls.work.material_energy_tolerance_j,
            ),
            (
                "elastic inertial support/dead-load work",
                energy.equilibrium_work_residual_j,
                controls.work.equilibrium_work_tolerance_j,
            ),
            (
                "elastic endpoint quadrature",
                constitutive_endpoint_defect_j,
                controls.work.maximum_endpoint_work_defect_j,
            ),
            (
                "elastic native kinetic work identity",
                inertia.balance_error_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "elastic retained kinetic storage",
                inertia.kinetic_before_j - self.inertia.kinetic_energy_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "elastic backward-Euler kinetic loss",
                inertia.numerical_loss_j,
                controls.maximum_numerical_kinetic_loss_j,
            ),
            (
                "complete elastic kinetic work balance",
                energy.total_balance_error_j,
                controls.work.material_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        Ok(energy)
    }
}

impl<'a> TensileTetMotionOnset<'a> {
    /// Source/plane diagnostics along the candidate STRAIN CHORD, not time events.
    pub fn events(&self) -> &[Option<ResolvedTensileOnset<'a>>] {
        &self.onset.events
    }
    pub fn initial_solution(&self) -> &SmallStrainTetSolution<()> {
        &self.onset.initial
    }
    /// Disposable original elastic response; no over-strength state is accepted.
    pub fn elastic_predictor(&self) -> &SmallStrainTetSolution<()> {
        &self.onset.predictor
    }
    pub fn requested_dt_s(&self) -> f64 {
        self.dt_s
    }
    pub fn requested_interval_s(&self) -> [f64; 2] {
        [self.inertia.time_s, self.inertia.time_s + self.dt_s]
    }
    pub fn prescribed_end_m(&self) -> &[(usize, f64)] {
        &self.onset.prescribed_end_m
    }
    pub fn nodal_forces_end_n(&self) -> &[(usize, f64)] {
        &self.onset.nodal_forces_end_n
    }
    pub fn solver_settings(&self) -> IncrementSettings {
        self.onset.solver_settings
    }

    /// Native projected widths, discovered planes and unchanged source bundles.
    pub fn resolve_bands(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Vec<ResolvedBrittleCrackBand>, TensileTetDynamicError> {
        self.onset.resolve_bands(cx).map_err(Into::into)
    }

    /// Transfer OLD accepted motion without replacing inertia or resetting time.
    ///
    /// Existing source-bound seeding proves strain/stress/storage continuity and
    /// refuses prior mixed-mode initiation. The exact original mass operator and
    /// velocity are retained, so neither momentum nor kinetic energy changes.
    pub fn handoff<'next>(
        &self,
        materials: &'next [ResolvedBrittleCrackBand],
        energy_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMotionHandoff<'next>, TensileTetDynamicError>
    where
        'a: 'next,
    {
        let handoff = self.onset.handoff(materials, energy_tolerance_j, cx)?;
        let old = &self.inertia;
        let motion = BrittleTetMotion::from_tensile_handoff(
            handoff.solid,
            old.velocity_mps.clone(),
            Arc::clone(&old.reference_density_kg_m3),
            Arc::clone(&old.mass),
            old.kinetic_energy_j,
            old.time_s,
            old.maximum_steps,
            old.accepted_steps,
            cx,
        )?;
        Ok(BrittleTetMotionHandoff {
            motion,
            elastic_storage_j: handoff.elastic_storage_j,
            storage_defect_j: handoff.storage_defect_j,
            maximum_stress_defect_pa: handoff.maximum_stress_defect_pa,
        })
    }

    /// Re-solve the exact original requested interval with native brittle motion.
    ///
    /// No predictor state, interpolation, new force ordering or synthetic
    /// support reaction enters the replay. Successful replay consumes exactly
    /// one original time step; any late refusal leaves this prediction reusable.
    pub fn replay<'next>(
        &self,
        materials: &'next [ResolvedBrittleCrackBand],
        handoff_energy_tolerance_j: f64,
        controls: BrittleTetDynamicControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMotionStep<'next>, TensileTetDynamicError>
    where
        'a: 'next,
    {
        let handoff = self.handoff(materials, handoff_energy_tolerance_j, cx)?;
        handoff
            .motion
            .advance(
                self.dt_s,
                &self.onset.prescribed_end_m,
                &self.onset.nodal_forces_end_n,
                self.onset.solver_settings,
                controls,
                cx,
            )
            .map_err(Into::into)
    }
}
