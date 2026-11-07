//! Retained physical motion with sourced thermal strain and tensile fracture.
//!
//! One original consistent reference mass is bound to the initial material
//! snapshots and mesh, then cached across every thermal/mechanical endpoint.
//! Successor density queries do not change that inventory or impart momentum.
//! Thermal strain, displacement, velocity and native fracture history advance
//! together through the existing P1 backward-Euler Newton owner.
//!
//! The static and dynamic paths share coefficient transfer, mechanical-strain
//! conversion and exact native work recovery. Coefficients change first at held
//! mechanical history; total displacement and eigenstrain then follow one
//! straight increment under the new coefficients. The caller resolves thermal
//! and mechanical reversals. No intermediate restrained damage state is accepted.
//!
//! Eigenstrain work is a mechanical port, not caloric heat. Physical fracture,
//! surviving-bond viscous loss,
//! endpoint constitutive quadrature and backward-Euler kinetic loss have
//! distinct accounts. The default step is first-order backward Euler; the
//! midpoint path retains centered momentum and analytic mean forces. Thermal
//! coefficient splitting, small strain and fixed planes remain explicit.
//! Completed-fracture continuation retains this same thermal state while the
//! original topology, consistent mass and normal-impact owners handle fragments.

/// Completed thermally driven cuts with conserved native fragment motion.
pub mod fracture;
/// Work-conjugate thermal fracture with centered native momentum.
pub mod midpoint;
/// Frozen native damage and current thermal stiffness with original reference inertia.
pub mod modes;

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;
use fs_material::brittle::{BrittleEnergy, BrittleState};

use super::{ThermalBrittleResponse, ThermalBrittleTetError, ThermalBrittleTetSolid};
use crate::brittle::solid::dynamics::{BrittleTetDynamicControls, BrittleTetDynamicError};
use crate::brittle::solid::temperature::BrittleTetTemperatureElementEnergy;
use crate::brittle::solid::{checkpoint, energy_gate};
use crate::inertia::P1ReferenceMass;
use crate::plastic3::dynamics::InertialEnergy;
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    SmallStrainTetSolution,
};

/// Accepted thermal solid, velocity and immutable original reference inventory.
#[derive(Debug)]
pub struct ThermalBrittleTetMotion<'a> {
    solid: ThermalBrittleTetSolid<'a>,
    velocity_mps: Vec<[f64; 3]>,
    /// Actual imposed isotropic free-strain rate at this accepted endpoint.
    /// One value per original element; zero-time impulses and cuts retain it.
    free_linear_strain_rates_per_s: Arc<[f64]>,
    reference_density_kg_m3: Arc<[f64]>,
    mass: Arc<P1ReferenceMass>,
    kinetic_energy_j: f64,
    time_s: f64,
    maximum_steps: usize,
    accepted_steps: usize,
}

/// Whole-solid physical work and separate numerical accounts, in joules.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalBrittleTetDynamicEnergy {
    pub inertia: InertialEnergy,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub parameter_work_j: f64,
    /// Exact mean stress contracted with the TOTAL displacement-strain change.
    pub mechanical_work_j: f64,
    /// Minus mean stress contracted with the FREE thermal-strain change.
    pub eigenstrain_work_j: f64,
    pub constitutive_strain_work_j: f64,
    pub fracture_dissipation_j: f64,
    /// Exact physical surviving-bond loss on the mechanical-strain path,
    /// distinct from intrinsic fracture, thermal work and numerical damping.
    pub viscous_dissipation_j: f64,
    pub parameter_fracture_dissipation_j: f64,
    pub strain_split_error_j: f64,
    pub material_balance_error_j: f64,
    pub internal_endpoint_work_j: f64,
    pub dead_load_endpoint_work_j: f64,
    /// Full native reactions, including consistent-mass support coupling.
    pub support_endpoint_work_j: f64,
    /// Applied/support work minus endpoint internal work and inertial work.
    pub equilibrium_work_residual_j: f64,
    /// Endpoint internal work minus exact TOTAL-strain mechanical work.
    /// This excludes physical eigenstrain work and numerical kinetic loss.
    pub constitutive_endpoint_defect_j: f64,
    /// External plus parameter/eigenstrain work minus storage, fracture/viscous loss,
    /// kinetic change and the two distinct numerical accounts.
    pub total_balance_error_j: f64,
}

impl ThermalBrittleTetDynamicEnergy {
    #[must_use]
    pub fn external_endpoint_work_j(&self) -> f64 {
        self.dead_load_endpoint_work_j + self.support_endpoint_work_j
    }
}

/// Complete accepted endpoint; every refusal preserves the original motion.
#[derive(Debug)]
pub struct ThermalBrittleTetMotionStep<'a> {
    pub motion: ThermalBrittleTetMotion<'a>,
    /// Intrinsic mechanical-strain work/storage/fracture per reference volume,
    /// J/m3. Complete work and separate physical viscosity are in `energy`.
    pub element_energy: Vec<BrittleEnergy>,
    pub element_transfer: Vec<BrittleTetTemperatureElementEnergy>,
    pub energy: ThermalBrittleTetDynamicEnergy,
}

impl<'a> ThermalBrittleTetMotion<'a> {
    /// The tensile predecessor already bound the original reference inventory.
    /// Its material handoff checks the source, expansion, strain and storage;
    /// this internal transfer retains the same inertia and physical clock.
    /// Current source density may differ from the frozen reference density.
    #[allow(clippy::too_many_arguments)]
    pub(in crate::brittle::solid::expansion) fn from_tensile_handoff(
        mut solid: ThermalBrittleTetSolid<'a>,
        velocity_mps: Vec<[f64; 3]>,
        reference_density_kg_m3: Arc<[f64]>,
        mass: Arc<P1ReferenceMass>,
        kinetic_energy_j: f64,
        time_s: f64,
        maximum_steps: usize,
        accepted_steps: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleTetDynamicError> {
        checkpoint(cx)?;
        if velocity_mps.len() != solid.geometry.nodes_m().len()
            || mass.ndofs() != 3 * velocity_mps.len()
            || reference_density_kg_m3.len() != solid.materials.len()
            || accepted_steps > maximum_steps
            || maximum_steps == 0
        {
            return Err(BrittleTetDynamicError::Invalid(
                "thermal tensile handoff must preserve original motion and inventory dimensions",
            ));
        }
        let free_linear_strain_rates_per_s: Arc<[f64]> = vec![0.0; solid.materials.len()].into();
        solid.refresh_rate_stresses(&velocity_mps, &free_linear_strain_rates_per_s, cx)?;
        Ok(Self {
            solid,
            velocity_mps,
            free_linear_strain_rates_per_s,
            reference_density_kg_m3,
            mass,
            kinetic_energy_j,
            time_s,
            maximum_steps,
            accepted_steps,
        })
    }

    /// Retain an accepted thermal solid and bind its sourced density snapshots
    /// once to the supplied original reference geometry.
    ///
    /// Use the thermal dynamic resolver to retain density in the SAME selected
    /// source bundle. The caller declares this initial reference inventory:
    /// neither the expansion reference nor a later current-density query implies
    /// a volume conversion. Starting motion preserves all prestress, damage and
    /// fracture loss; supplied velocity is actual initial kinetic storage.
    /// Initial temperature is held: its free-strain rate is zero. Use
    /// [`Self::try_from_solid_with_free_strain_rates`] for an already moving
    /// thermal boundary; no thermal derivative is inferred from temperature.
    pub fn try_from_solid(
        solid: &ThermalBrittleTetSolid<'a>,
        initial_velocity_mps: &[[f64; 3]],
        initial_time_s: f64,
        maximum_steps: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleTetDynamicError> {
        Self::from_solid_rates(
            solid,
            initial_velocity_mps,
            None,
            initial_time_s,
            maximum_steps,
            cx,
        )
    }

    /// Admit an explicitly supplied instantaneous isotropic free-strain rate
    /// per element together with actual nodal velocity. The source can convert
    /// a declared temperature rate using its current selected alpha(T).
    ///
    /// Initial total stress uses `B*v - free_rate*I`; neither damage, work nor
    /// physical time advances. Reference mass and history remain unchanged.
    pub fn try_from_solid_with_free_strain_rates(
        solid: &ThermalBrittleTetSolid<'a>,
        initial_velocity_mps: &[[f64; 3]],
        initial_free_linear_strain_rates_per_s: &[f64],
        initial_time_s: f64,
        maximum_steps: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleTetDynamicError> {
        Self::from_solid_rates(
            solid,
            initial_velocity_mps,
            Some(initial_free_linear_strain_rates_per_s),
            initial_time_s,
            maximum_steps,
            cx,
        )
    }

    fn from_solid_rates(
        solid: &ThermalBrittleTetSolid<'a>,
        initial_velocity_mps: &[[f64; 3]],
        initial_free_rates: Option<&[f64]>,
        initial_time_s: f64,
        maximum_steps: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleTetDynamicError> {
        checkpoint(cx)?;
        if !initial_time_s.is_finite()
            || initial_time_s < 0.0
            || maximum_steps == 0
            || initial_velocity_mps.len() != solid.geometry.nodes_m().len()
            || initial_velocity_mps
                .iter()
                .flatten()
                .any(|v| !v.is_finite())
            || initial_free_rates.is_some_and(|rates| {
                rates.len() != solid.materials.len() || rates.iter().any(|rate| !rate.is_finite())
            })
        {
            return Err(BrittleTetDynamicError::Invalid(
                "finite full initial velocity/free-strain rates, nonnegative time and positive step budget required",
            ));
        }
        let mut densities = Vec::with_capacity(solid.materials.len());
        for (element, material) in solid.materials.iter().enumerate() {
            checkpoint(cx)?;
            densities.push(
                material
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
        let mut kinetic = 0.0;
        for (dof, (&v, &p)) in initial_velocity_mps
            .as_flattened()
            .iter()
            .zip(&momentum)
            .enumerate()
        {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            kinetic = v.mul_add(p, kinetic);
        }
        kinetic *= 0.5;
        if !kinetic.is_finite() || kinetic < 0.0 {
            return Err(BrittleTetDynamicError::Invalid(
                "initial consistent kinetic storage must be finite and nonnegative",
            ));
        }
        let mut retained = ThermalBrittleTetSolid {
            geometry: solid.geometry,
            materials: solid.materials,
            bands: Arc::clone(&solid.bands),
            released_faces: Arc::clone(&solid.released_faces),
            budget: solid.budget,
            solution: solid.solution.clone(),
            stored_energy_j: solid.stored_energy_j,
            fracture_dissipation_j: solid.fracture_dissipation_j,
        };
        let free_linear_strain_rates_per_s: Arc<[f64]> =
            initial_free_rates.map_or_else(|| vec![0.0; solid.materials.len()].into(), Arc::from);
        retained.refresh_rate_stresses(
            initial_velocity_mps,
            &free_linear_strain_rates_per_s,
            cx,
        )?;
        checkpoint(cx)?;
        Ok(Self {
            solid: retained,
            velocity_mps: initial_velocity_mps.to_vec(),
            free_linear_strain_rates_per_s,
            reference_density_kg_m3: densities.into(),
            mass: Arc::new(mass),
            kinetic_energy_j: kinetic,
            time_s: initial_time_s,
            maximum_steps,
            accepted_steps: 0,
        })
    }

    #[must_use]
    pub const fn solid(&self) -> &ThermalBrittleTetSolid<'a> {
        &self.solid
    }
    #[must_use]
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }
    /// Accepted imposed isotropic free-strain rates in original element order.
    /// A physical step declares a straight free-strain path, giving
    /// `(free_next-free_old)/dt`. Initial rates are explicitly supplied or zero;
    /// zero-time contact and topology transfer preserve them unchanged.
    #[must_use]
    pub fn free_linear_strain_rates_per_s(&self) -> &[f64] {
        &self.free_linear_strain_rates_per_s
    }
    #[must_use]
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    #[must_use]
    pub fn reference_mass_kg(&self) -> f64 {
        self.mass.reference_mass_kg()
    }
    #[must_use]
    pub const fn kinetic_energy_j(&self) -> f64 {
        self.kinetic_energy_j
    }
    #[must_use]
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    #[must_use]
    pub const fn accepted_steps(&self) -> usize {
        self.accepted_steps
    }
    #[must_use]
    pub const fn remaining_steps(&self) -> usize {
        self.maximum_steps - self.accepted_steps
    }

    /// Solve thermal expansion and inertial fracture in one physical interval.
    ///
    /// Successor sources retain the original reference, selected claims and
    /// fixed projection. Loads and support displacements are absolute endpoints.
    /// Each Newton proposal uses the same coefficient-transferred mechanical
    /// history. The reference mass stays unchanged even if queried rho(T) changes.
    /// Total and source-integrated free strain follow straight paths over dt;
    /// this declares the free-strain rate, not an unseen temperature trajectory.
    /// Every source, solve, budget, cancellation or late work refusal leaves
    /// the complete old material state, inventory, velocity and clock retryable.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn advance<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetDynamicControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMotionStep<'next>, BrittleTetDynamicError>
    where
        'a: 'next,
    {
        let next_time = self.admit_step(
            dt_s,
            &[
                controls.work.material_energy_tolerance_j,
                controls.work.equilibrium_work_tolerance_j,
                controls.work.maximum_endpoint_work_defect_j,
                controls.kinetic_work_tolerance_j,
                controls.maximum_numerical_kinetic_loss_j,
            ],
            prescribed_m.len(),
            nodal_forces_n.len(),
            cx,
        )?;
        let prepared = self
            .solid
            .prepare_materials(
                next_materials,
                controls.work.material_energy_tolerance_j,
                cx,
            )
            .map_err(thermal_error)?;
        let free_rates = self.interval_free_strain_rates(next_materials, dt_s, cx)?;
        let mut laws = Vec::with_capacity(next_materials.len());
        for (element, material) in next_materials.iter().enumerate() {
            checkpoint(cx)?;
            laws.push(
                material
                    .small_strain_law_over_time(dt_s)
                    .map_err(|source| BrittleTetDynamicError::Source { element, source })?,
            );
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
        let mut step = problem.advance_backward_euler_with_mass(
            self.solid.solution.displacement_m.as_flattened(),
            &self.velocity_mps,
            &prepared.states,
            &self.mass,
            dt_s,
            cx,
        )?;
        self.solid
            .require_released_faces(&step.solution.displacement_m, cx)?;
        let recovered = self
            .solid
            .recover_material(
                next_materials,
                &prepared,
                &mut step.solution,
                ThermalBrittleResponse::Endpoint,
                Some(dt_s),
                controls.work.material_energy_tolerance_j,
                cx,
            )
            .map_err(thermal_error)?;
        let material = recovered.energy;
        let (dead, support) =
            self.external_work(&step.solution, prescribed_m, nodal_forces_n, cx)?;
        let external = dead + support;
        let storage = material.stored_after_j - material.stored_before_j;
        let constitutive = material.endpoint_work_defect_j;
        let equilibrium =
            (external - recovered.response_internal_work_j) - step.energy.inertial_work_j;
        let total = ((((((external + material.parameter_work_j + material.eigenstrain_work_j)
            - storage)
            - material.fracture_dissipation_j)
            - material.viscous_dissipation_j)
            - step.energy.kinetic_change_j)
            - step.energy.numerical_loss_j)
            - constitutive;
        let energy = ThermalBrittleTetDynamicEnergy {
            inertia: step.energy,
            stored_before_j: material.stored_before_j,
            stored_after_j: material.stored_after_j,
            parameter_work_j: material.parameter_work_j,
            mechanical_work_j: material.mechanical_work_j,
            eigenstrain_work_j: material.eigenstrain_work_j,
            constitutive_strain_work_j: material.constitutive_strain_work_j,
            fracture_dissipation_j: material.fracture_dissipation_j,
            viscous_dissipation_j: material.viscous_dissipation_j,
            parameter_fracture_dissipation_j: material.parameter_fracture_dissipation_j,
            strain_split_error_j: material.strain_split_error_j,
            material_balance_error_j: material.material_balance_error_j,
            internal_endpoint_work_j: material.internal_endpoint_work_j,
            dead_load_endpoint_work_j: dead,
            support_endpoint_work_j: support,
            equilibrium_work_residual_j: equilibrium,
            constitutive_endpoint_defect_j: constitutive,
            total_balance_error_j: total,
        };
        for (account, defect, tolerance) in [
            (
                "thermal dynamic equilibrium work",
                equilibrium,
                controls.work.equilibrium_work_tolerance_j,
            ),
            (
                "thermal dynamic mechanical endpoint defect",
                constitutive,
                controls.work.maximum_endpoint_work_defect_j,
            ),
            (
                "thermal dynamic native kinetic identity",
                step.energy.balance_error_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "thermal dynamic retained kinetic storage",
                step.energy.kinetic_before_j - self.kinetic_energy_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "thermal backward-Euler numerical kinetic loss",
                step.energy.numerical_loss_j,
                controls.maximum_numerical_kinetic_loss_j,
            ),
            (
                "thermal dynamic complete work balance",
                total,
                controls.work.material_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        let mut solid = ThermalBrittleTetSolid {
            geometry: self.solid.geometry,
            materials: next_materials,
            bands: Arc::clone(&self.solid.bands),
            released_faces: Arc::clone(&self.solid.released_faces),
            budget: self.solid.budget,
            solution: step.solution,
            stored_energy_j: material.stored_after_j,
            fracture_dissipation_j: recovered.cumulative_loss_j,
        };
        solid.refresh_rate_stresses(&step.velocity_mps, &free_rates, cx)?;
        checkpoint(cx)?;
        Ok(ThermalBrittleTetMotionStep {
            motion: self.successor(
                solid,
                step.velocity_mps,
                free_rates,
                step.energy.kinetic_after_j,
                next_time,
            ),
            element_energy: recovered.element_energy,
            element_transfer: prepared.element_transfer,
            energy,
        })
    }

    /// One declared affine free-strain path per original material element.
    fn interval_free_strain_rates(
        &self,
        next_materials: &[ResolvedThermalBrittle],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<Arc<[f64]>, BrittleTetDynamicError> {
        checkpoint(cx)?;
        if next_materials.len() != self.solid.materials.len() || !dt_s.is_finite() || dt_s <= 0.0 {
            return Err(BrittleTetDynamicError::Invalid(
                "complete thermal source field and positive physical interval required",
            ));
        }
        let mut rates = Vec::with_capacity(next_materials.len());
        for (old, next) in self.solid.materials.iter().zip(next_materials) {
            checkpoint(cx)?;
            let change = next.free_linear_strain() - old.free_linear_strain();
            let rate = change / dt_s;
            if !rate.is_finite() || (change != 0.0 && rate == 0.0) {
                return Err(BrittleTetDynamicError::Invalid(
                    "affine free-strain rate must be finite and representable",
                ));
            }
            rates.push(rate);
        }
        Ok(rates.into())
    }

    fn admit_step(
        &self,
        dt_s: f64,
        tolerances: &[f64],
        prescribed_entries: usize,
        force_entries: usize,
        cx: &Cx<'_>,
    ) -> Result<f64, BrittleTetDynamicError> {
        checkpoint(cx)?;
        if self.accepted_steps >= self.maximum_steps {
            return Err(BrittleTetDynamicError::StepBudget {
                accepted_steps: self.accepted_steps,
                maximum_steps: self.maximum_steps,
            });
        }
        let next_time = self.time_s + dt_s;
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || !next_time.is_finite()
            || next_time <= self.time_s
            || tolerances.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(BrittleTetDynamicError::Invalid(
                "representable positive time and finite nonnegative work gates required",
            ));
        }
        let dofs = self.solid.solution.displacement_m.as_flattened().len();
        if prescribed_entries > dofs || force_entries > dofs {
            return Err(BrittleTetDynamicError::Invalid(
                "load entries exceed original DOF count",
            ));
        }
        Ok(next_time)
    }

    /// Native solve admission has already validated load DOFs and reaction order.
    fn external_work(
        &self,
        solution: &SmallStrainTetSolution<BrittleState<'_>>,
        prescribed_m: &[(usize, f64)],
        forces_n: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<(f64, f64), BrittleTetDynamicError> {
        let old = self.solid.solution.displacement_m.as_flattened();
        let next = solution.displacement_m.as_flattened();
        let mut dead = 0.0;
        for &(dof, force) in forces_n {
            checkpoint(cx)?;
            dead = force.mul_add(next[dof] - old[dof], dead);
        }
        let mut support = 0.0;
        for (&(dof, _), &reaction) in prescribed_m.iter().zip(&solution.increments[0].reactions_n) {
            checkpoint(cx)?;
            support = reaction.mul_add(next[dof] - old[dof], support);
        }
        Ok((dead, support))
    }

    /// Only a completely gated endpoint reaches this immutable publication step.
    fn successor<'next>(
        &self,
        solid: ThermalBrittleTetSolid<'next>,
        velocity_mps: Vec<[f64; 3]>,
        free_linear_strain_rates_per_s: Arc<[f64]>,
        kinetic_energy_j: f64,
        time_s: f64,
    ) -> ThermalBrittleTetMotion<'next>
    where
        'a: 'next,
    {
        ThermalBrittleTetMotion {
            solid,
            velocity_mps,
            free_linear_strain_rates_per_s,
            reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
            mass: Arc::clone(&self.mass),
            kinetic_energy_j,
            time_s,
            maximum_steps: self.maximum_steps,
            accepted_steps: self.accepted_steps + 1,
        }
    }
}

impl ThermalBrittleTetSolid<'_> {
    /// Observe actual endpoint mechanical rate without evolving any history.
    /// The native P1 total-strain projection and the explicitly imposed thermal
    /// rate remain separate until the original source performs their subtraction.
    fn refresh_rate_stresses(
        &mut self,
        velocity_mps: &[[f64; 3]],
        free_rates_per_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(), BrittleTetDynamicError> {
        checkpoint(cx)?;
        if velocity_mps.len() != self.geometry.nodes_m().len()
            || free_rates_per_s.len() != self.materials.len()
            || velocity_mps
                .iter()
                .flatten()
                .any(|value| !value.is_finite())
            || free_rates_per_s.iter().any(|value| !value.is_finite())
        {
            return Err(BrittleTetDynamicError::Invalid(
                "finite complete endpoint velocity and imposed free-strain rates required",
            ));
        }
        if !self
            .materials
            .iter()
            .zip(&self.solution.element_states)
            .any(|(source, state)| source.band().law().viscosity().is_some() && !state.failed())
        {
            return Ok(());
        }
        // Only native geometry/type admission is used; strain_field performs
        // no untimed constitutive trial. Every coordinate is known here, so
        // no free unknown is introduced into the caller's solve budget.
        let laws: Vec<_> = self
            .materials
            .iter()
            .map(|source| source.small_strain_law())
            .collect();
        let prescribed: Vec<_> = velocity_mps
            .as_flattened()
            .iter()
            .copied()
            .enumerate()
            .collect();
        let problem = SmallStrainTetProblem {
            nodes_m: self.geometry.nodes_m(),
            tetrahedra: self.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m: &prescribed,
            nodal_forces_n: &[],
            load_path: &[1.0],
            budget: self.budget,
            settings: IncrementSettings::default(),
        };
        let total_rates = problem.strain_field(velocity_mps.as_flattened(), cx)?;
        let mut stresses = Vec::with_capacity(self.materials.len());
        for (element, source) in self.materials.iter().enumerate() {
            checkpoint(cx)?;
            stresses.push(
                source
                    .stress_at_rate(
                        &self.solution.element_states[element],
                        &total_rates[element],
                        free_rates_per_s[element],
                    )
                    .map_err(|source| BrittleTetDynamicError::Source { element, source })?,
            );
        }
        checkpoint(cx)?;
        self.solution.element_stress = stresses;
        Ok(())
    }
}

fn thermal_error(error: ThermalBrittleTetError) -> BrittleTetDynamicError {
    match error {
        ThermalBrittleTetError::Source { element, source } => {
            BrittleTetDynamicError::Source { element, source }
        }
        ThermalBrittleTetError::Solid(error) => error.into(),
    }
}
