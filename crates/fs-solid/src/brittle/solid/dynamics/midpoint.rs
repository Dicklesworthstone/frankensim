//! Work-conjugate brittle motion on the original mass, source and Newton owners.
//!
//! The force response is the material's analytic mean stress over the complete
//! straight strain path, including rupture within a step. Native midpoint
//! inertia supplies `v_next = 2*(u_next-u_old)/dt-v_old`. Their work is checked
//! against the actual kinetic/storage changes and irreversible fracture loss,
//! without an endpoint quadrature correction or numerical kinetic dissipation.
//!
//! Mean response stresses are distinct from the accepted physical endpoint
//! stresses. Every endpoint is replayed through the original constitutive owner
//! before publication; the returned solid retains its real endpoint history.
//! Nonsmooth paths may lack a classical derivative. An explicit solver control
//! admits the existing bounded numerical derivative of the same mean response;
//! it never changes material evolution or supplies a substitute endpoint law.
//!
//! Loads are interval mean forces conjugate to the actual displacement change;
//! prescribed coordinates are endpoint displacements, with the same linear
//! path used by material work. The caller resolves temporal forcing, reversals
//! and fracture events. Exact path work does not remove frequency warping or
//! establish a temporal error bound. Fixed planes, small rotations/strains,
//! source conditions and reference inventory remain unchanged during a step.

/// Original thermal source/history transfer followed by work-conjugate motion.
pub mod temperature;

/// Adaptive native half-step trajectories with mean-force work and error control.
pub mod time;

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::Voigt;
use fs_material::brittle::BrittleEnergy;

use super::super::{BrittleTetError, BrittleTetSolid, checkpoint, energy_gate};
use super::{BrittleTetDynamicError, BrittleTetMotion};
use crate::plastic3::dynamics::InertialEnergy;
use crate::plastic3::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem};

/// Explicit whole-mesh energy gates and an optional numerical Newton fallback.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetMidpointControls {
    /// Native work/storage/fracture, cumulative loss, and mean-stress work gates.
    pub material_energy_tolerance_j: f64,
    /// Applied/support mean work minus mean internal and inertial work.
    pub equilibrium_work_tolerance_j: f64,
    /// Native kinetic-work identity and retained kinetic-storage continuity.
    pub kinetic_work_tolerance_j: f64,
    /// Applied/support work minus kinetic/storage changes and fracture loss.
    pub total_energy_tolerance_j: f64,
    /// `None` requires every needed constitutive tangent to be analytic.
    /// `Some(h)` explicitly admits native numerical differences only where the
    /// mean tangent is unavailable, with maximum tensor-strain perturbation h.
    /// At most seven point-stress evaluations per such element/assembly;
    /// invalid, out-of-domain or unrepresentable stencils still refuse.
    pub difference_step_strain: Option<f64>,
}

/// Actual work and physical loss, with arithmetic/solve defects kept separate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetMidpointEnergy {
    /// Native consistent-mass receipt; numerical kinetic loss is zero.
    pub inertia: InertialEnergy,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub material_work_j: f64,
    pub fracture_dissipation_j: f64,
    pub material_balance_error_j: f64,
    /// Volume integral of mean stress double-contracted with actual strain change.
    pub mean_internal_work_j: f64,
    /// Mean-stress work minus the independently retained material work account.
    pub mean_work_error_j: f64,
    pub dead_load_mean_work_j: f64,
    /// Complete native mean reactions, including acceleration of constrained DOFs.
    pub support_mean_work_j: f64,
    pub equilibrium_work_residual_j: f64,
    /// External work minus elastic change, fracture loss and kinetic change.
    pub total_balance_error_j: f64,
}

impl BrittleTetMidpointEnergy {
    #[must_use]
    pub fn external_mean_work_j(&self) -> f64 {
        self.dead_load_mean_work_j + self.support_mean_work_j
    }
}

/// Complete accepted endpoint; no partial history or physical time is published.
#[derive(Debug)]
pub struct BrittleTetMidpointStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    /// Original per-element work, storage change and fracture loss in J/m³.
    pub element_energy: Vec<BrittleEnergy>,
    /// Work-conjugate response used by Newton, distinct from endpoint stresses
    /// available through `motion.solid().solution().element_stress`.
    pub mean_element_stress_pa: Vec<Voigt>,
    pub energy: BrittleTetMidpointEnergy,
    /// Number of element tangent evaluations using the explicitly admitted
    /// numerical fallback. Repeated Newton iterates count separately.
    pub numerical_tangent_evaluations: usize,
}

impl<'a> BrittleTetMotion<'a> {
    /// Advance using analytic material path work and native centered momentum.
    ///
    /// Material sources, frozen reference mass, crack planes and width remain
    /// the original owners. Prescribed displacements are absolute endpoints;
    /// supplied nodal loads are interval MEAN forces. Returned support reactions
    /// are means too, while displacement, velocity, stresses and material state
    /// in the retained solid are physical endpoints. Every trial uses the same
    /// accepted history. Cancellation or a late energy refusal leaves `self`
    /// unchanged, including its clock, step allowance and cumulative fracture.
    pub fn advance_midpoint(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMidpointStep<'a>, BrittleTetDynamicError> {
        let step = self.midpoint_candidate(
            dt_s,
            prescribed_m,
            mean_nodal_forces_n,
            settings,
            controls,
            cx,
        )?;
        self.solid
            .require_released_faces(&step.motion.solid.solution.displacement_m, cx)?;
        checkpoint(cx)?;
        Ok(step)
    }

    /// Unpublished response for the native sustained-contact force solve.
    /// Only released-face opening admission is deferred: signed force probes
    /// may overlap before complementarity is solved. Every material, source,
    /// history, inertia and work gate remains identical to public midpoint.
    /// The receiving contact owner must admit all openings before publication.
    pub(in crate::brittle::solid::dynamics) fn midpoint_candidate(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMidpointStep<'a>, BrittleTetDynamicError> {
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
            || [
                controls.material_energy_tolerance_j,
                controls.equilibrium_work_tolerance_j,
                controls.kinetic_work_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
            || controls
                .difference_step_strain
                .is_some_and(|h| !h.is_finite() || h <= 0.0)
        {
            return Err(BrittleTetDynamicError::Invalid(
                "representable positive time, finite work gates and admitted strain difference required",
            ));
        }
        let old_u = self.solid.solution.displacement_m.as_flattened();
        if prescribed_m.len() > old_u.len() || mean_nodal_forces_n.len() > old_u.len() {
            return Err(BrittleTetDynamicError::Invalid(
                "midpoint load entries exceed the original DOF count",
            ));
        }
        let mut laws = Vec::with_capacity(self.solid.materials.len());
        for material in self.solid.materials {
            checkpoint(cx)?;
            laws.push(material.law().mean_small_strain_law());
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
        let step = problem.advance_midpoint_with_mass(
            old_u,
            &self.velocity_mps,
            &self.solid.solution.element_states,
            &self.mass,
            dt_s,
            controls.difference_step_strain,
            cx,
        )?;
        let mut solution = step.solution;
        let mean_element_stress_pa = core::mem::take(&mut solution.element_stress);
        solution.element_stress = Vec::with_capacity(self.solid.materials.len());
        let mut element_energy = Vec::with_capacity(self.solid.materials.len());
        let mut stored_after = 0.0;
        let mut material_work = 0.0;
        let mut fracture_loss = 0.0;
        let mut cumulative_loss = 0.0;
        let mut mean_work = 0.0;
        for (element, material) in self.solid.materials.iter().enumerate() {
            checkpoint(cx)?;
            let old = &self.solid.solution.element_states[element];
            let path = material
                .law()
                .trial_with_mean_stress(&solution.element_strain[element], old)
                .map_err(|source| BrittleTetError::Material { element, source })?;
            let endpoint = path.endpoint;
            if path.mean_stress_pa != mean_element_stress_pa[element]
                || endpoint.state != solution.element_states[element]
                || endpoint.state.damage() < old.damage()
                || endpoint.state.dissipated_j_m3() < old.dissipated_j_m3()
                || endpoint.energy.dissipated < 0.0
            {
                return Err(BrittleTetDynamicError::Invalid(
                    "midpoint response, endpoint history or irreversibility mismatch",
                ));
            }
            let volume = self.solid.bands[element].reference_volume_m3();
            stored_after = volume.mul_add(endpoint.state.stored_energy_j_m3(), stored_after);
            material_work = volume.mul_add(endpoint.energy.work, material_work);
            fracture_loss = volume.mul_add(endpoint.energy.dissipated, fracture_loss);
            cumulative_loss = volume.mul_add(endpoint.state.dissipated_j_m3(), cumulative_loss);
            mean_work = volume.mul_add(path.path_work_j_m3, mean_work);
            solution.element_stress.push(endpoint.stress_pa);
            element_energy.push(endpoint.energy);
        }
        let new_u = solution.displacement_m.as_flattened();
        let mut dead_work = 0.0;
        for &(dof, force) in mean_nodal_forces_n {
            checkpoint(cx)?;
            dead_work = force.mul_add(new_u[dof] - old_u[dof], dead_work);
        }
        let mut support_work = 0.0;
        for (&(dof, _), &reaction) in prescribed_m.iter().zip(&solution.increments[0].reactions_n) {
            checkpoint(cx)?;
            support_work = reaction.mul_add(new_u[dof] - old_u[dof], support_work);
        }
        let external = dead_work + support_work;
        let storage_change = stored_after - self.solid.stored_energy_j;
        let energy = BrittleTetMidpointEnergy {
            inertia: step.energy,
            stored_before_j: self.solid.stored_energy_j,
            stored_after_j: stored_after,
            material_work_j: material_work,
            fracture_dissipation_j: fracture_loss,
            material_balance_error_j: (material_work - storage_change) - fracture_loss,
            mean_internal_work_j: mean_work,
            mean_work_error_j: mean_work - material_work,
            dead_load_mean_work_j: dead_work,
            support_mean_work_j: support_work,
            equilibrium_work_residual_j: (external - mean_work) - step.energy.inertial_work_j,
            total_balance_error_j: ((external - storage_change) - fracture_loss)
                - step.energy.kinetic_change_j,
        };
        for (account, defect, tolerance) in [
            (
                "midpoint material work/storage/fracture",
                energy.material_balance_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "midpoint cumulative fracture loss",
                (cumulative_loss - self.solid.fracture_dissipation_j) - fracture_loss,
                controls.material_energy_tolerance_j,
            ),
            (
                "analytic mean-stress work",
                energy.mean_work_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "midpoint support/dead-load equilibrium work",
                energy.equilibrium_work_residual_j,
                controls.equilibrium_work_tolerance_j,
            ),
            (
                "midpoint native kinetic work",
                step.energy.balance_error_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "midpoint retained kinetic storage",
                step.energy.kinetic_before_j - self.kinetic_energy_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "complete midpoint brittle work balance",
                energy.total_balance_error_j,
                controls.total_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        if step.energy.numerical_loss_j != 0.0
            || [stored_after, fracture_loss, cumulative_loss]
                .iter()
                .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(BrittleTetDynamicError::Invalid(
                "midpoint must retain physical nonnegative energy without numerical kinetic loss",
            ));
        }
        let solid = BrittleTetSolid {
            geometry: self.solid.geometry,
            materials: self.solid.materials,
            bands: Arc::clone(&self.solid.bands),
            released_faces: Arc::clone(&self.solid.released_faces),
            budget: self.solid.budget,
            solution,
            stored_energy_j: stored_after,
            fracture_dissipation_j: cumulative_loss,
        };
        checkpoint(cx)?;
        Ok(BrittleTetMidpointStep {
            motion: Self {
                solid,
                velocity_mps: step.velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::clone(&self.mass),
                kinetic_energy_j: step.energy.kinetic_after_j,
                time_s: next_time,
                maximum_steps: self.maximum_steps,
                accepted_steps: self.accepted_steps + 1,
            },
            element_energy,
            mean_element_stress_pa,
            energy,
            numerical_tangent_evaluations: step.numerical_tangent_evaluations,
        })
    }
}
