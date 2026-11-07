//! Work-conjugate thermal expansion and fracture with native centered inertia.
//!
//! The same source/history preparation and native recovery serve static,
//! backward-Euler and midpoint endpoints. The material's analytic path-mean
//! stress drives the existing P1 midpoint solver; physical endpoint stresses
//! are recovered before publication. Thermal eigenstrain work remains distinct
//! from displacement work, coefficient work and irreversible fracture loss.
//!
//! Mean external/support work plus parameter/eigenstrain work balances elastic
//! storage, physical fracture/viscous loss and native kinetic change, without endpoint
//! quadrature compensation or numerical kinetic loss. This work identity does
//! not make the coefficient-first thermal split second order. The caller must
//! resolve temperature changes, force means, thermal/mechanical reversals and
//! fracture transitions. Small-strain, fixed-plane/reference-inventory limits
//! remain; no heat, entropy or fragment-contact model is inferred.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::brittle::BrittleEnergy;
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;

use super::super::{ThermalBrittleResponse, ThermalBrittleTetSolid};
use super::{ThermalBrittleTetMotion, thermal_error};
use crate::brittle::solid::dynamics::BrittleTetDynamicError;
use crate::brittle::solid::dynamics::midpoint::BrittleTetMidpointControls;
use crate::brittle::solid::temperature::BrittleTetTemperatureElementEnergy;
use crate::brittle::solid::{checkpoint, energy_gate};
use crate::plastic3::dynamics::InertialEnergy;
use crate::plastic3::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem};

/// Physical mean-work accounts in joules; native numerical kinetic loss is zero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalBrittleTetMidpointEnergy {
    pub inertia: InertialEnergy,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub parameter_work_j: f64,
    /// Exact mean stress contracted with TOTAL displacement-strain change.
    pub mechanical_work_j: f64,
    /// Minus exact mean stress contracted with FREE thermal-strain change.
    pub eigenstrain_work_j: f64,
    pub constitutive_strain_work_j: f64,
    pub fracture_dissipation_j: f64,
    /// Exact surviving-bond physical loss on the mechanical-strain path.
    /// No part is counted as intrinsic fracture or numerical kinetic damping.
    pub viscous_dissipation_j: f64,
    pub parameter_fracture_dissipation_j: f64,
    pub strain_split_error_j: f64,
    pub material_balance_error_j: f64,
    /// Actual assembled mean response contracted with TOTAL strain change.
    pub mean_internal_work_j: f64,
    /// Solver mean work minus independently recovered total-strain work.
    pub mean_work_error_j: f64,
    pub dead_load_mean_work_j: f64,
    /// Complete native mean reactions, including support inertia coupling.
    pub support_mean_work_j: f64,
    pub equilibrium_work_residual_j: f64,
    /// External plus coefficient/eigenstrain work minus elastic storage change,
    /// physical fracture/viscous loss and kinetic change. No numerical loss is subtracted.
    pub total_balance_error_j: f64,
}

impl ThermalBrittleTetMidpointEnergy {
    #[must_use]
    pub fn external_mean_work_j(&self) -> f64 {
        self.dead_load_mean_work_j + self.support_mean_work_j
    }
}

/// Complete native thermal/inertial endpoint with unaveraged fracture history.
#[derive(Debug)]
pub struct ThermalBrittleTetMidpointStep<'a> {
    pub motion: ThermalBrittleTetMotion<'a>,
    /// Intrinsic mechanical-strain work/storage/fracture densities; total work
    /// and physical viscous loss are reported separately in `energy`.
    pub element_energy: Vec<BrittleEnergy>,
    pub element_transfer: Vec<BrittleTetTemperatureElementEnergy>,
    pub energy: ThermalBrittleTetMidpointEnergy,
    /// Original bounded numerical constitutive tangent evaluations, when
    /// explicitly enabled by the caller for unavailable classical derivatives.
    pub numerical_tangent_evaluations: usize,
}

impl<'a> ThermalBrittleTetMotion<'a> {
    /// Advance one coupled thermal/mechanical interval with work-conjugate forces.
    ///
    /// Applied forces are actual INTERVAL MEANS, conjugate to displacement
    /// change; prescribed displacements are absolute endpoints. Each native
    /// Newton trial retains the same coefficient-transferred mechanical history
    /// and the original reference mass. The mean material tangent differentiates
    /// the mean force with respect to the NEW total strain at fixed temperature.
    /// The optional native numerical fallback differentiates that SAME response.
    /// Total and free thermal strain share one affine path. Published endpoint
    /// stress uses the actual `B*v_next` minus that imposed free-strain rate;
    /// centered endpoint velocity generally differs from the path's mean rate.
    ///
    /// Every source, solver, time/budget, cancellation and physical-work gate
    /// passes before the new material owner, velocity and clock are published.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_midpoint<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMidpointStep<'next>, BrittleTetDynamicError>
    where
        'a: 'next,
    {
        let step = self.midpoint_candidate(
            next_materials,
            dt_s,
            prescribed_m,
            mean_forces_n,
            settings,
            controls,
            cx,
        )?;
        self.solid
            .require_released_faces(&step.motion.solid.solution.displacement_m, cx)?;
        checkpoint(cx)?;
        Ok(step)
    }

    /// Unpublished native response for the sustained-contact force solve.
    /// Only released-face opening admission is deferred: signed force probes
    /// may overlap while the original contact owner solves complementarity.
    /// Every thermal source, mechanical history, inertia and work gate remains
    /// identical to public midpoint, including an explicitly supplied source
    /// transition. The receiving contact owner must admit all actual openings
    /// and physical force/work accounts before publishing a successor.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::brittle::solid::expansion::dynamics) fn midpoint_candidate<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetMidpointStep<'next>, BrittleTetDynamicError>
    where
        'a: 'next,
    {
        let next_time = self.admit_step(
            dt_s,
            &[
                controls.material_energy_tolerance_j,
                controls.equilibrium_work_tolerance_j,
                controls.kinetic_work_tolerance_j,
                controls.total_energy_tolerance_j,
            ],
            prescribed_m.len(),
            mean_forces_n.len(),
            cx,
        )?;
        if controls
            .difference_step_strain
            .is_some_and(|h| !h.is_finite() || h <= 0.0)
        {
            return Err(BrittleTetDynamicError::Invalid(
                "positive finite numerical strain step required",
            ));
        }
        let prepared = self
            .solid
            .prepare_materials(next_materials, controls.material_energy_tolerance_j, cx)
            .map_err(thermal_error)?;
        let free_rates = self.interval_free_strain_rates(next_materials, dt_s, cx)?;
        let mut laws = Vec::with_capacity(next_materials.len());
        for (element, material) in next_materials.iter().enumerate() {
            checkpoint(cx)?;
            laws.push(
                material
                    .mean_small_strain_law_over_time(dt_s)
                    .map_err(|source| BrittleTetDynamicError::Source { element, source })?,
            );
        }
        let problem = SmallStrainTetProblem {
            nodes_m: self.solid.geometry.nodes_m(),
            tetrahedra: self.solid.geometry.tetrahedra(),
            laws: SmallStrainLawField::PerElement(&laws),
            prescribed_m,
            nodal_forces_n: mean_forces_n,
            load_path: &[1.0],
            budget: self.solid.budget,
            settings,
        };
        let mut step = problem.advance_midpoint_with_mass(
            self.solid.solution.displacement_m.as_flattened(),
            &self.velocity_mps,
            &prepared.states,
            &self.mass,
            dt_s,
            controls.difference_step_strain,
            cx,
        )?;
        let recovered = self
            .solid
            .recover_material(
                next_materials,
                &prepared,
                &mut step.solution,
                ThermalBrittleResponse::Mean,
                Some(dt_s),
                controls.material_energy_tolerance_j,
                cx,
            )
            .map_err(thermal_error)?;
        let material = recovered.energy;
        let (dead, support) =
            self.external_work(&step.solution, prescribed_m, mean_forces_n, cx)?;
        let external = dead + support;
        let storage = material.stored_after_j - material.stored_before_j;
        let equilibrium =
            (external - recovered.response_internal_work_j) - step.energy.inertial_work_j;
        let total = ((((external + material.parameter_work_j + material.eigenstrain_work_j)
            - storage)
            - material.fracture_dissipation_j)
            - material.viscous_dissipation_j)
            - step.energy.kinetic_change_j;
        let energy = ThermalBrittleTetMidpointEnergy {
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
            mean_internal_work_j: recovered.response_internal_work_j,
            mean_work_error_j: recovered.response_internal_work_j - material.mechanical_work_j,
            dead_load_mean_work_j: dead,
            support_mean_work_j: support,
            equilibrium_work_residual_j: equilibrium,
            total_balance_error_j: total,
        };
        for (account, defect, tolerance) in [
            (
                "thermal midpoint native mean work",
                energy.mean_work_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "thermal midpoint equilibrium work",
                equilibrium,
                controls.equilibrium_work_tolerance_j,
            ),
            (
                "thermal midpoint native kinetic identity",
                step.energy.balance_error_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "thermal midpoint retained kinetic storage",
                step.energy.kinetic_before_j - self.kinetic_energy_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "thermal midpoint complete work balance",
                total,
                controls.total_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        if step.energy.numerical_loss_j != 0.0 {
            return Err(BrittleTetDynamicError::Invalid(
                "thermal midpoint requires zero native numerical kinetic loss",
            ));
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
        Ok(ThermalBrittleTetMidpointStep {
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
            numerical_tangent_evaluations: step.numerical_tangent_evaluations,
        })
    }
}
