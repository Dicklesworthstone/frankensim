//! Elastic midpoint motion and source-bound onset without numerical kinetic loss.
//!
//! The original isotropic elastic backbone evaluated at the mean strain gives
//! the exact work-conjugate response on the supplied straight strain path. Its
//! endpoint derivative is C/2. Native centered inertia and the original Newton
//! owner solve the candidate; true ENDPOINT stresses are restored before the
//! shared tensile detector and elastic-storage observation inspect it.
//!
//! Inputs are interval MEAN forces and endpoint prescribed displacements. A
//! separate onset type preserves this force meaning through native brittle
//! midpoint replay. It transfers the same old accepted mass, velocity, kinetic
//! storage, clock and consumed step budget. No predictor is accepted as damaged
//! motion, no chord fraction becomes an event time, and no kinetic damping or
//! endpoint quadrature defect is introduced into the physical loss account.
//!
//! Exact supplied-path work does not resolve unseen reversals, remove frequency
//! warping or certify temporal error. The source-bound tensile-only, fixed-plane,
//! small-strain assumptions and explicit refinement responsibility remain.

use super::{
    BrittleTetMotionHandoff, RetainedTensileInertia, TensileElasticTetMotion,
    TensileTetDynamicError, TensileTetMotionOnset, TensileTetProbe,
};
use crate::brittle::solid::dynamics::BrittleTetDynamicError;
use crate::brittle::solid::dynamics::midpoint::{
    BrittleTetMidpointControls, BrittleTetMidpointEnergy, BrittleTetMidpointStep,
};
use crate::brittle::solid::{checkpoint, energy_gate};
use crate::inertia::InertialEnergy;
use crate::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, SmallStrainTetSolution,
};
use fs_evidence::ModelCard;
use fs_exec::Cx;
use fs_material::brittle::source::ResolvedBrittleCrackBand;
use fs_material::brittle::source::initiation::ResolvedTensileOnset;
use fs_material::tensor::contract;
use fs_material::{IsotropicElastic, MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use std::sync::Arc;

/// Complete elastic endpoint with a distinct interval-mean response and work.
#[derive(Debug)]
pub struct TensileElasticTetMidpointStep<'a> {
    pub motion: TensileElasticTetMotion<'a>,
    pub mean_element_stress_pa: Vec<Voigt>,
    pub energy: BrittleTetMidpointEnergy,
}

#[derive(Debug)]
pub enum TensileTetMidpointProbe<'a> {
    Elastic(TensileElasticTetMidpointStep<'a>),
    Onset(TensileTetMidpointOnset<'a>),
}

/// Disposable midpoint prediction with a retained MEAN-force replay request.
#[derive(Debug)]
pub struct TensileTetMidpointOnset<'a> {
    // Deliberately private: the endpoint-force replay method is not exposed.
    original: TensileTetMotionOnset<'a>,
}

/// A borrowed algebraic view of the original state-free linear elastic law.
/// No new coefficients, material history, source evidence or evolution law.
struct ElasticMean<'a> {
    owner: &'a IsotropicElastic,
    start: Voigt,
}

impl SmallStrainLaw for ElasticMean<'_> {
    type State = ();

    fn initial_state(&self) {}
    fn stress(&self, end: &Voigt, state: &()) -> Voigt {
        self.owner.stress(
            &core::array::from_fn(|c| 0.5 * self.start[c] + 0.5 * end[c]),
            state,
        )
    }
    fn tangent(&self, end: &Voigt, state: &()) -> Tangent6 {
        self.owner
            .tangent(end, state)
            .map(|row| row.map(|v| 0.5 * v))
    }
    fn update_state(&self, _end: &Voigt, _state: &()) {}
    fn admissibility(&self) -> MaterialAdmissibility {
        self.owner.admissibility()
    }
    fn card(&self) -> ModelCard {
        self.owner.card()
    }
}

impl<'a> TensileElasticTetMotion<'a> {
    /// Advance a native centered elastic candidate and inspect its strain chord.
    ///
    /// Loads are physical interval MEANS, including forces on moving supports;
    /// support coordinates are absolute endpoints. The original support set,
    /// sources, mass and initial history remain fixed. Elastic response always
    /// uses its analytic C/2 tangent; numerical fallback is only needed, when
    /// explicitly selected, by a subsequent nonsmooth brittle replay.
    ///
    /// A crossing returns the old motion and exact original mean-load request.
    /// A smaller physical interval must be solved again with its own supplied
    /// mean forces and support targets. This strain-chord detector is not a
    /// continuous-time event certificate or proof against unresolved reversals.
    #[allow(clippy::too_many_arguments)]
    pub fn probe_midpoint(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleTetMidpointControls,
        direction_hints: Option<&[Option<[f64; 3]>]>,
        cx: &Cx<'_>,
    ) -> Result<TensileTetMidpointProbe<'a>, TensileTetDynamicError> {
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
            return Err(TensileTetDynamicError::Invalid(
                "representable positive midpoint time and finite nonnegative work gates required",
            ));
        }
        let supports =
            self.solid
                .admit_endpoint(prescribed_m, mean_nodal_forces_n, direction_hints)?;
        let mut laws = Vec::with_capacity(self.solid.materials.len());
        for (source, &start) in self
            .solid
            .materials
            .iter()
            .zip(&self.solid.solution.element_strain)
        {
            checkpoint(cx)?;
            laws.push(ElasticMean {
                owner: source.law().elastic(),
                start,
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
            let strain = &solution.element_strain[element];
            if law.stress(strain, &()) != mean_element_stress_pa[element] {
                return Err(TensileTetDynamicError::Invalid(
                    "native elastic midpoint response mismatch",
                ));
            }
            solution.element_stress.push(law.owner.stress(strain, &()));
        }
        // The shared private request slots retain interval MEANS here. Neither
        // the enclosed static solid nor its quasistatic replay is exposed by
        // motion wrappers; an endpoint force must never be inferred from them.
        let (probe, stored_after_j) = self.solid.inspect_endpoint(
            solution,
            supports,
            prescribed_m,
            mean_nodal_forces_n,
            settings,
            self.initiation,
            direction_hints,
            cx,
        )?;
        let endpoint = match &probe {
            TensileTetProbe::Elastic(solid) => &solid.solution,
            TensileTetProbe::Onset(onset) => &onset.predictor,
        };
        let energy = self.admit_midpoint_energy(
            endpoint,
            &mean_element_stress_pa,
            stored_after_j,
            step.energy,
            prescribed_m,
            mean_nodal_forces_n,
            controls,
            cx,
        )?;
        checkpoint(cx)?;
        match probe {
            TensileTetProbe::Elastic(solid) => Ok(TensileTetMidpointProbe::Elastic(
                TensileElasticTetMidpointStep {
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
                    mean_element_stress_pa,
                    energy,
                },
            )),
            TensileTetProbe::Onset(onset) => {
                Ok(TensileTetMidpointProbe::Onset(TensileTetMidpointOnset {
                    original: TensileTetMotionOnset {
                        onset,
                        inertia: old.clone(),
                        dt_s,
                    },
                }))
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_midpoint_energy(
        &self,
        endpoint: &SmallStrainTetSolution<()>,
        mean_stress: &[Voigt],
        stored_after_j: f64,
        inertia: InertialEnergy,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMidpointEnergy, TensileTetDynamicError> {
        let mut material_work_j = 0.0;
        let mut mean_internal_work_j = 0.0;
        for (element, response) in mean_stress.iter().enumerate() {
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
            mean_internal_work_j =
                volume.mul_add(contract(response, &increment), mean_internal_work_j);
        }
        let old_u = self.solid.solution.displacement_m.as_flattened();
        let new_u = endpoint.displacement_m.as_flattened();
        let mut dead_load_mean_work_j = 0.0;
        for &(dof, force) in mean_nodal_forces_n {
            checkpoint(cx)?;
            dead_load_mean_work_j = force.mul_add(new_u[dof] - old_u[dof], dead_load_mean_work_j);
        }
        let mut support_mean_work_j = 0.0;
        for (&(dof, _), &reaction) in prescribed_m.iter().zip(&endpoint.increments[0].reactions_n) {
            checkpoint(cx)?;
            support_mean_work_j = reaction.mul_add(new_u[dof] - old_u[dof], support_mean_work_j);
        }
        let external = dead_load_mean_work_j + support_mean_work_j;
        let storage_change = stored_after_j - self.solid.stored_energy_j;
        let energy = BrittleTetMidpointEnergy {
            inertia,
            stored_before_j: self.solid.stored_energy_j,
            stored_after_j,
            material_work_j,
            fracture_dissipation_j: 0.0,
            viscous_dissipation_j: 0.0,
            material_balance_error_j: material_work_j - storage_change,
            mean_internal_work_j,
            mean_work_error_j: mean_internal_work_j - material_work_j,
            dead_load_mean_work_j,
            support_mean_work_j,
            equilibrium_work_residual_j: (external - mean_internal_work_j)
                - inertia.inertial_work_j,
            total_balance_error_j: (external - storage_change) - inertia.kinetic_change_j,
        };
        for (account, defect, tolerance) in [
            (
                "elastic midpoint work/storage",
                energy.material_balance_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "elastic midpoint mean response work",
                energy.mean_work_error_j,
                controls.material_energy_tolerance_j,
            ),
            (
                "elastic midpoint support/dead-load work",
                energy.equilibrium_work_residual_j,
                controls.equilibrium_work_tolerance_j,
            ),
            (
                "elastic midpoint native kinetic identity",
                inertia.balance_error_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "elastic midpoint retained kinetic storage",
                inertia.kinetic_before_j - self.inertia.kinetic_energy_j,
                controls.kinetic_work_tolerance_j,
            ),
            (
                "complete elastic midpoint work balance",
                energy.total_balance_error_j,
                controls.total_energy_tolerance_j,
            ),
        ] {
            energy_gate(account, defect, tolerance)?;
        }
        if inertia.numerical_loss_j != 0.0 {
            return Err(TensileTetDynamicError::Invalid(
                "midpoint may not introduce numerical kinetic loss",
            ));
        }
        Ok(energy)
    }
}

impl<'a> TensileTetMidpointOnset<'a> {
    /// Original source diagnostics on the strain chord, never physical event times.
    pub fn events(&self) -> &[Option<ResolvedTensileOnset<'a>>] {
        self.original.events()
    }
    pub fn initial_solution(&self) -> &SmallStrainTetSolution<()> {
        self.original.initial_solution()
    }
    /// Physical endpoint stresses recovered from the disposable elastic solve.
    pub fn elastic_predictor(&self) -> &SmallStrainTetSolution<()> {
        self.original.elastic_predictor()
    }
    pub fn requested_dt_s(&self) -> f64 {
        self.original.requested_dt_s()
    }
    pub fn requested_interval_s(&self) -> [f64; 2] {
        self.original.requested_interval_s()
    }
    pub fn prescribed_end_m(&self) -> &[(usize, f64)] {
        self.original.prescribed_end_m()
    }
    /// Exact ordered interval means, including duplicates and loads on supports.
    pub fn mean_nodal_forces_n(&self) -> &[(usize, f64)] {
        &self.original.onset.nodal_forces_end_n
    }
    pub fn solver_settings(&self) -> IncrementSettings {
        self.original.solver_settings()
    }
    pub fn resolve_bands(
        &self,
        cx: &Cx<'_>,
    ) -> Result<Vec<ResolvedBrittleCrackBand>, TensileTetDynamicError> {
        self.original.resolve_bands(cx)
    }

    /// Native source/stress/storage continuity with the SAME old motion inventory.
    pub fn handoff<'next>(
        &self,
        materials: &'next [ResolvedBrittleCrackBand],
        energy_tolerance_j: f64,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMotionHandoff<'next>, TensileTetDynamicError>
    where
        'a: 'next,
    {
        self.original.handoff(materials, energy_tolerance_j, cx)
    }

    /// Re-solve the original interval using native brittle MEAN-force motion.
    /// Source seeding and every work/history gate precede publication. A refused
    /// replay leaves the old motion and original request available unchanged.
    pub fn replay<'next>(
        &self,
        materials: &'next [ResolvedBrittleCrackBand],
        handoff_energy_tolerance_j: f64,
        controls: BrittleTetMidpointControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetMidpointStep<'next>, TensileTetDynamicError>
    where
        'a: 'next,
    {
        self.handoff(materials, handoff_energy_tolerance_j, cx)?
            .motion
            .advance_midpoint(
                self.requested_dt_s(),
                self.prescribed_end_m(),
                self.mean_nodal_forces_n(),
                self.solver_settings(),
                controls,
                cx,
            )
            .map_err(Into::into)
    }
}
