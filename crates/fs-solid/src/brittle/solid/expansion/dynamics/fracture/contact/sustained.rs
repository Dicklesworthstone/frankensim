//! Sustained normal contact for thermally fractured native components.
//!
//! The original nonlinear force/configuration solve is shared with nonthermal
//! fracture. Each disposable response uses the native thermal midpoint, its
//! original consistent mass and the same accepted mechanical history. Only the
//! final freshly replayed nonpenetrating response can advance physical time.
//!
//! Current-source dwell and an explicitly supplied thermal source transition
//! use this same force owner. Source coefficients, thermal strain and normal
//! forces are solved through the original native material path; no intermediate
//! restrained damage, gap projection, penalty spring or artificial impulse is
//! accepted. Eigenstrain and parameter work retain their mechanical meaning.

use fs_contact::sustained::solve::SolveReport;
use fs_exec::Cx;
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;

use crate::brittle::solid::dynamics::fracture::contact::BrittleFaceImpactRow;
pub use crate::brittle::solid::dynamics::fracture::contact::sustained::{
    BrittleFaceSustainedControls, BrittleFaceSustainedError,
};
use crate::brittle::solid::dynamics::fracture::contact::sustained::{
    NativeFaceSustainedTrial, SustainedScope, SustainedWork, solve_released_face_sustained,
};
use crate::brittle::solid::expansion::dynamics::ThermalBrittleTetMotion;
use crate::brittle::solid::expansion::dynamics::midpoint::{
    ThermalBrittleTetMidpointEnergy, ThermalBrittleTetMidpointStep,
};
use crate::plastic3::IncrementSettings;

/// Accepted thermal motion with native contact forces and separate work ports.
#[derive(Debug)]
pub struct ThermalBrittleFaceSustainedStep<'a> {
    pub motion: ThermalBrittleTetMotion<'a>,
    /// External dead-load and native support work, plus the original thermal
    /// parameter/eigenstrain ports. Internal normal work is gated separately.
    pub energy: ThermalBrittleTetMidpointEnergy,
    /// Original thermal midpoint receipt including applied contact nodal loads.
    pub native_energy: ThermalBrittleTetMidpointEnergy,
    pub solve: SolveReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub normal_contact_work_j: f64,
    pub absolute_contact_work_j: f64,
    pub contact_work_balance_error_j: f64,
    pub absolute_contact_couple_n_m: f64,
    pub numerical_tangent_evaluations: usize,
}

impl<'a> ThermalBrittleTetMotion<'a> {
    /// Resolve normal dwell or unloaded separation at the current thermal state.
    ///
    /// Every non-tied released pair participates. Loaded pairs must remain
    /// matched, active and normally at rest at both native endpoints. An
    /// approaching impact, open-to-loaded arrival, finite sliding or unresolved
    /// normal work refuses. The predecessor remains available for retry.
    ///
    /// Applied forces are interval means; prescribed total displacements are
    /// absolute endpoints. Original thermal sources, integrated expansion
    /// reference, tensile damage and frozen reference matter remain attached to
    /// the returned motion. Contact forces never supply a new material law or
    /// a conversion of dissipated fracture energy into heat.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_sustained(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceSustainedStep<'a>, BrittleFaceSustainedError> {
        self.advance_released_face_sustained_from(
            dt_s,
            prescribed_m,
            mean_nodal_forces_n,
            None,
            settings,
            controls,
            cx,
        )
    }

    /// Couple an explicitly resolved thermal endpoint to sustained normal force.
    ///
    /// The original thermal owner admits `next_materials`, transfers source
    /// coefficients at held mechanical history, and solves total displacement,
    /// eigenstrain and contact force together. Every signed-force trial replays
    /// this same predecessor and source successor; discarded probes cannot
    /// accept additional thermal damage or change the reference mass inventory.
    ///
    /// Native total and mechanical strains remain distinct. The complete
    /// external-only balance includes signed parameter and eigenstrain work,
    /// physical storage/fracture and kinetic change; contact work must pass its
    /// own numerical allowance. Temperature-dependent density cannot add matter.
    ///
    /// This retains the original coefficient-first, straight thermal/strain
    /// increment. Resolve temperature reversals and source knots explicitly;
    /// no automatic thermal subdivision, caloric heat, finite sliding or global
    /// continuous-collision guarantee is supplied by the force solve.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_sustained_at_materials<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceSustainedStep<'next>, BrittleFaceSustainedError>
    where
        'a: 'next,
    {
        self.sustained_at_materials_from(
            next_materials,
            dt_s,
            prescribed_m,
            mean_nodal_forces_n,
            None,
            settings,
            controls,
            cx,
        )
    }

    /// Previous normal forces seed the numerical solve, never material history.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_released_face_sustained_from(
        &self,
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        initial_normal_forces_n: Option<&[f64]>,
        settings: IncrementSettings,
        controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceSustainedStep<'a>, BrittleFaceSustainedError> {
        self.sustained_at_materials_from(
            self.solid.materials,
            dt_s,
            prescribed_m,
            mean_nodal_forces_n,
            initial_normal_forces_n,
            settings,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn sustained_at_materials_from<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        dt_s: f64,
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        initial_normal_forces_n: Option<&[f64]>,
        settings: IncrementSettings,
        controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceSustainedStep<'next>, BrittleFaceSustainedError>
    where
        'a: 'next,
    {
        let step = solve_released_face_sustained(
            self.released_face_motion_view(),
            prescribed_m,
            mean_nodal_forces_n,
            initial_normal_forces_n,
            SustainedScope::Complete,
            controls,
            cx,
            |total| {
                self.midpoint_candidate(
                    next_materials,
                    dt_s,
                    prescribed_m,
                    total,
                    settings,
                    controls.work,
                    cx,
                )
            },
        )?;
        Ok(ThermalBrittleFaceSustainedStep {
            motion: step.motion,
            energy: step.energy,
            native_energy: step.native_energy,
            solve: step.solve,
            rows: step.rows,
            normal_contact_work_j: step.normal_contact_work_j,
            absolute_contact_work_j: step.absolute_contact_work_j,
            contact_work_balance_error_j: step.contact_work_balance_error_j,
            absolute_contact_couple_n_m: step.absolute_contact_couple_n_m,
            numerical_tangent_evaluations: step.numerical_tangent_evaluations,
        })
    }
}

impl<'a> NativeFaceSustainedTrial for ThermalBrittleTetMidpointStep<'a> {
    type Motion = ThermalBrittleTetMotion<'a>;
    type Energy = ThermalBrittleTetMidpointEnergy;

    fn native_energy(&self) -> Self::Energy {
        self.energy
    }

    fn work(energy: &Self::Energy) -> SustainedWork {
        SustainedWork {
            dead_load_mean_work_j: energy.dead_load_mean_work_j,
            support_mean_work_j: energy.support_mean_work_j,
            mean_internal_work_j: energy.mean_internal_work_j,
            inertial_work_j: energy.inertia.inertial_work_j,
            storage_change_j: energy.stored_after_j - energy.stored_before_j,
            // The native fracture account already includes parameter loss.
            fracture_dissipation_j: energy.fracture_dissipation_j,
            kinetic_change_j: energy.inertia.kinetic_change_j,
            parameter_work_j: energy.parameter_work_j,
            eigenstrain_work_j: energy.eigenstrain_work_j,
        }
    }

    fn external_energy(
        native: Self::Energy,
        dead_work_j: f64,
        equilibrium_work_residual_j: f64,
        total_balance_error_j: f64,
    ) -> Self::Energy {
        ThermalBrittleTetMidpointEnergy {
            dead_load_mean_work_j: dead_work_j,
            equilibrium_work_residual_j,
            total_balance_error_j,
            ..native
        }
    }

    fn into_motion(self) -> Self::Motion {
        self.motion
    }
}
