//! Complete frictional motion of the original thermally fractured body.
//!
//! Current thermal sources are held over this call. The mechanical and thermal
//! owners share the same time loop, constrained-prefix locator, directional
//! contact laws and work admission. Every thermal midpoint still returns its
//! complete original total/mechanical/eigenstrain accounts and reference matter.

use fs_exec::Cx;
use fs_matdb::MaterialStateId;

use super::{
    BrittleFaceFriction, BrittleFaceFrictionControls, BrittleFaceFrictionError,
    BrittleFaceFrictionSource,
};
pub use crate::brittle::solid::dynamics::fracture::contact::friction::source::time::{
    BrittleFaceFrictionSourcePhase, BrittleFaceSourcedFrictionTimeError,
    BrittleFaceTimedFrictionSources,
};
use crate::brittle::solid::dynamics::fracture::contact::friction::source::time::{
    BrittleFaceSourcedFrictionTimeResult, SourcedFrictionTimeMotion, run_sourced_friction_time,
};
use crate::brittle::solid::dynamics::fracture::contact::friction::sustained::source::{
    BrittleFaceSourcedSustainedFrictionError, BrittleFaceSourcedSustainedFrictionResult,
    solve_released_face_sustained_sourced_frictional_from,
};
use crate::brittle::solid::dynamics::fracture::contact::friction::sustained::{
    BrittleFaceSustainedFrictionControls, BrittleFaceSustainedFrictionError,
    BrittleFaceSustainedFrictionResult, solve_released_face_sustained_frictional_from,
};
pub use crate::brittle::solid::dynamics::fracture::contact::friction::time::{
    BrittleFaceFrictionTimeError, BrittleFaceTimedFrictionImpact, BrittleFaceTimedSustainedFriction,
};
use crate::brittle::solid::dynamics::fracture::contact::friction::time::{
    BrittleFaceFrictionTimeResult, FrictionTimeMotion, FrictionTimePolicy,
};
use crate::brittle::solid::dynamics::fracture::contact::sustained::SustainedScope;
pub use crate::brittle::solid::dynamics::fracture::contact::time::{
    BrittleFaceTimeControls, BrittleFaceTimeEnergy,
};
use crate::brittle::solid::dynamics::fracture::contact::time::{
    ReleasedFaceTimeMotion, run_released_face_time_with,
};
use crate::brittle::solid::dynamics::fracture::contact::{
    BrittleFaceImpactRow, ReleasedFaceImpactBody,
};
use crate::brittle::solid::expansion::dynamics::ThermalBrittleTetMotion;
use crate::brittle::solid::expansion::dynamics::midpoint::ThermalBrittleTetMidpointEnergy;
use crate::plastic3::IncrementSettings;

pub type ThermalBrittleFaceFrictionTimeStep<'a> =
    BrittleFaceFrictionTimeResult<ThermalBrittleTetMotion<'a>, ThermalBrittleTetMidpointEnergy>;
pub type ThermalBrittleFaceSourcedFrictionTimeStep<'a> = BrittleFaceSourcedFrictionTimeResult<
    ThermalBrittleTetMotion<'a>,
    ThermalBrittleTetMidpointEnergy,
>;

impl<'a> ThermalBrittleTetMotion<'a> {
    /// Advance the complete duration through directional impacts, sticking
    /// load transfer and release, with the native current thermal sources held.
    ///
    /// Every original affine support/load prefix, native step, force response,
    /// impulse, material history and signed/absolute work gate remains active.
    /// Discarded contact-localization responses do not consume physical steps.
    /// Only complete-duration acceptance publishes a successor. Reference mass
    /// retains its original density even if current thermal sources differ.
    ///
    /// Temperature changes require explicit resolved source transitions between
    /// calls. This method does not invent intermediate coefficient/eigenstrain
    /// paths, contact heat, finite-sliding geometry or new fragment pairs.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_time_frictional(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        restitution: &[f64],
        friction: &[BrittleFaceFriction],
        origin_m: [f64; 3],
        settings: IncrementSettings,
        controls: BrittleFaceTimeControls<BrittleFaceFrictionControls>,
        sustained: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceFrictionTimeStep<'a>, BrittleFaceFrictionTimeError> {
        let mut policy = FrictionTimePolicy {
            friction,
            sustained,
        };
        let result = run_released_face_time_with(
            self,
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            restitution,
            origin_m,
            settings,
            controls,
            &mut policy,
            cx,
        )?;
        Ok(BrittleFaceFrictionTimeResult {
            motion: result.motion,
            intervals: result.intervals,
            impacts: result.impacts,
            sustained: result.sustained,
            trials: result.trials,
            numerical_tangent_evaluations: result.numerical_tangent_evaluations,
            energy: result.energy,
        })
    }

    /// Run the same native thermal duration with explicit original interface
    /// sources. Impacts use each caller-selected role; sticking uses Static.
    /// Actual incoming/resting speeds drive source queries and the last accepted
    /// source candidate is retained by complete paired-row identity. External
    /// contact conditions are held; slip/travel/finite-step temperature remain
    /// unchanged. All source records remain private until total time acceptance.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_time_sourced_frictional(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        restitution: &[f64],
        sources: &[BrittleFaceFrictionSource<'_>],
        origin_m: [f64; 3],
        settings: IncrementSettings,
        controls: BrittleFaceTimeControls<BrittleFaceFrictionControls>,
        sustained: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceSourcedFrictionTimeStep<'a>, BrittleFaceSourcedFrictionTimeError>
    {
        run_sourced_friction_time(
            self,
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            restitution,
            sources,
            origin_m,
            settings,
            controls,
            sustained,
            cx,
        )
    }
}

impl<'a> FrictionTimeMotion for ThermalBrittleTetMotion<'a> {
    fn frictional_impact_body(&self) -> ReleasedFaceImpactBody<'_> {
        ThermalBrittleTetMotion::frictional_impact_body(self)
    }
    fn friction_velocity_successor(
        &self,
        velocity_mps: Vec<[f64; 3]>,
        kinetic_energy_j: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleFaceFrictionError> {
        ThermalBrittleTetMotion::friction_velocity_successor(
            self,
            velocity_mps,
            kinetic_energy_j,
            cx,
        )
    }
    fn friction_force_response(
        &self,
        dt_s: f64,
        friction: &[BrittleFaceFriction],
        prescribed_m: &[(usize, f64)],
        mean_forces_n: &[(usize, f64)],
        seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
        scope: SustainedScope,
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<
        BrittleFaceSustainedFrictionResult<Self, Self::Energy>,
        BrittleFaceSustainedFrictionError,
    > {
        let directions: Vec<_> = friction.iter().map(|law| law.tangent_direction).collect();
        solve_released_face_sustained_frictional_from(
            self.time_view(),
            &directions,
            prescribed_m,
            mean_forces_n,
            seed,
            scope,
            controls,
            cx,
            |face, row, _, _, _| Ok(friction[face].vertices[row.vertex]),
            |total| {
                self.midpoint_candidate(
                    self.solid.materials,
                    dt_s,
                    prescribed_m,
                    total,
                    settings,
                    controls.contact.work,
                    cx,
                )
            },
        )
    }
}

impl<'a> SourcedFrictionTimeMotion for ThermalBrittleTetMotion<'a> {
    fn friction_material(&self, element: usize) -> Option<&MaterialStateId> {
        self.solid
            .materials
            .get(element)
            .map(|source| source.band().material().material())
    }
    fn sourced_friction_force_response(
        &self,
        dt_s: f64,
        sources: &[BrittleFaceFrictionSource<'_>],
        prescribed_m: &[(usize, f64)],
        mean_forces_n: &[(usize, f64)],
        seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
        scope: SustainedScope,
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<
        BrittleFaceSourcedSustainedFrictionResult<Self, Self::Energy>,
        BrittleFaceSourcedSustainedFrictionError,
    > {
        solve_released_face_sustained_sourced_frictional_from(
            self.time_view(),
            |element| self.friction_material(element),
            sources,
            prescribed_m,
            mean_forces_n,
            seed,
            scope,
            controls,
            cx,
            |total| {
                self.midpoint_candidate(
                    self.solid.materials,
                    dt_s,
                    prescribed_m,
                    total,
                    settings,
                    controls.contact.work,
                    cx,
                )
            },
        )
    }
}
