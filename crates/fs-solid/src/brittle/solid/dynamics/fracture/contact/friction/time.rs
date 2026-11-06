//! Complete physical durations with directional impacts, static dwell and release.
//!
//! The original time driver owns affine forcing, support paths, localization,
//! response/event budgets and whole-duration work. This policy supplies the
//! original coupled frictional impulse and constrained nonlinear midpoint;
//! no second integrator or event loop is introduced.

use fs_contact::frictional_multicontact::elliptic::{
    EllipticContactResiduals, EllipticFrictionalContactResult,
};
use fs_contact::sustained::elliptic::EllipticForceCoordinates;
use fs_contact::sustained::solve::SolveReport;
use fs_exec::Cx;

use super::sustained::approach::{
    BrittleFaceSustainedFrictionApproachError, approach_released_faces_with_friction,
};
use super::sustained::{
    BrittleFaceSustainedFrictionControls, BrittleFaceSustainedFrictionError,
    BrittleFaceSustainedFrictionResult, solve_released_face_sustained_frictional_from,
};
use super::{
    BrittleFaceFriction, BrittleFaceFrictionControls, BrittleFaceFrictionError,
    BrittleFaceFrictionFrame, ReleasedFaceFrictionResponse, released_face_frictional_impact_with,
};
use crate::brittle::solid::dynamics::BrittleTetMotion;
use crate::brittle::solid::dynamics::fracture::contact::approach::{
    BrittleFaceApproachControls, BrittleFaceContactLocation,
};
use crate::brittle::solid::dynamics::fracture::contact::sustained::{
    BrittleFaceSustainedControls, BrittleFaceSustainedError, SustainedScope,
};
pub use crate::brittle::solid::dynamics::fracture::contact::time::{
    BrittleFaceMotionRecord, BrittleFaceTimeControls, BrittleFaceTimeEnergy,
};
use crate::brittle::solid::dynamics::fracture::contact::time::{
    BrittleFaceTimeError, FaceTimeImpactReceipt, FaceTimeImpactWork, FaceTimeSustainedApproach,
    FaceTimeSustainedReceipt, FaceTimeSustainedWork, ReleasedFaceTimeMotion,
    ReleasedFaceTimePolicy, run_released_face_time_with,
};
use crate::brittle::solid::dynamics::fracture::contact::{
    BrittleFaceImpactRow, ReleasedFaceImpactBody, SolidImpulseReport,
};
use crate::brittle::solid::dynamics::midpoint::BrittleTetMidpointEnergy;
use crate::brittle::solid::{BrittleTetError, checkpoint};
use crate::plastic3::IncrementSettings;

pub type BrittleFaceFrictionTimeControls = BrittleFaceTimeControls<BrittleFaceFrictionControls>;

/// Original directional law and native impulse at one accepted event time.
#[derive(Debug)]
pub struct BrittleFaceTimedFrictionImpact {
    pub location: BrittleFaceContactLocation,
    pub contact: EllipticFrictionalContactResult,
    pub impulse: SolidImpulseReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub frames: Vec<BrittleFaceFrictionFrame>,
    pub measured_residuals: EllipticContactResiduals,
    pub maximum_velocity_residual_mps: f64,
    pub relative_work_j: f64,
    pub contact_work_error_j: f64,
    pub relative_work_error_j: f64,
    pub energy_balance_error_j: f64,
    pub absolute_contact_couple_kg_m2_s: f64,
}

impl BrittleFaceTimedFrictionImpact {
    pub(crate) fn from_response(
        location: BrittleFaceContactLocation,
        response: ReleasedFaceFrictionResponse,
    ) -> (Vec<[f64; 3]>, Self) {
        (
            response.velocity_mps,
            Self {
                location,
                contact: response.contact,
                impulse: response.impulse,
                rows: response.rows,
                frames: response.frames,
                measured_residuals: response.measured_residuals,
                maximum_velocity_residual_mps: response.maximum_velocity_residual_mps,
                relative_work_j: response.relative_work_j,
                contact_work_error_j: response.contact_work_error_j,
                relative_work_error_j: response.relative_work_error_j,
                energy_balance_error_j: response.energy_balance_error_j,
                absolute_contact_couple_kg_m2_s: response.absolute_contact_couple_kg_m2_s,
            },
        )
    }
}

impl FaceTimeImpactReceipt for BrittleFaceTimedFrictionImpact {
    fn impact_work(&self) -> FaceTimeImpactWork {
        FaceTimeImpactWork {
            time_s: self.location.time_s,
            positive_normal_impulse: self.contact.impulses.iter().any(|p| p[0] > 0.0),
            support_work_j: self.impulse.support_work_j,
            kinetic_change_j: self.impulse.kinetic_change_j,
            dissipated_energy_j: self.contact.dissipated_energy,
            balance_error_j: self.energy_balance_error_j,
        }
    }
}

/// A compact native loaded-dwell or release receipt, without retaining a body copy.
#[derive(Debug)]
pub struct BrittleFaceTimedSustainedFriction<Energy = BrittleTetMidpointEnergy> {
    pub start_time_s: f64,
    pub end_time_s: f64,
    /// Native solve in its explicitly transformed force/path coordinates.
    pub solve: SolveReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub frames: Vec<BrittleFaceFrictionFrame>,
    pub forces_n: Vec<[f64; 3]>,
    pub paths_m: Vec<[f64; 3]>,
    pub coordinates: Vec<EllipticForceCoordinates>,
    pub native_energy: Energy,
    pub normal_contact_work_j: f64,
    pub tangential_contact_work_j: f64,
    pub absolute_contact_work_j: f64,
    pub absolute_normal_contact_work_j: f64,
    pub absolute_tangential_contact_work_j: f64,
    pub coordinate_work_error_j: f64,
    pub contact_work_balance_error_j: f64,
    pub absolute_contact_couple_n_m: f64,
}

impl<Energy: Copy> BrittleFaceTimedSustainedFriction<Energy> {
    pub(crate) fn from_step<M: ReleasedFaceTimeMotion<Energy = Energy>>(
        start_time_s: f64,
        step: BrittleFaceSustainedFrictionResult<M, Energy>,
    ) -> (M, Energy, Self) {
        let end_time_s = step.motion.time_view().time_s;
        (
            step.motion,
            step.energy,
            Self {
                start_time_s,
                end_time_s,
                solve: step.solve,
                rows: step.rows,
                frames: step.frames,
                forces_n: step.forces_n,
                paths_m: step.paths_m,
                coordinates: step.coordinates,
                native_energy: step.native_energy,
                normal_contact_work_j: step.normal_contact_work_j,
                tangential_contact_work_j: step.tangential_contact_work_j,
                absolute_contact_work_j: step.absolute_contact_work_j,
                absolute_normal_contact_work_j: step.absolute_normal_contact_work_j,
                absolute_tangential_contact_work_j: step.absolute_tangential_contact_work_j,
                coordinate_work_error_j: step.coordinate_work_error_j,
                contact_work_balance_error_j: step.contact_work_balance_error_j,
                absolute_contact_couple_n_m: step.absolute_contact_couple_n_m,
            },
        )
    }
}

impl<Energy> FaceTimeSustainedReceipt for BrittleFaceTimedSustainedFriction<Energy> {
    fn sustained_work(&self) -> FaceTimeSustainedWork {
        FaceTimeSustainedWork {
            normal_contact_work_j: self.normal_contact_work_j,
            tangential_contact_work_j: self.tangential_contact_work_j,
            absolute_normal_contact_work_j: self.absolute_normal_contact_work_j,
            absolute_tangential_contact_work_j: self.absolute_tangential_contact_work_j,
        }
    }
    fn physical_seed(&self) -> (&[BrittleFaceImpactRow], &[[f64; 3]]) {
        (&self.rows, &self.forces_n)
    }
}

/// Only a complete duration returns an accepted motion and its original receipts.
#[derive(Debug)]
pub struct BrittleFaceFrictionTimeResult<Motion, Energy> {
    pub motion: Motion,
    pub intervals: Vec<BrittleFaceMotionRecord<Energy>>,
    pub impacts: Vec<BrittleFaceTimedFrictionImpact>,
    pub sustained: Vec<BrittleFaceTimedSustainedFriction<Energy>>,
    pub trials: usize,
    pub numerical_tangent_evaluations: usize,
    pub energy: BrittleFaceTimeEnergy,
}

pub type BrittleFaceFrictionTimeStep<'a> =
    BrittleFaceFrictionTimeResult<BrittleTetMotion<'a>, BrittleTetMidpointEnergy>;

#[derive(Debug)]
pub enum BrittleFaceFrictionTimeError {
    Time(BrittleFaceTimeError),
    Impact(BrittleFaceFrictionError),
    SustainedApproach(BrittleFaceSustainedFrictionApproachError),
}
impl From<BrittleFaceTimeError> for BrittleFaceFrictionTimeError {
    fn from(error: BrittleFaceTimeError) -> Self {
        Self::Time(error)
    }
}
impl From<BrittleFaceFrictionError> for BrittleFaceFrictionTimeError {
    fn from(error: BrittleFaceFrictionError) -> Self {
        Self::Impact(error)
    }
}
impl From<BrittleFaceSustainedFrictionApproachError> for BrittleFaceFrictionTimeError {
    fn from(error: BrittleFaceSustainedFrictionApproachError) -> Self {
        Self::SustainedApproach(error)
    }
}
impl From<BrittleTetError> for BrittleFaceFrictionTimeError {
    fn from(error: BrittleTetError) -> Self {
        Self::Time(error.into())
    }
}
impl core::fmt::Display for BrittleFaceFrictionTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "released-face frictional duration refused: {self:?}")
    }
}
impl std::error::Error for BrittleFaceFrictionTimeError {}

/// Actual native body entry points, shared by fixed-law and source policies.
pub(crate) trait FrictionTimeMotion: ReleasedFaceTimeMotion {
    fn frictional_impact_body(&self) -> ReleasedFaceImpactBody<'_>;
    fn friction_velocity_successor(
        &self,
        velocity_mps: Vec<[f64; 3]>,
        kinetic_after_j: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleFaceFrictionError>;
    #[allow(clippy::too_many_arguments)]
    fn friction_force_response(
        &self,
        dt_s: f64,
        friction: &[BrittleFaceFriction],
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
        scope: SustainedScope,
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<
        BrittleFaceSustainedFrictionResult<Self, Self::Energy>,
        BrittleFaceSustainedFrictionError,
    >;
}

impl FrictionTimeMotion for BrittleTetMotion<'_> {
    fn frictional_impact_body(&self) -> ReleasedFaceImpactBody<'_> {
        ReleasedFaceImpactBody {
            reference_nodes_m: self.solid.geometry.nodes_m(),
            displacement_m: &self.solid.solution.displacement_m,
            released_faces: &self.solid.released_faces,
            mass: &self.mass,
            velocity_mps: &self.velocity_mps,
            kinetic_energy_j: self.kinetic_energy_j,
        }
    }
    fn friction_velocity_successor(
        &self,
        velocity_mps: Vec<[f64; 3]>,
        kinetic_after_j: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleFaceFrictionError> {
        BrittleTetMotion::friction_velocity_successor(self, velocity_mps, kinetic_after_j, cx)
    }
    fn friction_force_response(
        &self,
        dt_s: f64,
        friction: &[BrittleFaceFriction],
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
        scope: SustainedScope,
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<
        BrittleFaceSustainedFrictionResult<Self, Self::Energy>,
        BrittleFaceSustainedFrictionError,
    > {
        if friction.len() != self.solid.released_faces.len() {
            return Err(BrittleFaceSustainedError::Invalid(
                "complete retained friction inputs required",
            )
            .into());
        }
        let directions: Vec<_> = friction.iter().map(|f| f.tangent_direction).collect();
        solve_released_face_sustained_frictional_from(
            self.released_face_motion_view(),
            &directions,
            prescribed_m,
            mean_nodal_forces_n,
            seed,
            scope,
            controls,
            cx,
            |face, row, _, _, _| Ok(friction[face].vertices[row.vertex]),
            |total| {
                self.midpoint_candidate(
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

pub(crate) struct FrictionTimePolicy<'f> {
    pub friction: &'f [BrittleFaceFriction],
    pub sustained: BrittleFaceSustainedFrictionControls,
}

impl<M: FrictionTimeMotion> ReleasedFaceTimePolicy<M> for FrictionTimePolicy<'_> {
    type ImpactControls = BrittleFaceFrictionControls;
    type Impact = BrittleFaceTimedFrictionImpact;
    type Sustained = BrittleFaceTimedSustainedFriction<M::Energy>;
    type Error = BrittleFaceFrictionTimeError;

    fn validate(&mut self, initial: &M, cx: &Cx<'_>) -> Result<(), Self::Error> {
        checkpoint(cx)?;
        if self.friction.len() != initial.time_view().released_faces.len() {
            return Err(BrittleFaceTimeError::Invalid(
                "complete retained friction inputs required",
            )
            .into());
        }
        Ok(())
    }
    fn impact_gates(controls: Self::ImpactControls) -> (f64, f64) {
        (
            controls.activation_tolerance_m,
            controls.velocity_tolerance_mps,
        )
    }
    fn sustained_controls(&self) -> Option<BrittleFaceSustainedControls> {
        Some(self.sustained.contact)
    }

    fn impact(
        &mut self,
        motion: &M,
        location: BrittleFaceContactLocation,
        restitution: &[f64],
        held: &[usize],
        origin: [f64; 3],
        controls: Self::ImpactControls,
        cx: &Cx<'_>,
    ) -> Result<(M, Self::Impact), Self::Error> {
        let directions: Vec<_> = self.friction.iter().map(|f| f.tangent_direction).collect();
        let response = released_face_frictional_impact_with::<BrittleFaceFrictionError>(
            motion.frictional_impact_body(),
            restitution,
            &directions,
            held,
            origin,
            controls,
            cx,
            |face, row, _, _| Ok(self.friction[face].vertices[row.vertex]),
        )?;
        let (velocity, receipt) = BrittleFaceTimedFrictionImpact::from_response(location, response);
        let next =
            motion.friction_velocity_successor(velocity, receipt.impulse.kinetic_after_j, cx)?;
        Ok((next, receipt))
    }

    fn sustained(
        &mut self,
        motion: &M,
        duration: f64,
        prescribed: &[(usize, f64)],
        forces: &[(usize, [f64; 2])],
        seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
        settings: IncrementSettings,
        approach: BrittleFaceApproachControls,
        controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<FaceTimeSustainedApproach<M, Self::Sustained>, Self::Error> {
        let prefix = approach_released_faces_with_friction(
            motion,
            duration,
            self.friction,
            prescribed,
            forces,
            seed,
            settings,
            approach,
            BrittleFaceSustainedFrictionControls {
                contact: controls,
                ..self.sustained
            },
            cx,
        )?;
        Ok(FaceTimeSustainedApproach {
            advance: prefix.advance.map(|step| {
                BrittleFaceTimedSustainedFriction::from_step(motion.time_view().time_s, step)
            }),
            contact: prefix.contact,
            response_evaluations: prefix.response_evaluations,
            numerical_tangent_evaluations: prefix.numerical_tangent_evaluations,
        })
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Continue one complete affine load duration through frictional impacts,
    /// static load transfer, mixed arrivals and unloaded release.
    ///
    /// Fixed per-face ellipses are explicit model inputs for both event and
    /// dwell laws. During a loaded prefix every participating pair must retain
    /// matched geometry, zero physical slip and resting endpoint velocities.
    /// The native locator keeps those constraints while checking all unforced
    /// pairs for the next arrival. Impact applies one simultaneous coupled
    /// impulse, then the next interval uses the actual native successor.
    ///
    /// All discarded original responses spend the same whole-duration budget;
    /// impulses spend no physical time. Physical force seeds retain row identity
    /// and are cleared after impact. Static normal/tangent work, support work,
    /// fracture loss and combined impact loss have separate whole-duration
    /// accounts and absolute gates. Refusal publishes no partial motion.
    ///
    /// This remains fixed-reference paired contact. Loaded finite sliding,
    /// contact-pair discovery, global CCD, static/kinetic coefficient switching
    /// and impact heat allocation are not supplied. Resolve forcing reversals
    /// and choose the maximum interval with the native temporal resolution.
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
        controls: BrittleFaceFrictionTimeControls,
        sustained_controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceFrictionTimeStep<'a>, BrittleFaceFrictionTimeError> {
        let step = run_released_face_time_with(
            self,
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            restitution,
            origin_m,
            settings,
            controls,
            &mut FrictionTimePolicy {
                friction,
                sustained: sustained_controls,
            },
            cx,
        )?;
        Ok(BrittleFaceFrictionTimeResult {
            motion: step.motion,
            intervals: step.intervals,
            impacts: step.impacts,
            sustained: step.sustained,
            trials: step.trials,
            numerical_tangent_evaluations: step.numerical_tangent_evaluations,
            energy: step.energy,
        })
    }
}
