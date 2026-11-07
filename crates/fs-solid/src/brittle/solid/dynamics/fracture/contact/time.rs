//! Complete physical intervals through native released-face contact modes.
//!
//! Motion uses the original work-conjugate midpoint and immutable material
//! history; contact uses the original simultaneous consistent-mass impulse.
//! Each actual prefix samples the original affine loads at its own times.
//! Only a complete duration publishes a successor. No penetrating candidate,
//! artificial gap correction, reset history or unresolved zero-time event is
//! allowed to escape into the retained trajectory.

use fs_contact::multicontact::MultiImpactResult;
use fs_contact::sustained::solve::SolveReport;
use fs_exec::Cx;

use super::approach::{
    BrittleFaceApproachControls, BrittleFaceApproachError, BrittleFaceContactLocation,
    ReleasedFaceMotionView,
};
use super::sustained::BrittleFaceSustainedControls;
use super::sustained::approach::BrittleFaceSustainedApproachError;
use super::{
    BrittleFaceImpactControls, BrittleFaceImpactError, BrittleFaceImpactRow, SolidImpulseReport,
};
use crate::brittle::solid::dynamics::midpoint::BrittleTetMidpointEnergy;
use crate::brittle::solid::dynamics::{BrittleTetDynamicError, BrittleTetMotion};
use crate::brittle::solid::path::{Sum, interpolate};
use crate::brittle::solid::{BrittleTetError, checkpoint, energy_gate};
use crate::plastic3::IncrementSettings;

/// Bounds on the complete duration, with the original motion/contact controls.
#[derive(Clone, Copy, Debug)]
pub struct BrittleFaceTimeControls<Impact = BrittleFaceImpactControls> {
    /// Largest free-flight or sustained-force interval proposed to native motion.
    /// This is a caller-selected temporal resolution, not an error estimate.
    pub maximum_interval_s: f64,
    pub maximum_intervals: usize,
    /// Includes instantaneous events at the original and final physical times.
    pub maximum_impacts: usize,
    /// All native motion solves, including discarded localization trials and
    /// every sustained-force response when that path is enabled.
    pub maximum_trials: usize,
    pub approach: BrittleFaceApproachControls,
    pub impact: Impact,
    /// Whole-duration signed balance AND sum of absolute local balance errors.
    pub total_energy_tolerance_j: f64,
}

/// One accepted physical midpoint interval, without a copy of its full fields.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceMotionRecord<Energy = BrittleTetMidpointEnergy> {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub energy: Energy,
}

/// Original coupled-contact and body-impulse receipts at one fixed configuration.
#[derive(Debug)]
pub struct BrittleFaceTimedImpact {
    pub location: BrittleFaceContactLocation,
    pub contact: MultiImpactResult,
    pub impulse: SolidImpulseReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub maximum_velocity_residual_mps: f64,
    pub relative_work_error_j: f64,
    pub energy_balance_error_j: f64,
    pub absolute_contact_couple_kg_m2_s: f64,
}

/// A loaded dwell or unloaded release on the original nonlinear body response.
#[derive(Debug)]
pub struct BrittleFaceTimedSustained<Energy = BrittleTetMidpointEnergy> {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub solve: SolveReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    /// Original total-load receipt; the matching interval uses external-only work.
    pub native_energy: Energy,
    pub normal_contact_work_j: f64,
    pub absolute_contact_work_j: f64,
    pub contact_work_balance_error_j: f64,
    pub absolute_contact_couple_n_m: f64,
}

/// Whole-duration physical work. Impact loss remains separate from fracture.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleFaceTimeEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Compensated sum of native polarized kinetic increments, including events.
    pub kinetic_change_j: f64,
    /// Mean material stress work on TOTAL displacement strain.
    pub material_work_j: f64,
    /// Signed work from source coefficients; zero for fixed nonthermal laws.
    pub parameter_work_j: f64,
    /// Minus mean stress contracted with the thermal free-strain change.
    pub eigenstrain_work_j: f64,
    pub fracture_dissipation_j: f64,
    /// Surviving-bond physical loss, kept separate from fracture and contact.
    pub viscous_dissipation_j: f64,
    /// Mechanical loss from the contact owner; no thermal allocation is inferred.
    pub contact_dissipation_j: f64,
    /// Bounded numerical work at sustained normal constraints, never a loss
    /// or an additional external-work port in the physical balance.
    pub normal_contact_work_j: f64,
    pub absolute_normal_contact_work_j: f64,
    /// Bounded static tangential work, separate from physical impact loss.
    pub tangential_contact_work_j: f64,
    pub absolute_tangential_contact_work_j: f64,
    pub dead_load_mean_work_j: f64,
    pub support_mean_work_j: f64,
    /// Signed reaction work at held, possibly moving velocity components.
    pub support_impulse_work_j: f64,
    pub material_balance_error_j: f64,
    pub mean_work_error_j: f64,
    pub absolute_mean_work_error_j: f64,
    /// External plus parameter/eigenstrain work minus elastic/kinetic changes
    /// and the separate fracture/viscous/contact losses.
    pub total_balance_error_j: f64,
    /// Prevents opposite local numerical defects from hiding in the total.
    pub absolute_local_balance_error_j: f64,
}

impl BrittleFaceTimeEnergy {
    #[must_use]
    pub fn external_work_j(&self) -> f64 {
        (self.dead_load_mean_work_j + self.support_mean_work_j) + self.support_impulse_work_j
    }
}

#[derive(Debug)]
pub struct BrittleFaceTimeStep<'a> {
    pub motion: BrittleTetMotion<'a>,
    pub intervals: Vec<BrittleFaceMotionRecord>,
    pub impacts: Vec<BrittleFaceTimedImpact>,
    /// Also represented in `intervals` with external-only physical work.
    pub sustained: Vec<BrittleFaceTimedSustained>,
    pub trials: usize,
    /// Completed native trial telemetry; a refused native solve has no receipt.
    pub numerical_tangent_evaluations: usize,
    pub energy: BrittleFaceTimeEnergy,
}

#[derive(Debug)]
pub enum BrittleFaceTimeError {
    Invalid(&'static str),
    Allocation,
    Motion(BrittleTetDynamicError),
    Approach(BrittleFaceApproachError),
    Impact(BrittleFaceImpactError),
    SustainedApproach(BrittleFaceSustainedApproachError),
    WorkLimit {
        intervals: usize,
        impacts: usize,
        trials: usize,
    },
    /// No two events may occur within one localization time tolerance.
    /// An instantaneous impulse cannot replace a sustained contact-force solve.
    UnresolvedContact {
        time_s: f64,
    },
}

impl core::fmt::Display for BrittleFaceTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "released-face physical interval refused: {self:?}")
    }
}
impl std::error::Error for BrittleFaceTimeError {}
impl From<BrittleTetError> for BrittleFaceTimeError {
    fn from(error: BrittleTetError) -> Self {
        Self::Motion(error.into())
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Advance the complete affine load interval, automatically resolving impacts.
    ///
    /// Supports follow an affine displacement path; their initial velocities
    /// must match its slope within the explicit velocity tolerance. The same
    /// support DOFs retain their velocities during each instantaneous impact.
    /// Applied-force pairs give values at this duration's original endpoints.
    /// Every shortened prefix and remainder recomputes its own mean forces.
    ///
    /// Each accepted midpoint consumes one native motion step. Impacts consume
    /// an event allowance and no physical timestep. Discarded localization
    /// trials consume only the solve budget. An exhausted budget, cancellation,
    /// unsupported geometry, physical-work refusal or unresolved contact returns
    /// no partially advanced motion; the original is available for retry.
    ///
    /// The native prefix locator admits nonpenetrating affine displacement
    /// traces on fixed paired faces. Its bracket is a local numerical event
    /// resolution, not a global continuous-collision certificate. Choose the
    /// maximum interval to resolve forcing and motion; split load reversals at
    /// their knots. Finite sliding, new fragment pairs, friction, persistent
    /// contact forces and thermal allocation of impact loss are outside this
    /// entry. Repeated events separated by no more than the localization time
    /// tolerance refuse instead of generating an unresolved impulse sequence.
    /// No temporal convergence order is asserted across impact events.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_time(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        restitution: &[f64],
        origin_m: [f64; 3],
        settings: IncrementSettings,
        controls: BrittleFaceTimeControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceTimeStep<'a>, BrittleFaceTimeError> {
        self.released_face_time(
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            restitution,
            origin_m,
            settings,
            controls,
            None,
            cx,
        )
    }

    /// Continue free motion, impacts, loaded contact and unloaded release.
    ///
    /// The same interval driver resolves approaching impacts first. At active,
    /// normally resting pairs it invokes the original sustained-force response
    /// instead of repeated zero-time impulses. Initially resting pairs retain
    /// their coupled force response while all released pairs are checked for
    /// new arrivals. A newly approaching pair shortens the constrained interval
    /// before the simultaneous impact includes every active pair. Every loaded
    /// pair retains its native dwell, complementarity, geometry and work gates.
    /// Released pairs return to the original free-flight path.
    ///
    /// Previous normal forces seed only the next numerical solve; they are
    /// cleared after impact and never become physical material history. The
    /// whole native-response budget includes every force trial and final replay.
    /// Original affine force means, prescribed support paths, immutable retry
    /// and complete-duration energy admission remain shared with impact-only
    /// time. Sustained normal work has an independent whole-duration absolute
    /// gate and is never attributed to contact loss, fracture or heat.
    ///
    /// The force and approach activation/normal-velocity tolerances must match.
    /// Fixed matching traces, small-strain source domains and caller-selected
    /// temporal resolution remain in force; this supplies no global CCD,
    /// finite sliding, friction or mixed-mode trajectory accuracy certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_time_with_sustained(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        restitution: &[f64],
        origin_m: [f64; 3],
        settings: IncrementSettings,
        controls: BrittleFaceTimeControls,
        sustained_controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceTimeStep<'a>, BrittleFaceTimeError> {
        self.released_face_time(
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            restitution,
            origin_m,
            settings,
            controls,
            Some(sustained_controls),
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn released_face_time(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        restitution: &[f64],
        origin_m: [f64; 3],
        settings: IncrementSettings,
        controls: BrittleFaceTimeControls,
        sustained_controls: Option<BrittleFaceSustainedControls>,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceTimeStep<'a>, BrittleFaceTimeError> {
        let step = run_released_face_time(
            self,
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            restitution,
            origin_m,
            settings,
            controls,
            sustained_controls.map(|controls| SustainedTimePath {
                controls,
                advance: advance_native_sustained,
            }),
            cx,
        )?;
        Ok(BrittleFaceTimeStep {
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

impl<'a> ReleasedFaceTimeMotion for BrittleTetMotion<'a> {
    type Energy = BrittleTetMidpointEnergy;

    fn time_view(&self) -> ReleasedFaceMotionView<'_> {
        self.released_face_motion_view()
    }

    fn time_inventory(&self) -> FaceTimeInventory {
        FaceTimeInventory {
            stored_energy_j: self.solid.stored_energy_j,
            fracture_dissipation_j: self.solid.fracture_dissipation_j,
            kinetic_energy_j: self.kinetic_energy_j,
        }
    }

    fn time_work(energy: &Self::Energy) -> FaceTimeWork {
        FaceTimeWork {
            material_work_j: energy.material_work_j,
            parameter_work_j: 0.0,
            eigenstrain_work_j: 0.0,
            fracture_dissipation_j: energy.fracture_dissipation_j,
            viscous_dissipation_j: energy.viscous_dissipation_j,
            dead_load_mean_work_j: energy.dead_load_mean_work_j,
            support_mean_work_j: energy.support_mean_work_j,
            kinetic_change_j: energy.inertia.kinetic_change_j,
            total_balance_error_j: energy.total_balance_error_j,
            mean_work_error_j: energy.mean_work_error_j,
            mean_internal_work_j: energy.mean_internal_work_j,
        }
    }

    fn time_approach(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleFaceApproachControls,
        cx: &Cx<'_>,
    ) -> Result<FaceTimeApproach<Self>, BrittleFaceApproachError> {
        let prefix = self.advance_toward_released_face_contact(
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            settings,
            controls,
            cx,
        )?;
        Ok(FaceTimeApproach {
            advance: prefix.advance.map(|step| (step.motion, step.energy)),
            contact: prefix.contact,
            trials: prefix.trials,
            numerical_tangent_evaluations: prefix.numerical_tangent_evaluations,
        })
    }

    fn time_impact(
        &self,
        location: BrittleFaceContactLocation,
        restitution: &[f64],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: BrittleFaceImpactControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, BrittleFaceTimedImpact), BrittleFaceImpactError> {
        let hit =
            self.impact_released_faces(restitution, held_velocity_dofs, origin_m, controls, cx)?;
        Ok((
            hit.motion,
            BrittleFaceTimedImpact {
                location,
                contact: hit.contact,
                impulse: hit.impulse,
                rows: hit.rows,
                maximum_velocity_residual_mps: hit.maximum_velocity_residual_mps,
                relative_work_error_j: hit.relative_work_error_j,
                energy_balance_error_j: hit.energy_balance_error_j,
                absolute_contact_couple_kg_m2_s: hit.absolute_contact_couple_kg_m2_s,
            },
        ))
    }
}

#[allow(clippy::too_many_arguments)]
fn advance_native_sustained<'a>(
    motion: &BrittleTetMotion<'a>,
    duration_s: f64,
    prescribed_end_m: &[(usize, f64)],
    force_intervals_n: &[(usize, [f64; 2])],
    seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
    settings: IncrementSettings,
    approach: BrittleFaceApproachControls,
    controls: BrittleFaceSustainedControls,
    cx: &Cx<'_>,
) -> Result<FaceTimeSustainedApproach<BrittleTetMotion<'a>>, BrittleFaceSustainedApproachError> {
    let prefix = motion.advance_toward_released_face_contact_with_sustained_from(
        duration_s,
        prescribed_end_m,
        force_intervals_n,
        seed,
        settings,
        approach,
        controls,
        cx,
    )?;
    let advance = prefix.advance.map(|step| {
        let record = BrittleFaceTimedSustained {
            start_time_s: motion.time_s,
            end_time_s: step.motion.time_s,
            solve: step.solve,
            rows: step.rows,
            native_energy: step.native_energy,
            normal_contact_work_j: step.normal_contact_work_j,
            absolute_contact_work_j: step.absolute_contact_work_j,
            contact_work_balance_error_j: step.contact_work_balance_error_j,
            absolute_contact_couple_n_m: step.absolute_contact_couple_n_m,
        };
        (step.motion, step.energy, record)
    });
    Ok(FaceTimeSustainedApproach {
        advance,
        contact: prefix.contact,
        response_evaluations: prefix.response_evaluations,
        numerical_tangent_evaluations: prefix.numerical_tangent_evaluations,
    })
}

/// Borrowed native motion and exact retained inventory; no constitutive cast.
pub(crate) struct FaceTimeInventory {
    pub stored_energy_j: f64,
    pub fracture_dissipation_j: f64,
    pub kinetic_energy_j: f64,
}

/// Physical work projection only; complete native receipts remain in records.
pub(crate) struct FaceTimeWork {
    pub material_work_j: f64,
    pub parameter_work_j: f64,
    pub eigenstrain_work_j: f64,
    pub fracture_dissipation_j: f64,
    pub viscous_dissipation_j: f64,
    pub dead_load_mean_work_j: f64,
    pub support_mean_work_j: f64,
    pub kinetic_change_j: f64,
    pub total_balance_error_j: f64,
    pub mean_work_error_j: f64,
    pub mean_internal_work_j: f64,
}

pub(crate) struct FaceTimeApproach<M: ReleasedFaceTimeMotion> {
    pub advance: Option<(M, M::Energy)>,
    pub contact: Option<BrittleFaceContactLocation>,
    pub trials: usize,
    pub numerical_tangent_evaluations: usize,
}

/// The two actual material owners provide native responses to one time loop.
/// No cloned substitute body, material interpolation or alternate solver enters.
pub(crate) trait ReleasedFaceTimeMotion: Sized {
    type Energy: Copy;

    fn time_view(&self) -> ReleasedFaceMotionView<'_>;
    fn time_inventory(&self) -> FaceTimeInventory;
    fn time_work(energy: &Self::Energy) -> FaceTimeWork;

    fn time_approach(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleFaceApproachControls,
        cx: &Cx<'_>,
    ) -> Result<FaceTimeApproach<Self>, BrittleFaceApproachError>;

    #[allow(clippy::too_many_arguments)]
    fn time_impact(
        &self,
        location: BrittleFaceContactLocation,
        restitution: &[f64],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: BrittleFaceImpactControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, BrittleFaceTimedImpact), BrittleFaceImpactError>;
}

/// Optional native constrained-prefix response, supplied only by an admitted owner.
/// The thermal impact-only entry supplies none; it has no stub force response.
pub(crate) struct SustainedTimePath<M: ReleasedFaceTimeMotion> {
    pub controls: BrittleFaceSustainedControls,
    pub advance: fn(
        &M,
        f64,
        &[(usize, f64)],
        &[(usize, [f64; 2])],
        Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
        IncrementSettings,
        BrittleFaceApproachControls,
        BrittleFaceSustainedControls,
        &Cx<'_>,
    ) -> Result<FaceTimeSustainedApproach<M>, BrittleFaceSustainedApproachError>,
}

pub(crate) struct FaceTimeSustainedApproach<
    M: ReleasedFaceTimeMotion,
    Receipt = BrittleFaceTimedSustained<<M as ReleasedFaceTimeMotion>::Energy>,
> {
    pub advance: Option<(M, M::Energy, Receipt)>,
    pub contact: Option<BrittleFaceContactLocation>,
    pub response_evaluations: usize,
    pub numerical_tangent_evaluations: usize,
}

pub(crate) struct NativeFaceTimeStep<
    M: ReleasedFaceTimeMotion,
    Impact = BrittleFaceTimedImpact,
    Sustained = BrittleFaceTimedSustained<<M as ReleasedFaceTimeMotion>::Energy>,
> {
    pub motion: M,
    pub intervals: Vec<BrittleFaceMotionRecord<M::Energy>>,
    pub impacts: Vec<Impact>,
    pub sustained: Vec<Sustained>,
    pub trials: usize,
    pub numerical_tangent_evaluations: usize,
    pub energy: BrittleFaceTimeEnergy,
}

/// Physical event accounts only; complete native law receipts remain typed.
pub(crate) struct FaceTimeImpactWork {
    pub time_s: f64,
    pub positive_normal_impulse: bool,
    pub support_work_j: f64,
    pub kinetic_change_j: f64,
    pub dissipated_energy_j: f64,
    pub balance_error_j: f64,
}

pub(crate) trait FaceTimeImpactReceipt {
    fn impact_work(&self) -> FaceTimeImpactWork;
}

pub(crate) struct FaceTimeSustainedWork {
    pub normal_contact_work_j: f64,
    pub tangential_contact_work_j: f64,
    pub absolute_normal_contact_work_j: f64,
    pub absolute_tangential_contact_work_j: f64,
}

pub(crate) trait FaceTimeSustainedReceipt {
    fn sustained_work(&self) -> FaceTimeSustainedWork;
    /// Physical forces in the retained contact frames, never mapped solver forces.
    fn physical_seed(&self) -> (&[BrittleFaceImpactRow], &[[f64; 3]]);
}

impl FaceTimeImpactReceipt for BrittleFaceTimedImpact {
    fn impact_work(&self) -> FaceTimeImpactWork {
        FaceTimeImpactWork {
            time_s: self.location.time_s,
            positive_normal_impulse: self.contact.normal_impulses.iter().any(|&p| p > 0.0),
            support_work_j: self.impulse.support_work_j,
            kinetic_change_j: self.impulse.kinetic_change_j,
            dissipated_energy_j: self.contact.dissipated_energy,
            balance_error_j: self.energy_balance_error_j,
        }
    }
}

impl<Energy> FaceTimeSustainedReceipt for BrittleFaceTimedSustained<Energy> {
    fn sustained_work(&self) -> FaceTimeSustainedWork {
        FaceTimeSustainedWork {
            normal_contact_work_j: self.normal_contact_work_j,
            tangential_contact_work_j: 0.0,
            absolute_normal_contact_work_j: self.absolute_contact_work_j,
            absolute_tangential_contact_work_j: 0.0,
        }
    }

    fn physical_seed(&self) -> (&[BrittleFaceImpactRow], &[[f64; 3]]) {
        (&self.rows, &self.solve.forces_n)
    }
}

/// A contact-law family supplies original impulses and constrained prefixes to
/// one physical-time driver. Source-bearing policies may stage candidates in
/// these callbacks; only complete successful time results are published.
pub(crate) trait ReleasedFaceTimePolicy<M: ReleasedFaceTimeMotion> {
    type ImpactControls: Copy;
    type Impact: FaceTimeImpactReceipt;
    type Sustained: FaceTimeSustainedReceipt;
    type Error: From<BrittleFaceTimeError>;

    fn validate(&mut self, initial: &M, cx: &Cx<'_>) -> Result<(), Self::Error>;
    fn impact_gates(controls: Self::ImpactControls) -> (f64, f64);
    fn sustained_controls(&self) -> Option<BrittleFaceSustainedControls>;

    #[allow(clippy::too_many_arguments)]
    fn impact(
        &mut self,
        motion: &M,
        location: BrittleFaceContactLocation,
        restitution: &[f64],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: Self::ImpactControls,
        cx: &Cx<'_>,
    ) -> Result<(M, Self::Impact), Self::Error>;

    #[allow(clippy::too_many_arguments)]
    fn sustained(
        &mut self,
        motion: &M,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
        settings: IncrementSettings,
        approach: BrittleFaceApproachControls,
        controls: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<FaceTimeSustainedApproach<M, Self::Sustained>, Self::Error>;
}

struct NormalTimePolicy<M: ReleasedFaceTimeMotion> {
    sustained: Option<SustainedTimePath<M>>,
}

impl<M: ReleasedFaceTimeMotion> ReleasedFaceTimePolicy<M> for NormalTimePolicy<M> {
    type ImpactControls = BrittleFaceImpactControls;
    type Impact = BrittleFaceTimedImpact;
    type Sustained = BrittleFaceTimedSustained<M::Energy>;
    type Error = BrittleFaceTimeError;

    fn validate(&mut self, _: &M, _: &Cx<'_>) -> Result<(), Self::Error> {
        Ok(())
    }

    fn impact_gates(controls: Self::ImpactControls) -> (f64, f64) {
        (
            controls.activation_tolerance_m,
            controls.velocity_tolerance_mps,
        )
    }

    fn sustained_controls(&self) -> Option<BrittleFaceSustainedControls> {
        self.sustained.as_ref().map(|path| path.controls)
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
        motion
            .time_impact(location, restitution, held, origin, controls, cx)
            .map_err(BrittleFaceTimeError::Impact)
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
        let path = self
            .sustained
            .as_ref()
            .ok_or(BrittleFaceTimeError::Invalid(
                "missing configured native sustained path",
            ))?;
        (path.advance)(
            motion, duration, prescribed, forces, seed, settings, approach, controls, cx,
        )
        .map_err(BrittleFaceTimeError::SustainedApproach)
    }
}

enum TimeRunError<E> {
    Time(BrittleFaceTimeError),
    Policy(E),
}
impl<E> From<BrittleFaceTimeError> for TimeRunError<E> {
    fn from(error: BrittleFaceTimeError) -> Self {
        Self::Time(error)
    }
}
impl<E> From<BrittleTetError> for TimeRunError<E> {
    fn from(error: BrittleTetError) -> Self {
        Self::Time(error.into())
    }
}

/// Complete native interval driver shared by thermal and nonthermal bodies.
/// Every accepted record retains its original owner's energy type.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(crate) fn run_released_face_time<M: ReleasedFaceTimeMotion>(
    initial: &M,
    duration_s: f64,
    prescribed_end_m: &[(usize, f64)],
    force_intervals_n: &[(usize, [f64; 2])],
    restitution: &[f64],
    origin_m: [f64; 3],
    settings: IncrementSettings,
    controls: BrittleFaceTimeControls,
    sustained_path: Option<SustainedTimePath<M>>,
    cx: &Cx<'_>,
) -> Result<NativeFaceTimeStep<M>, BrittleFaceTimeError> {
    run_released_face_time_with(
        initial,
        duration_s,
        prescribed_end_m,
        force_intervals_n,
        restitution,
        origin_m,
        settings,
        controls,
        &mut NormalTimePolicy {
            sustained: sustained_path,
        },
        cx,
    )
}

/// The original complete-duration loop, parameterized only by its contact law.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_released_face_time_with<
    M: ReleasedFaceTimeMotion,
    P: ReleasedFaceTimePolicy<M>,
>(
    initial: &M,
    duration_s: f64,
    prescribed_end_m: &[(usize, f64)],
    force_intervals_n: &[(usize, [f64; 2])],
    restitution: &[f64],
    origin_m: [f64; 3],
    settings: IncrementSettings,
    controls: BrittleFaceTimeControls<P::ImpactControls>,
    policy: &mut P,
    cx: &Cx<'_>,
) -> Result<NativeFaceTimeStep<M, P::Impact, P::Sustained>, P::Error> {
    run_time(
        initial,
        duration_s,
        prescribed_end_m,
        force_intervals_n,
        restitution,
        origin_m,
        settings,
        controls,
        policy,
        cx,
    )
    .map_err(|error| match error {
        TimeRunError::Time(error) => error.into(),
        TimeRunError::Policy(error) => error,
    })
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn run_time<M: ReleasedFaceTimeMotion, P: ReleasedFaceTimePolicy<M>>(
    initial: &M,
    duration_s: f64,
    prescribed_end_m: &[(usize, f64)],
    force_intervals_n: &[(usize, [f64; 2])],
    restitution: &[f64],
    origin_m: [f64; 3],
    settings: IncrementSettings,
    controls: BrittleFaceTimeControls<P::ImpactControls>,
    policy: &mut P,
    cx: &Cx<'_>,
) -> Result<NativeFaceTimeStep<M, P::Impact, P::Sustained>, TimeRunError<P::Error>> {
    use BrittleFaceTimeError as E;
    checkpoint(cx)?;
    policy.validate(initial, cx).map_err(TimeRunError::Policy)?;
    let initial_view = initial.time_view();
    let initial_inventory = initial.time_inventory();
    let start_time_s = initial_view.time_s;
    let target_time = start_time_s + duration_s;
    let sustained_controls = policy.sustained_controls();
    let (activation_tolerance_m, velocity_tolerance_mps) = P::impact_gates(controls.impact);
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || !target_time.is_finite()
        || target_time <= start_time_s
        || !controls.maximum_interval_s.is_finite()
        || controls.maximum_interval_s <= 0.0
        || controls.maximum_intervals == 0
        || controls.maximum_trials == 0
        || !controls.total_energy_tolerance_j.is_finite()
        || controls.total_energy_tolerance_j < 0.0
        || controls.approach.activation_tolerance_m != activation_tolerance_m
        || controls.approach.velocity_tolerance_mps < velocity_tolerance_mps
        || restitution.len() != initial_view.released_faces.len()
        || restitution
            .iter()
            .any(|e| !e.is_finite() || !(0.0..=1.0).contains(e))
        || origin_m.iter().any(|x| !x.is_finite())
        || sustained_controls.is_some_and(|c| {
            c.activation_tolerance_m != controls.approach.activation_tolerance_m
                || c.velocity_tolerance_mps != controls.approach.velocity_tolerance_mps
        })
    {
        return Err(E::Invalid(
            "positive bounded duration, compatible event gates and complete restitution required",
        )
        .into());
    }
    if sustained_controls.is_some() {
        super::approach::validate(
            duration_s,
            start_time_s,
            target_time,
            settings,
            controls.approach,
        )
        .map_err(E::Approach)?;
    }
    let duration = target_time - start_time_s;
    let initial_u = initial_view.displacement_m.as_flattened();
    let initial_v = initial_view.velocity_mps.as_flattened();
    if prescribed_end_m.len() > initial_u.len() || force_intervals_n.len() > initial_u.len() {
        return Err(E::Invalid("load entries exceed the original DOF count").into());
    }
    let mut supports = buffer(prescribed_end_m.len())?;
    let mut held = buffer(prescribed_end_m.len())?;
    for &(dof, end) in prescribed_end_m {
        checkpoint(cx)?;
        let Some(&start) = initial_u.get(dof) else {
            return Err(E::Invalid("prescribed DOF outside the retained mesh").into());
        };
        let speed = (end - start) / duration;
        if !end.is_finite()
            || !speed.is_finite()
            || (speed - initial_v[dof]).abs() > controls.approach.velocity_tolerance_mps
        {
            return Err(E::Invalid("initial support velocity must match its affine path").into());
        }
        supports.push((dof, start));
        held.push(dof);
    }
    let mut forces = buffer(force_intervals_n.len())?;
    for &(dof, values) in force_intervals_n {
        checkpoint(cx)?;
        if dof >= initial_u.len() || values.iter().any(|x| !x.is_finite()) {
            return Err(E::Invalid("finite applied forces on existing DOFs required").into());
        }
        forces.push((dof, values));
    }
    // Reuse the last accepted force receipt by paired-row identity. Changing
    // contact subsets cannot shift a numerical seed onto an unrelated row.
    let mut seed_record: Option<usize> = None;
    let mut current: Option<M> = None;
    let mut intervals = buffer(controls.maximum_intervals.min(32))?;
    let mut impacts: Vec<P::Impact> = buffer(controls.maximum_impacts.min(8))?;
    let mut sustained: Vec<P::Sustained> = buffer(if sustained_controls.is_some() {
        controls.maximum_intervals.min(8)
    } else {
        0
    })?;
    let mut sums = [Sum::default(); 13];
    let mut normal_work = Sum::default();
    let mut absolute_normal_work = Sum::default();
    let mut tangential_work = Sum::default();
    let mut absolute_tangential_work = Sum::default();
    let mut trials = 0usize;
    let mut tangent_evaluations = 0usize;
    loop {
        checkpoint(cx)?;
        let accepted = current.as_ref().unwrap_or(initial);
        let accepted_view = accepted.time_view();
        if accepted_view.time_s == target_time {
            break;
        }
        if trials >= controls.maximum_trials || intervals.len() == controls.maximum_intervals {
            return Err(E::WorkLimit {
                intervals: intervals.len(),
                impacts: impacts.len(),
                trials,
            }
            .into());
        }
        let end = (accepted_view.time_s + controls.maximum_interval_s).min(target_time);
        if !end.is_finite() || end <= accepted_view.time_s {
            return Err(E::Invalid("increasing native physical time must be representable").into());
        }
        let start_fraction = (accepted_view.time_s - start_time_s) / duration;
        let end_fraction = if end == target_time {
            1.0
        } else {
            (end - start_time_s) / duration
        };
        for (i, (entry, &(dof, target))) in supports.iter_mut().zip(prescribed_end_m).enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            entry.1 = interpolate(initial_u[dof], target, end_fraction);
        }
        for (i, (entry, &(_, [start, target]))) in
            forces.iter_mut().zip(force_intervals_n).enumerate()
        {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            entry.1 = [
                interpolate(start, target, start_fraction),
                interpolate(start, target, end_fraction),
            ];
        }
        let (advance, contact, calls, tangents) = if let Some(mut c) = sustained_controls
            && resting_contact(&accepted_view, controls.approach, cx)?
        {
            c.solve.maximum_response_evaluations = c
                .solve
                .maximum_response_evaluations
                .min(controls.maximum_trials - trials);
            let mut approach_controls = controls.approach;
            approach_controls.maximum_trials = approach_controls
                .maximum_trials
                .min(controls.maximum_trials - trials);
            let seed = seed_record.map(|index| sustained[index].physical_seed());
            let prefix = policy
                .sustained(
                    accepted,
                    end - accepted_view.time_s,
                    &supports,
                    &forces,
                    seed,
                    settings,
                    approach_controls,
                    c,
                    cx,
                )
                .map_err(TimeRunError::Policy)?;
            let advance = if let Some((motion, energy, record)) = prefix.advance {
                let work = record.sustained_work();
                normal_work.add(work.normal_contact_work_j);
                absolute_normal_work.add(work.absolute_normal_contact_work_j);
                tangential_work.add(work.tangential_contact_work_j);
                absolute_tangential_work.add(work.absolute_tangential_contact_work_j);
                if !normal_work.value().is_finite()
                    || !absolute_normal_work.value().is_finite()
                    || !tangential_work.value().is_finite()
                    || !absolute_tangential_work.value().is_finite()
                {
                    return Err(E::Invalid(
                        "complete normal and tangential contact work must be representable",
                    )
                    .into());
                }
                push(&mut sustained, record)?;
                seed_record = Some(sustained.len() - 1);
                Some((motion, energy))
            } else {
                None
            };
            (
                advance,
                prefix.contact,
                prefix.response_evaluations,
                prefix.numerical_tangent_evaluations,
            )
        } else {
            let mut approach_controls = controls.approach;
            approach_controls.maximum_trials = approach_controls
                .maximum_trials
                .min(controls.maximum_trials - trials);
            let prefix = accepted
                .time_approach(
                    end - accepted_view.time_s,
                    &supports,
                    &forces,
                    settings,
                    approach_controls,
                    cx,
                )
                .map_err(E::Approach)?;
            (
                prefix.advance,
                prefix.contact,
                prefix.trials,
                prefix.numerical_tangent_evaluations,
            )
        };
        trials = trials
            .checked_add(calls)
            .ok_or(E::Invalid("native trial count overflow"))?;
        tangent_evaluations = tangent_evaluations
            .checked_add(tangents)
            .ok_or(E::Invalid("native tangent count overflow"))?;
        if advance.is_none() && contact.is_none() {
            return Err(
                E::Invalid("prefix must advance physical time or resolve a contact").into(),
            );
        }
        if let Some((next, energy)) = advance {
            if intervals.len() == controls.maximum_intervals {
                return Err(E::WorkLimit {
                    intervals: intervals.len(),
                    impacts: impacts.len(),
                    trials,
                }
                .into());
            }
            let next_time_s = next.time_view().time_s;
            if next_time_s <= accepted_view.time_s || next_time_s > end {
                return Err(
                    E::Invalid("native prefix escaped its requested physical interval").into(),
                );
            }
            let work = M::time_work(&energy);
            accumulate(
                &mut sums,
                [
                    work.material_work_j,
                    work.fracture_dissipation_j,
                    work.dead_load_mean_work_j,
                    work.support_mean_work_j,
                    0.0,
                    work.kinetic_change_j,
                    0.0,
                    work.total_balance_error_j.abs(),
                    work.mean_work_error_j.abs(),
                    work.mean_internal_work_j,
                    work.parameter_work_j,
                    work.eigenstrain_work_j,
                    work.viscous_dissipation_j,
                ],
            )?;
            push(
                &mut intervals,
                BrittleFaceMotionRecord {
                    start_time_s: accepted_view.time_s,
                    end_time_s: next_time_s,
                    energy,
                },
            )?;
            current = Some(next);
        }
        if let Some(location) = contact {
            let at_contact = current.as_ref().unwrap_or(initial);
            let contact_time_s = at_contact.time_view().time_s;
            if location.time_s != contact_time_s
                || impacts.last().is_some_and(|last| {
                    contact_time_s - last.impact_work().time_s <= controls.approach.time_tolerance_s
                })
            {
                return Err(E::UnresolvedContact {
                    time_s: contact_time_s,
                }
                .into());
            }
            if impacts.len() == controls.maximum_impacts {
                return Err(E::WorkLimit {
                    intervals: intervals.len(),
                    impacts: impacts.len(),
                    trials,
                }
                .into());
            }
            let (motion, hit) = policy
                .impact(
                    at_contact,
                    location,
                    restitution,
                    &held,
                    origin_m,
                    controls.impact,
                    cx,
                )
                .map_err(TimeRunError::Policy)?;
            let event = hit.impact_work();
            if !event.positive_normal_impulse {
                return Err(E::UnresolvedContact {
                    time_s: contact_time_s,
                }
                .into());
            }
            accumulate(
                &mut sums,
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    event.support_work_j,
                    event.kinetic_change_j,
                    event.dissipated_energy_j,
                    event.balance_error_j.abs(),
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            )?;
            push(&mut impacts, hit)?;
            current = Some(motion);
            seed_record = None;
        }
    }
    let motion = current.ok_or(E::Invalid("no physical interval accepted"))?;
    let [
        material,
        fracture,
        dead,
        support,
        impulse_support,
        kinetic,
        contact,
        absolute_balance,
        absolute_mean,
        internal,
        parameter,
        eigenstrain,
        viscous,
    ] = sums.map(Sum::value);
    let final_inventory = motion.time_inventory();
    let storage_change = final_inventory.stored_energy_j - initial_inventory.stored_energy_j;
    let energy = BrittleFaceTimeEnergy {
        stored_before_j: initial_inventory.stored_energy_j,
        stored_after_j: final_inventory.stored_energy_j,
        kinetic_before_j: initial_inventory.kinetic_energy_j,
        kinetic_after_j: final_inventory.kinetic_energy_j,
        kinetic_change_j: kinetic,
        material_work_j: material,
        parameter_work_j: parameter,
        eigenstrain_work_j: eigenstrain,
        fracture_dissipation_j: fracture,
        viscous_dissipation_j: viscous,
        contact_dissipation_j: contact,
        normal_contact_work_j: normal_work.value(),
        absolute_normal_contact_work_j: absolute_normal_work.value(),
        tangential_contact_work_j: tangential_work.value(),
        absolute_tangential_contact_work_j: absolute_tangential_work.value(),
        dead_load_mean_work_j: dead,
        support_mean_work_j: support,
        support_impulse_work_j: impulse_support,
        material_balance_error_j: ((((material + parameter) + eigenstrain) - storage_change)
            - fracture)
            - viscous,
        mean_work_error_j: internal - material,
        absolute_mean_work_error_j: absolute_mean,
        total_balance_error_j: ((((((((dead + support) + impulse_support) + parameter)
            + eigenstrain)
            - storage_change)
            - fracture)
            - viscous)
            - contact)
            - kinetic,
        absolute_local_balance_error_j: absolute_balance,
    };
    for (account, defect, tolerance) in [
        (
            "released-face time material work/storage/fracture",
            energy.material_balance_error_j,
            controls.approach.work.material_energy_tolerance_j,
        ),
        (
            "released-face time cumulative fracture loss",
            (final_inventory.fracture_dissipation_j - initial_inventory.fracture_dissipation_j)
                - fracture,
            controls.approach.work.material_energy_tolerance_j,
        ),
        (
            "released-face time mean-stress work",
            energy.mean_work_error_j,
            controls.approach.work.material_energy_tolerance_j,
        ),
        (
            "released-face time absolute mean-stress work",
            absolute_mean,
            controls.approach.work.material_energy_tolerance_j,
        ),
        (
            "released-face time complete work balance",
            energy.total_balance_error_j,
            controls.total_energy_tolerance_j,
        ),
        (
            "released-face time absolute local work balance",
            absolute_balance,
            controls.total_energy_tolerance_j,
        ),
    ] {
        energy_gate(account, defect, tolerance)?;
    }
    if let Some(c) = sustained_controls {
        energy_gate(
            "released-face time normal contact work",
            energy.normal_contact_work_j,
            c.contact_work_tolerance_j,
        )?;
        energy_gate(
            "released-face time absolute normal contact work",
            energy.absolute_normal_contact_work_j,
            c.contact_work_tolerance_j,
        )?;
        energy_gate(
            "released-face time tangential contact work",
            energy.tangential_contact_work_j,
            c.contact_work_tolerance_j,
        )?;
        energy_gate(
            "released-face time absolute combined contact work",
            energy.absolute_normal_contact_work_j + energy.absolute_tangential_contact_work_j,
            c.contact_work_tolerance_j,
        )?;
    }
    checkpoint(cx)?;
    Ok(NativeFaceTimeStep {
        motion,
        intervals,
        impacts,
        sustained,
        trials,
        numerical_tangent_evaluations: tangent_evaluations,
        energy,
    })
}

// Inspect every paired trace: a later approach takes precedence over an
// earlier resting pair, and no later overlap may be hidden by mode selection.
fn resting_contact(
    motion: &ReleasedFaceMotionView<'_>,
    controls: BrittleFaceApproachControls,
    cx: &Cx<'_>,
) -> Result<bool, BrittleFaceTimeError> {
    let mut resting = false;
    let mut approaching = false;
    for face in motion.released_faces {
        checkpoint(cx)?;
        for vertex in 0..3 {
            let [a, b] = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            if a == b {
                continue;
            }
            let gap = super::dot(
                face.normal,
                core::array::from_fn(|axis| {
                    motion.displacement_m[b][axis] - motion.displacement_m[a][axis]
                }),
            );
            let speed = super::dot(
                face.normal,
                core::array::from_fn(|axis| {
                    motion.velocity_mps[b][axis] - motion.velocity_mps[a][axis]
                }),
            );
            if !gap.is_finite() || !speed.is_finite() {
                return Err(BrittleFaceTimeError::Invalid(
                    "finite paired geometry and speed required",
                ));
            }
            if gap < 0.0 {
                return Err(BrittleFaceTimeError::Approach(
                    BrittleFaceApproachError::Overlap {
                        face: face.id,
                        vertex,
                        gap_m: gap,
                    },
                ));
            }
            if gap <= controls.activation_tolerance_m {
                approaching |= speed < -controls.velocity_tolerance_mps;
                resting |= speed.abs() <= controls.velocity_tolerance_mps;
            }
        }
    }
    Ok(resting && !approaching)
}

fn buffer<T>(capacity: usize) -> Result<Vec<T>, BrittleFaceTimeError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(capacity)
        .map_err(|_| BrittleFaceTimeError::Allocation)?;
    Ok(result)
}
fn push<T>(values: &mut Vec<T>, value: T) -> Result<(), BrittleFaceTimeError> {
    if values.len() == values.capacity() {
        values
            .try_reserve(1)
            .map_err(|_| BrittleFaceTimeError::Allocation)?;
    }
    values.push(value);
    Ok(())
}
fn accumulate(sums: &mut [Sum; 13], values: [f64; 13]) -> Result<(), BrittleFaceTimeError> {
    for (sum, value) in sums.iter_mut().zip(values) {
        sum.add(value);
        if !value.is_finite() || !sum.value().is_finite() {
            return Err(BrittleFaceTimeError::Invalid(
                "complete physical work must be representable",
            ));
        }
    }
    Ok(())
}
