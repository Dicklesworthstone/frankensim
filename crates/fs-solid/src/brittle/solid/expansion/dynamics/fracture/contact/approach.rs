//! Native thermal-fragment motion up to an approaching released-face event.
//!
//! The same bounded contact locator serves both native material owners. Every
//! trial replays this thermal body's accepted state at its CURRENT source
//! snapshots and temperature, through its original midpoint force, inertia and
//! work recovery. Only load/support paths are sampled. No temperature, source
//! data, mechanical history or damaged state is interpolated between trials.
//!
//! The admitted trace is the native affine displacement path with nonnegative
//! paired openings. Its sampled endpoint bracket is not global continuous
//! collision detection. Temperature changes require separately resolved native
//! thermal endpoints; this entry supplies no thermal-path interpolation,
//! persistent contact force, finite sliding or new fragment-pair discovery.

use fs_exec::Cx;

pub use crate::brittle::solid::dynamics::fracture::contact::approach::{
    BrittleFaceApproachControls, BrittleFaceApproachError, BrittleFaceContactLocation,
};
use crate::brittle::solid::dynamics::fracture::contact::approach::{
    NativeFaceApproachTrial, ReleasedFaceMotionView, approach_released_faces,
};
use crate::brittle::solid::expansion::dynamics::ThermalBrittleTetMotion;
use crate::brittle::solid::expansion::dynamics::midpoint::ThermalBrittleTetMidpointStep;
use crate::plastic3::IncrementSettings;

/// Complete thermal prefix or an active approach at the original physical time.
#[derive(Debug)]
pub struct ThermalBrittleFaceApproachStep<'a> {
    /// `None` only for an initially active approaching pair; no step is spent.
    pub advance: Option<ThermalBrittleTetMidpointStep<'a>>,
    /// `None` means the complete requested native interval was accepted.
    pub contact: Option<BrittleFaceContactLocation>,
    pub trials: usize,
    /// Successful native trials, including discarded localization prefixes.
    /// Refused native trials expose no internal fallback count.
    pub numerical_tangent_evaluations: usize,
}

impl<'a> ThermalBrittleTetMotion<'a> {
    /// Locate an approaching contact while retaining the admitted temperature.
    ///
    /// Applied forces are affine endpoint pairs; each shortened trial uses
    /// their actual prefix mean. Prescribed total displacements follow the
    /// affine path from the accepted endpoint to the supplied target. Every
    /// native solve starts from this same accepted mechanical history and the
    /// same current thermal sources, free strain and frozen reference mass.
    ///
    /// Only a completed prefix is returned, with its original thermal midpoint
    /// work accounts and physical endpoint history. Initial contact spends no
    /// physical step; discarded solves spend only the bounded trial allowance.
    /// All negative gaps refuse, including those inside the solid's opening
    /// roundoff tolerance. Only explicit overlap narrows the event bracket;
    /// source, solver, work and cancellation errors remain terminal.
    ///
    /// The returned witness does not apply an impulse or change temperature.
    /// The subsequent simultaneous impact still gathers all active paired rows.
    /// Resolve temperature/load reversals and continuous-trajectory features
    /// through explicit physical intervals; no earliest-collision certificate
    /// or temperature-path convergence order is inferred from bisection.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_toward_released_face_contact(
        &self,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        controls: BrittleFaceApproachControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceApproachStep<'a>, BrittleFaceApproachError> {
        let materials = self.solid.materials;
        let result = approach_released_faces(
            self.released_face_motion_view(),
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            settings,
            controls,
            cx,
            |dt_s, prescribed, forces| {
                self.advance_midpoint(
                    materials,
                    dt_s,
                    prescribed,
                    forces,
                    settings,
                    controls.work,
                    cx,
                )
            },
        )?;
        Ok(ThermalBrittleFaceApproachStep {
            advance: result.advance,
            contact: result.contact,
            trials: result.trials,
            numerical_tangent_evaluations: result.numerical_tangent_evaluations,
        })
    }

    pub(crate) fn released_face_motion_view(&self) -> ReleasedFaceMotionView<'_> {
        ReleasedFaceMotionView {
            displacement_m: &self.solid.solution.displacement_m,
            velocity_mps: &self.velocity_mps,
            released_faces: &self.solid.released_faces,
            time_s: self.time_s,
        }
    }
}

impl NativeFaceApproachTrial for ThermalBrittleTetMidpointStep<'_> {
    fn face_motion(&self) -> ReleasedFaceMotionView<'_> {
        self.motion.released_face_motion_view()
    }

    fn tangent_evaluations(&self) -> usize {
        self.numerical_tangent_evaluations
    }
}
