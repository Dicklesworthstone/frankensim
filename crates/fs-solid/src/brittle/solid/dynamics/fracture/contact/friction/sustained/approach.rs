//! Native frictional dwell while another released pair approaches contact.
//!
//! The original constrained-prefix locator carries the full static force
//! response on initially eligible rows. Every other paired opening is checked
//! before a prefix can be accepted. Physical force seeds stay keyed by rows.

use fs_contact::sustained::solve::{MAX_CONTACTS, SolveError};
use fs_exec::Cx;

use super::super::time::FrictionTimeMotion;
use super::super::{BrittleFaceFriction, BrittleFaceFrictionFrame};
use super::{
    BrittleFaceSustainedFrictionControls, BrittleFaceSustainedFrictionError,
    BrittleFaceSustainedFrictionResult, BrittleFaceSustainedFrictionStep,
};
use crate::brittle::solid::dynamics::BrittleTetMotion;
use crate::brittle::solid::dynamics::fracture::contact::BrittleFaceImpactRow;
use crate::brittle::solid::dynamics::fracture::contact::approach::{
    BrittleFaceApproachControls, ReleasedFaceMotionView,
};
use crate::brittle::solid::dynamics::fracture::contact::sustained::approach::{
    BrittleFaceSustainedApproachError, BrittleFaceSustainedApproachResult,
    NativeSustainedApproachTrial, approach_released_faces_with_forces,
};
use crate::brittle::solid::dynamics::fracture::contact::sustained::{
    BrittleFaceSustainedError, SustainedScope, pair_view,
};
use crate::brittle::solid::dynamics::fracture::contact::time::ReleasedFaceTimeMotion;
use crate::brittle::solid::{BrittleTetError, checkpoint};
use crate::plastic3::IncrementSettings;

pub type BrittleFaceSustainedFrictionApproachStep<'a> =
    BrittleFaceSustainedApproachResult<BrittleFaceSustainedFrictionStep<'a>>;

#[derive(Debug)]
pub enum BrittleFaceSustainedFrictionApproachError {
    Approach(BrittleFaceSustainedApproachError),
    Contact(BrittleFaceSustainedFrictionError),
}
impl From<BrittleFaceSustainedApproachError> for BrittleFaceSustainedFrictionApproachError {
    fn from(error: BrittleFaceSustainedApproachError) -> Self {
        Self::Approach(error)
    }
}
impl From<BrittleFaceSustainedFrictionError> for BrittleFaceSustainedFrictionApproachError {
    fn from(error: BrittleFaceSustainedFrictionError) -> Self {
        Self::Contact(error)
    }
}
impl From<BrittleTetError> for BrittleFaceSustainedFrictionApproachError {
    fn from(error: BrittleTetError) -> Self {
        Self::Contact(error.into())
    }
}
impl core::fmt::Display for BrittleFaceSustainedFrictionApproachError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "frictional constrained approach refused: {self:?}")
    }
}
impl std::error::Error for BrittleFaceSustainedFrictionApproachError {}

impl<M: ReleasedFaceTimeMotion, Energy> NativeSustainedApproachTrial
    for BrittleFaceSustainedFrictionResult<M, Energy>
{
    fn face_motion(&self) -> ReleasedFaceMotionView<'_> {
        self.motion.time_view()
    }
    fn rows(&self) -> &[BrittleFaceImpactRow] {
        &self.rows
    }
    fn physical_forces_n(&self) -> &[[f64; 3]] {
        &self.forces_n
    }
    fn response_evaluations(&self) -> usize {
        self.solve.response_evaluations
    }
    fn tangent_evaluations(&self) -> usize {
        self.numerical_tangent_evaluations
    }
}

impl<'a> BrittleTetMotion<'a> {
    /// Locate the next arrival while initially static contacts carry their full
    /// coupled directional forces. Other pairs remain unforced during every
    /// private prefix. Only an independently nonpenetrating native prefix is
    /// returned; no source/force/mode refusal is replaced by free flight.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_toward_released_face_contact_with_sustained_frictional(
        &self,
        duration_s: f64,
        friction: &[BrittleFaceFriction],
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        settings: IncrementSettings,
        approach: BrittleFaceApproachControls,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<
        BrittleFaceSustainedFrictionApproachStep<'a>,
        BrittleFaceSustainedFrictionApproachError,
    > {
        approach_released_faces_with_friction(
            self,
            duration_s,
            friction,
            prescribed_end_m,
            force_intervals_n,
            None,
            settings,
            approach,
            controls,
            cx,
        )
    }
}

/// Shared fixed-law constrained approach for actual mechanical/thermal owners.
#[allow(clippy::too_many_arguments)]
pub(crate) fn approach_released_faces_with_friction<M: FrictionTimeMotion>(
    motion: &M,
    duration_s: f64,
    friction: &[BrittleFaceFriction],
    prescribed_end_m: &[(usize, f64)],
    force_intervals_n: &[(usize, [f64; 2])],
    seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
    settings: IncrementSettings,
    approach: BrittleFaceApproachControls,
    controls: BrittleFaceSustainedFrictionControls,
    cx: &Cx<'_>,
) -> Result<
    BrittleFaceSustainedApproachResult<BrittleFaceSustainedFrictionResult<M, M::Energy>>,
    BrittleFaceSustainedFrictionApproachError,
> {
    let directions: Vec<_> = friction.iter().map(|f| f.tangent_direction).collect();
    approach_released_faces_with_forces(
        motion.time_view(),
        duration_s,
        prescribed_end_m,
        force_intervals_n,
        seed,
        settings,
        approach,
        controls.contact,
        cx,
        || resting_friction_rows(motion.time_view(), &directions, controls, cx).map_err(Into::into),
        |dt, prescribed, forces, rows, seeds, remaining| {
            motion
                .friction_force_response(
                    dt,
                    friction,
                    prescribed,
                    forces,
                    Some((rows, seeds)),
                    SustainedScope::InitialResting,
                    settings,
                    BrittleFaceSustainedFrictionControls {
                        contact: remaining,
                        ..controls
                    },
                    cx,
                )
                .map_err(Into::into)
        },
        |error| {
            if let BrittleFaceSustainedFrictionApproachError::Contact(
                BrittleFaceSustainedFrictionError::Sustained(BrittleFaceSustainedError::Solve(
                    source,
                )),
            ) = error
                && let SolveError::ResponseBudget { evaluations, .. } = source.as_ref()
            {
                Some(*evaluations)
            } else {
                None
            }
        },
    )
}

/// The v10 loaded-eligibility rule, without querying a constitutive source.
pub(crate) fn resting_friction_rows(
    initial: ReleasedFaceMotionView<'_>,
    directions: &[[f64; 3]],
    controls: BrittleFaceSustainedFrictionControls,
    cx: &Cx<'_>,
) -> Result<Vec<BrittleFaceImpactRow>, BrittleFaceSustainedFrictionError> {
    if directions.len() != initial.released_faces.len() {
        return Err(BrittleFaceSustainedError::Invalid(
            "complete retained friction frames required",
        )
        .into());
    }
    let mut rows = Vec::new();
    for (face_index, face) in initial.released_faces.iter().enumerate() {
        checkpoint(cx)?;
        let normal = super::super::unit(face.normal)?;
        let mut frame = None;
        for vertex in 0..3 {
            let nodes = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            if nodes[0] == nodes[1] {
                continue;
            }
            let row = BrittleFaceImpactRow {
                face_id: face.id,
                vertex,
                nodes,
            };
            let geometry = pair_view(initial, row, normal)?;
            let frame = match frame {
                Some(frame) => frame,
                None => {
                    let value =
                        BrittleFaceFrictionFrame::from_unit_normal(normal, directions[face_index])?;
                    frame = Some(value);
                    value
                }
            };
            let velocity = super::super::velocity(initial.velocity_mps, row, frame)?;
            if geometry.gap <= controls.contact.activation_tolerance_m
                && geometry.tangent <= controls.contact.tangential_tolerance_m
                && velocity
                    .iter()
                    .all(|v| v.abs() <= controls.contact.velocity_tolerance_mps)
            {
                if rows.len() >= MAX_CONTACTS
                    || rows.len() >= controls.contact.solve.maximum_contacts
                {
                    return Err(BrittleFaceSustainedError::Invalid(
                        "all initially static friction rows must fit the force budget",
                    )
                    .into());
                }
                rows.push(row);
            }
        }
    }
    Ok(rows)
}
