//! Original static-friction sources held through an accepted sticking interval.
//!
//! The original coefficient sampler prepares the interval's fixed law at the
//! actual initial relative speed. Complete native force/motion acceptance is
//! required before the body and source candidates are returned together.

use super::super::source::{
    BrittleFaceFrictionSample, BrittleFaceFrictionSource, BrittleFaceSourcedFrictionError,
    sample_source, source_directions,
};
use super::*;
use crate::brittle::solid::dynamics::fracture::contact::approach::ReleasedFaceMotionView;
use crate::brittle::solid::dynamics::fracture::contact::sustained::{
    NativeFaceSustainedTrial, SustainedScope,
};
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_matdb::MaterialStateId;

/// A complete accepted sticking interval and its prepared original sources.
#[derive(Debug)]
pub struct BrittleFaceSourcedSustainedFrictionResult<Motion, Energy> {
    pub step: BrittleFaceSustainedFrictionResult<Motion, Energy>,
    /// Only initially active, matched, static pairs prepare source candidates.
    /// No accepted physical slip, travel or finite-step temperature is advanced.
    pub samples: Vec<BrittleFaceFrictionSample>,
}

pub type BrittleFaceSourcedSustainedFrictionStep<'a> =
    BrittleFaceSourcedSustainedFrictionResult<BrittleTetMotion<'a>, BrittleTetMidpointEnergy>;

#[derive(Debug)]
pub enum BrittleFaceSourcedSustainedFrictionError {
    Contact(BrittleFaceSustainedFrictionError),
    Source(BrittleFaceSourcedFrictionError),
}
impl From<BrittleFaceSustainedFrictionError> for BrittleFaceSourcedSustainedFrictionError {
    fn from(error: BrittleFaceSustainedFrictionError) -> Self {
        Self::Contact(error)
    }
}
impl From<BrittleFaceSourcedFrictionError> for BrittleFaceSourcedSustainedFrictionError {
    fn from(error: BrittleFaceSourcedFrictionError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for BrittleFaceSourcedSustainedFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "source-bound released-face sticking interval refused: {self:?}"
        )
    }
}
impl std::error::Error for BrittleFaceSourcedSustainedFrictionError {}

impl<'a> BrittleTetMotion<'a> {
    /// Resolve sustained static friction with an explicitly sourced fixed law.
    ///
    /// Bindings follow original face/vertex order and must select the Static
    /// role. Both ordered surface material IDs must match the incident elements.
    /// Each initially eligible pair samples its own actual initial speed through
    /// the original source owner; no source is inferred for an open, separating
    /// or moving pair. Later acquisition of load on such a pair is refused.
    ///
    /// The coefficients and declared external conditions are held throughout
    /// this interval. Loaded pairs must remain static and matched, so accepted
    /// path/slip/texture histories remain at their original values. This does
    /// not evolve the contact's finite-step temperature or resolve changing
    /// friction parameters within the interval. The caller must subdivide such
    /// changes explicitly. No pressure, travel, damping or heating is invented.
    ///
    /// Source queries and every native force trial are candidates. A late source,
    /// material, mode, budget or work refusal publishes neither history nor body.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_sustained_sourced_frictional(
        &self,
        dt_s: f64,
        sources: &[BrittleFaceFrictionSource<'_>],
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleFaceSourcedSustainedFrictionStep<'a>, BrittleFaceSourcedSustainedFrictionError>
    {
        solve_released_face_sustained_sourced_frictional(
            self.released_face_motion_view(),
            |element| {
                self.solid
                    .materials
                    .get(element)
                    .map(|source| source.material().material())
            },
            sources,
            prescribed_m,
            mean_nodal_forces_n,
            controls,
            cx,
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

/// Shared source/body transaction for native mechanical and thermal midpoints.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_released_face_sustained_sourced_frictional<
    'm,
    Step: NativeFaceSustainedTrial,
>(
    old_motion: ReleasedFaceMotionView<'_>,
    material: impl Fn(usize) -> Option<&'m MaterialStateId>,
    sources: &[BrittleFaceFrictionSource<'_>],
    prescribed_m: &[(usize, f64)],
    mean_nodal_forces_n: &[(usize, f64)],
    controls: BrittleFaceSustainedFrictionControls,
    cx: &Cx<'_>,
    candidate: impl FnMut(&[(usize, f64)]) -> Result<Step, BrittleTetDynamicError>,
) -> Result<
    BrittleFaceSourcedSustainedFrictionResult<Step::Motion, Step::Energy>,
    BrittleFaceSourcedSustainedFrictionError,
> {
    solve_released_face_sustained_sourced_frictional_from(
        old_motion,
        material,
        sources,
        prescribed_m,
        mean_nodal_forces_n,
        None,
        SustainedScope::Complete,
        controls,
        cx,
        candidate,
    )
}

/// The same original sampler and force owner, with physical numerical seeds
/// and an initially eligible subset for constrained contact localization.
#[allow(clippy::too_many_arguments)]
pub(crate) fn solve_released_face_sustained_sourced_frictional_from<
    'm,
    Step: NativeFaceSustainedTrial,
>(
    old_motion: ReleasedFaceMotionView<'_>,
    material: impl Fn(usize) -> Option<&'m MaterialStateId>,
    sources: &[BrittleFaceFrictionSource<'_>],
    prescribed_m: &[(usize, f64)],
    mean_nodal_forces_n: &[(usize, f64)],
    seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
    scope: SustainedScope,
    controls: BrittleFaceSustainedFrictionControls,
    cx: &Cx<'_>,
    candidate: impl FnMut(&[(usize, f64)]) -> Result<Step, BrittleTetDynamicError>,
) -> Result<
    BrittleFaceSourcedSustainedFrictionResult<Step::Motion, Step::Energy>,
    BrittleFaceSourcedSustainedFrictionError,
> {
    let directions = source_directions(old_motion.released_faces, material, sources, cx)?;
    for source in sources {
        cx.checkpoint()
            .map_err(|_| BrittleFaceSourcedFrictionError::Cancelled)?;
        if source
            .vertices
            .iter()
            .any(|vertex| vertex.coefficient_role != SourceImpulseCoefficientRole::Static)
        {
            return Err(BrittleFaceSourcedFrictionError::Invalid(
                "sustained sticking requires the explicit original Static coefficient role",
            )
            .into());
        }
    }
    let mut samples = Vec::new();
    let step = solve_released_face_sustained_frictional_from(
        old_motion,
        &directions,
        prescribed_m,
        mean_nodal_forces_n,
        seed,
        scope,
        controls,
        cx,
        |face_index, row, _frame, velocity, _gap| {
            let sample =
                sample_source(sources[face_index].vertices[row.vertex], *row, velocity, cx)?;
            let friction = EllipticFriction {
                coefficients_mu: sample.sample.coefficients_mu,
                principal_axis: sample.sample.state.constitutive_state().principal_axis(),
            };
            samples.push(sample);
            Ok::<_, BrittleFaceSourcedSustainedFrictionError>(friction)
        },
        candidate,
    )?;
    Ok(BrittleFaceSourcedSustainedFrictionResult { step, samples })
}
