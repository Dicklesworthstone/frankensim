//! Original interface sources carried through a complete native motion duration.
//!
//! Impact laws use the explicitly selected coefficient role at actual incoming
//! slip. Sustained sticking uses the same original sources' Static role. Source
//! candidates follow accepted physical events only; discarded localization
//! responses and a refused complete duration publish no source successor.

use super::super::sustained::source::{
    BrittleFaceSourcedSustainedFrictionError, BrittleFaceSourcedSustainedFrictionResult,
    solve_released_face_sustained_sourced_frictional_from,
};
use super::super::sustained::{
    BrittleFaceSustainedFrictionControls, BrittleFaceSustainedFrictionError, resting_friction_rows,
};
use super::super::time::{
    BrittleFaceFrictionTimeResult, BrittleFaceTimedFrictionImpact,
    BrittleFaceTimedSustainedFriction, FrictionTimeMotion,
};
use super::super::{BrittleFaceFrictionControls, BrittleFaceFrictionError};
use super::{
    BrittleFaceFrictionSample, BrittleFaceFrictionSource, BrittleFaceSourcedFrictionError,
};
use crate::brittle::solid::dynamics::BrittleTetMotion;
use crate::brittle::solid::dynamics::fracture::contact::BrittleFaceImpactRow;
use crate::brittle::solid::dynamics::fracture::contact::approach::{
    BrittleFaceApproachControls, BrittleFaceContactLocation, ReleasedFaceMotionView,
};
use crate::brittle::solid::dynamics::fracture::contact::sustained::approach::{
    BrittleFaceSustainedApproachError, NativeSustainedApproachTrial,
    approach_released_faces_with_forces,
};
use crate::brittle::solid::dynamics::fracture::contact::sustained::{
    BrittleFaceSustainedControls, SustainedScope,
};
use crate::brittle::solid::dynamics::fracture::contact::time::{
    BrittleFaceTimeControls, BrittleFaceTimeError, FaceTimeSustainedApproach,
    ReleasedFaceTimeMotion, ReleasedFaceTimePolicy, run_released_face_time_with,
};
use crate::brittle::solid::dynamics::midpoint::BrittleTetMidpointEnergy;
use crate::plastic3::IncrementSettings;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::SourceImpulseCoefficientRole;
use fs_contact::sustained::solve::SolveError;
use fs_exec::Cx;
use fs_matdb::MaterialStateId;

/// Complete physical motion and only the accepted original source queries.
#[derive(Debug)]
pub struct BrittleFaceSourcedFrictionTimeResult<Motion, Energy> {
    pub step: BrittleFaceFrictionTimeResult<Motion, Energy>,
    /// Ordered by accepted physical event/interval. The last sample for a row
    /// is its source successor; unqueried rows retain their original input.
    pub sources: Vec<BrittleFaceTimedFrictionSources>,
}

pub type BrittleFaceSourcedFrictionTimeStep<'a> =
    BrittleFaceSourcedFrictionTimeResult<BrittleTetMotion<'a>, BrittleTetMidpointEnergy>;

#[derive(Debug)]
pub enum BrittleFaceSourcedFrictionTimeError {
    Time(BrittleFaceTimeError),
    Approach(BrittleFaceSustainedApproachError),
    Contact(BrittleFaceSustainedFrictionError),
    Source(BrittleFaceSourcedFrictionError),
    Force(BrittleFaceSourcedSustainedFrictionError),
}
impl From<BrittleFaceTimeError> for BrittleFaceSourcedFrictionTimeError {
    fn from(error: BrittleFaceTimeError) -> Self {
        Self::Time(error)
    }
}
impl From<BrittleFaceSustainedApproachError> for BrittleFaceSourcedFrictionTimeError {
    fn from(error: BrittleFaceSustainedApproachError) -> Self {
        Self::Approach(error)
    }
}
impl From<BrittleFaceSustainedFrictionError> for BrittleFaceSourcedFrictionTimeError {
    fn from(error: BrittleFaceSustainedFrictionError) -> Self {
        Self::Contact(error)
    }
}
impl From<BrittleFaceSourcedFrictionError> for BrittleFaceSourcedFrictionTimeError {
    fn from(error: BrittleFaceSourcedFrictionError) -> Self {
        Self::Source(error)
    }
}
impl From<BrittleFaceSourcedSustainedFrictionError> for BrittleFaceSourcedFrictionTimeError {
    fn from(error: BrittleFaceSourcedSustainedFrictionError) -> Self {
        Self::Force(error)
    }
}
impl From<BrittleFaceFrictionError> for BrittleFaceSourcedFrictionTimeError {
    fn from(error: BrittleFaceFrictionError) -> Self {
        Self::Source(error.into())
    }
}
impl core::fmt::Display for BrittleFaceSourcedFrictionTimeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "source-bound released-face friction duration refused: {self:?}"
        )
    }
}
impl std::error::Error for BrittleFaceSourcedFrictionTimeError {}

impl<'a> BrittleTetMotion<'a> {
    /// Advance one complete duration through sourced impact, static load
    /// transfer and release. Each impact uses the caller's explicit role;
    /// sticking uses the same source's Static coefficients at actual rest.
    ///
    /// Complete external conditions, source identities and texture frames are
    /// held for this duration. Every source refresh starts from the last
    /// accepted sample for its complete paired row. Loaded paths must remain
    /// static, so physical travel, slip and finite-step contact temperature do
    /// not advance. Split changing external conditions and source knots between
    /// calls. This does not infer pressure from contact force or allocate impact
    /// loss to heat. Bulk material and reference matter remain with the native
    /// body throughout every accepted and discarded midpoint response.
    ///
    /// Source records are returned only after the entire requested duration
    /// passes the original time, physical mode and work gates. A late refusal
    /// leaves both the original body and all borrowed source histories intact.
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
    ) -> Result<BrittleFaceSourcedFrictionTimeStep<'a>, BrittleFaceSourcedFrictionTimeError> {
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

/// One source/time transaction for both original native material owners.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_sourced_friction_time<M: SourcedFrictionTimeMotion>(
    initial: &M,
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
) -> Result<BrittleFaceSourcedFrictionTimeResult<M, M::Energy>, BrittleFaceSourcedFrictionTimeError>
{
    let mut policy = SourcedFrictionTimePolicy {
        sources: FrictionTimeSources::new(sources),
        sustained,
    };
    let result = run_released_face_time_with(
        initial,
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
    Ok(BrittleFaceSourcedFrictionTimeResult {
        step: BrittleFaceFrictionTimeResult {
            motion: result.motion,
            intervals: result.intervals,
            impacts: result.impacts,
            sustained: result.sustained,
            trials: result.trials,
            numerical_tangent_evaluations: result.numerical_tangent_evaluations,
            energy: result.energy,
        },
        sources: policy.sources.into_records(),
    })
}

struct SourcedFrictionTimePolicy<'b, 's> {
    sources: FrictionTimeSources<'b, 's>,
    sustained: BrittleFaceSustainedFrictionControls,
}

impl<M: SourcedFrictionTimeMotion> ReleasedFaceTimePolicy<M> for SourcedFrictionTimePolicy<'_, '_> {
    type ImpactControls = BrittleFaceFrictionControls;
    type Impact = BrittleFaceTimedFrictionImpact;
    type Sustained = BrittleFaceTimedSustainedFriction<M::Energy>;
    type Error = BrittleFaceSourcedFrictionTimeError;

    fn validate(&mut self, initial: &M, cx: &Cx<'_>) -> Result<(), Self::Error> {
        super::source_directions(
            initial.time_view().released_faces,
            |element| initial.friction_material(element),
            self.sources.original,
            cx,
        )?;
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
        origin_m: [f64; 3],
        controls: Self::ImpactControls,
        cx: &Cx<'_>,
    ) -> Result<(M, Self::Impact), Self::Error> {
        let bindings = self
            .sources
            .bindings(BrittleFaceFrictionSourcePhase::Impact, cx)?;
        let (response, samples) = super::released_face_sourced_frictional_impact(
            motion.frictional_impact_body(),
            |element| motion.friction_material(element),
            restitution,
            &bindings,
            held,
            origin_m,
            controls,
            cx,
        )?;
        let (velocity, record) = BrittleFaceTimedFrictionImpact::from_response(location, response);
        let next =
            motion.friction_velocity_successor(velocity, record.impulse.kinetic_after_j, cx)?;
        self.sources.accept(
            location.time_s,
            location.time_s,
            BrittleFaceFrictionSourcePhase::Impact,
            samples,
            cx,
        )?;
        Ok((next, record))
    }
    fn sustained(
        &mut self,
        motion: &M,
        duration_s: f64,
        prescribed_end_m: &[(usize, f64)],
        force_intervals_n: &[(usize, [f64; 2])],
        seed: Option<(&[BrittleFaceImpactRow], &[[f64; 3]])>,
        settings: IncrementSettings,
        approach: BrittleFaceApproachControls,
        normal: BrittleFaceSustainedControls,
        cx: &Cx<'_>,
    ) -> Result<FaceTimeSustainedApproach<M, Self::Sustained>, Self::Error> {
        let controls = BrittleFaceSustainedFrictionControls {
            contact: normal,
            ..self.sustained
        };
        let bindings = self
            .sources
            .bindings(BrittleFaceFrictionSourcePhase::StaticInterval, cx)?;
        let initial = motion.time_view();
        let directions = super::source_directions(
            initial.released_faces,
            |element| motion.friction_material(element),
            &bindings,
            cx,
        )?;
        let prefix = approach_released_faces_with_forces(
            initial,
            duration_s,
            prescribed_end_m,
            force_intervals_n,
            seed,
            settings,
            approach,
            normal,
            cx,
            || {
                resting_friction_rows(motion.time_view(), &directions, controls, cx)
                    .map_err(Self::Error::from)
            },
            |dt, prescribed, mean, rows, seeds, remaining| {
                motion
                    .sourced_friction_force_response(
                        dt,
                        &bindings,
                        prescribed,
                        mean,
                        Some((rows, seeds)),
                        SustainedScope::InitialResting,
                        settings,
                        BrittleFaceSustainedFrictionControls {
                            contact: remaining,
                            ..controls
                        },
                        cx,
                    )
                    .map_err(Self::Error::from)
            },
            |error| {
                match error {
                BrittleFaceSourcedFrictionTimeError::Force(BrittleFaceSourcedSustainedFrictionError::Contact(
                    BrittleFaceSustainedFrictionError::Sustained(
                        crate::brittle::solid::dynamics::fracture::contact::sustained::BrittleFaceSustainedError::Solve(error)))) => {
                    if let SolveError::ResponseBudget { evaluations, .. } = error.as_ref() { Some(*evaluations) } else { None }
                },
                _ => None,
            }
            },
        )?;
        let advance = if let Some(accepted) = prefix.advance {
            let end_time_s = accepted.step.motion.time_view().time_s;
            self.sources.accept(
                initial.time_s,
                end_time_s,
                BrittleFaceFrictionSourcePhase::StaticInterval,
                accepted.samples,
                cx,
            )?;
            Some(BrittleFaceTimedSustainedFriction::from_step(
                initial.time_s,
                accepted.step,
            ))
        } else {
            None
        };
        Ok(FaceTimeSustainedApproach {
            advance,
            contact: prefix.contact,
            response_evaluations: prefix.response_evaluations,
            numerical_tangent_evaluations: prefix.numerical_tangent_evaluations,
        })
    }
}

/// Original material identity and original source-held midpoint adapter.
pub(crate) trait SourcedFrictionTimeMotion: FrictionTimeMotion {
    fn friction_material(&self, element: usize) -> Option<&MaterialStateId>;
    #[allow(clippy::too_many_arguments)]
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
    >;
}

impl<'a> SourcedFrictionTimeMotion for BrittleTetMotion<'a> {
    fn friction_material(&self, element: usize) -> Option<&MaterialStateId> {
        self.solid
            .materials
            .get(element)
            .map(|source| source.material().material())
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

impl<M: ReleasedFaceTimeMotion, E> NativeSustainedApproachTrial
    for BrittleFaceSourcedSustainedFrictionResult<M, E>
{
    fn face_motion(&self) -> ReleasedFaceMotionView<'_> {
        self.step.motion.time_view()
    }
    fn rows(&self) -> &[BrittleFaceImpactRow] {
        &self.step.rows
    }
    fn physical_forces_n(&self) -> &[[f64; 3]] {
        &self.step.forces_n
    }
    fn response_evaluations(&self) -> usize {
        self.step.solve.response_evaluations
    }
    fn tangent_evaluations(&self) -> usize {
        self.step.numerical_tangent_evaluations
    }
}

/// The physical mode whose accepted original source query is retained.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BrittleFaceFrictionSourcePhase {
    Impact,
    StaticInterval,
}

/// Exactly the active original samples for one accepted event or interval.
#[derive(Clone, Debug)]
pub struct BrittleFaceTimedFrictionSources {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub phase: BrittleFaceFrictionSourcePhase,
    pub samples: Vec<BrittleFaceFrictionSample>,
}

/// References to already retained samples avoid cloning full source bundles at
/// every force probe. The owning time policy is private until total acceptance.
pub(crate) struct FrictionTimeSources<'b, 's> {
    original: &'b [BrittleFaceFrictionSource<'s>],
    records: Vec<BrittleFaceTimedFrictionSources>,
    latest: Vec<(BrittleFaceImpactRow, usize, usize)>,
}

impl<'b, 's> FrictionTimeSources<'b, 's> {
    pub(crate) fn new(original: &'b [BrittleFaceFrictionSource<'s>]) -> Self {
        Self {
            original,
            records: Vec::new(),
            latest: Vec::new(),
        }
    }

    /// Fixed external conditions and source identity survive every stage;
    /// source state comes from the last accepted sample for this paired vertex.
    pub(crate) fn bindings(
        &self,
        phase: BrittleFaceFrictionSourcePhase,
        cx: &Cx<'_>,
    ) -> Result<Vec<BrittleFaceFrictionSource<'_>>, BrittleFaceSourcedFrictionError> {
        let mut result = Vec::new();
        result.try_reserve_exact(self.original.len()).map_err(|_| {
            BrittleFaceSourcedFrictionError::Invalid("source time binding allocation refused")
        })?;
        for source in self.original {
            cx.checkpoint()
                .map_err(|_| BrittleFaceSourcedFrictionError::Cancelled)?;
            let mut next = *source;
            for (vertex, binding) in next.vertices.iter_mut().enumerate() {
                if let Some(&(_, record, sample)) = self
                    .latest
                    .iter()
                    .find(|(row, _, _)| row.face_id == source.face_id && row.vertex == vertex)
                {
                    binding.state = &self.records[record].samples[sample].sample.state;
                }
                if phase == BrittleFaceFrictionSourcePhase::StaticInterval {
                    binding.coefficient_role = SourceImpulseCoefficientRole::Static;
                }
            }
            result.push(next);
        }
        Ok(result)
    }

    /// Accept only the selected physical prefix; all other candidates are
    /// dropped by the original native locator before this method is called.
    pub(crate) fn accept(
        &mut self,
        start_time_s: f64,
        end_time_s: f64,
        phase: BrittleFaceFrictionSourcePhase,
        samples: Vec<BrittleFaceFrictionSample>,
        cx: &Cx<'_>,
    ) -> Result<(), BrittleFaceSourcedFrictionError> {
        cx.checkpoint()
            .map_err(|_| BrittleFaceSourcedFrictionError::Cancelled)?;
        self.records.try_reserve(1).map_err(|_| {
            BrittleFaceSourcedFrictionError::Invalid("source time record allocation refused")
        })?;
        self.latest.try_reserve(samples.len()).map_err(|_| {
            BrittleFaceSourcedFrictionError::Invalid("source time index allocation refused")
        })?;
        let record = self.records.len();
        for (sample, value) in samples.iter().enumerate() {
            cx.checkpoint()
                .map_err(|_| BrittleFaceSourcedFrictionError::Cancelled)?;
            if let Some(entry) = self.latest.iter_mut().find(|(row, _, _)| {
                row.face_id == value.row.face_id && row.vertex == value.row.vertex
            }) {
                if entry.0 != value.row {
                    return Err(BrittleFaceSourcedFrictionError::Invalid(
                        "source time must retain the complete paired-row identity",
                    ));
                }
                *entry = (value.row, record, sample);
            } else {
                self.latest.push((value.row, record, sample));
            }
        }
        self.records.push(BrittleFaceTimedFrictionSources {
            start_time_s,
            end_time_s,
            phase,
            samples,
        });
        Ok(())
    }

    pub(crate) fn into_records(self) -> Vec<BrittleFaceTimedFrictionSources> {
        self.records
    }
}
