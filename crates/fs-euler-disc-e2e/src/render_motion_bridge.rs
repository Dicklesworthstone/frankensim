//! Event-aware bridge from admitted Euler trajectories to render transforms.
//!
//! This module is a visualization adapter. It reconstructs poses with the
//! existing [`TimelineResampler`], maps them into
//! [`fs_render::instances::RigidTransform`], and retains the source
//! [`RenderTrajectoryAuthority`] unchanged. It does not add mechanical
//! resolution, certify the reconstructed motion, or promote simulation
//! evidence into physical truth.

use core::fmt;

use fs_render::instances::{InstanceError, RigidTransform};
use fs_render::motion::{
    MotionTimeError, NormalizedShutterTime, ShotTimeBounds, ShutterConvention, ShutterDistribution,
    ShutterInterval, TimedRay,
};

use crate::render_trajectory::{RenderTrajectory, RenderTrajectoryAuthority};
use crate::timeline_resampling::{
    DeclaredTimelineDiscontinuity, EventEvaluationSide, ExposureEventPolicy, ExposurePartition,
    ResampledTimelineSample, TimelineResampler, TimelineResamplingError,
};

/// Version of the Euler-to-render motion mapping and authority semantics.
pub const EULER_RENDER_MOTION_BRIDGE_VERSION: u16 = 1;

/// How a prepared shutter is partitioned around known timeline events.
#[derive(Clone, Debug, PartialEq)]
pub enum EulerShutterPartition {
    /// A zero-width shutter has one exact time and crosses no open interval.
    Static {
        /// The sole shutter time [s].
        time_s: f64,
    },
    /// A positive-width shutter uses the timeline resampler's event partition.
    EventDelimited(ExposurePartition),
}

/// One event-delimited shutter segment prepared for deterministic sampling.
///
/// Segment-local normalized coordinates map only into this interval. At the
/// closing endpoint of every non-final segment, the bridge selects the event's
/// left limit; the next segment selects the right limit at its opening endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedEulerShutterSegment {
    index: usize,
    shutter: ShutterInterval,
    duration_weight: f64,
    closes_at_interior_event: bool,
}

impl PreparedEulerShutterSegment {
    /// Zero-based index in the prepared shutter's deterministic partition.
    #[must_use]
    pub const fn index(self) -> usize {
        self.index
    }

    /// Segment-local shutter used to construct timed rays.
    #[must_use]
    pub const fn shutter(self) -> ShutterInterval {
        self.shutter
    }

    /// Segment duration divided by the complete positive-width exposure.
    /// Static zero-width shutters carry weight one.
    #[must_use]
    pub const fn duration_weight(self) -> f64 {
        self.duration_weight
    }

    /// Deterministic one-sided policy for a segment-local coordinate.
    #[must_use]
    pub fn event_side_at(self, normalized_time: NormalizedShutterTime) -> EventEvaluationSide {
        if self.closes_at_interior_event && normalized_time.value().to_bits() == 1.0_f64.to_bits() {
            EventEvaluationSide::LeftLimit
        } else {
            EventEvaluationSide::RightLimit
        }
    }
}

/// A resolved shutter bound to one admitted trajectory, event model, and policy.
///
/// Fields are private so the prepared-shutter APIs cannot accidentally bypass
/// event-crossing admission or mix a shutter prepared for another trajectory
/// configuration. [`EulerRenderMotionBridge::sample_at_time`] remains the raw
/// expert API for an already-admitted absolute-time query.
#[derive(Clone, Debug)]
pub struct PreparedEulerShutter<'trajectory> {
    shutter: ShutterInterval,
    partition: EulerShutterPartition,
    event_policy: ExposureEventPolicy,
    declared_discontinuities: Vec<DeclaredTimelineDiscontinuity>,
    source_trajectory: &'trajectory RenderTrajectory,
    source_authority: RenderTrajectoryAuthority,
}

impl PartialEq for PreparedEulerShutter<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.shutter == other.shutter
            && self.partition == other.partition
            && self.event_policy == other.event_policy
            && self.declared_discontinuities == other.declared_discontinuities
            && core::ptr::eq(self.source_trajectory, other.source_trajectory)
            && self.source_authority == other.source_authority
    }
}

impl PreparedEulerShutter<'_> {
    /// Resolved render shutter semantics.
    #[must_use]
    pub const fn shutter(&self) -> ShutterInterval {
        self.shutter
    }

    /// Event-aware exposure partition, including the zero-width case.
    #[must_use]
    pub const fn partition(&self) -> &EulerShutterPartition {
        &self.partition
    }

    /// Explicit policy used while preparing the shutter.
    #[must_use]
    pub const fn event_policy(&self) -> ExposureEventPolicy {
        self.event_policy
    }

    /// Unchanged authority inherited from the admitted source trajectory.
    #[must_use]
    pub const fn source_authority(&self) -> RenderTrajectoryAuthority {
        self.source_authority
    }
}

/// One timeline sample mapped into the render transform convention.
#[derive(Clone, Debug, PartialEq)]
pub struct EulerRenderMotionSample {
    transform: RigidTransform,
    timeline_sample: ResampledTimelineSample,
    source_authority: RenderTrajectoryAuthority,
}

impl EulerRenderMotionSample {
    /// Proper-rigid body-to-world transform for the sampled time.
    #[must_use]
    pub const fn transform(&self) -> RigidTransform {
        self.transform
    }

    /// Timeline reconstruction and discrete event provenance.
    #[must_use]
    pub const fn timeline_sample(&self) -> &ResampledTimelineSample {
        &self.timeline_sample
    }

    /// Absolute trajectory/ray time [s].
    #[must_use]
    pub const fn absolute_time_s(&self) -> f64 {
        self.timeline_sample.time_s
    }

    /// Unchanged authority inherited from the admitted source trajectory.
    #[must_use]
    pub const fn source_authority(&self) -> RenderTrajectoryAuthority {
        self.source_authority
    }
}

/// Structured refusal from shutter preparation or pose-to-transform mapping.
#[derive(Clone, Debug, PartialEq)]
pub enum RenderMotionBridgeError {
    /// The render shutter or its shot bounds were invalid.
    MotionTime(MotionTimeError),
    /// Timeline reconstruction or event-crossing admission refused.
    Timeline(TimelineResamplingError),
    /// The reconstructed pose could not form a proper render transform.
    Transform(InstanceError),
    /// A prepared shutter belonged to a different trajectory configuration.
    PreparedShutterMismatch,
    /// A timed ray did not carry the time implied by the prepared shutter.
    TimedRayShutterMismatch,
    /// A subdivided shutter must be sampled through one explicit segment.
    ShutterSegmentSelectionRequired,
    /// A requested segment index was outside the prepared partition.
    InvalidShutterSegmentIndex {
        /// Requested zero-based segment index.
        index: usize,
        /// Number of segments available in the prepared shutter.
        segment_count: usize,
    },
    /// A timeline partition could not be represented by an exact segment shutter.
    UnrepresentableShutterSegment(usize),
    /// The timeline violated the one-query/one-sample bridge contract.
    TimelineCardinalityMismatch,
}

impl fmt::Display for RenderMotionBridgeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for RenderMotionBridgeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::MotionTime(error) => Some(error),
            Self::Timeline(error) => Some(error),
            Self::Transform(error) => Some(error),
            Self::PreparedShutterMismatch
            | Self::TimedRayShutterMismatch
            | Self::ShutterSegmentSelectionRequired
            | Self::InvalidShutterSegmentIndex { .. }
            | Self::UnrepresentableShutterSegment(_)
            | Self::TimelineCardinalityMismatch => None,
        }
    }
}

impl From<MotionTimeError> for RenderMotionBridgeError {
    fn from(error: MotionTimeError) -> Self {
        Self::MotionTime(error)
    }
}

impl From<TimelineResamplingError> for RenderMotionBridgeError {
    fn from(error: TimelineResamplingError) -> Self {
        Self::Timeline(error)
    }
}

impl From<InstanceError> for RenderMotionBridgeError {
    fn from(error: InstanceError) -> Self {
        Self::Transform(error)
    }
}

/// Binds one admitted Euler trajectory to deterministic render motion.
pub struct EulerRenderMotionBridge<'trajectory> {
    trajectory: &'trajectory RenderTrajectory,
    timeline: TimelineResampler<'trajectory>,
    declared_discontinuities: Vec<DeclaredTimelineDiscontinuity>,
}

impl<'trajectory> EulerRenderMotionBridge<'trajectory> {
    /// Bind an admitted trajectory with no additional continuation seams.
    #[must_use]
    pub fn new(trajectory: &'trajectory RenderTrajectory) -> Self {
        Self {
            trajectory,
            timeline: TimelineResampler::new(trajectory),
            declared_discontinuities: Vec::new(),
        }
    }

    /// Bind an admitted trajectory and source-sample-aligned declared seams.
    pub fn with_declared_discontinuities(
        trajectory: &'trajectory RenderTrajectory,
        declared_discontinuities: Vec<DeclaredTimelineDiscontinuity>,
    ) -> Result<Self, RenderMotionBridgeError> {
        let timeline = TimelineResampler::with_declared_discontinuities(
            trajectory,
            declared_discontinuities.clone(),
        )?;
        Ok(Self {
            trajectory,
            timeline,
            declared_discontinuities,
        })
    }

    /// Resolve a render shutter inside the admitted trajectory and apply the
    /// requested event-crossing policy before any ray samples are evaluated.
    pub fn resolve_shutter(
        &self,
        frame_time_s: f64,
        exposure_duration_s: f64,
        convention: ShutterConvention,
        distribution: ShutterDistribution,
        event_policy: ExposureEventPolicy,
    ) -> Result<PreparedEulerShutter<'trajectory>, RenderMotionBridgeError> {
        let samples = self.trajectory.samples();
        let first = samples[0].input().time_s;
        let last = samples[samples.len() - 1].input().time_s;
        let shot = ShotTimeBounds::try_new(first, last)?;
        let shutter = ShutterInterval::resolve(
            frame_time_s,
            exposure_duration_s,
            convention,
            distribution,
            shot,
        )?;
        let partition = if same_time(shutter.open_s(), shutter.close_s()) {
            EulerShutterPartition::Static {
                time_s: shutter.open_s(),
            }
        } else {
            EulerShutterPartition::EventDelimited(self.timeline.partition_exposure(
                shutter.open_s(),
                shutter.close_s(),
                event_policy,
            )?)
        };
        Ok(PreparedEulerShutter {
            shutter,
            partition,
            event_policy,
            declared_discontinuities: self.declared_discontinuities.clone(),
            source_trajectory: self.trajectory,
            source_authority: self.trajectory.metadata().authority,
        })
    }

    /// Sample one raw absolute trajectory time already admitted by the source.
    ///
    /// This expert resampling API does not apply shutter event policy. Use a
    /// prepared shutter API when evaluating an exposure.
    pub fn sample_at_time(
        &self,
        absolute_time_s: f64,
        event_side: EventEvaluationSide,
    ) -> Result<EulerRenderMotionSample, RenderMotionBridgeError> {
        let mut samples = self.timeline.resample(&[absolute_time_s], event_side)?;
        let timeline_sample = samples
            .pop()
            .ok_or(RenderMotionBridgeError::TimelineCardinalityMismatch)?;
        if !samples.is_empty() {
            return Err(RenderMotionBridgeError::TimelineCardinalityMismatch);
        }
        self.map_timeline_sample(timeline_sample)
    }

    /// Sample an explicit normalized coordinate from a prepared shutter.
    pub fn sample_shutter_coordinate(
        &self,
        prepared: &PreparedEulerShutter<'trajectory>,
        normalized_time: NormalizedShutterTime,
        event_side: EventEvaluationSide,
    ) -> Result<EulerRenderMotionSample, RenderMotionBridgeError> {
        self.validate_prepared(prepared)?;
        self.require_explicit_segment(prepared)?;
        self.sample_at_time(prepared.shutter.time_at(normalized_time), event_side)
    }

    /// Sample the pose associated with a render ray's explicit shutter time.
    ///
    /// The ray must have been constructed from the same resolved shutter. This
    /// prevents a shutter admitted under `Subdivide` or `Refuse` from being
    /// silently paired with a time from another exposure.
    pub fn sample_timed_ray<SpatialRay>(
        &self,
        prepared: &PreparedEulerShutter<'trajectory>,
        ray: &TimedRay<SpatialRay>,
        event_side: EventEvaluationSide,
    ) -> Result<EulerRenderMotionSample, RenderMotionBridgeError> {
        self.validate_prepared(prepared)?;
        self.require_explicit_segment(prepared)?;
        if ray.shutter() != prepared.shutter {
            return Err(RenderMotionBridgeError::TimedRayShutterMismatch);
        }
        let expected_time_s = prepared.shutter.time_at(ray.normalized_time());
        if !same_time(expected_time_s, ray.absolute_time_s()) {
            return Err(RenderMotionBridgeError::TimedRayShutterMismatch);
        }
        self.sample_at_time(ray.absolute_time_s(), event_side)
    }

    /// Resolve one event-delimited segment into a segment-local render shutter.
    ///
    /// Each positive-width segment weight is its duration divided by the full
    /// exposure duration; their mathematical sum is one (subject to binary64
    /// rounding). A static shutter has exactly one segment with weight one.
    pub fn shutter_segment(
        &self,
        prepared: &PreparedEulerShutter<'trajectory>,
        segment_index: usize,
    ) -> Result<PreparedEulerShutterSegment, RenderMotionBridgeError> {
        self.validate_prepared(prepared)?;
        match &prepared.partition {
            EulerShutterPartition::Static { .. } => {
                if segment_index != 0 {
                    return Err(RenderMotionBridgeError::InvalidShutterSegmentIndex {
                        index: segment_index,
                        segment_count: 1,
                    });
                }
                Ok(PreparedEulerShutterSegment {
                    index: 0,
                    shutter: prepared.shutter,
                    duration_weight: 1.0,
                    closes_at_interior_event: false,
                })
            }
            EulerShutterPartition::EventDelimited(partition) => {
                let segment_count = partition.segments.len();
                let segment = partition.segments.get(segment_index).ok_or(
                    RenderMotionBridgeError::InvalidShutterSegmentIndex {
                        index: segment_index,
                        segment_count,
                    },
                )?;
                let shutter = if segment_count == 1 {
                    prepared.shutter
                } else {
                    let duration_s = segment.end_s - segment.start_s;
                    let shot = ShotTimeBounds::try_new(
                        prepared.shutter.open_s(),
                        prepared.shutter.close_s(),
                    )?;
                    ShutterInterval::resolve(
                        segment.start_s,
                        duration_s,
                        ShutterConvention::FrontLoaded,
                        prepared.shutter.distribution(),
                        shot,
                    )?
                };
                if !same_time(shutter.open_s(), segment.start_s)
                    || !same_time(shutter.close_s(), segment.end_s)
                {
                    return Err(RenderMotionBridgeError::UnrepresentableShutterSegment(
                        segment_index,
                    ));
                }
                let duration_weight = shutter.duration_s() / prepared.shutter.duration_s();
                if !duration_weight.is_finite()
                    || !(0.0 < duration_weight && duration_weight <= 1.0)
                {
                    return Err(RenderMotionBridgeError::UnrepresentableShutterSegment(
                        segment_index,
                    ));
                }
                Ok(PreparedEulerShutterSegment {
                    index: segment_index,
                    shutter,
                    duration_weight,
                    closes_at_interior_event: segment_index + 1 < segment_count,
                })
            }
        }
    }

    /// Sample one segment-local coordinate with deterministic event-side rules.
    pub fn sample_segment_coordinate(
        &self,
        prepared: &PreparedEulerShutter<'trajectory>,
        segment_index: usize,
        normalized_time: NormalizedShutterTime,
    ) -> Result<EulerRenderMotionSample, RenderMotionBridgeError> {
        let segment = self.shutter_segment(prepared, segment_index)?;
        self.sample_at_time(
            segment.shutter.time_at(normalized_time),
            segment.event_side_at(normalized_time),
        )
    }

    /// Sample a timed ray constructed from one explicit segment shutter.
    pub fn sample_segment_timed_ray<SpatialRay>(
        &self,
        prepared: &PreparedEulerShutter<'trajectory>,
        segment_index: usize,
        ray: &TimedRay<SpatialRay>,
    ) -> Result<EulerRenderMotionSample, RenderMotionBridgeError> {
        let segment = self.shutter_segment(prepared, segment_index)?;
        if ray.shutter() != segment.shutter {
            return Err(RenderMotionBridgeError::TimedRayShutterMismatch);
        }
        let expected_time_s = segment.shutter.time_at(ray.normalized_time());
        if !same_time(expected_time_s, ray.absolute_time_s()) {
            return Err(RenderMotionBridgeError::TimedRayShutterMismatch);
        }
        self.sample_at_time(
            ray.absolute_time_s(),
            segment.event_side_at(ray.normalized_time()),
        )
    }

    fn validate_prepared(
        &self,
        prepared: &PreparedEulerShutter<'trajectory>,
    ) -> Result<(), RenderMotionBridgeError> {
        if !core::ptr::eq(prepared.source_trajectory, self.trajectory)
            || prepared.source_authority != self.trajectory.metadata().authority
            || prepared.declared_discontinuities != self.declared_discontinuities
        {
            return Err(RenderMotionBridgeError::PreparedShutterMismatch);
        }
        Ok(())
    }

    fn require_explicit_segment(
        &self,
        prepared: &PreparedEulerShutter<'trajectory>,
    ) -> Result<(), RenderMotionBridgeError> {
        if prepared.event_policy == ExposureEventPolicy::Subdivide
            && matches!(
                &prepared.partition,
                EulerShutterPartition::EventDelimited(partition) if partition.segments.len() > 1
            )
        {
            return Err(RenderMotionBridgeError::ShutterSegmentSelectionRequired);
        }
        Ok(())
    }

    fn map_timeline_sample(
        &self,
        timeline_sample: ResampledTimelineSample,
    ) -> Result<EulerRenderMotionSample, RenderMotionBridgeError> {
        let pose = timeline_sample.state.pose();
        let [w, x, y, z] = pose.orientation().components();
        let position = pose.position_world();
        let transform =
            RigidTransform::try_new([x, y, z, w], [position.x, position.y, position.z])?;
        Ok(EulerRenderMotionSample {
            transform,
            timeline_sample,
            source_authority: self.trajectory.metadata().authority,
        })
    }
}

fn same_time(first: f64, second: f64) -> bool {
    let first_bits = first.to_bits();
    let second_bits = second.to_bits();
    first_bits == second_bits || (first_bits << 1 == 0 && second_bits << 1 == 0)
}

#[cfg(test)]
mod tests {
    use core::f64::consts::FRAC_PI_3;

    use fs_blake3::{ContentHash, hash_domain};
    use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
    use fs_render::motion::{
        NormalizedShutterTime, ShutterConvention, ShutterDistribution, TimedRay,
    };

    use super::*;
    use crate::coupled_runner::{ChannelOwnership, ContactTransitionKind};
    use crate::render_trajectory::{
        DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, RenderBaseFrame,
        RenderBaseModeState, RenderChannelAvailability, RenderContactBranch, RenderContactGeometry,
        RenderContactTransition, RenderMassProperties, RenderSampleDisposition,
        RenderSupportFeature, RenderTrajectoryMetadata, RenderTrajectorySampleInput,
        RenderUnitSystem, RenderWorldFrame,
    };
    use crate::timeline_resampling::{
        DeclaredDiscontinuityKind, TimelineEvent, TimelineSampleSource,
    };

    fn identity(label: &str) -> ContentHash {
        hash_domain(
            "org.frankensim.test.render-motion-bridge.v1",
            label.as_bytes(),
        )
    }

    fn mass() -> MassProperties {
        MassProperties::new(1.0, Vec3::ZERO, Vec3::new(0.2, 0.2, 0.2)).unwrap()
    }

    fn orientation() -> UnitQuaternion {
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), FRAC_PI_3).unwrap()
    }

    fn assert_close(left: f64, right: f64) {
        let scale = 1.0_f64.max(left.abs()).max(right.abs());
        assert!((left - right).abs() <= 2.0e-12 * scale, "{left} != {right}");
    }

    fn state(time_s: f64, orientation: UnitQuaternion) -> RigidBodyState {
        RigidBodyState::new(
            Pose::new(Vec3::new(time_s, 2.0 * time_s, 3.0 * time_s), orientation).unwrap(),
            Vec3::new(1.0, 2.0, 3.0),
            Vec3::new(0.1, 0.0, 0.0),
        )
        .unwrap()
    }

    fn sample(
        time_s: f64,
        raw_orientation_wxyz: [f64; 4],
        branch: RenderContactBranch,
        disposition: RenderSampleDisposition,
    ) -> RenderTrajectorySampleInput {
        let orientation = UnitQuaternion::new(
            raw_orientation_wxyz[0],
            raw_orientation_wxyz[1],
            raw_orientation_wxyz[2],
            raw_orientation_wxyz[3],
        )
        .unwrap();
        let state = state(time_s, orientation);
        RenderTrajectorySampleInput {
            interval_start_time_s: 0.0,
            time_s,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            center_of_mass_world_m: state.pose().position_world(),
            orientation_body_to_world: raw_orientation_wxyz,
            linear_momentum_world_kg_m_per_s: state.linear_momentum_world(),
            angular_momentum_body_kg_m2_per_s: state.angular_momentum_body(),
            symmetry_axis_world: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
            contact_branch: branch,
            contact_geometry: (branch == RenderContactBranch::Closed).then_some(
                RenderContactGeometry {
                    point_world_m: Vec3::new(time_s, 2.0 * time_s, 0.0),
                    normal_world: Vec3::new(0.0, 0.0, 1.0),
                    support_feature: RenderSupportFeature::CylinderRim,
                },
            ),
            signed_gap_m: if branch == RenderContactBranch::Closed {
                0.0
            } else {
                1.0e-3
            },
            interval_contact_active: time_s > 0.0 && branch == RenderContactBranch::Closed,
            interval_normal_force_n: if time_s > 0.0 && branch == RenderContactBranch::Closed {
                1.0
            } else {
                0.0
            },
            contact_transitions: Vec::new(),
            base_mode: Some(RenderBaseModeState {
                displacement_m: time_s,
                velocity_m_per_s: 1.0,
            }),
            channels: ChannelOwnership::default(),
            mechanical_energy_j: 1.0,
            energy_defect_j: 0.0,
            qois: DerivedEulerQois::from_state(state, mass(), 0.0).unwrap(),
            disposition,
            terminal_event: None,
        }
    }

    fn metadata(first: &RenderTrajectorySampleInput) -> RenderTrajectoryMetadata {
        let raw = first.orientation_body_to_world;
        let orientation = UnitQuaternion::new(raw[0], raw[1], raw[2], raw[3]).unwrap();
        let initial_state = RigidBodyState::new(
            Pose::new(first.center_of_mass_world_m, orientation).unwrap(),
            first.linear_momentum_world_kg_m_per_s,
            first.angular_momentum_body_kg_m2_per_s,
        )
        .unwrap();
        RenderTrajectoryMetadata {
            schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            specimen_profile_identity: identity("profile"),
            specimen_chart_identity: identity("chart"),
            mass_properties: RenderMassProperties {
                identity: identity("mass"),
                properties: mass(),
            },
            initial_state,
            initial_base_mode: first.base_mode.unwrap(),
            base_model_identity: identity("base"),
            base_frame: RenderBaseFrame {
                origin_world_m: Vec3::ZERO,
                orientation_base_to_world: UnitQuaternion::IDENTITY,
            },
            model_identity: identity("model"),
            channel_availability: RenderChannelAvailability::ALL_AVAILABLE,
            configuration_identity: identity("configuration"),
            configuration_fingerprint: 0x6d6f_7469_6f6e_0001,
            timestep_s: 1.0,
            producer_version: "render-motion-bridge-test-v1".into(),
            applicability: "deterministic visualization reconstruction only".into(),
            no_claims: vec!["does not promote simulation evidence".into()],
            authority: RenderTrajectoryAuthority::SimulationEvidence,
        }
    }

    fn trajectory(raw_second_orientation: [f64; 4], with_event: bool) -> RenderTrajectory {
        let raw_first = orientation().components();
        let first_branch = if with_event {
            RenderContactBranch::Closed
        } else {
            RenderContactBranch::Open
        };
        let first = sample(
            0.0,
            raw_first,
            first_branch,
            RenderSampleDisposition::Continue,
        );
        let mut second = sample(
            1.0,
            raw_second_orientation,
            RenderContactBranch::Open,
            RenderSampleDisposition::HorizonCensored,
        );
        second.interval_start_time_s = first.time_s;
        if with_event {
            second.interval_contact_active = true;
        }
        if with_event {
            second.contact_transitions.push(RenderContactTransition {
                kind: ContactTransitionKind::Opening,
                time_s: 0.5,
                bracket_start_s: 0.49,
                bracket_end_s: 0.51,
            });
        }
        let metadata = metadata(&first);
        RenderTrajectory::try_new(metadata, vec![first, second]).unwrap()
    }

    fn analytic_trajectory() -> RenderTrajectory {
        trajectory(orientation().components(), false)
    }

    fn seam_source_trajectory() -> RenderTrajectory {
        let raw = orientation().components();
        let first = sample(
            0.0,
            raw,
            RenderContactBranch::Open,
            RenderSampleDisposition::Continue,
        );
        let middle = sample(
            0.5,
            raw,
            RenderContactBranch::Open,
            RenderSampleDisposition::Continue,
        );
        let mut last = sample(
            1.0,
            raw,
            RenderContactBranch::Open,
            RenderSampleDisposition::HorizonCensored,
        );
        last.interval_start_time_s = middle.time_s;
        let metadata = metadata(&first);
        RenderTrajectory::try_new(metadata, vec![first, middle, last]).unwrap()
    }

    fn single_sample_trajectory(time_s: f64) -> RenderTrajectory {
        let mut only = sample(
            time_s,
            orientation().components(),
            RenderContactBranch::Open,
            RenderSampleDisposition::HorizonCensored,
        );
        only.interval_start_time_s = time_s;
        let metadata = metadata(&only);
        RenderTrajectory::try_new(metadata, vec![only]).unwrap()
    }

    #[test]
    fn exact_endpoints_map_wxyz_to_xyzw_and_retain_authority() {
        let trajectory = analytic_trajectory();
        let bridge = EulerRenderMotionBridge::new(&trajectory);
        let first = bridge
            .sample_at_time(0.0, EventEvaluationSide::RightLimit)
            .unwrap();
        let last = bridge
            .sample_at_time(1.0, EventEvaluationSide::RightLimit)
            .unwrap();
        let [w, x, y, z] = orientation().components();

        for (observed, expected) in first
            .transform()
            .rotation_xyzw()
            .into_iter()
            .zip([x, y, z, w])
        {
            assert_close(observed, expected);
        }
        assert_eq!(first.transform().translation_m(), [0.0, 0.0, 0.0]);
        assert_eq!(last.transform().translation_m(), [1.0, 2.0, 3.0]);
        assert_eq!(
            first.timeline_sample().source,
            TimelineSampleSource::ExactSample { index: 0 }
        );
        assert_eq!(
            last.timeline_sample().source,
            TimelineSampleSource::ExactSample { index: 1 }
        );
        assert_eq!(
            last.source_authority(),
            RenderTrajectoryAuthority::SimulationEvidence
        );
    }

    #[test]
    fn quaternion_double_cover_has_one_render_transform() {
        let positive_raw = orientation().components();
        let negative_raw = positive_raw.map(|component| -component);
        let positive = trajectory(positive_raw, false);
        let negative = trajectory(negative_raw, false);
        let positive_transform = EulerRenderMotionBridge::new(&positive)
            .sample_at_time(1.0, EventEvaluationSide::RightLimit)
            .unwrap()
            .transform();
        let negative_transform = EulerRenderMotionBridge::new(&negative)
            .sample_at_time(1.0, EventEvaluationSide::RightLimit)
            .unwrap()
            .transform();

        assert_eq!(positive_transform, negative_transform);
    }

    #[test]
    fn event_crossing_shutter_subdivides_or_refuses_explicitly() {
        let trajectory = trajectory(orientation().components(), true);
        let bridge = EulerRenderMotionBridge::new(&trajectory);
        let prepared = bridge
            .resolve_shutter(
                0.5,
                0.5,
                ShutterConvention::Centered,
                ShutterDistribution::UniformCounterV1,
                ExposureEventPolicy::Subdivide,
            )
            .unwrap();
        let EulerShutterPartition::EventDelimited(partition) = prepared.partition() else {
            panic!("positive-width shutter must have an event partition");
        };
        assert_eq!(partition.segments.len(), 2);
        assert!(matches!(
            partition.interior_events.as_slice(),
            [TimelineEvent::Contact(RenderContactTransition {
                kind: ContactTransitionKind::Opening,
                ..
            })]
        ));
        assert_eq!(
            bridge.sample_shutter_coordinate(
                &prepared,
                NormalizedShutterTime::try_new(0.5).unwrap(),
                EventEvaluationSide::RightLimit,
            ),
            Err(RenderMotionBridgeError::ShutterSegmentSelectionRequired)
        );
        let global_ray = TimedRay::at_normalized(
            (),
            prepared.shutter(),
            NormalizedShutterTime::try_new(0.5).unwrap(),
        );
        assert_eq!(
            bridge.sample_timed_ray(&prepared, &global_ray, EventEvaluationSide::RightLimit),
            Err(RenderMotionBridgeError::ShutterSegmentSelectionRequired)
        );

        let before = bridge.shutter_segment(&prepared, 0).unwrap();
        let after = bridge.shutter_segment(&prepared, 1).unwrap();
        assert_eq!(before.index(), 0);
        assert_eq!(after.index(), 1);
        assert_eq!(before.shutter().open_s(), 0.25);
        assert_eq!(before.shutter().close_s(), 0.5);
        assert_eq!(after.shutter().open_s(), 0.5);
        assert_eq!(after.shutter().close_s(), 0.75);
        assert_eq!(before.duration_weight(), 0.5);
        assert_eq!(after.duration_weight(), 0.5);

        let close = NormalizedShutterTime::try_new(1.0).unwrap();
        let open = NormalizedShutterTime::try_new(0.0).unwrap();
        assert_eq!(before.event_side_at(close), EventEvaluationSide::LeftLimit);
        assert_eq!(after.event_side_at(open), EventEvaluationSide::RightLimit);
        let before_event = bridge
            .sample_segment_coordinate(&prepared, 0, close)
            .unwrap();
        let after_event = bridge
            .sample_segment_coordinate(&prepared, 1, open)
            .unwrap();
        assert_eq!(before_event.absolute_time_s(), 0.5);
        assert_eq!(after_event.absolute_time_s(), 0.5);
        assert_eq!(
            before_event.timeline_sample().contact_branch,
            RenderContactBranch::Closed
        );
        assert_eq!(
            after_event.timeline_sample().contact_branch,
            RenderContactBranch::Open
        );

        let before_ray = TimedRay::at_normalized((), before.shutter(), close);
        assert_eq!(
            bridge
                .sample_segment_timed_ray(&prepared, 0, &before_ray)
                .unwrap(),
            before_event
        );
        assert_eq!(
            bridge.sample_segment_timed_ray(&prepared, 1, &before_ray),
            Err(RenderMotionBridgeError::TimedRayShutterMismatch)
        );
        assert_eq!(
            bridge.shutter_segment(&prepared, 2),
            Err(RenderMotionBridgeError::InvalidShutterSegmentIndex {
                index: 2,
                segment_count: 2,
            })
        );
        assert_eq!(
            bridge.resolve_shutter(
                0.5,
                0.5,
                ShutterConvention::Centered,
                ShutterDistribution::UniformCounterV1,
                ExposureEventPolicy::Refuse,
            ),
            Err(RenderMotionBridgeError::Timeline(
                TimelineResamplingError::ExposureSpansEvent
            ))
        );
    }

    #[test]
    fn prepared_shutter_is_bound_to_declared_event_model() {
        let trajectory = seam_source_trajectory();
        let seamless = EulerRenderMotionBridge::new(&trajectory);
        let with_seam = EulerRenderMotionBridge::with_declared_discontinuities(
            &trajectory,
            vec![DeclaredTimelineDiscontinuity {
                time_s: 0.5,
                kind: DeclaredDiscontinuityKind::ContinuationSeam,
            }],
        )
        .unwrap();
        let prepared_without_seam = seamless
            .resolve_shutter(
                0.5,
                0.5,
                ShutterConvention::Centered,
                ShutterDistribution::UniformCounterV1,
                ExposureEventPolicy::Refuse,
            )
            .unwrap();

        assert_eq!(
            with_seam.sample_shutter_coordinate(
                &prepared_without_seam,
                NormalizedShutterTime::try_new(0.25).unwrap(),
                EventEvaluationSide::RightLimit,
            ),
            Err(RenderMotionBridgeError::PreparedShutterMismatch)
        );
    }

    #[test]
    fn prepared_shutter_equality_uses_exact_source_identity() {
        let first_trajectory = analytic_trajectory();
        let cloned_trajectory = first_trajectory.clone();
        assert_eq!(first_trajectory, cloned_trajectory);
        let first_bridge = EulerRenderMotionBridge::new(&first_trajectory);
        let cloned_bridge = EulerRenderMotionBridge::new(&cloned_trajectory);
        let first_prepared = first_bridge
            .resolve_shutter(
                0.5,
                0.25,
                ShutterConvention::Centered,
                ShutterDistribution::UniformCounterV1,
                ExposureEventPolicy::Refuse,
            )
            .unwrap();
        let cloned_prepared = cloned_bridge
            .resolve_shutter(
                0.5,
                0.25,
                ShutterConvention::Centered,
                ShutterDistribution::UniformCounterV1,
                ExposureEventPolicy::Refuse,
            )
            .unwrap();

        assert_ne!(first_prepared, cloned_prepared);
        assert_eq!(
            cloned_bridge.sample_shutter_coordinate(
                &first_prepared,
                NormalizedShutterTime::try_new(0.5).unwrap(),
                EventEvaluationSide::RightLimit,
            ),
            Err(RenderMotionBridgeError::PreparedShutterMismatch)
        );
    }

    #[test]
    fn out_of_range_times_and_shutters_refuse_structurally() {
        let trajectory = analytic_trajectory();
        let bridge = EulerRenderMotionBridge::new(&trajectory);
        assert_eq!(
            bridge.sample_at_time(-0.1, EventEvaluationSide::RightLimit),
            Err(RenderMotionBridgeError::Timeline(
                TimelineResamplingError::QueryOutOfRange {
                    index: 0,
                    time_s: -0.1,
                }
            ))
        );
        assert_eq!(
            bridge.resolve_shutter(
                0.95,
                0.2,
                ShutterConvention::Centered,
                ShutterDistribution::UniformCounterV1,
                ExposureEventPolicy::Subdivide,
            ),
            Err(RenderMotionBridgeError::MotionTime(
                MotionTimeError::ExposureOutsideShot
            ))
        );
    }

    #[test]
    fn positive_exposure_that_collapses_at_absolute_time_resolution_refuses() {
        let frame_time_s = 9_007_199_254_740_992.0;
        let trajectory = single_sample_trajectory(frame_time_s);
        let bridge = EulerRenderMotionBridge::new(&trajectory);

        assert_eq!(
            bridge.resolve_shutter(
                frame_time_s,
                0.5,
                ShutterConvention::FrontLoaded,
                ShutterDistribution::UniformCounterV1,
                ExposureEventPolicy::Refuse,
            ),
            Err(RenderMotionBridgeError::MotionTime(
                MotionTimeError::CollapsedExposure
            ))
        );
    }

    #[test]
    fn zero_width_shutter_reduces_all_coordinates_to_one_pose() {
        let trajectory = analytic_trajectory();
        let bridge = EulerRenderMotionBridge::new(&trajectory);
        let prepared = bridge
            .resolve_shutter(
                0.5,
                0.0,
                ShutterConvention::Centered,
                ShutterDistribution::UniformCounterV1,
                ExposureEventPolicy::Refuse,
            )
            .unwrap();
        assert_eq!(
            prepared.partition(),
            &EulerShutterPartition::Static { time_s: 0.5 }
        );
        let segment = bridge.shutter_segment(&prepared, 0).unwrap();
        assert_eq!(segment.shutter(), prepared.shutter());
        assert_eq!(segment.duration_weight(), 1.0);
        assert_eq!(
            segment.event_side_at(NormalizedShutterTime::try_new(1.0).unwrap()),
            EventEvaluationSide::RightLimit
        );
        let open = bridge
            .sample_shutter_coordinate(
                &prepared,
                NormalizedShutterTime::try_new(0.0).unwrap(),
                EventEvaluationSide::RightLimit,
            )
            .unwrap();
        let close = bridge
            .sample_shutter_coordinate(
                &prepared,
                NormalizedShutterTime::try_new(1.0).unwrap(),
                EventEvaluationSide::RightLimit,
            )
            .unwrap();
        assert_eq!(open, close);
        assert_eq!(open.absolute_time_s().to_bits(), 0.5_f64.to_bits());
    }

    #[test]
    fn timed_ray_mapping_replays_bit_identically() {
        let trajectory = analytic_trajectory();
        let bridge = EulerRenderMotionBridge::new(&trajectory);
        let prepared = bridge
            .resolve_shutter(
                0.5,
                0.8,
                ShutterConvention::Centered,
                ShutterDistribution::StratifiedCounterV1 { strata: 16 },
                ExposureEventPolicy::Subdivide,
            )
            .unwrap();
        let ray = TimedRay::from_sample((), prepared.shutter(), 41, 137);
        let first = bridge
            .sample_timed_ray(&prepared, &ray, EventEvaluationSide::RightLimit)
            .unwrap();
        let replay = bridge
            .sample_timed_ray(&prepared, &ray, EventEvaluationSide::RightLimit)
            .unwrap();

        assert_eq!(first, replay);
        assert_eq!(
            first.absolute_time_s().to_bits(),
            ray.absolute_time_s().to_bits()
        );
        assert_eq!(first.source_authority(), prepared.source_authority());

        let coincident_but_different_shutter = ShutterInterval::resolve(
            0.5,
            0.0,
            ShutterConvention::Centered,
            ShutterDistribution::UniformCounterV1,
            ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
        )
        .unwrap();
        let mismatched_ray = TimedRay::at_normalized(
            (),
            coincident_but_different_shutter,
            NormalizedShutterTime::try_new(0.5).unwrap(),
        );
        assert_eq!(mismatched_ray.absolute_time_s(), 0.5);
        assert_eq!(
            bridge.sample_timed_ray(&prepared, &mismatched_ray, EventEvaluationSide::RightLimit,),
            Err(RenderMotionBridgeError::TimedRayShutterMismatch)
        );
    }
}
