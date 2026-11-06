//! Deterministic, event-aware resampling of admitted Euler render trajectories.
//!
//! Pose and continuous base state are reconstructed between accepted samples.
//! Contact branches and terminal state are selected discretely; they are never
//! numerically blended. Reconstruction does not add mechanical resolution or
//! increase the authority of the source trajectory.

use core::fmt;

use fs_ivl::Interval;
use fs_mbd::{Pose, RigidBodyState, UnitQuaternion, Vec3};

mod kinematics;

use crate::render_trajectory::{
    MAX_RENDER_TRAJECTORY_SAMPLES, RenderBaseModeState, RenderContactBranch,
    RenderContactTransition, RenderNumericalRefusalReason, RenderSampleDisposition,
    RenderTerminalEvent, RenderTrajectory,
};

/// Version of the interpolation and event-side semantics in this module.
pub const EULER_TIMELINE_RESAMPLER_VERSION: u16 = 1;

/// Reconstruction method recorded by every interpolated sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TimelineInterpolationMethod {
    /// Cubic Hermite translation/base displacement plus shortest-arc SLERP.
    CubicHermiteSlerpV1,
}

/// Which one-sided value to use when a query exactly coincides with an event.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum EventEvaluationSide {
    /// Retain the discrete state immediately before the event.
    LeftLimit,
    /// Apply the event and retain the discrete state immediately after it.
    RightLimit,
}

/// A producer-declared discontinuity not already represented by contact or
/// terminal metadata.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DeclaredDiscontinuityKind {
    /// Boundary between separately executed continuation segments.
    ContinuationSeam,
    /// Other producer-known discontinuity with no stronger interpretation.
    ProducerDeclared,
}

/// Additional event boundary supplied by the trajectory composer.
///
/// Its time must coincide with an accepted source sample. Without that
/// one-sided state boundary, interpolating either side of a declared
/// discontinuity would invent continuity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeclaredTimelineDiscontinuity {
    /// Exact event time [s].
    pub time_s: f64,
    /// Stable discontinuity class.
    pub kind: DeclaredDiscontinuityKind,
}

/// Event metadata retained at its original time and uncertainty bracket.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TimelineEvent {
    /// Localized contact branch transition.
    Contact(RenderContactTransition),
    /// Localized terminal inclination threshold.
    TerminalInclination(RenderTerminalEvent),
    /// Composer-declared boundary with exact time but no inferred uncertainty.
    Declared(DeclaredTimelineDiscontinuity),
}

impl TimelineEvent {
    /// Retained event time [s].
    #[must_use]
    pub const fn time_s(self) -> f64 {
        match self {
            Self::Contact(event) => event.time_s,
            Self::TerminalInclination(event) => event.time_s,
            Self::Declared(event) => event.time_s,
        }
    }
}

/// Provenance of one reconstructed sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum TimelineSampleSource {
    /// The query reproduced an accepted source sample's continuous state.
    ExactSample {
        /// Accepted source sample index.
        index: usize,
    },
    /// The query was reconstructed between two accepted samples.
    Interpolated {
        /// Left accepted sample index.
        left_index: usize,
        /// Right accepted sample index.
        right_index: usize,
        /// Dimensionless normalized interval coordinate.
        alpha: f64,
        /// Explicit reconstruction method/version.
        method: TimelineInterpolationMethod,
    },
}

/// One deterministic timeline query result.
#[derive(Clone, Debug, PartialEq)]
pub struct ResampledTimelineSample {
    /// Requested time [s].
    pub time_s: f64,
    /// Reconstructed continuous rigid-body state.
    pub state: RigidBodyState,
    /// Reconstructed continuous one-mode base state.
    pub base_mode: RenderBaseModeState,
    /// One-sided, non-interpolated contact branch.
    pub contact_branch: RenderContactBranch,
    /// Terminal/censor disposition; non-source instants are `Continue`.
    pub disposition: RenderSampleDisposition,
    /// Source sample or interval that produced this result.
    pub source: TimelineSampleSource,
    /// All event boundaries in the source interval, with original brackets.
    pub interval_events: Vec<TimelineEvent>,
    /// Events exactly coincident with this query.
    pub events_at_query: Vec<TimelineEvent>,
}

/// Continuous upper speeds of the reconstructed pose, not its momentum channels.
///
/// These enclose the selected original interpolation span. They do not certify
/// the mechanical trajectory or add physical resolution to its accepted samples.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TimelinePoseSpeedBounds {
    /// Upper norm of the derivative of world center-of-mass position [m/s].
    pub maximum_translation_speed_m_s: f64,
    /// Upper angular speed of the normalized shortest-arc orientation [rad/s].
    pub maximum_angular_speed_rad_s: f64,
}

/// Policy for an exposure interval containing one or more event boundaries.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ExposureEventPolicy {
    /// Return event-delimited exposure segments.
    Subdivide,
    /// Refuse rather than selecting an implicit cross-event blur convention.
    Refuse,
    /// Permit contact-branch transitions while retaining hard discontinuities.
    ///
    /// This policy is for consumers, such as rigid-geometry rendering, whose
    /// inputs are the continuously reconstructed pose and base state rather
    /// than the discrete contact branch. Terminal events and producer-declared
    /// discontinuities still refuse the exposure.
    AllowContactTransitions,
}

/// One half-open exposure segment, except that the last end is inclusive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExposureSegment {
    /// Segment start [s].
    pub start_s: f64,
    /// Segment end [s].
    pub end_s: f64,
}

/// Event-aware partition of a shutter exposure.
#[derive(Clone, Debug, PartialEq)]
pub struct ExposurePartition {
    /// Contiguous segments delimited by the selected policy.
    pub segments: Vec<ExposureSegment>,
    /// All events strictly inside the requested exposure, in deterministic order.
    pub interior_events: Vec<TimelineEvent>,
}

/// Structured refusal from timeline reconstruction.
#[derive(Clone, Debug, PartialEq)]
pub enum TimelineResamplingError {
    /// No query times were provided.
    EmptyQueries,
    /// The query count exceeded the source trajectory's public resource ceiling.
    TooManyQueries(usize),
    /// A query was NaN or infinite.
    NonFiniteQuery(usize),
    /// Queries were duplicated or not strictly increasing.
    NonIncreasingQuery(usize),
    /// A query attempted to extrapolate beyond accepted source samples.
    QueryOutOfRange {
        /// Query index.
        index: usize,
        /// Refused query time [s].
        time_s: f64,
    },
    /// Declared discontinuities were invalid, duplicated, or out of range.
    InvalidDeclaredDiscontinuity(usize),
    /// An exposure was not finite, ordered, and within the source time range.
    InvalidExposure,
    /// Explicit policy refused a shutter interval that contains an event.
    ExposureSpansEvent,
    /// A reconstructed rigid-body state was not representable.
    InvalidReconstruction(String),
}

impl fmt::Display for TimelineResamplingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for TimelineResamplingError {}

/// Validated view over one admitted trajectory and optional declared seams.
pub struct TimelineResampler<'trajectory> {
    trajectory: &'trajectory RenderTrajectory,
    declared_discontinuities: Vec<DeclaredTimelineDiscontinuity>,
}

impl<'trajectory> TimelineResampler<'trajectory> {
    /// Construct a resampler with no additional producer-declared seams.
    #[must_use]
    pub const fn new(trajectory: &'trajectory RenderTrajectory) -> Self {
        Self {
            trajectory,
            declared_discontinuities: Vec::new(),
        }
    }

    /// Construct a resampler with strictly ordered, in-range declared seams.
    pub fn with_declared_discontinuities(
        trajectory: &'trajectory RenderTrajectory,
        declared_discontinuities: Vec<DeclaredTimelineDiscontinuity>,
    ) -> Result<Self, TimelineResamplingError> {
        let first = trajectory.samples()[0].input().time_s;
        let last = trajectory.samples()[trajectory.samples().len() - 1]
            .input()
            .time_s;
        let mut previous = None;
        for (index, event) in declared_discontinuities.iter().enumerate() {
            if !event.time_s.is_finite()
                || event.time_s < first
                || event.time_s > last
                || previous.is_some_and(|time| event.time_s <= time)
                || !trajectory
                    .samples()
                    .iter()
                    .any(|sample| same_time(sample.input().time_s, event.time_s))
            {
                return Err(TimelineResamplingError::InvalidDeclaredDiscontinuity(index));
            }
            previous = Some(event.time_s);
        }
        Ok(Self {
            trajectory,
            declared_discontinuities,
        })
    }

    /// Resample arbitrary strictly increasing query times without extrapolation.
    pub fn resample(
        &self,
        query_times_s: &[f64],
        event_side: EventEvaluationSide,
    ) -> Result<Vec<ResampledTimelineSample>, TimelineResamplingError> {
        validate_queries(self.trajectory, query_times_s)?;
        query_times_s
            .iter()
            .copied()
            .map(|time_s| self.sample_at(time_s, event_side))
            .collect()
    }

    /// Reconstruct one arbitrary time without allocating a query vector.
    ///
    /// This is the streaming seam for high-rate downstream consumers such as
    /// acoustic observer evaluation. It has exactly the same event-side and
    /// no-extrapolation semantics as [`Self::resample`].
    pub fn sample(
        &self,
        time_s: f64,
        event_side: EventEvaluationSide,
    ) -> Result<ResampledTimelineSample, TimelineResamplingError> {
        validate_queries(self.trajectory, core::slice::from_ref(&time_s))?;
        self.sample_at(time_s, event_side)
    }

    /// Read the original continuous center-of-mass position without event work.
    ///
    /// Uses exactly the translation interpolant of `sample`, including original
    /// knot values and the same finite domain. No quaternion, momentum, base-mode
    /// or event reconstruction is performed, and successful queries allocate
    /// nothing. This is the geometry-only seam for retarded arrival residuals.
    pub fn position(&self, time_s: f64) -> Result<Vec3, TimelineResamplingError> {
        validate_queries(self.trajectory, core::slice::from_ref(&time_s))?;
        let samples = self.trajectory.samples();
        match self.sample_source(time_s)? {
            TimelineSampleSource::ExactSample { index } => {
                Ok(samples[index].state().pose().position_world())
            }
            TimelineSampleSource::Interpolated {
                left_index,
                right_index,
                alpha,
                ..
            } => {
                let position = interpolate_position(
                    samples[left_index].state(),
                    samples[right_index].state(),
                    self.trajectory.metadata().mass_properties.properties.mass(),
                    samples[right_index].input().time_s - samples[left_index].input().time_s,
                    alpha,
                );
                if [position.x, position.y, position.z]
                    .iter()
                    .any(|value| !value.is_finite())
                {
                    return Err(TimelineResamplingError::InvalidReconstruction(
                        "position interpolation produced non-finite geometry".into(),
                    ));
                }
                Ok(position)
            }
        }
    }

    /// Read the original continuous pose without reconstructing other channels.
    ///
    /// Pose is independent of the one-sided discrete contact/event convention.
    /// Original knots and Hermite/SLERP arithmetic are shared with `sample`;
    /// accepted queries allocate nothing and never extrapolate.
    pub fn pose(&self, time_s: f64) -> Result<Pose, TimelineResamplingError> {
        validate_queries(self.trajectory, core::slice::from_ref(&time_s))?;
        let samples = self.trajectory.samples();
        match self.sample_source(time_s)? {
            TimelineSampleSource::ExactSample { index } => Ok(samples[index].state().pose()),
            TimelineSampleSource::Interpolated {
                left_index,
                right_index,
                alpha,
                ..
            } => interpolate_pose(
                samples[left_index].state(),
                samples[right_index].state(),
                self.trajectory.metadata().mass_properties.properties.mass(),
                samples[right_index].input().time_s - samples[left_index].input().time_s,
                alpha,
            ),
        }
    }

    /// Bound pose speed throughout a closed part of one original sample span.
    ///
    /// The left sample index selects the same Hermite/SLERP interpolant as
    /// sampling. The interval must lie inside that span; no extrapolation or
    /// endpoint-momentum approximation is permitted. Work and storage are fixed.
    pub fn pose_speed_bounds(
        &self,
        left_index: usize,
        interval_s: [f64; 2],
    ) -> Result<TimelinePoseSpeedBounds, TimelineResamplingError> {
        let invalid = || TimelineResamplingError::InvalidReconstruction(
            "pose-speed bound requires a representable interval inside one original span".into(),
        );
        let samples = self.trajectory.samples()
            .get(left_index..)
            .and_then(|samples| samples.get(..2))
            .ok_or_else(invalid)?;
        let first_s = samples[0].input().time_s;
        let last_s = samples[1].input().time_s;
        let duration_s = last_s - first_s;
        if interval_s.iter().any(|time| !time.is_finite())
            || interval_s[0] < first_s
            || interval_s[1] > last_s
            || interval_s[0] > interval_s[1]
            || !(duration_s > 0.0 && duration_s.is_finite())
        {
            return Err(invalid());
        }
        let unit = Interval::new(0.0, 1.0);
        let alpha = interval_s.map(|time_s| {
            if same_time(time_s, first_s) {
                Interval::point(0.0)
            } else if same_time(time_s, last_s) {
                Interval::point(1.0)
            } else {
                ((Interval::point(time_s) - Interval::point(first_s))
                    / Interval::point(duration_s))
                    .intersect(unit)
                    .unwrap_or(unit)
            }
        });
        let states = [samples[0].state(), samples[1].state()];
        let inverse_mass = self.trajectory.metadata().mass_properties.properties.mass().recip();
        let components = |value: Vec3| [value.x, value.y, value.z];
        let position = states.map(|state| components(state.pose().position_world()));
        // These are exactly the endpoint tangents used by interpolate_state.
        let velocity = states.map(|state| {
            components(state.linear_momentum_world().scale(inverse_mass))
        });
        let orientation = states.map(|state| state.pose().orientation().components());
        let [maximum_translation_speed_m_s, maximum_angular_speed_rad_s] =
            kinematics::pose_speed_bounds(position, velocity, orientation, duration_s, alpha)
                .ok_or_else(invalid)?;
        Ok(TimelinePoseSpeedBounds {
            maximum_translation_speed_m_s,
            maximum_angular_speed_rad_s,
        })
    }

    /// Partition a shutter interval at every known event boundary.
    pub fn partition_exposure(
        &self,
        open_s: f64,
        close_s: f64,
        policy: ExposureEventPolicy,
    ) -> Result<ExposurePartition, TimelineResamplingError> {
        let samples = self.trajectory.samples();
        let first = samples[0].input().time_s;
        let last = samples[samples.len() - 1].input().time_s;
        if !open_s.is_finite()
            || !close_s.is_finite()
            || open_s >= close_s
            || open_s < first
            || close_s > last
        {
            return Err(TimelineResamplingError::InvalidExposure);
        }
        let interior_events: Vec<_> = self
            .all_events()
            .into_iter()
            .filter(|event| open_s < event.time_s() && event.time_s() < close_s)
            .collect();
        let contains_hard_discontinuity = interior_events
            .iter()
            .any(|event| !matches!(event, TimelineEvent::Contact(_)));
        if (policy == ExposureEventPolicy::Refuse && !interior_events.is_empty())
            || (policy == ExposureEventPolicy::AllowContactTransitions
                && contains_hard_discontinuity)
        {
            return Err(TimelineResamplingError::ExposureSpansEvent);
        }
        let partition_events = if policy == ExposureEventPolicy::Subdivide {
            interior_events.as_slice()
        } else {
            &[]
        };
        let mut segments = Vec::with_capacity(partition_events.len() + 1);
        let mut start_s = open_s;
        for event in partition_events {
            let end_s = event.time_s();
            if end_s > start_s {
                segments.push(ExposureSegment { start_s, end_s });
                start_s = end_s;
            }
        }
        segments.push(ExposureSegment {
            start_s,
            end_s: close_s,
        });
        Ok(ExposurePartition {
            segments,
            interior_events,
        })
    }

    #[allow(
        clippy::too_many_lines,
        reason = "exact and interpolated paths share one event-side decision boundary"
    )]
    fn sample_at(
        &self,
        time_s: f64,
        event_side: EventEvaluationSide,
    ) -> Result<ResampledTimelineSample, TimelineResamplingError> {
        let samples = self.trajectory.samples();
        match self.sample_source(time_s)? {
            TimelineSampleSource::ExactSample { index } => {
                let sample = &samples[index];
                let input = sample.input();
                let interval_events = self.interval_events(index.saturating_sub(1), index);
                let events_at_query = events_at_time(&interval_events, time_s);
                let contact_branch = if events_at_query.is_empty()
                    || event_side == EventEvaluationSide::RightLimit
                {
                    input.contact_branch
                } else {
                    branch_before_events(input.contact_branch, &events_at_query)
                };
                let disposition =
                    disposition_at_event(input.disposition, &events_at_query, event_side, true);
                Ok(ResampledTimelineSample {
                    time_s,
                    state: sample.state(),
                    base_mode: input.base_mode.ok_or_else(|| {
                        TimelineResamplingError::InvalidReconstruction(
                            "admitted source sample lost its base state".into(),
                        )
                    })?,
                    contact_branch,
                    disposition,
                    source: TimelineSampleSource::ExactSample { index },
                    interval_events,
                    events_at_query,
                })
            }
            TimelineSampleSource::Interpolated {
                left_index,
                right_index,
                alpha,
                ..
            } => {
                let left = &samples[left_index];
                let right = &samples[right_index];
                let duration_s = right.input().time_s - left.input().time_s;
                let interval_events = self.interval_events(left_index, right_index);
                let events_at_query = events_at_time(&interval_events, time_s);
                let state = interpolate_state(
                    left.state(),
                    right.state(),
                    self.trajectory.metadata().mass_properties.properties.mass(),
                    duration_s,
                    alpha,
                )?;
                let left_base = left.input().base_mode.ok_or_else(|| {
                    TimelineResamplingError::InvalidReconstruction(
                        "admitted left sample lost its base state".into(),
                    )
                })?;
                let right_base = right.input().base_mode.ok_or_else(|| {
                    TimelineResamplingError::InvalidReconstruction(
                        "admitted right sample lost its base state".into(),
                    )
                })?;
                let base_mode = interpolate_base(left_base, right_base, duration_s, alpha)?;
                let contact_branch = branch_at_time(
                    left.input().contact_branch,
                    &interval_events,
                    time_s,
                    event_side,
                );
                Ok(ResampledTimelineSample {
                    time_s,
                    state,
                    base_mode,
                    contact_branch,
                    disposition: disposition_at_event(
                        right.input().disposition,
                        &events_at_query,
                        event_side,
                        false,
                    ),
                    source: TimelineSampleSource::Interpolated {
                        left_index,
                        right_index,
                        alpha,
                        method: TimelineInterpolationMethod::CubicHermiteSlerpV1,
                    },
                    interval_events,
                    events_at_query,
                })
            }
        }
    }

    // All callers have admitted the finite query domain. Share the exact knot
    // search and normalized interval arithmetic across full and geometry-only
    // queries, including signed zero and unrepresentable interior fractions.
    fn sample_source(&self, time_s: f64) -> Result<TimelineSampleSource, TimelineResamplingError> {
        let samples = self.trajectory.samples();
        match samples.binary_search_by(|sample| {
            let sample_time = sample.input().time_s;
            if same_time(sample_time, time_s) {
                core::cmp::Ordering::Equal
            } else {
                sample_time.total_cmp(&time_s)
            }
        }) {
            Ok(index) => Ok(TimelineSampleSource::ExactSample { index }),
            Err(right_index) => {
                let left_index = right_index - 1;
                let start_s = samples[left_index].input().time_s;
                let duration_s = samples[right_index].input().time_s - start_s;
                let alpha = (time_s - start_s) / duration_s;
                if !(duration_s.is_finite()
                    && duration_s > 0.0
                    && alpha.is_finite()
                    && 0.0 < alpha
                    && alpha < 1.0)
                {
                    return Err(TimelineResamplingError::InvalidReconstruction(
                        "source interval is not representable".into(),
                    ));
                }
                Ok(TimelineSampleSource::Interpolated {
                    left_index,
                    right_index,
                    alpha,
                    method: TimelineInterpolationMethod::CubicHermiteSlerpV1,
                })
            }
        }
    }

    fn interval_events(&self, left_index: usize, right_index: usize) -> Vec<TimelineEvent> {
        let samples = self.trajectory.samples();
        let start_s = samples[left_index].input().time_s;
        let end_s = samples[right_index].input().time_s;
        let mut events: Vec<_> = samples[right_index]
            .input()
            .contact_transitions
            .iter()
            .copied()
            .map(TimelineEvent::Contact)
            .collect();
        if let Some(event) = samples[right_index].input().terminal_event {
            events.push(TimelineEvent::TerminalInclination(event));
        }
        events.extend(
            self.declared_discontinuities
                .iter()
                .copied()
                .filter(|event| start_s <= event.time_s && event.time_s <= end_s)
                .map(TimelineEvent::Declared),
        );
        events.sort_by(|first, second| first.time_s().total_cmp(&second.time_s()));
        events
    }

    fn all_events(&self) -> Vec<TimelineEvent> {
        let mut events = Vec::new();
        for sample in self.trajectory.samples() {
            events.extend(
                sample
                    .input()
                    .contact_transitions
                    .iter()
                    .copied()
                    .map(TimelineEvent::Contact),
            );
            if let Some(event) = sample.input().terminal_event {
                events.push(TimelineEvent::TerminalInclination(event));
            }
        }
        events.extend(
            self.declared_discontinuities
                .iter()
                .copied()
                .map(TimelineEvent::Declared),
        );
        events.sort_by(|first, second| first.time_s().total_cmp(&second.time_s()));
        events
    }
}

fn validate_queries(
    trajectory: &RenderTrajectory,
    query_times_s: &[f64],
) -> Result<(), TimelineResamplingError> {
    if query_times_s.is_empty() {
        return Err(TimelineResamplingError::EmptyQueries);
    }
    if query_times_s.len() > MAX_RENDER_TRAJECTORY_SAMPLES {
        return Err(TimelineResamplingError::TooManyQueries(query_times_s.len()));
    }
    let samples = trajectory.samples();
    let first = samples[0].input().time_s;
    let last = samples[samples.len() - 1].input().time_s;
    for (index, &time_s) in query_times_s.iter().enumerate() {
        if !time_s.is_finite() {
            return Err(TimelineResamplingError::NonFiniteQuery(index));
        }
        if index > 0 && time_s <= query_times_s[index - 1] {
            return Err(TimelineResamplingError::NonIncreasingQuery(index));
        }
        if time_s < first || time_s > last {
            return Err(TimelineResamplingError::QueryOutOfRange { index, time_s });
        }
    }
    Ok(())
}

fn interpolate_position(
    left: RigidBodyState,
    right: RigidBodyState,
    mass_kg: f64,
    duration_s: f64,
    alpha: f64,
) -> Vec3 {
    let left_velocity = left.linear_momentum_world().scale(mass_kg.recip());
    let right_velocity = right.linear_momentum_world().scale(mass_kg.recip());
    hermite(
        left.pose().position_world(),
        left_velocity,
        right.pose().position_world(),
        right_velocity,
        duration_s,
        alpha,
    )
}

fn interpolate_pose(
    left: RigidBodyState,
    right: RigidBodyState,
    mass_kg: f64,
    duration_s: f64,
    alpha: f64,
) -> Result<Pose, TimelineResamplingError> {
    let position = interpolate_position(left, right, mass_kg, duration_s, alpha);
    let orientation = slerp_shortest(left.pose().orientation(), right.pose().orientation(), alpha)?;
    Pose::new(position, orientation)
        .map_err(|error| TimelineResamplingError::InvalidReconstruction(error.to_string()))
}

fn interpolate_state(
    left: RigidBodyState,
    right: RigidBodyState,
    mass_kg: f64,
    duration_s: f64,
    alpha: f64,
) -> Result<RigidBodyState, TimelineResamplingError> {
    let pose = interpolate_pose(left, right, mass_kg, duration_s, alpha)?;
    let linear_momentum = lerp_vec(
        left.linear_momentum_world(),
        right.linear_momentum_world(),
        alpha,
    );
    let angular_momentum = lerp_vec(
        left.angular_momentum_body(),
        right.angular_momentum_body(),
        alpha,
    );
    RigidBodyState::new(pose, linear_momentum, angular_momentum)
        .map_err(|error| TimelineResamplingError::InvalidReconstruction(error.to_string()))
}

fn interpolate_base(
    left: RenderBaseModeState,
    right: RenderBaseModeState,
    duration_s: f64,
    alpha: f64,
) -> Result<RenderBaseModeState, TimelineResamplingError> {
    let displacement_m = hermite_scalar(
        left.displacement_m,
        left.velocity_m_per_s,
        right.displacement_m,
        right.velocity_m_per_s,
        duration_s,
        alpha,
    );
    let velocity_m_per_s = hermite_scalar_derivative(
        left.displacement_m,
        left.velocity_m_per_s,
        right.displacement_m,
        right.velocity_m_per_s,
        duration_s,
        alpha,
    );
    if !displacement_m.is_finite() || !velocity_m_per_s.is_finite() {
        return Err(TimelineResamplingError::InvalidReconstruction(
            "base interpolation produced non-finite state".into(),
        ));
    }
    Ok(RenderBaseModeState {
        displacement_m,
        velocity_m_per_s,
    })
}

fn slerp_shortest(
    left: UnitQuaternion,
    right: UnitQuaternion,
    alpha: f64,
) -> Result<UnitQuaternion, TimelineResamplingError> {
    let first = left.components();
    let (second, dot) = kinematics::shortest_rotation(first, right.components());
    let components: [f64; 4] = if dot > kinematics::LINEAR_ROTATION_DOT {
        core::array::from_fn(|index| first[index] + alpha * (second[index] - first[index]))
    } else {
        let angle = dot.acos();
        let denominator = angle.sin();
        let left_weight = ((1.0 - alpha) * angle).sin() / denominator;
        let right_weight = (alpha * angle).sin() / denominator;
        core::array::from_fn(|index| {
            left_weight.mul_add(first[index], right_weight * second[index])
        })
    };
    UnitQuaternion::new(components[0], components[1], components[2], components[3])
        .map_err(|error| TimelineResamplingError::InvalidReconstruction(error.to_string()))
}

fn hermite(
    left_position: Vec3,
    left_velocity: Vec3,
    right_position: Vec3,
    right_velocity: Vec3,
    duration_s: f64,
    alpha: f64,
) -> Vec3 {
    Vec3::new(
        hermite_scalar(
            left_position.x,
            left_velocity.x,
            right_position.x,
            right_velocity.x,
            duration_s,
            alpha,
        ),
        hermite_scalar(
            left_position.y,
            left_velocity.y,
            right_position.y,
            right_velocity.y,
            duration_s,
            alpha,
        ),
        hermite_scalar(
            left_position.z,
            left_velocity.z,
            right_position.z,
            right_velocity.z,
            duration_s,
            alpha,
        ),
    )
}

fn hermite_scalar(
    left_position: f64,
    left_velocity: f64,
    right_position: f64,
    right_velocity: f64,
    duration_s: f64,
    alpha: f64,
) -> f64 {
    let alpha2 = alpha * alpha;
    let alpha3 = alpha2 * alpha;
    let h00 = 2.0 * alpha3 - 3.0 * alpha2 + 1.0;
    let h10 = alpha3 - 2.0 * alpha2 + alpha;
    let h01 = -2.0 * alpha3 + 3.0 * alpha2;
    let h11 = alpha3 - alpha2;
    h00.mul_add(
        left_position,
        h10.mul_add(
            duration_s * left_velocity,
            h01.mul_add(right_position, h11 * duration_s * right_velocity),
        ),
    )
}

fn hermite_scalar_derivative(
    left_position: f64,
    left_velocity: f64,
    right_position: f64,
    right_velocity: f64,
    duration_s: f64,
    alpha: f64,
) -> f64 {
    let alpha2 = alpha * alpha;
    let dh00 = 6.0 * alpha2 - 6.0 * alpha;
    let dh10 = 3.0 * alpha2 - 4.0 * alpha + 1.0;
    let dh01 = -dh00;
    let dh11 = 3.0 * alpha2 - 2.0 * alpha;
    dh00.mul_add(
        left_position,
        dh10.mul_add(
            duration_s * left_velocity,
            dh01.mul_add(right_position, dh11 * duration_s * right_velocity),
        ),
    ) / duration_s
}

fn lerp_vec(left: Vec3, right: Vec3, alpha: f64) -> Vec3 {
    left.add(right.sub(left).scale(alpha))
}

fn events_at_time(events: &[TimelineEvent], time_s: f64) -> Vec<TimelineEvent> {
    events
        .iter()
        .copied()
        .filter(|event| same_time(event.time_s(), time_s))
        .collect()
}

fn branch_at_time(
    mut branch: RenderContactBranch,
    events: &[TimelineEvent],
    time_s: f64,
    event_side: EventEvaluationSide,
) -> RenderContactBranch {
    for event in events {
        let TimelineEvent::Contact(event) = event else {
            continue;
        };
        if event.time_s < time_s
            || (same_time(event.time_s, time_s) && event_side == EventEvaluationSide::RightLimit)
        {
            branch = transition_result(event.kind);
        }
    }
    branch
}

fn same_time(first: f64, second: f64) -> bool {
    let first_bits = first.to_bits();
    let second_bits = second.to_bits();
    first_bits == second_bits || (first_bits << 1 == 0 && second_bits << 1 == 0)
}

fn disposition_at_event(
    source_disposition: RenderSampleDisposition,
    events_at_query: &[TimelineEvent],
    event_side: EventEvaluationSide,
    exact_source_sample: bool,
) -> RenderSampleDisposition {
    let terminal_event = events_at_query
        .iter()
        .any(|event| matches!(event, TimelineEvent::TerminalInclination(_)));
    let refusing_reimpact = source_disposition
        == RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::ReimpactLimitExceeded,
        )
        && events_at_query.iter().any(|event| {
            matches!(
                event,
                TimelineEvent::Contact(RenderContactTransition {
                    kind: crate::coupled_runner::ContactTransitionKind::Reimpact,
                    ..
                })
            )
        });
    if event_side == EventEvaluationSide::LeftLimit && (terminal_event || refusing_reimpact) {
        RenderSampleDisposition::Continue
    } else if event_side == EventEvaluationSide::RightLimit && terminal_event {
        RenderSampleDisposition::TerminalInclination
    } else if refusing_reimpact || exact_source_sample {
        source_disposition
    } else {
        RenderSampleDisposition::Continue
    }
}

fn branch_before_events(
    mut branch_after: RenderContactBranch,
    events: &[TimelineEvent],
) -> RenderContactBranch {
    for event in events.iter().rev() {
        if let TimelineEvent::Contact(event) = event {
            branch_after = match event.kind {
                crate::coupled_runner::ContactTransitionKind::Opening => {
                    RenderContactBranch::Closed
                }
                crate::coupled_runner::ContactTransitionKind::Reimpact => RenderContactBranch::Open,
            };
        }
    }
    branch_after
}

fn transition_result(kind: crate::coupled_runner::ContactTransitionKind) -> RenderContactBranch {
    match kind {
        crate::coupled_runner::ContactTransitionKind::Opening => RenderContactBranch::Open,
        crate::coupled_runner::ContactTransitionKind::Reimpact => RenderContactBranch::Closed,
    }
}
