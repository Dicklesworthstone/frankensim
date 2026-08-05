//! Source-bound visualization and audio controls for Euler-disc trajectories.
//!
//! Endpoint pose/velocity controls remain point sampled. Wrench, work, and
//! power controls remain attached to their exact accepted intervals. The only
//! downsampler in this module is a whole-interval boxcar: it integrates before
//! decimating and treats contact-event intervals as barriers. It never invents
//! bandwidth, event impulses, or acoustic-frequency authority.

use core::{fmt, num::NonZeroUsize};

use fs_exec::Cx;
use fs_mbd::{Pose, UnitQuaternion, Vec3};

use crate::{
    DerivedEulerQois, RenderBaseModeState, RenderChannelAvailability, RenderContactBranch,
    RenderContactTransition, RenderSampleDisposition, RenderSupportFeature, RenderTrajectory,
    RenderTrajectoryAuthority,
    coupled_runner::{ChannelOwnership, ChannelWrench, ContactTransitionKind},
};

/// Schema version for the in-memory raw control semantics.
pub const EULER_CONTROL_STREAM_SCHEMA_VERSION: u16 = 1;

/// Time semantics of the only downsampling filter implemented here.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum AudioControlFilter {
    /// Duration-weighted box integration over complete accepted intervals,
    /// with eventful intervals emitted alone.
    WholeIntervalBoxcarV1,
}

/// Why a contact event has no amplitude in this control schema.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ContactEventMeasure {
    /// V1 retains only class, time, and a localization bracket. It does not
    /// retain an event-specific impulse or resolved force history.
    TimingOnly,
}

/// One localized contact event on the source trajectory clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ControlContactEvent {
    /// Source interval that owns this event.
    pub source_sample_index: usize,
    /// Opening or reimpact class.
    pub kind: ContactTransitionKind,
    /// Localized event time [s].
    pub time_s: f64,
    /// Inclusive localization bracket start [s].
    pub bracket_start_s: f64,
    /// Inclusive localization bracket end [s].
    pub bracket_end_s: f64,
    /// Width of the retained timing-uncertainty bracket [s].
    pub localization_width_s: f64,
    /// Explicit absence of an admitted impulse magnitude.
    pub measure: ContactEventMeasure,
}

/// Contact coordinates and point velocity at one exact closed endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactFrameCoordinates {
    /// Contact position in the frozen world frame [m].
    pub point_world_m: Vec3,
    /// Contact position relative to the disc center of mass in body axes [m].
    pub point_body_m: Vec3,
    /// Contact position in the displaced reduced-base frame [m].
    pub point_base_m: Vec3,
    /// Base-to-disc contact normal in world axes.
    pub normal_world: Vec3,
    /// Base-to-disc contact normal in body axes.
    pub normal_body: Vec3,
    /// Base-to-disc contact normal in reduced-base axes.
    pub normal_base: Vec3,
    /// Resolved support feature at the exact endpoint.
    pub support_feature: RenderSupportFeature,
    /// Disc material-point velocity at the endpoint contact location [m/s].
    pub disc_point_velocity_world_m_per_s: Vec3,
    /// Reduced-base material-point velocity at the contact location [m/s].
    pub base_point_velocity_world_m_per_s: Vec3,
    /// Disc-minus-base point velocity [m/s].
    pub relative_point_velocity_world_m_per_s: Vec3,
}

/// Exact endpoint controls intended for rendering and scientific overlays.
///
/// These fields are raw SI quantities. Artistic normalization, clamping,
/// exposure, glow, and color mapping belong to a separate consumer-owned style
/// object and are deliberately absent here.
#[derive(Clone, Debug, PartialEq)]
pub struct VisualizationControlPoint {
    /// Index of the exact source sample.
    pub source_sample_index: usize,
    /// Exact accepted endpoint time [s].
    pub time_s: f64,
    /// Exact accepted disc pose.
    pub disc_pose: Pose,
    /// Center-of-mass velocity in world axes [m/s].
    pub center_of_mass_velocity_world_m_per_s: Vec3,
    /// Angular velocity in principal body axes [rad/s].
    pub angular_velocity_body_rad_per_s: Vec3,
    /// Angular velocity in world axes [rad/s].
    pub angular_velocity_world_rad_per_s: Vec3,
    /// Disc symmetry axis in world coordinates.
    pub symmetry_axis_world: Vec3,
    /// Exact reduced-base modal state.
    pub base_mode: RenderBaseModeState,
    /// Origin of the displaced reduced-base frame in world coordinates [m].
    pub displaced_base_origin_world_m: Vec3,
    /// Base-frame orientation, unchanged by the one-mode displacement.
    pub orientation_base_to_world: UnitQuaternion,
    /// Reduced-base point velocity in world coordinates [m/s].
    pub base_velocity_world_m_per_s: Vec3,
    /// Exact contact geometry, available only at a closed retained endpoint.
    pub contact: Option<ContactFrameCoordinates>,
    /// Signed endpoint support gap [m].
    pub signed_gap_m: f64,
    /// Post-interval unilateral branch.
    pub contact_branch: RenderContactBranch,
    /// Audited intrinsic Euler-disc quantities.
    pub qois: DerivedEulerQois,
    /// Terminal/censor/refusal state at this endpoint.
    pub disposition: RenderSampleDisposition,
    /// Audio interval ending at this point, or `None` for a zero-duration
    /// initial point.
    pub preceding_audio_interval_index: Option<usize>,
}

/// Available mean-wrench and signed-work data for one accepted interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AvailableChannelControl {
    /// Duration-weighted mean force in world axes [N].
    pub mean_force_world_n: Vec3,
    /// Duration-weighted mean torque about the source channel's declared
    /// reference point, in world axes [N m].
    pub mean_torque_world_nm: Vec3,
    /// Exact retained signed work [J]. Positive means work into the disc/body
    /// under the source channel convention.
    pub signed_work_j: f64,
    /// Signed interval-mean work rate [W].
    pub signed_mean_work_rate_w: f64,
    /// Integral of mean force over the accepted interval [N s]. This is an
    /// aggregate interval measure, never an event-specific impulse.
    pub force_time_measure_world_n_s: Vec3,
    /// Integral of mean torque over the accepted interval [N m s].
    pub torque_time_measure_world_nm_s: Vec3,
}

/// A channel is either explicitly available or explicitly absent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ChannelControl {
    /// The source declared and retained this channel.
    Available(AvailableChannelControl),
    /// The source explicitly declared the channel unavailable. Numerical zero
    /// is never used to infer this state.
    Unavailable,
}

impl ChannelControl {
    /// Returns the available payload without conflating absence with zero.
    #[must_use]
    pub const fn available(self) -> Option<AvailableChannelControl> {
        match self {
            Self::Available(value) => Some(value),
            Self::Unavailable => None,
        }
    }
}

/// Fixed channel vocabulary shared by raw and coarsened audio controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChannelControlSet {
    /// Gravity body work, retained for accounting rather than sound synthesis.
    pub gravity: ChannelControl,
    /// Aggregate contact body wrench/work. This is not tangential-only work.
    pub contact: ChannelControl,
    /// Reduced rolling-resistance channel.
    pub rolling: ChannelControl,
    /// Reduced-base damping channel, not total contact work into the base.
    pub base: ChannelControl,
    /// Exterior-gas body wrench/work, not relative gas dissipation.
    pub gas: ChannelControl,
}

/// One exact accepted interval of raw audio/control excitation.
#[derive(Clone, Debug, PartialEq)]
pub struct AudioControlInterval {
    /// Source sample whose endpoint closes this interval.
    pub source_sample_index: usize,
    /// Exact interval start [s].
    pub start_time_s: f64,
    /// Exact interval end [s].
    pub end_time_s: f64,
    /// Exact positive duration [s].
    pub duration_s: f64,
    /// True when any accepted subinterval used the closed branch.
    pub interval_contact_active: bool,
    /// `+z_base` component of the duration-weighted mean contact force [N], or
    /// `None` when the contact channel is unavailable.
    pub mean_base_normal_contact_force_n: Option<f64>,
    /// Producer diagnostic evaluated at the midpoint of only the first
    /// accepted subinterval [N]. It is deliberately not named an interval mean.
    pub first_subinterval_midpoint_normal_force_n: f64,
    /// Exact raw channel controls.
    pub channels: ChannelControlSet,
    /// Exact endpoint center-of-mass velocity [m/s].
    pub endpoint_center_of_mass_velocity_world_m_per_s: Vec3,
    /// Exact endpoint angular velocity [rad/s].
    pub endpoint_angular_velocity_world_rad_per_s: Vec3,
    /// Exact endpoint reduced-base velocity in world axes [m/s].
    pub endpoint_base_velocity_world_m_per_s: Vec3,
    /// Localized events retained inside this interval. Their measures are
    /// timing-only; no impulse magnitude is synthesized.
    pub events: Vec<ControlContactEvent>,
}

/// A source-bound pair of point-sampled visualization controls and
/// interval-sampled audio controls.
#[derive(Debug)]
pub struct EulerControlStream<'trajectory> {
    source: &'trajectory RenderTrajectory,
    visualization: Vec<VisualizationControlPoint>,
    audio: Vec<AudioControlInterval>,
    reconciliation: ChannelWorkIntegralChecks,
}

impl<'trajectory> EulerControlStream<'trajectory> {
    /// Derives all controls transactionally from one admitted trajectory.
    ///
    /// Cancellation is polled before every source sample and before
    /// publication. An error returns no partial stream.
    pub fn try_derive(
        source: &'trajectory RenderTrajectory,
        cx: &Cx<'_>,
    ) -> Result<Self, ControlStreamError> {
        cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
        let metadata = source.metadata();
        let mass = metadata.mass_properties.properties;
        let base_orientation = metadata.base_frame.orientation_base_to_world;
        let base_axis_world = base_orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        let mut visualization = Vec::with_capacity(source.samples().len());
        let mut audio = Vec::with_capacity(source.samples().len());
        let mut work_accumulators = ChannelWorkAccumulators::new(metadata.channel_availability);

        for (sample_index, sample) in source.samples().iter().enumerate() {
            cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
            let input = sample.input();
            let state = sample.state();
            let pose = state.pose();
            let center_velocity = state.center_of_mass_velocity_world(mass).map_err(|error| {
                ControlStreamError::DerivedState {
                    sample: sample_index,
                    message: error.to_string(),
                }
            })?;
            let angular_body = mass
                .angular_velocity_body_checked(state.angular_momentum_body())
                .map_err(|error| ControlStreamError::DerivedState {
                    sample: sample_index,
                    message: error.to_string(),
                })?;
            let angular_world = pose.orientation().rotate_body_to_world(angular_body);
            finite_vec(
                angular_world,
                sample_index,
                "angular_velocity_world_rad_per_s",
            )?;
            let base_mode = input
                .base_mode
                .ok_or(ControlStreamError::MissingBaseState(sample_index))?;
            let displaced_base_origin = metadata
                .base_frame
                .origin_world_m
                .add(base_axis_world.scale(base_mode.displacement_m));
            let base_velocity_world = base_axis_world.scale(base_mode.velocity_m_per_s);
            finite_vec(
                displaced_base_origin,
                sample_index,
                "displaced_base_origin_world_m",
            )?;
            finite_vec(
                base_velocity_world,
                sample_index,
                "base_velocity_world_m_per_s",
            )?;

            let contact = input
                .contact_geometry
                .map(|geometry| {
                    let arm_world = geometry.point_world_m.sub(pose.position_world());
                    let point_body = pose.orientation().rotate_world_to_body(arm_world);
                    let point_base = base_orientation
                        .rotate_world_to_body(geometry.point_world_m.sub(displaced_base_origin));
                    let normal_body = pose
                        .orientation()
                        .rotate_world_to_body(geometry.normal_world);
                    let normal_base = base_orientation.rotate_world_to_body(geometry.normal_world);
                    let disc_point_velocity = center_velocity.add(angular_world.cross(arm_world));
                    let relative_velocity = disc_point_velocity.sub(base_velocity_world);
                    for (field, value) in [
                        ("contact.point_body_m", point_body),
                        ("contact.point_base_m", point_base),
                        ("contact.normal_body", normal_body),
                        ("contact.normal_base", normal_base),
                        (
                            "contact.disc_point_velocity_world_m_per_s",
                            disc_point_velocity,
                        ),
                        (
                            "contact.relative_point_velocity_world_m_per_s",
                            relative_velocity,
                        ),
                    ] {
                        finite_vec(value, sample_index, field)?;
                    }
                    Ok(ContactFrameCoordinates {
                        point_world_m: geometry.point_world_m,
                        point_body_m: point_body,
                        point_base_m: point_base,
                        normal_world: geometry.normal_world,
                        normal_body,
                        normal_base,
                        support_feature: geometry.support_feature,
                        disc_point_velocity_world_m_per_s: disc_point_velocity,
                        base_point_velocity_world_m_per_s: base_velocity_world,
                        relative_point_velocity_world_m_per_s: relative_velocity,
                    })
                })
                .transpose()?;

            let duration_s = input.time_s - input.interval_start_time_s;
            let preceding_audio_interval_index = if duration_s > 0.0 {
                let channels = derive_channel_set(
                    input.channels,
                    metadata.channel_availability,
                    duration_s,
                    sample_index,
                )?;
                let mean_base_normal_contact_force_n = channels
                    .contact
                    .available()
                    .map(|channel| channel.mean_force_world_n.dot(base_axis_world));
                if mean_base_normal_contact_force_n.is_some_and(|value| !value.is_finite()) {
                    return Err(ControlStreamError::NonFiniteDerived {
                        sample: sample_index,
                        field: "mean_base_normal_contact_force_n",
                    });
                }
                let events = input
                    .contact_transitions
                    .iter()
                    .map(|event| control_event(sample_index, *event))
                    .collect::<Result<Vec<_>, _>>()?;
                work_accumulators.accumulate(channels, duration_s, sample_index)?;
                let index = audio.len();
                audio.push(AudioControlInterval {
                    source_sample_index: sample_index,
                    start_time_s: input.interval_start_time_s,
                    end_time_s: input.time_s,
                    duration_s,
                    interval_contact_active: input.interval_contact_active,
                    mean_base_normal_contact_force_n,
                    first_subinterval_midpoint_normal_force_n: input.interval_normal_force_n,
                    channels,
                    endpoint_center_of_mass_velocity_world_m_per_s: center_velocity,
                    endpoint_angular_velocity_world_rad_per_s: angular_world,
                    endpoint_base_velocity_world_m_per_s: base_velocity_world,
                    events,
                });
                Some(index)
            } else {
                None
            };

            visualization.push(VisualizationControlPoint {
                source_sample_index: sample_index,
                time_s: input.time_s,
                disc_pose: pose,
                center_of_mass_velocity_world_m_per_s: center_velocity,
                angular_velocity_body_rad_per_s: angular_body,
                angular_velocity_world_rad_per_s: angular_world,
                symmetry_axis_world: input.symmetry_axis_world,
                base_mode,
                displaced_base_origin_world_m: displaced_base_origin,
                orientation_base_to_world: base_orientation,
                base_velocity_world_m_per_s: base_velocity_world,
                contact,
                signed_gap_m: input.signed_gap_m,
                contact_branch: input.contact_branch,
                qois: input.qois,
                disposition: input.disposition,
                preceding_audio_interval_index,
            });
        }
        cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
        let reconciliation = work_accumulators.finish(source.samples().len().saturating_sub(1))?;
        Ok(Self {
            source,
            visualization,
            audio,
            reconciliation,
        })
    }

    /// Exact source trajectory. Durable content identity belongs to the codec
    /// layer rather than this borrowed in-memory binding.
    #[must_use]
    pub const fn source(&self) -> &'trajectory RenderTrajectory {
        self.source
    }

    /// Tests pointer identity against the exact admitted in-memory source.
    #[must_use]
    pub fn is_bound_to(&self, source: &RenderTrajectory) -> bool {
        core::ptr::eq(self.source, source)
    }

    /// Inherited authority ceiling; derivation never promotes it.
    #[must_use]
    pub const fn authority(&self) -> RenderTrajectoryAuthority {
        self.source.metadata().authority
    }

    /// Exact point-sampled rendering controls.
    #[must_use]
    pub fn visualization(&self) -> &[VisualizationControlPoint] {
        &self.visualization
    }

    /// Exact interval-sampled audio/control excitations.
    #[must_use]
    pub fn audio(&self) -> &[AudioControlInterval] {
        &self.audio
    }

    /// Checks that raw mean work rates integrate back to the exact retained
    /// per-channel work.
    #[must_use]
    pub const fn work_integral_checks(&self) -> ChannelWorkIntegralChecks {
        self.reconciliation
    }

    /// Applies duration-weighted whole-interval box filtering before reducing
    /// temporal resolution. Eventful source intervals are output alone and are
    /// never blended across a contact transition.
    pub fn boxcar_coarsen(
        &self,
        intervals_per_bin: NonZeroUsize,
        cx: &Cx<'_>,
    ) -> Result<CoarsenedAudioControls<'trajectory>, ControlStreamError> {
        cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
        if self.audio.is_empty() {
            return Err(ControlStreamError::NoPositiveDurationIntervals);
        }
        let mut bins = Vec::with_capacity(self.audio.len().div_ceil(intervals_per_bin.get()));
        let mut cursor = 0;
        while cursor < self.audio.len() {
            cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
            if !self.audio[cursor].events.is_empty() {
                bins.push(coarsen_group(&self.audio[cursor..=cursor], true, cx)?);
                cursor += 1;
                continue;
            }
            let start = cursor;
            while cursor < self.audio.len()
                && cursor - start < intervals_per_bin.get()
                && self.audio[cursor].events.is_empty()
            {
                cursor += 1;
            }
            bins.push(coarsen_group(&self.audio[start..cursor], false, cx)?);
        }
        cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
        let represented =
            coarsened_work_checks(&bins, self.source.metadata().channel_availability, cx)?;
        let reconciliation = reconcile_against_raw(self.reconciliation, represented)?;
        Ok(CoarsenedAudioControls {
            source: self.source,
            filter: AudioControlFilter::WholeIntervalBoxcarV1,
            intervals_per_bin,
            bins,
            reconciliation,
        })
    }
}

/// One anti-aliased output control bin.
#[derive(Clone, Debug, PartialEq)]
pub struct CoarsenedAudioBin {
    /// First source interval included in this bin.
    pub first_source_sample_index: usize,
    /// Last source interval included in this bin, inclusive.
    pub last_source_sample_index: usize,
    /// Exact bin start [s].
    pub start_time_s: f64,
    /// Exact bin end [s].
    pub end_time_s: f64,
    /// Positive bin duration [s].
    pub duration_s: f64,
    /// True when any source interval used the closed contact branch.
    pub interval_contact_active: bool,
    /// Duration-weighted mean normal contact force [N], or unavailable.
    pub mean_base_normal_contact_force_n: Option<f64>,
    /// Conservatively combined channel controls.
    pub channels: ChannelControlSet,
    /// Events retained without interpolation. Eventful bins contain exactly
    /// one source interval.
    pub events: Vec<ControlContactEvent>,
    /// True when the bin was isolated as an event barrier.
    pub event_barrier: bool,
}

/// Result of deterministic whole-interval anti-alias filtering.
#[derive(Debug)]
pub struct CoarsenedAudioControls<'trajectory> {
    source: &'trajectory RenderTrajectory,
    filter: AudioControlFilter,
    intervals_per_bin: NonZeroUsize,
    bins: Vec<CoarsenedAudioBin>,
    reconciliation: ChannelWorkIntegralChecks,
}

impl<'trajectory> CoarsenedAudioControls<'trajectory> {
    /// Exact trajectory to which these controls remain bound.
    #[must_use]
    pub const fn source(&self) -> &'trajectory RenderTrajectory {
        self.source
    }

    /// Tests pointer identity against the exact admitted in-memory source.
    #[must_use]
    pub fn is_bound_to(&self, source: &RenderTrajectory) -> bool {
        core::ptr::eq(self.source, source)
    }

    /// Declared anti-alias filter.
    #[must_use]
    pub const fn filter(&self) -> AudioControlFilter {
        self.filter
    }

    /// Maximum event-free source intervals requested per bin.
    #[must_use]
    pub const fn intervals_per_bin(&self) -> NonZeroUsize {
        self.intervals_per_bin
    }

    /// Coarsened bins in source-clock order.
    #[must_use]
    pub fn bins(&self) -> &[CoarsenedAudioBin] {
        &self.bins
    }

    /// Reconciliation of coarsened signed work against raw retained work.
    #[must_use]
    pub const fn work_integral_checks(&self) -> ChannelWorkIntegralChecks {
        self.reconciliation
    }
}

/// One signed-work integration comparison [J].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WorkIntegralCheck {
    /// Exact work summed from the reference representation [J].
    pub retained_work_j: f64,
    /// Work reconstructed by integrating the represented mean rates [J].
    pub integrated_work_j: f64,
    /// `integrated_work_j - retained_work_j` [J].
    pub residual_j: f64,
}

/// Per-channel work reconciliation; `None` means explicitly unavailable.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChannelWorkIntegralChecks {
    /// Gravity channel check.
    pub gravity: Option<WorkIntegralCheck>,
    /// Aggregate contact channel check.
    pub contact: Option<WorkIntegralCheck>,
    /// Rolling channel check.
    pub rolling: Option<WorkIntegralCheck>,
    /// Base damping channel check.
    pub base: Option<WorkIntegralCheck>,
    /// Exterior gas body-work check.
    pub gas: Option<WorkIntegralCheck>,
}

impl ChannelWorkIntegralChecks {
    /// Returns whether every available residual is within an explicit finite,
    /// nonnegative absolute tolerance [J].
    pub fn within_tolerance(self, tolerance_j: f64) -> Result<bool, ControlStreamError> {
        if !tolerance_j.is_finite() || tolerance_j < 0.0 {
            return Err(ControlStreamError::InvalidWorkTolerance);
        }
        Ok([
            self.gravity,
            self.contact,
            self.rolling,
            self.base,
            self.gas,
        ]
        .into_iter()
        .flatten()
        .all(|check| check.residual_j.abs() <= tolerance_j))
    }
}

/// Typed refusal from control derivation or anti-aliased coarsening.
#[derive(Clone, Debug, PartialEq)]
pub enum ControlStreamError {
    /// Cancellation was observed before atomic publication.
    Cancelled,
    /// An admitted sample unexpectedly lacked its mandatory reduced-base state.
    MissingBaseState(usize),
    /// A checked rigid-body derivation refused.
    DerivedState {
        /// Source sample index.
        sample: usize,
        /// Upstream diagnostic.
        message: String,
    },
    /// Finite source values produced a non-representable derived control.
    NonFiniteDerived {
        /// Source sample index.
        sample: usize,
        /// Stable semantic field.
        field: &'static str,
    },
    /// Coarsening was requested for a point-only trajectory.
    NoPositiveDurationIntervals,
    /// A work reconciliation tolerance was negative or non-finite.
    InvalidWorkTolerance,
    /// Internally adjacent source intervals were not clock-contiguous.
    NonContiguousIntervals {
        /// Source sample index at which continuity failed.
        interval: usize,
    },
    /// An internal control representation contradicted source availability.
    ChannelAvailabilityMismatch {
        /// Source sample index at which the contradiction was observed.
        sample: usize,
        /// Stable channel name.
        channel: &'static str,
    },
}

impl fmt::Display for ControlStreamError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ControlStreamError {}

fn derive_channel_set(
    channels: ChannelOwnership,
    availability: RenderChannelAvailability,
    duration_s: f64,
    sample: usize,
) -> Result<ChannelControlSet, ControlStreamError> {
    Ok(ChannelControlSet {
        gravity: derive_channel(channels.gravity, availability.gravity, duration_s, sample)?,
        contact: derive_channel(channels.contact, availability.contact, duration_s, sample)?,
        rolling: derive_channel(channels.rolling, availability.rolling, duration_s, sample)?,
        base: derive_channel(channels.base, availability.base, duration_s, sample)?,
        gas: derive_channel(channels.gas, availability.gas, duration_s, sample)?,
    })
}

fn derive_channel(
    channel: ChannelWrench,
    available: bool,
    duration_s: f64,
    sample: usize,
) -> Result<ChannelControl, ControlStreamError> {
    if !available {
        return Ok(ChannelControl::Unavailable);
    }
    let signed_mean_work_rate_w = channel.work_j / duration_s;
    let force_time_measure_world_n_s = channel.force_world_n.scale(duration_s);
    let torque_time_measure_world_nm_s = channel.torque_world_nm.scale(duration_s);
    if !signed_mean_work_rate_w.is_finite() {
        return Err(ControlStreamError::NonFiniteDerived {
            sample,
            field: "signed_mean_work_rate_w",
        });
    }
    finite_vec(
        force_time_measure_world_n_s,
        sample,
        "force_time_measure_world_n_s",
    )?;
    finite_vec(
        torque_time_measure_world_nm_s,
        sample,
        "torque_time_measure_world_nm_s",
    )?;
    Ok(ChannelControl::Available(AvailableChannelControl {
        mean_force_world_n: channel.force_world_n,
        mean_torque_world_nm: channel.torque_world_nm,
        signed_work_j: channel.work_j,
        signed_mean_work_rate_w,
        force_time_measure_world_n_s,
        torque_time_measure_world_nm_s,
    }))
}

fn control_event(
    source_sample_index: usize,
    event: RenderContactTransition,
) -> Result<ControlContactEvent, ControlStreamError> {
    let localization_width_s = event.bracket_end_s - event.bracket_start_s;
    if !localization_width_s.is_finite() {
        return Err(ControlStreamError::NonFiniteDerived {
            sample: source_sample_index,
            field: "contact_event.localization_width_s",
        });
    }
    Ok(ControlContactEvent {
        source_sample_index,
        kind: event.kind,
        time_s: event.time_s,
        bracket_start_s: event.bracket_start_s,
        bracket_end_s: event.bracket_end_s,
        localization_width_s,
        measure: ContactEventMeasure::TimingOnly,
    })
}

fn coarsen_group(
    intervals: &[AudioControlInterval],
    event_barrier: bool,
    cx: &Cx<'_>,
) -> Result<CoarsenedAudioBin, ControlStreamError> {
    let first = &intervals[0];
    let last = &intervals[intervals.len() - 1];
    for pair in intervals.windows(2) {
        cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
        if pair[0].end_time_s.to_bits() != pair[1].start_time_s.to_bits() {
            return Err(ControlStreamError::NonContiguousIntervals {
                interval: pair[1].source_sample_index,
            });
        }
    }
    let duration_s = last.end_time_s - first.start_time_s;
    let channels = aggregate_channel_sets(intervals, duration_s, cx)?;
    let mut weighted_normal_force = 0.0;
    let mut normal_available = true;
    let mut interval_contact_active = false;
    let mut events = Vec::new();
    for interval in intervals {
        cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
        interval_contact_active |= interval.interval_contact_active;
        if let Some(value) = interval.mean_base_normal_contact_force_n {
            weighted_normal_force = value.mul_add(interval.duration_s, weighted_normal_force);
        } else {
            normal_available = false;
        }
        events.extend_from_slice(&interval.events);
    }
    let mean_base_normal_contact_force_n = if normal_available {
        let value = weighted_normal_force / duration_s;
        if !value.is_finite() {
            return Err(ControlStreamError::NonFiniteDerived {
                sample: first.source_sample_index,
                field: "coarsened.mean_base_normal_contact_force_n",
            });
        }
        Some(value)
    } else {
        None
    };
    Ok(CoarsenedAudioBin {
        first_source_sample_index: first.source_sample_index,
        last_source_sample_index: last.source_sample_index,
        start_time_s: first.start_time_s,
        end_time_s: last.end_time_s,
        duration_s,
        interval_contact_active,
        mean_base_normal_contact_force_n,
        channels,
        events,
        event_barrier,
    })
}

fn aggregate_channel_sets(
    intervals: &[AudioControlInterval],
    duration_s: f64,
    cx: &Cx<'_>,
) -> Result<ChannelControlSet, ControlStreamError> {
    Ok(ChannelControlSet {
        gravity: aggregate_channel(intervals, duration_s, |set| set.gravity, cx)?,
        contact: aggregate_channel(intervals, duration_s, |set| set.contact, cx)?,
        rolling: aggregate_channel(intervals, duration_s, |set| set.rolling, cx)?,
        base: aggregate_channel(intervals, duration_s, |set| set.base, cx)?,
        gas: aggregate_channel(intervals, duration_s, |set| set.gas, cx)?,
    })
}

fn aggregate_channel(
    intervals: &[AudioControlInterval],
    duration_s: f64,
    select: impl Fn(ChannelControlSet) -> ChannelControl,
    cx: &Cx<'_>,
) -> Result<ChannelControl, ControlStreamError> {
    if matches!(select(intervals[0].channels), ChannelControl::Unavailable) {
        return Ok(ChannelControl::Unavailable);
    }
    if intervals.len() == 1 {
        return Ok(select(intervals[0].channels));
    }
    let mut force_time = Vec3::ZERO;
    let mut torque_time = Vec3::ZERO;
    let mut signed_work_j = 0.0;
    for interval in intervals {
        cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
        let Some(channel) = select(interval.channels).available() else {
            return Ok(ChannelControl::Unavailable);
        };
        force_time = force_time.add(channel.force_time_measure_world_n_s);
        torque_time = torque_time.add(channel.torque_time_measure_world_nm_s);
        signed_work_j += channel.signed_work_j;
    }
    let sample = intervals[0].source_sample_index;
    let mean_force = force_time.scale(duration_s.recip());
    let mean_torque = torque_time.scale(duration_s.recip());
    let mean_rate = signed_work_j / duration_s;
    if !signed_work_j.is_finite() || !mean_rate.is_finite() {
        return Err(ControlStreamError::NonFiniteDerived {
            sample,
            field: "coarsened.signed_work_or_rate",
        });
    }
    finite_vec(mean_force, sample, "coarsened.mean_force_world_n")?;
    finite_vec(mean_torque, sample, "coarsened.mean_torque_world_nm")?;
    finite_vec(force_time, sample, "coarsened.force_time_measure_world_n_s")?;
    finite_vec(
        torque_time,
        sample,
        "coarsened.torque_time_measure_world_nm_s",
    )?;
    Ok(ChannelControl::Available(AvailableChannelControl {
        mean_force_world_n: mean_force,
        mean_torque_world_nm: mean_torque,
        signed_work_j,
        signed_mean_work_rate_w: mean_rate,
        force_time_measure_world_n_s: force_time,
        torque_time_measure_world_nm_s: torque_time,
    }))
}

fn coarsened_work_checks(
    bins: &[CoarsenedAudioBin],
    availability: RenderChannelAvailability,
    cx: &Cx<'_>,
) -> Result<ChannelWorkIntegralChecks, ControlStreamError> {
    let mut accumulators = ChannelWorkAccumulators::new(availability);
    for bin in bins {
        cx.checkpoint().map_err(|_| ControlStreamError::Cancelled)?;
        accumulators.accumulate(bin.channels, bin.duration_s, bin.last_source_sample_index)?;
    }
    accumulators.finish(bins.last().map_or(0, |bin| bin.last_source_sample_index))
}

#[derive(Clone, Copy, Debug, Default)]
struct WorkAccumulator {
    retained_work_j: f64,
    integrated_work_j: f64,
}

#[derive(Clone, Copy, Debug)]
struct ChannelWorkAccumulators {
    gravity: Option<WorkAccumulator>,
    contact: Option<WorkAccumulator>,
    rolling: Option<WorkAccumulator>,
    base: Option<WorkAccumulator>,
    gas: Option<WorkAccumulator>,
}

impl ChannelWorkAccumulators {
    fn new(availability: RenderChannelAvailability) -> Self {
        Self {
            gravity: availability.gravity.then_some(WorkAccumulator {
                retained_work_j: 0.0,
                integrated_work_j: 0.0,
            }),
            contact: availability.contact.then_some(WorkAccumulator {
                retained_work_j: 0.0,
                integrated_work_j: 0.0,
            }),
            rolling: availability.rolling.then_some(WorkAccumulator {
                retained_work_j: 0.0,
                integrated_work_j: 0.0,
            }),
            base: availability.base.then_some(WorkAccumulator {
                retained_work_j: 0.0,
                integrated_work_j: 0.0,
            }),
            gas: availability.gas.then_some(WorkAccumulator {
                retained_work_j: 0.0,
                integrated_work_j: 0.0,
            }),
        }
    }

    fn accumulate(
        &mut self,
        channels: ChannelControlSet,
        duration_s: f64,
        sample: usize,
    ) -> Result<(), ControlStreamError> {
        accumulate_channel_work(
            &mut self.gravity,
            channels.gravity,
            duration_s,
            sample,
            "gravity",
        )?;
        accumulate_channel_work(
            &mut self.contact,
            channels.contact,
            duration_s,
            sample,
            "contact",
        )?;
        accumulate_channel_work(
            &mut self.rolling,
            channels.rolling,
            duration_s,
            sample,
            "rolling",
        )?;
        accumulate_channel_work(
            &mut self.base,
            channels.base,
            duration_s,
            sample,
            "base",
        )?;
        accumulate_channel_work(
            &mut self.gas,
            channels.gas,
            duration_s,
            sample,
            "gas",
        )
    }

    fn finish(self, sample: usize) -> Result<ChannelWorkIntegralChecks, ControlStreamError> {
        Ok(ChannelWorkIntegralChecks {
            gravity: finish_channel_work(self.gravity, sample)?,
            contact: finish_channel_work(self.contact, sample)?,
            rolling: finish_channel_work(self.rolling, sample)?,
            base: finish_channel_work(self.base, sample)?,
            gas: finish_channel_work(self.gas, sample)?,
        })
    }
}

fn accumulate_channel_work(
    accumulator: &mut Option<WorkAccumulator>,
    channel: ChannelControl,
    duration_s: f64,
    sample: usize,
    channel_name: &'static str,
) -> Result<(), ControlStreamError> {
    let (Some(accumulator), Some(channel)) = (accumulator, channel.available()) else {
        if accumulator.is_none() && channel.available().is_none() {
            return Ok(());
        }
        return Err(ControlStreamError::ChannelAvailabilityMismatch {
            sample,
            channel: channel_name,
        });
    };
    accumulator.retained_work_j += channel.signed_work_j;
    accumulator.integrated_work_j = channel
        .signed_mean_work_rate_w
        .mul_add(duration_s, accumulator.integrated_work_j);
    if !accumulator.retained_work_j.is_finite() || !accumulator.integrated_work_j.is_finite() {
        return Err(ControlStreamError::NonFiniteDerived {
            sample,
            field: "work_integral_accumulator",
        });
    }
    Ok(())
}

fn finish_channel_work(
    accumulator: Option<WorkAccumulator>,
    sample: usize,
) -> Result<Option<WorkIntegralCheck>, ControlStreamError> {
    accumulator
        .map(|accumulator| {
            let residual_j = accumulator.integrated_work_j - accumulator.retained_work_j;
            if !residual_j.is_finite() {
                return Err(ControlStreamError::NonFiniteDerived {
                    sample,
                    field: "work_integral_residual_j",
                });
            }
            Ok(WorkIntegralCheck {
                retained_work_j: accumulator.retained_work_j,
                integrated_work_j: accumulator.integrated_work_j,
                residual_j,
            })
        })
        .transpose()
}

fn reconcile_against_raw(
    raw: ChannelWorkIntegralChecks,
    represented: ChannelWorkIntegralChecks,
) -> Result<ChannelWorkIntegralChecks, ControlStreamError> {
    Ok(ChannelWorkIntegralChecks {
        gravity: reconcile_channel(raw.gravity, represented.gravity)?,
        contact: reconcile_channel(raw.contact, represented.contact)?,
        rolling: reconcile_channel(raw.rolling, represented.rolling)?,
        base: reconcile_channel(raw.base, represented.base)?,
        gas: reconcile_channel(raw.gas, represented.gas)?,
    })
}

fn reconcile_channel(
    raw: Option<WorkIntegralCheck>,
    represented: Option<WorkIntegralCheck>,
) -> Result<Option<WorkIntegralCheck>, ControlStreamError> {
    match (raw, represented) {
        (Some(raw), Some(represented)) => {
            let residual_j = represented.integrated_work_j - raw.retained_work_j;
            if !residual_j.is_finite() {
                return Err(ControlStreamError::NonFiniteDerived {
                    sample: 0,
                    field: "coarsened.work_integral_residual_j",
                });
            }
            Ok(Some(WorkIntegralCheck {
                retained_work_j: raw.retained_work_j,
                integrated_work_j: represented.integrated_work_j,
                residual_j,
            }))
        }
        _ => Ok(None),
    }
}

fn finite_vec(value: Vec3, sample: usize, field: &'static str) -> Result<(), ControlStreamError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(ControlStreamError::NonFiniteDerived { sample, field })
    }
}
