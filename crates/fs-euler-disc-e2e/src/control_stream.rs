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
    RenderContactTransition, RenderNormalForceSampling, RenderSampleDisposition,
    RenderSupportFeature, RenderTrajectory, RenderTrajectoryAuthority,
    coupled_runner::{ChannelOwnership, ChannelWrench, ContactTransitionKind},
};

/// Schema version for the in-memory raw control semantics.
pub const EULER_CONTROL_STREAM_SCHEMA_VERSION: u16 = 3;

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
    /// V3 retains only class, time, and a localization bracket. It does not
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

/// Exact visualization endpoints available for an audio/control interval.
///
/// `start_visualization_index == None` identifies retained preroll: the
/// interval's closing state is available, but its opening rigid/base/contact
/// state was not retained by the source trajectory. Consumers must not animate
/// that interval by silently substituting configuration metadata.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct AudioVisualCoverage {
    /// Visualization point at the interval start, when exactly retained.
    pub start_visualization_index: Option<usize>,
    /// Visualization point at the interval end.
    pub end_visualization_index: usize,
}

impl AudioVisualCoverage {
    /// Whether both endpoints needed for synchronized motion are retained.
    #[must_use]
    pub const fn is_fully_bracketed(self) -> bool {
        self.start_visualization_index.is_some()
    }
}

/// Positive-duration clock range covered by both visualization and audio.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AudioVisualHorizon {
    /// First fully bracketed audio interval start [s].
    pub start_time_s: f64,
    /// Last fully bracketed audio interval end [s].
    pub end_time_s: f64,
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
    /// Whether both visual endpoint states exist for this interval.
    pub visual_coverage: AudioVisualCoverage,
    /// True when any accepted subinterval used the closed branch.
    pub interval_contact_active: bool,
    /// `+z_base` component of the duration-weighted mean contact force [N].
    /// This comes from the full contact channel when available, otherwise from
    /// a separately declared normal-load-only source; `None` means neither has
    /// authority.
    pub mean_base_normal_contact_force_n: Option<f64>,
    /// Sampling rule of [`Self::declared_normal_force_n`].
    pub normal_force_sampling: RenderNormalForceSampling,
    /// Producer-declared normal-load scalar [N]. Its exact meaning is carried
    /// by [`Self::normal_force_sampling`]; only a declared interval mean or
    /// applied zero-order hold enters the mean-force sound seam.
    pub declared_normal_force_n: f64,
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
        let requested_controls = source.samples().len();
        let mut visualization = Vec::new();
        visualization
            .try_reserve_exact(requested_controls)
            .map_err(|_| ControlStreamError::Capacity {
                artifact: "visualization controls",
                requested: requested_controls,
            })?;
        let mut audio = Vec::new();
        audio
            .try_reserve_exact(requested_controls)
            .map_err(|_| ControlStreamError::Capacity {
                artifact: "audio controls",
                requested: requested_controls,
            })?;
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
                    &input.channels,
                    metadata.channel_availability,
                    duration_s,
                    sample_index,
                )?;
                let normal_force_sampling = metadata.channel_availability.normal_force_sampling;
                let mean_base_normal_contact_force_n = channels
                    .contact
                    .available()
                    .map(|contact| contact.mean_force_world_n.dot(base_axis_world))
                    .or_else(|| {
                        matches!(
                            normal_force_sampling,
                            RenderNormalForceSampling::IntervalMean
                                | RenderNormalForceSampling::AppliedSubstepZeroOrderHold
                        )
                        .then_some(input.interval_normal_force_n)
                    });
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
                work_accumulators.accumulate(&channels, duration_s, sample_index)?;
                let index = audio.len();
                audio.push(AudioControlInterval {
                    source_sample_index: sample_index,
                    start_time_s: input.interval_start_time_s,
                    end_time_s: input.time_s,
                    duration_s,
                    visual_coverage: AudioVisualCoverage {
                        start_visualization_index: sample_index.checked_sub(1),
                        end_visualization_index: sample_index,
                    },
                    interval_contact_active: input.interval_contact_active,
                    mean_base_normal_contact_force_n,
                    normal_force_sampling,
                    declared_normal_force_n: input.interval_normal_force_n,
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

    /// Audio intervals whose opening and closing visualization states are both
    /// exact retained points. At most the first raw interval is omitted as
    /// endpoint-only preroll.
    #[must_use]
    pub fn fully_bracketed_audio(&self) -> &[AudioControlInterval] {
        let first = usize::from(
            self.audio
                .first()
                .is_some_and(|interval| !interval.visual_coverage.is_fully_bracketed()),
        );
        &self.audio[first..]
    }

    /// Common positive-duration visualization/audio clock range, or `None`
    /// when the source contains only a point or endpoint-only preroll.
    #[must_use]
    pub fn audio_visual_horizon(&self) -> Option<AudioVisualHorizon> {
        let synchronized = self.fully_bracketed_audio();
        Some(AudioVisualHorizon {
            start_time_s: synchronized.first()?.start_time_s,
            end_time_s: synchronized.last()?.end_time_s,
        })
    }

    /// Checks that raw mean work rates integrate back to the exact retained
    /// per-channel work.
    #[must_use]
    pub const fn work_integral_checks(&self) -> ChannelWorkIntegralChecks {
        self.reconciliation
    }

    /// Applies duration-weighted whole-interval box filtering before reducing
    /// temporal resolution. Eventful source intervals are output alone and are
    /// never blended across a contact transition. This prefilter mitigates
    /// aliasing but does not claim band-limited sample-rate conversion.
    pub fn boxcar_coarsen(
        &self,
        intervals_per_bin: NonZeroUsize,
        cx: &Cx<'_>,
    ) -> Result<CoarsenedAudioControls<'trajectory>, ControlStreamError> {
        let mut checkpoint = || cx.checkpoint().map_err(|_| ControlStreamError::Cancelled);
        checkpoint()?;
        if self.audio.is_empty() {
            return Err(ControlStreamError::NoPositiveDurationIntervals);
        }
        let requested_bins = self.audio.len().div_ceil(intervals_per_bin.get());
        let mut bins = Vec::new();
        bins.try_reserve_exact(requested_bins)
            .map_err(|_| ControlStreamError::Capacity {
                artifact: "coarsened audio bins",
                requested: requested_bins,
            })?;
        let mut cursor = 0;
        while cursor < self.audio.len() {
            checkpoint()?;
            if !self.audio[cursor].visual_coverage.is_fully_bracketed() {
                let event_barrier = !self.audio[cursor].events.is_empty();
                push_coarsened_bin(
                    &mut bins,
                    coarsen_group(&self.audio[cursor..=cursor], event_barrier, &mut checkpoint)?,
                )?;
                cursor += 1;
                continue;
            }
            if !self.audio[cursor].events.is_empty() {
                push_coarsened_bin(
                    &mut bins,
                    coarsen_group(&self.audio[cursor..=cursor], true, &mut checkpoint)?,
                )?;
                cursor += 1;
                continue;
            }
            let start = cursor;
            while cursor < self.audio.len()
                && cursor - start < intervals_per_bin.get()
                && self.audio[cursor].events.is_empty()
            {
                checkpoint()?;
                cursor += 1;
            }
            push_coarsened_bin(
                &mut bins,
                coarsen_group(&self.audio[start..cursor], false, &mut checkpoint)?,
            )?;
        }
        checkpoint()?;
        let represented = coarsened_work_checks(
            &bins,
            self.source.metadata().channel_availability,
            &mut checkpoint,
        )?;
        let last_sample = bins.last().map_or(0, |bin| bin.last_source_sample_index);
        let reconciliation = reconcile_against_raw(self.reconciliation, represented, last_sample)?;
        Ok(CoarsenedAudioControls {
            source: self.source,
            filter: AudioControlFilter::WholeIntervalBoxcarV1,
            intervals_per_bin,
            bins,
            reconciliation,
        })
    }
}

/// One whole-interval boxcar-coarsened output control bin.
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
    /// Visualization endpoint coverage inherited from the source intervals.
    pub visual_coverage: AudioVisualCoverage,
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

/// Result of deterministic whole-interval boxcar prefiltering and coarsening.
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

    /// Declared boxcar prefilter/coarsening rule.
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

    /// Coarsened bins whose opening and closing visualization states are both
    /// exactly retained.
    #[must_use]
    pub fn fully_bracketed_bins(&self) -> &[CoarsenedAudioBin] {
        let first = usize::from(
            self.bins
                .first()
                .is_some_and(|bin| !bin.visual_coverage.is_fully_bracketed()),
        );
        &self.bins[first..]
    }

    /// Common positive-duration visualization/audio clock range after
    /// coarsening, or `None` when every bin includes endpoint-only preroll.
    #[must_use]
    pub fn audio_visual_horizon(&self) -> Option<AudioVisualHorizon> {
        let synchronized = self.fully_bracketed_bins();
        Some(AudioVisualHorizon {
            start_time_s: synchronized.first()?.start_time_s,
            end_time_s: synchronized.last()?.end_time_s,
        })
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

/// Typed refusal from control derivation or boxcar coarsening.
#[derive(Clone, Debug, PartialEq)]
pub enum ControlStreamError {
    /// Cancellation was observed before atomic publication.
    Cancelled,
    /// The allocator refused an explicitly preflighted output capacity.
    Capacity {
        /// Output collection being constructed.
        artifact: &'static str,
        /// Requested element capacity.
        requested: usize,
    },
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
    channels: &ChannelOwnership,
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
    checkpoint: &mut impl FnMut() -> Result<(), ControlStreamError>,
) -> Result<CoarsenedAudioBin, ControlStreamError> {
    let first = &intervals[0];
    let last = &intervals[intervals.len() - 1];
    for pair in intervals.windows(2) {
        checkpoint()?;
        if pair[0].end_time_s.to_bits() != pair[1].start_time_s.to_bits() {
            return Err(ControlStreamError::NonContiguousIntervals {
                interval: pair[1].source_sample_index,
            });
        }
    }
    let duration_s = last.end_time_s - first.start_time_s;
    let channels = aggregate_channel_sets(intervals, duration_s, &mut *checkpoint)?;
    let mut weighted_normal_force = 0.0;
    let mut normal_available = true;
    let mut interval_contact_active = false;
    let mut events = Vec::new();
    for interval in intervals {
        checkpoint()?;
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
        visual_coverage: AudioVisualCoverage {
            start_visualization_index: first.visual_coverage.start_visualization_index,
            end_visualization_index: last.visual_coverage.end_visualization_index,
        },
        interval_contact_active,
        mean_base_normal_contact_force_n,
        channels,
        events,
        event_barrier,
    })
}

fn push_coarsened_bin(
    bins: &mut Vec<CoarsenedAudioBin>,
    bin: CoarsenedAudioBin,
) -> Result<(), ControlStreamError> {
    if bins.len() == bins.capacity() {
        bins.try_reserve(1)
            .map_err(|_| ControlStreamError::Capacity {
                artifact: "coarsened audio bins",
                requested: bins.len().saturating_add(1),
            })?;
    }
    bins.push(bin);
    Ok(())
}

fn aggregate_channel_sets(
    intervals: &[AudioControlInterval],
    duration_s: f64,
    checkpoint: &mut impl FnMut() -> Result<(), ControlStreamError>,
) -> Result<ChannelControlSet, ControlStreamError> {
    Ok(ChannelControlSet {
        gravity: aggregate_channel(intervals, duration_s, |set| set.gravity, &mut *checkpoint)?,
        contact: aggregate_channel(intervals, duration_s, |set| set.contact, &mut *checkpoint)?,
        rolling: aggregate_channel(intervals, duration_s, |set| set.rolling, &mut *checkpoint)?,
        base: aggregate_channel(intervals, duration_s, |set| set.base, &mut *checkpoint)?,
        gas: aggregate_channel(intervals, duration_s, |set| set.gas, &mut *checkpoint)?,
    })
}

fn aggregate_channel(
    intervals: &[AudioControlInterval],
    duration_s: f64,
    select: impl Fn(ChannelControlSet) -> ChannelControl,
    checkpoint: &mut impl FnMut() -> Result<(), ControlStreamError>,
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
        checkpoint()?;
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
    checkpoint: &mut impl FnMut() -> Result<(), ControlStreamError>,
) -> Result<ChannelWorkIntegralChecks, ControlStreamError> {
    let mut accumulators = ChannelWorkAccumulators::new(availability);
    for bin in bins {
        checkpoint()?;
        accumulators.accumulate(&bin.channels, bin.duration_s, bin.last_source_sample_index)?;
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
        channels: &ChannelControlSet,
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
        accumulate_channel_work(&mut self.base, channels.base, duration_s, sample, "base")?;
        accumulate_channel_work(&mut self.gas, channels.gas, duration_s, sample, "gas")
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
    let (accumulator, channel) = match (accumulator.as_mut(), channel) {
        (Some(accumulator), ChannelControl::Available(channel)) => (accumulator, channel),
        (None, ChannelControl::Unavailable) => return Ok(()),
        _ => {
            return Err(ControlStreamError::ChannelAvailabilityMismatch {
                sample,
                channel: channel_name,
            });
        }
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
    sample: usize,
) -> Result<ChannelWorkIntegralChecks, ControlStreamError> {
    Ok(ChannelWorkIntegralChecks {
        gravity: reconcile_channel(raw.gravity, represented.gravity, sample, "gravity")?,
        contact: reconcile_channel(raw.contact, represented.contact, sample, "contact")?,
        rolling: reconcile_channel(raw.rolling, represented.rolling, sample, "rolling")?,
        base: reconcile_channel(raw.base, represented.base, sample, "base")?,
        gas: reconcile_channel(raw.gas, represented.gas, sample, "gas")?,
    })
}

fn reconcile_channel(
    raw: Option<WorkIntegralCheck>,
    represented: Option<WorkIntegralCheck>,
    sample: usize,
    channel_name: &'static str,
) -> Result<Option<WorkIntegralCheck>, ControlStreamError> {
    match (raw, represented) {
        (Some(raw), Some(represented)) => {
            let residual_j = represented.integrated_work_j - raw.retained_work_j;
            if !residual_j.is_finite() {
                return Err(ControlStreamError::NonFiniteDerived {
                    sample,
                    field: "coarsened.work_integral_residual_j",
                });
            }
            Ok(Some(WorkIntegralCheck {
                retained_work_j: raw.retained_work_j,
                integrated_work_j: represented.integrated_work_j,
                residual_j,
            }))
        }
        (None, None) => Ok(None),
        _ => Err(ControlStreamError::ChannelAvailabilityMismatch {
            sample,
            channel: channel_name,
        }),
    }
}

fn finite_vec(value: Vec3, sample: usize, field: &'static str) -> Result<(), ControlStreamError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(ControlStreamError::NonFiniteDerived { sample, field })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn available_channel() -> ChannelControl {
        ChannelControl::Available(AvailableChannelControl {
            mean_force_world_n: Vec3::new(0.0, 0.0, 1.0),
            mean_torque_world_nm: Vec3::new(0.0, 0.0, -0.25),
            signed_work_j: -0.5,
            signed_mean_work_rate_w: -0.5,
            force_time_measure_world_n_s: Vec3::new(0.0, 0.0, 1.0),
            torque_time_measure_world_nm_s: Vec3::new(0.0, 0.0, -0.25),
        })
    }

    fn interval(index: usize) -> AudioControlInterval {
        let start_time_s = index as f64;
        let channel = available_channel();
        AudioControlInterval {
            source_sample_index: index + 1,
            start_time_s,
            end_time_s: start_time_s + 1.0,
            duration_s: 1.0,
            visual_coverage: AudioVisualCoverage {
                start_visualization_index: Some(index),
                end_visualization_index: index + 1,
            },
            interval_contact_active: true,
            mean_base_normal_contact_force_n: Some(1.0),
            normal_force_sampling: RenderNormalForceSampling::IntervalMean,
            declared_normal_force_n: 1.0,
            channels: ChannelControlSet {
                gravity: channel,
                contact: channel,
                rolling: channel,
                base: channel,
                gas: channel,
            },
            endpoint_center_of_mass_velocity_world_m_per_s: Vec3::ZERO,
            endpoint_angular_velocity_world_rad_per_s: Vec3::ZERO,
            endpoint_base_velocity_world_m_per_s: Vec3::ZERO,
            events: Vec::new(),
        }
    }

    #[test]
    fn g4_inner_coarsening_loops_observe_deterministic_cancellation() {
        let intervals = [interval(0), interval(1), interval(2)];

        // The selected checkpoints land respectively in continuity scanning,
        // channel aggregation, and the post-aggregation interval fold.
        for cancel_at in [1, 3, 18] {
            let mut observed = 0;
            let mut checkpoint = || {
                observed += 1;
                if observed == cancel_at {
                    Err(ControlStreamError::Cancelled)
                } else {
                    Ok(())
                }
            };
            assert_eq!(
                coarsen_group(&intervals, false, &mut checkpoint).unwrap_err(),
                ControlStreamError::Cancelled
            );
            assert_eq!(observed, cancel_at);
        }

        let bin = coarsen_group(&intervals, false, &mut || Ok(())).unwrap();
        let mut observed = 0;
        let mut checkpoint = || {
            observed += 1;
            Err(ControlStreamError::Cancelled)
        };
        assert_eq!(
            coarsened_work_checks(
                &[bin],
                RenderChannelAvailability::ALL_AVAILABLE,
                &mut checkpoint,
            )
            .unwrap_err(),
            ControlStreamError::Cancelled
        );
        assert_eq!(observed, 1);
    }
}
