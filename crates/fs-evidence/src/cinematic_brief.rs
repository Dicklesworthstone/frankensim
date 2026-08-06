//! Typed creative brief and exact A/V timeline for the reference Euler film.
//!
//! This module is acceptance input shared by rendering, sound, CLI, and report
//! layers. It does not implement a camera or modify simulation state.

use core::f64::consts::TAU;
use core::fmt;

use fs_blake3::{ContentHash, hash_domain};

use crate::cinematic::CinematicDeliverableContract;

/// Version of the canonical creative-brief identity preimage.
pub const CINEMATIC_BRIEF_IDENTITY_VERSION: u16 = 1;
/// Domain separating creative briefs from every referenced trajectory,
/// renderer, image, and audio artifact.
pub const CINEMATIC_BRIEF_IDENTITY_DOMAIN: &str = "org.frankensim.cinematic-brief.identity.v1";

const CINEMATIC_BRIEF_MAGIC: &[u8; 8] = b"FSCBRF01";

/// One half-open range of video frames.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FrameRange {
    start: u32,
    end_exclusive: u32,
}

impl FrameRange {
    /// Construct a nonempty half-open range.
    pub fn try_new(start: u32, end_exclusive: u32) -> Result<Self, CinematicBriefError> {
        if end_exclusive <= start {
            return Err(CinematicBriefError::EmptyFrameRange);
        }
        Ok(Self {
            start,
            end_exclusive,
        })
    }

    /// Inclusive first frame.
    #[must_use]
    pub const fn start(self) -> u32 {
        self.start
    }

    /// Exclusive final frame.
    #[must_use]
    pub const fn end_exclusive(self) -> u32 {
        self.end_exclusive
    }

    /// Number of frames.
    #[must_use]
    pub const fn len(self) -> u32 {
        self.end_exclusive - self.start
    }

    /// Whether the range is empty. Admitted ranges always return false.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.start == self.end_exclusive
    }

    const fn contains(self, frame: u32) -> bool {
        frame >= self.start && frame < self.end_exclusive
    }
}

/// A point or direction in the brief's declared studio frame, in metres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BriefVec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl BriefVec3 {
    fn finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    fn norm_squared(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
}

/// Camera keyframe; interpolation occurs in the declared studio frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BriefCameraKeyframe {
    /// Video frame carrying the keyframe.
    pub frame: u32,
    /// Camera optical center in metres.
    pub eye_m: BriefVec3,
    /// Look-at target in metres.
    pub target_m: BriefVec3,
    /// Up reference; need not be normalized but must not be singular.
    pub up: BriefVec3,
    /// Focus distance in metres.
    pub focus_distance_m: f64,
    /// Semantic focus target used to diagnose missing/incorrect focus pulls.
    pub focus_target: FocusTarget,
}

/// Semantic focus target.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FocusTarget {
    /// Disc center of mass/visual center.
    DiscCenter,
    /// Disc/base contact neighborhood.
    ContactPoint,
    /// Readable outer rim/precession silhouette.
    DiscRim,
}

/// Stable narrative role of a shot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReferenceShotRole {
    /// Studio product establishing view.
    EstablishingProduct,
    /// Inclination and contact-orbit explanation view.
    InclinationAndContactOrbit,
    /// Macro precession and wobble view.
    MacroPrecession,
    /// Terminal close-up.
    TerminalCloseUp,
}

/// Lens and exposure intent; the renderer remains free to implement it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BriefOptics {
    /// Focal length in micrometres (50 mm = 50,000 µm).
    pub focal_length_um: u32,
    /// F-number multiplied by 1,000.
    pub f_number_milli: u32,
    /// Shutter-open offset from frame center in millionths of one frame.
    pub shutter_open_microframes: i32,
    /// Shutter-close offset from frame center in millionths of one frame.
    pub shutter_close_microframes: i32,
    /// Exposure intent.
    pub exposure: ExposureIntent,
}

/// Exposure intent, separate from a renderer-specific EV implementation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ExposureIntent {
    /// Preserve specular detail on the metal disc.
    ProtectMetalHighlights,
    /// Preserve glass caustic/highlight detail in a close view.
    ProtectGlassHighlights,
}

/// Reference lighting preset.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LightingPreset {
    /// Large soft key, controlled rim, dark studio environment.
    StudioSoftboxRim,
    /// Raking macro light that reveals contact and brushing.
    RakingContactMacro,
    /// Controlled terminal glint without clipped highlights.
    TerminalGlint,
}

/// Reference visible disc material.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DiscMaterialPreset {
    /// Brushed tungsten.
    BrushedTungsten,
    /// Brushed stainless steel.
    BrushedStainlessSteel,
}

/// Background intent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BackgroundPreset {
    /// Near-black neutral studio sweep without a false horizon.
    NeutralBlackSweep,
}

/// Audio listening perspective for a shot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AudioPerspective {
    /// Human-scale studio listening position.
    StudioObserver,
    /// Close material/contact microphone perspective.
    ContactMacro,
    /// Close terminal perspective with controlled gain taper.
    TerminalDetail,
}

/// How a cut is formed. The v1 reference uses only hard cuts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ShotTransition {
    /// No blended frames; shutter support is clipped at the cut.
    HardCut,
}

/// Mapping from presentation frames to simulation ticks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ShotTimeMapping {
    /// Physically linear use of accepted simulation state.
    PhysicalLinear {
        /// Inclusive first simulation tick.
        start_tick: u64,
        /// Exclusive terminal simulation tick.
        end_tick_exclusive: u64,
    },
    /// Explicit visualization-only hold for censored input.
    VisualizationHold {
        /// Last accepted source tick held on screen.
        source_tick: u64,
        /// Stable reason label; a hold can never be visually silent.
        label: VisualizationHoldLabel,
    },
}

/// Human/machine label for a visualization-only hold.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VisualizationHoldLabel {
    /// The accepted trajectory ended before the requested presentation.
    CensoredTrajectory,
}

/// Required apparent-motion/readability result for a shot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ApparentRotationIntent {
    /// Product view must show forward rotation without obscuring inclination.
    ForwardCueAndInclinationReadable,
    /// Contact orbit is primary; false spin reversal/flicker are forbidden.
    ContactOrbitPrimaryNoFalseReversal,
    /// Wobble/precession is primary; cue must not freeze or dominate it.
    WobblePrimaryNoFrozenCue,
    /// Terminal motion may slow but may not acquire alias-driven reversal.
    TerminalSlowdownNoAliasReversal,
}

/// One editable shot specification.
#[derive(Debug, Clone, PartialEq)]
pub struct ReferenceShotInput {
    /// Narrative role.
    pub role: ReferenceShotRole,
    /// Contiguous frame range.
    pub frames: FrameRange,
    /// Camera keyframes within `frames`.
    pub camera: Vec<BriefCameraKeyframe>,
    /// Lens, aperture, shutter, and exposure intent.
    pub optics: BriefOptics,
    /// Lighting intent.
    pub lighting: LightingPreset,
    /// Disc appearance.
    pub disc_material: DiscMaterialPreset,
    /// Background intent.
    pub background: BackgroundPreset,
    /// Whether the glass plate must remain visible.
    pub glass_visible: bool,
    /// Whether the base/housing must remain visible.
    pub base_visible: bool,
    /// Audio perspective.
    pub audio_perspective: AudioPerspective,
    /// Transition entering this shot.
    pub transition: ShotTransition,
    /// Physical or explicitly visualization-only time mapping.
    pub time_mapping: ShotTimeMapping,
    /// Shot-specific apparent-motion acceptance intent.
    pub apparent_rotation_intent: ApparentRotationIntent,
}

/// Opaque admitted shot.
#[derive(Debug, Clone, PartialEq)]
pub struct ReferenceShot(ReferenceShotInput);

impl ReferenceShot {
    /// Narrative role.
    #[must_use]
    pub const fn role(&self) -> ReferenceShotRole {
        self.0.role
    }

    /// Frame range.
    #[must_use]
    pub const fn frames(&self) -> FrameRange {
        self.0.frames
    }

    /// Camera keyframes.
    #[must_use]
    pub fn camera(&self) -> &[BriefCameraKeyframe] {
        &self.0.camera
    }

    /// Optics/exposure intent.
    #[must_use]
    pub const fn optics(&self) -> BriefOptics {
        self.0.optics
    }

    /// Lighting intent.
    #[must_use]
    pub const fn lighting(&self) -> LightingPreset {
        self.0.lighting
    }

    /// Disc appearance preset.
    #[must_use]
    pub const fn disc_material(&self) -> DiscMaterialPreset {
        self.0.disc_material
    }

    /// Background intent.
    #[must_use]
    pub const fn background(&self) -> BackgroundPreset {
        self.0.background
    }

    /// Glass visibility intent.
    #[must_use]
    pub const fn glass_visible(&self) -> bool {
        self.0.glass_visible
    }

    /// Base/housing visibility intent.
    #[must_use]
    pub const fn base_visible(&self) -> bool {
        self.0.base_visible
    }

    /// Audio listening perspective.
    #[must_use]
    pub const fn audio_perspective(&self) -> AudioPerspective {
        self.0.audio_perspective
    }

    /// Incoming transition.
    #[must_use]
    pub const fn transition(&self) -> ShotTransition {
        self.0.transition
    }

    /// Time mapping.
    #[must_use]
    pub const fn time_mapping(&self) -> ShotTimeMapping {
        self.0.time_mapping
    }

    /// Required apparent-motion/readability behavior.
    #[must_use]
    pub const fn apparent_rotation_intent(&self) -> ApparentRotationIntent {
        self.0.apparent_rotation_intent
    }
}

/// Safe-area insets in thousandths of frame width/height.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SafeAreaInsetsPermille {
    /// Action-safe inset.
    pub action: u16,
    /// Title-safe inset.
    pub title: u16,
}

/// Visualization-only spin cue. It has deliberately no mass/contact fields.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SpinCue {
    /// Rotational repetition count of subtle anisotropic brushing.
    pub brushing_marking_frequency: u16,
    /// Whether a tiny engraved radial mark is enabled.
    pub engraved_radial_mark: bool,
}

/// Handling of a trajectory that ends before the requested presentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CensoredTrajectoryPolicy {
    /// Refuse an unlabeled continuation. A preview may hold the last accepted
    /// state for at most this many frames with an explicit censored label and
    /// audio taper.
    LabeledHoldAndAudioTaper {
        /// Maximum visualization-only hold.
        maximum_hold_frames: u32,
    },
}

/// How presentation-frame instants are placed on the master clock.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FrameTimeConvention {
    /// Frame `k` is centered at `k / 24` seconds; the master endpoint is the
    /// exclusive `total_frames / 24` boundary.
    IntegerFrameCentersHalfOpenMaster,
}

/// Camera and scalar interpolation semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BriefInterpolationDomain {
    /// Rigid camera pose in the studio frame plus independently interpolated
    /// focus/optics scalars; no interpolation modifies simulation state.
    StudioFrameRigidPoseAndScalarOptics,
}

/// How shutter support behaves at cuts and master boundaries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ShutterBoundaryPolicy {
    /// Clip support to the current shot and master; never integrate two shots.
    ClipAtCutsAndMaster,
}

/// Complete caller input for an editable brief.
#[derive(Debug, Clone, PartialEq)]
pub struct CinematicBriefInput {
    /// Total video frames.
    pub total_frames: u32,
    /// Audio sample frames; stereo samples at one instant count as one frame.
    pub total_audio_sample_frames: u64,
    /// Simulation clock rate in ticks per second.
    pub simulation_ticks_per_second: u32,
    /// Inclusive available simulation tick.
    pub trajectory_start_tick: u64,
    /// Exclusive available simulation tick.
    pub trajectory_end_tick_exclusive: u64,
    /// Ordered shot list.
    pub shots: Vec<ReferenceShotInput>,
    /// Safe action/title regions.
    pub safe_areas: SafeAreaInsetsPermille,
    /// Visualization-only rotation cue.
    pub spin_cue: SpinCue,
    /// Censor/early-termination display policy.
    pub censored_policy: CensoredTrajectoryPolicy,
    /// Frame-center and first/last convention.
    pub frame_time_convention: FrameTimeConvention,
    /// Camera/keyframe interpolation domain.
    pub interpolation_domain: BriefInterpolationDomain,
    /// Cut/master shutter behavior.
    pub shutter_boundary_policy: ShutterBoundaryPolicy,
    /// Signed audio lead (positive) or lag (negative), in 48 kHz samples.
    pub audio_lead_samples: i64,
    /// Muted-audio review is mandatory.
    pub muted_review_required: bool,
    /// Essential scientific context may not exist only in sound.
    pub essential_context_in_audio_only: bool,
}

/// Opaque, validated reference or user-edited cinematic brief.
#[derive(Debug, Clone, PartialEq)]
pub struct CinematicBrief {
    total_frames: u32,
    total_audio_sample_frames: u64,
    simulation_ticks_per_second: u32,
    trajectory_start_tick: u64,
    trajectory_end_tick_exclusive: u64,
    shots: Vec<ReferenceShot>,
    safe_areas: SafeAreaInsetsPermille,
    spin_cue: SpinCue,
    censored_policy: CensoredTrajectoryPolicy,
    frame_time_convention: FrameTimeConvention,
    interpolation_domain: BriefInterpolationDomain,
    shutter_boundary_policy: ShutterBoundaryPolicy,
    audio_lead_samples: i64,
}

impl CinematicBrief {
    /// Validate an editable brief against the frozen v1 delivery envelope.
    pub fn try_new(input: CinematicBriefInput) -> Result<Self, CinematicBriefError> {
        let deliverable = CinematicDeliverableContract::euler_disc_v1();
        deliverable
            .validate_timeline(input.total_frames, input.total_audio_sample_frames)
            .map_err(|_| CinematicBriefError::InvalidMasterTimeline)?;
        if input.simulation_ticks_per_second == 0
            || input.trajectory_end_tick_exclusive <= input.trajectory_start_tick
        {
            return Err(CinematicBriefError::InvalidTrajectoryClock);
        }
        if input.shots.is_empty() {
            return Err(CinematicBriefError::MissingShots);
        }
        if input.shots.len() > input.total_frames as usize {
            return Err(CinematicBriefError::TooManyShots {
                maximum: input.total_frames,
                got: u64::try_from(input.shots.len()).unwrap_or(u64::MAX),
            });
        }
        if input.safe_areas.action > 500
            || input.safe_areas.title > 500
            || input.safe_areas.title < input.safe_areas.action
        {
            return Err(CinematicBriefError::InvalidSafeAreas);
        }
        if input.spin_cue.brushing_marking_frequency == 0 {
            return Err(CinematicBriefError::InvalidSpinCue);
        }
        let maximum_hold_frames = match input.censored_policy {
            CensoredTrajectoryPolicy::LabeledHoldAndAudioTaper {
                maximum_hold_frames,
            } if maximum_hold_frames > 0 => maximum_hold_frames,
            CensoredTrajectoryPolicy::LabeledHoldAndAudioTaper { .. } => {
                return Err(CinematicBriefError::InvalidCensoredPolicy);
            }
        };
        // The reference master deliberately freezes sync at zero offset. A
        // future non-zero policy must first define boundary padding/truncation
        // rather than silently shifting samples outside the master.
        if input.audio_lead_samples != 0 {
            return Err(CinematicBriefError::UnsupportedAudioLeadLag);
        }
        if !input.muted_review_required || input.essential_context_in_audio_only {
            return Err(CinematicBriefError::AudioOnlyMeaning);
        }

        let mut expected_start = 0;
        let mut held_frames = 0_u32;
        let mut saw_hold = false;
        let mut shots = Vec::with_capacity(input.shots.len());
        for shot in input.shots {
            if shot.camera.len() > shot.frames.len() as usize {
                return Err(CinematicBriefError::TooManyCameraKeyframes {
                    maximum: shot.frames.len(),
                    got: u64::try_from(shot.camera.len()).unwrap_or(u64::MAX),
                });
            }
            if shot.frames.start != expected_start {
                return Err(CinematicBriefError::ShotGapOrOverlap {
                    expected_start,
                    got: shot.frames.start,
                });
            }
            if shot.frames.end_exclusive > input.total_frames {
                return Err(CinematicBriefError::ShotOutsideMaster);
            }
            validate_shot(
                &shot,
                input.trajectory_start_tick,
                input.trajectory_end_tick_exclusive,
            )?;
            match shot.time_mapping {
                ShotTimeMapping::PhysicalLinear { .. } if saw_hold => {
                    return Err(CinematicBriefError::PhysicalShotAfterCensoredHold);
                }
                ShotTimeMapping::PhysicalLinear { .. } => {}
                ShotTimeMapping::VisualizationHold { source_tick, .. } => {
                    saw_hold = true;
                    if source_tick + 1 != input.trajectory_end_tick_exclusive {
                        return Err(CinematicBriefError::HoldIsNotTerminalState);
                    }
                    held_frames = held_frames.saturating_add(shot.frames.len());
                    if held_frames > maximum_hold_frames {
                        return Err(CinematicBriefError::CensoredHoldTooLong {
                            maximum: maximum_hold_frames,
                            got: held_frames,
                        });
                    }
                }
            }
            expected_start = shot.frames.end_exclusive;
            shots.push(ReferenceShot(shot));
        }
        if expected_start != input.total_frames {
            return Err(CinematicBriefError::ShotGapOrOverlap {
                expected_start,
                got: input.total_frames,
            });
        }

        Ok(Self {
            total_frames: input.total_frames,
            total_audio_sample_frames: input.total_audio_sample_frames,
            simulation_ticks_per_second: input.simulation_ticks_per_second,
            trajectory_start_tick: input.trajectory_start_tick,
            trajectory_end_tick_exclusive: input.trajectory_end_tick_exclusive,
            shots,
            safe_areas: input.safe_areas,
            spin_cue: input.spin_cue,
            censored_policy: input.censored_policy,
            frame_time_convention: input.frame_time_convention,
            interpolation_domain: input.interpolation_domain,
            shutter_boundary_policy: input.shutter_boundary_policy,
            audio_lead_samples: input.audio_lead_samples,
        })
    }

    /// Frozen ten-second, four-shot Euler reference brief.
    pub fn euler_disc_v1() -> Result<Self, CinematicBriefError> {
        Self::try_new(reference_input())
    }

    /// Canonical, locator-free encoding of every admitted brief semantic.
    ///
    /// This is an identity preimage rather than a general interchange codec.
    /// It lets independent render/audio finalizers bind the exact shot, clock,
    /// cut, camera, material, and no-audio-only-meaning contract they consumed.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(CINEMATIC_BRIEF_MAGIC);
        bytes.extend_from_slice(&CINEMATIC_BRIEF_IDENTITY_VERSION.to_le_bytes());
        bytes.extend_from_slice(&self.total_frames.to_le_bytes());
        bytes.extend_from_slice(&self.total_audio_sample_frames.to_le_bytes());
        bytes.extend_from_slice(&self.simulation_ticks_per_second.to_le_bytes());
        bytes.extend_from_slice(&self.trajectory_start_tick.to_le_bytes());
        bytes.extend_from_slice(&self.trajectory_end_tick_exclusive.to_le_bytes());
        bytes.extend_from_slice(
            &u32::try_from(self.shots.len())
                .expect("admission bounds shots by the u32 frame count")
                .to_le_bytes(),
        );
        for shot in &self.shots {
            push_shot(&mut bytes, shot);
        }
        bytes.extend_from_slice(&self.safe_areas.action.to_le_bytes());
        bytes.extend_from_slice(&self.safe_areas.title.to_le_bytes());
        bytes.extend_from_slice(&self.spin_cue.brushing_marking_frequency.to_le_bytes());
        bytes.push(u8::from(self.spin_cue.engraved_radial_mark));
        match self.censored_policy {
            CensoredTrajectoryPolicy::LabeledHoldAndAudioTaper {
                maximum_hold_frames,
            } => {
                bytes.push(1);
                bytes.extend_from_slice(&maximum_hold_frames.to_le_bytes());
            }
        }
        bytes.push(frame_time_convention_tag(self.frame_time_convention));
        bytes.push(interpolation_domain_tag(self.interpolation_domain));
        bytes.push(shutter_boundary_policy_tag(self.shutter_boundary_policy));
        bytes.extend_from_slice(&self.audio_lead_samples.to_le_bytes());
        // These values are invariants of every admitted brief even though the
        // struct need not retain redundant booleans after admission.
        bytes.push(1); // muted_review_required
        bytes.push(0); // essential_context_in_audio_only
        bytes
    }

    /// Domain-separated identity of [`Self::canonical_bytes`].
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        hash_domain(CINEMATIC_BRIEF_IDENTITY_DOMAIN, &self.canonical_bytes())
    }

    /// Total video frames.
    #[must_use]
    pub const fn total_frames(&self) -> u32 {
        self.total_frames
    }

    /// Total audio sample frames.
    #[must_use]
    pub const fn total_audio_sample_frames(&self) -> u64 {
        self.total_audio_sample_frames
    }

    /// Exact simulation clock rate.
    #[must_use]
    pub const fn simulation_ticks_per_second(&self) -> u32 {
        self.simulation_ticks_per_second
    }

    /// Ordered shots.
    #[must_use]
    pub fn shots(&self) -> &[ReferenceShot] {
        &self.shots
    }

    /// Safe regions.
    #[must_use]
    pub const fn safe_areas(&self) -> SafeAreaInsetsPermille {
        self.safe_areas
    }

    /// Visualization-only spin cue.
    #[must_use]
    pub const fn spin_cue(&self) -> SpinCue {
        self.spin_cue
    }

    /// Censored input policy.
    #[must_use]
    pub const fn censored_policy(&self) -> CensoredTrajectoryPolicy {
        self.censored_policy
    }

    /// Frame-center and first/last convention.
    #[must_use]
    pub const fn frame_time_convention(&self) -> FrameTimeConvention {
        self.frame_time_convention
    }

    /// Interpolation domain.
    #[must_use]
    pub const fn interpolation_domain(&self) -> BriefInterpolationDomain {
        self.interpolation_domain
    }

    /// Shutter support at cuts/master boundaries.
    #[must_use]
    pub const fn shutter_boundary_policy(&self) -> ShutterBoundaryPolicy {
        self.shutter_boundary_policy
    }

    /// Signed audio lead/lag in sample frames.
    #[must_use]
    pub const fn audio_lead_samples(&self) -> i64 {
        self.audio_lead_samples
    }

    /// Exact start sample for a video-frame boundary.
    pub fn audio_sample_for_frame(&self, frame: u32) -> Result<u64, CinematicBriefError> {
        if frame > self.total_frames {
            return Err(CinematicBriefError::FrameOutsideMaster(frame));
        }
        let numerator = u128::from(frame) * u128::from(self.total_audio_sample_frames);
        let denominator = u128::from(self.total_frames);
        debug_assert_eq!(numerator % denominator, 0);
        Ok(u64::try_from(numerator / denominator).unwrap_or(self.total_audio_sample_frames))
    }

    /// Exact rational simulation tick for one audio sample frame.
    ///
    /// The audio sample frame is mapped continuously within the same shot-time
    /// interval as video. This avoids quantizing sound synthesis to video
    /// frames while keeping cuts sample-exact.
    pub fn simulation_tick_for_audio_sample(
        &self,
        sample: u64,
    ) -> Result<RationalSimulationTick, CinematicBriefError> {
        if sample >= self.total_audio_sample_frames {
            return Err(CinematicBriefError::AudioSampleOutsideMaster(sample));
        }
        let samples_per_frame = self.total_audio_sample_frames / u64::from(self.total_frames);
        let frame = u32::try_from(sample / samples_per_frame)
            .map_err(|_| CinematicBriefError::AudioSampleOutsideMaster(sample))?;
        let shot = self.shot_for_frame(frame)?;
        match shot.time_mapping() {
            ShotTimeMapping::PhysicalLinear {
                start_tick,
                end_tick_exclusive,
            } => {
                let shot_start_sample = u64::from(shot.frames().start()) * samples_per_frame;
                let shot_sample_count = u64::from(shot.frames().len()) * samples_per_frame;
                let local_sample = sample - shot_start_sample;
                Ok(RationalSimulationTick {
                    numerator: u128::from(start_tick) * u128::from(shot_sample_count)
                        + u128::from(end_tick_exclusive - start_tick) * u128::from(local_sample),
                    denominator: shot_sample_count,
                    visualization_only: false,
                })
            }
            ShotTimeMapping::VisualizationHold { source_tick, .. } => Ok(RationalSimulationTick {
                numerator: u128::from(source_tick),
                denominator: 1,
                visualization_only: true,
            }),
        }
    }

    /// Shutter support for one frame after the declared cut/master clipping.
    pub fn effective_shutter_window(
        &self,
        frame: u32,
    ) -> Result<EffectiveShutterWindow, CinematicBriefError> {
        const MICROFRAMES_PER_FRAME: i64 = 1_000_000;
        let shot = self.shot_for_frame(frame)?;
        let center = i64::from(frame) * MICROFRAMES_PER_FRAME;
        let requested_start = center + i64::from(shot.optics().shutter_open_microframes);
        let requested_end = center + i64::from(shot.optics().shutter_close_microframes);
        let shot_start = i64::from(shot.frames().start()) * MICROFRAMES_PER_FRAME;
        let shot_end = i64::from(shot.frames().end_exclusive()) * MICROFRAMES_PER_FRAME;
        let start = requested_start.max(shot_start).max(0);
        let end = requested_end
            .min(shot_end)
            .min(i64::from(self.total_frames) * MICROFRAMES_PER_FRAME);
        Ok(EffectiveShutterWindow {
            start_microframes: start,
            end_microframes: end,
            clipped_at_boundary: start != requested_start || end != requested_end,
        })
    }

    /// Shot containing a presentation frame.
    pub fn shot_for_frame(&self, frame: u32) -> Result<&ReferenceShot, CinematicBriefError> {
        self.shots
            .iter()
            .find(|shot| shot.frames().contains(frame))
            .ok_or(CinematicBriefError::FrameOutsideMaster(frame))
    }

    /// Exact rational simulation tick for the center of a physical frame.
    pub fn simulation_tick_for_frame(
        &self,
        frame: u32,
    ) -> Result<RationalSimulationTick, CinematicBriefError> {
        let shot = self.shot_for_frame(frame)?;
        match shot.time_mapping() {
            ShotTimeMapping::PhysicalLinear {
                start_tick,
                end_tick_exclusive,
            } => {
                let local = u64::from(frame - shot.frames().start());
                let length = u64::from(shot.frames().len());
                Ok(RationalSimulationTick {
                    numerator: u128::from(start_tick) * u128::from(length)
                        + u128::from(end_tick_exclusive - start_tick) * u128::from(local),
                    denominator: length,
                    visualization_only: false,
                })
            }
            ShotTimeMapping::VisualizationHold {
                source_tick,
                label: _,
            } => Ok(RationalSimulationTick {
                numerator: u128::from(source_tick),
                denominator: 1,
                visualization_only: true,
            }),
        }
    }

    /// Generate the complete low-cost, render-independent storyboard proxy.
    pub fn storyboard_proxy(&self) -> Result<Vec<StoryboardFrame>, CinematicBriefError> {
        let mut frames = Vec::new();
        for frame in 0..self.total_frames {
            let shot = self.shot_for_frame(frame)?;
            frames.push(StoryboardFrame {
                frame,
                audio_sample: self.audio_sample_for_frame(frame)?,
                shot_role: shot.role(),
                simulation_tick: self.simulation_tick_for_frame(frame)?,
                audio_muted_for_review: true,
            });
        }
        Ok(frames)
    }

    /// Bounded review metadata for the muted proxy artifact.
    #[must_use]
    pub fn muted_review_manifest_json(&self) -> String {
        format!(
            "{{\"schema\":\"euler-cinematic-muted-storyboard-v1\",\"total_frames\":{},\"audio_muted\":true,\"essential_context_in_audio_only\":false,\"shot_count\":{}}}",
            self.total_frames,
            self.shots.len(),
        )
    }
}

fn push_shot(bytes: &mut Vec<u8>, shot: &ReferenceShot) {
    bytes.push(reference_shot_role_tag(shot.role()));
    bytes.extend_from_slice(&shot.frames().start().to_le_bytes());
    bytes.extend_from_slice(&shot.frames().end_exclusive().to_le_bytes());
    bytes.extend_from_slice(
        &u32::try_from(shot.camera().len())
            .expect("admission bounds keyframes by the u32 shot length")
            .to_le_bytes(),
    );
    for keyframe in shot.camera() {
        bytes.extend_from_slice(&keyframe.frame.to_le_bytes());
        for vector in [keyframe.eye_m, keyframe.target_m, keyframe.up] {
            bytes.extend_from_slice(&canonical_f64_bits(vector.x).to_le_bytes());
            bytes.extend_from_slice(&canonical_f64_bits(vector.y).to_le_bytes());
            bytes.extend_from_slice(&canonical_f64_bits(vector.z).to_le_bytes());
        }
        bytes.extend_from_slice(&canonical_f64_bits(keyframe.focus_distance_m).to_le_bytes());
        bytes.push(focus_target_tag(keyframe.focus_target));
    }
    let optics = shot.optics();
    bytes.extend_from_slice(&optics.focal_length_um.to_le_bytes());
    bytes.extend_from_slice(&optics.f_number_milli.to_le_bytes());
    bytes.extend_from_slice(&optics.shutter_open_microframes.to_le_bytes());
    bytes.extend_from_slice(&optics.shutter_close_microframes.to_le_bytes());
    bytes.push(exposure_intent_tag(optics.exposure));
    bytes.push(lighting_preset_tag(shot.lighting()));
    bytes.push(disc_material_tag(shot.disc_material()));
    bytes.push(background_preset_tag(shot.background()));
    bytes.push(u8::from(shot.glass_visible()));
    bytes.push(u8::from(shot.base_visible()));
    bytes.push(audio_perspective_tag(shot.audio_perspective()));
    bytes.push(shot_transition_tag(shot.transition()));
    match shot.time_mapping() {
        ShotTimeMapping::PhysicalLinear {
            start_tick,
            end_tick_exclusive,
        } => {
            bytes.push(1);
            bytes.extend_from_slice(&start_tick.to_le_bytes());
            bytes.extend_from_slice(&end_tick_exclusive.to_le_bytes());
        }
        ShotTimeMapping::VisualizationHold { source_tick, label } => {
            bytes.push(2);
            bytes.extend_from_slice(&source_tick.to_le_bytes());
            bytes.push(visualization_hold_label_tag(label));
        }
    }
    bytes.push(apparent_rotation_intent_tag(
        shot.apparent_rotation_intent(),
    ));
}

const fn canonical_f64_bits(value: f64) -> u64 {
    if value == 0.0 { 0 } else { value.to_bits() }
}

const fn focus_target_tag(value: FocusTarget) -> u8 {
    match value {
        FocusTarget::DiscCenter => 1,
        FocusTarget::ContactPoint => 2,
        FocusTarget::DiscRim => 3,
    }
}

const fn reference_shot_role_tag(value: ReferenceShotRole) -> u8 {
    match value {
        ReferenceShotRole::EstablishingProduct => 1,
        ReferenceShotRole::InclinationAndContactOrbit => 2,
        ReferenceShotRole::MacroPrecession => 3,
        ReferenceShotRole::TerminalCloseUp => 4,
    }
}

const fn exposure_intent_tag(value: ExposureIntent) -> u8 {
    match value {
        ExposureIntent::ProtectMetalHighlights => 1,
        ExposureIntent::ProtectGlassHighlights => 2,
    }
}

const fn lighting_preset_tag(value: LightingPreset) -> u8 {
    match value {
        LightingPreset::StudioSoftboxRim => 1,
        LightingPreset::RakingContactMacro => 2,
        LightingPreset::TerminalGlint => 3,
    }
}

const fn disc_material_tag(value: DiscMaterialPreset) -> u8 {
    match value {
        DiscMaterialPreset::BrushedTungsten => 1,
        DiscMaterialPreset::BrushedStainlessSteel => 2,
    }
}

const fn background_preset_tag(value: BackgroundPreset) -> u8 {
    match value {
        BackgroundPreset::NeutralBlackSweep => 1,
    }
}

const fn audio_perspective_tag(value: AudioPerspective) -> u8 {
    match value {
        AudioPerspective::StudioObserver => 1,
        AudioPerspective::ContactMacro => 2,
        AudioPerspective::TerminalDetail => 3,
    }
}

const fn shot_transition_tag(value: ShotTransition) -> u8 {
    match value {
        ShotTransition::HardCut => 1,
    }
}

const fn visualization_hold_label_tag(value: VisualizationHoldLabel) -> u8 {
    match value {
        VisualizationHoldLabel::CensoredTrajectory => 1,
    }
}

const fn apparent_rotation_intent_tag(value: ApparentRotationIntent) -> u8 {
    match value {
        ApparentRotationIntent::ForwardCueAndInclinationReadable => 1,
        ApparentRotationIntent::ContactOrbitPrimaryNoFalseReversal => 2,
        ApparentRotationIntent::WobblePrimaryNoFrozenCue => 3,
        ApparentRotationIntent::TerminalSlowdownNoAliasReversal => 4,
    }
}

const fn frame_time_convention_tag(value: FrameTimeConvention) -> u8 {
    match value {
        FrameTimeConvention::IntegerFrameCentersHalfOpenMaster => 1,
    }
}

const fn interpolation_domain_tag(value: BriefInterpolationDomain) -> u8 {
    match value {
        BriefInterpolationDomain::StudioFrameRigidPoseAndScalarOptics => 1,
    }
}

const fn shutter_boundary_policy_tag(value: ShutterBoundaryPolicy) -> u8 {
    match value {
        ShutterBoundaryPolicy::ClipAtCutsAndMaster => 1,
    }
}

/// One render-independent storyboard frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StoryboardFrame {
    /// Video frame number.
    pub frame: u32,
    /// Exact aligned audio sample frame.
    pub audio_sample: u64,
    /// Narrative shot role.
    pub shot_role: ReferenceShotRole,
    /// Exact source/hold tick.
    pub simulation_tick: RationalSimulationTick,
    /// Always true for the required muted review artifact.
    pub audio_muted_for_review: bool,
}

/// Exact simulation-tick mapping for one frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RationalSimulationTick {
    /// Numerator in simulation ticks.
    pub numerator: u128,
    /// Positive denominator.
    pub denominator: u64,
    /// True only for an explicitly labeled presentation hold.
    pub visualization_only: bool,
}

/// Effective half-open shutter support on the integer-frame master clock.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EffectiveShutterWindow {
    /// Inclusive start in millionths of a frame from master time zero.
    pub start_microframes: i64,
    /// Exclusive end in millionths of a frame from master time zero.
    pub end_microframes: i64,
    /// True when support was clipped to avoid integrating across a cut/master.
    pub clipped_at_boundary: bool,
}

/// Apparent-rotation classification at one constant angular velocity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApparentRotationAssessment {
    /// Directional aliasing classification.
    pub direction: ApparentDirection,
    /// Exposure/readability classification.
    pub exposure_legibility: ExposureLegibility,
    /// Aliased signed marking cycles per frame in `[-0.5, 0.5]`.
    pub apparent_cycles_per_frame: f64,
}

/// Directional result of sampling the rotational cue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ApparentDirection {
    /// Marking motion retains the physical direction and remains visible.
    Faithful,
    /// Sampling reverses the apparent sign.
    FalseReversal,
    /// The marking is nearly stationary between frames.
    FrozenMarking,
}

/// Shutter/readability result for the rotational cue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ExposureLegibility {
    /// No named shutter or near-Nyquist hazard was detected.
    Resolved,
    /// Near-Nyquist phase jumps with short exposure may flicker.
    ObjectionableFlickerRisk,
    /// One shutter integrates at least one full marking period.
    ShutterSmear,
}

/// Evaluate temporal aliasing of a visualization-only rotational cue.
pub fn assess_apparent_rotation(
    angular_velocity_rad_s: f64,
    marking_frequency: u16,
    shutter_duration_microframes: u32,
) -> Result<ApparentRotationAssessment, CinematicBriefError> {
    assess_apparent_rotation_at_rate(
        angular_velocity_rad_s,
        marking_frequency,
        shutter_duration_microframes,
        24,
        1,
    )
}

/// Evaluate temporal aliasing at an explicit rational presentation rate.
pub fn assess_apparent_rotation_at_rate(
    angular_velocity_rad_s: f64,
    marking_frequency: u16,
    shutter_duration_microframes: u32,
    frames_per_second_numerator: u32,
    frames_per_second_denominator: u32,
) -> Result<ApparentRotationAssessment, CinematicBriefError> {
    if !angular_velocity_rad_s.is_finite()
        || marking_frequency == 0
        || shutter_duration_microframes > 1_000_000
        || frames_per_second_numerator == 0
        || frames_per_second_denominator == 0
    {
        return Err(CinematicBriefError::InvalidApparentRotationInput);
    }
    let cycles_per_frame = angular_velocity_rad_s
        * f64::from(marking_frequency)
        * f64::from(frames_per_second_denominator)
        / (TAU * f64::from(frames_per_second_numerator));
    let apparent = cycles_per_frame - cycles_per_frame.round();
    let actual_sign = angular_velocity_rad_s.signum();
    let apparent_sign = apparent.signum();
    let shutter_cycles =
        cycles_per_frame.abs() * f64::from(shutter_duration_microframes) / 1_000_000.0;
    let direction = if apparent.abs() <= 1.0e-6 {
        ApparentDirection::FrozenMarking
    } else if actual_sign * apparent_sign < 0.0 {
        ApparentDirection::FalseReversal
    } else {
        ApparentDirection::Faithful
    };
    let exposure_legibility = if shutter_cycles >= 1.0 {
        ExposureLegibility::ShutterSmear
    } else if apparent.abs() >= 0.4 && shutter_cycles < 0.5 {
        ExposureLegibility::ObjectionableFlickerRisk
    } else {
        ExposureLegibility::Resolved
    };
    Ok(ApparentRotationAssessment {
        direction,
        exposure_legibility,
        apparent_cycles_per_frame: apparent,
    })
}

fn validate_shot(
    shot: &ReferenceShotInput,
    trajectory_start: u64,
    trajectory_end: u64,
) -> Result<(), CinematicBriefError> {
    if shot.camera.is_empty() {
        return Err(CinematicBriefError::MissingCameraKeyframe);
    }
    if shot.optics.focal_length_um == 0
        || shot.optics.f_number_milli == 0
        || shot.optics.shutter_open_microframes < -500_000
        || shot.optics.shutter_close_microframes > 500_000
        || shot.optics.shutter_close_microframes <= shot.optics.shutter_open_microframes
    {
        return Err(CinematicBriefError::InvalidOptics);
    }
    let mut previous = None;
    for keyframe in &shot.camera {
        if !shot.frames.contains(keyframe.frame)
            || previous.is_some_and(|frame| keyframe.frame <= frame)
        {
            return Err(CinematicBriefError::InvalidCameraKeyframeOrder);
        }
        if !keyframe.eye_m.finite()
            || !keyframe.target_m.finite()
            || !keyframe.up.finite()
            || !keyframe.focus_distance_m.is_finite()
            || keyframe.focus_distance_m <= 0.0
        {
            return Err(CinematicBriefError::InvalidCameraValue);
        }
        let forward = keyframe.target_m.sub(keyframe.eye_m);
        if forward.norm_squared() <= 1.0e-18
            || keyframe.up.norm_squared() <= 1.0e-18
            || forward.cross(keyframe.up).norm_squared() <= 1.0e-18
        {
            return Err(CinematicBriefError::CameraSingularity);
        }
        previous = Some(keyframe.frame);
    }
    match shot.time_mapping {
        ShotTimeMapping::PhysicalLinear {
            start_tick,
            end_tick_exclusive,
        } if start_tick >= trajectory_start
            && end_tick_exclusive <= trajectory_end
            && end_tick_exclusive > start_tick => {}
        ShotTimeMapping::VisualizationHold {
            source_tick,
            label: _,
        } if source_tick >= trajectory_start && source_tick < trajectory_end => {}
        _ => return Err(CinematicBriefError::TimeOutsideTrajectory),
    }
    Ok(())
}

fn reference_input() -> CinematicBriefInput {
    let shots = [
        reference_shot(ReferenceShotRole::EstablishingProduct, 0, 60, 0, 2_500),
        reference_shot(
            ReferenceShotRole::InclinationAndContactOrbit,
            60,
            120,
            2_500,
            5_000,
        ),
        reference_shot(ReferenceShotRole::MacroPrecession, 120, 192, 5_000, 8_000),
        reference_shot(ReferenceShotRole::TerminalCloseUp, 192, 240, 8_000, 10_000),
    ];
    CinematicBriefInput {
        total_frames: 240,
        total_audio_sample_frames: 480_000,
        simulation_ticks_per_second: 1_000,
        trajectory_start_tick: 0,
        trajectory_end_tick_exclusive: 10_000,
        shots: shots.into_iter().collect(),
        safe_areas: SafeAreaInsetsPermille {
            action: 50,
            title: 100,
        },
        spin_cue: SpinCue {
            brushing_marking_frequency: 2,
            engraved_radial_mark: false,
        },
        censored_policy: CensoredTrajectoryPolicy::LabeledHoldAndAudioTaper {
            maximum_hold_frames: 12,
        },
        frame_time_convention: FrameTimeConvention::IntegerFrameCentersHalfOpenMaster,
        interpolation_domain: BriefInterpolationDomain::StudioFrameRigidPoseAndScalarOptics,
        shutter_boundary_policy: ShutterBoundaryPolicy::ClipAtCutsAndMaster,
        audio_lead_samples: 0,
        muted_review_required: true,
        essential_context_in_audio_only: false,
    }
}

fn reference_shot(
    role: ReferenceShotRole,
    start: u32,
    end: u32,
    simulation_start: u64,
    simulation_end: u64,
) -> ReferenceShotInput {
    let (eye, lighting, audio, focal_length_um, exposure) = match role {
        ReferenceShotRole::EstablishingProduct => (
            BriefVec3 {
                x: 0.20,
                y: -0.24,
                z: 0.16,
            },
            LightingPreset::StudioSoftboxRim,
            AudioPerspective::StudioObserver,
            50_000,
            ExposureIntent::ProtectMetalHighlights,
        ),
        ReferenceShotRole::InclinationAndContactOrbit => (
            BriefVec3 {
                x: 0.13,
                y: -0.16,
                z: 0.075,
            },
            LightingPreset::RakingContactMacro,
            AudioPerspective::ContactMacro,
            85_000,
            ExposureIntent::ProtectGlassHighlights,
        ),
        ReferenceShotRole::MacroPrecession => (
            BriefVec3 {
                x: 0.085,
                y: -0.11,
                z: 0.055,
            },
            LightingPreset::RakingContactMacro,
            AudioPerspective::ContactMacro,
            100_000,
            ExposureIntent::ProtectMetalHighlights,
        ),
        ReferenceShotRole::TerminalCloseUp => (
            BriefVec3 {
                x: 0.060,
                y: -0.075,
                z: 0.040,
            },
            LightingPreset::TerminalGlint,
            AudioPerspective::TerminalDetail,
            120_000,
            ExposureIntent::ProtectGlassHighlights,
        ),
    };
    ReferenceShotInput {
        role,
        frames: FrameRange {
            start,
            end_exclusive: end,
        },
        camera: reference_camera(role, start, end, eye),
        optics: BriefOptics {
            focal_length_um,
            f_number_milli: 4_000,
            shutter_open_microframes: -180_000,
            shutter_close_microframes: 180_000,
            exposure,
        },
        lighting,
        disc_material: DiscMaterialPreset::BrushedTungsten,
        background: BackgroundPreset::NeutralBlackSweep,
        glass_visible: true,
        base_visible: role != ReferenceShotRole::MacroPrecession,
        audio_perspective: audio,
        transition: ShotTransition::HardCut,
        time_mapping: ShotTimeMapping::PhysicalLinear {
            start_tick: simulation_start,
            end_tick_exclusive: simulation_end,
        },
        apparent_rotation_intent: match role {
            ReferenceShotRole::EstablishingProduct => {
                ApparentRotationIntent::ForwardCueAndInclinationReadable
            }
            ReferenceShotRole::InclinationAndContactOrbit => {
                ApparentRotationIntent::ContactOrbitPrimaryNoFalseReversal
            }
            ReferenceShotRole::MacroPrecession => ApparentRotationIntent::WobblePrimaryNoFrozenCue,
            ReferenceShotRole::TerminalCloseUp => {
                ApparentRotationIntent::TerminalSlowdownNoAliasReversal
            }
        },
    }
}

fn reference_camera(
    role: ReferenceShotRole,
    start: u32,
    end: u32,
    eye: BriefVec3,
) -> Vec<BriefCameraKeyframe> {
    let focus_target = match role {
        ReferenceShotRole::EstablishingProduct => FocusTarget::DiscCenter,
        ReferenceShotRole::InclinationAndContactOrbit => FocusTarget::ContactPoint,
        ReferenceShotRole::MacroPrecession | ReferenceShotRole::TerminalCloseUp => {
            FocusTarget::DiscRim
        }
    };
    vec![
        BriefCameraKeyframe {
            frame: start,
            eye_m: eye,
            target_m: BriefVec3 {
                x: 0.0,
                y: 0.0,
                z: 0.015,
            },
            up: BriefVec3 {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
            focus_distance_m: 0.25,
            focus_target,
        },
        BriefCameraKeyframe {
            frame: end - 1,
            eye_m: BriefVec3 {
                x: eye.x * 0.96,
                y: eye.y,
                z: eye.z * 0.98,
            },
            target_m: BriefVec3 {
                x: 0.0,
                y: 0.0,
                z: 0.012,
            },
            up: BriefVec3 {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
            focus_distance_m: 0.24,
            focus_target,
        },
    ]
}

/// Stable, actionable brief refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CinematicBriefError {
    /// Frame range was empty or reversed.
    EmptyFrameRange,
    /// Video/audio duration disagreed with the deliverable.
    InvalidMasterTimeline,
    /// Simulation clock or available range was invalid.
    InvalidTrajectoryClock,
    /// No shots were supplied.
    MissingShots,
    /// More shots than nonempty master-frame ranges were supplied.
    TooManyShots {
        /// Maximum possible number of nonempty shots.
        maximum: u32,
        /// Supplied shot count.
        got: u64,
    },
    /// Ordered shots overlap or leave a gap.
    ShotGapOrOverlap {
        /// Required next start.
        expected_start: u32,
        /// Supplied start/end boundary.
        got: u32,
    },
    /// Shot extends beyond the master.
    ShotOutsideMaster,
    /// No camera keyframe was supplied.
    MissingCameraKeyframe,
    /// A shot supplied more unique keyframes than frames in its range.
    TooManyCameraKeyframes {
        /// Maximum possible number of strictly ordered keyframes.
        maximum: u32,
        /// Supplied keyframe count.
        got: u64,
    },
    /// Keyframe frames were outside the shot or not strictly increasing.
    InvalidCameraKeyframeOrder,
    /// Camera coordinates or focus distance were non-finite/invalid.
    InvalidCameraValue,
    /// Eye/target/up geometry was singular.
    CameraSingularity,
    /// Lens, aperture, or shutter was invalid.
    InvalidOptics,
    /// Physical/hold mapping escaped the available trajectory.
    TimeOutsideTrajectory,
    /// Action/title safe regions were invalid.
    InvalidSafeAreas,
    /// Rotation cue had no visual frequency.
    InvalidSpinCue,
    /// The visual story depended on audio-only meaning.
    AudioOnlyMeaning,
    /// The censored-input policy allowed no labeled hold frames.
    InvalidCensoredPolicy,
    /// The frozen v1 master does not admit a non-zero A/V lead or lag.
    UnsupportedAudioLeadLag,
    /// A physical shot followed a visualization-only censored hold.
    PhysicalShotAfterCensoredHold,
    /// A censored hold did not use the last available trajectory state.
    HoldIsNotTerminalState,
    /// The aggregate visualization-only hold exceeded its declared bound.
    CensoredHoldTooLong {
        /// Declared maximum hold length.
        maximum: u32,
        /// Aggregate requested hold length.
        got: u32,
    },
    /// Frame does not belong to the master.
    FrameOutsideMaster(u32),
    /// Audio sample frame does not belong to the half-open master.
    AudioSampleOutsideMaster(u64),
    /// Apparent-rotation inputs were non-finite or out of bounds.
    InvalidApparentRotationInput,
}

impl CinematicBriefError {
    /// Stable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::EmptyFrameRange => "cinematic-brief-empty-frame-range",
            Self::InvalidMasterTimeline => "cinematic-brief-invalid-master-timeline",
            Self::InvalidTrajectoryClock => "cinematic-brief-invalid-trajectory-clock",
            Self::MissingShots => "cinematic-brief-missing-shots",
            Self::TooManyShots { .. } => "cinematic-brief-too-many-shots",
            Self::ShotGapOrOverlap { .. } => "cinematic-brief-shot-gap-or-overlap",
            Self::ShotOutsideMaster => "cinematic-brief-shot-outside-master",
            Self::MissingCameraKeyframe => "cinematic-brief-missing-camera-keyframe",
            Self::TooManyCameraKeyframes { .. } => "cinematic-brief-too-many-camera-keyframes",
            Self::InvalidCameraKeyframeOrder => "cinematic-brief-camera-keyframe-order",
            Self::InvalidCameraValue => "cinematic-brief-invalid-camera-value",
            Self::CameraSingularity => "cinematic-brief-camera-singularity",
            Self::InvalidOptics => "cinematic-brief-invalid-optics",
            Self::TimeOutsideTrajectory => "cinematic-brief-time-outside-trajectory",
            Self::InvalidSafeAreas => "cinematic-brief-invalid-safe-areas",
            Self::InvalidSpinCue => "cinematic-brief-invalid-spin-cue",
            Self::AudioOnlyMeaning => "cinematic-brief-audio-only-meaning",
            Self::InvalidCensoredPolicy => "cinematic-brief-invalid-censored-policy",
            Self::UnsupportedAudioLeadLag => "cinematic-brief-unsupported-audio-lead-lag",
            Self::PhysicalShotAfterCensoredHold => {
                "cinematic-brief-physical-shot-after-censored-hold"
            }
            Self::HoldIsNotTerminalState => "cinematic-brief-hold-not-terminal-state",
            Self::CensoredHoldTooLong { .. } => "cinematic-brief-censored-hold-too-long",
            Self::FrameOutsideMaster(_) => "cinematic-brief-frame-outside-master",
            Self::AudioSampleOutsideMaster(_) => "cinematic-brief-audio-outside-master",
            Self::InvalidApparentRotationInput => "cinematic-brief-invalid-rotation-input",
        }
    }
}

impl fmt::Display for CinematicBriefError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {self:?}", self.code())
    }
}

impl std::error::Error for CinematicBriefError {}
