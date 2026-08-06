//! Independent, read-only verification of cinematic frame and audio bundles.
//!
//! The producer's success flag and its manifest are deliberately not treated
//! as an oracle.  A [`CinematicFinalizationPlan`] reconstructs the expected
//! frame/segment/role inventory from the admitted composition, quality
//! profile, render plan, and original prepared exposures.  Verification then
//! consumes persisted bytes through independently pinned codecs and publishes
//! only a deterministic report.
//!
//! This gate proves artifact completeness and internal compatibility.  It
//! does not authenticate caller-supplied hashes, judge visual taste, prove
//! that a biased pixel transform was executed by the named implementation, or
//! promote the underlying mechanics/acoustics to experimental authority.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use fs_blake3::{Blake3, ContentHash, DomainHasher, hash_domain};
use fs_evidence::cinematic::{
    CinematicArtifactKind, CinematicAuthorityClass, CinematicAuthorityError,
    CinematicAuthorityRecord, CinematicClock, CinematicClockDomain, CinematicNoClaim,
    CinematicTransformDisposition, CinematicUnitContract, DeclaredAcousticCalibrationReceipt,
    SoundAuthority, required_no_claims,
};
use fs_evidence::cinematic_brief::{CINEMATIC_BRIEF_IDENTITY_VERSION, CinematicBrief};
use fs_evidence::cinematic_budget::{
    AovPreset, CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION, CinematicQualityProfile, DenoisePolicy,
};
use fs_evidence::cinematic_config::{CINEMATIC_CONFIG_SCHEMA_VERSION, CinematicConfig};
use fs_evidence::cinematic_sound::{
    SOUND_SYNTHESIS_SCHEMA_VERSION, SoundSynthesisConfig, SoundSynthesisReceipt,
};
use fs_exec::Cx;
use fs_img::{
    ExpectedFrameArtifact, ExrInspectLimits, ExrInspection, ExrRawFrameSemanticLimits,
    FrameArtifactDescriptor, FrameArtifactFormat, FrameArtifactKey, FrameArtifactRole,
    FrameChannel, FrameChannelType, FrameSamplingStats, FrameSequenceContext, FrameSequenceError,
    FrameSequenceLimits, FrameSequenceManifest, FrameSequenceState, ImgError, PixelType, PngColor,
    PngInspectLimits, SOURCE_ARTIFACT_HASH_ATTRIBUTE, inspect_exr_with_poll, inspect_png_with_poll,
    validate_exr_raw_frame_payload_against_with_poll, verify_exr_float_channel_constant_with_poll,
};
use fs_math::det;
use fs_render::aov::{
    CINEMATIC_AOV_CHANNEL_SEMANTICS, CINEMATIC_AOV_INVALID_SEMANTICS,
    CINEMATIC_AOV_PALETTE_ZERO_SEMANTICS, CINEMATIC_AOV_SEMANTICS_VERSION, CinematicAovConfig,
    CinematicAovError, CinematicAovLimits, CinematicAovPalette, CinematicAovProfile,
    CinematicAovProvenance, cinematic_export_metadata_payload_bound,
    cinematic_render_semantics_versions, encode_material_palette_with_poll,
    encode_object_palette_with_poll,
};
use fs_render::camera::CutSide;
use fs_render::motion::{ShutterConvention, ShutterDistribution};
use fs_render::tracer::{DirectStrategy, Sampler, TracerError};

use crate::audio_artifact::{
    AUDIO_ARTIFACT_SCHEMA_VERSION, AudioArtifactBudget, AudioArtifactError, AudioArtifactManifest,
    AudioArtifactRole, AudioMeters, AudioSignalPath, WavSampleEncoding, decode_stereo_wav,
    measure_audio,
};
use crate::audio_resampling::{
    AudioVideoAlignment, AudioVideoSyncMarker, ResampledAudioEvent, validate_resampled_audio_event,
};
use crate::render_checkpoint::euler_render_checkpoint_frame_identity;
use crate::render_scene_bridge::EulerCinematicScene;
use crate::render_sharding::{EulerRenderFrameInput, EulerUniformRenderPlan};

/// Canonical report identity schema.
pub const CINEMATIC_FINALIZATION_REPORT_VERSION: u16 = 1;
/// Canonical independently reconstructed plan identity schema.
pub const CINEMATIC_FINALIZATION_PLAN_VERSION: u16 = 1;
/// Domain separating finalization reports from every inspected child artifact.
pub const CINEMATIC_FINALIZATION_REPORT_DOMAIN: &str =
    "org.frankensim.euler-cinematic.finalization-report.v1";
const CINEMATIC_FINALIZATION_PLAN_DOMAIN: &str =
    "org.frankensim.euler-cinematic.finalization-plan.v1";

const REPORT_MAGIC: &[u8; 8] = b"FSFINL01";
const ALIGNMENT_RECEIPT_MAGIC: &[u8; 8] = b"FSAVSYN1";
const EVENT_RECEIPT_MAGIC: &[u8; 8] = b"FSAEVT01";
const AUDIO_MANIFEST_RECEIPT_MAGIC: &[u8; 8] = b"FSAMNF01";
const RECEIPT_CODEC_VERSION: u16 = 1;
const ALIGNMENT_RECEIPT_DOMAIN: &str = "org.frankensim.euler-cinematic.av-alignment.v1";
const EVENT_RECEIPT_DOMAIN: &str = "org.frankensim.euler-cinematic.audio-events.v1";
const AUDIO_MANIFEST_RECEIPT_DOMAIN: &str =
    "org.frankensim.euler-cinematic.audio-manifest-receipt.v1";
const AUDIO_MANIFEST_IDENTITY_DOMAIN: &str = "org.frankensim.euler-cinematic.audio-manifest.v1";
const CHANNEL_RECEIPT_IDENTITY_DOMAIN: &str =
    "org.frankensim.euler-cinematic.audio-channel-receipt.v1";
const HASH_POLL_BYTES: usize = 64 * 1024;
const MAX_AUTHORITY_RECORD_WIRE_BYTES: u64 = 64 * 1024;
const MAX_EXACT_AOV_PALETTE_INDEX: u32 = 1 << 24;
const ALLOWED_AOV_VALIDITY_BITS: u32 = fs_render::aov::validity::PRIMARY
    | fs_render::aov::validity::ALBEDO
    | fs_render::aov::validity::PREVIOUS_MOTION
    | fs_render::aov::validity::OBJECT_ID
    | fs_render::aov::validity::MATERIAL_ID
    | fs_render::aov::validity::CONTRIBUTION_SPLIT;

/// Maximum encoded bytes reserved for each image role.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicFrameArtifactCeilings {
    /// One raw EXR master.
    pub raw_master_bytes: u64,
    /// One separately labeled denoised EXR.
    pub denoised_intermediate_bytes: u64,
    /// One display-referred PNG preview.
    pub display_preview_bytes: u64,
}

impl CinematicFrameArtifactCeilings {
    fn validate(self) -> Result<(), CinematicFinalizationPlanError> {
        if self.raw_master_bytes == 0
            || self.denoised_intermediate_bytes == 0
            || self.display_preview_bytes == 0
        {
            Err(CinematicFinalizationPlanError::InvalidLimit(
                "frame artifact ceiling",
            ))
        } else {
            Ok(())
        }
    }
}

/// Independent verification ceilings.  These cap bytes before parsing or
/// image-sized allocation; the inspectors retain metadata only.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicFinalizationLimits {
    /// Admission limits for the persisted frame manifest.
    pub frame_sequence: FrameSequenceLimits,
    /// Structural EXR limits.
    pub exr: ExrInspectLimits,
    /// Structural PNG limits.
    pub png: PngInspectLimits,
    /// Strict WAV decode/meter limits.
    pub audio: AudioArtifactBudget,
    /// Maximum canonical audio-manifest receipt bytes.
    pub max_audio_manifest_bytes: u64,
    /// Maximum retained A/V sync markers inspected.
    pub max_sync_markers: u32,
    /// Maximum exact resampled event receipts inspected.
    pub max_audio_events: u32,
    /// Caller ceiling for any authority envelope. Verification also applies
    /// a fixed 64-KiB schema-parser ceiling.
    pub max_authority_record_bytes: u64,
    /// Maximum aggregate payload bytes inspected in one call.
    pub max_bundle_bytes: u64,
}

impl CinematicFinalizationLimits {
    fn validate(self) -> Result<(), CinematicFinalizationPlanError> {
        if self.max_sync_markers == 0
            || self.max_audio_events == 0
            || self.max_authority_record_bytes == 0
            || self.max_audio_manifest_bytes == 0
            || self.max_bundle_bytes == 0
            || self.exr.max_input_bytes == 0
            || self.exr.max_header_bytes == 0
            || self.exr.max_decoded_bytes == 0
            || self.exr.max_metadata_bytes == 0
            || self.png.max_input_bytes == 0
            || self.png.max_decoded_bytes == 0
        {
            return Err(CinematicFinalizationPlanError::InvalidLimit(
                "finalization limit",
            ));
        }
        Ok(())
    }
}

/// Whether a successful integrity check is eligible to satisfy the frozen
/// final-4K delivery contract.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CinematicFinalizationTarget {
    /// A deliberately small cross-codec fixture.  It can pass integrity but
    /// can never be represented as final delivery.
    IntegrityFixture,
    /// Integrity-only production target. This includes non-final quality tiers
    /// and a nominal final raster produced by a path that lacks a required
    /// delivery-completion receipt.
    NonFinal,
    /// Complete frozen 4K master target.
    Final4k,
}

impl CinematicFinalizationTarget {
    const fn tag(self) -> u8 {
        match self {
            Self::IntegrityFixture => 1,
            Self::NonFinal => 2,
            Self::Final4k => 3,
        }
    }
}

/// Inputs used to reconstruct a production Euler-disc artifact inventory.
pub struct EulerCinematicFinalizationPlanInput<'a, 'scene, 'frame> {
    /// Complete admitted composition.
    pub configuration: &'a CinematicConfig,
    /// Exact admitted image/audio resource profile.
    pub quality_profile: &'a CinematicQualityProfile,
    /// Exact creative/timeline brief.
    pub brief: &'a CinematicBrief,
    /// Immutable uniform render partition.
    pub render_plan: &'a EulerUniformRenderPlan,
    /// Exact admitted scene used to prepare and render the exposures. Camera
    /// shot identities and AOV palettes are derived from this object rather
    /// than accepted as caller assertions.
    pub scene: &'a EulerCinematicScene<'scene>,
    /// Original prepared exposures, needed because the render plan omits cut
    /// and shutter metadata.
    pub render_frames: &'a [EulerRenderFrameInput<'frame>],
    /// AOV resource/configuration limits used by the renderer.
    pub aov_limits: CinematicAovLimits,
    /// Sequence resource ceilings.
    pub sequence_limits: FrameSequenceLimits,
    /// Per-role output reservations.
    pub artifact_ceilings: CinematicFrameArtifactCeilings,
    /// Independently supplied build identity.
    pub build_identity: ContentHash,
    /// Complete admitted sound configuration.
    pub sound_configuration: &'a SoundSynthesisConfig,
    /// Expected identity of the exact pre-encoding stereo samples, retained
    /// independently from the produced audio manifest.
    pub expected_audio_source_signal_identity: ContentHash,
    /// Expected canonical dry/spatialized channel-path identity, retained
    /// independently from the produced audio manifest.
    pub expected_audio_channel_layout_identity: ContentHash,
    /// Expected dry-mix identity. This is `None` for a separately spatialized
    /// stereo source and `Some` for the canonical dry-stem mixer.
    pub expected_audio_mix_identity: Option<ContentHash>,
    /// Independently retained expected event-placement receipts.  Exact
    /// reproduction requires the source-grid origin, which is not retained by
    /// an `AudioExcitationEvent` alone.
    pub expected_audio_events: &'a [ResampledAudioEvent],
    /// Maximum event receipts copied into the immutable finalization plan.
    pub max_audio_events: u32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct RawFrameExpectation {
    authority_source_identity: ContentHash,
    authority_transform_identity: ContentHash,
    required_attributes: BTreeMap<String, Arc<[u8]>>,
    object_palette_entries: u32,
    material_palette_entries: u32,
    expected_uniform_spp: Option<u32>,
}

/// Immutable, independently reconstructed finalization oracle.
#[derive(Clone, Debug, PartialEq)]
pub struct CinematicFinalizationPlan {
    target: CinematicFinalizationTarget,
    identity: ContentHash,
    configuration_identity: ContentHash,
    image_configuration_identity: ContentHash,
    build_identity: ContentHash,
    profile_identity: ContentHash,
    brief_identity: ContentHash,
    render_plan_identity: ContentHash,
    image_pipeline_identity: ContentHash,
    expected_sequence: FrameSequenceManifest,
    raw_frames: BTreeMap<FrameArtifactKey, RawFrameExpectation>,
    sound_receipt: SoundSynthesisReceipt,
    expected_audio_source_signal_identity: ContentHash,
    expected_audio_channel_layout_identity: ContentHash,
    expected_audio_mix_identity: Option<ContentHash>,
    expected_acoustic_calibration: Option<DeclaredAcousticCalibrationReceipt>,
    expected_audio_events: Vec<ResampledAudioEvent>,
    total_video_frames: u32,
    total_audio_sample_frames: u64,
    frames_per_second: u32,
    audio_sample_rate_hz: u32,
    cut_frame_boundaries: Vec<u32>,
}

impl CinematicFinalizationPlan {
    /// Reconstruct the expected artifact set without consulting producer
    /// manifest contents.
    pub fn try_from_euler_disc(
        input: EulerCinematicFinalizationPlanInput<'_, '_, '_>,
        cx: &Cx<'_>,
    ) -> Result<Self, CinematicFinalizationPlanError> {
        checkpoint_plan(cx)?;
        input.artifact_ceilings.validate()?;
        if input.max_audio_events == 0 {
            return Err(CinematicFinalizationPlanError::InvalidLimit(
                "maximum audio events",
            ));
        }
        let event_count = u32::try_from(input.expected_audio_events.len()).map_err(|_| {
            CinematicFinalizationPlanError::ResourceLimit {
                resource: "expected audio events",
                requested: u64::MAX,
                limit: u64::from(input.max_audio_events),
            }
        })?;
        if event_count > input.max_audio_events {
            return Err(CinematicFinalizationPlanError::ResourceLimit {
                resource: "expected audio events",
                requested: u64::from(event_count),
                limit: u64::from(input.max_audio_events),
            });
        }
        if is_zero(input.build_identity) {
            return Err(CinematicFinalizationPlanError::MissingIdentity(
                "build identity",
            ));
        }
        let profile = input.quality_profile.input();
        if profile.frames_per_second != 24 {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "quality-profile frame rate",
            ));
        }
        let aov_profile = aov_profile(profile.aov_preset);
        let settings = input.render_plan.settings();
        let tile_layout = input
            .render_plan
            .tile_layout()
            .map_err(|_| CinematicFinalizationPlanError::Incompatible("render tile layout"))?;
        if settings.width != profile.width_pixels
            || settings.height != profile.height_pixels
            || settings.spp < profile.spp_floor
            || settings.spp > profile.spp_ceiling
            || settings.max_depth != u32::from(profile.max_path_depth)
            || tile_layout.tile_width() != u32::from(profile.tile_width)
            || tile_layout.tile_height() != u32::from(profile.tile_height)
        {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "render settings and quality profile",
            ));
        }
        if profile.first_frame != 0 || profile.frame_count != input.brief.total_frames() {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "partial-range A/V finalization requires an explicit range clock",
            ));
        }
        let expected_artifact_count = preflight_expected_inventory(
            input.render_plan.segments().len(),
            aov_profile,
            profile.denoise_policy,
            profile.width_pixels,
            profile.height_pixels,
            profile.output_ceiling_bytes,
            input.artifact_ceilings,
            input.sequence_limits,
            input.aov_limits,
        )?;
        let profile_identity = input.quality_profile.identity();
        let brief_identity = input.brief.identity();
        let config_input = input.configuration.input();
        let sound_input = input.sound_configuration.input();
        if config_input.render_budget_profile.is_none_or(|reference| {
            reference.identity() != profile_identity
                || reference.version() != u32::from(CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION)
        }) {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "render budget profile binding",
            ));
        }
        if config_input.audio_budget_profile.is_none_or(|reference| {
            reference.identity() != profile_identity
                || reference.version() != u32::from(CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION)
        }) {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "audio budget profile binding",
            ));
        }
        if config_input.trajectory != sound_input.trajectory
            || config_input.timeline != sound_input.timeline
            || config_input.timeline.identity() != brief_identity
            || config_input.timeline.version() != u32::from(CINEMATIC_BRIEF_IDENTITY_VERSION)
            || config_input.audio_excitation != sound_input.excitation
            || config_input.sound_model != sound_input.sound_model
            || config_input.microphone != sound_input.microphone
            || config_input.room != sound_input.room
            || input.render_plan.sequence_identity() != brief_identity
            || config_input.trajectory.identity() != input.render_plan.source_trajectory_identity()
            || input.sound_configuration.receipt().trajectory_identity
                != input.render_plan.source_trajectory_identity()
            || input.render_plan.source_configuration_identity()
                != input.scene.source_configuration_identity()
            || input.render_plan.source_trajectory_identity()
                != input.scene.source_trajectory_identity()
            || input.render_plan.scene_identity() != input.scene.scene_identity()
        {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "composition source identities",
            ));
        }
        let expected_video_clock = CinematicClock::try_new(
            CinematicClockDomain::Video,
            profile.frames_per_second,
            1,
            0,
            i64::from(input.brief.total_frames()),
        )
        .map_err(|_| CinematicFinalizationPlanError::Incompatible("brief video clock"))?;
        let expected_audio_clock = CinematicClock::try_new(
            CinematicClockDomain::Audio,
            48_000,
            1,
            0,
            i64::try_from(input.brief.total_audio_sample_frames()).map_err(|_| {
                CinematicFinalizationPlanError::ArithmeticOverflow("brief audio clock")
            })?,
        )
        .map_err(|_| CinematicFinalizationPlanError::Incompatible("brief audio clock"))?;
        if sound_input.video_clock != expected_video_clock
            || sound_input.audio_clock != expected_audio_clock
            || input.brief.audio_lead_samples() != 0
        {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "sound and brief master clocks",
            ));
        }
        if config_input.seed != Some(settings.seed) {
            return Err(CinematicFinalizationPlanError::Incompatible("render seed"));
        }
        let first = u64::from(profile.first_frame);
        let end = first.checked_add(u64::from(profile.frame_count)).ok_or(
            CinematicFinalizationPlanError::ArithmeticOverflow("profile frame range"),
        )?;
        let profile_frame_count = usize::try_from(profile.frame_count).map_err(|_| {
            CinematicFinalizationPlanError::ArithmeticOverflow("profile frame count")
        })?;
        if input.render_plan.frames().len() != profile_frame_count
            || input.render_frames.len() != input.render_plan.frames().len()
            || end > u64::from(input.brief.total_frames())
        {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "render frame inventory",
            ));
        }
        for (index, frame) in input.render_plan.frames().iter().enumerate() {
            if index.is_multiple_of(1_024) {
                checkpoint_plan(cx)?;
            }
            let expected_ordinal = first
                .checked_add(u64::try_from(index).map_err(|_| {
                    CinematicFinalizationPlanError::ArithmeticOverflow("render frame position")
                })?)
                .ok_or(CinematicFinalizationPlanError::ArithmeticOverflow(
                    "render frame ordinal",
                ))?;
            if frame.frame_ordinal() != expected_ordinal {
                return Err(CinematicFinalizationPlanError::Incompatible(
                    "render frame inventory",
                ));
            }
        }

        let frame_inputs = index_frame_inputs(input.render_frames, input.render_plan, cx)?;
        validate_prepared_frames_against_brief(
            &frame_inputs,
            input.render_plan,
            input.brief,
            profile.frames_per_second,
            cx,
        )?;
        let context = FrameSequenceContext::try_new(
            brief_identity,
            input.render_plan.source_trajectory_identity(),
            input.configuration.image_identity(),
            input.render_plan.scene_identity(),
            input.build_identity,
            profile_identity,
        )
        .map_err(CinematicFinalizationPlanError::Sequence)?;
        let palette = CinematicAovPalette::try_from_scene(
            input.scene.scene(),
            input.aov_limits,
            matches!(aov_profile, CinematicAovProfile::FinalDiagnostic),
            cx,
        )
        .map_err(map_aov_plan_error)?;
        let object_palette_entries = u32::try_from(palette.object_ids().len()).map_err(|_| {
            CinematicFinalizationPlanError::ArithmeticOverflow("object palette entries")
        })?;
        let material_palette_entries =
            u32::try_from(palette.material_identities().len()).map_err(|_| {
                CinematicFinalizationPlanError::ArithmeticOverflow("material palette entries")
            })?;
        let metadata_bound = cinematic_export_metadata_payload_bound(
            palette.object_ids().len(),
            palette.material_identities().len(),
        )
        .map_err(map_aov_plan_error)?;
        if metadata_bound > input.aov_limits.max_export_metadata_bytes() {
            return Err(CinematicFinalizationPlanError::ResourceLimit {
                resource: "AOV export metadata bytes",
                requested: metadata_bound,
                limit: input.aov_limits.max_export_metadata_bytes(),
            });
        }
        let object_palette_metadata: Arc<[u8]> =
            encode_object_palette_with_poll(palette.object_ids(), || cx.checkpoint().is_ok())
                .map_err(map_aov_plan_error)?
                .into_bytes()
                .into();
        let material_palette_metadata: Arc<[u8]> =
            encode_material_palette_with_poll(palette.material_identities(), || {
                cx.checkpoint().is_ok()
            })
            .map_err(map_aov_plan_error)?
            .into_bytes()
            .into();
        let render_semantics_versions: Arc<[u8]> =
            cinematic_render_semantics_versions().into_bytes().into();
        let raw_channels = aov_channels(aov_profile)?;
        let rgb_float = rgb_channels(FrameChannelType::Float32)?;
        let rgb_u16 = rgb_channels(FrameChannelType::Uint16)?;
        let mut expected = Vec::new();
        expected
            .try_reserve_exact(expected_artifact_count)
            .map_err(|_| CinematicFinalizationPlanError::Capacity("expected frame artifacts"))?;
        let mut raw_frames = BTreeMap::new();
        for segment in input.render_plan.segments() {
            checkpoint_plan(cx)?;
            let prepared = frame_inputs.get(&segment.frame_ordinal()).ok_or(
                CinematicFinalizationPlanError::Incompatible("prepared frame inventory"),
            )?;
            let segment_index = usize::try_from(segment.segment_index()).map_err(|_| {
                CinematicFinalizationPlanError::ArithmeticOverflow("prepared segment index")
            })?;
            let (shutter, cut_side) = input
                .scene
                .prepared_segment_shard_binding(prepared.prepared(), segment_index)
                .map_err(|_| {
                    CinematicFinalizationPlanError::Incompatible("prepared segment scene binding")
                })?;
            let exposure = input
                .scene
                .camera()
                .admit_shutter(cx, shutter, cut_side)
                .map_err(|_| {
                    CinematicFinalizationPlanError::Incompatible("prepared segment camera shot")
                })?;
            let recomputed =
                euler_render_checkpoint_frame_identity(prepared.prepared(), segment_index)
                    .map_err(|_| {
                        CinematicFinalizationPlanError::Incompatible("prepared frame identity")
                    })?;
            if recomputed != segment.frame_identity()
                || prepared.prepared().scene_identity() != input.render_plan.scene_identity()
            {
                return Err(CinematicFinalizationPlanError::Incompatible(
                    "prepared frame binding",
                ));
            }
            let segment_index_u32 = u32::try_from(segment.segment_index()).map_err(|_| {
                CinematicFinalizationPlanError::ArithmeticOverflow("frame segment index")
            })?;
            let frame_time_s =
                segment.frame_ordinal() as f64 / f64::from(profile.frames_per_second);
            let raw_descriptor = FrameArtifactDescriptor::try_new(
                segment.frame_ordinal(),
                segment_index_u32,
                FrameArtifactRole::RawMaster,
                frame_time_s,
                FrameArtifactFormat::OpenExr,
                profile.width_pixels,
                profile.height_pixels,
                raw_channels.clone(),
                FrameSamplingStats::Uniform { spp: settings.spp },
            )
            .map_err(CinematicFinalizationPlanError::Sequence)?;
            let raw_key = raw_descriptor.key();
            expected.push(
                ExpectedFrameArtifact::try_new(
                    raw_descriptor,
                    input.artifact_ceilings.raw_master_bytes,
                    None,
                )
                .map_err(CinematicFinalizationPlanError::Sequence)?,
            );
            let provenance = frame_provenance(
                segment.frame_ordinal(),
                input.brief.total_frames(),
                profile.frames_per_second,
                input.render_plan.source_trajectory_identity(),
                input.render_plan.scene_identity(),
                input.configuration.composition_identity(),
            )?;
            let aov_config = CinematicAovConfig::new(aov_profile, provenance, input.aov_limits);
            raw_frames.insert(
                raw_key,
                RawFrameExpectation {
                    authority_source_identity: segment.frame_identity(),
                    authority_transform_identity: aov_config.identity(),
                    required_attributes: raw_attributes(
                        aov_config,
                        settings,
                        exposure.shot_id(),
                        cut_side,
                        shutter,
                        &object_palette_metadata,
                        &material_palette_metadata,
                        &render_semantics_versions,
                    ),
                    object_palette_entries,
                    material_palette_entries,
                    expected_uniform_spp: matches!(
                        aov_profile,
                        CinematicAovProfile::FinalDiagnostic
                    )
                    .then_some(settings.spp),
                },
            );
            let mut preview_source = raw_key;
            if profile.denoise_policy == DenoisePolicy::SeparateBiasedDerivative {
                let descriptor = FrameArtifactDescriptor::try_new(
                    segment.frame_ordinal(),
                    segment_index_u32,
                    FrameArtifactRole::DenoisedIntermediate,
                    frame_time_s,
                    FrameArtifactFormat::OpenExr,
                    profile.width_pixels,
                    profile.height_pixels,
                    rgb_float.clone(),
                    FrameSamplingStats::Uniform { spp: settings.spp },
                )
                .map_err(CinematicFinalizationPlanError::Sequence)?;
                preview_source = descriptor.key();
                expected.push(
                    ExpectedFrameArtifact::try_new(
                        descriptor,
                        input.artifact_ceilings.denoised_intermediate_bytes,
                        Some(raw_key),
                    )
                    .map_err(CinematicFinalizationPlanError::Sequence)?,
                );
            }
            let preview = FrameArtifactDescriptor::try_new(
                segment.frame_ordinal(),
                segment_index_u32,
                FrameArtifactRole::DisplayPreview,
                frame_time_s,
                FrameArtifactFormat::Png16,
                profile.width_pixels,
                profile.height_pixels,
                rgb_u16.clone(),
                FrameSamplingStats::Uniform { spp: settings.spp },
            )
            .map_err(CinematicFinalizationPlanError::Sequence)?;
            expected.push(
                ExpectedFrameArtifact::try_new(
                    preview,
                    input.artifact_ceilings.display_preview_bytes,
                    Some(preview_source),
                )
                .map_err(CinematicFinalizationPlanError::Sequence)?,
            );
        }
        let expected_sequence = FrameSequenceManifest::try_new_with_poll(
            context,
            expected,
            input.sequence_limits,
            input.sequence_limits.max_output_bytes(),
            || cx.checkpoint().is_ok(),
        )
        .map_err(map_sequence_plan_error)?;
        // A uniform shard plan can prove byte integrity at 4K, but it carries
        // neither the profile's adaptive stopping receipt nor its explicit
        // temporal-sample completion contract.  Treating resolution/tier alone
        // as final-delivery authority would silently discard those semantics.
        let target = CinematicFinalizationTarget::NonFinal;
        let cut_count = input.brief.shots().len().saturating_sub(1);
        let mut cut_frame_boundaries = Vec::new();
        cut_frame_boundaries
            .try_reserve_exact(cut_count)
            .map_err(|_| CinematicFinalizationPlanError::Capacity("brief cut boundaries"))?;
        for (index, shot) in input.brief.shots().iter().take(cut_count).enumerate() {
            if index.is_multiple_of(64) {
                checkpoint_plan(cx)?;
            }
            cut_frame_boundaries.push(shot.frames().end_exclusive());
        }
        let mut expected_audio_events = Vec::new();
        expected_audio_events
            .try_reserve_exact(input.expected_audio_events.len())
            .map_err(|_| CinematicFinalizationPlanError::Capacity("expected audio events"))?;
        for (index, event) in input.expected_audio_events.iter().enumerate() {
            if index.is_multiple_of(1_024) {
                checkpoint_plan(cx)?;
            }
            expected_audio_events.push(event.clone());
        }
        Self::from_parts(
            target,
            input.configuration.composition_identity(),
            input.configuration.image_identity(),
            input.build_identity,
            profile_identity,
            input.brief.identity(),
            input.render_plan.plan_identity(),
            config_input.image_pipeline.identity(),
            input.scene.source_configuration_identity(),
            expected_sequence,
            raw_frames,
            input.sound_configuration.receipt(),
            input.expected_audio_source_signal_identity,
            input.expected_audio_channel_layout_identity,
            input.expected_audio_mix_identity,
            sound_input.calibration,
            expected_audio_events,
            input.brief.total_frames(),
            input.brief.total_audio_sample_frames(),
            profile.frames_per_second,
            48_000,
            cut_frame_boundaries,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn from_parts(
        target: CinematicFinalizationTarget,
        configuration_identity: ContentHash,
        image_configuration_identity: ContentHash,
        build_identity: ContentHash,
        profile_identity: ContentHash,
        brief_identity: ContentHash,
        render_plan_identity: ContentHash,
        image_pipeline_identity: ContentHash,
        scene_configuration_identity: ContentHash,
        expected_sequence: FrameSequenceManifest,
        raw_frames: BTreeMap<FrameArtifactKey, RawFrameExpectation>,
        sound_receipt: SoundSynthesisReceipt,
        expected_audio_source_signal_identity: ContentHash,
        expected_audio_channel_layout_identity: ContentHash,
        expected_audio_mix_identity: Option<ContentHash>,
        expected_acoustic_calibration: Option<DeclaredAcousticCalibrationReceipt>,
        expected_audio_events: Vec<ResampledAudioEvent>,
        total_video_frames: u32,
        total_audio_sample_frames: u64,
        frames_per_second: u32,
        audio_sample_rate_hz: u32,
        cut_frame_boundaries: Vec<u32>,
        cx: &Cx<'_>,
    ) -> Result<Self, CinematicFinalizationPlanError> {
        checkpoint_plan(cx)?;
        for (index, (name, identity)) in [
            ("configuration identity", configuration_identity),
            ("image configuration identity", image_configuration_identity),
            ("build identity", build_identity),
            ("profile identity", profile_identity),
            ("brief identity", brief_identity),
            ("render plan identity", render_plan_identity),
            ("image pipeline identity", image_pipeline_identity),
            ("scene configuration identity", scene_configuration_identity),
            (
                "sound configuration identity",
                sound_receipt.configuration_identity,
            ),
            (
                "sound trajectory identity",
                sound_receipt.trajectory_identity,
            ),
            (
                "sound excitation identity",
                sound_receipt.excitation_identity,
            ),
            ("sound model identity", sound_receipt.sound_model_identity),
            ("sound timeline identity", sound_receipt.timeline_identity),
            (
                "audio source signal identity",
                expected_audio_source_signal_identity,
            ),
            (
                "audio channel layout identity",
                expected_audio_channel_layout_identity,
            ),
        ]
        .into_iter()
        .enumerate()
        {
            if index.is_multiple_of(16) {
                checkpoint_plan(cx)?;
            }
            if is_zero(identity) {
                return Err(CinematicFinalizationPlanError::MissingIdentity(name));
            }
        }
        if sound_receipt.schema_version != SOUND_SYNTHESIS_SCHEMA_VERSION {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "sound synthesis schema version",
            ));
        }
        if expected_audio_mix_identity.is_some_and(is_zero) {
            return Err(CinematicFinalizationPlanError::MissingIdentity(
                "audio mix identity",
            ));
        }
        let audio_frames_per_video_frame = audio_sample_rate_hz
            .checked_div(frames_per_second)
            .filter(|_| audio_sample_rate_hz % frames_per_second == 0)
            .ok_or(CinematicFinalizationPlanError::Incompatible(
                "master clocks",
            ))?;
        let expected_audio_sample_frames = u64::from(total_video_frames)
            .checked_mul(u64::from(audio_frames_per_video_frame))
            .ok_or(CinematicFinalizationPlanError::ArithmeticOverflow(
                "master clock duration",
            ))?;
        if total_video_frames == 0
            || total_audio_sample_frames == 0
            || frames_per_second == 0
            || audio_sample_rate_hz == 0
            || total_audio_sample_frames != expected_audio_sample_frames
        {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "master clocks",
            ));
        }
        let marker_count = total_video_frames.checked_add(1).ok_or(
            CinematicFinalizationPlanError::ArithmeticOverflow("sync marker count"),
        )?;
        let _ = marker_count;
        u32::try_from(expected_audio_events.len())
            .map_err(|_| CinematicFinalizationPlanError::ArithmeticOverflow("audio event count"))?;
        validate_event_sequence(&expected_audio_events, total_audio_sample_frames, cx)?;
        let mut previous_boundary = None;
        for (index, boundary) in cut_frame_boundaries.iter().copied().enumerate() {
            if index.is_multiple_of(1_024) {
                checkpoint_plan(cx)?;
            }
            if boundary == 0
                || boundary >= total_video_frames
                || previous_boundary.is_some_and(|previous| previous >= boundary)
            {
                return Err(CinematicFinalizationPlanError::Incompatible(
                    "cut frame boundaries",
                ));
            }
            previous_boundary = Some(boundary);
        }
        let snapshot = expected_sequence
            .snapshot_with_poll(|| cx.checkpoint().is_ok())
            .map_err(map_sequence_plan_error)?;
        validate_raw_expectations(&expected_sequence, &raw_frames, cx)?;
        let mut identity_writer = PlanIdentityHasher::new(cx)?;
        identity_writer.update(&CINEMATIC_FINALIZATION_PLAN_VERSION.to_le_bytes())?;
        identity_writer.byte(target.tag())?;
        for identity in [
            configuration_identity,
            image_configuration_identity,
            build_identity,
            profile_identity,
            brief_identity,
            render_plan_identity,
            image_pipeline_identity,
            scene_configuration_identity,
            snapshot.identity(),
            sound_receipt.configuration_identity,
            sound_receipt.trajectory_identity,
            sound_receipt.excitation_identity,
            sound_receipt.sound_model_identity,
            sound_receipt.timeline_identity,
            expected_audio_source_signal_identity,
            expected_audio_channel_layout_identity,
        ] {
            identity_writer.update(identity.as_bytes())?;
        }
        match expected_audio_mix_identity {
            None => identity_writer.byte(0)?,
            Some(identity) => {
                identity_writer.byte(1)?;
                identity_writer.update(identity.as_bytes())?;
            }
        }
        identity_writer.update(&total_video_frames.to_le_bytes())?;
        identity_writer.update(&total_audio_sample_frames.to_le_bytes())?;
        identity_writer.update(&frames_per_second.to_le_bytes())?;
        identity_writer.update(&audio_sample_rate_hz.to_le_bytes())?;
        identity_writer.update(&sound_receipt.schema_version.to_le_bytes())?;
        identity_writer.string(sound_receipt.authority.code())?;
        match expected_acoustic_calibration {
            None => identity_writer.byte(0)?,
            Some(calibration) => {
                identity_writer.byte(1)?;
                identity_writer.update(calibration.dataset_identity().as_bytes())?;
                identity_writer.update(calibration.method_identity().as_bytes())?;
                identity_writer.update(calibration.validity_identity().as_bytes())?;
                identity_writer.update(&calibration.version().to_le_bytes())?;
            }
        }
        hash_raw_expectations(&mut identity_writer, &raw_frames)?;
        hash_resampled_events(&mut identity_writer, &expected_audio_events)?;
        identity_writer.update(
            &u32::try_from(cut_frame_boundaries.len())
                .expect("cut boundaries are indexed by admitted video frames")
                .to_le_bytes(),
        )?;
        for boundary in &cut_frame_boundaries {
            identity_writer.update(&boundary.to_le_bytes())?;
        }
        let identity = identity_writer.finish()?;
        Ok(Self {
            target,
            identity,
            configuration_identity,
            image_configuration_identity,
            build_identity,
            profile_identity,
            brief_identity,
            render_plan_identity,
            image_pipeline_identity,
            expected_sequence,
            raw_frames,
            sound_receipt,
            expected_audio_source_signal_identity,
            expected_audio_channel_layout_identity,
            expected_audio_mix_identity,
            expected_acoustic_calibration,
            expected_audio_events,
            total_video_frames,
            total_audio_sample_frames,
            frames_per_second,
            audio_sample_rate_hz,
            cut_frame_boundaries,
        })
    }

    /// Content identity of every expectation consumed by the verifier.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Delivery-authority class reconstructed for this plan. The current
    /// uniform sharding constructor returns [`CinematicFinalizationTarget::NonFinal`]
    /// even for a 4K integrity check because it has no adaptive/temporal
    /// completion receipt.
    #[must_use]
    pub const fn target(&self) -> CinematicFinalizationTarget {
        self.target
    }

    /// Independently reconstructed expected frame inventory.
    #[must_use]
    pub const fn expected_sequence(&self) -> &FrameSequenceManifest {
        &self.expected_sequence
    }
}

/// Planning/refusal failures detected before artifact inspection begins.
#[derive(Clone, Debug, PartialEq)]
pub enum CinematicFinalizationPlanError {
    /// Cancellation was observed during expected-list reconstruction.
    Cancelled,
    /// A nonzero identity was required.
    MissingIdentity(&'static str),
    /// A declared resource ceiling was zero.
    InvalidLimit(&'static str),
    /// A declared plan-construction resource ceiling was exceeded.
    ResourceLimit {
        /// Stable resource name.
        resource: &'static str,
        /// Requested logical count.
        requested: u64,
        /// Admitted logical count.
        limit: u64,
    },
    /// A fallible plan-construction allocation was refused.
    Capacity(&'static str),
    /// Admitted inputs disagree semantically.
    Incompatible(&'static str),
    /// Checked integer conversion/arithmetic failed.
    ArithmeticOverflow(&'static str),
    /// Frame-sequence admission failed.
    Sequence(fs_img::FrameSequenceError),
}

/// One supplied image payload plus an externally pinned authority envelope.
#[derive(Clone, Copy, Debug)]
pub struct CinematicFrameArtifact<'a> {
    /// Canonical relative path named by the sequence manifest.
    pub relative_path: &'a str,
    /// Exact image bytes.
    pub bytes: &'a [u8],
    /// Canonical authority-record bytes.
    pub authority_bytes: &'a [u8],
    /// Authority identity retained outside those bytes.
    pub authority_identity: ContentHash,
}

/// Supplied audio payload and cross-clock receipts.
pub struct CinematicAudioArtifact<'a> {
    /// Exact persisted canonical binary snapshot of the audio manifest.  The
    /// verifier decodes these bytes rather than accepting an in-memory
    /// producer object as its oracle.
    pub manifest_bytes: &'a [u8],
    /// Externally pinned identity of `manifest_bytes`.
    pub manifest_identity: ContentHash,
    /// Exact authoritative float WAV bytes.
    pub wav_bytes: &'a [u8],
    /// WAV content identity retained independently from the manifest bytes.
    pub wav_identity: ContentHash,
    /// Canonical authority-record bytes for the WAV.
    pub authority_bytes: &'a [u8],
    /// Externally pinned audio-authority identity.
    pub authority_identity: ContentHash,
    /// Canonical persisted A/V marker receipt bytes.
    pub alignment_bytes: &'a [u8],
    /// Alignment receipt identity retained outside those bytes.
    pub alignment_identity: ContentHash,
    /// Canonical persisted resampled-event receipt bytes.
    pub event_bytes: &'a [u8],
    /// Event receipt identity retained outside those bytes.
    pub event_identity: ContentHash,
}

/// Owned canonical bytes plus an identity that must be persisted through an
/// independent channel. Each receipt family uses a distinct hash domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicCanonicalReceipt {
    bytes: Vec<u8>,
    identity: ContentHash,
}

impl CinematicCanonicalReceipt {
    /// Exact canonical receipt bytes.
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Domain-separated identity of the canonical bytes.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}

/// Refusal while constructing a canonical persisted receipt.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CinematicReceiptError {
    /// A caller ceiling was zero.
    InvalidLimit,
    /// The caller-owned record or byte ceiling was exceeded.
    BudgetExceeded,
    /// Cancellation was observed at a bounded record or hash boundary.
    Cancelled,
    /// The public collection cannot be represented by the u32 wire count.
    TooManyRecords,
    /// A source sample index cannot be represented by the u64 wire field.
    SourceIndexOverflow,
    /// Checked canonical byte-size arithmetic overflowed.
    SizeOverflow,
    /// The allocator refused the exact canonical receipt reservation.
    Capacity,
}

/// Canonically encode the exact A/V marker table.
pub fn encode_audio_video_alignment_receipt(
    alignment: &AudioVideoAlignment,
    max_markers: u32,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<CinematicCanonicalReceipt, CinematicReceiptError> {
    if max_markers == 0 || max_bytes == 0 {
        return Err(CinematicReceiptError::InvalidLimit);
    }
    cx.checkpoint()
        .map_err(|_| CinematicReceiptError::Cancelled)?;
    let count = u32::try_from(alignment.markers.len())
        .map_err(|_| CinematicReceiptError::TooManyRecords)?;
    if count > max_markers {
        return Err(CinematicReceiptError::BudgetExceeded);
    }
    let mut bytes = Vec::new();
    let encoded_len = 26_usize
        .checked_add(
            alignment
                .markers
                .len()
                .checked_mul(24)
                .ok_or(CinematicReceiptError::SizeOverflow)?,
        )
        .ok_or(CinematicReceiptError::SizeOverflow)?;
    if u64::try_from(encoded_len).map_err(|_| CinematicReceiptError::SizeOverflow)? > max_bytes {
        return Err(CinematicReceiptError::BudgetExceeded);
    }
    bytes
        .try_reserve_exact(encoded_len)
        .map_err(|_| CinematicReceiptError::Capacity)?;
    bytes.extend_from_slice(ALIGNMENT_RECEIPT_MAGIC);
    bytes.extend_from_slice(&RECEIPT_CODEC_VERSION.to_le_bytes());
    bytes.extend_from_slice(&alignment.audio_frames_per_video_frame.to_le_bytes());
    bytes.extend_from_slice(&alignment.endpoint_drift_audio_frames.to_le_bytes());
    bytes.extend_from_slice(&count.to_le_bytes());
    for (index, marker) in alignment.markers.iter().enumerate() {
        if index % 1_024 == 0 {
            cx.checkpoint()
                .map_err(|_| CinematicReceiptError::Cancelled)?;
        }
        bytes.extend_from_slice(&marker.video_tick.to_le_bytes());
        bytes.extend_from_slice(&marker.audio_tick.to_le_bytes());
        bytes.extend_from_slice(&marker.audio_frame_offset.to_le_bytes());
    }
    let identity = hash_domain_with_cancellation(ALIGNMENT_RECEIPT_DOMAIN, &bytes, cx)
        .ok_or(CinematicReceiptError::Cancelled)?;
    Ok(CinematicCanonicalReceipt { identity, bytes })
}

/// Canonically encode every exact resampled-event placement receipt.
pub fn encode_resampled_audio_event_receipt(
    events: &[ResampledAudioEvent],
    max_events: u32,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<CinematicCanonicalReceipt, CinematicReceiptError> {
    if max_events == 0 || max_bytes == 0 {
        return Err(CinematicReceiptError::InvalidLimit);
    }
    cx.checkpoint()
        .map_err(|_| CinematicReceiptError::Cancelled)?;
    let count = u32::try_from(events.len()).map_err(|_| CinematicReceiptError::TooManyRecords)?;
    if count > max_events {
        return Err(CinematicReceiptError::BudgetExceeded);
    }
    let mut encoded_len = 14_usize;
    for (index, event) in events.iter().enumerate() {
        if index % 1_024 == 0 {
            cx.checkpoint()
                .map_err(|_| CinematicReceiptError::Cancelled)?;
        }
        u64::try_from(event.source.source_sample_index)
            .map_err(|_| CinematicReceiptError::SourceIndexOverflow)?;
        let event_bytes = 109_usize
            .checked_add(event.source.artistic.map_or(0, |_| 56))
            .and_then(|bytes| bytes.checked_add(event.left_frame_offset.map_or(0, |_| 8)))
            .and_then(|bytes| bytes.checked_add(event.right_frame_offset.map_or(0, |_| 8)))
            .ok_or(CinematicReceiptError::SizeOverflow)?;
        encoded_len = encoded_len
            .checked_add(event_bytes)
            .ok_or(CinematicReceiptError::SizeOverflow)?;
    }
    let mut bytes = Vec::new();
    if u64::try_from(encoded_len).map_err(|_| CinematicReceiptError::SizeOverflow)? > max_bytes {
        return Err(CinematicReceiptError::BudgetExceeded);
    }
    bytes
        .try_reserve_exact(encoded_len)
        .map_err(|_| CinematicReceiptError::Capacity)?;
    bytes.extend_from_slice(EVENT_RECEIPT_MAGIC);
    bytes.extend_from_slice(&RECEIPT_CODEC_VERSION.to_le_bytes());
    bytes.extend_from_slice(&count.to_le_bytes());
    for (index, event) in events.iter().enumerate() {
        if index % 1_024 == 0 {
            cx.checkpoint()
                .map_err(|_| CinematicReceiptError::Cancelled)?;
        }
        push_resampled_event(&mut bytes, event);
    }
    let identity = hash_domain_with_cancellation(EVENT_RECEIPT_DOMAIN, &bytes, cx)
        .ok_or(CinematicReceiptError::Cancelled)?;
    Ok(CinematicCanonicalReceipt { identity, bytes })
}

/// Canonically snapshot every typed field required to verify an audio
/// artifact without retaining the producer's in-memory manifest object.
///
/// The external identity returned with the bytes is deliberately distinct
/// from [`AudioArtifactManifest::identity`]: the former authenticates this
/// persistence codec, while the latter is carried inside and independently
/// recomputed by the verifier.
pub fn encode_audio_artifact_manifest_receipt(
    manifest: &AudioArtifactManifest,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<CinematicCanonicalReceipt, CinematicReceiptError> {
    if max_bytes == 0 {
        return Err(CinematicReceiptError::InvalidLimit);
    }
    cx.checkpoint()
        .map_err(|_| CinematicReceiptError::Cancelled)?;

    let path = manifest.channel_layout().path();
    let encoded_len = 8_usize // magic
        + 2 // receipt codec version
        + 2 // audio-artifact schema version
        + 32 // original manifest identity
        + 2 // synthesis schema
        + 32 // synthesis configuration
        + 1 // synthesis authority
        + 32 * 4 // trajectory, excitation, model, timeline
        + 1 // manifest authority
        + 32 // channel-layout identity
        + 1 // signal path
        + usize::from(matches!(path, AudioSignalPath::SpatializedStereo { .. })) * 32
        + 32 // source-signal identity
        + 1 // optional mix tag
        + usize::from(manifest.mix_identity().is_some()) * 32
        + 32 * 2 // WAV and metadata identities
        + 8 * 2 // WAV byte and sample-frame counts
        + 4 // sample rate
        + 2 // sample encoding
        + 1 // artifact role
        + 8 * 5 // scalar meters
        + 1 // optional loudness tag
        + usize::from(manifest.meters().integrated_loudness_lufs.is_some()) * 8
        + 8 * 3 // meter block counts
        + 8 * 4 // video/audio clock endpoints
        + 4 // audio frames per video frame
        + 8; // admitted headroom
    if u64::try_from(encoded_len).map_err(|_| CinematicReceiptError::SizeOverflow)? > max_bytes {
        return Err(CinematicReceiptError::BudgetExceeded);
    }

    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(encoded_len)
        .map_err(|_| CinematicReceiptError::Capacity)?;
    bytes.extend_from_slice(AUDIO_MANIFEST_RECEIPT_MAGIC);
    bytes.extend_from_slice(&RECEIPT_CODEC_VERSION.to_le_bytes());
    bytes.extend_from_slice(&AUDIO_ARTIFACT_SCHEMA_VERSION.to_le_bytes());
    bytes.extend_from_slice(manifest.identity().as_bytes());
    let synthesis = manifest.synthesis();
    bytes.extend_from_slice(&synthesis.schema_version.to_le_bytes());
    bytes.extend_from_slice(synthesis.configuration_identity.as_bytes());
    bytes.push(sound_authority_tag(synthesis.authority));
    bytes.extend_from_slice(synthesis.trajectory_identity.as_bytes());
    bytes.extend_from_slice(synthesis.excitation_identity.as_bytes());
    bytes.extend_from_slice(synthesis.sound_model_identity.as_bytes());
    bytes.extend_from_slice(synthesis.timeline_identity.as_bytes());
    bytes.push(sound_authority_tag(manifest.authority()));
    bytes.extend_from_slice(manifest.channel_layout().identity().as_bytes());
    match path {
        AudioSignalPath::CanonicalDryStereo => bytes.push(1),
        AudioSignalPath::SpatializedStereo {
            spatialization_identity,
        } => {
            bytes.push(2);
            bytes.extend_from_slice(spatialization_identity.as_bytes());
        }
    }
    bytes.extend_from_slice(manifest.source_signal_identity().as_bytes());
    match manifest.mix_identity() {
        None => bytes.push(0),
        Some(identity) => {
            bytes.push(1);
            bytes.extend_from_slice(identity.as_bytes());
        }
    }
    let wav = manifest.wav();
    bytes.extend_from_slice(wav.wav_identity().as_bytes());
    bytes.extend_from_slice(wav.metadata_identity().as_bytes());
    bytes.extend_from_slice(&wav.byte_len().to_le_bytes());
    bytes.extend_from_slice(&wav.sample_frame_count().to_le_bytes());
    bytes.extend_from_slice(&wav.sample_rate_hz().to_le_bytes());
    bytes.extend_from_slice(&(wav.encoding() as u16).to_le_bytes());
    bytes.push(audio_artifact_role_tag(manifest.role()));
    push_audio_meters(&mut bytes, manifest.meters());
    let (video_start, video_end) = manifest.video_ticks();
    let (audio_start, audio_end) = manifest.audio_ticks();
    for tick in [video_start, video_end, audio_start, audio_end] {
        bytes.extend_from_slice(&tick.to_le_bytes());
    }
    bytes.extend_from_slice(&manifest.audio_frames_per_video_frame().to_le_bytes());
    bytes.extend_from_slice(&manifest.admitted_headroom_db().to_bits().to_le_bytes());
    debug_assert_eq!(bytes.len(), encoded_len);
    let identity = hash_domain_with_cancellation(AUDIO_MANIFEST_RECEIPT_DOMAIN, &bytes, cx)
        .ok_or(CinematicReceiptError::Cancelled)?;
    Ok(CinematicCanonicalReceipt { identity, bytes })
}

/// Complete read-only input bundle.
pub struct CinematicBundle<'a> {
    /// Canonical finalized sequence bytes.
    pub sequence_bytes: &'a [u8],
    /// Sequence identity retained independently from the bytes.
    pub sequence_identity: ContentHash,
    /// Exact image artifacts.  Input order is irrelevant; paths are unique.
    pub frames: &'a [CinematicFrameArtifact<'a>],
    /// Authoritative audio master and synchronization evidence.
    pub audio: CinematicAudioArtifact<'a>,
}

/// Stable top-level outcome.  `Pass` still does not imply physical truth or
/// aesthetic approval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CinematicFinalizationDisposition {
    /// Every expected byte-level and semantic check passed.
    Pass,
    /// Required persisted evidence is absent or the sequence is unfinished.
    Incomplete,
    /// A budget or cancellation boundary refused inspection.
    Refused,
    /// Bytes are malformed, noncanonical, truncated, or hash-inconsistent.
    Corrupt,
    /// Well-formed artifacts disagree with the independent plan.
    Incompatible,
}

impl CinematicFinalizationDisposition {
    const fn tag(self) -> u8 {
        match self {
            Self::Pass => 1,
            Self::Incomplete => 2,
            Self::Refused => 3,
            Self::Corrupt => 4,
            Self::Incompatible => 5,
        }
    }
}

/// Stable coordinate of the first divergence in verification order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CinematicFinalizationCoordinate {
    /// Plan or global budget admission.
    Bundle,
    /// Persisted frame-sequence snapshot.
    Sequence,
    /// One canonical frame artifact.
    Frame {
        /// Stable artifact key.
        key: FrameArtifactKey,
        /// Canonical relative path, when known.
        relative_path: String,
    },
    /// WAV or audio manifest.
    Audio,
    /// A/V marker index including the endpoint.
    SyncMarker(u32),
    /// Resampled event receipt index.
    AudioEvent(u32),
}

/// Stable defect class for automation and hostile-fixture assertions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum CinematicFinalizationDivergenceCode {
    /// Execution scope requested cancellation.
    Cancelled = 1,
    /// Aggregate verification byte/count ceiling was exceeded.
    BundleBudgetExceeded = 2,
    /// Frame-sequence bytes were not a valid canonical snapshot.
    SequenceDecode = 3,
    /// Frame sequence is structurally valid but unfinished.
    SequenceIncomplete = 4,
    /// Sequence rows disagree with the independent expected inventory.
    SequenceInventory = 5,
    /// A required artifact was not supplied.
    MissingArtifact = 6,
    /// The bundle supplied one artifact path more than once.
    DuplicateArtifact = 7,
    /// The bundle supplied a path absent from the expected inventory.
    UnexpectedArtifact = 8,
    /// Exact artifact bytes disagree with their independently retained hash.
    ArtifactHash = 9,
    /// Image codec structure is malformed or noncanonical.
    ImageStructure = 10,
    /// Image raster dimensions disagree with the plan.
    ImageDimensions = 11,
    /// Image channel names or sample types disagree with the plan.
    ImageChannels = 12,
    /// Image metadata disagrees with reconstructed renderer semantics.
    ImageMetadata = 13,
    /// Required authority bytes were absent.
    AuthorityMissing = 14,
    /// Authority bytes failed strict canonical decoding.
    AuthorityCodec = 15,
    /// Authority bytes disagree with their external identity pin.
    AuthorityIdentity = 16,
    /// Authority kind, role, transform, units, or disclosures disagree.
    AuthoritySemantics = 17,
    /// Composition/image configuration identity disagrees.
    ConfigurationIdentity = 18,
    /// Producer build identity disagrees.
    BuildIdentity = 19,
    /// Quality-profile identity disagrees.
    ProfileIdentity = 20,
    /// Source trajectory, scene, or derivation lineage disagrees.
    SourceIdentity = 21,
    /// Audio-manifest bytes disagree with an external or internal identity.
    AudioManifestIdentity = 22,
    /// Decoded audio-manifest roles, paths, clocks, or lineage disagree.
    AudioManifestSemantics = 23,
    /// WAV header, payload, meters, or independent content pin disagrees.
    WavStructure = 24,
    /// Exact audio/video master durations disagree.
    AudioVideoDuration = 25,
    /// One exact cross-clock synchronization marker disagrees.
    SyncMarker = 26,
    /// One exact shot-cut boundary marker disagrees.
    CutMarker = 27,
    /// One ordered resampled audio event receipt disagrees.
    AudioEvent = 28,
    /// Raw EXR per-pixel sample-count evidence disagrees.
    ImageSampleCount = 29,
    /// Raw EXR pixel payload is non-finite or semantically invalid.
    ImagePayload = 30,
}

impl CinematicFinalizationDivergenceCode {
    const fn tag(self) -> u8 {
        self as u8
    }

    /// Stable machine-readable code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Cancelled => "finalize-cancelled",
            Self::BundleBudgetExceeded => "finalize-bundle-budget",
            Self::SequenceDecode => "finalize-sequence-decode",
            Self::SequenceIncomplete => "finalize-sequence-incomplete",
            Self::SequenceInventory => "finalize-sequence-inventory",
            Self::MissingArtifact => "finalize-missing-artifact",
            Self::DuplicateArtifact => "finalize-duplicate-artifact",
            Self::UnexpectedArtifact => "finalize-unexpected-artifact",
            Self::ArtifactHash => "finalize-artifact-hash",
            Self::ImageStructure => "finalize-image-structure",
            Self::ImageDimensions => "finalize-image-dimensions",
            Self::ImageChannels => "finalize-image-channels",
            Self::ImageMetadata => "finalize-image-metadata",
            Self::AuthorityMissing => "finalize-authority-missing",
            Self::AuthorityCodec => "finalize-authority-codec",
            Self::AuthorityIdentity => "finalize-authority-identity",
            Self::AuthoritySemantics => "finalize-authority-semantics",
            Self::ConfigurationIdentity => "finalize-configuration-identity",
            Self::BuildIdentity => "finalize-build-identity",
            Self::ProfileIdentity => "finalize-profile-identity",
            Self::SourceIdentity => "finalize-source-identity",
            Self::AudioManifestIdentity => "finalize-audio-manifest-identity",
            Self::AudioManifestSemantics => "finalize-audio-manifest-semantics",
            Self::WavStructure => "finalize-wav-structure",
            Self::AudioVideoDuration => "finalize-av-duration",
            Self::SyncMarker => "finalize-sync-marker",
            Self::CutMarker => "finalize-cut-marker",
            Self::AudioEvent => "finalize-audio-event",
            Self::ImageSampleCount => "finalize-image-sample-count",
            Self::ImagePayload => "finalize-image-payload",
        }
    }
}

/// First observed mismatch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicFinalizationDivergence {
    /// Stable defect class.
    pub code: CinematicFinalizationDivergenceCode,
    /// Exact verification coordinate.
    pub coordinate: CinematicFinalizationCoordinate,
}

/// Deterministically ranked repair advice.  Repairs never mutate artifacts.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[repr(u8)]
pub enum CinematicFinalizationRepair {
    /// Resume the producer from retained state or finish missing work.
    CompleteOrResumeProduction = 1,
    /// Restore the artifact named by the expected inventory.
    RestoreExpectedArtifact = 2,
    /// Re-run the producer from the exact pinned inputs.
    RegenerateArtifactFromPinnedInputs = 3,
    /// Restore the canonical manifest/receipt bytes.
    RestoreCanonicalManifest = 4,
    /// Supply the independently retained pin for the exact bytes.
    SupplyCorrectExternalIdentityPin = 5,
    /// Rebuild using the configuration admitted by the plan.
    RebuildWithAdmittedConfiguration = 6,
    /// Increase an explicit verifier ceiling and retry.
    IncreaseExplicitVerificationBudget = 7,
    /// Retry under a live, non-cancelled execution scope.
    RetryInLiveExecutionScope = 8,
}

impl CinematicFinalizationRepair {
    const fn tag(self) -> u8 {
        self as u8
    }
}

/// Deterministic finalization report.  The identity binds child pins, outcome,
/// first divergence, repairs, and the union of verified no-claims.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CinematicFinalizationReport {
    identity: ContentHash,
    plan_identity: ContentHash,
    sequence_identity: ContentHash,
    audio_manifest_identity: ContentHash,
    wav_identity: ContentHash,
    alignment_identity: ContentHash,
    event_identity: ContentHash,
    target: CinematicFinalizationTarget,
    disposition: CinematicFinalizationDisposition,
    first_divergence: Option<CinematicFinalizationDivergence>,
    repairs: Vec<CinematicFinalizationRepair>,
    verified_frame_artifacts: u32,
    verified_sync_markers: u32,
    verified_audio_events: u32,
    no_claims: Vec<CinematicNoClaim>,
}

impl CinematicFinalizationReport {
    /// Delivery-authority class inherited from the independent plan.
    #[must_use]
    pub const fn target(&self) -> CinematicFinalizationTarget {
        self.target
    }

    /// Overall stable outcome.
    #[must_use]
    pub const fn disposition(&self) -> CinematicFinalizationDisposition {
        self.disposition
    }

    /// Earliest deterministic mismatch, if any.
    #[must_use]
    pub const fn first_divergence(&self) -> Option<&CinematicFinalizationDivergence> {
        self.first_divergence.as_ref()
    }

    /// Ranked repair list.
    #[must_use]
    pub fn repairs(&self) -> &[CinematicFinalizationRepair] {
        &self.repairs
    }

    /// Union of disclosures verified on child authority records.
    #[must_use]
    pub fn no_claims(&self) -> &[CinematicNoClaim] {
        &self.no_claims
    }

    /// Only a passing complete `Final4k` plan is delivery eligible.
    #[must_use]
    pub const fn final_delivery_eligible(&self) -> bool {
        matches!(self.target, CinematicFinalizationTarget::Final4k)
            && matches!(self.disposition, CinematicFinalizationDisposition::Pass)
    }

    /// Domain-separated report identity.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Canonical compact report preimage.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        report_bytes(
            self.plan_identity,
            self.sequence_identity,
            self.audio_manifest_identity,
            self.wav_identity,
            self.alignment_identity,
            self.event_identity,
            self.target,
            self.disposition,
            self.first_divergence.as_ref(),
            &self.repairs,
            self.verified_frame_artifacts,
            self.verified_sync_markers,
            self.verified_audio_events,
            &self.no_claims,
        )
    }
}

struct ReportBuilder<'a> {
    plan: &'a CinematicFinalizationPlan,
    sequence_identity: ContentHash,
    audio_manifest_identity: ContentHash,
    wav_identity: ContentHash,
    alignment_identity: ContentHash,
    event_identity: ContentHash,
    verified_frames: u32,
    verified_sync_markers: u32,
    verified_audio_events: u32,
    no_claims: BTreeSet<CinematicNoClaim>,
}

impl<'a> ReportBuilder<'a> {
    fn fail(
        &self,
        disposition: CinematicFinalizationDisposition,
        code: CinematicFinalizationDivergenceCode,
        coordinate: CinematicFinalizationCoordinate,
    ) -> CinematicFinalizationReport {
        let divergence = CinematicFinalizationDivergence { code, coordinate };
        self.finish(
            disposition,
            Some(divergence),
            repairs_for(disposition, code),
        )
    }

    fn pass(&self) -> CinematicFinalizationReport {
        self.finish(CinematicFinalizationDisposition::Pass, None, Vec::new())
    }

    fn finish(
        &self,
        disposition: CinematicFinalizationDisposition,
        first_divergence: Option<CinematicFinalizationDivergence>,
        repairs: Vec<CinematicFinalizationRepair>,
    ) -> CinematicFinalizationReport {
        let no_claims: Vec<_> = self.no_claims.iter().copied().collect();
        let bytes = report_bytes(
            self.plan.identity,
            self.sequence_identity,
            self.audio_manifest_identity,
            self.wav_identity,
            self.alignment_identity,
            self.event_identity,
            self.plan.target,
            disposition,
            first_divergence.as_ref(),
            &repairs,
            self.verified_frames,
            self.verified_sync_markers,
            self.verified_audio_events,
            &no_claims,
        );
        CinematicFinalizationReport {
            identity: hash_domain(CINEMATIC_FINALIZATION_REPORT_DOMAIN, &bytes),
            plan_identity: self.plan.identity,
            sequence_identity: self.sequence_identity,
            audio_manifest_identity: self.audio_manifest_identity,
            wav_identity: self.wav_identity,
            alignment_identity: self.alignment_identity,
            event_identity: self.event_identity,
            target: self.plan.target,
            disposition,
            first_divergence,
            repairs,
            verified_frame_artifacts: self.verified_frames,
            verified_sync_markers: self.verified_sync_markers,
            verified_audio_events: self.verified_audio_events,
            no_claims,
        }
    }
}

/// Independently verify all persisted frame/audio bytes and synchronization
/// receipts.  No supplied object is mutated and no filesystem path is opened.
pub fn verify_cinematic_bundle(
    plan: &CinematicFinalizationPlan,
    bundle: &CinematicBundle<'_>,
    limits: CinematicFinalizationLimits,
    cx: &Cx<'_>,
) -> CinematicFinalizationReport {
    let mut report = ReportBuilder {
        plan,
        sequence_identity: bundle.sequence_identity,
        audio_manifest_identity: bundle.manifest_identity(),
        wav_identity: bundle.audio.wav_identity,
        alignment_identity: bundle.audio.alignment_identity,
        event_identity: bundle.audio.event_identity,
        verified_frames: 0,
        verified_sync_markers: 0,
        verified_audio_events: 0,
        no_claims: BTreeSet::new(),
    };
    if limits.validate().is_err() {
        return report.fail(
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
            CinematicFinalizationCoordinate::Bundle,
        );
    }
    if cx.checkpoint().is_err() {
        return report.fail(
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
            CinematicFinalizationCoordinate::Bundle,
        );
    }
    if bundle.frames.len() > plan.expected_sequence.entries().len() {
        return report.fail(
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::UnexpectedArtifact,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    let bundle_bytes = match bundle_size(bundle, cx) {
        Ok(bytes) => bytes,
        Err(BundleSizeError::Overflow) => u64::MAX,
        Err(BundleSizeError::Cancelled) => {
            return report.fail(
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::Cancelled,
                CinematicFinalizationCoordinate::Bundle,
            );
        }
    };
    if bundle_bytes > limits.max_bundle_bytes
        || u32::try_from(plan.expected_audio_events.len())
            .map_or(true, |count| count > limits.max_audio_events)
    {
        return report.fail(
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
            CinematicFinalizationCoordinate::Bundle,
        );
    }
    if bundle.sequence_bytes.is_empty() {
        return report.fail(
            CinematicFinalizationDisposition::Incomplete,
            CinematicFinalizationDivergenceCode::SequenceIncomplete,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    let actual_sequence = match FrameSequenceManifest::decode_snapshot_with_poll(
        bundle.sequence_bytes,
        bundle.sequence_identity,
        limits.frame_sequence,
        limits.frame_sequence.max_output_bytes(),
        || cx.checkpoint().is_ok(),
    ) {
        Ok(manifest) => manifest,
        Err(error) => {
            let (disposition, code) = map_sequence_decode_error(&error);
            return report.fail(disposition, code, CinematicFinalizationCoordinate::Sequence);
        }
    };
    if cx.checkpoint().is_err() {
        return report.fail(
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    if actual_sequence.state() != FrameSequenceState::Finalized {
        return report.fail(
            CinematicFinalizationDisposition::Incomplete,
            CinematicFinalizationDivergenceCode::SequenceIncomplete,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    let expected_context = plan.expected_sequence.context();
    let actual_context = actual_sequence.context();
    if actual_context.build_id() != plan.build_identity {
        return report.fail(
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::BuildIdentity,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    if actual_context.render_config_id() != plan.image_configuration_identity {
        return report.fail(
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::ConfigurationIdentity,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    if actual_context.profile_id() != plan.profile_identity {
        return report.fail(
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::ProfileIdentity,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    if actual_context != expected_context
        || actual_sequence.limits() != plan.expected_sequence.limits()
    {
        return report.fail(
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::SourceIdentity,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    if actual_sequence.entries().len() != plan.expected_sequence.entries().len() {
        return report.fail(
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::SequenceInventory,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    for (index, (actual, expected)) in actual_sequence
        .entries()
        .iter()
        .zip(plan.expected_sequence.entries())
        .enumerate()
    {
        if index.is_multiple_of(256) && cx.checkpoint().is_err() {
            return report.fail(
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::Cancelled,
                CinematicFinalizationCoordinate::Sequence,
            );
        }
        if actual.relative_path() != expected.relative_path()
            || actual.descriptor() != expected.descriptor()
            || actual.max_bytes() != expected.max_bytes()
            || actual.source() != expected.source()
        {
            return report.fail(
                CinematicFinalizationDisposition::Incompatible,
                CinematicFinalizationDivergenceCode::SequenceInventory,
                CinematicFinalizationCoordinate::Sequence,
            );
        }
    }

    // Bundle order is intentionally irrelevant, but verifier-owned indexing
    // must remain fallible and bounded.  Sorting borrowed rows avoids one
    // allocation per BTree node and never duplicates artifact payload bytes.
    let mut artifacts = Vec::new();
    if artifacts.try_reserve_exact(bundle.frames.len()).is_err() {
        return report.fail(
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    for (index, artifact) in bundle.frames.iter().enumerate() {
        if index % 256 == 0 && cx.checkpoint().is_err() {
            return report.fail(
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::Cancelled,
                CinematicFinalizationCoordinate::Sequence,
            );
        }
        artifacts.push(*artifact);
    }
    if cx.checkpoint().is_err() {
        return report.fail(
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    artifacts.sort_unstable_by(|left, right| left.relative_path.cmp(right.relative_path));
    if cx.checkpoint().is_err() {
        return report.fail(
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    for (index, pair) in artifacts.windows(2).enumerate() {
        if index.is_multiple_of(256) && cx.checkpoint().is_err() {
            return report.fail(
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::Cancelled,
                CinematicFinalizationCoordinate::Sequence,
            );
        }
        if pair[0].relative_path == pair[1].relative_path {
            return report.fail(
                CinematicFinalizationDisposition::Corrupt,
                CinematicFinalizationDivergenceCode::DuplicateArtifact,
                CinematicFinalizationCoordinate::Sequence,
            );
        }
    }
    // The sequence's canonical order is by artifact key rather than path, so
    // use a separate fallible borrowed-path index to classify an unexpected
    // artifact before the corresponding missing expected artifact.
    let mut expected_paths = Vec::new();
    if expected_paths
        .try_reserve_exact(actual_sequence.entries().len())
        .is_err()
    {
        return report.fail(
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
            CinematicFinalizationCoordinate::Sequence,
        );
    }
    for (index, entry) in actual_sequence.entries().iter().enumerate() {
        if index.is_multiple_of(256) && cx.checkpoint().is_err() {
            return report.fail(
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::Cancelled,
                CinematicFinalizationCoordinate::Sequence,
            );
        }
        expected_paths.push(entry.relative_path());
    }
    expected_paths.sort_unstable();
    for (index, artifact) in artifacts.iter().enumerate() {
        if index.is_multiple_of(256) && cx.checkpoint().is_err() {
            return report.fail(
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::Cancelled,
                CinematicFinalizationCoordinate::Sequence,
            );
        }
        if expected_paths
            .binary_search(&artifact.relative_path)
            .is_err()
        {
            return report.fail(
                CinematicFinalizationDisposition::Corrupt,
                CinematicFinalizationDivergenceCode::UnexpectedArtifact,
                CinematicFinalizationCoordinate::Sequence,
            );
        }
    }
    for entry in actual_sequence.entries() {
        let coordinate = || CinematicFinalizationCoordinate::Frame {
            key: entry.descriptor().key(),
            relative_path: entry.relative_path().to_owned(),
        };
        let Ok(artifact_index) = artifacts
            .binary_search_by(|artifact| artifact.relative_path.cmp(entry.relative_path()))
        else {
            return report.fail(
                CinematicFinalizationDisposition::Incomplete,
                CinematicFinalizationDivergenceCode::MissingArtifact,
                coordinate(),
            );
        };
        let artifact = artifacts[artifact_index];
        let Some(file_state) = entry.file_state() else {
            return report.fail(
                CinematicFinalizationDisposition::Incomplete,
                CinematicFinalizationDivergenceCode::SequenceIncomplete,
                coordinate(),
            );
        };
        let actual_hash = match hash_with_cancellation(artifact.bytes, cx) {
            Some(hash) => hash,
            None => {
                return report.fail(
                    CinematicFinalizationDisposition::Refused,
                    CinematicFinalizationDivergenceCode::Cancelled,
                    coordinate(),
                );
            }
        };
        if actual_hash != file_state.content_hash()
            || u64::try_from(artifact.bytes.len()).ok() != Some(file_state.byte_size())
        {
            return report.fail(
                CinematicFinalizationDisposition::Corrupt,
                CinematicFinalizationDivergenceCode::ArtifactHash,
                coordinate(),
            );
        }
        let inspection = match inspect_frame(entry.descriptor(), artifact.bytes, limits, cx) {
            Ok(inspection) => inspection,
            Err((disposition, code)) => return report.fail(disposition, code, coordinate()),
        };
        let authority = match decode_authority(
            artifact.authority_bytes,
            artifact.authority_identity,
            limits.max_authority_record_bytes,
            cx,
        ) {
            Ok(record) => record,
            Err((disposition, code)) => return report.fail(disposition, code, coordinate()),
        };
        if let Err(code) = verify_frame_authority(plan, entry, actual_hash, &authority) {
            return report.fail(
                CinematicFinalizationDisposition::Incompatible,
                code,
                coordinate(),
            );
        }
        if entry.descriptor().key().role() == FrameArtifactRole::RawMaster {
            let expected_raw = match plan.raw_frames.get(&entry.descriptor().key()) {
                Some(expected) => expected,
                None => {
                    return report.fail(
                        CinematicFinalizationDisposition::Incompatible,
                        CinematicFinalizationDivergenceCode::SequenceInventory,
                        coordinate(),
                    );
                }
            };
            if let Err(code) = verify_raw_metadata(&inspection, expected_raw) {
                return report.fail(
                    CinematicFinalizationDisposition::Incompatible,
                    code,
                    coordinate(),
                );
            }
            if let Err((disposition, code)) =
                verify_raw_payload(artifact.bytes, expected_raw, limits.exr, cx)
            {
                return report.fail(disposition, code, coordinate());
            }
            if let Some(expected_spp) = expected_raw.expected_uniform_spp
                && let Err((disposition, code)) =
                    verify_raw_sample_count(artifact.bytes, expected_spp, limits.exr, cx)
            {
                return report.fail(disposition, code, coordinate());
            }
        } else if entry.descriptor().format() == FrameArtifactFormat::OpenExr {
            let Some(source_hash) = entry.source_content_hash() else {
                return report.fail(
                    CinematicFinalizationDisposition::Incompatible,
                    CinematicFinalizationDivergenceCode::SourceIdentity,
                    coordinate(),
                );
            };
            if let Err(code) = verify_derived_source(&inspection, source_hash) {
                return report.fail(
                    CinematicFinalizationDisposition::Incompatible,
                    code,
                    coordinate(),
                );
            }
        }
        report.no_claims.extend(authority.no_claims());
        report.verified_frames = report
            .verified_frames
            .checked_add(1)
            .expect("decoded manifest limits frame-artifact count to u32");
    }

    let audio = &bundle.audio;
    let audio_manifest = match decode_audio_manifest_receipt(
        audio.manifest_bytes,
        audio.manifest_identity,
        limits.max_audio_manifest_bytes,
        cx,
    ) {
        Ok(manifest) => manifest,
        Err(error) => {
            let (disposition, code) = map_receipt_decode_error(
                error,
                CinematicFinalizationDivergenceCode::AudioManifestIdentity,
            );
            return report.fail(disposition, code, CinematicFinalizationCoordinate::Audio);
        }
    };
    if audio_manifest.wav.wav_identity != audio.wav_identity {
        return report.fail(
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::WavStructure,
            CinematicFinalizationCoordinate::Audio,
        );
    }
    if audio_manifest.synthesis != plan.sound_receipt
        || audio_manifest.role != AudioArtifactRole::AuthoritativeFloat32Master
        || audio_manifest.wav.encoding != WavSampleEncoding::Float32
        || audio_manifest.wav.sample_rate_hz != plan.audio_sample_rate_hz
        || audio_manifest.wav.sample_frame_count != plan.total_audio_sample_frames
        || (
            audio_manifest.video_start_tick,
            audio_manifest.video_end_tick_exclusive,
        ) != (0, i64::from(plan.total_video_frames))
        || (
            audio_manifest.audio_start_tick,
            audio_manifest.audio_end_tick_exclusive,
        ) != (
            0,
            i64::try_from(plan.total_audio_sample_frames).unwrap_or(i64::MAX),
        )
        || audio_manifest.audio_frames_per_video_frame
            != plan.audio_sample_rate_hz / plan.frames_per_second
        || audio_manifest.source_signal_identity != plan.expected_audio_source_signal_identity
        || audio_manifest.channel_layout_identity != plan.expected_audio_channel_layout_identity
        || audio_manifest.mix_identity != plan.expected_audio_mix_identity
    {
        return report.fail(
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::AudioManifestSemantics,
            CinematicFinalizationCoordinate::Audio,
        );
    }
    if audio.wav_bytes.is_empty() {
        return report.fail(
            CinematicFinalizationDisposition::Incomplete,
            CinematicFinalizationDivergenceCode::WavStructure,
            CinematicFinalizationCoordinate::Audio,
        );
    }
    if let Err((disposition, code)) =
        verify_wav_against_snapshot(&audio_manifest, audio.wav_bytes, limits.audio, cx)
    {
        return report.fail(disposition, code, CinematicFinalizationCoordinate::Audio);
    }
    let audio_authority = match decode_authority(
        audio.authority_bytes,
        audio.authority_identity,
        limits.max_authority_record_bytes,
        cx,
    ) {
        Ok(record) => record,
        Err((disposition, code)) => {
            return report.fail(disposition, code, CinematicFinalizationCoordinate::Audio);
        }
    };
    if let Err(code) = verify_audio_authority(plan, &audio_manifest, &audio_authority) {
        return report.fail(
            CinematicFinalizationDisposition::Incompatible,
            code,
            CinematicFinalizationCoordinate::Audio,
        );
    }
    report.no_claims.extend(audio_authority.no_claims());
    let alignment = match decode_alignment_receipt(
        audio.alignment_bytes,
        audio.alignment_identity,
        limits.max_sync_markers,
        cx,
    ) {
        Ok(alignment) => alignment,
        Err(error) => {
            let (disposition, code) =
                map_receipt_decode_error(error, CinematicFinalizationDivergenceCode::SyncMarker);
            return report.fail(disposition, code, CinematicFinalizationCoordinate::Audio);
        }
    };
    match verify_alignment(plan, &alignment, cx) {
        Err(()) => {
            return report.fail(
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::Cancelled,
                CinematicFinalizationCoordinate::Audio,
            );
        }
        Ok(Some((index, cut))) => {
            return report.fail(
                CinematicFinalizationDisposition::Incompatible,
                if cut {
                    CinematicFinalizationDivergenceCode::CutMarker
                } else {
                    CinematicFinalizationDivergenceCode::SyncMarker
                },
                CinematicFinalizationCoordinate::SyncMarker(index),
            );
        }
        Ok(None) => {}
    }
    report.verified_sync_markers = u32::try_from(alignment.markers.len())
        .expect("alignment length was checked against admitted u32 marker count");
    let events = match decode_event_receipt(
        audio.event_bytes,
        audio.event_identity,
        limits.max_audio_events,
        cx,
    ) {
        Ok(events) => events,
        Err(error) => {
            let (disposition, code) =
                map_receipt_decode_error(error, CinematicFinalizationDivergenceCode::AudioEvent);
            return report.fail(disposition, code, CinematicFinalizationCoordinate::Audio);
        }
    };
    if events.len() != plan.expected_audio_events.len() {
        return report.fail(
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::AudioEvent,
            CinematicFinalizationCoordinate::AudioEvent(
                u32::try_from(events.len().min(plan.expected_audio_events.len()))
                    .expect("receipt and plan event counts were admitted as u32"),
            ),
        );
    }
    let mut previous_event_time_s = None;
    for (index, (actual, expected)) in events.iter().zip(&plan.expected_audio_events).enumerate() {
        if index % 1_024 == 0 && cx.checkpoint().is_err() {
            return report.fail(
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::Cancelled,
                CinematicFinalizationCoordinate::AudioEvent(
                    u32::try_from(index).expect("plan admission bounds event count to u32"),
                ),
            );
        }
        if validate_resampled_audio_event(actual, plan.total_audio_sample_frames).is_err()
            || previous_event_time_s.is_some_and(|previous| actual.source.time_s <= previous)
            || !resampled_event_eq(actual, expected)
        {
            return report.fail(
                CinematicFinalizationDisposition::Incompatible,
                CinematicFinalizationDivergenceCode::AudioEvent,
                CinematicFinalizationCoordinate::AudioEvent(
                    u32::try_from(index).unwrap_or(u32::MAX),
                ),
            );
        }
        previous_event_time_s = Some(actual.source.time_s);
        report.verified_audio_events = report
            .verified_audio_events
            .checked_add(1)
            .expect("plan admission limits event count to u32");
    }
    report.pass()
}

impl CinematicBundle<'_> {
    const fn manifest_identity(&self) -> ContentHash {
        self.audio.manifest_identity
    }
}

fn checkpoint_plan(cx: &Cx<'_>) -> Result<(), CinematicFinalizationPlanError> {
    cx.checkpoint()
        .map_err(|_| CinematicFinalizationPlanError::Cancelled)
}

fn map_sequence_plan_error(error: FrameSequenceError) -> CinematicFinalizationPlanError {
    match error {
        FrameSequenceError::Cancelled => CinematicFinalizationPlanError::Cancelled,
        FrameSequenceError::ResourceLimit {
            resource,
            requested,
            limit,
        } => CinematicFinalizationPlanError::ResourceLimit {
            resource,
            requested,
            limit,
        },
        FrameSequenceError::AllocationRefused { resource, .. } => {
            CinematicFinalizationPlanError::Capacity(resource)
        }
        FrameSequenceError::SizeOverflow { context } => {
            CinematicFinalizationPlanError::ArithmeticOverflow(context)
        }
        error => CinematicFinalizationPlanError::Sequence(error),
    }
}

fn is_zero(identity: ContentHash) -> bool {
    identity.as_bytes().iter().all(|byte| *byte == 0)
}

fn index_frame_inputs<'frame>(
    inputs: &[EulerRenderFrameInput<'frame>],
    plan: &EulerUniformRenderPlan,
    cx: &Cx<'_>,
) -> Result<BTreeMap<u64, EulerRenderFrameInput<'frame>>, CinematicFinalizationPlanError> {
    let mut indexed = BTreeMap::new();
    for input in inputs {
        checkpoint_plan(cx)?;
        if indexed.insert(input.frame_ordinal(), *input).is_some() {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "duplicate prepared frame ordinal",
            ));
        }
    }
    if indexed.len() != plan.frames().len() {
        return Err(CinematicFinalizationPlanError::Incompatible(
            "prepared frame count",
        ));
    }
    for frame in plan.frames() {
        let input = indexed.get(&frame.frame_ordinal()).ok_or(
            CinematicFinalizationPlanError::Incompatible("missing prepared frame"),
        )?;
        if u64::try_from(input.prepared().segments().len()).ok() != Some(frame.segment_count()) {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "prepared frame segment count",
            ));
        }
    }
    Ok(indexed)
}

fn validate_prepared_frames_against_brief(
    inputs: &BTreeMap<u64, EulerRenderFrameInput<'_>>,
    plan: &EulerUniformRenderPlan,
    brief: &CinematicBrief,
    frames_per_second: u32,
    cx: &Cx<'_>,
) -> Result<(), CinematicFinalizationPlanError> {
    let cut_starts: BTreeSet<_> = brief
        .shots()
        .iter()
        .skip(1)
        .map(|shot| shot.frames().start())
        .collect();
    for frame in plan.frames() {
        checkpoint_plan(cx)?;
        let frame_ordinal = u32::try_from(frame.frame_ordinal()).map_err(|_| {
            CinematicFinalizationPlanError::ArithmeticOverflow("brief frame ordinal")
        })?;
        let prepared = inputs
            .get(&frame.frame_ordinal())
            .ok_or(CinematicFinalizationPlanError::Incompatible(
                "prepared frame inventory",
            ))?
            .prepared();
        let segments = prepared.segments();
        let (Some(first), Some(last)) = (segments.first(), segments.last()) else {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "prepared shutter inventory",
            ));
        };
        let window = brief
            .effective_shutter_window(frame_ordinal)
            .map_err(|_| CinematicFinalizationPlanError::Incompatible("brief shutter window"))?;
        let scale = 1_000_000.0 * f64::from(frames_per_second);
        let expected_open_s = window.start_microframes as f64 / scale;
        let expected_close_s = window.end_microframes as f64 / scale;
        let mut duration_weight_sum = 0.0;
        let mut previous_close_s = None;
        let mut segment_discontinuity = false;
        for (index, segment) in segments.iter().enumerate() {
            if index.is_multiple_of(1_024) {
                checkpoint_plan(cx)?;
            }
            if previous_close_s
                .is_some_and(|close_s| !same_clock_second(close_s, segment.shutter().open_s()))
            {
                segment_discontinuity = true;
            }
            previous_close_s = Some(segment.shutter().close_s());
            duration_weight_sum += segment.duration_weight();
        }
        if !same_clock_second(first.shutter().open_s(), expected_open_s)
            || !same_clock_second(last.shutter().close_s(), expected_close_s)
            || segment_discontinuity
            || !same_clock_second(duration_weight_sum, 1.0)
        {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "prepared shutter and brief",
            ));
        }
        if cut_starts.contains(&frame_ordinal) && prepared.cut_side() != CutSide::After {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "prepared cut ownership",
            ));
        }
    }
    Ok(())
}

fn same_clock_second(actual: f64, expected: f64) -> bool {
    actual.is_finite()
        && expected.is_finite()
        && (actual - expected).abs()
            <= 8.0 * f64::EPSILON * actual.abs().max(expected.abs()).max(1.0)
}

fn validate_raw_expectations(
    sequence: &FrameSequenceManifest,
    raw_frames: &BTreeMap<FrameArtifactKey, RawFrameExpectation>,
    cx: &Cx<'_>,
) -> Result<(), CinematicFinalizationPlanError> {
    let mut expected_raw_keys = BTreeSet::new();
    for (index, entry) in sequence.entries().iter().enumerate() {
        if index.is_multiple_of(1_024) {
            checkpoint_plan(cx)?;
        }
        if entry.descriptor().key().role() == FrameArtifactRole::RawMaster {
            expected_raw_keys.insert(entry.descriptor().key());
        }
    }
    if raw_frames.len() != expected_raw_keys.len()
        || raw_frames
            .keys()
            .any(|key| !expected_raw_keys.contains(key))
    {
        return Err(CinematicFinalizationPlanError::Incompatible(
            "raw frame expectations",
        ));
    }
    for (index, expectation) in raw_frames.values().enumerate() {
        if index.is_multiple_of(1_024) {
            checkpoint_plan(cx)?;
        }
        if is_zero(expectation.authority_source_identity)
            || is_zero(expectation.authority_transform_identity)
        {
            return Err(CinematicFinalizationPlanError::MissingIdentity(
                "raw frame authority identity",
            ));
        }
        if expectation.expected_uniform_spp == Some(0) {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "raw sample-count expectation",
            ));
        }
        if expectation.object_palette_entries >= MAX_EXACT_AOV_PALETTE_INDEX
            || expectation.material_palette_entries >= MAX_EXACT_AOV_PALETTE_INDEX
        {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "raw palette cardinality expectation",
            ));
        }
        u32::try_from(expectation.required_attributes.len()).map_err(|_| {
            CinematicFinalizationPlanError::ArithmeticOverflow("raw attribute count")
        })?;
        for (attribute_index, (name, value)) in expectation.required_attributes.iter().enumerate() {
            if attribute_index.is_multiple_of(64) {
                checkpoint_plan(cx)?;
            }
            u32::try_from(name.len()).map_err(|_| {
                CinematicFinalizationPlanError::ArithmeticOverflow("raw attribute name")
            })?;
            u32::try_from(value.len()).map_err(|_| {
                CinematicFinalizationPlanError::ArithmeticOverflow("raw attribute value")
            })?;
        }
    }
    Ok(())
}

const fn aov_profile(preset: AovPreset) -> CinematicAovProfile {
    match preset {
        AovPreset::BeautyXyz => CinematicAovProfile::BeautyOnly,
        AovPreset::DailyCore => CinematicAovProfile::DailyCore,
        AovPreset::FinalDiagnostic => CinematicAovProfile::FinalDiagnostic,
    }
}

#[allow(clippy::too_many_arguments)]
fn preflight_expected_inventory(
    segment_count: usize,
    aov_profile: CinematicAovProfile,
    denoise_policy: DenoisePolicy,
    width: u32,
    height: u32,
    profile_output_ceiling: u64,
    artifact_ceilings: CinematicFrameArtifactCeilings,
    sequence_limits: FrameSequenceLimits,
    aov_limits: CinematicAovLimits,
) -> Result<usize, CinematicFinalizationPlanError> {
    let segment_count = u64::try_from(segment_count)
        .map_err(|_| CinematicFinalizationPlanError::ArithmeticOverflow("render segment count"))?;
    if segment_count == 0 {
        return Err(CinematicFinalizationPlanError::Incompatible(
            "empty render segment inventory",
        ));
    }
    let artifacts_per_segment = if denoise_policy == DenoisePolicy::SeparateBiasedDerivative {
        3_u64
    } else {
        2_u64
    };
    let artifact_count = segment_count.checked_mul(artifacts_per_segment).ok_or(
        CinematicFinalizationPlanError::ArithmeticOverflow("expected frame artifact count"),
    )?;
    if artifact_count > u64::from(sequence_limits.max_artifacts()) {
        return Err(CinematicFinalizationPlanError::ResourceLimit {
            resource: "expected frame artifact count",
            requested: artifact_count,
            limit: u64::from(sequence_limits.max_artifacts()),
        });
    }
    let raw_channels = u64::try_from(aov_profile.exr_channel_layout().len())
        .map_err(|_| CinematicFinalizationPlanError::ArithmeticOverflow("raw AOV channel count"))?;
    if raw_channels > u64::from(sequence_limits.max_channels_per_artifact()) {
        return Err(CinematicFinalizationPlanError::ResourceLimit {
            resource: "raw AOV channel count",
            requested: raw_channels,
            limit: u64::from(sequence_limits.max_channels_per_artifact()),
        });
    }
    let pixels = u64::from(width).checked_mul(u64::from(height)).ok_or(
        CinematicFinalizationPlanError::ArithmeticOverflow("AOV raster pixels"),
    )?;
    if pixels > aov_limits.max_pixels() {
        return Err(CinematicFinalizationPlanError::ResourceLimit {
            resource: "AOV pixels",
            requested: pixels,
            limit: aov_limits.max_pixels(),
        });
    }
    let mut reserved_per_segment = artifact_ceilings
        .raw_master_bytes
        .checked_add(artifact_ceilings.display_preview_bytes)
        .ok_or(CinematicFinalizationPlanError::ArithmeticOverflow(
            "per-segment output reservation",
        ))?;
    if denoise_policy == DenoisePolicy::SeparateBiasedDerivative {
        reserved_per_segment = reserved_per_segment
            .checked_add(artifact_ceilings.denoised_intermediate_bytes)
            .ok_or(CinematicFinalizationPlanError::ArithmeticOverflow(
                "per-segment output reservation",
            ))?;
    }
    let reserved_output_bytes = segment_count.checked_mul(reserved_per_segment).ok_or(
        CinematicFinalizationPlanError::ArithmeticOverflow("sequence output reservation"),
    )?;
    let output_limit = sequence_limits
        .max_output_bytes()
        .min(profile_output_ceiling);
    if reserved_output_bytes > output_limit {
        return Err(CinematicFinalizationPlanError::ResourceLimit {
            resource: "sequence output reservation",
            requested: reserved_output_bytes,
            limit: output_limit,
        });
    }
    usize::try_from(artifact_count).map_err(|_| {
        CinematicFinalizationPlanError::ArithmeticOverflow("expected frame artifact capacity")
    })
}

fn map_aov_plan_error(error: CinematicAovError) -> CinematicFinalizationPlanError {
    match error {
        CinematicAovError::Tracer(TracerError::Cancelled) => {
            CinematicFinalizationPlanError::Cancelled
        }
        CinematicAovError::PaletteLimit {
            requested, limit, ..
        } => CinematicFinalizationPlanError::ResourceLimit {
            resource: "AOV palette entries",
            requested,
            limit,
        },
        CinematicAovError::PixelLimit { requested, limit } => {
            CinematicFinalizationPlanError::ResourceLimit {
                resource: "AOV pixels",
                requested,
                limit,
            }
        }
        CinematicAovError::AllocationRefused => {
            CinematicFinalizationPlanError::Capacity("scene AOV palette")
        }
        CinematicAovError::SizeOverflow => {
            CinematicFinalizationPlanError::ArithmeticOverflow("scene AOV palette")
        }
        _ => CinematicFinalizationPlanError::Incompatible("scene AOV palette"),
    }
}

fn validate_event_sequence(
    events: &[ResampledAudioEvent],
    total_audio_sample_frames: u64,
    cx: &Cx<'_>,
) -> Result<(), CinematicFinalizationPlanError> {
    let mut previous_time_s = None;
    for (index, event) in events.iter().enumerate() {
        if index.is_multiple_of(1_024) {
            checkpoint_plan(cx)?;
        }
        u64::try_from(event.source.source_sample_index).map_err(|_| {
            CinematicFinalizationPlanError::ArithmeticOverflow("audio source sample index")
        })?;
        validate_resampled_audio_event(event, total_audio_sample_frames).map_err(|_| {
            CinematicFinalizationPlanError::Incompatible("audio event placement receipt")
        })?;
        if previous_time_s.is_some_and(|previous| event.source.time_s <= previous) {
            return Err(CinematicFinalizationPlanError::Incompatible(
                "strictly increasing audio event time",
            ));
        }
        previous_time_s = Some(event.source.time_s);
    }
    Ok(())
}

fn aov_channels(
    profile: CinematicAovProfile,
) -> Result<Vec<FrameChannel>, CinematicFinalizationPlanError> {
    profile
        .exr_channel_layout()
        .iter()
        .map(|(name, pixel_type)| {
            FrameChannel::try_new(
                *name,
                match pixel_type {
                    PixelType::Half => FrameChannelType::Float16,
                    PixelType::Float => FrameChannelType::Float32,
                },
            )
            .map_err(CinematicFinalizationPlanError::Sequence)
        })
        .collect()
}

fn rgb_channels(
    sample_type: FrameChannelType,
) -> Result<Vec<FrameChannel>, CinematicFinalizationPlanError> {
    ["R", "G", "B"]
        .into_iter()
        .map(|name| {
            FrameChannel::try_new(name, sample_type)
                .map_err(CinematicFinalizationPlanError::Sequence)
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn frame_provenance(
    frame: u64,
    total_frames: u32,
    frames_per_second: u32,
    trajectory: ContentHash,
    scene: ContentHash,
    composition: ContentHash,
) -> Result<CinematicAovProvenance, CinematicFinalizationPlanError> {
    let fps = f64::from(frames_per_second);
    let current = frame as f64 / fps;
    let previous = frame.saturating_sub(1) as f64 / fps;
    let next = frame.saturating_add(1).min(u64::from(total_frames)) as f64 / fps;
    CinematicAovProvenance::try_new(
        frame,
        current,
        previous,
        next,
        trajectory,
        scene,
        composition,
    )
    .map_err(|_| CinematicFinalizationPlanError::Incompatible("AOV provenance"))
}

fn raw_attributes(
    config: CinematicAovConfig,
    settings: fs_render::tracer::Settings,
    shot_id: u64,
    cut_side: CutSide,
    shutter: fs_render::motion::ShutterInterval,
    object_palette: &Arc<[u8]>,
    material_palette: &Arc<[u8]>,
    render_semantics_versions: &Arc<[u8]>,
) -> BTreeMap<String, Arc<[u8]>> {
    let provenance = config.provenance();
    let mut attributes = BTreeMap::new();
    let mut insert = |name: &str, value: String| {
        attributes.insert(name.to_owned(), Arc::from(value.into_bytes()));
    };
    insert("frankensim.aov.authority", "raw-estimate".to_owned());
    insert(
        "frankensim.aov.schemaVersion",
        CINEMATIC_AOV_SEMANTICS_VERSION.to_string(),
    );
    insert("frankensim.aov.profile", config.profile().code().to_owned());
    insert("frankensim.aov.configHash", config.identity().to_hex());
    insert(
        "frankensim.aov.channelSemantics",
        CINEMATIC_AOV_CHANNEL_SEMANTICS.to_owned(),
    );
    insert(
        "frankensim.aov.invalidSemantics",
        CINEMATIC_AOV_INVALID_SEMANTICS.to_owned(),
    );
    insert(
        "frankensim.aov.materialDomain",
        fs_render::tracer::MATERIAL_CONTENT_IDENTITY_DOMAIN.to_owned(),
    );
    insert(
        "frankensim.frame.index",
        provenance.frame_index().to_string(),
    );
    insert(
        "frankensim.frame.timeSeconds",
        f64_bits_string(provenance.frame_time_s()),
    );
    insert(
        "frankensim.frame.previousTimeS",
        f64_bits_string(provenance.previous_frame_time_s()),
    );
    insert(
        "frankensim.frame.nextTimeS",
        f64_bits_string(provenance.next_frame_time_s()),
    );
    insert(
        "frankensim.source.trajectory",
        provenance.source_trajectory_identity().to_hex(),
    );
    insert(
        "frankensim.source.sceneHash",
        provenance.scene_identity().to_hex(),
    );
    insert(
        "frankensim.source.composition",
        provenance.composition_identity().to_hex(),
    );
    insert("frankensim.render.seed", settings.seed.to_string());
    insert(
        "frankensim.render.sampler",
        match settings.sampler {
            Sampler::Iid => "iid-philox",
            Sampler::OwenSobol => "owen-sobol",
        }
        .to_owned(),
    );
    insert(
        "frankensim.render.strategy",
        match settings.strategy {
            DirectStrategy::NeeOnly => "nee-only",
            DirectStrategy::BsdfOnly => "bsdf-only",
            DirectStrategy::Mis => "mis",
        }
        .to_owned(),
    );
    insert("frankensim.render.maxDepth", settings.max_depth.to_string());
    insert("frankensim.render.sampleMode", "uniform".to_owned());
    insert("frankensim.render.spp", settings.spp.to_string());
    insert("frankensim.render.sppCeiling", settings.spp.to_string());
    insert("frankensim.render.shotId", shot_id.to_string());
    insert(
        "frankensim.render.cutSide",
        match cut_side {
            CutSide::Before => "before",
            CutSide::After => "after",
        }
        .to_owned(),
    );
    insert(
        "frankensim.render.shutterOpenS",
        f64_bits_string(shutter.open_s()),
    );
    insert(
        "frankensim.render.shutterCloseS",
        f64_bits_string(shutter.close_s()),
    );
    let convention = match shutter.convention() {
        ShutterConvention::Centered => "centered",
        ShutterConvention::FrontLoaded => "front-loaded",
        ShutterConvention::BackLoaded => "back-loaded",
    };
    let (distribution, strata) = match shutter.distribution() {
        ShutterDistribution::UniformCounterV1 => ("uniform-counter-v1", 0),
        ShutterDistribution::StratifiedCounterV1 { strata } => ("stratified-counter-v1", strata),
    };
    insert(
        "frankensim.render.shutter",
        format!("convention={convention};distribution={distribution};strata={strata}"),
    );
    drop(insert);
    attributes.insert(
        "frankensim.render.versions".to_owned(),
        Arc::clone(render_semantics_versions),
    );
    attributes.insert(
        "frankensim.aov.objectPalette".to_owned(),
        Arc::clone(object_palette),
    );
    attributes.insert(
        "frankensim.aov.materialPalette".to_owned(),
        Arc::clone(material_palette),
    );
    let mut insert = |name: &str, value: String| {
        attributes.insert(name.to_owned(), Arc::from(value.into_bytes()));
    };
    insert(
        "frankensim.aov.paletteZero",
        CINEMATIC_AOV_PALETTE_ZERO_SEMANTICS.to_owned(),
    );
    attributes
}

fn f64_bits_string(value: f64) -> String {
    let value = if value == 0.0 { 0.0 } else { value };
    format!("{}@0x{:016x}", value, value.to_bits())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BundleSizeError {
    Cancelled,
    Overflow,
}

fn bundle_size(bundle: &CinematicBundle<'_>, cx: &Cx<'_>) -> Result<u64, BundleSizeError> {
    let mut total =
        u64::try_from(bundle.sequence_bytes.len()).map_err(|_| BundleSizeError::Overflow)?;
    for (index, frame) in bundle.frames.iter().enumerate() {
        if index % 256 == 0 {
            cx.checkpoint().map_err(|_| BundleSizeError::Cancelled)?;
        }
        total = total
            .checked_add(u64::try_from(frame.bytes.len()).map_err(|_| BundleSizeError::Overflow)?)
            .ok_or(BundleSizeError::Overflow)?;
        total = total
            .checked_add(
                u64::try_from(frame.authority_bytes.len())
                    .map_err(|_| BundleSizeError::Overflow)?,
            )
            .ok_or(BundleSizeError::Overflow)?;
    }
    cx.checkpoint().map_err(|_| BundleSizeError::Cancelled)?;
    for bytes in [
        bundle.audio.wav_bytes,
        bundle.audio.manifest_bytes,
        bundle.audio.authority_bytes,
        bundle.audio.alignment_bytes,
        bundle.audio.event_bytes,
    ] {
        total = total
            .checked_add(u64::try_from(bytes.len()).map_err(|_| BundleSizeError::Overflow)?)
            .ok_or(BundleSizeError::Overflow)?;
    }
    Ok(total)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ReceiptDecodeError {
    Missing,
    Cancelled,
    Budget,
    Corrupt,
    Incompatible,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct AudioManifestWavSnapshot {
    wav_identity: ContentHash,
    metadata_identity: ContentHash,
    byte_len: u64,
    sample_frame_count: u64,
    sample_rate_hz: u32,
    encoding: WavSampleEncoding,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct DecodedAudioManifest {
    manifest_identity: ContentHash,
    synthesis: SoundSynthesisReceipt,
    authority: SoundAuthority,
    channel_layout_identity: ContentHash,
    signal_path: AudioSignalPath,
    source_signal_identity: ContentHash,
    mix_identity: Option<ContentHash>,
    wav: AudioManifestWavSnapshot,
    role: AudioArtifactRole,
    meters: AudioMeters,
    video_start_tick: i64,
    video_end_tick_exclusive: i64,
    audio_start_tick: i64,
    audio_end_tick_exclusive: i64,
    audio_frames_per_video_frame: u32,
    admitted_headroom_db: f64,
}

fn decode_audio_manifest_receipt(
    bytes: &[u8],
    expected_identity: ContentHash,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<DecodedAudioManifest, ReceiptDecodeError> {
    if bytes.is_empty() {
        return Err(ReceiptDecodeError::Missing);
    }
    if u64::try_from(bytes.len()).map_err(|_| ReceiptDecodeError::Budget)? > max_bytes {
        return Err(ReceiptDecodeError::Budget);
    }
    let identity = hash_domain_with_cancellation(AUDIO_MANIFEST_RECEIPT_DOMAIN, bytes, cx)
        .ok_or(ReceiptDecodeError::Cancelled)?;
    if identity != expected_identity {
        return Err(ReceiptDecodeError::Corrupt);
    }
    let mut reader = ReceiptReader::new(bytes);
    if &reader.array::<8>()? != AUDIO_MANIFEST_RECEIPT_MAGIC {
        return Err(ReceiptDecodeError::Corrupt);
    }
    if reader.u16()? != RECEIPT_CODEC_VERSION {
        return Err(ReceiptDecodeError::Incompatible);
    }
    if reader.u16()? != AUDIO_ARTIFACT_SCHEMA_VERSION {
        return Err(ReceiptDecodeError::Incompatible);
    }
    let manifest_identity = reader.hash()?;
    let synthesis = SoundSynthesisReceipt {
        schema_version: reader.u16()?,
        configuration_identity: reader.hash()?,
        authority: decode_sound_authority(reader.u8()?)?,
        trajectory_identity: reader.hash()?,
        excitation_identity: reader.hash()?,
        sound_model_identity: reader.hash()?,
        timeline_identity: reader.hash()?,
    };
    if synthesis.schema_version != SOUND_SYNTHESIS_SCHEMA_VERSION {
        return Err(ReceiptDecodeError::Incompatible);
    }
    let authority = decode_sound_authority(reader.u8()?)?;
    let channel_layout_identity = reader.hash()?;
    let signal_path = match reader.u8()? {
        1 => AudioSignalPath::CanonicalDryStereo,
        2 => AudioSignalPath::SpatializedStereo {
            spatialization_identity: reader.hash()?,
        },
        _ => return Err(ReceiptDecodeError::Corrupt),
    };
    let source_signal_identity = reader.hash()?;
    let mix_identity = match reader.u8()? {
        0 => None,
        1 => Some(reader.hash()?),
        _ => return Err(ReceiptDecodeError::Corrupt),
    };
    let wav = AudioManifestWavSnapshot {
        wav_identity: reader.hash()?,
        metadata_identity: reader.hash()?,
        byte_len: reader.u64()?,
        sample_frame_count: reader.u64()?,
        sample_rate_hz: reader.u32()?,
        encoding: match reader.u16()? {
            1 => WavSampleEncoding::Pcm24,
            3 => WavSampleEncoding::Float32,
            _ => return Err(ReceiptDecodeError::Corrupt),
        },
    };
    let role = match reader.u8()? {
        1 => AudioArtifactRole::AuthoritativeFloat32Master,
        2 => AudioArtifactRole::QuantizedPcm24Derivative,
        _ => return Err(ReceiptDecodeError::Corrupt),
    };
    let meters = decode_audio_meters(&mut reader)?;
    let manifest = DecodedAudioManifest {
        manifest_identity,
        synthesis,
        authority,
        channel_layout_identity,
        signal_path,
        source_signal_identity,
        mix_identity,
        wav,
        role,
        meters,
        video_start_tick: reader.i64()?,
        video_end_tick_exclusive: reader.i64()?,
        audio_start_tick: reader.i64()?,
        audio_end_tick_exclusive: reader.i64()?,
        audio_frames_per_video_frame: reader.u32()?,
        admitted_headroom_db: reader.f64()?,
    };
    if !reader.finished() {
        return Err(ReceiptDecodeError::Corrupt);
    }
    validate_decoded_audio_manifest(&manifest)?;
    Ok(manifest)
}

fn decode_sound_authority(tag: u8) -> Result<SoundAuthority, ReceiptDecodeError> {
    match tag {
        1 => Ok(SoundAuthority::Artistic),
        2 => Ok(SoundAuthority::PhysicallyInformed),
        3 => Ok(SoundAuthority::Calibrated),
        _ => Err(ReceiptDecodeError::Corrupt),
    }
}

fn decode_audio_meters(reader: &mut ReceiptReader<'_>) -> Result<AudioMeters, ReceiptDecodeError> {
    let sample_peak_fs = reader.f64()?;
    let true_peak_estimate_fs = reader.f64()?;
    let stereo_rms_fs = reader.f64()?;
    let dc_left_fs = reader.f64()?;
    let dc_right_fs = reader.f64()?;
    let integrated_loudness_lufs = match reader.u8()? {
        0 => None,
        1 => Some(reader.f64()?),
        _ => return Err(ReceiptDecodeError::Corrupt),
    };
    Ok(AudioMeters {
        sample_peak_fs,
        true_peak_estimate_fs,
        stereo_rms_fs,
        dc_left_fs,
        dc_right_fs,
        integrated_loudness_lufs,
        loudness_block_count: reader.u64()?,
        absolute_gated_block_count: reader.u64()?,
        relative_gated_block_count: reader.u64()?,
    })
}

fn validate_decoded_audio_manifest(
    manifest: &DecodedAudioManifest,
) -> Result<(), ReceiptDecodeError> {
    let identities = [
        manifest.manifest_identity,
        manifest.synthesis.configuration_identity,
        manifest.synthesis.trajectory_identity,
        manifest.synthesis.excitation_identity,
        manifest.synthesis.sound_model_identity,
        manifest.synthesis.timeline_identity,
        manifest.channel_layout_identity,
        manifest.source_signal_identity,
        manifest.wav.wav_identity,
        manifest.wav.metadata_identity,
    ];
    if identities.into_iter().any(is_zero)
        || manifest.mix_identity.is_some_and(is_zero)
        || matches!(
            manifest.signal_path,
            AudioSignalPath::SpatializedStereo {
                spatialization_identity
            } if is_zero(spatialization_identity)
        )
        || manifest.authority != manifest.synthesis.authority
        || manifest.channel_layout_identity != channel_layout_identity(manifest.signal_path)
        || manifest.manifest_identity != decoded_audio_manifest_identity(manifest)
        || manifest.wav.byte_len == 0
        || manifest.wav.sample_frame_count == 0
        || manifest.wav.sample_rate_hz == 0
        || manifest.video_start_tick >= manifest.video_end_tick_exclusive
        || manifest.audio_start_tick >= manifest.audio_end_tick_exclusive
        || manifest.audio_frames_per_video_frame == 0
        || !manifest.admitted_headroom_db.is_finite()
        || manifest.admitted_headroom_db < 0.0
        || !valid_audio_meters(manifest.meters)
        || !matches!(
            (manifest.wav.encoding, manifest.role),
            (
                WavSampleEncoding::Float32,
                AudioArtifactRole::AuthoritativeFloat32Master
            ) | (
                WavSampleEncoding::Pcm24,
                AudioArtifactRole::QuantizedPcm24Derivative
            )
        )
        || !matches!(
            (manifest.signal_path, manifest.mix_identity),
            (AudioSignalPath::CanonicalDryStereo, Some(_))
                | (AudioSignalPath::SpatializedStereo { .. }, None)
        )
    {
        return Err(ReceiptDecodeError::Corrupt);
    }
    Ok(())
}

fn valid_audio_meters(meters: AudioMeters) -> bool {
    [
        meters.sample_peak_fs,
        meters.true_peak_estimate_fs,
        meters.stereo_rms_fs,
        meters.dc_left_fs,
        meters.dc_right_fs,
    ]
    .into_iter()
    .all(f64::is_finite)
        && meters.sample_peak_fs >= 0.0
        && meters.true_peak_estimate_fs >= 0.0
        && meters.stereo_rms_fs >= 0.0
        && meters.integrated_loudness_lufs.is_none_or(f64::is_finite)
        && meters.relative_gated_block_count <= meters.absolute_gated_block_count
        && meters.absolute_gated_block_count <= meters.loudness_block_count
}

const fn sound_authority_tag(authority: SoundAuthority) -> u8 {
    match authority {
        SoundAuthority::Artistic => 1,
        SoundAuthority::PhysicallyInformed => 2,
        SoundAuthority::Calibrated => 3,
    }
}

const fn audio_artifact_role_tag(role: AudioArtifactRole) -> u8 {
    match role {
        AudioArtifactRole::AuthoritativeFloat32Master => 1,
        AudioArtifactRole::QuantizedPcm24Derivative => 2,
    }
}

fn push_audio_meters(bytes: &mut Vec<u8>, meters: AudioMeters) {
    for value in [
        meters.sample_peak_fs,
        meters.true_peak_estimate_fs,
        meters.stereo_rms_fs,
        meters.dc_left_fs,
        meters.dc_right_fs,
    ] {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    match meters.integrated_loudness_lufs {
        None => bytes.push(0),
        Some(value) => {
            bytes.push(1);
            bytes.extend_from_slice(&value.to_bits().to_le_bytes());
        }
    }
    bytes.extend_from_slice(&meters.loudness_block_count.to_le_bytes());
    bytes.extend_from_slice(&meters.absolute_gated_block_count.to_le_bytes());
    bytes.extend_from_slice(&meters.relative_gated_block_count.to_le_bytes());
}

fn hash_audio_meters(hasher: &mut DomainHasher, meters: AudioMeters) {
    for value in [
        meters.sample_peak_fs,
        meters.true_peak_estimate_fs,
        meters.stereo_rms_fs,
        meters.dc_left_fs,
        meters.dc_right_fs,
    ] {
        hasher.update(&value.to_bits().to_le_bytes());
    }
    match meters.integrated_loudness_lufs {
        None => hasher.update(&[0]),
        Some(value) => {
            hasher.update(&[1]);
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    hasher.update(&meters.loudness_block_count.to_le_bytes());
    hasher.update(&meters.absolute_gated_block_count.to_le_bytes());
    hasher.update(&meters.relative_gated_block_count.to_le_bytes());
}

fn channel_layout_identity(path: AudioSignalPath) -> ContentHash {
    let mut hasher = DomainHasher::new(CHANNEL_RECEIPT_IDENTITY_DOMAIN);
    hasher.update(&[2]);
    match path {
        AudioSignalPath::CanonicalDryStereo => hasher.update(&[1]),
        AudioSignalPath::SpatializedStereo {
            spatialization_identity,
        } => {
            hasher.update(&[2]);
            hasher.update(spatialization_identity.as_bytes());
        }
    }
    hasher.finalize()
}

fn decoded_audio_manifest_identity(manifest: &DecodedAudioManifest) -> ContentHash {
    let mut hasher = DomainHasher::new(AUDIO_MANIFEST_IDENTITY_DOMAIN);
    hasher.update(&AUDIO_ARTIFACT_SCHEMA_VERSION.to_le_bytes());
    let synthesis = manifest.synthesis;
    hasher.update(&synthesis.schema_version.to_le_bytes());
    hasher.update(synthesis.configuration_identity.as_bytes());
    hasher.update(synthesis.authority.code().as_bytes());
    hasher.update(synthesis.trajectory_identity.as_bytes());
    hasher.update(synthesis.excitation_identity.as_bytes());
    hasher.update(synthesis.sound_model_identity.as_bytes());
    hasher.update(synthesis.timeline_identity.as_bytes());
    hasher.update(manifest.authority.code().as_bytes());
    hasher.update(manifest.channel_layout_identity.as_bytes());
    hasher.update(manifest.source_signal_identity.as_bytes());
    match manifest.mix_identity {
        None => hasher.update(&[0]),
        Some(identity) => {
            hasher.update(&[1]);
            hasher.update(identity.as_bytes());
        }
    }
    hasher.update(manifest.wav.wav_identity.as_bytes());
    hasher.update(manifest.wav.metadata_identity.as_bytes());
    hasher.update(&manifest.wav.byte_len.to_le_bytes());
    hasher.update(&manifest.wav.sample_frame_count.to_le_bytes());
    hasher.update(&manifest.wav.sample_rate_hz.to_le_bytes());
    hasher.update(&(manifest.wav.encoding as u16).to_le_bytes());
    hasher.update(manifest.role.code().as_bytes());
    hash_audio_meters(&mut hasher, manifest.meters);
    hasher.update(&manifest.video_start_tick.to_le_bytes());
    hasher.update(&manifest.video_end_tick_exclusive.to_le_bytes());
    hasher.update(&manifest.audio_start_tick.to_le_bytes());
    hasher.update(&manifest.audio_end_tick_exclusive.to_le_bytes());
    hasher.update(&manifest.audio_frames_per_video_frame.to_le_bytes());
    hasher.update(&manifest.admitted_headroom_db.to_bits().to_le_bytes());
    hasher.finalize()
}

fn map_receipt_decode_error(
    error: ReceiptDecodeError,
    semantic_code: CinematicFinalizationDivergenceCode,
) -> (
    CinematicFinalizationDisposition,
    CinematicFinalizationDivergenceCode,
) {
    match error {
        ReceiptDecodeError::Missing => {
            (CinematicFinalizationDisposition::Incomplete, semantic_code)
        }
        ReceiptDecodeError::Cancelled => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ),
        ReceiptDecodeError::Budget => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
        ),
        ReceiptDecodeError::Corrupt => (CinematicFinalizationDisposition::Corrupt, semantic_code),
        ReceiptDecodeError::Incompatible => (
            CinematicFinalizationDisposition::Incompatible,
            semantic_code,
        ),
    }
}

fn decode_alignment_receipt(
    bytes: &[u8],
    expected_identity: ContentHash,
    max_markers: u32,
    cx: &Cx<'_>,
) -> Result<AudioVideoAlignment, ReceiptDecodeError> {
    if bytes.is_empty() {
        return Err(ReceiptDecodeError::Missing);
    }
    let identity = hash_domain_with_cancellation(ALIGNMENT_RECEIPT_DOMAIN, bytes, cx)
        .ok_or(ReceiptDecodeError::Cancelled)?;
    if identity != expected_identity {
        return Err(ReceiptDecodeError::Corrupt);
    }
    let mut reader = ReceiptReader::new(bytes);
    if &reader.array::<8>()? != ALIGNMENT_RECEIPT_MAGIC {
        return Err(ReceiptDecodeError::Corrupt);
    }
    if reader.u16()? != RECEIPT_CODEC_VERSION {
        return Err(ReceiptDecodeError::Incompatible);
    }
    let ratio = reader.u32()?;
    let drift = reader.i64()?;
    let count = reader.u32()?;
    let exact_marker_bytes = u64::from(count)
        .checked_mul(24)
        .ok_or(ReceiptDecodeError::Corrupt)?;
    if u64::try_from(reader.remaining()).map_err(|_| ReceiptDecodeError::Corrupt)?
        != exact_marker_bytes
    {
        return Err(ReceiptDecodeError::Corrupt);
    }
    if count > max_markers {
        return Err(ReceiptDecodeError::Budget);
    }
    let count_usize = usize::try_from(count).map_err(|_| ReceiptDecodeError::Budget)?;
    let mut markers = Vec::new();
    markers
        .try_reserve_exact(count_usize)
        .map_err(|_| ReceiptDecodeError::Budget)?;
    for index in 0..count_usize {
        if index % 1_024 == 0 && cx.checkpoint().is_err() {
            return Err(ReceiptDecodeError::Cancelled);
        }
        markers.push(AudioVideoSyncMarker {
            video_tick: reader.i64()?,
            audio_tick: reader.i64()?,
            audio_frame_offset: reader.u64()?,
        });
    }
    if !reader.finished() {
        return Err(ReceiptDecodeError::Corrupt);
    }
    let alignment = AudioVideoAlignment {
        audio_frames_per_video_frame: ratio,
        markers,
        endpoint_drift_audio_frames: drift,
    };
    Ok(alignment)
}

fn decode_event_receipt(
    bytes: &[u8],
    expected_identity: ContentHash,
    max_events: u32,
    cx: &Cx<'_>,
) -> Result<Vec<ResampledAudioEvent>, ReceiptDecodeError> {
    const MIN_EVENT_BYTES: u64 = 109;
    if bytes.is_empty() {
        return Err(ReceiptDecodeError::Missing);
    }
    let identity = hash_domain_with_cancellation(EVENT_RECEIPT_DOMAIN, bytes, cx)
        .ok_or(ReceiptDecodeError::Cancelled)?;
    if identity != expected_identity {
        return Err(ReceiptDecodeError::Corrupt);
    }
    let mut reader = ReceiptReader::new(bytes);
    if &reader.array::<8>()? != EVENT_RECEIPT_MAGIC {
        return Err(ReceiptDecodeError::Corrupt);
    }
    if reader.u16()? != RECEIPT_CODEC_VERSION {
        return Err(ReceiptDecodeError::Incompatible);
    }
    let count = reader.u32()?;
    let minimum_bytes = u64::from(count)
        .checked_mul(MIN_EVENT_BYTES)
        .ok_or(ReceiptDecodeError::Corrupt)?;
    if u64::try_from(reader.remaining()).map_err(|_| ReceiptDecodeError::Corrupt)? < minimum_bytes {
        return Err(ReceiptDecodeError::Corrupt);
    }
    if count > max_events {
        return Err(ReceiptDecodeError::Budget);
    }
    let count_usize = usize::try_from(count).map_err(|_| ReceiptDecodeError::Budget)?;
    let mut events = Vec::new();
    events
        .try_reserve_exact(count_usize)
        .map_err(|_| ReceiptDecodeError::Budget)?;
    for index in 0..count_usize {
        if index % 1_024 == 0 && cx.checkpoint().is_err() {
            return Err(ReceiptDecodeError::Cancelled);
        }
        events.push(decode_resampled_event(&mut reader)?);
    }
    if !reader.finished() {
        return Err(ReceiptDecodeError::Corrupt);
    }
    Ok(events)
}

fn decode_resampled_event(
    reader: &mut ReceiptReader<'_>,
) -> Result<ResampledAudioEvent, ReceiptDecodeError> {
    let source_sample_index =
        usize::try_from(reader.u64()?).map_err(|_| ReceiptDecodeError::Incompatible)?;
    let kind = match reader.u8()? {
        1 => crate::coupled_runner::ContactTransitionKind::Opening,
        2 => crate::coupled_runner::ContactTransitionKind::Reimpact,
        _ => return Err(ReceiptDecodeError::Corrupt),
    };
    let measure = match reader.u8()? {
        1 => crate::control_stream::ContactEventMeasure::TimingOnly,
        _ => return Err(ReceiptDecodeError::Corrupt),
    };
    let time_s = reader.f64()?;
    let bracket_start_s = reader.f64()?;
    let bracket_end_s = reader.f64()?;
    let requested_sample_position = reader.f64()?;
    let left_weight = reader.f64()?;
    let right_weight = reader.f64()?;
    let centroid_error_frames = reader.f64()?;
    let bracket_start_sample_position = reader.f64()?;
    let bracket_end_sample_position = reader.f64()?;
    let physical_impulse_n_s = reader.modal_values()?;
    let artistic = match reader.u8()? {
        0 => None,
        1 => Some(crate::audio_excitation::ArtisticEventExcitation {
            stream_identity: ContentHash(reader.array::<32>()?),
            impulse_n_s: reader.modal_values()?,
        }),
        _ => return Err(ReceiptDecodeError::Corrupt),
    };
    let left_frame_offset = reader.optional_u64()?;
    let right_frame_offset = reader.optional_u64()?;
    Ok(ResampledAudioEvent {
        source: crate::audio_excitation::AudioExcitationEvent {
            source_sample_index,
            kind,
            time_s,
            bracket_start_s,
            bracket_end_s,
            measure,
            physical_impulse_n_s,
            artistic,
        },
        requested_sample_position,
        left_frame_offset,
        right_frame_offset,
        left_weight,
        right_weight,
        centroid_error_frames,
        bracket_start_sample_position,
        bracket_end_sample_position,
    })
}

struct ReceiptReader<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> ReceiptReader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, position: 0 }
    }

    fn array<const N: usize>(&mut self) -> Result<[u8; N], ReceiptDecodeError> {
        let end = self
            .position
            .checked_add(N)
            .ok_or(ReceiptDecodeError::Corrupt)?;
        let value = self
            .bytes
            .get(self.position..end)
            .ok_or(ReceiptDecodeError::Corrupt)?;
        self.position = end;
        value.try_into().map_err(|_| ReceiptDecodeError::Corrupt)
    }

    fn u8(&mut self) -> Result<u8, ReceiptDecodeError> {
        Ok(self.array::<1>()?[0])
    }

    fn u16(&mut self) -> Result<u16, ReceiptDecodeError> {
        Ok(u16::from_le_bytes(self.array()?))
    }

    fn u32(&mut self) -> Result<u32, ReceiptDecodeError> {
        Ok(u32::from_le_bytes(self.array()?))
    }

    fn u64(&mut self) -> Result<u64, ReceiptDecodeError> {
        Ok(u64::from_le_bytes(self.array()?))
    }

    fn i64(&mut self) -> Result<i64, ReceiptDecodeError> {
        Ok(i64::from_le_bytes(self.array()?))
    }

    fn f64(&mut self) -> Result<f64, ReceiptDecodeError> {
        Ok(f64::from_bits(self.u64()?))
    }

    fn hash(&mut self) -> Result<ContentHash, ReceiptDecodeError> {
        Ok(ContentHash(self.array::<32>()?))
    }

    fn modal_values(
        &mut self,
    ) -> Result<crate::modal_synthesis::ModalComponentValues, ReceiptDecodeError> {
        Ok(crate::modal_synthesis::ModalComponentValues {
            disc: self.f64()?,
            glass_plate: self.f64()?,
            base_assembly: self.f64()?,
        })
    }

    fn optional_u64(&mut self) -> Result<Option<u64>, ReceiptDecodeError> {
        match self.u8()? {
            0 => Ok(None),
            1 => Ok(Some(self.u64()?)),
            _ => Err(ReceiptDecodeError::Corrupt),
        }
    }

    fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.position)
    }

    fn finished(&self) -> bool {
        self.position == self.bytes.len()
    }
}

fn hash_with_cancellation(bytes: &[u8], cx: &Cx<'_>) -> Option<ContentHash> {
    let mut hasher = Blake3::new();
    for chunk in bytes.chunks(HASH_POLL_BYTES) {
        cx.checkpoint().ok()?;
        hasher.update(chunk);
    }
    cx.checkpoint().ok()?;
    Some(hasher.finalize())
}

fn hash_domain_with_cancellation(domain: &str, bytes: &[u8], cx: &Cx<'_>) -> Option<ContentHash> {
    let mut hasher = DomainHasher::new(domain);
    for chunk in bytes.chunks(HASH_POLL_BYTES) {
        cx.checkpoint().ok()?;
        hasher.update(chunk);
    }
    cx.checkpoint().ok()?;
    Some(hasher.finalize())
}

struct PlanIdentityHasher<'cx, 'scope> {
    inner: DomainHasher,
    cx: &'cx Cx<'scope>,
}

impl<'cx, 'scope> PlanIdentityHasher<'cx, 'scope> {
    fn new(cx: &'cx Cx<'scope>) -> Result<Self, CinematicFinalizationPlanError> {
        checkpoint_plan(cx)?;
        Ok(Self {
            inner: DomainHasher::new(CINEMATIC_FINALIZATION_PLAN_DOMAIN),
            cx,
        })
    }

    fn update(&mut self, bytes: &[u8]) -> Result<(), CinematicFinalizationPlanError> {
        if bytes.is_empty() {
            checkpoint_plan(self.cx)?;
            return Ok(());
        }
        for chunk in bytes.chunks(HASH_POLL_BYTES) {
            checkpoint_plan(self.cx)?;
            self.inner.update(chunk);
        }
        Ok(())
    }

    fn byte(&mut self, value: u8) -> Result<(), CinematicFinalizationPlanError> {
        self.update(&[value])
    }

    fn string(&mut self, value: &str) -> Result<(), CinematicFinalizationPlanError> {
        let length = u32::try_from(value.len()).map_err(|_| {
            CinematicFinalizationPlanError::ArithmeticOverflow("finalization plan string")
        })?;
        self.update(&length.to_le_bytes())?;
        self.update(value.as_bytes())
    }

    fn finish(self) -> Result<ContentHash, CinematicFinalizationPlanError> {
        checkpoint_plan(self.cx)?;
        Ok(self.inner.finalize())
    }
}

fn inspect_frame<'a>(
    descriptor: &FrameArtifactDescriptor,
    bytes: &'a [u8],
    limits: CinematicFinalizationLimits,
    cx: &Cx<'_>,
) -> Result<
    FrameInspection<'a>,
    (
        CinematicFinalizationDisposition,
        CinematicFinalizationDivergenceCode,
    ),
> {
    match descriptor.format() {
        FrameArtifactFormat::OpenExr => {
            let inspection = inspect_exr_with_poll(bytes, limits.exr, || cx.checkpoint().is_ok())
                .map_err(map_image_error)?;
            if inspection.width != descriptor.width() || inspection.height != descriptor.height() {
                return Err((
                    CinematicFinalizationDisposition::Incompatible,
                    CinematicFinalizationDivergenceCode::ImageDimensions,
                ));
            }
            let channels_match = inspection.channels.len() == descriptor.channels().len()
                && inspection.channels.iter().zip(descriptor.channels()).all(
                    |(actual, expected)| {
                        actual.name == expected.name()
                            && match actual.ty {
                                PixelType::Half => {
                                    expected.sample_type() == FrameChannelType::Float16
                                }
                                PixelType::Float => {
                                    expected.sample_type() == FrameChannelType::Float32
                                }
                            }
                    },
                );
            if !channels_match {
                return Err((
                    CinematicFinalizationDisposition::Incompatible,
                    CinematicFinalizationDivergenceCode::ImageChannels,
                ));
            }
            Ok(FrameInspection::Exr(inspection))
        }
        FrameArtifactFormat::Png8 | FrameArtifactFormat::Png16 => {
            let inspection = inspect_png_with_poll(bytes, limits.png, || cx.checkpoint().is_ok())
                .map_err(map_image_error)?;
            if inspection.width != descriptor.width() || inspection.height != descriptor.height() {
                return Err((
                    CinematicFinalizationDisposition::Incompatible,
                    CinematicFinalizationDivergenceCode::ImageDimensions,
                ));
            }
            let expected_depth = match descriptor.format() {
                FrameArtifactFormat::Png8 => 8,
                FrameArtifactFormat::Png16 => 16,
                FrameArtifactFormat::OpenExr => unreachable!(),
            };
            let expected_color = match descriptor.channels().len() {
                1 => PngColor::Gray,
                3 => PngColor::Rgb,
                4 => PngColor::Rgba,
                _ => {
                    return Err((
                        CinematicFinalizationDisposition::Incompatible,
                        CinematicFinalizationDivergenceCode::ImageChannels,
                    ));
                }
            };
            if inspection.depth != expected_depth || inspection.color != expected_color {
                return Err((
                    CinematicFinalizationDisposition::Incompatible,
                    CinematicFinalizationDivergenceCode::ImageChannels,
                ));
            }
            Ok(FrameInspection::Png)
        }
    }
}

enum FrameInspection<'a> {
    Exr(ExrInspection<'a>),
    Png,
}

fn map_image_error(
    error: ImgError,
) -> (
    CinematicFinalizationDisposition,
    CinematicFinalizationDivergenceCode,
) {
    match error {
        ImgError::Cancelled { .. } => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ),
        ImgError::ResourceLimit { .. } | ImgError::AllocationRefused { .. } => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
        ),
        ImgError::Unsupported { .. } => (
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::ImageStructure,
        ),
        ImgError::SizeOverflow { .. } | ImgError::Shape { .. } | ImgError::Malformed { .. } => (
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::ImageStructure,
        ),
    }
}

fn map_sequence_decode_error(
    error: &FrameSequenceError,
) -> (
    CinematicFinalizationDisposition,
    CinematicFinalizationDivergenceCode,
) {
    match error {
        FrameSequenceError::ResourceLimit { .. } | FrameSequenceError::AllocationRefused { .. } => {
            (
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
            )
        }
        FrameSequenceError::Cancelled => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ),
        FrameSequenceError::UnsupportedVersion { .. } => (
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::SequenceDecode,
        ),
        _ => (
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::SequenceDecode,
        ),
    }
}

fn map_audio_artifact_error(
    error: &AudioArtifactError,
) -> (
    CinematicFinalizationDisposition,
    CinematicFinalizationDivergenceCode,
) {
    match error {
        AudioArtifactError::Cancelled => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ),
        AudioArtifactError::InvalidBudget(_)
        | AudioArtifactError::BudgetExceeded { .. }
        | AudioArtifactError::Capacity { .. } => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
        ),
        AudioArtifactError::ManifestIdentityMismatch => (
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::AudioManifestIdentity,
        ),
        AudioArtifactError::UnsupportedWav(reason) if is_noncanonical_wav_structure(reason) => (
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::WavStructure,
        ),
        AudioArtifactError::UnsupportedWav(_) | AudioArtifactError::InvalidSampleRate(_) => (
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::WavStructure,
        ),
        _ => (
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::WavStructure,
        ),
    }
}

fn is_noncanonical_wav_structure(reason: &str) -> bool {
    matches!(
        reason,
        "noncanonical chunk before fmt"
            | "noncanonical chunk before float fact"
            | "noncanonical or unknown chunk"
            | "chunk after data"
            | "non-INFO LIST chunk"
            | "noncanonical INFO metadata"
    )
}

fn verify_wav_against_snapshot(
    manifest: &DecodedAudioManifest,
    wav_bytes: &[u8],
    budget: AudioArtifactBudget,
    cx: &Cx<'_>,
) -> Result<
    (),
    (
        CinematicFinalizationDisposition,
        CinematicFinalizationDivergenceCode,
    ),
> {
    if u64::try_from(wav_bytes.len()).ok() != Some(manifest.wav.byte_len) {
        return Err((
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::WavStructure,
        ));
    }
    let decoded = decode_stereo_wav(wav_bytes, budget, cx)
        .map_err(|error| map_audio_artifact_error(&error))?;
    let receipt = decoded.receipt;
    if receipt.wav_identity() != manifest.wav.wav_identity
        || receipt.metadata_identity() != manifest.wav.metadata_identity
        || receipt.byte_len() != manifest.wav.byte_len
        || receipt.sample_frame_count() != manifest.wav.sample_frame_count
        || receipt.sample_rate_hz() != manifest.wav.sample_rate_hz
        || receipt.encoding() != manifest.wav.encoding
    {
        return Err((
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::WavStructure,
        ));
    }
    let meters = measure_audio(&decoded.samples, budget, cx)
        .map_err(|error| map_audio_artifact_error(&error))?;
    if meters != manifest.meters {
        return Err((
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::WavStructure,
        ));
    }
    let allowed_peak_fs =
        det::exp(-manifest.admitted_headroom_db * core::f64::consts::LN_10 / 20.0);
    if meters.sample_peak_fs.max(meters.true_peak_estimate_fs) > allowed_peak_fs {
        return Err((
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::WavStructure,
        ));
    }
    Ok(())
}

fn verify_raw_metadata(
    inspection: &FrameInspection<'_>,
    expected: &RawFrameExpectation,
) -> Result<(), CinematicFinalizationDivergenceCode> {
    let FrameInspection::Exr(inspection) = inspection else {
        return Err(CinematicFinalizationDivergenceCode::ImageMetadata);
    };
    if inspection.attributes.len() != expected.required_attributes.len() {
        return Err(CinematicFinalizationDivergenceCode::ImageMetadata);
    }
    for (actual, (expected_name, expected_value)) in inspection
        .attributes
        .iter()
        .zip(&expected.required_attributes)
    {
        if actual.name != expected_name
            || actual.ty != "string"
            || actual.value != expected_value.as_ref()
        {
            return Err(CinematicFinalizationDivergenceCode::ImageMetadata);
        }
    }
    Ok(())
}

fn verify_raw_payload(
    bytes: &[u8],
    expected: &RawFrameExpectation,
    limits: ExrInspectLimits,
    cx: &Cx<'_>,
) -> Result<
    (),
    (
        CinematicFinalizationDisposition,
        CinematicFinalizationDivergenceCode,
    ),
> {
    let semantic_limits = ExrRawFrameSemanticLimits::try_new(
        ALLOWED_AOV_VALIDITY_BITS,
        expected.object_palette_entries,
        expected.material_palette_entries,
    )
    .map_err(map_raw_payload_error)?;
    validate_exr_raw_frame_payload_against_with_poll(bytes, semantic_limits, limits, || {
        cx.checkpoint().is_ok()
    })
    .map_err(map_raw_payload_error)
}

fn map_raw_payload_error(
    error: ImgError,
) -> (
    CinematicFinalizationDisposition,
    CinematicFinalizationDivergenceCode,
) {
    match error {
        ImgError::Cancelled { .. } => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ),
        ImgError::ResourceLimit { .. } | ImgError::AllocationRefused { .. } => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
        ),
        ImgError::Unsupported { .. } => (
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::ImagePayload,
        ),
        ImgError::SizeOverflow { .. } | ImgError::Shape { .. } | ImgError::Malformed { .. } => (
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::ImagePayload,
        ),
    }
}

fn verify_raw_sample_count(
    bytes: &[u8],
    expected_spp: u32,
    limits: ExrInspectLimits,
    cx: &Cx<'_>,
) -> Result<
    (),
    (
        CinematicFinalizationDisposition,
        CinematicFinalizationDivergenceCode,
    ),
> {
    verify_exr_float_channel_constant_with_poll(
        bytes,
        "samples",
        expected_spp as f32,
        limits,
        || cx.checkpoint().is_ok(),
    )
    .map_err(|error| match error {
        ImgError::Cancelled { .. } => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ),
        ImgError::ResourceLimit { .. } | ImgError::AllocationRefused { .. } => (
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
        ),
        ImgError::Unsupported { .. } | ImgError::Shape { .. } | ImgError::Malformed { .. } => (
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::ImageSampleCount,
        ),
        ImgError::SizeOverflow { .. } => (
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::ImageSampleCount,
        ),
    })
}

fn verify_derived_source(
    inspection: &FrameInspection<'_>,
    expected_source: ContentHash,
) -> Result<(), CinematicFinalizationDivergenceCode> {
    let FrameInspection::Exr(inspection) = inspection else {
        return Err(CinematicFinalizationDivergenceCode::ImageMetadata);
    };
    let expected_source_hex = expected_source.to_hex();
    if inspection
        .attributes
        .binary_search_by(|attribute| attribute.name.cmp(SOURCE_ARTIFACT_HASH_ATTRIBUTE))
        .ok()
        .map(|index| inspection.attributes[index])
        .is_none_or(|attribute| {
            attribute.ty != "string" || attribute.value != expected_source_hex.as_bytes()
        })
    {
        return Err(CinematicFinalizationDivergenceCode::SourceIdentity);
    }
    Ok(())
}

fn decode_authority(
    bytes: &[u8],
    expected_identity: ContentHash,
    max_bytes: u64,
    cx: &Cx<'_>,
) -> Result<
    CinematicAuthorityRecord,
    (
        CinematicFinalizationDisposition,
        CinematicFinalizationDivergenceCode,
    ),
> {
    if cx.checkpoint().is_err() {
        return Err((
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ));
    }
    if bytes.is_empty() {
        return Err((
            CinematicFinalizationDisposition::Incomplete,
            CinematicFinalizationDivergenceCode::AuthorityMissing,
        ));
    }
    if u64::try_from(bytes.len()).unwrap_or(u64::MAX)
        > max_bytes.min(MAX_AUTHORITY_RECORD_WIRE_BYTES)
    {
        return Err((
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
        ));
    }
    let record = CinematicAuthorityRecord::from_canonical_bytes(bytes).map_err(|error| {
        (
            if matches!(error, CinematicAuthorityError::UnsupportedSchemaVersion(_)) {
                CinematicFinalizationDisposition::Incompatible
            } else {
                CinematicFinalizationDisposition::Corrupt
            },
            CinematicFinalizationDivergenceCode::AuthorityCodec,
        )
    })?;
    if cx.checkpoint().is_err() {
        return Err((
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ));
    }
    let canonical = record.canonical_bytes();
    if cx.checkpoint().is_err() {
        return Err((
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ));
    }
    let actual_identity = record.identity();
    if cx.checkpoint().is_err() {
        return Err((
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        ));
    }
    if canonical != bytes || actual_identity != expected_identity {
        return Err((
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::AuthorityIdentity,
        ));
    }
    Ok(record)
}

fn verify_frame_authority(
    plan: &CinematicFinalizationPlan,
    entry: &fs_img::FrameArtifactEntry,
    artifact_identity: ContentHash,
    record: &CinematicAuthorityRecord,
) -> Result<(), CinematicFinalizationDivergenceCode> {
    let key = entry.descriptor().key();
    let (kind, class, unit, source, transform) = match key.role() {
        FrameArtifactRole::RawMaster => {
            let expected = plan
                .raw_frames
                .get(&key)
                .ok_or(CinematicFinalizationDivergenceCode::SequenceInventory)?;
            (
                CinematicArtifactKind::RenderEstimate,
                CinematicAuthorityClass::MonteCarloRender,
                CinematicUnitContract::SpectralRadianceSi,
                expected.authority_source_identity,
                expected.authority_transform_identity,
            )
        }
        FrameArtifactRole::DenoisedIntermediate => (
            CinematicArtifactKind::Visualization,
            CinematicAuthorityClass::VisualizationDerivative,
            CinematicUnitContract::SpectralRadianceSi,
            entry
                .source_content_hash()
                .ok_or(CinematicFinalizationDivergenceCode::SourceIdentity)?,
            plan.image_pipeline_identity,
        ),
        FrameArtifactRole::DisplayPreview | FrameArtifactRole::ScientificOverlay => (
            CinematicArtifactKind::Visualization,
            CinematicAuthorityClass::VisualizationDerivative,
            CinematicUnitContract::DisplayEncoded,
            entry
                .source_content_hash()
                .ok_or(CinematicFinalizationDivergenceCode::SourceIdentity)?,
            plan.image_pipeline_identity,
        ),
    };
    if record.artifact_identity() != artifact_identity {
        return Err(CinematicFinalizationDivergenceCode::AuthorityIdentity);
    }
    if record.source_identity() != source {
        return Err(CinematicFinalizationDivergenceCode::SourceIdentity);
    }
    if record.configuration_identity() != plan.configuration_identity
        || record.configuration_version() != u32::from(CINEMATIC_CONFIG_SCHEMA_VERSION)
    {
        return Err(CinematicFinalizationDivergenceCode::ConfigurationIdentity);
    }
    if record.artifact_kind() != kind
        || record.authority_class() != class
        || record.unit_contract() != unit
        || record.transform_identity() != transform
    {
        return Err(CinematicFinalizationDivergenceCode::AuthoritySemantics);
    }
    let expected_clock = CinematicClock::try_new(
        CinematicClockDomain::Video,
        plan.frames_per_second,
        1,
        i64::try_from(key.frame_index()).unwrap_or(i64::MAX),
        i64::try_from(key.frame_index().saturating_add(1)).unwrap_or(i64::MAX),
    )
    .map_err(|_| CinematicFinalizationDivergenceCode::AuthoritySemantics)?;
    if record.clock() != expected_clock {
        return Err(CinematicFinalizationDivergenceCode::AudioVideoDuration);
    }
    let disposition_matches = match (key.role(), record.transform_disposition()) {
        (FrameArtifactRole::RawMaster, CinematicTransformDisposition::MonteCarloEstimator) => true,
        (
            FrameArtifactRole::DenoisedIntermediate
            | FrameArtifactRole::DisplayPreview
            | FrameArtifactRole::ScientificOverlay,
            CinematicTransformDisposition::BiasedVisualization(_),
        ) => true,
        _ => false,
    };
    if !disposition_matches
        || required_no_claims(class)
            .iter()
            .any(|claim| record.no_claims().binary_search(claim).is_err())
    {
        return Err(CinematicFinalizationDivergenceCode::AuthoritySemantics);
    }
    Ok(())
}

fn verify_audio_authority(
    plan: &CinematicFinalizationPlan,
    manifest: &DecodedAudioManifest,
    record: &CinematicAuthorityRecord,
) -> Result<(), CinematicFinalizationDivergenceCode> {
    let class = CinematicAuthorityClass::Sound(manifest.authority);
    let expected_clock = CinematicClock::try_new(
        CinematicClockDomain::Audio,
        plan.audio_sample_rate_hz,
        1,
        0,
        i64::try_from(plan.total_audio_sample_frames).unwrap_or(i64::MAX),
    )
    .map_err(|_| CinematicFinalizationDivergenceCode::AuthoritySemantics)?;
    if record.artifact_identity() != manifest.wav.wav_identity {
        return Err(CinematicFinalizationDivergenceCode::AuthorityIdentity);
    }
    if record.source_identity() != manifest.source_signal_identity {
        return Err(CinematicFinalizationDivergenceCode::SourceIdentity);
    }
    if record.configuration_identity() != plan.configuration_identity
        || record.configuration_version() != u32::from(CINEMATIC_CONFIG_SCHEMA_VERSION)
    {
        return Err(CinematicFinalizationDivergenceCode::ConfigurationIdentity);
    }
    if record.artifact_kind() != CinematicArtifactKind::Audio
        || record.authority_class() != class
        || record.transform_identity() != plan.sound_receipt.configuration_identity
        || record.unit_contract() != CinematicUnitContract::DigitalAudioFullScale
        || record.clock() != expected_clock
        || record.acoustic_calibration() != plan.expected_acoustic_calibration
        || !matches!(
            record.transform_disposition(),
            CinematicTransformDisposition::SoundSynthesis(_)
        )
        || required_no_claims(class)
            .iter()
            .any(|claim| record.no_claims().binary_search(claim).is_err())
    {
        return Err(CinematicFinalizationDivergenceCode::AuthoritySemantics);
    }
    Ok(())
}

/// Return the first bad marker and whether that marker is a shot cut.
fn verify_alignment(
    plan: &CinematicFinalizationPlan,
    alignment: &AudioVideoAlignment,
    cx: &Cx<'_>,
) -> Result<Option<(u32, bool)>, ()> {
    let ratio = plan.audio_sample_rate_hz / plan.frames_per_second;
    let Some(expected_markers) = plan.total_video_frames.checked_add(1) else {
        return Ok(Some((0, false)));
    };
    if alignment.audio_frames_per_video_frame != ratio
        || alignment.endpoint_drift_audio_frames != 0
        || u32::try_from(alignment.markers.len()).ok() != Some(expected_markers)
    {
        return Ok(Some((0, false)));
    }
    for frame in 0..=plan.total_video_frames {
        if frame % 1_024 == 0 && cx.checkpoint().is_err() {
            return Err(());
        }
        let Some(marker) = alignment.markers.get(frame as usize) else {
            return Ok(Some((
                frame,
                plan.cut_frame_boundaries.binary_search(&frame).is_ok(),
            )));
        };
        let Some(audio_offset) = u64::from(frame).checked_mul(u64::from(ratio)) else {
            return Ok(Some((
                frame,
                plan.cut_frame_boundaries.binary_search(&frame).is_ok(),
            )));
        };
        let Ok(audio_tick) = i64::try_from(audio_offset) else {
            return Ok(Some((
                frame,
                plan.cut_frame_boundaries.binary_search(&frame).is_ok(),
            )));
        };
        let valid = marker.video_tick == i64::from(frame)
            && marker.audio_tick == audio_tick
            && marker.audio_frame_offset == audio_offset;
        if !valid {
            return Ok(Some((
                frame,
                plan.cut_frame_boundaries.binary_search(&frame).is_ok(),
            )));
        }
    }
    Ok(None)
}

fn resampled_event_eq(actual: &ResampledAudioEvent, expected: &ResampledAudioEvent) -> bool {
    audio_excitation_event_eq(&actual.source, &expected.source)
        && actual.requested_sample_position.to_bits()
            == expected.requested_sample_position.to_bits()
        && actual.left_frame_offset == expected.left_frame_offset
        && actual.right_frame_offset == expected.right_frame_offset
        && actual.left_weight.to_bits() == expected.left_weight.to_bits()
        && actual.right_weight.to_bits() == expected.right_weight.to_bits()
        && actual.centroid_error_frames.to_bits() == expected.centroid_error_frames.to_bits()
        && actual.bracket_start_sample_position.to_bits()
            == expected.bracket_start_sample_position.to_bits()
        && actual.bracket_end_sample_position.to_bits()
            == expected.bracket_end_sample_position.to_bits()
}

fn audio_excitation_event_eq(
    actual: &crate::audio_excitation::AudioExcitationEvent,
    expected: &crate::audio_excitation::AudioExcitationEvent,
) -> bool {
    actual.source_sample_index == expected.source_sample_index
        && actual.kind == expected.kind
        && actual.time_s.to_bits() == expected.time_s.to_bits()
        && actual.bracket_start_s.to_bits() == expected.bracket_start_s.to_bits()
        && actual.bracket_end_s.to_bits() == expected.bracket_end_s.to_bits()
        && actual.measure == expected.measure
        && modal_values_eq(actual.physical_impulse_n_s, expected.physical_impulse_n_s)
        && match (actual.artistic, expected.artistic) {
            (None, None) => true,
            (Some(actual), Some(expected)) => {
                actual.stream_identity == expected.stream_identity
                    && modal_values_eq(actual.impulse_n_s, expected.impulse_n_s)
            }
            _ => false,
        }
}

fn modal_values_eq(
    actual: crate::modal_synthesis::ModalComponentValues,
    expected: crate::modal_synthesis::ModalComponentValues,
) -> bool {
    actual.disc.to_bits() == expected.disc.to_bits()
        && actual.glass_plate.to_bits() == expected.glass_plate.to_bits()
        && actual.base_assembly.to_bits() == expected.base_assembly.to_bits()
}

fn repairs_for(
    disposition: CinematicFinalizationDisposition,
    code: CinematicFinalizationDivergenceCode,
) -> Vec<CinematicFinalizationRepair> {
    use CinematicFinalizationDivergenceCode as Code;
    use CinematicFinalizationRepair as Repair;
    if disposition == CinematicFinalizationDisposition::Incomplete {
        return vec![
            Repair::CompleteOrResumeProduction,
            Repair::RestoreExpectedArtifact,
        ];
    }
    match code {
        Code::Cancelled => vec![Repair::RetryInLiveExecutionScope],
        Code::BundleBudgetExceeded => vec![Repair::IncreaseExplicitVerificationBudget],
        Code::SequenceIncomplete | Code::MissingArtifact | Code::AuthorityMissing => {
            vec![
                Repair::CompleteOrResumeProduction,
                Repair::RestoreExpectedArtifact,
            ]
        }
        Code::SequenceDecode | Code::DuplicateArtifact | Code::UnexpectedArtifact => {
            vec![Repair::RestoreCanonicalManifest]
        }
        Code::ArtifactHash
        | Code::ImageStructure
        | Code::ImageDimensions
        | Code::ImageChannels
        | Code::ImageMetadata
        | Code::ImageSampleCount
        | Code::ImagePayload
        | Code::AuthorityCodec
        | Code::WavStructure => vec![Repair::RegenerateArtifactFromPinnedInputs],
        Code::AuthorityIdentity | Code::AudioManifestIdentity => vec![
            Repair::SupplyCorrectExternalIdentityPin,
            Repair::RegenerateArtifactFromPinnedInputs,
        ],
        Code::SequenceInventory
        | Code::AuthoritySemantics
        | Code::ConfigurationIdentity
        | Code::BuildIdentity
        | Code::ProfileIdentity
        | Code::SourceIdentity
        | Code::AudioManifestSemantics
        | Code::AudioVideoDuration
        | Code::SyncMarker
        | Code::CutMarker
        | Code::AudioEvent => vec![Repair::RebuildWithAdmittedConfiguration],
    }
}

#[allow(clippy::too_many_arguments)]
fn report_bytes(
    plan_identity: ContentHash,
    sequence_identity: ContentHash,
    audio_manifest_identity: ContentHash,
    wav_identity: ContentHash,
    alignment_identity: ContentHash,
    event_identity: ContentHash,
    target: CinematicFinalizationTarget,
    disposition: CinematicFinalizationDisposition,
    divergence: Option<&CinematicFinalizationDivergence>,
    repairs: &[CinematicFinalizationRepair],
    verified_frames: u32,
    verified_sync_markers: u32,
    verified_audio_events: u32,
    no_claims: &[CinematicNoClaim],
) -> Vec<u8> {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(REPORT_MAGIC);
    bytes.extend_from_slice(&CINEMATIC_FINALIZATION_REPORT_VERSION.to_le_bytes());
    for identity in [
        plan_identity,
        sequence_identity,
        audio_manifest_identity,
        wav_identity,
        alignment_identity,
        event_identity,
    ] {
        bytes.extend_from_slice(identity.as_bytes());
    }
    bytes.push(target.tag());
    bytes.push(disposition.tag());
    match divergence {
        None => bytes.push(0),
        Some(divergence) => {
            bytes.push(1);
            bytes.push(divergence.code.tag());
            push_coordinate(&mut bytes, &divergence.coordinate);
        }
    }
    bytes.extend_from_slice(
        &u32::try_from(repairs.len())
            .expect("repair vocabulary is bounded by this module")
            .to_le_bytes(),
    );
    bytes.extend(repairs.iter().map(|repair| repair.tag()));
    bytes.extend_from_slice(&verified_frames.to_le_bytes());
    bytes.extend_from_slice(&verified_sync_markers.to_le_bytes());
    bytes.extend_from_slice(&verified_audio_events.to_le_bytes());
    bytes.extend_from_slice(
        &u32::try_from(no_claims.len())
            .expect("no-claim vocabulary is bounded by fs-evidence")
            .to_le_bytes(),
    );
    for claim in no_claims {
        push_string(&mut bytes, claim.code());
    }
    bytes
}

fn push_coordinate(bytes: &mut Vec<u8>, coordinate: &CinematicFinalizationCoordinate) {
    match coordinate {
        CinematicFinalizationCoordinate::Bundle => bytes.push(1),
        CinematicFinalizationCoordinate::Sequence => bytes.push(2),
        CinematicFinalizationCoordinate::Frame { key, relative_path } => {
            bytes.push(3);
            bytes.extend_from_slice(&key.frame_index().to_le_bytes());
            bytes.extend_from_slice(&key.segment_index().to_le_bytes());
            bytes.push(match key.role() {
                FrameArtifactRole::RawMaster => 1,
                FrameArtifactRole::DenoisedIntermediate => 2,
                FrameArtifactRole::DisplayPreview => 3,
                FrameArtifactRole::ScientificOverlay => 4,
            });
            push_string(bytes, relative_path);
        }
        CinematicFinalizationCoordinate::Audio => bytes.push(4),
        CinematicFinalizationCoordinate::SyncMarker(index) => {
            bytes.push(5);
            bytes.extend_from_slice(&index.to_le_bytes());
        }
        CinematicFinalizationCoordinate::AudioEvent(index) => {
            bytes.push(6);
            bytes.extend_from_slice(&index.to_le_bytes());
        }
    }
}

fn push_string(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(
        &u32::try_from(value.len())
            .expect("canonical finalization strings are schema-bounded")
            .to_le_bytes(),
    );
    bytes.extend_from_slice(value.as_bytes());
}

fn hash_raw_expectations(
    writer: &mut PlanIdentityHasher<'_, '_>,
    expectations: &BTreeMap<FrameArtifactKey, RawFrameExpectation>,
) -> Result<(), CinematicFinalizationPlanError> {
    writer.update(
        &u32::try_from(expectations.len())
            .expect("raw expectation count was admitted before identity construction")
            .to_le_bytes(),
    )?;
    for (key, expectation) in expectations {
        writer.update(&key.frame_index().to_le_bytes())?;
        writer.update(&key.segment_index().to_le_bytes())?;
        debug_assert_eq!(key.role(), FrameArtifactRole::RawMaster);
        writer.update(expectation.authority_source_identity.as_bytes())?;
        writer.update(expectation.authority_transform_identity.as_bytes())?;
        writer.update(&expectation.object_palette_entries.to_le_bytes())?;
        writer.update(&expectation.material_palette_entries.to_le_bytes())?;
        match expectation.expected_uniform_spp {
            None => writer.byte(0)?,
            Some(spp) => {
                writer.byte(1)?;
                writer.update(&spp.to_le_bytes())?;
            }
        }
        writer.update(
            &u32::try_from(expectation.required_attributes.len())
                .expect("raw attribute count was admitted before identity construction")
                .to_le_bytes(),
        )?;
        for (name, value) in &expectation.required_attributes {
            writer.string(name)?;
            writer.update(
                &u32::try_from(value.len())
                    .expect("raw attribute bytes were admitted before identity construction")
                    .to_le_bytes(),
            )?;
            writer.update(value)?;
        }
    }
    Ok(())
}

fn hash_resampled_events(
    writer: &mut PlanIdentityHasher<'_, '_>,
    events: &[ResampledAudioEvent],
) -> Result<(), CinematicFinalizationPlanError> {
    writer.update(
        &u32::try_from(events.len())
            .expect("plan admission bounds event count to u32")
            .to_le_bytes(),
    )?;
    for event in events {
        write_resampled_event(event, &mut |chunk| writer.update(chunk))?;
    }
    Ok(())
}

fn push_resampled_event(bytes: &mut Vec<u8>, event: &ResampledAudioEvent) {
    let result: Result<(), core::convert::Infallible> =
        write_resampled_event(event, &mut |chunk| {
            bytes.extend_from_slice(chunk);
            Ok(())
        });
    match result {
        Ok(()) => {}
        Err(never) => match never {},
    }
}

fn write_resampled_event<E>(
    event: &ResampledAudioEvent,
    write: &mut impl FnMut(&[u8]) -> Result<(), E>,
) -> Result<(), E> {
    write(
        &u64::try_from(event.source.source_sample_index)
            .expect("plan admission bounds source sample indices to u64")
            .to_le_bytes(),
    )?;
    write(&[match event.source.kind {
        crate::coupled_runner::ContactTransitionKind::Opening => 1,
        crate::coupled_runner::ContactTransitionKind::Reimpact => 2,
    }])?;
    write(&[match event.source.measure {
        crate::control_stream::ContactEventMeasure::TimingOnly => 1,
    }])?;
    for value in [
        event.source.time_s,
        event.source.bracket_start_s,
        event.source.bracket_end_s,
        event.requested_sample_position,
        event.left_weight,
        event.right_weight,
        event.centroid_error_frames,
        event.bracket_start_sample_position,
        event.bracket_end_sample_position,
    ] {
        write(&value.to_bits().to_le_bytes())?;
    }
    write_modal_values(write, event.source.physical_impulse_n_s)?;
    match event.source.artistic {
        None => write(&[0])?,
        Some(artistic) => {
            write(&[1])?;
            write(artistic.stream_identity.as_bytes())?;
            write_modal_values(write, artistic.impulse_n_s)?;
        }
    }
    for offset in [event.left_frame_offset, event.right_frame_offset] {
        match offset {
            Some(offset) => {
                write(&[1])?;
                write(&offset.to_le_bytes())?;
            }
            None => write(&[0])?,
        }
    }
    Ok(())
}

fn write_modal_values<E>(
    write: &mut impl FnMut(&[u8]) -> Result<(), E>,
    values: crate::modal_synthesis::ModalComponentValues,
) -> Result<(), E> {
    for value in [values.disc, values.glass_plate, values.base_assembly] {
        write(&value.to_bits().to_le_bytes())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_evidence::{
        cinematic::{CINEMATIC_AUTHORITY_SCHEMA_VERSION, CinematicAuthorityInput, SoundAuthority},
        cinematic_budget::CinematicQualityTier,
        cinematic_config::{
            CinematicArtifactRoot, CinematicAssetBinding, CinematicAssetInterpretation,
            CinematicCapabilities, CinematicComponentRef, CinematicComponentRole,
            CinematicConfigInput, CinematicConfigUnits, CinematicMuxRequest,
        },
        cinematic_sound::{
            ListenerFrame, ListenerPose, SOUND_MASTER_SAMPLE_RATE_HZ, SoundAmplitudeReference,
            SoundChannelLayout, SoundExcitationChannel, SoundExcitationControl,
            SoundModalComponent, SoundMode, SoundModeParticipation, SoundModelAssumption,
            SoundRoomResponse, SoundSynthesisInput, SoundTerminalPolicy,
            SoundTrajectoryDisposition,
        },
    };
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_geom::{Point3, Vec3 as GeomVec3};
    use fs_img::{
        Channel, ExrAttribute, FrameArtifactFileState, RegistrationOutcome,
        write_exr_with_attributes, write_png16,
    };
    use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3 as MbdVec3};
    use fs_render::camera::{
        AnimatedCamera, Aperture, CameraKeyframe, CameraProjection, CameraShot, PhysicalCamera,
    };
    use fs_rep_frep::SquatDiscEdgeTreatment;

    use crate::{
        DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, EulerRenderTrajectoryArtifact,
        ExposureEventPolicy, RenderBaseFrame, RenderBaseModeState, RenderChannelAvailability,
        RenderContactBranch, RenderMassProperties, RenderSampleDisposition, RenderTrajectory,
        RenderTrajectoryAuthority, RenderTrajectoryCodecBudget, RenderTrajectoryMetadata,
        RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame,
        audio_artifact::{AudioMasterSource, SoundWavArtifact, StereoSample, WavMetadata},
        audio_excitation::{ArtisticEventExcitation, AudioExcitationEvent},
        control_stream::ContactEventMeasure,
        coupled_runner::{ChannelOwnership, ContactTransitionKind},
        modal_synthesis::ModalComponentValues,
        render_scene_bridge::{
            EulerFrameRequest, EulerPreparedFrame, EulerSceneConfig, EulerTessellationConfig,
            euler_scene_smoke_settings,
        },
        render_sharding::EulerRenderShardLimits,
        specimen::{DiscProfileSpec, ResolvedDiscProfile},
    };

    const VIDEO_FRAMES: u32 = 192;
    const AUDIO_FRAMES: u64 = VIDEO_FRAMES as u64 * 2_000;
    const PRODUCTION_VIDEO_FRAMES: u32 = 240;
    const PRODUCTION_AUDIO_FRAMES: u64 = PRODUCTION_VIDEO_FRAMES as u64 * 2_000;

    fn with_cx<R>(cancelled: bool, operation: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        if cancelled {
            gate.request();
        }
        let pool = ArenaPool::new(ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x4649_4e41_4c49_5a45,
                    kernel_id: 0x4255_4e44,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            operation(&cx)
        })
    }

    fn test_identity(label: &str) -> ContentHash {
        hash_domain(
            "org.frankensim.test.cinematic-finalizer.v1",
            label.as_bytes(),
        )
    }

    fn component(role: CinematicComponentRole, label: &str) -> CinematicComponentRef {
        CinematicComponentRef::try_new(role, test_identity(label), 1).expect("test component")
    }

    #[derive(Clone, Copy)]
    struct ProductionBindings {
        trajectory: CinematicComponentRef,
        timeline: CinematicComponentRef,
        excitation: CinematicComponentRef,
        sound_model: CinematicComponentRef,
        microphone: CinematicComponentRef,
        room: CinematicComponentRef,
    }

    impl ProductionBindings {
        fn with_timeline(self, timeline: CinematicComponentRef) -> Self {
            Self { timeline, ..self }
        }
    }

    fn production_specimen(cx: &Cx<'_>) -> ResolvedDiscProfile {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
        }
        .resolve(7_800.0, cx)
        .expect("production-constructor specimen")
    }

    fn production_mass(specimen: &ResolvedDiscProfile) -> MassProperties {
        MassProperties::new(
            specimen.mass_properties.mass,
            MbdVec3::ZERO,
            MbdVec3::new(
                specimen.mass_properties.principal_inertia.transverse,
                specimen.mass_properties.principal_inertia.transverse,
                specimen.mass_properties.principal_inertia.axial,
            ),
        )
        .expect("resolved production-constructor mass")
    }

    fn production_state() -> RigidBodyState {
        let orientation = UnitQuaternion::from_axis_angle(MbdVec3::new(1.0, 0.0, 0.0), 1.0)
            .expect("tilted production-constructor pose");
        RigidBodyState::new(
            Pose::new(MbdVec3::new(0.0, 0.0, 0.045), orientation)
                .expect("production-constructor pose"),
            MbdVec3::ZERO,
            MbdVec3::ZERO,
        )
        .expect("stationary production-constructor state")
    }

    fn production_sample(
        time_s: f64,
        disposition: RenderSampleDisposition,
        specimen: &ResolvedDiscProfile,
        mass: MassProperties,
        cx: &Cx<'_>,
    ) -> RenderTrajectorySampleInput {
        let state = production_state();
        let orientation = state.pose().orientation();
        let contact = crate::profile_contact_geometry(
            &specimen.chart,
            specimen.mass_properties,
            state.pose(),
            cx,
        )
        .expect("open production-constructor support geometry");
        RenderTrajectorySampleInput {
            interval_start_time_s: 0.0,
            time_s,
            world_frame: RenderWorldFrame::RightHandedZUp,
            units: RenderUnitSystem::SiRadians,
            center_of_mass_world_m: state.pose().position_world(),
            orientation_body_to_world: orientation.components(),
            linear_momentum_world_kg_m_per_s: state.linear_momentum_world(),
            angular_momentum_body_kg_m2_per_s: state.angular_momentum_body(),
            symmetry_axis_world: orientation.rotate_body_to_world(MbdVec3::new(0.0, 0.0, 1.0)),
            contact_branch: RenderContactBranch::Open,
            contact_geometry: None,
            signed_gap_m: contact.contact.gap_m,
            interval_contact_active: false,
            interval_normal_force_n: 0.0,
            contact_transitions: Vec::new(),
            base_mode: Some(RenderBaseModeState {
                displacement_m: 0.0,
                velocity_m_per_s: 0.0,
            }),
            channels: ChannelOwnership::default(),
            mechanical_energy_j: 1.0,
            energy_defect_j: 0.0,
            qois: DerivedEulerQois::from_state(state, mass, 0.0)
                .expect("finite production-constructor QoIs"),
            disposition,
            terminal_event: None,
        }
    }

    fn production_trajectory_artifact(
        specimen: &ResolvedDiscProfile,
        cx: &Cx<'_>,
    ) -> EulerRenderTrajectoryArtifact {
        let mass = production_mass(specimen);
        let first = production_sample(0.0, RenderSampleDisposition::Continue, specimen, mass, cx);
        let last = production_sample(
            10.0,
            RenderSampleDisposition::HorizonCensored,
            specimen,
            mass,
            cx,
        );
        let identities = specimen.content_identities();
        let trajectory = RenderTrajectory::try_new(
            RenderTrajectoryMetadata {
                schema_version: EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
                world_frame: RenderWorldFrame::RightHandedZUp,
                units: RenderUnitSystem::SiRadians,
                specimen_profile_identity: identities.profile,
                specimen_chart_identity: identities.chart,
                mass_properties: RenderMassProperties {
                    identity: identities.mass_properties,
                    properties: mass,
                },
                initial_state: production_state(),
                initial_base_mode: first.base_mode.expect("fixture base mode"),
                base_model_identity: test_identity("production-base"),
                base_frame: RenderBaseFrame {
                    origin_world_m: MbdVec3::ZERO,
                    orientation_base_to_world: UnitQuaternion::IDENTITY,
                },
                model_identity: test_identity("production-model"),
                channel_availability: RenderChannelAvailability::NONE_AVAILABLE,
                configuration_identity: test_identity("production-trajectory-configuration"),
                configuration_fingerprint: 0x4649_4e41_4c34_4b54,
                timestep_s: 10.0,
                producer_version: "cinematic-finalization-constructor-test-v1".to_owned(),
                applicability: "deterministic production-constructor binding fixture".to_owned(),
                no_claims: vec!["fixture does not validate Euler-disc mechanics".to_owned()],
                authority: RenderTrajectoryAuthority::SimulationEvidence,
            },
            vec![first, last],
        )
        .expect("production-constructor trajectory");
        EulerRenderTrajectoryArtifact::try_from_trajectory(
            test_identity("production-campaign"),
            trajectory,
            Vec::new(),
            RenderTrajectoryCodecBudget::DEFAULT,
            cx,
        )
        .expect("production-constructor trajectory artifact")
    }

    fn production_camera() -> AnimatedCamera {
        let eye = Point3::new(0.24, -0.30, 0.18);
        let target = Point3::new(0.0, 0.0, 0.025);
        let physical = PhysicalCamera::try_look_at(
            eye,
            target,
            GeomVec3::new(0.0, 0.0, 1.0),
            CameraProjection::try_half_tangent(0.48).expect("fixture projection"),
            target.delta_from(eye).norm(),
            Aperture::try_circular(0.0).expect("fixture pinhole"),
        )
        .expect("fixture physical camera");
        let bounds = [
            (101, 0.0, 2.5),
            (202, 2.5, 5.0),
            (303, 5.0, 8.0),
            (404, 8.0, 10.0),
        ];
        let shots = bounds
            .into_iter()
            .map(|(shot_id, start_s, end_s)| {
                CameraShot::try_new(
                    shot_id,
                    start_s,
                    end_s,
                    vec![
                        CameraKeyframe::try_new(start_s, physical.clone())
                            .expect("fixture camera keyframe"),
                    ],
                )
                .expect("fixture camera shot")
            })
            .collect();
        AnimatedCamera::try_new(shots).expect("four-cut fixture camera")
    }

    fn production_scene_config() -> EulerSceneConfig {
        let mut config = EulerSceneConfig::reference(production_camera());
        config.tessellation = EulerTessellationConfig {
            azimuthal_segments: 16,
            arc_subdivisions_per_arc: 4,
        };
        config
    }

    fn production_prepared_frames(
        scene: &EulerCinematicScene<'_>,
        brief: &CinematicBrief,
    ) -> Vec<EulerPreparedFrame> {
        (0..brief.total_frames())
            .map(|frame| {
                let window = brief
                    .effective_shutter_window(frame)
                    .expect("admitted fixture shutter window");
                let scale = 1_000_000.0 * 24.0;
                let open_s = window.start_microframes as f64 / scale;
                let close_s = window.end_microframes as f64 / scale;
                let at_shot_start = brief
                    .shot_for_frame(frame)
                    .expect("fixture frame belongs to a shot")
                    .frames()
                    .start()
                    == frame;
                scene
                    .prepare_frame(EulerFrameRequest {
                        frame_time_s: if at_shot_start { open_s } else { close_s },
                        exposure_duration_s: close_s - open_s,
                        convention: if at_shot_start {
                            ShutterConvention::FrontLoaded
                        } else {
                            ShutterConvention::BackLoaded
                        },
                        distribution: ShutterDistribution::UniformCounterV1,
                        event_policy: ExposureEventPolicy::Refuse,
                        cut_side: CutSide::After,
                    })
                    .expect("brief-matched prepared exposure")
            })
            .collect()
    }

    fn production_settings(profile: &CinematicQualityProfile) -> fs_render::tracer::Settings {
        let profile = profile.input();
        let mut settings = euler_scene_smoke_settings(profile.width_pixels, profile.height_pixels);
        settings.spp = profile.spp_floor;
        settings.max_depth = u32::from(profile.max_path_depth);
        settings
    }

    fn production_render_plan(
        scene: &EulerCinematicScene<'_>,
        brief: &CinematicBrief,
        profile: &CinematicQualityProfile,
        inputs: &[EulerRenderFrameInput<'_>],
        cx: &Cx<'_>,
    ) -> EulerUniformRenderPlan {
        let profile_input = profile.input();
        let settings = production_settings(profile);
        let tile_count = u64::from(
            profile_input
                .width_pixels
                .div_ceil(u32::from(profile_input.tile_width)),
        ) * u64::from(
            profile_input
                .height_pixels
                .div_ceil(u32::from(profile_input.tile_height)),
        );
        let maximum_paths = u64::from(profile_input.width_pixels)
            * u64::from(profile_input.height_pixels)
            * u64::from(settings.spp);
        let limits = EulerRenderShardLimits::try_new(
            u64::from(PRODUCTION_VIDEO_FRAMES),
            512,
            1 << 20,
            maximum_paths,
            2 << 30,
            1 << 40,
        )
        .expect("production-constructor sharding limits");
        EulerUniformRenderPlan::try_new(
            scene,
            brief.identity(),
            inputs,
            settings,
            u32::from(profile_input.tile_width),
            u32::from(profile_input.tile_height),
            tile_count,
            settings.spp,
            1,
            limits,
            cx,
        )
        .expect("production-constructor render plan")
    }

    fn production_bindings(
        trajectory_identity: ContentHash,
        brief: &CinematicBrief,
    ) -> ProductionBindings {
        ProductionBindings {
            trajectory: CinematicComponentRef::try_new(
                CinematicComponentRole::Trajectory,
                trajectory_identity,
                1,
            )
            .expect("trajectory binding"),
            timeline: CinematicComponentRef::try_new(
                CinematicComponentRole::Timeline,
                brief.identity(),
                u32::from(CINEMATIC_BRIEF_IDENTITY_VERSION),
            )
            .expect("timeline binding"),
            excitation: component(
                CinematicComponentRole::AudioExcitation,
                "production-excitation",
            ),
            sound_model: component(CinematicComponentRole::SoundModel, "production-sound-model"),
            microphone: component(CinematicComponentRole::Microphone, "production-microphone"),
            room: component(CinematicComponentRole::Room, "production-room"),
        }
    }

    fn production_sound_configuration(
        bindings: ProductionBindings,
        brief: &CinematicBrief,
    ) -> SoundSynthesisConfig {
        let video_clock = CinematicClock::try_new(
            CinematicClockDomain::Video,
            24,
            1,
            0,
            i64::from(brief.total_frames()),
        )
        .expect("production video clock");
        let audio_clock = CinematicClock::try_new(
            CinematicClockDomain::Audio,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            1,
            0,
            i64::try_from(brief.total_audio_sample_frames()).expect("fixture audio clock fits i64"),
        )
        .expect("production audio clock");
        SoundSynthesisConfig::try_admit(SoundSynthesisInput {
            schema_version: SOUND_SYNTHESIS_SCHEMA_VERSION,
            authority: SoundAuthority::PhysicallyInformed,
            trajectory: bindings.trajectory,
            excitation: bindings.excitation,
            sound_model: bindings.sound_model,
            microphone: bindings.microphone,
            room: bindings.room,
            timeline: bindings.timeline,
            video_clock,
            audio_clock,
            channel_layout: SoundChannelLayout::Stereo,
            listener: ListenerPose {
                frame: ListenerFrame::AnimatedCamera,
                position_m: [0.0, 0.0, 0.0],
                forward: [0.0, 0.0, -1.0],
                up: [0.0, 1.0, 0.0],
            },
            excitation_controls: vec![SoundExcitationControl {
                channel: SoundExcitationChannel::ContactNormalForce,
                target_component: SoundModalComponent::Disc,
                source_scale: 1.0,
            }],
            modes: vec![SoundMode {
                mode_id: 1,
                component: SoundModalComponent::Disc,
                frequency_hz: 440.0,
                damping_ratio: 0.02,
                modal_mass_kg: 0.2,
                source_participation: SoundModeParticipation {
                    disc: 1.0,
                    glass_plate: 0.0,
                    base_assembly: 0.0,
                },
                radiation_gain_fs_s_per_m: 0.1,
                material_identity: test_identity("production-sound-material"),
                base_identity: test_identity("production-sound-base"),
            }],
            room_response: SoundRoomResponse::Dry,
            amplitude_reference: SoundAmplitudeReference::DigitalFullScale { headroom_db: 6.0 },
            trajectory_disposition: SoundTrajectoryDisposition::HorizonCensored,
            terminal_policy: SoundTerminalPolicy::FadeAtLastAccepted {
                fade_sample_frames: 240,
            },
            resampler_identity: test_identity("production-resampler"),
            resampler_version: 1,
            filter_identity: test_identity("production-filter"),
            filter_version: 1,
            assumptions: vec![
                SoundModelAssumption::LinearModalSuperposition,
                SoundModelAssumption::TimeInvariantDamping,
                SoundModelAssumption::DeclaredExcitationCompleteness,
                SoundModelAssumption::DeclaredRoomResponse,
            ],
            calibration: None,
        })
        .expect("production sound configuration")
    }

    fn production_asset(
        label: &str,
        interpretation: CinematicAssetInterpretation,
    ) -> CinematicAssetBinding {
        CinematicAssetBinding::from_bytes(
            label.as_bytes(),
            interpretation,
            1,
            format!("/fixture/{label}"),
        )
        .expect("production fixture asset")
    }

    fn production_configuration(
        profile: &CinematicQualityProfile,
        bindings: ProductionBindings,
        seed: u64,
    ) -> CinematicConfig {
        CinematicConfig::try_new(CinematicConfigInput {
            schema_version: CINEMATIC_CONFIG_SCHEMA_VERSION,
            units: Some(CinematicConfigUnits::SiMetersKilogramsSecondsRadians),
            seed: Some(seed),
            capabilities: Some(
                CinematicCapabilities::try_new(
                    CinematicCapabilities::RENDER | CinematicCapabilities::AUDIO,
                )
                .expect("production fixture capabilities"),
            ),
            render_budget_profile: Some(
                CinematicComponentRef::try_new(
                    CinematicComponentRole::RenderBudgetProfile,
                    profile.identity(),
                    u32::from(CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION),
                )
                .expect("render budget binding"),
            ),
            audio_budget_profile: Some(
                CinematicComponentRef::try_new(
                    CinematicComponentRole::AudioBudgetProfile,
                    profile.identity(),
                    u32::from(CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION),
                )
                .expect("audio budget binding"),
            ),
            trajectory: bindings.trajectory,
            timeline: bindings.timeline,
            camera: component(CinematicComponentRole::Camera, "production-camera"),
            scene_geometry: component(
                CinematicComponentRole::SceneGeometry,
                "production-scene-geometry",
            ),
            instance_mapping: component(
                CinematicComponentRole::InstanceMapping,
                "production-instance-mapping",
            ),
            renderer: component(CinematicComponentRole::Renderer, "production-renderer"),
            image_pipeline: component(
                CinematicComponentRole::ImagePipeline,
                "production-image-pipeline",
            ),
            audio_excitation: bindings.excitation,
            sound_model: bindings.sound_model,
            microphone: bindings.microphone,
            room: bindings.room,
            material_assets: vec![production_asset(
                "production-metal",
                CinematicAssetInterpretation::SpectralReflectance,
            )],
            light_assets: vec![production_asset(
                "production-softbox",
                CinematicAssetInterpretation::SpectralEmission,
            )],
            environment_asset: production_asset(
                "production-environment",
                CinematicAssetInterpretation::SpectralEmission,
            ),
            artifact_root: CinematicArtifactRoot::try_new(
                "fixture/euler-finalization".to_owned(),
                "/fixture/output".to_owned(),
            )
            .expect("production fixture artifact root"),
            mux_request: CinematicMuxRequest::None,
        })
        .expect("production cinematic configuration")
    }

    fn production_finalization_input<'a, 'scene, 'frame>(
        configuration: &'a CinematicConfig,
        profile: &'a CinematicQualityProfile,
        brief: &'a CinematicBrief,
        render_plan: &'a EulerUniformRenderPlan,
        scene: &'a EulerCinematicScene<'scene>,
        render_frames: &'a [EulerRenderFrameInput<'frame>],
        sound_configuration: &'a SoundSynthesisConfig,
    ) -> EulerCinematicFinalizationPlanInput<'a, 'scene, 'frame> {
        EulerCinematicFinalizationPlanInput {
            configuration,
            quality_profile: profile,
            brief,
            render_plan,
            scene,
            render_frames,
            aov_limits: CinematicAovLimits::default(),
            sequence_limits: FrameSequenceLimits::default(),
            artifact_ceilings: CinematicFrameArtifactCeilings {
                raw_master_bytes: 800 * 1024 * 1024,
                denoised_intermediate_bytes: 100 * 1024 * 1024,
                display_preview_bytes: 40 * 1024 * 1024,
            },
            build_identity: test_identity("production-build"),
            sound_configuration,
            expected_audio_source_signal_identity: test_identity("production-audio-source"),
            expected_audio_channel_layout_identity: test_identity(
                "production-audio-channel-layout",
            ),
            expected_audio_mix_identity: None,
            expected_audio_events: &[],
            max_audio_events: 1,
        }
    }

    fn sound_configuration() -> SoundSynthesisConfig {
        let video_clock = CinematicClock::try_new(
            CinematicClockDomain::Video,
            24,
            1,
            0,
            i64::from(VIDEO_FRAMES),
        )
        .expect("video clock");
        let audio_clock = CinematicClock::try_new(
            CinematicClockDomain::Audio,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            1,
            0,
            i64::try_from(AUDIO_FRAMES).expect("fixture audio clock fits i64"),
        )
        .expect("audio clock");
        SoundSynthesisConfig::try_admit(SoundSynthesisInput {
            schema_version: SOUND_SYNTHESIS_SCHEMA_VERSION,
            authority: SoundAuthority::PhysicallyInformed,
            trajectory: component(CinematicComponentRole::Trajectory, "trajectory"),
            excitation: component(CinematicComponentRole::AudioExcitation, "excitation"),
            sound_model: component(CinematicComponentRole::SoundModel, "sound-model"),
            microphone: component(CinematicComponentRole::Microphone, "microphone"),
            room: component(CinematicComponentRole::Room, "room"),
            timeline: component(CinematicComponentRole::Timeline, "timeline"),
            video_clock,
            audio_clock,
            channel_layout: SoundChannelLayout::Stereo,
            listener: ListenerPose {
                frame: ListenerFrame::AnimatedCamera,
                position_m: [0.0, 0.0, 0.0],
                forward: [0.0, 0.0, -1.0],
                up: [0.0, 1.0, 0.0],
            },
            excitation_controls: vec![SoundExcitationControl {
                channel: SoundExcitationChannel::ContactNormalForce,
                target_component: SoundModalComponent::Disc,
                source_scale: 1.0,
            }],
            modes: vec![SoundMode {
                mode_id: 1,
                component: SoundModalComponent::Disc,
                frequency_hz: 440.0,
                damping_ratio: 0.02,
                modal_mass_kg: 0.2,
                source_participation: SoundModeParticipation {
                    disc: 1.0,
                    glass_plate: 0.0,
                    base_assembly: 0.0,
                },
                radiation_gain_fs_s_per_m: 0.1,
                material_identity: test_identity("material"),
                base_identity: test_identity("base"),
            }],
            room_response: SoundRoomResponse::Dry,
            amplitude_reference: SoundAmplitudeReference::DigitalFullScale { headroom_db: 6.0 },
            trajectory_disposition: SoundTrajectoryDisposition::HorizonCensored,
            terminal_policy: SoundTerminalPolicy::FadeAtLastAccepted {
                fade_sample_frames: 240,
            },
            resampler_identity: test_identity("resampler"),
            resampler_version: 1,
            filter_identity: test_identity("filter"),
            filter_version: 1,
            assumptions: vec![
                SoundModelAssumption::LinearModalSuperposition,
                SoundModelAssumption::TimeInvariantDamping,
                SoundModelAssumption::DeclaredExcitationCompleteness,
                SoundModelAssumption::DeclaredRoomResponse,
            ],
            calibration: None,
        })
        .expect("sound fixture")
    }

    fn event_receipt() -> ResampledAudioEvent {
        ResampledAudioEvent {
            source: AudioExcitationEvent {
                source_sample_index: 7,
                kind: ContactTransitionKind::Reimpact,
                time_s: 1.0,
                bracket_start_s: 0.999,
                bracket_end_s: 1.001,
                measure: ContactEventMeasure::TimingOnly,
                physical_impulse_n_s: ModalComponentValues::ZERO,
                artistic: Some(ArtisticEventExcitation {
                    stream_identity: test_identity("event-stream"),
                    impulse_n_s: ModalComponentValues {
                        disc: 0.01,
                        glass_plate: 0.0,
                        base_assembly: 0.0,
                    },
                }),
            },
            requested_sample_position: 48_000.0,
            left_frame_offset: Some(48_000),
            right_frame_offset: None,
            left_weight: 1.0,
            right_weight: 0.0,
            centroid_error_frames: 0.0,
            bracket_start_sample_position: 47_952.0,
            bracket_end_sample_position: 48_048.0,
        }
    }

    fn alignment() -> AudioVideoAlignment {
        AudioVideoAlignment {
            audio_frames_per_video_frame: 2_000,
            markers: (0..=VIDEO_FRAMES)
                .map(|frame| AudioVideoSyncMarker {
                    video_tick: i64::from(frame),
                    audio_tick: i64::from(frame) * 2_000,
                    audio_frame_offset: u64::from(frame) * 2_000,
                })
                .collect(),
            endpoint_drift_audio_frames: 0,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn authority(
        artifact_kind: CinematicArtifactKind,
        authority_class: CinematicAuthorityClass,
        artifact_identity: ContentHash,
        source_identity: ContentHash,
        transform_identity: ContentHash,
        transform_name: &str,
        configuration_identity: ContentHash,
        unit_contract: CinematicUnitContract,
        clock: CinematicClock,
        transform_disposition: CinematicTransformDisposition,
    ) -> (Vec<u8>, ContentHash) {
        let record = CinematicAuthorityRecord::try_new(CinematicAuthorityInput {
            schema_version: CINEMATIC_AUTHORITY_SCHEMA_VERSION,
            artifact_kind,
            authority_class,
            artifact_identity,
            source_identity,
            transform_identity,
            transform_name: transform_name.to_owned(),
            configuration_identity,
            configuration_version: u32::from(CINEMATIC_CONFIG_SCHEMA_VERSION),
            unit_contract,
            clock,
            transform_disposition,
            no_claims: required_no_claims(authority_class).to_vec(),
            acoustic_calibration: None,
        })
        .expect("authority fixture");
        (record.canonical_bytes(), record.identity())
    }

    fn rgb_channels(sample_type: FrameChannelType) -> Vec<FrameChannel> {
        ["R", "G", "B"]
            .into_iter()
            .map(|name| FrameChannel::try_new(name, sample_type).expect("RGB channel"))
            .collect()
    }

    #[derive(Clone)]
    struct OwnedFrame {
        relative_path: String,
        bytes: Vec<u8>,
        authority_bytes: Vec<u8>,
        authority_identity: ContentHash,
    }

    #[derive(Clone)]
    struct BundleData {
        sequence_bytes: Vec<u8>,
        sequence_identity: ContentHash,
        frames: Vec<OwnedFrame>,
        manifest_bytes: Vec<u8>,
        manifest_identity: ContentHash,
        wav_bytes: Vec<u8>,
        wav_identity: ContentHash,
        audio_authority_bytes: Vec<u8>,
        audio_authority_identity: ContentHash,
        alignment_bytes: Vec<u8>,
        alignment_identity: ContentHash,
        event_bytes: Vec<u8>,
        event_identity: ContentHash,
    }

    struct Fixture {
        plan: CinematicFinalizationPlan,
        data: BundleData,
        alignment: AudioVideoAlignment,
        events: Vec<ResampledAudioEvent>,
        limits: CinematicFinalizationLimits,
        scene_configuration_identity: ContentHash,
    }

    fn build_fixture(cx: &Cx<'_>) -> Fixture {
        let configuration_identity = test_identity("composition");
        let image_configuration_identity = test_identity("image-config");
        let build_identity = test_identity("build");
        let profile_identity = test_identity("profile");
        let brief_identity = test_identity("brief");
        let render_plan_identity = test_identity("render-plan");
        let image_pipeline_identity = test_identity("image-pipeline");
        let scene_configuration_identity = test_identity("scene-config");
        let raw_source_identity = test_identity("raw-source");
        let raw_transform_identity = test_identity("raw-transform");

        let raw_attribute = ExrAttribute {
            name: "fixture.raw".to_owned(),
            ty: "string".to_owned(),
            value: b"independent".to_vec(),
        };
        let raw_bytes = write_exr_with_attributes(
            1,
            1,
            &[
                Channel {
                    name: "R".to_owned(),
                    ty: PixelType::Float,
                    data: vec![0.25],
                },
                Channel {
                    name: "G".to_owned(),
                    ty: PixelType::Float,
                    data: vec![0.5],
                },
                Channel {
                    name: "B".to_owned(),
                    ty: PixelType::Float,
                    data: vec![0.75],
                },
            ],
            std::slice::from_ref(&raw_attribute),
        )
        .expect("raw EXR");
        let raw_hash = FrameArtifactFileState::from_bytes(&raw_bytes)
            .expect("raw state")
            .content_hash();
        let denoised_bytes = write_exr_with_attributes(
            1,
            1,
            &[
                Channel {
                    name: "R".to_owned(),
                    ty: PixelType::Float,
                    data: vec![0.3],
                },
                Channel {
                    name: "G".to_owned(),
                    ty: PixelType::Float,
                    data: vec![0.5],
                },
                Channel {
                    name: "B".to_owned(),
                    ty: PixelType::Float,
                    data: vec![0.7],
                },
            ],
            &[ExrAttribute {
                name: SOURCE_ARTIFACT_HASH_ATTRIBUTE.to_owned(),
                ty: "string".to_owned(),
                value: raw_hash.to_hex().into_bytes(),
            }],
        )
        .expect("denoised EXR");
        let denoised_hash = FrameArtifactFileState::from_bytes(&denoised_bytes)
            .expect("denoised state")
            .content_hash();
        let preview_bytes =
            write_png16(1, 1, PngColor::Rgb, &[1_000, 2_000, 3_000]).expect("preview PNG");
        let preview_hash = FrameArtifactFileState::from_bytes(&preview_bytes)
            .expect("preview state")
            .content_hash();

        let raw_descriptor = FrameArtifactDescriptor::try_new(
            0,
            0,
            FrameArtifactRole::RawMaster,
            0.0,
            FrameArtifactFormat::OpenExr,
            1,
            1,
            rgb_channels(FrameChannelType::Float32),
            FrameSamplingStats::Uniform { spp: 4 },
        )
        .expect("raw descriptor");
        let raw_key = raw_descriptor.key();
        let denoised_descriptor = FrameArtifactDescriptor::try_new(
            0,
            0,
            FrameArtifactRole::DenoisedIntermediate,
            0.0,
            FrameArtifactFormat::OpenExr,
            1,
            1,
            rgb_channels(FrameChannelType::Float32),
            FrameSamplingStats::Uniform { spp: 4 },
        )
        .expect("denoised descriptor");
        let denoised_key = denoised_descriptor.key();
        let preview_descriptor = FrameArtifactDescriptor::try_new(
            0,
            0,
            FrameArtifactRole::DisplayPreview,
            0.0,
            FrameArtifactFormat::Png16,
            1,
            1,
            rgb_channels(FrameChannelType::Uint16),
            FrameSamplingStats::Uniform { spp: 4 },
        )
        .expect("preview descriptor");
        let expected = vec![
            ExpectedFrameArtifact::try_new(raw_descriptor.clone(), 1 << 20, None)
                .expect("raw expectation"),
            ExpectedFrameArtifact::try_new(
                denoised_descriptor.clone(),
                denoised_bytes.len() as u64,
                Some(raw_key),
            )
            .expect("denoised expectation"),
            ExpectedFrameArtifact::try_new(
                preview_descriptor.clone(),
                preview_bytes.len() as u64,
                Some(denoised_key),
            )
            .expect("preview expectation"),
        ];
        let output_bytes = expected.iter().map(ExpectedFrameArtifact::max_bytes).sum();
        let sequence_limits = FrameSequenceLimits::try_new(3, 3, 512, 1 << 20, output_bytes)
            .expect("sequence limits");
        let context = FrameSequenceContext::try_new(
            brief_identity,
            test_identity("trajectory"),
            image_configuration_identity,
            test_identity("scene"),
            build_identity,
            profile_identity,
        )
        .expect("sequence context");
        let expected_sequence =
            FrameSequenceManifest::try_new(context, expected, sequence_limits, output_bytes)
                .expect("expected sequence");
        let path = |key| {
            expected_sequence
                .entries()
                .iter()
                .find(|entry| entry.descriptor().key() == key)
                .expect("fixture path")
                .relative_path()
                .to_owned()
        };
        let raw_path = path(raw_key);
        let denoised_path = path(denoised_key);
        let preview_path = path(preview_descriptor.key());

        let frame_clock =
            CinematicClock::try_new(CinematicClockDomain::Video, 24, 1, 0, 1).expect("frame clock");
        let (raw_authority, raw_authority_identity) = authority(
            CinematicArtifactKind::RenderEstimate,
            CinematicAuthorityClass::MonteCarloRender,
            raw_hash,
            raw_source_identity,
            raw_transform_identity,
            "fixture-render",
            configuration_identity,
            CinematicUnitContract::SpectralRadianceSi,
            frame_clock,
            CinematicTransformDisposition::MonteCarloEstimator,
        );
        let (denoised_authority, denoised_authority_identity) = authority(
            CinematicArtifactKind::Visualization,
            CinematicAuthorityClass::VisualizationDerivative,
            denoised_hash,
            raw_hash,
            image_pipeline_identity,
            "fixture-denoise",
            configuration_identity,
            CinematicUnitContract::SpectralRadianceSi,
            frame_clock,
            CinematicTransformDisposition::BiasedVisualization("fixture-denoise".to_owned()),
        );
        let (preview_authority, preview_authority_identity) = authority(
            CinematicArtifactKind::Visualization,
            CinematicAuthorityClass::VisualizationDerivative,
            preview_hash,
            denoised_hash,
            image_pipeline_identity,
            "fixture-display",
            configuration_identity,
            CinematicUnitContract::DisplayEncoded,
            frame_clock,
            CinematicTransformDisposition::BiasedVisualization("fixture-display".to_owned()),
        );
        let frames = vec![
            OwnedFrame {
                relative_path: raw_path.clone(),
                bytes: raw_bytes.clone(),
                authority_bytes: raw_authority,
                authority_identity: raw_authority_identity,
            },
            OwnedFrame {
                relative_path: denoised_path.clone(),
                bytes: denoised_bytes.clone(),
                authority_bytes: denoised_authority,
                authority_identity: denoised_authority_identity,
            },
            OwnedFrame {
                relative_path: preview_path.clone(),
                bytes: preview_bytes.clone(),
                authority_bytes: preview_authority,
                authority_identity: preview_authority_identity,
            },
        ];

        let mut actual_sequence = expected_sequence.clone();
        for (path, descriptor, bytes, source_hash) in [
            (&raw_path, raw_descriptor, &raw_bytes, None),
            (
                &denoised_path,
                denoised_descriptor,
                &denoised_bytes,
                Some(raw_hash),
            ),
            (
                &preview_path,
                preview_descriptor,
                &preview_bytes,
                Some(denoised_hash),
            ),
        ] {
            assert_eq!(
                actual_sequence
                    .register_artifact_bytes(
                        path,
                        descriptor,
                        profile_identity,
                        bytes,
                        source_hash,
                    )
                    .expect("register fixture frame"),
                RegistrationOutcome::Recorded,
            );
        }
        let observed: BTreeMap<_, _> = frames
            .iter()
            .map(|frame| {
                (
                    frame.relative_path.clone(),
                    FrameArtifactFileState::from_bytes(&frame.bytes).expect("frame state"),
                )
            })
            .collect();
        let seal = actual_sequence
            .finalize_with(|| true, |path| observed.get(path).copied())
            .expect("finalize fixture sequence");

        let sound_configuration = sound_configuration();
        let samples = vec![StereoSample::default(); AUDIO_FRAMES as usize];
        let sound_artifact = SoundWavArtifact::try_build(
            &sound_configuration,
            AudioMasterSource::SpatializedStereo {
                frames: &samples,
                spatialization_identity: test_identity("spatialization"),
                source_synthesis: sound_configuration.receipt(),
            },
            WavSampleEncoding::Float32,
            WavMetadata::default(),
            AudioArtifactBudget::DEFAULT,
            cx,
        )
        .expect("float WAV fixture");
        let manifest = sound_artifact.manifest().clone();
        let wav_bytes = sound_artifact.wav_bytes().to_vec();
        let audio_clock = CinematicClock::try_new(
            CinematicClockDomain::Audio,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            1,
            0,
            i64::try_from(AUDIO_FRAMES).expect("audio frames fit i64"),
        )
        .expect("audio authority clock");
        let audio_class = CinematicAuthorityClass::Sound(manifest.authority());
        let (audio_authority_bytes, audio_authority_identity) = authority(
            CinematicArtifactKind::Audio,
            audio_class,
            manifest.wav().wav_identity(),
            manifest.source_signal_identity(),
            sound_configuration.receipt().configuration_identity,
            "fixture-synthesis",
            configuration_identity,
            CinematicUnitContract::DigitalAudioFullScale,
            audio_clock,
            CinematicTransformDisposition::SoundSynthesis("fixture-synthesis".to_owned()),
        );
        let alignment = alignment();
        let events = vec![event_receipt()];
        let alignment_receipt =
            encode_audio_video_alignment_receipt(&alignment, 1_000, 1 << 20, cx)
                .expect("alignment receipt");
        let event_snapshot =
            encode_resampled_audio_event_receipt(&events, 32, 1 << 20, cx).expect("event receipt");
        let mut required_attributes = BTreeMap::new();
        required_attributes.insert(raw_attribute.name, raw_attribute.value.into());
        let mut raw_frames = BTreeMap::new();
        raw_frames.insert(
            raw_key,
            RawFrameExpectation {
                authority_source_identity: raw_source_identity,
                authority_transform_identity: raw_transform_identity,
                required_attributes,
                object_palette_entries: 0,
                material_palette_entries: 0,
                expected_uniform_spp: None,
            },
        );
        let plan = CinematicFinalizationPlan::from_parts(
            CinematicFinalizationTarget::IntegrityFixture,
            configuration_identity,
            image_configuration_identity,
            build_identity,
            profile_identity,
            brief_identity,
            render_plan_identity,
            image_pipeline_identity,
            scene_configuration_identity,
            expected_sequence,
            raw_frames,
            sound_configuration.receipt(),
            manifest.source_signal_identity(),
            manifest.channel_layout().identity(),
            manifest.mix_identity(),
            None,
            events.clone(),
            VIDEO_FRAMES,
            AUDIO_FRAMES,
            24,
            SOUND_MASTER_SAMPLE_RATE_HZ,
            vec![96],
            cx,
        )
        .expect("finalization fixture plan");
        let manifest_receipt = encode_audio_artifact_manifest_receipt(&manifest, 1 << 16, cx)
            .expect("audio manifest receipt");
        let data = BundleData {
            sequence_bytes: seal.bytes().to_vec(),
            sequence_identity: seal.identity(),
            frames,
            manifest_bytes: manifest_receipt.bytes().to_vec(),
            manifest_identity: manifest_receipt.identity(),
            wav_bytes,
            wav_identity: manifest.wav().wav_identity(),
            audio_authority_bytes,
            audio_authority_identity,
            alignment_bytes: alignment_receipt.bytes().to_vec(),
            alignment_identity: alignment_receipt.identity(),
            event_bytes: event_snapshot.bytes().to_vec(),
            event_identity: event_snapshot.identity(),
        };
        Fixture {
            plan,
            data,
            alignment,
            events,
            limits: CinematicFinalizationLimits {
                frame_sequence: sequence_limits,
                exr: ExrInspectLimits {
                    max_input_bytes: 1 << 20,
                    max_header_bytes: 1 << 16,
                    max_decoded_bytes: 1 << 20,
                    max_metadata_bytes: 1 << 16,
                },
                png: PngInspectLimits {
                    max_input_bytes: 1 << 20,
                    max_decoded_bytes: 1 << 20,
                },
                audio: AudioArtifactBudget::DEFAULT,
                max_audio_manifest_bytes: 1 << 16,
                max_sync_markers: 1_000,
                max_audio_events: 32,
                max_authority_record_bytes: 1 << 16,
                max_bundle_bytes: 16 << 20,
            },
            scene_configuration_identity,
        }
    }

    fn verify_data(
        fixture: &Fixture,
        data: &BundleData,
        limits: CinematicFinalizationLimits,
        cx: &Cx<'_>,
    ) -> CinematicFinalizationReport {
        verify_plan_data(&fixture.plan, data, limits, cx)
    }

    fn verify_plan_data(
        plan: &CinematicFinalizationPlan,
        data: &BundleData,
        limits: CinematicFinalizationLimits,
        cx: &Cx<'_>,
    ) -> CinematicFinalizationReport {
        let frame_views: Vec<_> = data
            .frames
            .iter()
            .map(|frame| CinematicFrameArtifact {
                relative_path: &frame.relative_path,
                bytes: &frame.bytes,
                authority_bytes: &frame.authority_bytes,
                authority_identity: frame.authority_identity,
            })
            .collect();
        verify_cinematic_bundle(
            plan,
            &CinematicBundle {
                sequence_bytes: &data.sequence_bytes,
                sequence_identity: data.sequence_identity,
                frames: &frame_views,
                audio: CinematicAudioArtifact {
                    manifest_bytes: &data.manifest_bytes,
                    manifest_identity: data.manifest_identity,
                    wav_bytes: &data.wav_bytes,
                    wav_identity: data.wav_identity,
                    authority_bytes: &data.audio_authority_bytes,
                    authority_identity: data.audio_authority_identity,
                    alignment_bytes: &data.alignment_bytes,
                    alignment_identity: data.alignment_identity,
                    event_bytes: &data.event_bytes,
                    event_identity: data.event_identity,
                },
            },
            limits,
            cx,
        )
    }

    fn assert_failure(
        report: &CinematicFinalizationReport,
        disposition: CinematicFinalizationDisposition,
        code: CinematicFinalizationDivergenceCode,
    ) {
        assert_eq!(report.disposition(), disposition);
        assert_eq!(
            report
                .first_divergence()
                .expect("failed report has divergence")
                .code,
            code,
        );
        assert!(!report.repairs().is_empty());
    }

    fn reseal_sequence(plan: &CinematicFinalizationPlan, data: &mut BundleData) {
        let states: BTreeMap<_, _> = plan
            .expected_sequence
            .entries()
            .iter()
            .map(|entry| {
                let frame = data
                    .frames
                    .iter()
                    .find(|frame| frame.relative_path == entry.relative_path())
                    .expect("every expected frame has fixture bytes");
                (
                    entry.descriptor().key(),
                    FrameArtifactFileState::from_bytes(&frame.bytes).expect("frame state"),
                )
            })
            .collect();
        let mut sequence = plan.expected_sequence.clone();
        for entry in plan.expected_sequence.entries() {
            let frame = data
                .frames
                .iter()
                .find(|frame| frame.relative_path == entry.relative_path())
                .expect("every expected frame has fixture bytes");
            let source_hash = entry.source().map(|source| {
                states
                    .get(&source)
                    .expect("source precedes derivative")
                    .content_hash()
            });
            assert_eq!(
                sequence
                    .register_artifact_bytes(
                        entry.relative_path(),
                        entry.descriptor().clone(),
                        plan.profile_identity,
                        &frame.bytes,
                        source_hash,
                    )
                    .expect("register hostile frame"),
                RegistrationOutcome::Recorded,
            );
        }
        let observed: BTreeMap<_, _> = plan
            .expected_sequence
            .entries()
            .iter()
            .map(|entry| {
                (
                    entry.relative_path().to_owned(),
                    *states
                        .get(&entry.descriptor().key())
                        .expect("registered state"),
                )
            })
            .collect();
        let seal = sequence
            .finalize_with(|| true, |path| observed.get(path).copied())
            .expect("finalize hostile sequence");
        data.sequence_bytes = seal.bytes().to_vec();
        data.sequence_identity = seal.identity();
    }

    fn replace_raw_frame(fixture: &Fixture, data: &mut BundleData, bytes: Vec<u8>) {
        let entry = fixture
            .plan
            .expected_sequence
            .entries()
            .iter()
            .find(|entry| entry.descriptor().key().role() == FrameArtifactRole::RawMaster)
            .expect("raw fixture entry");
        let expected = fixture
            .plan
            .raw_frames
            .get(&entry.descriptor().key())
            .expect("raw expectation");
        let state = FrameArtifactFileState::from_bytes(&bytes).expect("replacement raw state");
        let clock = CinematicClock::try_new(
            CinematicClockDomain::Video,
            fixture.plan.frames_per_second,
            1,
            i64::try_from(entry.descriptor().key().frame_index()).expect("frame index"),
            i64::try_from(entry.descriptor().key().frame_index() + 1).expect("frame endpoint"),
        )
        .expect("replacement frame clock");
        let (authority_bytes, authority_identity) = authority(
            CinematicArtifactKind::RenderEstimate,
            CinematicAuthorityClass::MonteCarloRender,
            state.content_hash(),
            expected.authority_source_identity,
            expected.authority_transform_identity,
            "hostile-fixture-render",
            fixture.plan.configuration_identity,
            CinematicUnitContract::SpectralRadianceSi,
            clock,
            CinematicTransformDisposition::MonteCarloEstimator,
        );
        let frame = data
            .frames
            .iter_mut()
            .find(|frame| frame.relative_path == entry.relative_path())
            .expect("raw fixture bytes");
        frame.bytes = bytes;
        frame.authority_bytes = authority_bytes;
        frame.authority_identity = authority_identity;
        reseal_sequence(&fixture.plan, data);
    }

    #[test]
    fn g0_production_constructor_derives_cut_shots_and_refuses_cross_bound_inputs() {
        with_cx(false, |cx| {
            let brief = CinematicBrief::euler_disc_v1().expect("reference production brief");
            assert_eq!(brief.total_frames(), PRODUCTION_VIDEO_FRAMES);
            assert_eq!(brief.total_audio_sample_frames(), PRODUCTION_AUDIO_FRAMES);
            let profile = CinematicQualityProfile::canonical(CinematicQualityTier::Final4k)
                .expect("canonical Final4K profile");
            let specimen = production_specimen(cx);
            let artifact = production_trajectory_artifact(&specimen, cx);
            let scene =
                EulerCinematicScene::try_build(&artifact, &specimen, production_scene_config(), cx)
                    .expect("production-constructor scene");
            let prepared = production_prepared_frames(&scene, &brief);
            let render_frames: Vec<_> = prepared
                .iter()
                .enumerate()
                .map(|(frame, prepared)| {
                    EulerRenderFrameInput::new(
                        u64::try_from(frame).expect("production frame ordinal"),
                        prepared,
                    )
                })
                .collect();
            let render_plan = production_render_plan(&scene, &brief, &profile, &render_frames, cx);
            let bindings = production_bindings(scene.source_trajectory_identity(), &brief);
            let sound = production_sound_configuration(bindings, &brief);
            let configuration =
                production_configuration(&profile, bindings, render_plan.settings().seed);

            let plan = CinematicFinalizationPlan::try_from_euler_disc(
                production_finalization_input(
                    &configuration,
                    &profile,
                    &brief,
                    &render_plan,
                    &scene,
                    &render_frames,
                    &sound,
                ),
                cx,
            )
            .expect("coherent production constructor inputs");

            assert_eq!(plan.target(), CinematicFinalizationTarget::NonFinal);
            assert_eq!(plan.cut_frame_boundaries, vec![60, 120, 192]);
            assert_eq!(plan.raw_frames.len(), PRODUCTION_VIDEO_FRAMES as usize);
            assert_eq!(
                plan.expected_sequence().entries().len(),
                3 * PRODUCTION_VIDEO_FRAMES as usize,
                "Final4K uniform output retains raw, denoised, and preview roles",
            );
            let shot_id = |frame_index: u64| {
                let expectation = plan
                    .raw_frames
                    .iter()
                    .find(|(key, _)| key.frame_index() == frame_index)
                    .map(|(_, expectation)| expectation)
                    .expect("raw expectation for production frame");
                let encoded = expectation
                    .required_attributes
                    .get("frankensim.render.shotId")
                    .expect("scene-derived shot attribute");
                std::str::from_utf8(encoded)
                    .expect("shot identity is UTF-8")
                    .parse::<u64>()
                    .expect("shot identity is a canonical integer")
            };
            assert_eq!(
                [
                    shot_id(0),
                    shot_id(59),
                    shot_id(60),
                    shot_id(119),
                    shot_id(120),
                    shot_id(191),
                    shot_id(192),
                    shot_id(239),
                ],
                [101, 101, 202, 202, 303, 303, 404, 404],
                "hard-cut shot ownership must come from each scene-bound exposure",
            );

            let mut partial_input =
                CinematicQualityProfile::canonical(CinematicQualityTier::Qualification4kFrame)
                    .expect("qualification profile")
                    .input()
                    .clone();
            partial_input.max_path_depth = profile.input().max_path_depth;
            let partial_profile = CinematicQualityProfile::try_new(partial_input)
                .expect("settings-compatible partial profile");
            assert_eq!(
                CinematicFinalizationPlan::try_from_euler_disc(
                    production_finalization_input(
                        &configuration,
                        &partial_profile,
                        &brief,
                        &render_plan,
                        &scene,
                        &render_frames,
                        &sound,
                    ),
                    cx,
                ),
                Err(CinematicFinalizationPlanError::Incompatible(
                    "partial-range A/V finalization requires an explicit range clock",
                )),
                "a one-frame profile cannot borrow the full-master audio clock",
            );

            let mut stale_profile_input = configuration.input().clone();
            stale_profile_input.render_budget_profile = Some(
                CinematicComponentRef::try_new(
                    CinematicComponentRole::RenderBudgetProfile,
                    profile.identity(),
                    u32::from(CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION) + 1,
                )
                .expect("stale-version profile reference"),
            );
            let stale_profile_configuration = CinematicConfig::try_new(stale_profile_input)
                .expect("structurally admitted stale profile reference");
            assert_eq!(
                CinematicFinalizationPlan::try_from_euler_disc(
                    production_finalization_input(
                        &stale_profile_configuration,
                        &profile,
                        &brief,
                        &render_plan,
                        &scene,
                        &render_frames,
                        &sound,
                    ),
                    cx,
                ),
                Err(CinematicFinalizationPlanError::Incompatible(
                    "render budget profile binding",
                )),
            );

            let stale_timeline = CinematicComponentRef::try_new(
                CinematicComponentRole::Timeline,
                brief.identity(),
                u32::from(CINEMATIC_BRIEF_IDENTITY_VERSION) + 1,
            )
            .expect("stale-version timeline reference");
            let stale_timeline_bindings = bindings.with_timeline(stale_timeline);
            let stale_timeline_sound =
                production_sound_configuration(stale_timeline_bindings, &brief);
            let mut stale_timeline_input = configuration.input().clone();
            stale_timeline_input.timeline = stale_timeline;
            let stale_timeline_configuration = CinematicConfig::try_new(stale_timeline_input)
                .expect("structurally admitted stale timeline reference");
            assert_eq!(
                CinematicFinalizationPlan::try_from_euler_disc(
                    production_finalization_input(
                        &stale_timeline_configuration,
                        &profile,
                        &brief,
                        &render_plan,
                        &scene,
                        &render_frames,
                        &stale_timeline_sound,
                    ),
                    cx,
                ),
                Err(CinematicFinalizationPlanError::Incompatible(
                    "composition source identities",
                )),
                "matching config/sound references with a stale brief version still refuse",
            );
        });
    }

    #[test]
    fn g2_real_cross_codec_bundle_passes_and_hostile_twins_fail_at_first_divergence() {
        let fixture = with_cx(false, build_fixture);
        let first = with_cx(false, |cx| {
            verify_data(&fixture, &fixture.data, fixture.limits, cx)
        });
        let second = with_cx(false, |cx| {
            verify_data(&fixture, &fixture.data, fixture.limits, cx)
        });
        assert_eq!(first.disposition(), CinematicFinalizationDisposition::Pass);
        assert_eq!(first.identity(), second.identity());
        assert_eq!(first.canonical_bytes(), second.canonical_bytes());
        assert_eq!(first.verified_frame_artifacts, 3);
        assert_eq!(first.verified_sync_markers, VIDEO_FRAMES + 1);
        assert_eq!(first.verified_audio_events, 1);
        assert!(!first.final_delivery_eligible());

        let mut data = fixture.data.clone();
        data.sequence_identity = test_identity("wrong sequence pin");
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::SequenceDecode,
        );

        let mut data = fixture.data.clone();
        data.sequence_bytes.clear();
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incomplete,
            CinematicFinalizationDivergenceCode::SequenceIncomplete,
        );
        assert_eq!(
            report.repairs().first(),
            Some(&CinematicFinalizationRepair::CompleteOrResumeProduction),
        );

        let snapshot = fixture
            .plan
            .expected_sequence
            .snapshot()
            .expect("incomplete canonical snapshot");
        let mut data = fixture.data.clone();
        data.sequence_bytes = snapshot.bytes().to_vec();
        data.sequence_identity = snapshot.identity();
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incomplete,
            CinematicFinalizationDivergenceCode::SequenceIncomplete,
        );

        let mut data = fixture.data.clone();
        data.frames.pop();
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incomplete,
            CinematicFinalizationDivergenceCode::MissingArtifact,
        );

        let mut data = fixture.data.clone();
        data.frames[0].bytes[0] ^= 1;
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::ArtifactHash,
        );

        let mut data = fixture.data.clone();
        replace_raw_frame(
            &fixture,
            &mut data,
            write_exr_with_attributes(
                1,
                1,
                &[
                    Channel {
                        name: "R".to_owned(),
                        ty: PixelType::Float,
                        data: vec![0.25],
                    },
                    Channel {
                        name: "G".to_owned(),
                        ty: PixelType::Float,
                        data: vec![0.5],
                    },
                    Channel {
                        name: "B".to_owned(),
                        ty: PixelType::Float,
                        data: vec![0.75],
                    },
                ],
                &[
                    ExrAttribute {
                        name: "fixture.raw".to_owned(),
                        ty: "string".to_owned(),
                        value: b"independent".to_vec(),
                    },
                    ExrAttribute {
                        name: "vendor.untrusted".to_owned(),
                        ty: "string".to_owned(),
                        value: b"must-not-bypass-the-exact-oracle".to_vec(),
                    },
                ],
            )
            .expect("hash-consistent raw EXR with an unexpected attribute"),
        );
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::ImageMetadata,
        );

        let mut data = fixture.data.clone();
        replace_raw_frame(
            &fixture,
            &mut data,
            write_exr_with_attributes(
                1,
                1,
                &[
                    Channel {
                        name: "R".to_owned(),
                        ty: PixelType::Float,
                        data: vec![f32::NAN],
                    },
                    Channel {
                        name: "G".to_owned(),
                        ty: PixelType::Float,
                        data: vec![0.5],
                    },
                    Channel {
                        name: "B".to_owned(),
                        ty: PixelType::Float,
                        data: vec![0.75],
                    },
                ],
                &[ExrAttribute {
                    name: "fixture.raw".to_owned(),
                    ty: "string".to_owned(),
                    value: b"independent".to_vec(),
                }],
            )
            .expect("hash-consistent nonfinite raw EXR"),
        );
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::ImagePayload,
        );

        let mut data = fixture.data.clone();
        replace_raw_frame(
            &fixture,
            &mut data,
            write_exr_with_attributes(
                2,
                1,
                &[
                    Channel {
                        name: "R".to_owned(),
                        ty: PixelType::Float,
                        data: vec![0.25, 0.25],
                    },
                    Channel {
                        name: "G".to_owned(),
                        ty: PixelType::Float,
                        data: vec![0.5, 0.5],
                    },
                    Channel {
                        name: "B".to_owned(),
                        ty: PixelType::Float,
                        data: vec![0.75, 0.75],
                    },
                ],
                &[ExrAttribute {
                    name: "fixture.raw".to_owned(),
                    ty: "string".to_owned(),
                    value: b"independent".to_vec(),
                }],
            )
            .expect("hash-consistent wrong-size raw EXR"),
        );
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::ImageDimensions,
        );

        let mut data = fixture.data.clone();
        data.frames[0].authority_identity = test_identity("wrong authority pin");
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::AuthorityIdentity,
        );

        let mut data = fixture.data.clone();
        let original =
            CinematicAuthorityRecord::from_canonical_bytes(&data.frames[1].authority_bytes)
                .expect("derived authority fixture");
        assert!(
            !original.no_claims().is_empty(),
            "hostile twin must actually remove a required disclosure",
        );
        let claim_count = original.no_claims().len();
        let count_offset = data.frames[1]
            .authority_bytes
            .len()
            .checked_sub(claim_count + 2)
            .expect("authority claim suffix");
        data.frames[1].authority_bytes.truncate(count_offset);
        data.frames[1]
            .authority_bytes
            .extend_from_slice(&0_u16.to_le_bytes());
        data.frames[1].authority_identity = test_identity("authority missing no-claims");
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::AuthorityCodec,
        );

        let mut data = fixture.data.clone();
        data.manifest_bytes.clear();
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incomplete,
            CinematicFinalizationDivergenceCode::AudioManifestIdentity,
        );

        let mut data = fixture.data.clone();
        let last = data.manifest_bytes.len() - 1;
        data.manifest_bytes[last] ^= 1;
        data.manifest_identity = hash_domain(AUDIO_MANIFEST_RECEIPT_DOMAIN, &data.manifest_bytes);
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::AudioManifestIdentity,
        );

        let mut data = fixture.data.clone();
        data.wav_bytes.clear();
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incomplete,
            CinematicFinalizationDivergenceCode::WavStructure,
        );
        assert_eq!(
            report.repairs().first(),
            Some(&CinematicFinalizationRepair::CompleteOrResumeProduction),
        );

        let mut data = fixture.data.clone();
        data.wav_identity = test_identity("wrong WAV pin");
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::WavStructure,
        );

        let mut data = fixture.data.clone();
        data.wav_bytes[22..24].copy_from_slice(&1_u16.to_le_bytes());
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::WavStructure,
        );

        let mut data = fixture.data.clone();
        data.wav_bytes[24..28].copy_from_slice(&44_100_u32.to_le_bytes());
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::WavStructure,
        );

        let mut data = fixture.data.clone();
        let last = data.wav_bytes.len() - 1;
        data.wav_bytes[last] ^= 1;
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Corrupt,
            CinematicFinalizationDivergenceCode::WavStructure,
        );

        let mut bad_alignment = fixture.alignment.clone();
        bad_alignment.markers[96].audio_tick += 1;
        let alignment_receipt = with_cx(false, |cx| {
            encode_audio_video_alignment_receipt(&bad_alignment, 1_000, 1 << 20, cx)
                .expect("hostile alignment receipt")
        });
        let mut data = fixture.data.clone();
        data.alignment_bytes = alignment_receipt.bytes().to_vec();
        data.alignment_identity = alignment_receipt.identity();
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::CutMarker,
        );

        let mut bad_events = fixture.events.clone();
        bad_events[0].requested_sample_position += 1.0;
        let event_snapshot = with_cx(false, |cx| {
            encode_resampled_audio_event_receipt(&bad_events, 32, 1 << 20, cx)
                .expect("hostile event receipt")
        });
        let mut data = fixture.data.clone();
        data.event_bytes = event_snapshot.bytes().to_vec();
        data.event_identity = event_snapshot.identity();
        let report = with_cx(false, |cx| verify_data(&fixture, &data, fixture.limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::AudioEvent,
        );

        let mut limits = fixture.limits;
        limits.max_bundle_bytes = 1;
        let report = with_cx(false, |cx| verify_data(&fixture, &fixture.data, limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
        );
        let mut limits = fixture.limits;
        limits.max_audio_manifest_bytes = 1;
        let report = with_cx(false, |cx| verify_data(&fixture, &fixture.data, limits, cx));
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
        );
        let report = with_cx(true, |cx| {
            verify_data(&fixture, &fixture.data, fixture.limits, cx)
        });
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Refused,
            CinematicFinalizationDivergenceCode::Cancelled,
        );

        let mut raw_frames = fixture.plan.raw_frames.clone();
        raw_frames
            .values_mut()
            .next()
            .expect("raw expectation")
            .required_attributes
            .insert("fixture.extra".to_owned(), Arc::from(&b"bound"[..]));
        let changed_plan = with_cx(false, |cx| {
            CinematicFinalizationPlan::from_parts(
                fixture.plan.target,
                fixture.plan.configuration_identity,
                fixture.plan.image_configuration_identity,
                fixture.plan.build_identity,
                fixture.plan.profile_identity,
                fixture.plan.brief_identity,
                fixture.plan.render_plan_identity,
                fixture.plan.image_pipeline_identity,
                fixture.scene_configuration_identity,
                fixture.plan.expected_sequence.clone(),
                raw_frames,
                fixture.plan.sound_receipt,
                fixture.plan.expected_audio_source_signal_identity,
                fixture.plan.expected_audio_channel_layout_identity,
                fixture.plan.expected_audio_mix_identity,
                fixture.plan.expected_acoustic_calibration,
                fixture.plan.expected_audio_events.clone(),
                fixture.plan.total_video_frames,
                fixture.plan.total_audio_sample_frames,
                fixture.plan.frames_per_second,
                fixture.plan.audio_sample_rate_hz,
                fixture.plan.cut_frame_boundaries.clone(),
                cx,
            )
        })
        .expect("changed raw oracle remains structurally valid");
        assert_ne!(fixture.plan.identity(), changed_plan.identity());
        let report = with_cx(false, |cx| {
            verify_plan_data(&changed_plan, &fixture.data, fixture.limits, cx)
        });
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::ImageMetadata,
        );

        let changed_audio_plan = with_cx(false, |cx| {
            CinematicFinalizationPlan::from_parts(
                fixture.plan.target,
                fixture.plan.configuration_identity,
                fixture.plan.image_configuration_identity,
                fixture.plan.build_identity,
                fixture.plan.profile_identity,
                fixture.plan.brief_identity,
                fixture.plan.render_plan_identity,
                fixture.plan.image_pipeline_identity,
                fixture.scene_configuration_identity,
                fixture.plan.expected_sequence.clone(),
                fixture.plan.raw_frames.clone(),
                fixture.plan.sound_receipt,
                test_identity("wrong independently retained audio source"),
                fixture.plan.expected_audio_channel_layout_identity,
                fixture.plan.expected_audio_mix_identity,
                fixture.plan.expected_acoustic_calibration,
                fixture.plan.expected_audio_events.clone(),
                fixture.plan.total_video_frames,
                fixture.plan.total_audio_sample_frames,
                fixture.plan.frames_per_second,
                fixture.plan.audio_sample_rate_hz,
                fixture.plan.cut_frame_boundaries.clone(),
                cx,
            )
        })
        .expect("alternate audio oracle remains structurally valid");
        let report = with_cx(false, |cx| {
            verify_plan_data(&changed_audio_plan, &fixture.data, fixture.limits, cx)
        });
        assert_failure(
            &report,
            CinematicFinalizationDisposition::Incompatible,
            CinematicFinalizationDivergenceCode::AudioManifestSemantics,
        );
    }

    #[test]
    fn g0_g3_receipt_codecs_are_bounded_cancellable_and_reject_every_truncation() {
        let alignment = alignment();
        let events = vec![event_receipt()];
        let alignment_receipt = with_cx(false, |cx| {
            encode_audio_video_alignment_receipt(&alignment, 1_000, 1 << 20, cx).unwrap()
        });
        let event_receipt = with_cx(false, |cx| {
            encode_resampled_audio_event_receipt(&events, 32, 1 << 20, cx).unwrap()
        });
        with_cx(false, |cx| {
            assert_eq!(
                decode_alignment_receipt(
                    alignment_receipt.bytes(),
                    alignment_receipt.identity(),
                    1_000,
                    cx,
                )
                .unwrap(),
                alignment,
            );
            assert!(resampled_event_eq(
                &decode_event_receipt(event_receipt.bytes(), event_receipt.identity(), 32, cx,)
                    .unwrap()[0],
                &events[0],
            ));
        });
        for end in 0..event_receipt.bytes().len() {
            let bytes = &event_receipt.bytes()[..end];
            let identity = hash_domain(EVENT_RECEIPT_DOMAIN, bytes);
            let result = with_cx(false, |cx| decode_event_receipt(bytes, identity, 32, cx));
            assert!(
                matches!(
                    result,
                    Err(ReceiptDecodeError::Missing | ReceiptDecodeError::Corrupt)
                ),
                "truncated event prefix {end} was accepted: {result:?}",
            );
        }
        assert_eq!(
            with_cx(false, |cx| {
                encode_resampled_audio_event_receipt(&events, 0, 1 << 20, cx)
            }),
            Err(CinematicReceiptError::InvalidLimit),
        );
        assert_eq!(
            with_cx(false, |cx| {
                encode_resampled_audio_event_receipt(&events, 32, 1, cx)
            }),
            Err(CinematicReceiptError::BudgetExceeded),
        );
        assert_eq!(
            with_cx(true, |cx| {
                encode_resampled_audio_event_receipt(&events, 32, 1 << 20, cx)
            }),
            Err(CinematicReceiptError::Cancelled),
        );
    }

    #[test]
    fn g0_g3_audio_manifest_receipt_is_independent_bounded_and_strict() {
        let fixture = with_cx(false, build_fixture);
        let decoded = with_cx(false, |cx| {
            decode_audio_manifest_receipt(
                &fixture.data.manifest_bytes,
                fixture.data.manifest_identity,
                fixture.limits.max_audio_manifest_bytes,
                cx,
            )
            .expect("valid persisted audio manifest receipt")
        });
        assert_eq!(
            decoded.source_signal_identity,
            fixture.plan.expected_audio_source_signal_identity,
        );
        assert_eq!(decoded.wav.wav_identity, fixture.data.wav_identity);

        for end in 0..fixture.data.manifest_bytes.len() {
            let bytes = &fixture.data.manifest_bytes[..end];
            let identity = hash_domain(AUDIO_MANIFEST_RECEIPT_DOMAIN, bytes);
            let result = with_cx(false, |cx| {
                decode_audio_manifest_receipt(
                    bytes,
                    identity,
                    fixture.limits.max_audio_manifest_bytes,
                    cx,
                )
            });
            assert!(
                matches!(
                    result,
                    Err(ReceiptDecodeError::Missing | ReceiptDecodeError::Corrupt)
                ),
                "truncated audio-manifest prefix {end} was accepted: {result:?}",
            );
        }
        assert_eq!(
            with_cx(false, |cx| {
                decode_audio_manifest_receipt(
                    &fixture.data.manifest_bytes,
                    fixture.data.manifest_identity,
                    1,
                    cx,
                )
            }),
            Err(ReceiptDecodeError::Budget),
        );
        assert_eq!(
            with_cx(true, |cx| {
                decode_audio_manifest_receipt(
                    &fixture.data.manifest_bytes,
                    fixture.data.manifest_identity,
                    fixture.limits.max_audio_manifest_bytes,
                    cx,
                )
            }),
            Err(ReceiptDecodeError::Cancelled),
        );
        let mut future = fixture.data.manifest_bytes.clone();
        future[8..10].copy_from_slice(&2_u16.to_le_bytes());
        let future_identity = hash_domain(AUDIO_MANIFEST_RECEIPT_DOMAIN, &future);
        assert_eq!(
            with_cx(false, |cx| {
                decode_audio_manifest_receipt(
                    &future,
                    future_identity,
                    fixture.limits.max_audio_manifest_bytes,
                    cx,
                )
            }),
            Err(ReceiptDecodeError::Incompatible),
        );
    }

    #[test]
    fn g0_raw_sample_count_mismatch_is_incompatible_not_a_structural_crash() {
        let bytes = write_exr_with_attributes(
            1,
            1,
            &[Channel {
                name: "samples".to_owned(),
                ty: PixelType::Float,
                data: vec![3.0],
            }],
            &[],
        )
        .unwrap();
        let result = with_cx(false, |cx| {
            verify_raw_sample_count(&bytes, 4, ExrInspectLimits::UNBOUNDED, cx)
        });
        assert_eq!(
            result,
            Err((
                CinematicFinalizationDisposition::Incompatible,
                CinematicFinalizationDivergenceCode::ImageSampleCount,
            )),
        );
    }

    #[test]
    fn g0_wav_structure_classification_separates_noncanonical_from_unsupported() {
        for reason in [
            "noncanonical chunk before fmt",
            "noncanonical chunk before float fact",
            "noncanonical or unknown chunk",
            "chunk after data",
            "non-INFO LIST chunk",
            "noncanonical INFO metadata",
        ] {
            assert_eq!(
                map_audio_artifact_error(&AudioArtifactError::UnsupportedWav(reason)),
                (
                    CinematicFinalizationDisposition::Corrupt,
                    CinematicFinalizationDivergenceCode::WavStructure,
                ),
                "canonical-subset violation {reason:?} is corrupt",
            );
        }
        for reason in [
            "non-RIFF container",
            "non-WAVE RIFF form",
            "audio format tag",
            "fmt extension",
            "non-stereo channel layout",
        ] {
            assert_eq!(
                map_audio_artifact_error(&AudioArtifactError::UnsupportedWav(reason)),
                (
                    CinematicFinalizationDisposition::Incompatible,
                    CinematicFinalizationDivergenceCode::WavStructure,
                ),
                "well-formed unsupported feature {reason:?} is incompatible",
            );
        }
    }

    #[test]
    fn g4_authority_decode_has_a_hard_byte_ceiling_and_cancellation_boundaries() {
        let fixture = with_cx(false, build_fixture);
        let frame = &fixture.data.frames[0];
        with_cx(false, |cx| {
            decode_authority(
                &frame.authority_bytes,
                frame.authority_identity,
                u64::MAX,
                cx,
            )
            .expect("bounded canonical authority record");
        });

        let oversized = vec![0_u8; MAX_AUTHORITY_RECORD_WIRE_BYTES as usize + 1];
        assert_eq!(
            with_cx(false, |cx| {
                decode_authority(&oversized, test_identity("unused pin"), u64::MAX, cx)
            }),
            Err((
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::BundleBudgetExceeded,
            )),
        );
        assert_eq!(
            with_cx(true, |cx| {
                decode_authority(
                    &frame.authority_bytes,
                    frame.authority_identity,
                    u64::MAX,
                    cx,
                )
            }),
            Err((
                CinematicFinalizationDisposition::Refused,
                CinematicFinalizationDivergenceCode::Cancelled,
            )),
        );
    }
}
