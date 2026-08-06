//! Static, checked resource admission for Euler-disc cinematic work.
//!
//! These estimates are conservative planning bounds, not performance claims.
//! Wall-time admission requires an explicit host measurement and every final
//! quality reduction remains a new, visible configuration decision.

use core::fmt;
use fs_blake3::{ContentHash, hash_domain};

/// Version of the canonical quality-profile identity preimage.
pub const CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION: u16 = 1;
/// Domain separating quality profiles from every other content identity.
pub const CINEMATIC_QUALITY_PROFILE_IDENTITY_DOMAIN: &str =
    "org.frankensim.cinematic-quality-profile.v1";

const MAX_DIMENSION: u32 = 7_680;
const MAX_FRAMES: u32 = 12 * 24;
const MAX_SPP: u32 = 4_096;
const MAX_PATH_DEPTH: u16 = 64;
const MAX_SHUTTER_SAMPLES: u16 = 1_024;
const MAX_TILE_EDGE: u16 = 1_024;
const MAX_WORKERS: u16 = 1_024;
const AUDIO_RATE_HZ: u64 = 48_000;
const AUDIO_CHANNELS: u64 = 2;
const AUDIO_BYTES_PER_SAMPLE: u64 = 4;

/// Stable quality/use tier. Qualification is deliberately distinct from final.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum CinematicQualityTier {
    /// Cheap story/camera/synchronization smoke preview.
    StoryboardSmoke = 0,
    /// Full-timeline 1080p creative daily.
    Daily1080p = 1,
    /// One representative 4K frame used to qualify the final configuration.
    Qualification4kFrame = 2,
    /// Full 4K image-master sequence.
    Final4k = 3,
}

/// Denoising disposition. Final raw estimates are never overwritten.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum DenoisePolicy {
    /// No denoising.
    Disabled = 0,
    /// Preview-only denoising; output is a biased visualization.
    PreviewBiased = 1,
    /// Preserve raw masters and emit denoising only as a separate derivative.
    SeparateBiasedDerivative = 2,
}

/// Frozen AOV bundles with a deterministic float-channel count.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum AovPreset {
    /// XYZ beauty only.
    BeautyXyz = 0,
    /// Beauty, albedo, normal, depth, primary coverage, variance, and motion.
    DailyCore = 1,
    /// Daily core plus geometric normal, direct, indirect, emission, object and
    /// material IDs, sample count, and validity.
    FinalDiagnostic = 2,
}

impl AovPreset {
    const fn float_channels(self) -> u64 {
        match self {
            Self::BeautyXyz => 3,
            Self::DailyCore => 14,
            Self::FinalDiagnostic => 30,
        }
    }
}

/// Editable input. Use [`CinematicQualityProfile::canonical`] for frozen tiers.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CinematicQualityProfileInput {
    /// Intended use and anti-downgrade rule set.
    pub tier: CinematicQualityTier,
    /// Raster width.
    pub width_pixels: u32,
    /// Raster height.
    pub height_pixels: u32,
    /// Integer master rate; v1 requires 24.
    pub frames_per_second: u32,
    /// First master-frame index rendered.
    pub first_frame: u32,
    /// Number of contiguous frames rendered.
    pub frame_count: u32,
    /// Minimum accepted samples per pixel.
    pub spp_floor: u32,
    /// Maximum budgeted samples per pixel.
    pub spp_ceiling: u32,
    /// Maximum scattering-event depth.
    pub max_path_depth: u16,
    /// Relative adaptive standard-error target in parts per million.
    pub adaptive_error_ppm: u32,
    /// Denoising and raw-master preservation policy.
    pub denoise_policy: DenoisePolicy,
    /// Required auxiliary-output bundle.
    pub aov_preset: AovPreset,
    /// Explicit temporal samples per pixel sample.
    pub shutter_samples: u16,
    /// Horizontal scheduling tile edge.
    pub tile_width: u16,
    /// Vertical scheduling tile edge.
    pub tile_height: u16,
    /// Maximum participating workers.
    pub worker_limit: u16,
    /// Progressive checkpoint interval.
    pub checkpoint_cadence_spp: u32,
    /// Profile-level peak live-memory ceiling.
    pub memory_ceiling_bytes: u64,
    /// Profile-level per-frame wall-time ceiling.
    pub per_frame_wall_time_ceiling_s: u64,
    /// Profile-level sequence wall-time ceiling.
    pub sequence_wall_time_ceiling_s: u64,
    /// Profile-level aggregate working-storage ceiling.
    pub output_ceiling_bytes: u64,
    /// Disk space that must remain free after admission.
    pub minimum_free_space_reserve_bytes: u64,
}

/// Opaque profile admitted against structural/tier rules.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CinematicQualityProfile(CinematicQualityProfileInput);

impl CinematicQualityProfile {
    /// Validate dimensions, bounds, and tier-specific anti-downgrade rules.
    pub fn try_new(input: CinematicQualityProfileInput) -> Result<Self, CinematicBudgetError> {
        if input.width_pixels == 0
            || input.height_pixels == 0
            || input.width_pixels > MAX_DIMENSION
            || input.height_pixels > MAX_DIMENSION
        {
            return Err(CinematicBudgetError::InvalidDimensions);
        }
        if input.frames_per_second != 24
            || input.frame_count == 0
            || input.frame_count > MAX_FRAMES
            || input
                .first_frame
                .checked_add(input.frame_count)
                .is_none_or(|end| end > MAX_FRAMES)
        {
            return Err(CinematicBudgetError::InvalidFrameRange);
        }
        if input.spp_floor == 0
            || input.spp_ceiling < input.spp_floor
            || input.spp_ceiling > MAX_SPP
        {
            return Err(CinematicBudgetError::InvalidSppRange);
        }
        if input.max_path_depth == 0 || input.max_path_depth > MAX_PATH_DEPTH {
            return Err(CinematicBudgetError::InvalidPathDepth);
        }
        if input.adaptive_error_ppm == 0 || input.adaptive_error_ppm > 1_000_000 {
            return Err(CinematicBudgetError::InvalidAdaptiveTarget);
        }
        if input.shutter_samples == 0 || input.shutter_samples > MAX_SHUTTER_SAMPLES {
            return Err(CinematicBudgetError::InvalidShutterSamples);
        }
        if input.tile_width == 0
            || input.tile_height == 0
            || input.tile_width > MAX_TILE_EDGE
            || input.tile_height > MAX_TILE_EDGE
        {
            return Err(CinematicBudgetError::InvalidTileShape);
        }
        if input.worker_limit == 0 || input.worker_limit > MAX_WORKERS {
            return Err(CinematicBudgetError::InvalidWorkerLimit);
        }
        if input.checkpoint_cadence_spp == 0 || input.checkpoint_cadence_spp > input.spp_ceiling {
            return Err(CinematicBudgetError::InvalidCheckpointCadence);
        }
        if input.memory_ceiling_bytes == 0
            || input.per_frame_wall_time_ceiling_s == 0
            || input.sequence_wall_time_ceiling_s < input.per_frame_wall_time_ceiling_s
            || input.output_ceiling_bytes == 0
            || input.minimum_free_space_reserve_bytes == 0
        {
            return Err(CinematicBudgetError::InvalidResourceEnvelope);
        }
        validate_tier(&input)?;
        Ok(Self(input))
    }

    /// Frozen reproducible reference profile.
    pub fn canonical(tier: CinematicQualityTier) -> Result<Self, CinematicBudgetError> {
        Self::try_new(canonical_input(tier))
    }

    /// Validated fields for configuration identity and admission.
    #[must_use]
    pub const fn input(&self) -> &CinematicQualityProfileInput {
        &self.0
    }

    /// Canonical fixed-width encoding of every admitted profile field.
    ///
    /// This is the byte-level bridge between static budget admission and the
    /// `RenderBudgetProfile` / `AudioBudgetProfile` references carried by a
    /// [`crate::cinematic_config::CinematicConfig`]. Without it, a composition
    /// could name unrelated bytes while the CLI admitted a different profile.
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let input = self.input();
        let mut bytes = Vec::with_capacity(91);
        bytes.extend_from_slice(&CINEMATIC_QUALITY_PROFILE_IDENTITY_VERSION.to_le_bytes());
        bytes.push(input.tier as u8);
        bytes.extend_from_slice(&input.width_pixels.to_le_bytes());
        bytes.extend_from_slice(&input.height_pixels.to_le_bytes());
        bytes.extend_from_slice(&input.frames_per_second.to_le_bytes());
        bytes.extend_from_slice(&input.first_frame.to_le_bytes());
        bytes.extend_from_slice(&input.frame_count.to_le_bytes());
        bytes.extend_from_slice(&input.spp_floor.to_le_bytes());
        bytes.extend_from_slice(&input.spp_ceiling.to_le_bytes());
        bytes.extend_from_slice(&input.max_path_depth.to_le_bytes());
        bytes.extend_from_slice(&input.adaptive_error_ppm.to_le_bytes());
        bytes.push(input.denoise_policy as u8);
        bytes.push(input.aov_preset as u8);
        bytes.extend_from_slice(&input.shutter_samples.to_le_bytes());
        bytes.extend_from_slice(&input.tile_width.to_le_bytes());
        bytes.extend_from_slice(&input.tile_height.to_le_bytes());
        bytes.extend_from_slice(&input.worker_limit.to_le_bytes());
        bytes.extend_from_slice(&input.checkpoint_cadence_spp.to_le_bytes());
        bytes.extend_from_slice(&input.memory_ceiling_bytes.to_le_bytes());
        bytes.extend_from_slice(&input.per_frame_wall_time_ceiling_s.to_le_bytes());
        bytes.extend_from_slice(&input.sequence_wall_time_ceiling_s.to_le_bytes());
        bytes.extend_from_slice(&input.output_ceiling_bytes.to_le_bytes());
        bytes.extend_from_slice(&input.minimum_free_space_reserve_bytes.to_le_bytes());
        bytes
    }

    /// Domain-separated identity consumed by cinematic composition references.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        hash_domain(
            CINEMATIC_QUALITY_PROFILE_IDENTITY_DOMAIN,
            &self.canonical_bytes(),
        )
    }
}

fn validate_tier(input: &CinematicQualityProfileInput) -> Result<(), CinematicBudgetError> {
    let valid = match input.tier {
        CinematicQualityTier::StoryboardSmoke => {
            input.width_pixels <= 960
                && input.height_pixels <= 540
                && input.frame_count <= 48
                && input.spp_ceiling <= 4
                && input.aov_preset == AovPreset::BeautyXyz
                && input.denoise_policy == DenoisePolicy::Disabled
        }
        CinematicQualityTier::Daily1080p => {
            input.width_pixels == 1_920
                && input.height_pixels == 1_080
                && input.first_frame == 0
                && input.frame_count == 240
                && input.spp_floor >= 16
                && input.spp_ceiling <= 64
                && input.aov_preset == AovPreset::DailyCore
        }
        CinematicQualityTier::Qualification4kFrame => {
            input.width_pixels == 3_840
                && input.height_pixels == 2_160
                && input.frame_count == 1
                && input.spp_floor >= 64
                && input.spp_ceiling <= 256
                && input.aov_preset == AovPreset::FinalDiagnostic
        }
        CinematicQualityTier::Final4k => {
            input.width_pixels == 3_840
                && input.height_pixels == 2_160
                && (192..=288).contains(&input.frame_count)
                && input.spp_floor >= 256
                && input.spp_ceiling >= input.spp_floor
                && input.max_path_depth >= 12
                && input.shutter_samples >= 8
                && input.aov_preset == AovPreset::FinalDiagnostic
                && input.denoise_policy == DenoisePolicy::SeparateBiasedDerivative
        }
    };
    if valid {
        Ok(())
    } else {
        Err(CinematicBudgetError::UnsupportedTierCombination(input.tier))
    }
}

/// Host/resource facts supplied to dry admission.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CinematicResourceAvailability {
    /// Memory available to this job.
    pub memory_bytes: u64,
    /// Currently free storage before preserving the profile reserve.
    pub free_storage_bytes: u64,
    /// Wall time granted to the job.
    pub wall_time_available_s: u64,
    /// Workers the host/session can supply.
    pub worker_capacity: u16,
    /// Measured camera-path throughput for this host/configuration.
    pub measured_camera_paths_per_second: u64,
}

/// Checked conservative estimate for one profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CinematicResourceEstimate {
    /// Raster pixels in one frame.
    pub pixels_per_frame: u64,
    /// Upper-bound primary camera paths including shutter sampling.
    pub camera_paths: u64,
    /// Peak film, staging, AOV, and worker-scratch bytes.
    pub live_memory_bytes: u64,
    /// One f64 XYZ accumulation film.
    pub film_bytes: u64,
    /// Transactional/staging copy of the f64 XYZ film.
    pub staging_bytes: u64,
    /// One frame of float AOV storage.
    pub aov_bytes: u64,
    /// One restart checkpoint.
    pub checkpoint_bytes: u64,
    /// Float EXR/AOV master sequence bytes.
    pub exr_sequence_bytes: u64,
    /// Conservative uncompressed RGBA display/PNG sequence bound.
    pub png_sequence_bytes: u64,
    /// Combined master and display image sequence bytes.
    pub image_sequence_bytes: u64,
    /// Stereo 48 kHz float WAV bytes, excluding a negligible header.
    pub wav_bytes: u64,
    /// Temporary display sequence plus WAV used for muxing.
    pub temporary_mux_bytes: u64,
    /// Aggregate masters, two checkpoints, and temporary mux inputs.
    pub total_storage_bytes: u64,
    /// Host-measurement-derived upper-bound time for one frame.
    pub per_frame_wall_time_s: u64,
    /// Host-measurement-derived upper-bound time for the selected range.
    pub sequence_wall_time_s: u64,
}

/// Source of a resource ceiling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResourceLimitSource {
    /// Limit declared by the versioned quality profile.
    ProfileEnvelope,
    /// Limit supplied from current host/session facts.
    HostAvailability,
}

/// Resource whose exact bound was exceeded.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CinematicResourceKind {
    /// Peak resident job memory.
    LiveMemoryBytes,
    /// Aggregate working storage after preserving the free-space reserve.
    StorageBytes,
    /// Time for the slowest estimated frame.
    PerFrameWallTimeSeconds,
    /// Time for all selected frames.
    SequenceWallTimeSeconds,
    /// Concurrent worker capacity.
    Workers,
}

/// One deterministic resource deficit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CinematicResourceDeficit {
    /// Resource that failed admission.
    pub kind: CinematicResourceKind,
    /// Whether the profile or host supplied the failing limit.
    pub source: ResourceLimitSource,
    /// Estimated requirement.
    pub required: u64,
    /// Admitted limit.
    pub available: u64,
}

/// Ranked, explicit repair. Advice never mutates a profile.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum CinematicBudgetRepair {
    /// Supply more host memory.
    IncreaseHostMemory,
    /// Free or supply more storage while preserving the reserve.
    IncreaseFreeStorage,
    /// Grant a longer wall-time envelope.
    ExtendWallTime,
    /// Supply the requested worker capacity.
    IncreaseWorkerCapacity,
    /// Create a new non-final configuration with lower preview SPP.
    LowerPreviewSppWithNewConfiguration,
    /// Create a new non-final configuration with fewer preview AOVs.
    ReducePreviewAovsWithNewConfiguration,
    /// Create a new non-final configuration with a shorter range.
    ShortenRangeWithNewConfiguration,
    /// Version and re-admit a changed profile ceiling.
    RaiseProfileEnvelopeWithNewConfiguration,
}

/// Successful static admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AdmittedCinematicBudget {
    tier: CinematicQualityTier,
    estimate: CinematicResourceEstimate,
}

impl AdmittedCinematicBudget {
    /// Checked resource estimate that passed all profile and host limits.
    #[must_use]
    pub const fn estimate(&self) -> CinematicResourceEstimate {
        self.estimate
    }

    /// Bounded deterministic dry-admission summary for CLI/report consumers.
    #[must_use]
    pub fn summary_json(&self) -> String {
        format!(
            "{{\"schema\":\"cinematic-budget-admission-v1\",\"tier\":\"{}\",\"camera_paths\":{},\"live_memory_bytes\":{},\"total_storage_bytes\":{},\"sequence_wall_time_s\":{},\"verdict\":\"admitted\"}}",
            tier_name(self.tier),
            self.estimate.camera_paths,
            self.estimate.live_memory_bytes,
            self.estimate.total_storage_bytes,
            self.estimate.sequence_wall_time_s,
        )
    }
}

/// Estimate and admit before allocation or render work begins.
pub fn admit_cinematic_budget(
    profile: &CinematicQualityProfile,
    available: CinematicResourceAvailability,
) -> Result<AdmittedCinematicBudget, CinematicBudgetError> {
    if available.measured_camera_paths_per_second == 0 {
        return Err(CinematicBudgetError::MissingThroughputMeasurement);
    }
    let p = profile.input();
    let estimate = estimate(p, available.measured_camera_paths_per_second)?;
    let mut deficits = Vec::new();
    deficit(
        &mut deficits,
        CinematicResourceKind::LiveMemoryBytes,
        ResourceLimitSource::ProfileEnvelope,
        estimate.live_memory_bytes,
        p.memory_ceiling_bytes,
    );
    deficit(
        &mut deficits,
        CinematicResourceKind::LiveMemoryBytes,
        ResourceLimitSource::HostAvailability,
        estimate.live_memory_bytes,
        available.memory_bytes,
    );
    deficit(
        &mut deficits,
        CinematicResourceKind::StorageBytes,
        ResourceLimitSource::ProfileEnvelope,
        estimate.total_storage_bytes,
        p.output_ceiling_bytes,
    );
    let usable_storage = available
        .free_storage_bytes
        .saturating_sub(p.minimum_free_space_reserve_bytes);
    deficit(
        &mut deficits,
        CinematicResourceKind::StorageBytes,
        ResourceLimitSource::HostAvailability,
        estimate.total_storage_bytes,
        usable_storage,
    );
    deficit(
        &mut deficits,
        CinematicResourceKind::PerFrameWallTimeSeconds,
        ResourceLimitSource::ProfileEnvelope,
        estimate.per_frame_wall_time_s,
        p.per_frame_wall_time_ceiling_s,
    );
    deficit(
        &mut deficits,
        CinematicResourceKind::SequenceWallTimeSeconds,
        ResourceLimitSource::ProfileEnvelope,
        estimate.sequence_wall_time_s,
        p.sequence_wall_time_ceiling_s,
    );
    deficit(
        &mut deficits,
        CinematicResourceKind::SequenceWallTimeSeconds,
        ResourceLimitSource::HostAvailability,
        estimate.sequence_wall_time_s,
        available.wall_time_available_s,
    );
    deficit(
        &mut deficits,
        CinematicResourceKind::Workers,
        ResourceLimitSource::HostAvailability,
        u64::from(p.worker_limit),
        u64::from(available.worker_capacity),
    );
    if deficits.is_empty() {
        Ok(AdmittedCinematicBudget {
            tier: p.tier,
            estimate,
        })
    } else {
        Err(CinematicBudgetError::InsufficientResources {
            repairs: ranked_repairs(p.tier, &deficits),
            deficits,
            estimate: Box::new(estimate),
        })
    }
}

fn estimate(
    p: &CinematicQualityProfileInput,
    paths_per_second: u64,
) -> Result<CinematicResourceEstimate, CinematicBudgetError> {
    let pixels = checked_mul(
        u64::from(p.width_pixels),
        u64::from(p.height_pixels),
        "pixels",
    )?;
    let frame_count = u64::from(p.frame_count);
    let paths_per_frame = checked_mul(
        checked_mul(pixels, u64::from(p.spp_ceiling), "paths per frame")?,
        u64::from(p.shutter_samples),
        "shutter paths per frame",
    )?;
    let camera_paths = checked_mul(paths_per_frame, frame_count, "sequence paths")?;
    let film_bytes = checked_mul(pixels, 3 * 8, "f64 XYZ film")?;
    let aov_bytes = checked_mul(
        pixels,
        checked_mul(p.aov_preset.float_channels(), 4, "AOV stride")?,
        "AOV frame",
    )?;
    let tile_pixels = checked_mul(u64::from(p.tile_width), u64::from(p.tile_height), "tile")?;
    let worker_scratch = checked_mul(
        checked_mul(tile_pixels, u64::from(p.worker_limit), "worker tiles")?,
        64,
        "worker scratch",
    )?;
    let staging_bytes = film_bytes;
    let live_memory_bytes = checked_add(
        checked_add(film_bytes, staging_bytes, "transactional film")?,
        checked_add(aov_bytes, worker_scratch, "AOV+worker scratch")?,
        "live memory",
    )?;
    let checkpoint_bytes = checked_add(film_bytes, aov_bytes, "checkpoint")?;
    let display_bytes = checked_mul(pixels, 4, "display frame")?;
    let exr_sequence_bytes = checked_mul(aov_bytes, frame_count, "EXR sequence")?;
    let png_sequence_bytes = checked_mul(display_bytes, frame_count, "PNG sequence")?;
    let image_sequence_bytes = checked_add(
        exr_sequence_bytes,
        png_sequence_bytes,
        "master and display sequences",
    )?;
    let audio_sample_frames = checked_mul(frame_count, AUDIO_RATE_HZ / 24, "audio samples")?;
    let wav_bytes = checked_mul(
        checked_mul(audio_sample_frames, AUDIO_CHANNELS, "audio channels")?,
        AUDIO_BYTES_PER_SAMPLE,
        "WAV bytes",
    )?;
    let temporary_mux_bytes = checked_add(
        checked_mul(display_bytes, frame_count, "temporary display sequence")?,
        wav_bytes,
        "temporary mux",
    )?;
    let total_storage_bytes = checked_add(
        checked_add(image_sequence_bytes, wav_bytes, "masters")?,
        checked_add(
            checked_mul(checkpoint_bytes, 2, "checkpoint pair")?,
            temporary_mux_bytes,
            "working storage",
        )?,
        "total storage",
    )?;
    let per_frame_wall_time_s = ceil_div(paths_per_frame, paths_per_second);
    let sequence_wall_time_s = ceil_div(camera_paths, paths_per_second);
    Ok(CinematicResourceEstimate {
        pixels_per_frame: pixels,
        camera_paths,
        live_memory_bytes,
        film_bytes,
        staging_bytes,
        aov_bytes,
        checkpoint_bytes,
        exr_sequence_bytes,
        png_sequence_bytes,
        image_sequence_bytes,
        wav_bytes,
        temporary_mux_bytes,
        total_storage_bytes,
        per_frame_wall_time_s,
        sequence_wall_time_s,
    })
}

fn checked_mul(left: u64, right: u64, term: &'static str) -> Result<u64, CinematicBudgetError> {
    left.checked_mul(right)
        .ok_or(CinematicBudgetError::ArithmeticOverflow(term))
}

fn checked_add(left: u64, right: u64, term: &'static str) -> Result<u64, CinematicBudgetError> {
    left.checked_add(right)
        .ok_or(CinematicBudgetError::ArithmeticOverflow(term))
}

fn ceil_div(value: u64, divisor: u64) -> u64 {
    value / divisor + u64::from(!value.is_multiple_of(divisor))
}

fn deficit(
    out: &mut Vec<CinematicResourceDeficit>,
    kind: CinematicResourceKind,
    source: ResourceLimitSource,
    required: u64,
    available: u64,
) {
    if required > available {
        out.push(CinematicResourceDeficit {
            kind,
            source,
            required,
            available,
        });
    }
}

fn ranked_repairs(
    tier: CinematicQualityTier,
    deficits: &[CinematicResourceDeficit],
) -> Vec<CinematicBudgetRepair> {
    let mut repairs = Vec::new();
    for deficit in deficits {
        let repair = match (deficit.kind, deficit.source) {
            (CinematicResourceKind::LiveMemoryBytes, ResourceLimitSource::HostAvailability) => {
                CinematicBudgetRepair::IncreaseHostMemory
            }
            (CinematicResourceKind::StorageBytes, ResourceLimitSource::HostAvailability) => {
                CinematicBudgetRepair::IncreaseFreeStorage
            }
            (CinematicResourceKind::Workers, _) => CinematicBudgetRepair::IncreaseWorkerCapacity,
            (CinematicResourceKind::PerFrameWallTimeSeconds, _)
            | (
                CinematicResourceKind::SequenceWallTimeSeconds,
                ResourceLimitSource::HostAvailability,
            ) => CinematicBudgetRepair::ExtendWallTime,
            (_, ResourceLimitSource::ProfileEnvelope) => {
                CinematicBudgetRepair::RaiseProfileEnvelopeWithNewConfiguration
            }
        };
        repairs.push(repair);
    }
    if tier != CinematicQualityTier::Final4k {
        repairs.push(CinematicBudgetRepair::LowerPreviewSppWithNewConfiguration);
        repairs.push(CinematicBudgetRepair::ReducePreviewAovsWithNewConfiguration);
        repairs.push(CinematicBudgetRepair::ShortenRangeWithNewConfiguration);
    }
    repairs.sort_unstable();
    repairs.dedup();
    repairs
}

const fn gib(value: u64) -> u64 {
    value * 1024 * 1024 * 1024
}

fn canonical_input(tier: CinematicQualityTier) -> CinematicQualityProfileInput {
    match tier {
        CinematicQualityTier::StoryboardSmoke => CinematicQualityProfileInput {
            tier,
            width_pixels: 960,
            height_pixels: 540,
            frames_per_second: 24,
            first_frame: 0,
            frame_count: 24,
            spp_floor: 1,
            spp_ceiling: 4,
            max_path_depth: 4,
            adaptive_error_ppm: 100_000,
            denoise_policy: DenoisePolicy::Disabled,
            aov_preset: AovPreset::BeautyXyz,
            shutter_samples: 1,
            tile_width: 32,
            tile_height: 32,
            worker_limit: 4,
            checkpoint_cadence_spp: 4,
            memory_ceiling_bytes: gib(1),
            per_frame_wall_time_ceiling_s: 30,
            sequence_wall_time_ceiling_s: 600,
            output_ceiling_bytes: gib(1),
            minimum_free_space_reserve_bytes: gib(1),
        },
        CinematicQualityTier::Daily1080p => CinematicQualityProfileInput {
            tier,
            width_pixels: 1_920,
            height_pixels: 1_080,
            frames_per_second: 24,
            first_frame: 0,
            frame_count: 240,
            spp_floor: 16,
            spp_ceiling: 64,
            max_path_depth: 8,
            adaptive_error_ppm: 20_000,
            denoise_policy: DenoisePolicy::PreviewBiased,
            aov_preset: AovPreset::DailyCore,
            shutter_samples: 4,
            tile_width: 32,
            tile_height: 32,
            worker_limit: 16,
            checkpoint_cadence_spp: 16,
            memory_ceiling_bytes: gib(2),
            per_frame_wall_time_ceiling_s: 900,
            sequence_wall_time_ceiling_s: 86_400,
            output_ceiling_bytes: gib(32),
            minimum_free_space_reserve_bytes: gib(20),
        },
        CinematicQualityTier::Qualification4kFrame => CinematicQualityProfileInput {
            tier,
            width_pixels: 3_840,
            height_pixels: 2_160,
            frames_per_second: 24,
            first_frame: 120,
            frame_count: 1,
            spp_floor: 64,
            spp_ceiling: 256,
            max_path_depth: 12,
            adaptive_error_ppm: 5_000,
            denoise_policy: DenoisePolicy::PreviewBiased,
            aov_preset: AovPreset::FinalDiagnostic,
            shutter_samples: 8,
            tile_width: 32,
            tile_height: 32,
            worker_limit: 32,
            checkpoint_cadence_spp: 32,
            memory_ceiling_bytes: gib(4),
            per_frame_wall_time_ceiling_s: 21_600,
            sequence_wall_time_ceiling_s: 21_600,
            output_ceiling_bytes: gib(4),
            minimum_free_space_reserve_bytes: gib(20),
        },
        CinematicQualityTier::Final4k => CinematicQualityProfileInput {
            tier,
            width_pixels: 3_840,
            height_pixels: 2_160,
            frames_per_second: 24,
            first_frame: 0,
            frame_count: 240,
            spp_floor: 256,
            spp_ceiling: 1_024,
            max_path_depth: 16,
            adaptive_error_ppm: 1_000,
            denoise_policy: DenoisePolicy::SeparateBiasedDerivative,
            aov_preset: AovPreset::FinalDiagnostic,
            shutter_samples: 16,
            tile_width: 32,
            tile_height: 32,
            worker_limit: 64,
            checkpoint_cadence_spp: 64,
            memory_ceiling_bytes: gib(8),
            per_frame_wall_time_ceiling_s: 86_400,
            sequence_wall_time_ceiling_s: 60 * 86_400,
            output_ceiling_bytes: gib(256),
            minimum_free_space_reserve_bytes: gib(50),
        },
    }
}

const fn tier_name(tier: CinematicQualityTier) -> &'static str {
    match tier {
        CinematicQualityTier::StoryboardSmoke => "storyboard-smoke",
        CinematicQualityTier::Daily1080p => "daily-1080p",
        CinematicQualityTier::Qualification4kFrame => "qualification-4k-frame",
        CinematicQualityTier::Final4k => "final-4k",
    }
}

/// Stable, actionable admission refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CinematicBudgetError {
    /// A dimension was zero or exceeded the structural maximum.
    InvalidDimensions,
    /// Frame rate/range was empty, overflowing, or outside the v1 envelope.
    InvalidFrameRange,
    /// SPP floor/ceiling was zero, reversed, or excessive.
    InvalidSppRange,
    /// Path depth was zero or excessive.
    InvalidPathDepth,
    /// Adaptive target was zero or exceeded one whole unit.
    InvalidAdaptiveTarget,
    /// No shutter sample was requested.
    InvalidShutterSamples,
    /// Tile edge was zero or excessive.
    InvalidTileShape,
    /// Worker limit was zero or excessive.
    InvalidWorkerLimit,
    /// Checkpoint cadence was zero or beyond the SPP ceiling.
    InvalidCheckpointCadence,
    /// A byte/time/reserve ceiling was zero or contradictory.
    InvalidResourceEnvelope,
    /// Fields violated the named tier's anti-downgrade contract.
    UnsupportedTierCombination(CinematicQualityTier),
    /// No host-specific path-throughput measurement was supplied.
    MissingThroughputMeasurement,
    /// A named checked estimate term exceeded `u64`.
    ArithmeticOverflow(&'static str),
    /// Estimate exceeded one or more explicit profile/host limits.
    InsufficientResources {
        /// Deterministically ordered exact failures.
        deficits: Vec<CinematicResourceDeficit>,
        /// Deterministically ranked explicit configuration/resource changes.
        repairs: Vec<CinematicBudgetRepair>,
        /// Full estimate retained for actionable diagnostics.
        estimate: Box<CinematicResourceEstimate>,
    },
}

impl CinematicBudgetError {
    /// Stable machine-readable diagnostic code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::InvalidDimensions => "cinematic-budget-invalid-dimensions",
            Self::InvalidFrameRange => "cinematic-budget-invalid-frame-range",
            Self::InvalidSppRange => "cinematic-budget-invalid-spp-range",
            Self::InvalidPathDepth => "cinematic-budget-invalid-path-depth",
            Self::InvalidAdaptiveTarget => "cinematic-budget-invalid-adaptive-target",
            Self::InvalidShutterSamples => "cinematic-budget-invalid-shutter-samples",
            Self::InvalidTileShape => "cinematic-budget-invalid-tile-shape",
            Self::InvalidWorkerLimit => "cinematic-budget-invalid-worker-limit",
            Self::InvalidCheckpointCadence => "cinematic-budget-invalid-checkpoint-cadence",
            Self::InvalidResourceEnvelope => "cinematic-budget-invalid-resource-envelope",
            Self::UnsupportedTierCombination(_) => "cinematic-budget-unsupported-tier-combination",
            Self::MissingThroughputMeasurement => "cinematic-budget-missing-throughput",
            Self::ArithmeticOverflow(_) => "cinematic-budget-arithmetic-overflow",
            Self::InsufficientResources { .. } => "cinematic-budget-insufficient-resources",
        }
    }
}

impl fmt::Display for CinematicBudgetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}: {self:?}", self.code())
    }
}

impl std::error::Error for CinematicBudgetError {}
