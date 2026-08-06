//! Deterministic, content-addressed image-sequence manifests.
//!
//! This module owns the L5 artifact grammar and an in-memory transactional
//! state machine. It deliberately owns no filesystem policy: callers choose
//! where relative paths live and how a completed manifest is persisted.

use core::fmt;
use core::fmt::Write as _;

pub use fs_blake3::ContentHash;
use fs_blake3::{DomainHasher, hash_bytes};

/// Version of the exact binary frame-sequence manifest grammar.
pub const FRAME_SEQUENCE_MANIFEST_VERSION: u16 = 1;

const MAGIC: &[u8; 8] = b"FSIMSEQ1";
const MANIFEST_IDENTITY_DOMAIN: &str = "org.frankensim.fs-img.frame-sequence-manifest.v1";
const IDENTITY_POLL_BYTES: usize = 64 * 1024;
const ZERO_HASH: ContentHash = ContentHash([0; 32]);

/// Hard limits for construction, decoding, and storage admission.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrameSequenceLimits {
    max_artifacts: u32,
    max_channels_per_artifact: u16,
    max_relative_path_bytes: u16,
    max_manifest_bytes: u64,
    max_output_bytes: u64,
}

impl FrameSequenceLimits {
    /// Construct a nonzero set of independent resource ceilings.
    ///
    /// # Errors
    /// Returns [`FrameSequenceError::InvalidLimit`] for a zero ceiling.
    pub fn try_new(
        max_artifacts: u32,
        max_channels_per_artifact: u16,
        max_relative_path_bytes: u16,
        max_manifest_bytes: u64,
        max_output_bytes: u64,
    ) -> Result<Self, FrameSequenceError> {
        for (field, value) in [
            ("max_artifacts", u64::from(max_artifacts)),
            (
                "max_channels_per_artifact",
                u64::from(max_channels_per_artifact),
            ),
            (
                "max_relative_path_bytes",
                u64::from(max_relative_path_bytes),
            ),
            ("max_manifest_bytes", max_manifest_bytes),
            ("max_output_bytes", max_output_bytes),
        ] {
            if value == 0 {
                return Err(FrameSequenceError::InvalidLimit { field });
            }
        }
        Ok(Self {
            max_artifacts,
            max_channels_per_artifact,
            max_relative_path_bytes,
            max_manifest_bytes,
            max_output_bytes,
        })
    }

    /// Maximum number of expected artifacts.
    #[must_use]
    pub const fn max_artifacts(self) -> u32 {
        self.max_artifacts
    }

    /// Maximum channels declared by one artifact.
    #[must_use]
    pub const fn max_channels_per_artifact(self) -> u16 {
        self.max_channels_per_artifact
    }

    /// Maximum UTF-8 byte length of one canonical relative path.
    #[must_use]
    pub const fn max_relative_path_bytes(self) -> u16 {
        self.max_relative_path_bytes
    }

    /// Maximum canonical manifest/snapshot length.
    #[must_use]
    pub const fn max_manifest_bytes(self) -> u64 {
        self.max_manifest_bytes
    }

    /// Maximum reserved output bytes for the complete sequence.
    #[must_use]
    pub const fn max_output_bytes(self) -> u64 {
        self.max_output_bytes
    }
}

impl Default for FrameSequenceLimits {
    fn default() -> Self {
        Self {
            max_artifacts: 100_000,
            max_channels_per_artifact: 64,
            max_relative_path_bytes: 512,
            max_manifest_bytes: 256 * 1024 * 1024,
            max_output_bytes: 16 * 1024 * 1024 * 1024 * 1024,
        }
    }
}

/// Immutable source and configuration identities shared by a sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrameSequenceContext {
    shot_id: ContentHash,
    trajectory_id: ContentHash,
    render_config_id: ContentHash,
    scene_id: ContentHash,
    build_id: ContentHash,
    profile_id: ContentHash,
}

impl FrameSequenceContext {
    /// Bind a sequence to six non-placeholder identities.
    ///
    /// # Errors
    /// Returns [`FrameSequenceError::PlaceholderIdentity`] for an all-zero
    /// identity. This is a placeholder check, not authenticity validation.
    pub fn try_new(
        shot_id: ContentHash,
        trajectory_id: ContentHash,
        render_config_id: ContentHash,
        scene_id: ContentHash,
        build_id: ContentHash,
        profile_id: ContentHash,
    ) -> Result<Self, FrameSequenceError> {
        for (field, value) in [
            ("shot_id", shot_id),
            ("trajectory_id", trajectory_id),
            ("render_config_id", render_config_id),
            ("scene_id", scene_id),
            ("build_id", build_id),
            ("profile_id", profile_id),
        ] {
            if value == ZERO_HASH {
                return Err(FrameSequenceError::PlaceholderIdentity { field });
            }
        }
        Ok(Self {
            shot_id,
            trajectory_id,
            render_config_id,
            scene_id,
            build_id,
            profile_id,
        })
    }

    /// Shot identity used in every canonical relative name.
    #[must_use]
    pub const fn shot_id(self) -> ContentHash {
        self.shot_id
    }

    /// Source trajectory identity.
    #[must_use]
    pub const fn trajectory_id(self) -> ContentHash {
        self.trajectory_id
    }

    /// Render configuration identity.
    #[must_use]
    pub const fn render_config_id(self) -> ContentHash {
        self.render_config_id
    }

    /// Scene identity.
    #[must_use]
    pub const fn scene_id(self) -> ContentHash {
        self.scene_id
    }

    /// Build identity asserted by the producer.
    #[must_use]
    pub const fn build_id(self) -> ContentHash {
        self.build_id
    }

    /// Exact image-profile identity.
    #[must_use]
    pub const fn profile_id(self) -> ContentHash {
        self.profile_id
    }
}

/// Semantic role of one frame artifact.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FrameArtifactRole {
    /// Raw scene-linear master, normally OpenEXR.
    RawMaster,
    /// Bias-labeled denoised scene-linear intermediate.
    DenoisedIntermediate,
    /// Display-referred PNG preview.
    DisplayPreview,
    /// Optional scientific overlay kept separate from beauty masters.
    ScientificOverlay,
}

impl FrameArtifactRole {
    const fn tag(self) -> u8 {
        match self {
            Self::RawMaster => 1,
            Self::DenoisedIntermediate => 2,
            Self::DisplayPreview => 3,
            Self::ScientificOverlay => 4,
        }
    }

    const fn directory(self) -> &'static str {
        match self {
            Self::RawMaster => "raw-masters",
            Self::DenoisedIntermediate => "denoised-intermediates",
            Self::DisplayPreview => "display-previews",
            Self::ScientificOverlay => "scientific-overlays",
        }
    }

    const fn rank(self) -> u8 {
        self.tag()
    }

    fn from_tag(tag: u8) -> Result<Self, FrameSequenceError> {
        match tag {
            1 => Ok(Self::RawMaster),
            2 => Ok(Self::DenoisedIntermediate),
            3 => Ok(Self::DisplayPreview),
            4 => Ok(Self::ScientificOverlay),
            _ => Err(FrameSequenceError::Malformed {
                field: "artifact role",
            }),
        }
    }

    /// Whether this role is a raw master rather than a derived artifact.
    #[must_use]
    pub const fn is_raw(self) -> bool {
        matches!(self, Self::RawMaster)
    }
}

/// On-disk image encoding of one artifact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameArtifactFormat {
    /// Single-part OpenEXR.
    OpenExr,
    /// Eight-bit PNG.
    Png8,
    /// Sixteen-bit PNG.
    Png16,
}

impl FrameArtifactFormat {
    const fn tag(self) -> u8 {
        match self {
            Self::OpenExr => 1,
            Self::Png8 => 2,
            Self::Png16 => 3,
        }
    }

    const fn extension(self) -> &'static str {
        match self {
            Self::OpenExr => "exr",
            Self::Png8 | Self::Png16 => "png",
        }
    }

    fn from_tag(tag: u8) -> Result<Self, FrameSequenceError> {
        match tag {
            1 => Ok(Self::OpenExr),
            2 => Ok(Self::Png8),
            3 => Ok(Self::Png16),
            _ => Err(FrameSequenceError::Malformed {
                field: "artifact format",
            }),
        }
    }
}

/// Scalar storage type of one named image channel.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameChannelType {
    /// IEEE binary16.
    Float16,
    /// IEEE binary32.
    Float32,
    /// Unsigned eight-bit sample.
    Uint8,
    /// Unsigned sixteen-bit sample.
    Uint16,
}

impl FrameChannelType {
    const fn tag(self) -> u8 {
        match self {
            Self::Float16 => 1,
            Self::Float32 => 2,
            Self::Uint8 => 3,
            Self::Uint16 => 4,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, FrameSequenceError> {
        match tag {
            1 => Ok(Self::Float16),
            2 => Ok(Self::Float32),
            3 => Ok(Self::Uint8),
            4 => Ok(Self::Uint16),
            _ => Err(FrameSequenceError::Malformed {
                field: "channel type",
            }),
        }
    }
}

/// One named, typed channel in exact stored order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameChannel {
    name: String,
    sample_type: FrameChannelType,
}

impl FrameChannel {
    /// Construct one nonempty NUL-free short channel name.
    ///
    /// # Errors
    /// Names outside the OpenEXR-v2-compatible 1..=31-byte subset refuse.
    pub fn try_new(
        name: impl Into<String>,
        sample_type: FrameChannelType,
    ) -> Result<Self, FrameSequenceError> {
        let name = name.into();
        if name.is_empty() || name.len() > 31 || name.as_bytes().contains(&0) {
            return Err(FrameSequenceError::InvalidChannelName);
        }
        Ok(Self { name, sample_type })
    }

    /// Channel name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Channel scalar storage type.
    #[must_use]
    pub const fn sample_type(&self) -> FrameChannelType {
        self.sample_type
    }
}

/// Exact per-frame sample-count summary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameSamplingStats {
    /// Every pixel received exactly `spp` samples.
    Uniform {
        /// Samples per pixel.
        spp: u32,
    },
    /// Per-pixel sample counts varied within an exact closed range.
    Adaptive {
        /// Minimum samples at any pixel.
        min_spp: u32,
        /// Maximum samples at any pixel.
        max_spp: u32,
        /// Exact sum of all per-pixel sample counts.
        total_samples: u64,
        /// Pixels that terminated before reaching `max_spp`.
        converged_pixels: u64,
        /// Pixels that reached `max_spp`.
        maximum_sample_pixels: u64,
    },
}

impl FrameSamplingStats {
    fn validate(self, width: u32, height: u32) -> Result<(), FrameSequenceError> {
        let pixels = u64::from(width)
            .checked_mul(u64::from(height))
            .ok_or(FrameSequenceError::SizeOverflow {
                context: "pixel count",
            })?;
        match self {
            Self::Uniform { spp } if spp != 0 => {
                pixels
                    .checked_mul(u64::from(spp))
                    .ok_or(FrameSequenceError::SizeOverflow {
                        context: "uniform sample count",
                    })?;
                Ok(())
            }
            Self::Adaptive {
                min_spp,
                max_spp,
                total_samples,
                converged_pixels,
                maximum_sample_pixels,
            } if min_spp != 0 && min_spp <= max_spp => {
                if converged_pixels
                    .checked_add(maximum_sample_pixels)
                    != Some(pixels)
                {
                    return Err(FrameSequenceError::InvalidSampling);
                }
                if min_spp == max_spp && converged_pixels != 0 {
                    return Err(FrameSequenceError::InvalidSampling);
                }
                let converged_minimum = converged_pixels
                    .checked_mul(u64::from(min_spp))
                    .ok_or(FrameSequenceError::SizeOverflow {
                        context: "adaptive converged minimum samples",
                    })?;
                let maximum_sample_total = maximum_sample_pixels
                    .checked_mul(u64::from(max_spp))
                    .ok_or(FrameSequenceError::SizeOverflow {
                        context: "adaptive maximum-sample pixels",
                    })?;
                let minimum = converged_minimum.checked_add(maximum_sample_total).ok_or(
                    FrameSequenceError::SizeOverflow {
                        context: "adaptive minimum samples",
                    },
                )?;
                let converged_maximum_spp = max_spp.saturating_sub(1);
                let converged_maximum = converged_pixels
                    .checked_mul(u64::from(converged_maximum_spp))
                    .ok_or(FrameSequenceError::SizeOverflow {
                        context: "adaptive converged maximum samples",
                    })?;
                let maximum = converged_maximum.checked_add(maximum_sample_total).ok_or(
                    FrameSequenceError::SizeOverflow {
                        context: "adaptive maximum samples",
                    },
                )?;
                if (minimum..=maximum).contains(&total_samples) {
                    Ok(())
                } else {
                    Err(FrameSequenceError::InvalidSampling)
                }
            }
            _ => Err(FrameSequenceError::InvalidSampling),
        }
    }
}

/// Stable key of one expected artifact within a sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FrameArtifactKey {
    frame_index: u64,
    segment_index: u32,
    role: FrameArtifactRole,
}

impl FrameArtifactKey {
    /// Construct a key from logical frame/segment indices and artifact role.
    #[must_use]
    pub const fn new(
        frame_index: u64,
        segment_index: u32,
        role: FrameArtifactRole,
    ) -> Self {
        Self {
            frame_index,
            segment_index,
            role,
        }
    }

    /// Logical frame index.
    #[must_use]
    pub const fn frame_index(self) -> u64 {
        self.frame_index
    }

    /// Event-delimited segment within the logical presentation frame.
    #[must_use]
    pub const fn segment_index(self) -> u32 {
        self.segment_index
    }

    /// Artifact role.
    #[must_use]
    pub const fn role(self) -> FrameArtifactRole {
        self.role
    }
}

/// Exact expected image metadata, independent of file content.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameArtifactDescriptor {
    key: FrameArtifactKey,
    frame_time_bits: u64,
    format: FrameArtifactFormat,
    width: u32,
    height: u32,
    channels: Vec<FrameChannel>,
    sampling: FrameSamplingStats,
}

impl FrameArtifactDescriptor {
    /// Construct and validate one expected frame artifact descriptor.
    ///
    /// Time is retained as exact finite binary64 bits with both signed zeros
    /// normalized to positive zero. EXR channels are sorted by name to match
    /// the writer; PNG channels are normalized to their standard packed order.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        frame_index: u64,
        segment_index: u32,
        role: FrameArtifactRole,
        frame_time_s: f64,
        format: FrameArtifactFormat,
        width: u32,
        height: u32,
        mut channels: Vec<FrameChannel>,
        sampling: FrameSamplingStats,
    ) -> Result<Self, FrameSequenceError> {
        if !frame_time_s.is_finite() {
            return Err(FrameSequenceError::InvalidFrameTime { frame_index });
        }
        if width == 0 || height == 0 {
            return Err(FrameSequenceError::InvalidDimensions { frame_index });
        }
        if channels.is_empty() {
            return Err(FrameSequenceError::InvalidChannelSet { frame_index });
        }
        channels.sort_unstable_by(|left, right| match format {
            FrameArtifactFormat::OpenExr => left.name.cmp(&right.name),
            FrameArtifactFormat::Png8 | FrameArtifactFormat::Png16 => {
                png_channel_rank(&left.name)
                    .cmp(&png_channel_rank(&right.name))
                    .then_with(|| left.name.cmp(&right.name))
            }
        });
        if channels
            .windows(2)
            .any(|pair| pair[0].name == pair[1].name)
        {
            return Err(FrameSequenceError::DuplicateChannel { frame_index });
        }
        validate_format_channels(role, format, &channels, frame_index)?;
        sampling.validate(width, height)?;
        let canonical_time = if frame_time_s == 0.0 {
            0.0
        } else {
            frame_time_s
        };
        Ok(Self {
            key: FrameArtifactKey::new(frame_index, segment_index, role),
            frame_time_bits: canonical_time.to_bits(),
            format,
            width,
            height,
            channels,
            sampling,
        })
    }

    /// Stable artifact key.
    #[must_use]
    pub const fn key(&self) -> FrameArtifactKey {
        self.key
    }

    /// Exact normalized frame time.
    #[must_use]
    pub fn frame_time_s(&self) -> f64 {
        f64::from_bits(self.frame_time_bits)
    }

    /// Exact normalized frame-time bit pattern.
    #[must_use]
    pub const fn frame_time_bits(&self) -> u64 {
        self.frame_time_bits
    }

    /// File format.
    #[must_use]
    pub const fn format(&self) -> FrameArtifactFormat {
        self.format
    }

    /// Raster width.
    #[must_use]
    pub const fn width(&self) -> u32 {
        self.width
    }

    /// Raster height.
    #[must_use]
    pub const fn height(&self) -> u32 {
        self.height
    }

    /// Exact channel schema in stored order.
    #[must_use]
    pub fn channels(&self) -> &[FrameChannel] {
        &self.channels
    }

    /// Sampling summary.
    #[must_use]
    pub const fn sampling(&self) -> FrameSamplingStats {
        self.sampling
    }
}

fn png_channel_rank(name: &str) -> u8 {
    match name {
        "Y" => 0,
        "R" => 1,
        "G" => 2,
        "B" => 3,
        "A" => 4,
        _ => u8::MAX,
    }
}

fn validate_format_channels(
    role: FrameArtifactRole,
    format: FrameArtifactFormat,
    channels: &[FrameChannel],
    frame_index: u64,
) -> Result<(), FrameSequenceError> {
    let types_match = match format {
        FrameArtifactFormat::OpenExr => channels.iter().all(|channel| {
            matches!(
                channel.sample_type,
                FrameChannelType::Float16 | FrameChannelType::Float32
            )
        }),
        FrameArtifactFormat::Png8 => png_channels_match(channels, FrameChannelType::Uint8),
        FrameArtifactFormat::Png16 => png_channels_match(channels, FrameChannelType::Uint16),
    };
    let role_matches = match role {
        FrameArtifactRole::RawMaster | FrameArtifactRole::DenoisedIntermediate => {
            format == FrameArtifactFormat::OpenExr
        }
        FrameArtifactRole::DisplayPreview => {
            matches!(
                format,
                FrameArtifactFormat::Png8 | FrameArtifactFormat::Png16
            )
        }
        FrameArtifactRole::ScientificOverlay => true,
    };
    if types_match && role_matches {
        Ok(())
    } else {
        Err(FrameSequenceError::InvalidChannelSet { frame_index })
    }
}

fn png_channels_match(channels: &[FrameChannel], sample_type: FrameChannelType) -> bool {
    let expected_names: &[&str] = match channels.len() {
        1 => &["Y"],
        3 => &["R", "G", "B"],
        4 => &["R", "G", "B", "A"],
        _ => return false,
    };
    channels
        .iter()
        .zip(expected_names)
        .all(|(channel, expected)| channel.name == *expected && channel.sample_type == sample_type)
}

/// Plan row supplied before rendering begins.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExpectedFrameArtifact {
    descriptor: FrameArtifactDescriptor,
    max_bytes: u64,
    source: Option<FrameArtifactKey>,
}

impl ExpectedFrameArtifact {
    /// Declare one output reservation and optional same-frame source link.
    ///
    /// Derived artifacts must point strictly backward in the closed role
    /// order; raw masters must not name a source.
    pub fn try_new(
        descriptor: FrameArtifactDescriptor,
        max_bytes: u64,
        source: Option<FrameArtifactKey>,
    ) -> Result<Self, FrameSequenceError> {
        if max_bytes == 0 {
            return Err(FrameSequenceError::InvalidArtifactLimit {
                key: descriptor.key,
            });
        }
        match (descriptor.key.role, source) {
            (FrameArtifactRole::RawMaster, None) => {}
            (FrameArtifactRole::RawMaster, Some(_)) | (_, None) => {
                return Err(FrameSequenceError::InvalidSource {
                    key: descriptor.key,
                });
            }
            (_, Some(source))
                if source.frame_index == descriptor.key.frame_index
                    && source.segment_index == descriptor.key.segment_index
                    && source.role.rank() < descriptor.key.role.rank() => {}
            _ => {
                return Err(FrameSequenceError::InvalidSource {
                    key: descriptor.key,
                });
            }
        }
        Ok(Self {
            descriptor,
            max_bytes,
            source,
        })
    }

    /// Expected descriptor.
    #[must_use]
    pub const fn descriptor(&self) -> &FrameArtifactDescriptor {
        &self.descriptor
    }

    /// Per-artifact byte reservation.
    #[must_use]
    pub const fn max_bytes(&self) -> u64 {
        self.max_bytes
    }

    /// Same-frame source artifact for a derived output.
    #[must_use]
    pub const fn source(&self) -> Option<FrameArtifactKey> {
        self.source
    }
}

/// Hash and exact byte size observed for a file, without decoding it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrameArtifactFileState {
    content_hash: ContentHash,
    byte_size: u64,
}

impl FrameArtifactFileState {
    /// Construct from a precomputed exact-byte observation.
    #[must_use]
    pub const fn new(content_hash: ContentHash, byte_size: u64) -> Self {
        Self {
            content_hash,
            byte_size,
        }
    }

    /// Hash an in-memory artifact with plain BLAKE3.
    ///
    /// # Errors
    /// Returns size overflow on platforms whose slice length exceeds `u64`.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, FrameSequenceError> {
        Ok(Self {
            content_hash: hash_bytes(bytes),
            byte_size: u64::try_from(bytes.len()).map_err(|_| {
                FrameSequenceError::SizeOverflow {
                    context: "artifact byte size",
                }
            })?,
        })
    }

    /// Plain BLAKE3 of exact file bytes.
    #[must_use]
    pub const fn content_hash(self) -> ContentHash {
        self.content_hash
    }

    /// Exact file length.
    #[must_use]
    pub const fn byte_size(self) -> u64 {
        self.byte_size
    }
}

/// Producer observation registered against one expected relative path.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameArtifactObservation {
    descriptor: FrameArtifactDescriptor,
    profile_id: ContentHash,
    file: FrameArtifactFileState,
    source_content_hash: Option<ContentHash>,
}

impl FrameArtifactObservation {
    /// Bind exact producer metadata to an already hashed file.
    #[must_use]
    pub const fn new(
        descriptor: FrameArtifactDescriptor,
        profile_id: ContentHash,
        file: FrameArtifactFileState,
        source_content_hash: Option<ContentHash>,
    ) -> Self {
        Self {
            descriptor,
            profile_id,
            file,
            source_content_hash,
        }
    }

    /// Hash exact bytes and construct an observation.
    pub fn from_bytes(
        descriptor: FrameArtifactDescriptor,
        profile_id: ContentHash,
        bytes: &[u8],
        source_content_hash: Option<ContentHash>,
    ) -> Result<Self, FrameSequenceError> {
        Ok(Self::new(
            descriptor,
            profile_id,
            FrameArtifactFileState::from_bytes(bytes)?,
            source_content_hash,
        ))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct FrameArtifactCompletion {
    file: FrameArtifactFileState,
    source_content_hash: Option<ContentHash>,
}

/// One canonical manifest row.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameArtifactEntry {
    relative_path: String,
    descriptor: FrameArtifactDescriptor,
    max_bytes: u64,
    source: Option<FrameArtifactKey>,
    completion: Option<FrameArtifactCompletion>,
}

impl FrameArtifactEntry {
    /// Canonical, root-independent relative path.
    #[must_use]
    pub fn relative_path(&self) -> &str {
        &self.relative_path
    }

    /// Expected descriptor.
    #[must_use]
    pub const fn descriptor(&self) -> &FrameArtifactDescriptor {
        &self.descriptor
    }

    /// Reserved maximum file bytes.
    #[must_use]
    pub const fn max_bytes(&self) -> u64 {
        self.max_bytes
    }

    /// Same-frame source key for a derived artifact.
    #[must_use]
    pub const fn source(&self) -> Option<FrameArtifactKey> {
        self.source
    }

    /// Registered file state, if complete.
    #[must_use]
    pub const fn file_state(&self) -> Option<FrameArtifactFileState> {
        match self.completion {
            Some(completion) => Some(completion.file),
            None => None,
        }
    }

    /// Registered source hash for a derived artifact.
    #[must_use]
    pub const fn source_content_hash(&self) -> Option<ContentHash> {
        match self.completion {
            Some(completion) => completion.source_content_hash,
            None => None,
        }
    }
}

/// Lifecycle state encoded in every snapshot.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameSequenceState {
    /// Some or all expected outputs are not yet sealed.
    Incomplete,
    /// Every output and source link passed independent file observation.
    Finalized,
}

impl FrameSequenceState {
    const fn tag(self) -> u8 {
        match self {
            Self::Incomplete => 0,
            Self::Finalized => 1,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, FrameSequenceError> {
        match tag {
            0 => Ok(Self::Incomplete),
            1 => Ok(Self::Finalized),
            _ => Err(FrameSequenceError::Malformed {
                field: "sequence state",
            }),
        }
    }
}

/// Result of idempotent artifact registration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RegistrationOutcome {
    /// The expected row transitioned from pending to complete.
    Recorded,
    /// An exact retry observed the already recorded row.
    AlreadyRecorded,
}

/// Immutable canonical bytes for resumable or finalized sequence state.
///
/// The identity is domain-separated from ordinary file-byte hashes. Callers
/// must persist or transmit it through an independent trusted channel before
/// using it as the pin supplied to [`FrameSequenceManifest::decode_snapshot`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameSequenceSnapshot {
    identity: ContentHash,
    bytes: Vec<u8>,
    state: FrameSequenceState,
}

impl FrameSequenceSnapshot {
    /// Domain-separated identity of the exact canonical snapshot bytes.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Exact canonical snapshot bytes.
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Lifecycle state encoded in the snapshot.
    #[must_use]
    pub const fn state(&self) -> FrameSequenceState {
        self.state
    }

    /// Consume the snapshot and return its canonical bytes.
    #[must_use]
    pub fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

/// Immutable finalized canonical bytes plus their typed identity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameSequenceSeal {
    identity: ContentHash,
    bytes: Vec<u8>,
    artifact_count: u32,
    output_bytes: u64,
}

impl FrameSequenceSeal {
    /// Domain-separated identity of the exact canonical bytes.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Canonical finalized manifest bytes.
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    /// Number of completed artifacts.
    #[must_use]
    pub const fn artifact_count(&self) -> u32 {
        self.artifact_count
    }

    /// Exact sum of completed artifact bytes.
    #[must_use]
    pub const fn output_bytes(&self) -> u64 {
        self.output_bytes
    }

    /// Consume the seal and return its canonical bytes.
    #[must_use]
    pub fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

/// Mutable, resumable sequence manifest.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrameSequenceManifest {
    context: FrameSequenceContext,
    limits: FrameSequenceLimits,
    entries: Vec<FrameArtifactEntry>,
    state: FrameSequenceState,
    completed_bytes: u64,
}

impl FrameSequenceManifest {
    /// Admit a complete expected sequence before rendering begins.
    ///
    /// `available_output_bytes` is a caller observation used only for this
    /// admission. It is intentionally not serialized, because free space is
    /// location- and time-dependent.
    ///
    /// # Errors
    /// Refuses duplicate/missing source rows, unsafe paths, arithmetic
    /// overflow, or a storage/manifest ceiling before publishing a builder.
    pub fn try_new(
        context: FrameSequenceContext,
        expected: Vec<ExpectedFrameArtifact>,
        limits: FrameSequenceLimits,
        available_output_bytes: u64,
    ) -> Result<Self, FrameSequenceError> {
        validate_limits(limits)?;
        if expected.is_empty() {
            return Err(FrameSequenceError::EmptySequence);
        }
        let count = u32::try_from(expected.len()).map_err(|_| {
            FrameSequenceError::ResourceLimit {
                resource: "artifact count",
                requested: u64::MAX,
                limit: u64::from(limits.max_artifacts),
            }
        })?;
        if count > limits.max_artifacts {
            return Err(FrameSequenceError::ResourceLimit {
                resource: "artifact count",
                requested: u64::from(count),
                limit: u64::from(limits.max_artifacts),
            });
        }
        let mut expected = expected;
        expected.sort_unstable_by_key(|artifact| artifact.descriptor.key);
        if expected.windows(2).any(|pair| {
            pair[0].descriptor.key == pair[1].descriptor.key
        }) {
            return Err(FrameSequenceError::DuplicateExpectedArtifact);
        }

        let mut entries = Vec::new();
        entries
            .try_reserve_exact(expected.len())
            .map_err(|_| FrameSequenceError::AllocationRefused {
                resource: "sequence entries",
                requested: u64::from(count),
            })?;
        let mut reserved_bytes = 0_u64;
        for artifact in expected {
            let channel_count = u16::try_from(artifact.descriptor.channels.len()).map_err(|_| {
                FrameSequenceError::ResourceLimit {
                    resource: "artifact channels",
                    requested: u64::MAX,
                    limit: u64::from(limits.max_channels_per_artifact),
                }
            })?;
            if channel_count > limits.max_channels_per_artifact {
                return Err(FrameSequenceError::ResourceLimit {
                    resource: "artifact channels",
                    requested: u64::from(channel_count),
                    limit: u64::from(limits.max_channels_per_artifact),
                });
            }
            let relative_path = canonical_relative_path(context, &artifact.descriptor)?;
            if relative_path.len() > usize::from(limits.max_relative_path_bytes) {
                return Err(FrameSequenceError::ResourceLimit {
                    resource: "relative path bytes",
                    requested: u64::try_from(relative_path.len()).unwrap_or(u64::MAX),
                    limit: u64::from(limits.max_relative_path_bytes),
                });
            }
            reserved_bytes = reserved_bytes.checked_add(artifact.max_bytes).ok_or(
                FrameSequenceError::SizeOverflow {
                    context: "reserved output bytes",
                },
            )?;
            entries.push(FrameArtifactEntry {
                relative_path,
                descriptor: artifact.descriptor,
                max_bytes: artifact.max_bytes,
                source: artifact.source,
                completion: None,
            });
        }
        validate_sources(&entries, false)?;
        admit_storage(reserved_bytes, limits.max_output_bytes, available_output_bytes)?;
        let manifest = Self {
            context,
            limits,
            entries,
            state: FrameSequenceState::Incomplete,
            completed_bytes: 0,
        };
        manifest.encoded_len()?;
        Ok(manifest)
    }

    /// Decode a strict canonical incomplete or finalized snapshot and
    /// re-admit its remaining reservations at the current location.
    ///
    /// # Errors
    /// Malformed, noncanonical, over-budget, stale-version, or structurally
    /// inconsistent bytes refuse without publishing partial state.
    pub fn decode_snapshot(
        bytes: &[u8],
        expected_identity: ContentHash,
        admission_limits: FrameSequenceLimits,
        available_output_bytes: u64,
    ) -> Result<Self, FrameSequenceError> {
        validate_limits(admission_limits)?;
        let byte_len = u64::try_from(bytes.len()).map_err(|_| {
            FrameSequenceError::SizeOverflow {
                context: "manifest input length",
            }
        })?;
        if byte_len > admission_limits.max_manifest_bytes {
            return Err(FrameSequenceError::ResourceLimit {
                resource: "manifest bytes",
                requested: byte_len,
                limit: admission_limits.max_manifest_bytes,
            });
        }
        let mut never_cancel = || true;
        let actual_identity = identity_with_poll(bytes, &mut never_cancel)?;
        if actual_identity != expected_identity {
            return Err(FrameSequenceError::IdentityMismatch {
                expected: expected_identity,
                actual: actual_identity,
            });
        }
        let mut reader = Reader::new(bytes);
        if reader.take(8)? != MAGIC {
            return Err(FrameSequenceError::Malformed { field: "magic" });
        }
        let version = reader.u16()?;
        if version != FRAME_SEQUENCE_MANIFEST_VERSION {
            return Err(FrameSequenceError::UnsupportedVersion { version });
        }
        let state = FrameSequenceState::from_tag(reader.u8()?)?;
        let context = FrameSequenceContext::try_new(
            reader.hash()?,
            reader.hash()?,
            reader.hash()?,
            reader.hash()?,
            reader.hash()?,
            reader.hash()?,
        )?;
        let limits = FrameSequenceLimits::try_new(
            reader.u32()?,
            reader.u16()?,
            reader.u16()?,
            reader.u64()?,
            reader.u64()?,
        )?;
        validate_nested_limits(limits, admission_limits)?;
        let entry_count = reader.u32()?;
        if entry_count > limits.max_artifacts {
            return Err(FrameSequenceError::ResourceLimit {
                resource: "artifact count",
                requested: u64::from(entry_count),
                limit: u64::from(limits.max_artifacts),
            });
        }
        let encoded_completed_bytes = reader.u64()?;
        let entry_count_usize = usize::try_from(entry_count).map_err(|_| {
            FrameSequenceError::SizeOverflow {
                context: "artifact count on this platform",
            }
        })?;
        let mut entries = Vec::new();
        entries.try_reserve_exact(entry_count_usize).map_err(|_| {
            FrameSequenceError::AllocationRefused {
                resource: "decoded sequence entries",
                requested: u64::from(entry_count),
            }
        })?;
        for _ in 0..entry_count_usize {
            entries.push(decode_entry(&mut reader, context, limits)?);
        }
        if !reader.is_empty() {
            return Err(FrameSequenceError::TrailingBytes);
        }
        let mut manifest = Self {
            context,
            limits,
            entries,
            state,
            completed_bytes: 0,
        };
        manifest.validate_structure()?;
        if manifest.completed_bytes != encoded_completed_bytes {
            return Err(FrameSequenceError::Malformed {
                field: "completed byte total",
            });
        }
        manifest.admit_remaining_storage(available_output_bytes)?;
        if manifest.encoded_len()? != byte_len {
            return Err(FrameSequenceError::NonCanonical);
        }
        let canonical = encode_manifest(&manifest, manifest.state, &mut never_cancel)?;
        if canonical.as_slice() != bytes {
            return Err(FrameSequenceError::NonCanonical);
        }
        Ok(manifest)
    }

    /// Immutable sequence context.
    #[must_use]
    pub const fn context(&self) -> FrameSequenceContext {
        self.context
    }

    /// Frozen resource limits encoded in snapshots.
    #[must_use]
    pub const fn limits(&self) -> FrameSequenceLimits {
        self.limits
    }

    /// Lifecycle state.
    #[must_use]
    pub const fn state(&self) -> FrameSequenceState {
        self.state
    }

    /// Canonically key-sorted expected rows.
    #[must_use]
    pub fn entries(&self) -> &[FrameArtifactEntry] {
        &self.entries
    }

    /// Number of completed rows.
    #[must_use]
    pub fn completed_artifacts(&self) -> u32 {
        u32::try_from(
            self.entries
                .iter()
                .filter(|entry| entry.completion.is_some())
                .count(),
        )
        .unwrap_or(u32::MAX)
    }

    /// Exact completed file bytes.
    #[must_use]
    pub const fn completed_bytes(&self) -> u64 {
        self.completed_bytes
    }

    /// Worst-case bytes still reserved for pending rows.
    #[must_use]
    pub fn remaining_reserved_bytes(&self) -> Result<u64, FrameSequenceError> {
        self.entries
            .iter()
            .filter(|entry| entry.completion.is_none())
            .try_fold(0_u64, |total, entry| {
                total.checked_add(entry.max_bytes).ok_or(
                    FrameSequenceError::SizeOverflow {
                        context: "remaining reserved bytes",
                    },
                )
            })
    }

    /// Re-admit the pending reservation after restart or path relocation.
    pub fn admit_remaining_storage(
        &self,
        available_output_bytes: u64,
    ) -> Result<(), FrameSequenceError> {
        let remaining = self.remaining_reserved_bytes()?;
        if remaining > available_output_bytes {
            return Err(FrameSequenceError::ResourceLimit {
                resource: "available output bytes",
                requested: remaining,
                limit: available_output_bytes,
            });
        }
        Ok(())
    }

    /// Register one exact producer observation transactionally.
    ///
    /// Exact retry is idempotent; conflicting duplicate content refuses.
    pub fn register_artifact(
        &mut self,
        relative_path: &str,
        observation: FrameArtifactObservation,
    ) -> Result<RegistrationOutcome, FrameSequenceError> {
        if self.state == FrameSequenceState::Finalized {
            return Err(FrameSequenceError::AlreadyFinalized);
        }
        let index = self
            .entries
            .iter()
            .position(|entry| entry.relative_path == relative_path)
            .ok_or_else(|| FrameSequenceError::UnexpectedArtifact {
                path: relative_path.to_owned(),
            })?;
        let entry = &self.entries[index];
        validate_observation(self.context, entry, &observation)?;
        if observation.file.byte_size == 0 {
            return Err(FrameSequenceError::EmptyArtifact {
                path: relative_path.to_owned(),
            });
        }
        if observation.file.byte_size > entry.max_bytes {
            return Err(FrameSequenceError::ResourceLimit {
                resource: "artifact bytes",
                requested: observation.file.byte_size,
                limit: entry.max_bytes,
            });
        }
        let completion = FrameArtifactCompletion {
            file: observation.file,
            source_content_hash: observation.source_content_hash,
        };
        if let Some(existing) = entry.completion {
            return if existing == completion {
                Ok(RegistrationOutcome::AlreadyRecorded)
            } else {
                Err(FrameSequenceError::ConflictingDuplicate {
                    path: relative_path.to_owned(),
                })
            };
        }
        let completed_bytes = self
            .completed_bytes
            .checked_add(observation.file.byte_size)
            .ok_or(FrameSequenceError::SizeOverflow {
                context: "completed artifact bytes",
            })?;
        if completed_bytes > self.limits.max_output_bytes {
            return Err(FrameSequenceError::ResourceLimit {
                resource: "completed output bytes",
                requested: completed_bytes,
                limit: self.limits.max_output_bytes,
            });
        }
        self.entries[index].completion = Some(completion);
        self.completed_bytes = completed_bytes;
        Ok(RegistrationOutcome::Recorded)
    }

    /// Hash bytes and register them against an expected row.
    pub fn register_artifact_bytes(
        &mut self,
        relative_path: &str,
        descriptor: FrameArtifactDescriptor,
        profile_id: ContentHash,
        bytes: &[u8],
        source_content_hash: Option<ContentHash>,
    ) -> Result<RegistrationOutcome, FrameSequenceError> {
        let observation = FrameArtifactObservation::from_bytes(
            descriptor,
            profile_id,
            bytes,
            source_content_hash,
        )?;
        self.register_artifact(relative_path, observation)
    }

    /// Encode a deterministic incomplete or finalized resumable snapshot.
    ///
    /// This bounded convenience path never observes cancellation. Use
    /// [`Self::snapshot_with_poll`] when a caller needs cancellation polling.
    pub fn snapshot(&self) -> Result<FrameSequenceSnapshot, FrameSequenceError> {
        self.snapshot_with_poll(|| true)
    }

    /// Encode and identify a deterministic snapshot while polling at artifact
    /// boundaries and fixed-size hash chunks.
    pub fn snapshot_with_poll(
        &self,
        mut poll: impl FnMut() -> bool,
    ) -> Result<FrameSequenceSnapshot, FrameSequenceError> {
        let bytes = encode_manifest(self, self.state, &mut poll)?;
        let identity = identity_with_poll(&bytes, &mut poll)?;
        Ok(FrameSequenceSnapshot {
            identity,
            bytes,
            state: self.state,
        })
    }

    /// Independently re-observe all files named by a finalized manifest.
    /// The callback resolves the root-independent relative names and need not
    /// decode any image.
    pub fn audit_with(
        &self,
        mut poll: impl FnMut() -> bool,
        mut observe: impl FnMut(&str) -> Option<FrameArtifactFileState>,
    ) -> Result<(), FrameSequenceError> {
        if self.state != FrameSequenceState::Finalized {
            return Err(FrameSequenceError::NotFinalized);
        }
        self.verify_complete(&mut poll, &mut observe)
    }

    /// Verify every expected file and atomically transition to finalized.
    /// Cancellation, missing/stale content, or allocation failure leaves the
    /// manifest unchanged and resumable.
    pub fn finalize_with(
        &mut self,
        mut poll: impl FnMut() -> bool,
        mut observe: impl FnMut(&str) -> Option<FrameArtifactFileState>,
    ) -> Result<FrameSequenceSeal, FrameSequenceError> {
        self.verify_complete(&mut poll, &mut observe)?;
        let bytes = encode_manifest(self, FrameSequenceState::Finalized, &mut poll)?;
        let identity = identity_with_poll(&bytes, &mut poll)?;
        self.state = FrameSequenceState::Finalized;
        Ok(FrameSequenceSeal {
            identity,
            bytes,
            artifact_count: u32::try_from(self.entries.len()).unwrap_or(u32::MAX),
            output_bytes: self.completed_bytes,
        })
    }

    fn verify_complete(
        &self,
        poll: &mut impl FnMut() -> bool,
        observe: &mut impl FnMut(&str) -> Option<FrameArtifactFileState>,
    ) -> Result<(), FrameSequenceError> {
        if !poll() {
            return Err(FrameSequenceError::Cancelled);
        }
        validate_sources(&self.entries, true)?;
        for entry in &self.entries {
            if !poll() {
                return Err(FrameSequenceError::Cancelled);
            }
            let completion = entry.completion.ok_or_else(|| {
                FrameSequenceError::MissingArtifact {
                    path: entry.relative_path.clone(),
                }
            })?;
            let actual = observe(&entry.relative_path).ok_or_else(|| {
                FrameSequenceError::MissingArtifact {
                    path: entry.relative_path.clone(),
                }
            })?;
            if actual != completion.file {
                return Err(FrameSequenceError::StaleArtifact {
                    path: entry.relative_path.clone(),
                    expected: completion.file,
                    actual,
                });
            }
        }
        Ok(())
    }

    fn validate_structure(&mut self) -> Result<(), FrameSequenceError> {
        if self.entries.is_empty() {
            return Err(FrameSequenceError::EmptySequence);
        }
        if self.entries.len() > usize::try_from(self.limits.max_artifacts).unwrap_or(usize::MAX) {
            return Err(FrameSequenceError::ResourceLimit {
                resource: "artifact count",
                requested: u64::try_from(self.entries.len()).unwrap_or(u64::MAX),
                limit: u64::from(self.limits.max_artifacts),
            });
        }
        if self
            .entries
            .windows(2)
            .any(|pair| pair[0].descriptor.key >= pair[1].descriptor.key)
        {
            return Err(FrameSequenceError::NonCanonical);
        }
        let mut reserved = 0_u64;
        let mut completed = 0_u64;
        for entry in &self.entries {
            if entry.descriptor.channels.len()
                > usize::from(self.limits.max_channels_per_artifact)
            {
                return Err(FrameSequenceError::ResourceLimit {
                    resource: "artifact channels",
                    requested: u64::try_from(entry.descriptor.channels.len()).unwrap_or(u64::MAX),
                    limit: u64::from(self.limits.max_channels_per_artifact),
                });
            }
            if entry.relative_path.len() > usize::from(self.limits.max_relative_path_bytes) {
                return Err(FrameSequenceError::ResourceLimit {
                    resource: "relative path bytes",
                    requested: u64::try_from(entry.relative_path.len()).unwrap_or(u64::MAX),
                    limit: u64::from(self.limits.max_relative_path_bytes),
                });
            }
            if entry.relative_path != canonical_relative_path(self.context, &entry.descriptor)? {
                return Err(FrameSequenceError::NonCanonical);
            }
            reserved = reserved.checked_add(entry.max_bytes).ok_or(
                FrameSequenceError::SizeOverflow {
                    context: "reserved output bytes",
                },
            )?;
            if let Some(completion) = entry.completion {
                if completion.file.byte_size == 0 || completion.file.byte_size > entry.max_bytes {
                    return Err(FrameSequenceError::Malformed {
                        field: "completed artifact size",
                    });
                }
                completed = completed.checked_add(completion.file.byte_size).ok_or(
                    FrameSequenceError::SizeOverflow {
                        context: "completed output bytes",
                    },
                )?;
            }
        }
        if reserved > self.limits.max_output_bytes {
            return Err(FrameSequenceError::ResourceLimit {
                resource: "reserved output bytes",
                requested: reserved,
                limit: self.limits.max_output_bytes,
            });
        }
        self.completed_bytes = completed;
        validate_sources(
            &self.entries,
            self.state == FrameSequenceState::Finalized,
        )?;
        if self.state == FrameSequenceState::Finalized
            && self.entries.iter().any(|entry| entry.completion.is_none())
        {
            return Err(FrameSequenceError::Malformed {
                field: "finalized completeness",
            });
        }
        Ok(())
    }

    fn encoded_len(&self) -> Result<u64, FrameSequenceError> {
        manifest_encoded_len(self)
    }
}

fn validate_limits(limits: FrameSequenceLimits) -> Result<(), FrameSequenceError> {
    FrameSequenceLimits::try_new(
        limits.max_artifacts,
        limits.max_channels_per_artifact,
        limits.max_relative_path_bytes,
        limits.max_manifest_bytes,
        limits.max_output_bytes,
    )
    .map(|_| ())
}

fn validate_nested_limits(
    nested: FrameSequenceLimits,
    admission: FrameSequenceLimits,
) -> Result<(), FrameSequenceError> {
    for (resource, requested, limit) in [
        (
            "artifact count",
            u64::from(nested.max_artifacts),
            u64::from(admission.max_artifacts),
        ),
        (
            "artifact channels",
            u64::from(nested.max_channels_per_artifact),
            u64::from(admission.max_channels_per_artifact),
        ),
        (
            "relative path bytes",
            u64::from(nested.max_relative_path_bytes),
            u64::from(admission.max_relative_path_bytes),
        ),
        (
            "manifest bytes",
            nested.max_manifest_bytes,
            admission.max_manifest_bytes,
        ),
        (
            "reserved output bytes",
            nested.max_output_bytes,
            admission.max_output_bytes,
        ),
    ] {
        if requested > limit {
            return Err(FrameSequenceError::ResourceLimit {
                resource,
                requested,
                limit,
            });
        }
    }
    Ok(())
}

fn admit_storage(
    reserved: u64,
    output_limit: u64,
    available: u64,
) -> Result<(), FrameSequenceError> {
    if reserved > output_limit {
        return Err(FrameSequenceError::ResourceLimit {
            resource: "reserved output bytes",
            requested: reserved,
            limit: output_limit,
        });
    }
    if reserved > available {
        return Err(FrameSequenceError::ResourceLimit {
            resource: "available output bytes",
            requested: reserved,
            limit: available,
        });
    }
    Ok(())
}

fn validate_sources(
    entries: &[FrameArtifactEntry],
    require_complete: bool,
) -> Result<(), FrameSequenceError> {
    for entry in entries {
        match (entry.source, entry.completion) {
            (None, Some(completion)) if completion.source_content_hash.is_some() => {
                return Err(FrameSequenceError::InvalidSource {
                    key: entry.descriptor.key,
                });
            }
            (Some(source), completion) => {
                let source_entry = entries
                    .binary_search_by_key(&source, |candidate| candidate.descriptor.key)
                    .ok()
                    .map(|index| &entries[index])
                    .ok_or(FrameSequenceError::InvalidSource {
                        key: entry.descriptor.key,
                    })?;
                if let Some(completion) = completion {
                    let declared = completion.source_content_hash.ok_or(
                        FrameSequenceError::InvalidSource {
                            key: entry.descriptor.key,
                        },
                    )?;
                    match source_entry.completion {
                        Some(source_completion) if source_completion.file.content_hash == declared => {}
                        Some(source_completion) => {
                            return Err(FrameSequenceError::SourceHashMismatch {
                                path: entry.relative_path.clone(),
                                expected: source_completion.file.content_hash,
                                actual: declared,
                            });
                        }
                        None if require_complete => {
                            return Err(FrameSequenceError::MissingArtifact {
                                path: source_entry.relative_path.clone(),
                            });
                        }
                        None => {}
                    }
                } else if require_complete {
                    return Err(FrameSequenceError::MissingArtifact {
                        path: entry.relative_path.clone(),
                    });
                }
            }
            _ => {}
        }
    }
    Ok(())
}

fn validate_observation(
    context: FrameSequenceContext,
    entry: &FrameArtifactEntry,
    observation: &FrameArtifactObservation,
) -> Result<(), FrameSequenceError> {
    let path = || entry.relative_path.clone();
    if observation.profile_id != context.profile_id {
        return Err(FrameSequenceError::DescriptorMismatch {
            path: path(),
            field: "profile identity",
        });
    }
    let expected = &entry.descriptor;
    let actual = &observation.descriptor;
    for (mismatch, field) in [
        (actual.key != expected.key, "artifact key"),
        (
            actual.frame_time_bits != expected.frame_time_bits,
            "frame time",
        ),
        (actual.format != expected.format, "format"),
        (
            actual.width != expected.width || actual.height != expected.height,
            "dimensions",
        ),
        (actual.channels != expected.channels, "channels"),
        (actual.sampling != expected.sampling, "sampling statistics"),
    ] {
        if mismatch {
            return Err(FrameSequenceError::DescriptorMismatch {
                path: path(),
                field,
            });
        }
    }
    match (entry.source, observation.source_content_hash) {
        (None, None) | (Some(_), Some(_)) => Ok(()),
        _ => Err(FrameSequenceError::DescriptorMismatch {
            path: path(),
            field: "source content hash",
        }),
    }
}

fn canonical_relative_path(
    context: FrameSequenceContext,
    descriptor: &FrameArtifactDescriptor,
) -> Result<String, FrameSequenceError> {
    let directory = descriptor.key.role.directory();
    let extension = descriptor.format.extension();
    let required = directory
        .len()
        .checked_add(1 + 64)
        .and_then(|value| value.checked_add(1 + 6 + 20))
        .and_then(|value| value.checked_add(6 + 16))
        .and_then(|value| value.checked_add(9 + 64))
        .and_then(|value| value.checked_add(1 + extension.len()))
        .ok_or(FrameSequenceError::SizeOverflow {
            context: "relative path length",
        })?;
    let mut path = String::new();
    path.try_reserve_exact(required)
        .map_err(|_| FrameSequenceError::AllocationRefused {
            resource: "relative path",
            requested: u64::try_from(required).unwrap_or(u64::MAX),
        })?;
    path.push_str(directory);
    path.push('/');
    push_hash_hex(&mut path, context.shot_id);
    write!(
        &mut path,
        "/frame-{:020}-time-{:016x}-profile-",
        descriptor.key.frame_index, descriptor.frame_time_bits
    )
    .map_err(|_| FrameSequenceError::AllocationRefused {
        resource: "relative path",
        requested: u64::try_from(required).unwrap_or(u64::MAX),
    })?;
    push_hash_hex(&mut path, context.profile_id);
    path.push('.');
    path.push_str(extension);
    debug_assert_eq!(path.len(), required);
    Ok(path)
}

fn push_hash_hex(output: &mut String, hash: ContentHash) {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    for byte in hash.as_bytes() {
        output.push(char::from(HEX[usize::from(byte >> 4)]));
        output.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
}

fn manifest_encoded_len(manifest: &FrameSequenceManifest) -> Result<u64, FrameSequenceError> {
    let mut bytes = 8_u64
        .checked_add(2 + 1)
        .and_then(|value| value.checked_add(6 * 32))
        .and_then(|value| value.checked_add(4 + 2 + 2 + 8 + 8))
        .and_then(|value| value.checked_add(4 + 8))
        .ok_or(FrameSequenceError::SizeOverflow {
            context: "manifest fixed header",
        })?;
    for entry in &manifest.entries {
        bytes = bytes
            .checked_add(2)
            .and_then(|value| {
                value.checked_add(u64::try_from(entry.relative_path.len()).unwrap_or(u64::MAX))
            })
            .and_then(|value| value.checked_add(8 + 1 + 8 + 1 + 4 + 4 + 2))
            .ok_or(FrameSequenceError::SizeOverflow {
                context: "manifest entry header",
            })?;
        for channel in &entry.descriptor.channels {
            bytes = bytes
                .checked_add(1)
                .and_then(|value| {
                    value.checked_add(u64::try_from(channel.name.len()).unwrap_or(u64::MAX))
                })
                .and_then(|value| value.checked_add(1))
                .ok_or(FrameSequenceError::SizeOverflow {
                    context: "manifest channel row",
                })?;
        }
        bytes = bytes
            .checked_add(match entry.descriptor.sampling {
                FrameSamplingStats::Uniform { .. } => 1 + 4,
                FrameSamplingStats::Adaptive { .. } => 1 + 4 + 4 + 8,
            })
            .and_then(|value| value.checked_add(8 + 1))
            .and_then(|value| {
                value.checked_add(if entry.source.is_some() { 8 + 1 } else { 0 })
            })
            .and_then(|value| value.checked_add(1))
            .and_then(|value| {
                value.checked_add(match entry.completion {
                    None => 0,
                    Some(completion) => {
                        32 + 8 + 1 + if completion.source_content_hash.is_some() {
                            32
                        } else {
                            0
                        }
                    }
                })
            })
            .ok_or(FrameSequenceError::SizeOverflow {
                context: "manifest entry tail",
            })?;
    }
    if bytes > manifest.limits.max_manifest_bytes {
        return Err(FrameSequenceError::ResourceLimit {
            resource: "manifest bytes",
            requested: bytes,
            limit: manifest.limits.max_manifest_bytes,
        });
    }
    Ok(bytes)
}

fn encode_manifest(
    manifest: &FrameSequenceManifest,
    state: FrameSequenceState,
    poll: &mut impl FnMut() -> bool,
) -> Result<Vec<u8>, FrameSequenceError> {
    if state == FrameSequenceState::Finalized {
        validate_sources(&manifest.entries, true)?;
        if manifest.entries.iter().any(|entry| entry.completion.is_none()) {
            return Err(FrameSequenceError::NotFinalized);
        }
    }
    let required = manifest_encoded_len(manifest)?;
    let capacity = usize::try_from(required).map_err(|_| FrameSequenceError::SizeOverflow {
        context: "manifest bytes on this platform",
    })?;
    let mut output = Vec::new();
    output.try_reserve_exact(capacity).map_err(|_| {
        FrameSequenceError::AllocationRefused {
            resource: "manifest bytes",
            requested: required,
        }
    })?;
    output.extend_from_slice(MAGIC);
    push_u16(&mut output, FRAME_SEQUENCE_MANIFEST_VERSION);
    output.push(state.tag());
    for hash in [
        manifest.context.shot_id,
        manifest.context.trajectory_id,
        manifest.context.render_config_id,
        manifest.context.scene_id,
        manifest.context.build_id,
        manifest.context.profile_id,
    ] {
        output.extend_from_slice(hash.as_bytes());
    }
    push_u32(&mut output, manifest.limits.max_artifacts);
    push_u16(&mut output, manifest.limits.max_channels_per_artifact);
    push_u16(&mut output, manifest.limits.max_relative_path_bytes);
    push_u64(&mut output, manifest.limits.max_manifest_bytes);
    push_u64(&mut output, manifest.limits.max_output_bytes);
    push_u32(
        &mut output,
        u32::try_from(manifest.entries.len()).unwrap_or(u32::MAX),
    );
    push_u64(&mut output, manifest.completed_bytes);
    for entry in &manifest.entries {
        if !poll() {
            return Err(FrameSequenceError::Cancelled);
        }
        push_short_string(&mut output, &entry.relative_path)?;
        push_u64(&mut output, entry.descriptor.key.frame_index);
        output.push(entry.descriptor.key.role.tag());
        push_u64(&mut output, entry.descriptor.frame_time_bits);
        output.push(entry.descriptor.format.tag());
        push_u32(&mut output, entry.descriptor.width);
        push_u32(&mut output, entry.descriptor.height);
        push_u16(
            &mut output,
            u16::try_from(entry.descriptor.channels.len()).unwrap_or(u16::MAX),
        );
        for channel in &entry.descriptor.channels {
            output.push(u8::try_from(channel.name.len()).unwrap_or(u8::MAX));
            output.extend_from_slice(channel.name.as_bytes());
            output.push(channel.sample_type.tag());
        }
        match entry.descriptor.sampling {
            FrameSamplingStats::Uniform { spp } => {
                output.push(1);
                push_u32(&mut output, spp);
            }
            FrameSamplingStats::Adaptive {
                min_spp,
                max_spp,
                total_samples,
            } => {
                output.push(2);
                push_u32(&mut output, min_spp);
                push_u32(&mut output, max_spp);
                push_u64(&mut output, total_samples);
            }
        }
        push_u64(&mut output, entry.max_bytes);
        match entry.source {
            None => output.push(0),
            Some(source) => {
                output.push(1);
                push_u64(&mut output, source.frame_index);
                output.push(source.role.tag());
            }
        }
        match entry.completion {
            None => output.push(0),
            Some(completion) => {
                output.push(1);
                output.extend_from_slice(completion.file.content_hash.as_bytes());
                push_u64(&mut output, completion.file.byte_size);
                match completion.source_content_hash {
                    None => output.push(0),
                    Some(source_hash) => {
                        output.push(1);
                        output.extend_from_slice(source_hash.as_bytes());
                    }
                }
            }
        }
    }
    debug_assert_eq!(output.len(), capacity);
    Ok(output)
}

fn identity_with_poll(
    bytes: &[u8],
    poll: &mut impl FnMut() -> bool,
) -> Result<ContentHash, FrameSequenceError> {
    let mut hasher = DomainHasher::new(MANIFEST_IDENTITY_DOMAIN);
    for chunk in bytes.chunks(IDENTITY_POLL_BYTES) {
        if !poll() {
            return Err(FrameSequenceError::Cancelled);
        }
        hasher.update(chunk);
    }
    Ok(hasher.finalize())
}

fn decode_entry(
    reader: &mut Reader<'_>,
    context: FrameSequenceContext,
    limits: FrameSequenceLimits,
) -> Result<FrameArtifactEntry, FrameSequenceError> {
    let relative_path = reader.short_string(usize::from(limits.max_relative_path_bytes))?;
    let frame_index = reader.u64()?;
    let role = FrameArtifactRole::from_tag(reader.u8()?)?;
    let frame_time_bits = reader.u64()?;
    let frame_time_s = f64::from_bits(frame_time_bits);
    if !frame_time_s.is_finite() || (frame_time_s == 0.0 && frame_time_bits != 0) {
        return Err(FrameSequenceError::InvalidFrameTime { frame_index });
    }
    let format = FrameArtifactFormat::from_tag(reader.u8()?)?;
    let width = reader.u32()?;
    let height = reader.u32()?;
    let channel_count = reader.u16()?;
    if channel_count == 0 || channel_count > limits.max_channels_per_artifact {
        return Err(FrameSequenceError::ResourceLimit {
            resource: "artifact channels",
            requested: u64::from(channel_count),
            limit: u64::from(limits.max_channels_per_artifact),
        });
    }
    let mut channels = Vec::new();
    channels
        .try_reserve_exact(usize::from(channel_count))
        .map_err(|_| FrameSequenceError::AllocationRefused {
            resource: "decoded channels",
            requested: u64::from(channel_count),
        })?;
    for _ in 0..channel_count {
        let name_len = usize::from(reader.u8()?);
        let name = reader.string_exact(name_len)?;
        channels.push(FrameChannel::try_new(
            name,
            FrameChannelType::from_tag(reader.u8()?)?,
        )?);
    }
    let sampling = match reader.u8()? {
        1 => FrameSamplingStats::Uniform { spp: reader.u32()? },
        2 => FrameSamplingStats::Adaptive {
            min_spp: reader.u32()?,
            max_spp: reader.u32()?,
            total_samples: reader.u64()?,
        },
        _ => {
            return Err(FrameSequenceError::Malformed {
                field: "sampling mode",
            });
        }
    };
    let descriptor = FrameArtifactDescriptor::try_new(
        frame_index,
        role,
        frame_time_s,
        format,
        width,
        height,
        channels,
        sampling,
    )?;
    let max_bytes = reader.u64()?;
    let source = match reader.u8()? {
        0 => None,
        1 => Some(FrameArtifactKey::new(
            reader.u64()?,
            FrameArtifactRole::from_tag(reader.u8()?)?,
        )),
        _ => {
            return Err(FrameSequenceError::Malformed {
                field: "source presence",
            });
        }
    };
    ExpectedFrameArtifact::try_new(descriptor.clone(), max_bytes, source)?;
    let completion = match reader.u8()? {
        0 => None,
        1 => {
            let file = FrameArtifactFileState::new(reader.hash()?, reader.u64()?);
            let source_content_hash = match reader.u8()? {
                0 => None,
                1 => Some(reader.hash()?),
                _ => {
                    return Err(FrameSequenceError::Malformed {
                        field: "source hash presence",
                    });
                }
            };
            Some(FrameArtifactCompletion {
                file,
                source_content_hash,
            })
        }
        _ => {
            return Err(FrameSequenceError::Malformed {
                field: "completion presence",
            });
        }
    };
    if relative_path != canonical_relative_path(context, &descriptor)? {
        return Err(FrameSequenceError::NonCanonical);
    }
    Ok(FrameArtifactEntry {
        relative_path,
        descriptor,
        max_bytes,
        source,
        completion,
    })
}

fn push_short_string(output: &mut Vec<u8>, value: &str) -> Result<(), FrameSequenceError> {
    let len = u16::try_from(value.len()).map_err(|_| FrameSequenceError::SizeOverflow {
        context: "short string length",
    })?;
    push_u16(output, len);
    output.extend_from_slice(value.as_bytes());
    Ok(())
}

fn push_u16(output: &mut Vec<u8>, value: u16) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn push_u32(output: &mut Vec<u8>, value: u32) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(output: &mut Vec<u8>, value: u64) {
    output.extend_from_slice(&value.to_le_bytes());
}

struct Reader<'a> {
    bytes: &'a [u8],
    position: usize,
}

impl<'a> Reader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, position: 0 }
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], FrameSequenceError> {
        let end = self
            .position
            .checked_add(len)
            .ok_or(FrameSequenceError::Truncated)?;
        let value = self
            .bytes
            .get(self.position..end)
            .ok_or(FrameSequenceError::Truncated)?;
        self.position = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, FrameSequenceError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, FrameSequenceError> {
        let bytes: [u8; 2] = self
            .take(2)?
            .try_into()
            .map_err(|_| FrameSequenceError::Truncated)?;
        Ok(u16::from_le_bytes(bytes))
    }

    fn u32(&mut self) -> Result<u32, FrameSequenceError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| FrameSequenceError::Truncated)?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, FrameSequenceError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| FrameSequenceError::Truncated)?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn hash(&mut self) -> Result<ContentHash, FrameSequenceError> {
        ContentHash::from_slice(self.take(32)?).ok_or(FrameSequenceError::Truncated)
    }

    fn short_string(&mut self, max_len: usize) -> Result<String, FrameSequenceError> {
        let len = usize::from(self.u16()?);
        if len > max_len {
            return Err(FrameSequenceError::ResourceLimit {
                resource: "relative path bytes",
                requested: u64::try_from(len).unwrap_or(u64::MAX),
                limit: u64::try_from(max_len).unwrap_or(u64::MAX),
            });
        }
        self.string_exact(len)
    }

    fn string_exact(&mut self, len: usize) -> Result<String, FrameSequenceError> {
        let text = core::str::from_utf8(self.take(len)?).map_err(|_| {
            FrameSequenceError::Malformed {
                field: "UTF-8 string",
            }
        })?;
        let mut output = String::new();
        output
            .try_reserve_exact(len)
            .map_err(|_| FrameSequenceError::AllocationRefused {
                resource: "decoded string",
                requested: u64::try_from(len).unwrap_or(u64::MAX),
            })?;
        output.push_str(text);
        Ok(output)
    }

    fn is_empty(&self) -> bool {
        self.position == self.bytes.len()
    }
}

/// Structured frame-sequence refusal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FrameSequenceError {
    /// A resource ceiling was zero.
    InvalidLimit {
        /// Stable field name.
        field: &'static str,
    },
    /// One context identity was the all-zero placeholder.
    PlaceholderIdentity {
        /// Stable identity field.
        field: &'static str,
    },
    /// Frame time was NaN or infinite.
    InvalidFrameTime {
        /// Logical frame index.
        frame_index: u64,
    },
    /// Raster dimensions included zero.
    InvalidDimensions {
        /// Logical frame index.
        frame_index: u64,
    },
    /// Channel name was empty, contained NUL, or exceeded 31 bytes.
    InvalidChannelName,
    /// Channel names were duplicated.
    DuplicateChannel {
        /// Logical frame index.
        frame_index: u64,
    },
    /// Channel types/count or format disagreed with the artifact role.
    InvalidChannelSet {
        /// Logical frame index.
        frame_index: u64,
    },
    /// Sample statistics were empty, reversed, or inconsistent with pixels.
    InvalidSampling,
    /// A per-artifact reservation was zero.
    InvalidArtifactLimit {
        /// Affected row.
        key: FrameArtifactKey,
    },
    /// A raw/derived source relationship was missing, cyclic, or mismatched.
    InvalidSource {
        /// Affected row.
        key: FrameArtifactKey,
    },
    /// Expected rows repeated one key.
    DuplicateExpectedArtifact,
    /// An exact caller ceiling was exceeded.
    ResourceLimit {
        /// Stable resource name.
        resource: &'static str,
        /// Exact required value.
        requested: u64,
        /// Supplied limit.
        limit: u64,
    },
    /// Checked size arithmetic overflowed.
    SizeOverflow {
        /// Stable quantity name.
        context: &'static str,
    },
    /// A fallible allocation was refused after admission.
    AllocationRefused {
        /// Stable allocation name.
        resource: &'static str,
        /// Logical requested units or bytes.
        requested: u64,
    },
    /// Registration named no expected relative path.
    UnexpectedArtifact {
        /// Supplied path.
        path: String,
    },
    /// Producer metadata differed from the expected row.
    DescriptorMismatch {
        /// Canonical expected path.
        path: String,
        /// First mismatching field.
        field: &'static str,
    },
    /// A completed file had zero bytes.
    EmptyArtifact {
        /// Canonical path.
        path: String,
    },
    /// A retry disagreed with an already completed row.
    ConflictingDuplicate {
        /// Canonical path.
        path: String,
    },
    /// A derived artifact named a stale source content hash.
    SourceHashMismatch {
        /// Canonical derived path.
        path: String,
        /// Actual registered source identity.
        expected: ContentHash,
        /// Identity claimed by the derived artifact.
        actual: ContentHash,
    },
    /// One expected or independently observed file was absent.
    MissingArtifact {
        /// Canonical path.
        path: String,
    },
    /// Independent file observation disagreed with the manifest.
    StaleArtifact {
        /// Canonical path.
        path: String,
        /// Manifest state.
        expected: FrameArtifactFileState,
        /// Fresh observation.
        actual: FrameArtifactFileState,
    },
    /// Registration was attempted after finalization.
    AlreadyFinalized,
    /// An operation requiring a final manifest received an incomplete one.
    NotFinalized,
    /// Caller cancellation was observed at a bounded poll point.
    Cancelled,
    /// Snapshot ended before a complete field.
    Truncated,
    /// Snapshot contained bytes after the complete canonical record.
    TrailingBytes,
    /// Snapshot used an unsupported schema version.
    UnsupportedVersion {
        /// Observed version.
        version: u16,
    },
    /// Snapshot field had an invalid closed-union tag or encoding.
    Malformed {
        /// Stable field name.
        field: &'static str,
    },
    /// Snapshot was structurally valid-looking but not canonical.
    NonCanonical,
}

impl fmt::Display for FrameSequenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLimit { field } => write!(formatter, "sequence limit {field} must be nonzero"),
            Self::PlaceholderIdentity { field } => write!(formatter, "sequence identity {field} is an all-zero placeholder"),
            Self::InvalidFrameTime { frame_index } => write!(formatter, "frame {frame_index} has a non-finite time"),
            Self::InvalidDimensions { frame_index } => write!(formatter, "frame {frame_index} has zero dimensions"),
            Self::InvalidChannelName => formatter.write_str("channel name must be NUL-free and 1..=31 UTF-8 bytes"),
            Self::DuplicateChannel { frame_index } => write!(formatter, "frame {frame_index} repeats a channel name"),
            Self::InvalidChannelSet { frame_index } => write!(formatter, "frame {frame_index} has channels incompatible with its role or format"),
            Self::InvalidSampling => formatter.write_str("sample statistics are inconsistent with the raster"),
            Self::InvalidArtifactLimit { key } => write!(formatter, "artifact {key:?} has a zero byte reservation"),
            Self::InvalidSource { key } => write!(formatter, "artifact {key:?} has an invalid source relationship"),
            Self::DuplicateExpectedArtifact => formatter.write_str("sequence repeats an expected artifact key"),
            Self::ResourceLimit { resource, requested, limit } => write!(formatter, "{resource} requires {requested}, above limit {limit}"),
            Self::SizeOverflow { context } => write!(formatter, "size arithmetic overflowed for {context}"),
            Self::AllocationRefused { resource, requested } => write!(formatter, "allocator refused {requested} units for {resource}"),
            Self::UnexpectedArtifact { path } => write!(formatter, "artifact path {path:?} is not expected"),
            Self::DescriptorMismatch { path, field } => write!(formatter, "artifact {path:?} mismatches expected {field}"),
            Self::EmptyArtifact { path } => write!(formatter, "artifact {path:?} is empty"),
            Self::ConflictingDuplicate { path } => write!(formatter, "artifact {path:?} conflicts with its completed row"),
            Self::SourceHashMismatch { path, expected, actual } => write!(formatter, "artifact {path:?} names stale source {actual}, expected {expected}"),
            Self::MissingArtifact { path } => write!(formatter, "artifact {path:?} is incomplete or missing"),
            Self::StaleArtifact { path, expected, actual } => write!(formatter, "artifact {path:?} changed from {expected:?} to {actual:?}"),
            Self::AlreadyFinalized => formatter.write_str("sequence is already finalized"),
            Self::NotFinalized => formatter.write_str("sequence is not complete and finalized"),
            Self::Cancelled => formatter.write_str("sequence operation was cancelled"),
            Self::Truncated => formatter.write_str("sequence snapshot is truncated"),
            Self::TrailingBytes => formatter.write_str("sequence snapshot has trailing bytes"),
            Self::UnsupportedVersion { version } => write!(formatter, "unsupported frame-sequence manifest version {version}"),
            Self::Malformed { field } => write!(formatter, "malformed frame-sequence field {field}"),
            Self::NonCanonical => formatter.write_str("frame-sequence snapshot is not canonical"),
        }
    }
}

impl std::error::Error for FrameSequenceError {}
