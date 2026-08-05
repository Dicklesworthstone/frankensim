//! Deterministic, bounded fixed-SPP work shards for local multi-process
//! rendering.
//!
//! A shard owns a rectangular block of the logical `(tile, sample)` grid and
//! returns immutable raw XYZ partial sums.  The random stream remains keyed by
//! absolute pixel and sample identities, so workers neither share mutable
//! sampler state nor depend on arrival order.  Full-SPP tile-only shards retain
//! the legacy serial film's per-pixel addition order.  When samples are split,
//! merge order is frozen by the plan's ascending sample ranges; that result is
//! bit-stable for every execution of the same plan, but floating-point
//! non-associativity means it is not claimed bit-identical to the legacy
//! monolithic accumulator.

use super::{
    CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION, CameraPath,
    DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION, DirectStrategy, Film, FilmTimeMode,
    LIGHTING_TRACER_BIT_SEMANTICS_VERSION, MOTION_TRACER_BIT_SEMANTICS_VERSION, RenderTileLayout,
    Sampler, Scene, Settings, TRACER_BIT_SEMANTICS_VERSION, TracerError, checked_pixel_len,
    requested_time_mode, trace_pixel_sample, validate_scene,
};
use crate::camera::{AnimatedCamera, CutSide};
use crate::motion::{ShutterConvention, ShutterDistribution, ShutterInterval};
use crate::spectral::y_integral;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_rand::qmc::Sobol;
use std::collections::BTreeMap;

/// Canonical shard-result wire schema.
pub const RENDER_SHARD_SCHEMA_VERSION: u16 = 1;
/// Bit-affecting fixed-SPP shard accumulation and merge semantics.
pub const RENDER_SHARD_SEMANTICS_VERSION: u32 = 1;
/// Domain for one fully bound logical renderer shard.
pub const RENDER_SHARD_IDENTITY_DOMAIN: &str = "org.frankensim.fs-render.uniform-render-shard.v1";
/// Domain for the canonical semantic result body and its integrity trailer.
pub const RENDER_SHARD_ARTIFACT_DOMAIN: &str =
    "org.frankensim.fs-render.uniform-render-shard-result.v1";

const RENDER_SHARD_EXECUTION_ENVIRONMENT_DOMAIN: &str =
    "org.frankensim.fs-render.uniform-render-shard-environment.v1";
const RESULT_MAGIC: &[u8; 8] = b"FSRSHD01";
const TIME_MODE_BYTES: u64 = 48;
const SPEC_HEADER_BYTES: u64 = 320;
const PAYLOAD_COUNT_BYTES: u64 = 8;
const RESULT_SEAL_BYTES: u64 = 32;
const RESULT_FIXED_BYTES: u64 = SPEC_HEADER_BYTES + PAYLOAD_COUNT_BYTES + RESULT_SEAL_BYTES;

/// Per-worker admission caps frozen into a shard identity.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RenderShardLimits {
    max_paths: u64,
    max_result_bytes: u64,
}

impl RenderShardLimits {
    /// Admit positive traced-path and encoded-result ceilings.
    pub fn try_new(max_paths: u64, max_result_bytes: u64) -> Result<Self, RenderShardError> {
        if max_paths == 0 {
            return Err(RenderShardError::InvalidLimit("max_paths"));
        }
        if max_result_bytes == 0 {
            return Err(RenderShardError::InvalidLimit("max_result_bytes"));
        }
        Ok(Self {
            max_paths,
            max_result_bytes,
        })
    }

    /// Maximum absolute pixel/sample paths evaluated by this shard.
    #[must_use]
    pub const fn max_paths(self) -> u64 {
        self.max_paths
    }

    /// Maximum canonical result bytes retained or emitted by this shard.
    #[must_use]
    pub const fn max_result_bytes(self) -> u64 {
        self.max_result_bytes
    }

    /// Admit the canonical result size for a shard covering `pixel_count`
    /// pixels and return that exact encoded byte count.
    pub fn admit_result_pixels(self, pixel_count: u64) -> Result<u64, RenderShardError> {
        let observed = encoded_result_bytes(pixel_count)?;
        if observed > self.max_result_bytes {
            return Err(RenderShardError::ResultByteLimit {
                limit: self.max_result_bytes,
                observed,
            });
        }
        Ok(observed)
    }
}

/// Aggregate caps for validating and publishing one complete film.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RenderShardMergeLimits {
    max_input_bytes: u64,
    max_output_bytes: u64,
}

impl RenderShardMergeLimits {
    /// Admit positive aggregate encoded-input and raw-film ceilings.
    pub fn try_new(max_input_bytes: u64, max_output_bytes: u64) -> Result<Self, RenderShardError> {
        if max_input_bytes == 0 {
            return Err(RenderShardError::InvalidLimit("max_input_bytes"));
        }
        if max_output_bytes == 0 {
            return Err(RenderShardError::InvalidLimit("max_output_bytes"));
        }
        Ok(Self {
            max_input_bytes,
            max_output_bytes,
        })
    }

    /// Aggregate canonical bytes accepted from all submitted results.
    #[must_use]
    pub const fn max_input_bytes(self) -> u64 {
        self.max_input_bytes
    }

    /// Maximum raw XYZ film payload published after complete validation.
    #[must_use]
    pub const fn max_output_bytes(self) -> u64 {
        self.max_output_bytes
    }
}

/// Fully bound rectangular fixed-SPP work item.
#[derive(Clone, Copy, Debug)]
pub struct UniformRenderShardSpec {
    plan_identity: ContentHash,
    frame_identity: ContentHash,
    frame_ordinal: u64,
    settings: Settings,
    time_mode: FilmTimeMode,
    layout: RenderTileLayout,
    tile_start: u64,
    tile_end: u64,
    sample_start: u32,
    sample_end: u32,
    limits: RenderShardLimits,
    execution_environment_identity: ContentHash,
    shard_identity: ContentHash,
    path_count: u64,
    payload_pixel_count: u64,
    encoded_result_bytes: u64,
}

impl UniformRenderShardSpec {
    /// Validate and freeze one contiguous tile/sample rectangle.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        plan_identity: ContentHash,
        frame_identity: ContentHash,
        frame_ordinal: u64,
        settings: Settings,
        time_mode: FilmTimeMode,
        layout: RenderTileLayout,
        tile_start: u64,
        tile_end: u64,
        sample_start: u32,
        sample_end: u32,
        limits: RenderShardLimits,
    ) -> Result<Self, RenderShardError> {
        require_nonzero("plan_identity", plan_identity)?;
        require_nonzero("frame_identity", frame_identity)?;
        validate_settings(settings)?;
        if layout.image_width() != settings.width || layout.image_height() != settings.height {
            return Err(RenderShardError::SpecMismatch("layout dimensions"));
        }
        validate_time_mode(time_mode, settings.seed)?;
        if tile_start >= tile_end || tile_end > layout.tile_count() {
            return Err(RenderShardError::InvalidTileRange {
                start: tile_start,
                end: tile_end,
                tile_count: layout.tile_count(),
            });
        }
        if sample_start >= sample_end || sample_end > settings.spp {
            return Err(RenderShardError::InvalidSampleRange {
                start: sample_start,
                end: sample_end,
                spp: settings.spp,
            });
        }
        let payload_pixel_count = tile_pixel_count(layout, tile_start, tile_end)?;
        let path_count = payload_pixel_count
            .checked_mul(u64::from(sample_end - sample_start))
            .ok_or(RenderShardError::ArithmeticOverflow("path_count"))?;
        if path_count > limits.max_paths {
            return Err(RenderShardError::PathLimit {
                limit: limits.max_paths,
                observed: path_count,
            });
        }
        let encoded_result_bytes = limits.admit_result_pixels(payload_pixel_count)?;
        let execution_environment_identity = execution_environment_identity();
        let mut spec = Self {
            plan_identity,
            frame_identity,
            frame_ordinal,
            settings,
            time_mode,
            layout,
            tile_start,
            tile_end,
            sample_start,
            sample_end,
            limits,
            execution_environment_identity,
            shard_identity: ContentHash([0; 32]),
            path_count,
            payload_pixel_count,
            encoded_result_bytes,
        };
        spec.shard_identity = spec_identity(&spec);
        Ok(spec)
    }

    /// External plan identity that authorized this work.
    #[must_use]
    pub const fn plan_identity(self) -> ContentHash {
        self.plan_identity
    }
    /// Exact logical frame identity.
    #[must_use]
    pub const fn frame_identity(self) -> ContentHash {
        self.frame_identity
    }
    /// Stable ordinal of the frame within its sequence.
    #[must_use]
    pub const fn frame_ordinal(self) -> u64 {
        self.frame_ordinal
    }
    /// Complete fixed-SPP tracer settings.
    #[must_use]
    pub const fn settings(self) -> Settings {
        self.settings
    }
    /// Static, motion, or cinematic time semantics.
    #[must_use]
    pub const fn time_mode(self) -> FilmTimeMode {
        self.time_mode
    }
    /// Canonical row-major tile layout.
    #[must_use]
    pub const fn layout(self) -> RenderTileLayout {
        self.layout
    }
    /// Inclusive logical tile index.
    #[must_use]
    pub const fn tile_start(self) -> u64 {
        self.tile_start
    }
    /// Exclusive logical tile index.
    #[must_use]
    pub const fn tile_end(self) -> u64 {
        self.tile_end
    }
    /// Inclusive absolute sample index.
    #[must_use]
    pub const fn sample_start(self) -> u32 {
        self.sample_start
    }
    /// Exclusive absolute sample index.
    #[must_use]
    pub const fn sample_end(self) -> u32 {
        self.sample_end
    }
    /// Resource caps frozen into this work identity.
    #[must_use]
    pub const fn limits(self) -> RenderShardLimits {
        self.limits
    }
    /// Domain-separated identity of the fully bound work item.
    #[must_use]
    pub const fn shard_identity(self) -> ContentHash {
        self.shard_identity
    }
    /// Exact number of pixel/sample paths traced by this shard.
    #[must_use]
    pub const fn path_count(self) -> u64 {
        self.path_count
    }
    /// Number of pixel partial sums in the result payload.
    #[must_use]
    pub const fn payload_pixel_count(self) -> u64 {
        self.payload_pixel_count
    }
    /// Exact canonical result length in bytes.
    #[must_use]
    pub const fn encoded_result_bytes(self) -> u64 {
        self.encoded_result_bytes
    }
}

impl PartialEq for UniformRenderShardSpec {
    fn eq(&self, other: &Self) -> bool {
        self.plan_identity == other.plan_identity
            && self.frame_identity == other.frame_identity
            && self.frame_ordinal == other.frame_ordinal
            && self.settings == other.settings
            && time_mode_bits_eq(self.time_mode, other.time_mode)
            && self.layout == other.layout
            && self.tile_start == other.tile_start
            && self.tile_end == other.tile_end
            && self.sample_start == other.sample_start
            && self.sample_end == other.sample_end
            && self.limits == other.limits
            && self.execution_environment_identity == other.execution_environment_identity
            && self.shard_identity == other.shard_identity
            && self.path_count == other.path_count
            && self.payload_pixel_count == other.payload_pixel_count
            && self.encoded_result_bytes == other.encoded_result_bytes
    }
}

impl Eq for UniformRenderShardSpec {}

/// Immutable complete result for one [`UniformRenderShardSpec`].
#[derive(Clone, Debug, PartialEq)]
pub struct UniformRenderShardResult {
    spec: UniformRenderShardSpec,
    xyz: Vec<[f64; 3]>,
    result_identity: ContentHash,
}

impl UniformRenderShardResult {
    /// Fully bound specification that produced this result.
    #[must_use]
    pub const fn spec(&self) -> &UniformRenderShardSpec {
        &self.spec
    }
    /// Identity of the work item that produced this result.
    #[must_use]
    pub const fn shard_identity(&self) -> ContentHash {
        self.spec.shard_identity
    }
    /// Domain-separated digest of the canonical result body.
    #[must_use]
    pub const fn result_identity(&self) -> ContentHash {
        self.result_identity
    }
    /// Raw XYZ sums in canonical tile/pixel order.
    #[must_use]
    pub fn xyz_sums(&self) -> &[[f64; 3]] {
        &self.xyz
    }
    /// Exact canonical result length in bytes.
    #[must_use]
    pub const fn encoded_result_bytes(&self) -> u64 {
        self.spec.encoded_result_bytes
    }

    /// Emit canonical bytes under both the spec's frozen cap and this call's
    /// explicit cap. The result already exists immutably; cancellation only
    /// prevents the returned buffer from escaping.
    pub fn encode_canonical(
        &self,
        max_bytes: u64,
        cx: &Cx<'_>,
    ) -> Result<Vec<u8>, RenderShardError> {
        let limit = max_bytes.min(self.spec.limits.max_result_bytes);
        if self.spec.encoded_result_bytes > limit {
            return Err(RenderShardError::ResultByteLimit {
                limit,
                observed: self.spec.encoded_result_bytes,
            });
        }
        checkpoint(cx)?;
        let capacity = usize::try_from(self.spec.encoded_result_bytes)
            .map_err(|_| RenderShardError::ArithmeticOverflow("encoded_result_bytes"))?;
        let mut bytes = Vec::new();
        bytes
            .try_reserve_exact(capacity)
            .map_err(|_| RenderShardError::Capacity("encoded result"))?;
        encode_spec_header(&mut bytes, &self.spec);
        bytes.extend_from_slice(&self.spec.payload_pixel_count.to_le_bytes());
        for (index, xyz) in self.xyz.iter().enumerate() {
            if index.is_multiple_of(4096) {
                checkpoint(cx)?;
            }
            for value in xyz {
                bytes.extend_from_slice(&value.to_bits().to_le_bytes());
            }
        }
        let actual = hash_result_body(&bytes);
        if actual != self.result_identity {
            return Err(RenderShardError::Integrity);
        }
        bytes.extend_from_slice(actual.as_bytes());
        if bytes.len() != capacity {
            return Err(RenderShardError::ArithmeticOverflow(
                "encoded result length",
            ));
        }
        checkpoint(cx)?;
        Ok(bytes)
    }

    /// Decode only against an externally pinned expected spec and identities.
    /// The wire bytes cannot mint their own plan or work authority.
    #[allow(clippy::too_many_arguments)]
    pub fn decode_canonical(
        bytes: &[u8],
        max_bytes: u64,
        expected_spec: &UniformRenderShardSpec,
        expected_plan_pin: ContentHash,
        expected_shard_pin: ContentHash,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderShardError> {
        if expected_spec.plan_identity != expected_plan_pin {
            return Err(RenderShardError::PlanIdentityMismatch {
                expected: expected_plan_pin,
                actual: expected_spec.plan_identity,
            });
        }
        if expected_spec.shard_identity != expected_shard_pin {
            return Err(RenderShardError::ShardIdentityMismatch {
                expected: expected_shard_pin,
                actual: expected_spec.shard_identity,
            });
        }
        let observed = u64::try_from(bytes.len())
            .map_err(|_| RenderShardError::ArithmeticOverflow("input length"))?;
        let limit = max_bytes.min(expected_spec.limits.max_result_bytes);
        if observed > limit || expected_spec.encoded_result_bytes > limit {
            return Err(RenderShardError::ResultByteLimit {
                limit,
                observed: observed.max(expected_spec.encoded_result_bytes),
            });
        }
        if observed < expected_spec.encoded_result_bytes {
            return Err(RenderShardError::Truncated);
        }
        if observed > expected_spec.encoded_result_bytes {
            return Err(RenderShardError::TrailingBytes);
        }
        checkpoint(cx)?;
        let header_len = usize::try_from(SPEC_HEADER_BYTES).expect("fixed header fits usize");
        let mut expected_header = Vec::with_capacity(header_len);
        encode_spec_header(&mut expected_header, expected_spec);
        if bytes.get(..header_len) != Some(expected_header.as_slice()) {
            return Err(RenderShardError::NonCanonical("result header"));
        }
        let count_end = header_len + 8;
        let count_bytes: [u8; 8] = bytes
            .get(header_len..count_end)
            .and_then(|raw| raw.try_into().ok())
            .ok_or(RenderShardError::Truncated)?;
        let count = u64::from_le_bytes(count_bytes);
        if count != expected_spec.payload_pixel_count {
            return Err(RenderShardError::NonCanonical("payload count"));
        }
        let seal_start = bytes.len() - 32;
        let actual = hash_result_body(&bytes[..seal_start]);
        let stored = ContentHash(
            bytes[seal_start..]
                .try_into()
                .map_err(|_| RenderShardError::Truncated)?,
        );
        if actual != stored {
            return Err(RenderShardError::Integrity);
        }
        let count_usize = usize::try_from(count)
            .map_err(|_| RenderShardError::ArithmeticOverflow("payload count"))?;
        let payload = &bytes[count_end..seal_start];
        let payload_len = count_usize
            .checked_mul(24)
            .ok_or(RenderShardError::ArithmeticOverflow("payload byte length"))?;
        if payload.len() != payload_len {
            return Err(RenderShardError::NonCanonical("payload length"));
        }
        let mut xyz = Vec::new();
        xyz.try_reserve_exact(count_usize)
            .map_err(|_| RenderShardError::Capacity("decoded payload"))?;
        for (index, raw) in payload.chunks_exact(24).enumerate() {
            if index.is_multiple_of(4096) {
                checkpoint(cx)?;
            }
            let mut sample = [0.0; 3];
            for channel in 0..3 {
                let start = channel * 8;
                let bits = u64::from_le_bytes(
                    raw[start..start + 8]
                        .try_into()
                        .map_err(|_| RenderShardError::Truncated)?,
                );
                sample[channel] = f64::from_bits(bits);
            }
            if sample.iter().any(|value| !value.is_finite()) {
                return Err(RenderShardError::NonFinitePayload);
            }
            xyz.push(sample);
        }
        checkpoint(cx)?;
        Ok(Self {
            spec: *expected_spec,
            xyz,
            result_identity: actual,
        })
    }
}

/// Fail-closed refusal from shard admission, execution, codec, or merge.
#[derive(Debug)]
pub enum RenderShardError {
    /// A named resource ceiling was zero.
    InvalidLimit(&'static str),
    /// A mandatory authority identity was the all-zero sentinel.
    ZeroIdentity(&'static str),
    /// A named fixed-SPP setting was zero or otherwise invalid.
    InvalidSettings(&'static str),
    /// The film time mode was uninitialized or inconsistent with the seed.
    InvalidTimeMode,
    /// The logical tile interval was empty or outside the layout.
    InvalidTileRange {
        /// Inclusive requested start.
        start: u64,
        /// Exclusive requested end.
        end: u64,
        /// Number of tiles in the bound layout.
        tile_count: u64,
    },
    /// The absolute sample interval was empty or outside the fixed SPP.
    InvalidSampleRange {
        /// Inclusive requested start.
        start: u32,
        /// Exclusive requested end.
        end: u32,
        /// Samples per pixel in the bound settings.
        spp: u32,
    },
    /// A named semantic field did not match the expected shard contract.
    SpecMismatch(&'static str),
    /// Checked integer arithmetic could not represent a named quantity.
    ArithmeticOverflow(&'static str),
    /// One work item exceeded its traced-path cap.
    PathLimit {
        /// Admitted maximum paths.
        limit: u64,
        /// Required paths.
        observed: u64,
    },
    /// One canonical result exceeded its encoded-byte cap.
    ResultByteLimit {
        /// Admitted maximum bytes.
        limit: u64,
        /// Required or supplied bytes.
        observed: u64,
    },
    /// Unique complete results exceeded the aggregate input cap.
    AggregateInputLimit {
        /// Admitted maximum bytes.
        limit: u64,
        /// Required bytes.
        observed: u64,
    },
    /// The completed raw XYZ film exceeded its output cap.
    OutputByteLimit {
        /// Admitted maximum bytes.
        limit: u64,
        /// Required bytes.
        observed: u64,
    },
    /// A named fallible retained allocation was refused.
    Capacity(&'static str),
    /// The execution scope requested cancellation.
    Cancelled,
    /// The underlying path tracer refused the work.
    Tracer(TracerError),
    /// A worker payload or completed accumulation was non-finite.
    NonFinitePayload,
    /// A canonical body digest or in-memory result invariant failed.
    Integrity,
    /// Canonical bytes ended before the expected boundary.
    Truncated,
    /// Canonical bytes continued past the expected boundary.
    TrailingBytes,
    /// A named wire field did not use its unique canonical encoding.
    NonCanonical(&'static str),
    /// The externally trusted plan pin did not match the expected spec.
    PlanIdentityMismatch {
        /// Trusted external plan identity.
        expected: ContentHash,
        /// Plan identity carried by the expected spec.
        actual: ContentHash,
    },
    /// The externally trusted shard pin did not match the expected spec.
    ShardIdentityMismatch {
        /// Trusted external shard identity.
        expected: ContentHash,
        /// Shard identity carried by the expected spec.
        actual: ContentHash,
    },
    /// A result belongs to another plan.
    ForeignPlan(ContentHash),
    /// A result belongs to another frame.
    ForeignFrame(ContentHash),
    /// A result names no shard in the expected complete work set.
    UnexpectedShard(ContentHash),
    /// Two valid, different results claimed the same shard identity.
    ConflictingDuplicate(ContentHash),
    /// No valid result was supplied for an expected shard.
    MissingShard(ContentHash),
    /// The expected shard rectangles leave a logical cell uncovered.
    CoverageGap {
        /// First uncovered tile.
        tile: u64,
        /// First uncovered absolute sample within that tile block.
        sample: u32,
    },
    /// The expected shard rectangles cover a logical cell more than once.
    CoverageOverlap {
        /// First multiply covered tile.
        tile: u64,
        /// First multiply covered absolute sample within that tile block.
        sample: u32,
    },
}

impl core::fmt::Display for RenderShardError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl core::error::Error for RenderShardError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Tracer(error) => Some(error),
            _ => None,
        }
    }
}

impl From<TracerError> for RenderShardError {
    fn from(error: TracerError) -> Self {
        if matches!(error, TracerError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Tracer(error)
        }
    }
}

/// Render one static fixed-SPP shard transactionally.
pub fn render_static_shard(
    scene: &Scene,
    cx: &Cx<'_>,
    spec: &UniformRenderShardSpec,
) -> Result<UniformRenderShardResult, RenderShardError> {
    render_shard_impl(scene, cx, spec, None, CameraPath::Legacy)
}

/// Render one legacy-camera motion shard transactionally.
pub fn render_motion_shard(
    scene: &Scene,
    cx: &Cx<'_>,
    spec: &UniformRenderShardSpec,
    shutter: ShutterInterval,
) -> Result<UniformRenderShardResult, RenderShardError> {
    render_shard_impl(scene, cx, spec, Some(shutter), CameraPath::Legacy)
}

/// Render one cinematic-camera shard transactionally.
pub fn render_cinematic_shard(
    scene: &Scene,
    camera: &AnimatedCamera,
    cut_side: CutSide,
    cx: &Cx<'_>,
    spec: &UniformRenderShardSpec,
    shutter: ShutterInterval,
) -> Result<UniformRenderShardResult, RenderShardError> {
    let exposure = camera
        .admit_shutter(cx, shutter, cut_side)
        .map_err(TracerError::from)?;
    render_shard_impl(
        scene,
        cx,
        spec,
        Some(shutter),
        CameraPath::Cinematic { camera, exposure },
    )
}

fn render_shard_impl(
    scene: &Scene,
    cx: &Cx<'_>,
    spec: &UniformRenderShardSpec,
    shutter: Option<ShutterInterval>,
    camera_path: CameraPath<'_>,
) -> Result<UniformRenderShardResult, RenderShardError> {
    checkpoint(cx)?;
    if spec.execution_environment_identity != execution_environment_identity() {
        return Err(RenderShardError::SpecMismatch("execution environment"));
    }
    let requested = requested_time_mode(shutter, spec.settings.seed, camera_path)?;
    if !time_mode_bits_eq(requested, spec.time_mode) {
        return Err(RenderShardError::SpecMismatch("time mode"));
    }
    let _ = checked_pixel_len(spec.settings.width, spec.settings.height)?;
    let lighting = validate_scene(scene, shutter)?;
    let payload_len = usize::try_from(spec.payload_pixel_count)
        .map_err(|_| RenderShardError::ArithmeticOverflow("payload_pixel_count"))?;
    let mut xyz = Vec::new();
    xyz.try_reserve_exact(payload_len)
        .map_err(|_| RenderShardError::Capacity("shard XYZ payload"))?;
    let key = [
        (spec.settings.seed & 0xffff_ffff) as u32,
        (spec.settings.seed >> 32) as u32,
    ];
    let sobol = (spec.settings.sampler == Sampler::OwenSobol)
        .then(|| Sobol::scrambled(3, spec.settings.seed));
    let kn = 1.0 / y_integral();
    for tile in spec.tile_start..spec.tile_end {
        checkpoint(cx)?;
        let bounds = spec
            .layout
            .bounds(tile)
            .ok_or(RenderShardError::InvalidTileRange {
                start: spec.tile_start,
                end: spec.tile_end,
                tile_count: spec.layout.tile_count(),
            })?;
        for py in bounds.y..bounds.y + bounds.height {
            checkpoint(cx)?;
            for px in bounds.x..bounds.x + bounds.width {
                let pixel = py
                    .checked_mul(spec.settings.width)
                    .and_then(|row| row.checked_add(px))
                    .ok_or(RenderShardError::ArithmeticOverflow("pixel identity"))?;
                let mut sum = [0.0; 3];
                for sample in spec.sample_start..spec.sample_end {
                    checkpoint(cx)?;
                    let value = trace_pixel_sample(
                        scene,
                        &lighting,
                        cx,
                        &spec.settings,
                        kn,
                        sobol.as_ref(),
                        key,
                        pixel,
                        sample,
                        shutter,
                        camera_path,
                    )?;
                    for channel in 0..3 {
                        sum[channel] += value[channel];
                    }
                }
                if sum.iter().any(|value| !value.is_finite()) {
                    return Err(RenderShardError::NonFinitePayload);
                }
                xyz.push(sum);
            }
        }
    }
    if xyz.len() != payload_len {
        return Err(RenderShardError::ArithmeticOverflow(
            "rendered payload length",
        ));
    }
    checkpoint(cx)?;
    let result_identity = result_identity(spec, &xyz);
    Ok(UniformRenderShardResult {
        spec: *spec,
        xyz,
        result_identity,
    })
}

/// Validate a complete result set privately, then publish one raw film.
///
/// Exact duplicates are ignored. A different valid result for one shard ID is
/// a conflict. Submitted order is nonsemantic: arithmetic follows ascending
/// tile ranges and, within each tile range, ascending sample ranges.
pub fn merge_uniform_shards(
    expected_specs: &[UniformRenderShardSpec],
    results: &[UniformRenderShardResult],
    limits: RenderShardMergeLimits,
    cx: &Cx<'_>,
) -> Result<Film, RenderShardError> {
    checkpoint(cx)?;
    let first = expected_specs
        .first()
        .ok_or(RenderShardError::SpecMismatch("empty expected shard set"))?;
    let output_bytes = u64::from(first.settings.width)
        .checked_mul(u64::from(first.settings.height))
        .and_then(|pixels| pixels.checked_mul(24))
        .ok_or(RenderShardError::ArithmeticOverflow("output film bytes"))?;
    if output_bytes > limits.max_output_bytes {
        return Err(RenderShardError::OutputByteLimit {
            limit: limits.max_output_bytes,
            observed: output_bytes,
        });
    }
    let mut ordered_specs: Vec<&UniformRenderShardSpec> = Vec::new();
    ordered_specs
        .try_reserve_exact(expected_specs.len())
        .map_err(|_| RenderShardError::Capacity("ordered specs"))?;
    ordered_specs.extend(expected_specs.iter());
    ordered_specs.sort_by_key(|spec| {
        (
            spec.tile_start,
            spec.tile_end,
            spec.sample_start,
            spec.sample_end,
            spec.shard_identity,
        )
    });
    validate_expected_specs(&ordered_specs, first)?;

    let mut expected_by_id = BTreeMap::new();
    for spec in &ordered_specs {
        if expected_by_id.insert(spec.shard_identity, **spec).is_some() {
            return Err(RenderShardError::CoverageOverlap {
                tile: spec.tile_start,
                sample: spec.sample_start,
            });
        }
    }
    let mut aggregate_input = 0_u64;
    let mut result_by_id: BTreeMap<ContentHash, &UniformRenderShardResult> = BTreeMap::new();
    for result in results {
        checkpoint(cx)?;
        let identity = result.shard_identity();
        let Some(expected) = expected_by_id.get(&identity) else {
            if result.spec.plan_identity != first.plan_identity {
                return Err(RenderShardError::ForeignPlan(result.spec.plan_identity));
            }
            if result.spec.frame_identity != first.frame_identity {
                return Err(RenderShardError::ForeignFrame(result.spec.frame_identity));
            }
            return Err(RenderShardError::UnexpectedShard(identity));
        };
        if result.spec != *expected {
            return Err(RenderShardError::UnexpectedShard(identity));
        }
        if result.xyz.len()
            != usize::try_from(result.spec.payload_pixel_count)
                .map_err(|_| RenderShardError::ArithmeticOverflow("payload count"))?
            || result.xyz.iter().flatten().any(|value| !value.is_finite())
            || result_identity(&result.spec, &result.xyz) != result.result_identity
        {
            return Err(RenderShardError::Integrity);
        }
        if let Some(prior) = result_by_id.get(&identity) {
            if prior.result_identity != result.result_identity {
                return Err(RenderShardError::ConflictingDuplicate(identity));
            }
        } else {
            aggregate_input = aggregate_input
                .checked_add(result.encoded_result_bytes())
                .ok_or(RenderShardError::ArithmeticOverflow(
                    "aggregate result bytes",
                ))?;
            if aggregate_input > limits.max_input_bytes {
                return Err(RenderShardError::AggregateInputLimit {
                    limit: limits.max_input_bytes,
                    observed: aggregate_input,
                });
            }
            result_by_id.insert(identity, result);
        }
    }
    for spec in &ordered_specs {
        if !result_by_id.contains_key(&spec.shard_identity) {
            return Err(RenderShardError::MissingShard(spec.shard_identity));
        }
    }

    let mut film = Film::try_new(first.settings.width, first.settings.height)?;
    for spec in ordered_specs {
        checkpoint(cx)?;
        let result = result_by_id
            .get(&spec.shard_identity)
            .ok_or(RenderShardError::MissingShard(spec.shard_identity))?;
        merge_result_into_film(&mut film, spec, result, cx)?;
    }
    checkpoint(cx)?;
    film.spp_done = first.settings.spp;
    film.time_mode = first.time_mode;
    Ok(film)
}

fn validate_expected_specs(
    ordered: &[&UniformRenderShardSpec],
    first: &UniformRenderShardSpec,
) -> Result<(), RenderShardError> {
    let mut tile_cursor = 0_u64;
    let mut index = 0usize;
    while index < ordered.len() {
        let tile_start = ordered[index].tile_start;
        let tile_end = ordered[index].tile_end;
        if tile_start > tile_cursor {
            return Err(RenderShardError::CoverageGap {
                tile: tile_cursor,
                sample: 0,
            });
        }
        if tile_start < tile_cursor {
            return Err(RenderShardError::CoverageOverlap {
                tile: tile_start,
                sample: 0,
            });
        }
        let mut sample_cursor = 0_u32;
        while index < ordered.len()
            && ordered[index].tile_start == tile_start
            && ordered[index].tile_end == tile_end
        {
            let spec = ordered[index];
            validate_common_spec(spec, first)?;
            if spec.sample_start > sample_cursor {
                return Err(RenderShardError::CoverageGap {
                    tile: tile_start,
                    sample: sample_cursor,
                });
            }
            if spec.sample_start < sample_cursor {
                return Err(RenderShardError::CoverageOverlap {
                    tile: tile_start,
                    sample: spec.sample_start,
                });
            }
            sample_cursor = spec.sample_end;
            index += 1;
        }
        if sample_cursor != first.settings.spp {
            return Err(RenderShardError::CoverageGap {
                tile: tile_start,
                sample: sample_cursor,
            });
        }
        tile_cursor = tile_end;
    }
    if tile_cursor != first.layout.tile_count() {
        return Err(RenderShardError::CoverageGap {
            tile: tile_cursor,
            sample: 0,
        });
    }
    Ok(())
}

fn validate_common_spec(
    spec: &UniformRenderShardSpec,
    first: &UniformRenderShardSpec,
) -> Result<(), RenderShardError> {
    if spec.plan_identity != first.plan_identity {
        return Err(RenderShardError::ForeignPlan(spec.plan_identity));
    }
    if spec.frame_identity != first.frame_identity || spec.frame_ordinal != first.frame_ordinal {
        return Err(RenderShardError::ForeignFrame(spec.frame_identity));
    }
    if spec.settings != first.settings
        || !time_mode_bits_eq(spec.time_mode, first.time_mode)
        || spec.layout != first.layout
        || spec.execution_environment_identity != first.execution_environment_identity
    {
        return Err(RenderShardError::SpecMismatch("frame render semantics"));
    }
    Ok(())
}

fn merge_result_into_film(
    film: &mut Film,
    spec: &UniformRenderShardSpec,
    result: &UniformRenderShardResult,
    cx: &Cx<'_>,
) -> Result<(), RenderShardError> {
    let mut source = 0usize;
    for tile in spec.tile_start..spec.tile_end {
        checkpoint(cx)?;
        let bounds = spec
            .layout
            .bounds(tile)
            .ok_or(RenderShardError::InvalidTileRange {
                start: spec.tile_start,
                end: spec.tile_end,
                tile_count: spec.layout.tile_count(),
            })?;
        for py in bounds.y..bounds.y + bounds.height {
            let destination = py as usize * film.width as usize + bounds.x as usize;
            for column in 0..bounds.width as usize {
                let value = *result.xyz.get(source).ok_or(RenderShardError::Integrity)?;
                let slot = &mut film.xyz[destination + column];
                if spec.sample_start == 0 {
                    *slot = value;
                } else {
                    for channel in 0..3 {
                        slot[channel] += value[channel];
                        if !slot[channel].is_finite() {
                            return Err(RenderShardError::NonFinitePayload);
                        }
                    }
                }
                source += 1;
            }
        }
    }
    if source != result.xyz.len() {
        return Err(RenderShardError::Integrity);
    }
    Ok(())
}

fn require_nonzero(field: &'static str, identity: ContentHash) -> Result<(), RenderShardError> {
    if identity.as_bytes().iter().all(|byte| *byte == 0) {
        Err(RenderShardError::ZeroIdentity(field))
    } else {
        Ok(())
    }
}

fn validate_settings(settings: Settings) -> Result<(), RenderShardError> {
    if settings.width == 0 {
        return Err(RenderShardError::InvalidSettings("width"));
    }
    if settings.height == 0 {
        return Err(RenderShardError::InvalidSettings("height"));
    }
    if settings.spp == 0 {
        return Err(RenderShardError::InvalidSettings("spp"));
    }
    if settings.max_depth == 0 {
        return Err(RenderShardError::InvalidSettings("max_depth"));
    }
    Ok(())
}

fn validate_time_mode(mode: FilmTimeMode, seed: u64) -> Result<(), RenderShardError> {
    match mode {
        FilmTimeMode::Uninitialized => Err(RenderShardError::InvalidTimeMode),
        FilmTimeMode::Static => Ok(()),
        FilmTimeMode::Motion {
            stream_identity, ..
        }
        | FilmTimeMode::Cinematic {
            stream_identity, ..
        } if stream_identity != seed => Err(RenderShardError::InvalidTimeMode),
        FilmTimeMode::Cinematic { shot_id: 0, .. } => Err(RenderShardError::InvalidTimeMode),
        FilmTimeMode::Motion { .. } | FilmTimeMode::Cinematic { .. } => Ok(()),
    }
}

fn tile_pixel_count(
    layout: RenderTileLayout,
    tile_start: u64,
    tile_end: u64,
) -> Result<u64, RenderShardError> {
    layout
        .pixel_count_in_range(tile_start, tile_end)
        .ok_or(RenderShardError::InvalidTileRange {
            start: tile_start,
            end: tile_end,
            tile_count: layout.tile_count(),
        })
}

fn encoded_result_bytes(payload_pixel_count: u64) -> Result<u64, RenderShardError> {
    RESULT_FIXED_BYTES
        .checked_add(
            payload_pixel_count
                .checked_mul(24)
                .ok_or(RenderShardError::ArithmeticOverflow("payload bytes"))?,
        )
        .ok_or(RenderShardError::ArithmeticOverflow("result bytes"))
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), RenderShardError> {
    cx.checkpoint().map_err(|_| RenderShardError::Cancelled)
}

fn execution_environment_identity() -> ContentHash {
    let probe = fs_substrate::CapabilityProbe::topology_only();
    let mut hasher = DomainHasher::new(RENDER_SHARD_EXECUTION_ENVIRONMENT_DOMAIN);
    let isa = match probe.isa {
        fs_substrate::Isa::Aarch64Apple => 0,
        fs_substrate::Isa::Aarch64Other => 1,
        fs_substrate::Isa::X86_64 => 2,
        fs_substrate::Isa::Other => 3,
    };
    hasher.update(&[isa]);
    hasher.update(&(probe.features.len() as u64).to_le_bytes());
    for feature in probe.features {
        hasher.update(&(feature.len() as u64).to_le_bytes());
        hasher.update(feature.as_bytes());
    }
    hasher.finalize()
}

fn spec_identity(spec: &UniformRenderShardSpec) -> ContentHash {
    let mut hasher = DomainHasher::new(RENDER_SHARD_IDENTITY_DOMAIN);
    hash_spec_fields(&mut hasher, spec);
    hasher.finalize()
}

fn hash_spec_fields(hasher: &mut DomainHasher, spec: &UniformRenderShardSpec) {
    hasher.update(&RENDER_SHARD_SCHEMA_VERSION.to_le_bytes());
    hasher.update(&RENDER_SHARD_SEMANTICS_VERSION.to_le_bytes());
    for version in [
        TRACER_BIT_SEMANTICS_VERSION,
        MOTION_TRACER_BIT_SEMANTICS_VERSION,
        CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION,
        DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION,
        LIGHTING_TRACER_BIT_SEMANTICS_VERSION,
    ] {
        hasher.update(&version.to_le_bytes());
    }
    hasher.update(spec.plan_identity.as_bytes());
    hasher.update(spec.frame_identity.as_bytes());
    hasher.update(spec.execution_environment_identity.as_bytes());
    hasher.update(&spec.frame_ordinal.to_le_bytes());
    hash_settings(hasher, spec.settings);
    hash_layout(hasher, spec.layout);
    hasher.update(&spec.tile_start.to_le_bytes());
    hasher.update(&spec.tile_end.to_le_bytes());
    hasher.update(&spec.sample_start.to_le_bytes());
    hasher.update(&spec.sample_end.to_le_bytes());
    hasher.update(&spec.limits.max_paths.to_le_bytes());
    hasher.update(&spec.limits.max_result_bytes.to_le_bytes());
    hasher.update(&spec.path_count.to_le_bytes());
    hasher.update(&spec.payload_pixel_count.to_le_bytes());
    hasher.update(&spec.encoded_result_bytes.to_le_bytes());
    hash_time_mode(hasher, spec.time_mode);
}

fn hash_settings(hasher: &mut DomainHasher, settings: Settings) {
    hasher.update(&settings.width.to_le_bytes());
    hasher.update(&settings.height.to_le_bytes());
    hasher.update(&settings.spp.to_le_bytes());
    hasher.update(&settings.max_depth.to_le_bytes());
    hasher.update(&[sampler_tag(settings.sampler)]);
    hasher.update(&[strategy_tag(settings.strategy)]);
    hasher.update(&settings.seed.to_le_bytes());
}

fn hash_layout(hasher: &mut DomainHasher, layout: RenderTileLayout) {
    hasher.update(&layout.image_width().to_le_bytes());
    hasher.update(&layout.image_height().to_le_bytes());
    hasher.update(&layout.tile_width().to_le_bytes());
    hasher.update(&layout.tile_height().to_le_bytes());
    hasher.update(&layout.tiles_x().to_le_bytes());
    hasher.update(&layout.tiles_y().to_le_bytes());
    hasher.update(&layout.tile_count().to_le_bytes());
}

fn hash_time_mode(hasher: &mut DomainHasher, mode: FilmTimeMode) {
    let mut bytes = Vec::with_capacity(TIME_MODE_BYTES as usize);
    encode_time_mode(&mut bytes, mode);
    hasher.update(&bytes);
}

fn result_identity(spec: &UniformRenderShardSpec, xyz: &[[f64; 3]]) -> ContentHash {
    let mut bytes = Vec::with_capacity(SPEC_HEADER_BYTES as usize + 8);
    encode_spec_header(&mut bytes, spec);
    bytes.extend_from_slice(&spec.payload_pixel_count.to_le_bytes());
    let mut hasher = DomainHasher::new(RENDER_SHARD_ARTIFACT_DOMAIN);
    hasher.update(&bytes);
    for sample in xyz {
        for value in sample {
            hasher.update(&value.to_bits().to_le_bytes());
        }
    }
    hasher.finalize()
}

fn hash_result_body(body: &[u8]) -> ContentHash {
    let mut hasher = DomainHasher::new(RENDER_SHARD_ARTIFACT_DOMAIN);
    hasher.update(body);
    hasher.finalize()
}

fn encode_spec_header(bytes: &mut Vec<u8>, spec: &UniformRenderShardSpec) {
    let start = bytes.len();
    bytes.extend_from_slice(RESULT_MAGIC);
    bytes.extend_from_slice(&RENDER_SHARD_SCHEMA_VERSION.to_le_bytes());
    bytes.extend_from_slice(&0_u16.to_le_bytes());
    bytes.extend_from_slice(&RENDER_SHARD_SEMANTICS_VERSION.to_le_bytes());
    bytes.extend_from_slice(spec.plan_identity.as_bytes());
    bytes.extend_from_slice(spec.shard_identity.as_bytes());
    bytes.extend_from_slice(spec.frame_identity.as_bytes());
    bytes.extend_from_slice(spec.execution_environment_identity.as_bytes());
    bytes.extend_from_slice(&spec.frame_ordinal.to_le_bytes());
    bytes.extend_from_slice(&spec.settings.width.to_le_bytes());
    bytes.extend_from_slice(&spec.settings.height.to_le_bytes());
    bytes.extend_from_slice(&spec.settings.spp.to_le_bytes());
    bytes.extend_from_slice(&spec.settings.max_depth.to_le_bytes());
    bytes.push(sampler_tag(spec.settings.sampler));
    bytes.push(strategy_tag(spec.settings.strategy));
    bytes.extend_from_slice(&[0; 6]);
    bytes.extend_from_slice(&spec.settings.seed.to_le_bytes());
    bytes.extend_from_slice(&spec.layout.image_width().to_le_bytes());
    bytes.extend_from_slice(&spec.layout.image_height().to_le_bytes());
    bytes.extend_from_slice(&spec.layout.tile_width().to_le_bytes());
    bytes.extend_from_slice(&spec.layout.tile_height().to_le_bytes());
    bytes.extend_from_slice(&spec.layout.tiles_x().to_le_bytes());
    bytes.extend_from_slice(&spec.layout.tiles_y().to_le_bytes());
    bytes.extend_from_slice(&spec.layout.tile_count().to_le_bytes());
    bytes.extend_from_slice(&spec.tile_start.to_le_bytes());
    bytes.extend_from_slice(&spec.tile_end.to_le_bytes());
    bytes.extend_from_slice(&spec.sample_start.to_le_bytes());
    bytes.extend_from_slice(&spec.sample_end.to_le_bytes());
    bytes.extend_from_slice(&spec.limits.max_paths.to_le_bytes());
    bytes.extend_from_slice(&spec.limits.max_result_bytes.to_le_bytes());
    bytes.extend_from_slice(&spec.path_count.to_le_bytes());
    bytes.extend_from_slice(&spec.encoded_result_bytes.to_le_bytes());
    encode_time_mode(bytes, spec.time_mode);
    debug_assert_eq!(bytes.len() - start, SPEC_HEADER_BYTES as usize);
}

fn encode_time_mode(bytes: &mut Vec<u8>, mode: FilmTimeMode) {
    let start = bytes.len();
    let (tag, shutter, stream_identity, shot_id) = match mode {
        FilmTimeMode::Uninitialized => (0, None, 0, 0),
        FilmTimeMode::Static => (1, None, 0, 0),
        FilmTimeMode::Motion {
            shutter,
            stream_identity,
        } => (2, Some(shutter), stream_identity, 0),
        FilmTimeMode::Cinematic {
            shutter,
            stream_identity,
            shot_id,
        } => (3, Some(shutter), stream_identity, shot_id),
    };
    let (convention, distribution, strata, open_bits, close_bits) =
        shutter.map_or((0, 0, 0, 0, 0), |shutter| {
            let convention = match shutter.convention() {
                ShutterConvention::Centered => 0,
                ShutterConvention::FrontLoaded => 1,
                ShutterConvention::BackLoaded => 2,
            };
            let (distribution, strata) = match shutter.distribution() {
                ShutterDistribution::UniformCounterV1 => (0, 0),
                ShutterDistribution::StratifiedCounterV1 { strata } => (1, strata),
            };
            (
                convention,
                distribution,
                strata,
                shutter.open_s().to_bits(),
                shutter.close_s().to_bits(),
            )
        });
    bytes.push(tag);
    bytes.push(convention);
    bytes.push(distribution);
    bytes.push(0);
    bytes.extend_from_slice(&strata.to_le_bytes());
    bytes.extend_from_slice(&open_bits.to_le_bytes());
    bytes.extend_from_slice(&close_bits.to_le_bytes());
    bytes.extend_from_slice(&stream_identity.to_le_bytes());
    bytes.extend_from_slice(&shot_id.to_le_bytes());
    bytes.extend_from_slice(&0_u64.to_le_bytes());
    debug_assert_eq!(bytes.len() - start, TIME_MODE_BYTES as usize);
}

fn time_mode_bits_eq(left: FilmTimeMode, right: FilmTimeMode) -> bool {
    match (left, right) {
        (FilmTimeMode::Uninitialized, FilmTimeMode::Uninitialized)
        | (FilmTimeMode::Static, FilmTimeMode::Static) => true,
        (
            FilmTimeMode::Motion {
                shutter: left_shutter,
                stream_identity: left_stream,
            },
            FilmTimeMode::Motion {
                shutter: right_shutter,
                stream_identity: right_stream,
            },
        ) => shutter_bits_eq(left_shutter, right_shutter) && left_stream == right_stream,
        (
            FilmTimeMode::Cinematic {
                shutter: left_shutter,
                stream_identity: left_stream,
                shot_id: left_shot,
            },
            FilmTimeMode::Cinematic {
                shutter: right_shutter,
                stream_identity: right_stream,
                shot_id: right_shot,
            },
        ) => {
            shutter_bits_eq(left_shutter, right_shutter)
                && left_stream == right_stream
                && left_shot == right_shot
        }
        _ => false,
    }
}

fn shutter_bits_eq(left: ShutterInterval, right: ShutterInterval) -> bool {
    left.open_s().to_bits() == right.open_s().to_bits()
        && left.close_s().to_bits() == right.close_s().to_bits()
        && left.convention() == right.convention()
        && left.distribution() == right.distribution()
}

const fn sampler_tag(sampler: Sampler) -> u8 {
    match sampler {
        Sampler::Iid => 0,
        Sampler::OwenSobol => 1,
    }
}

const fn strategy_tag(strategy: DirectStrategy) -> u8 {
    match strategy {
        DirectStrategy::NeeOnly => 0,
        DirectStrategy::BsdfOnly => 1,
        DirectStrategy::Mis => 2,
    }
}
