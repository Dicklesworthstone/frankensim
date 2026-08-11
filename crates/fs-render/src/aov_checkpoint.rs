//! Canonical, bounded checkpoints for progressive uniform cinematic AOV films.
//!
//! This is deliberately a separate wire format from the legacy tracer
//! checkpoint. It snapshots the exact raw beauty and AOV accumulators owned by
//! [`CinematicAovFilm`] after a committed sample prefix, so restore followed by
//! [`crate::tracer::render_cinematic_range_with_aovs`] preserves sequential
//! floating-point order. The codec owns no filesystem policy: callers stream
//! the sealed bytes into a transactional artifact store.

#[allow(clippy::wildcard_imports)]
// private codec is intentionally coupled to its parent state
use super::*;
use crate::motion::{ShutterConvention, ShutterDistribution};
use crate::tracer::{FilmTimeMode, Sampler};
use fs_exec::Cx;

/// Canonical cinematic-AOV checkpoint schema.
pub const CINEMATIC_AOV_CHECKPOINT_SCHEMA_VERSION: u16 = 1;
/// Domain-separated digest over every checkpoint byte except its final seal.
pub const CINEMATIC_AOV_CHECKPOINT_CONTENT_DOMAIN: &str =
    "org.frankensim.fs-render.cinematic-aov-checkpoint.v1";

const MAGIC: &[u8; 8] = b"FSRAOVC1";
const SEAL_BYTES: u64 = 32;
const STREAM_CHUNK_BYTES: usize = 64 * 1024;
const RESTORE_POLL_PIXELS: usize = 1_024;
const FIXED_HEADER_BYTES: u64 = 276;
const BOUND_HEADER_BYTES: u64 = 97;
const BEAUTY_PIXEL_BYTES: u64 = 3 * 8;
const COMMON_PIXEL_BYTES: u64 = 11 * 8 + 5 * 4;
const FINAL_PIXEL_BYTES: u64 = 12 * 8 + 1 + 8 + 4 + 8 + 4 + 4;

#[derive(Clone, Copy, Debug, PartialEq)]
struct ExpectedRenderBinding {
    settings: Settings,
    shutter: ShutterInterval,
    shot_id: u64,
    cut_side: CutSide,
    committed_samples_per_pixel: u32,
}

/// Caller-owned identity and render binding expected from a checkpoint.
///
/// A self-seal detects byte corruption but cannot decide which valid job the
/// caller intended to resume. Requiring this independent expectation prevents
/// a well-formed checkpoint for another frame, scene assertion, shutter, or
/// sample stream from being accepted merely because its bytes are coherent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CinematicAovCheckpointExpectation {
    config: CinematicAovConfig,
    binding: Option<ExpectedRenderBinding>,
}

impl CinematicAovCheckpointExpectation {
    /// Expect a pristine, unbound AOV film carrying exactly `config`.
    #[must_use]
    pub const fn unbound(config: CinematicAovConfig) -> Self {
        Self {
            config,
            binding: None,
        }
    }

    /// Expect a committed uniform sample prefix for one exact render binding.
    pub fn bound(
        config: CinematicAovConfig,
        settings: Settings,
        shutter: ShutterInterval,
        shot_id: u64,
        cut_side: CutSide,
        committed_samples_per_pixel: u32,
    ) -> Result<Self, CinematicAovCheckpointError> {
        if settings.width == 0
            || settings.height == 0
            || settings.spp == 0
            || shot_id == 0
            || committed_samples_per_pixel == 0
            || committed_samples_per_pixel > settings.spp
        {
            return Err(CinematicAovCheckpointError::InvalidState {
                field: "expected render binding",
            });
        }
        validate_reference_times(config, shutter)?;
        Ok(Self {
            config,
            binding: Some(ExpectedRenderBinding {
                settings,
                shutter,
                shot_id,
                cut_side,
                committed_samples_per_pixel,
            }),
        })
    }

    /// Expected complete AOV configuration identity.
    #[must_use]
    pub const fn config_identity(self) -> ContentHash {
        self.config.identity
    }

    /// Whether a committed render binding is required.
    #[must_use]
    pub const fn expects_bound_prefix(self) -> bool {
        self.binding.is_some()
    }

    /// Exact committed uniform sample prefix required from the checkpoint.
    #[must_use]
    pub const fn committed_samples_per_pixel(self) -> u32 {
        match self.binding {
            Some(binding) => binding.committed_samples_per_pixel,
            None => 0,
        }
    }
}

/// Evidence returned only after a complete checkpoint seal is emitted or
/// verified.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CinematicAovCheckpointReceipt {
    content_hash: ContentHash,
    byte_len: u64,
    samples_per_pixel: u32,
}

impl CinematicAovCheckpointReceipt {
    /// Domain-separated identity of the complete pre-seal payload.
    #[must_use]
    pub const fn content_hash(self) -> ContentHash {
        self.content_hash
    }

    /// Exact encoded byte length, including the final 32-byte seal.
    #[must_use]
    pub const fn byte_len(self) -> u64 {
        self.byte_len
    }

    /// Uniform committed sample prefix retained by the checkpoint.
    #[must_use]
    pub const fn samples_per_pixel(self) -> u32 {
        self.samples_per_pixel
    }
}

/// Fail-closed checkpoint validation, resource, codec, or cancellation error.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CinematicAovCheckpointError {
    /// The source or reconstructed AOV film violated its typed contract.
    Aov(CinematicAovError),
    /// An internal film field was inconsistent with its committed prefix.
    InvalidState {
        /// Stable name of the inconsistent state family.
        field: &'static str,
    },
    /// The canonical artifact exceeded the caller's explicit byte budget.
    ByteLimit {
        /// Exact required or supplied bytes.
        required: u64,
        /// Caller-declared ceiling.
        limit: u64,
    },
    /// Checked byte arithmetic overflowed.
    LengthOverflow,
    /// A bounded allocation was refused.
    AllocationRefused,
    /// The operation observed cancellation before publishing a receipt.
    Cancelled,
    /// The input ended before a complete canonical record was available.
    Truncated,
    /// Magic, tags, values, or canonical ordering were invalid.
    InvalidEnvelope,
    /// The encoded schema is not implemented by this binary.
    UnsupportedSchema {
        /// Rejected schema number.
        found: u16,
    },
    /// A bit-affecting semantic version differed from this binary.
    SemanticsMismatch,
    /// The coherent checkpoint did not describe the independently expected
    /// configuration or render binding.
    ExpectedBindingMismatch {
        /// Stable name of the mismatched expectation family.
        field: &'static str,
    },
    /// The final domain-separated content seal did not verify.
    SealMismatch,
    /// Bytes remained after the exact profile-specific payload.
    TrailingBytes,
}

impl core::fmt::Display for CinematicAovCheckpointError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Aov(error) => write!(formatter, "cinematic AOV checkpoint refused: {error}"),
            Self::InvalidState { field } => {
                write!(
                    formatter,
                    "cinematic AOV checkpoint state {field} is invalid"
                )
            }
            Self::ByteLimit { required, limit } => write!(
                formatter,
                "cinematic AOV checkpoint needs {required} bytes above limit {limit}"
            ),
            Self::LengthOverflow => {
                formatter.write_str("cinematic AOV checkpoint length overflowed")
            }
            Self::AllocationRefused => {
                formatter.write_str("cinematic AOV checkpoint allocation was refused")
            }
            Self::Cancelled => formatter.write_str("cinematic AOV checkpoint was cancelled"),
            Self::Truncated => formatter.write_str("cinematic AOV checkpoint is truncated"),
            Self::InvalidEnvelope => {
                formatter.write_str("cinematic AOV checkpoint envelope is invalid")
            }
            Self::UnsupportedSchema { found } => {
                write!(
                    formatter,
                    "unsupported cinematic AOV checkpoint schema {found}"
                )
            }
            Self::SemanticsMismatch => formatter
                .write_str("cinematic AOV checkpoint semantic versions do not match this binary"),
            Self::ExpectedBindingMismatch { field } => write!(
                formatter,
                "cinematic AOV checkpoint does not match expected {field}"
            ),
            Self::SealMismatch => {
                formatter.write_str("cinematic AOV checkpoint content seal does not match")
            }
            Self::TrailingBytes => {
                formatter.write_str("cinematic AOV checkpoint has trailing bytes")
            }
        }
    }
}

impl core::error::Error for CinematicAovCheckpointError {}

impl From<CinematicAovError> for CinematicAovCheckpointError {
    fn from(error: CinematicAovError) -> Self {
        Self::Aov(error)
    }
}

/// Streaming checkpoint failure, separating renderer refusal from sink I/O.
#[derive(Debug)]
pub enum CinematicAovCheckpointWriteError<E> {
    /// Renderer-side validation, budget, allocation, or cancellation refusal.
    Checkpoint(CinematicAovCheckpointError),
    /// Caller-provided transactional sink refused a chunk.
    Sink(E),
}

impl<E: core::fmt::Display> core::fmt::Display for CinematicAovCheckpointWriteError<E> {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Checkpoint(error) => error.fmt(formatter),
            Self::Sink(error) => {
                write!(formatter, "cinematic AOV checkpoint sink refused: {error}")
            }
        }
    }
}

impl<E: core::fmt::Debug + core::fmt::Display> core::error::Error
    for CinematicAovCheckpointWriteError<E>
{
}

impl CinematicAovFilm {
    /// Exact canonical checkpoint length for the currently committed prefix.
    pub fn checkpoint_byte_len(&self) -> Result<u64, CinematicAovCheckpointError> {
        validate_checkpoint_structure(self)?;
        checkpoint_byte_len(self)
    }

    /// Stream a canonical self-sealed checkpoint in chunks no larger than
    /// 64 KiB. Callers must discard any already-emitted chunks if this method
    /// returns an error; the receipt is the publication boundary.
    pub fn write_checkpoint<E>(
        &self,
        max_bytes: u64,
        cx: &Cx<'_>,
        mut emit: impl FnMut(&[u8]) -> Result<(), E>,
    ) -> Result<CinematicAovCheckpointReceipt, CinematicAovCheckpointWriteError<E>> {
        cx.checkpoint().map_err(|_| {
            CinematicAovCheckpointWriteError::Checkpoint(CinematicAovCheckpointError::Cancelled)
        })?;
        validate_checkpoint_structure(self)
            .map_err(CinematicAovCheckpointWriteError::Checkpoint)?;
        let required =
            checkpoint_byte_len(self).map_err(CinematicAovCheckpointWriteError::Checkpoint)?;
        if required > max_bytes {
            return Err(CinematicAovCheckpointWriteError::Checkpoint(
                CinematicAovCheckpointError::ByteLimit {
                    required,
                    limit: max_bytes,
                },
            ));
        }
        validate_checkpoint_state(self, cx)
            .map_err(CinematicAovCheckpointWriteError::Checkpoint)?;
        let mut writer = StreamWriter::try_new(cx, &mut emit)?;
        encode_film(&mut writer, self)?;
        let (content_hash, byte_len) = writer.finish()?;
        if byte_len != required {
            return Err(CinematicAovCheckpointWriteError::Checkpoint(
                CinematicAovCheckpointError::InvalidState {
                    field: "encoded byte length",
                },
            ));
        }
        Ok(CinematicAovCheckpointReceipt {
            content_hash,
            byte_len,
            samples_per_pixel: self.beauty.spp_done,
        })
    }

    /// Restore a complete canonical checkpoint transactionally. No partially
    /// decoded film escapes on refusal. A later progressive append still
    /// revalidates the current scene/camera palette and continuity guard.
    #[allow(clippy::too_many_lines)] // one linear decoder mirrors the canonical wire schema
    pub fn restore_checkpoint(
        expected: CinematicAovCheckpointExpectation,
        bytes: &[u8],
        max_bytes: u64,
        cx: &Cx<'_>,
    ) -> Result<(Self, CinematicAovCheckpointReceipt), CinematicAovCheckpointError> {
        cx.checkpoint()
            .map_err(|_| CinematicAovCheckpointError::Cancelled)?;
        let byte_len =
            u64::try_from(bytes.len()).map_err(|_| CinematicAovCheckpointError::LengthOverflow)?;
        if byte_len > max_bytes {
            return Err(CinematicAovCheckpointError::ByteLimit {
                required: byte_len,
                limit: max_bytes,
            });
        }
        if bytes.len() < MAGIC.len() + SEAL_BYTES as usize {
            return Err(CinematicAovCheckpointError::Truncated);
        }
        let payload_len = bytes.len() - SEAL_BYTES as usize;
        let (payload, seal) = bytes.split_at(payload_len);
        let mut hasher = DomainHasher::new(CINEMATIC_AOV_CHECKPOINT_CONTENT_DOMAIN);
        for chunk in payload.chunks(STREAM_CHUNK_BYTES) {
            cx.checkpoint()
                .map_err(|_| CinematicAovCheckpointError::Cancelled)?;
            hasher.update(chunk);
        }
        let content_hash = hasher.finalize();
        if seal != content_hash.as_bytes() {
            return Err(CinematicAovCheckpointError::SealMismatch);
        }

        let mut reader = Reader::new(payload);
        if reader.take(MAGIC.len())? != MAGIC {
            return Err(CinematicAovCheckpointError::InvalidEnvelope);
        }
        let schema = reader.u16()?;
        if schema != CINEMATIC_AOV_CHECKPOINT_SCHEMA_VERSION {
            return Err(CinematicAovCheckpointError::UnsupportedSchema { found: schema });
        }
        let profile = decode_profile(reader.u8()?)?;
        if reader.u32()? != CINEMATIC_AOV_SEMANTICS_VERSION
            || reader.u32()? != CINEMATIC_AOV_CATEGORY_SEMANTICS_VERSION
            || reader.u32()? != CINEMATIC_AOV_ALBEDO_SEMANTICS_VERSION
            || reader.u32()? != TRACER_BIT_SEMANTICS_VERSION
            || reader.u32()? != MOTION_TRACER_BIT_SEMANTICS_VERSION
            || reader.u32()? != CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION
            || reader.u32()? != DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION
            || reader.u32()? != LIGHTING_TRACER_BIT_SEMANTICS_VERSION
            || reader.u32()? != MOTION_VECTOR_SEMANTICS_VERSION
            || reader.u32()? != crate::charts::CHART_BACKEND_BIT_SEMANTICS_VERSION
        {
            return Err(CinematicAovCheckpointError::SemanticsMismatch);
        }
        let width = reader.u32()?;
        let height = reader.u32()?;
        let spp_done = reader.u32()?;
        let provenance = CinematicAovProvenance::try_new(
            reader.u64()?,
            reader.canonical_f64()?,
            reader.canonical_f64()?,
            reader.canonical_f64()?,
            reader.hash()?,
            reader.hash()?,
            reader.hash()?,
        )?;
        let limits = CinematicAovLimits::try_new(
            reader.u64()?,
            reader.u64()?,
            reader.u64()?,
            reader.u64()?,
            reader.u64()?,
            reader.u64()?,
            reader.u32()?,
        )?;
        let config = CinematicAovConfig::new(profile, provenance, limits);
        if reader.hash()? != config.identity() {
            return Err(CinematicAovCheckpointError::SemanticsMismatch);
        }
        if config != expected.config {
            return Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
                field: "AOV configuration",
            });
        }

        let binding = match reader.u8()? {
            0 => {
                let expected = checkpoint_byte_len_fields(width, height, profile, None, 0, 0)?;
                require_exact_length(byte_len, expected)?;
                None
            }
            1 => {
                let settings = Settings {
                    width: reader.u32()?,
                    height: reader.u32()?,
                    spp: reader.u32()?,
                    max_depth: reader.u32()?,
                    sampler: decode_sampler(reader.u8()?)?,
                    strategy: decode_strategy(reader.u8()?)?,
                    seed: reader.u64()?,
                };
                let shutter = decode_shutter(&mut reader)?;
                let shot_id = reader.u64()?;
                let cut_side = match reader.u8()? {
                    0 => CutSide::Before,
                    1 => CutSide::After,
                    _ => return Err(CinematicAovCheckpointError::InvalidEnvelope),
                };
                let continuity_fingerprint = reader.hash()?;
                let object_count = reader.u32()?;
                let material_count = reader.u32()?;
                validate_palette_counts(object_count, material_count, limits)?;
                let expected = checkpoint_byte_len_fields(
                    width,
                    height,
                    profile,
                    Some(settings),
                    object_count,
                    material_count,
                )?;
                require_exact_length(byte_len, expected)?;
                let object_ids = reader.ascending_nonzero_u64s(object_count, cx)?;
                let material_identities = reader.ascending_nonzero_hashes(material_count, cx)?;
                Some(CinematicAovRenderBinding {
                    settings,
                    shutter,
                    shot_id,
                    cut_side,
                    palette: CinematicAovPalette {
                        object_ids,
                        material_identities,
                    },
                    continuity_fingerprint,
                    adaptive_policy: None,
                })
            }
            _ => return Err(CinematicAovCheckpointError::InvalidEnvelope),
        };
        require_expected_binding(binding.as_ref(), spp_done, expected.binding)?;

        let mut film = CinematicAovFilm::try_new(width, height, config)?;
        film.beauty.spp_done = spp_done;
        film.beauty.time_mode = binding
            .as_ref()
            .map_or(FilmTimeMode::Uninitialized, |binding| {
                FilmTimeMode::Cinematic {
                    shutter: binding.shutter,
                    stream_identity: binding.settings.seed,
                    shot_id: binding.shot_id,
                }
            });
        film.binding = binding;

        for pixel in 0..film.beauty.xyz.len() {
            if pixel.is_multiple_of(RESTORE_POLL_PIXELS) {
                cx.checkpoint()
                    .map_err(|_| CinematicAovCheckpointError::Cancelled)?;
            }
            film.beauty.xyz[pixel] = reader.f64_array()?;
            if let Some(common) = &mut film.common {
                common[pixel] = decode_common(&mut reader)?;
            }
            if let Some(final_diagnostic) = &mut film.final_diagnostic {
                final_diagnostic[pixel] = decode_final(&mut reader)?;
            }
        }
        if !reader.is_done() {
            return Err(CinematicAovCheckpointError::TrailingBytes);
        }
        validate_checkpoint_state(&film, cx)?;
        cx.checkpoint()
            .map_err(|_| CinematicAovCheckpointError::Cancelled)?;
        Ok((
            film,
            CinematicAovCheckpointReceipt {
                content_hash,
                byte_len,
                samples_per_pixel: spp_done,
            },
        ))
    }
}

struct StreamWriter<'cx, 'scope, 'emit, E, F>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    cx: &'cx Cx<'scope>,
    emit: &'emit mut F,
    buffer: Vec<u8>,
    hasher: DomainHasher,
    emitted: u64,
}

impl<'cx, 'scope, 'emit, E, F> StreamWriter<'cx, 'scope, 'emit, E, F>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    fn try_new(
        cx: &'cx Cx<'scope>,
        emit: &'emit mut F,
    ) -> Result<Self, CinematicAovCheckpointWriteError<E>> {
        let mut buffer = Vec::new();
        buffer.try_reserve_exact(STREAM_CHUNK_BYTES).map_err(|_| {
            CinematicAovCheckpointWriteError::Checkpoint(
                CinematicAovCheckpointError::AllocationRefused,
            )
        })?;
        Ok(Self {
            cx,
            emit,
            buffer,
            hasher: DomainHasher::new(CINEMATIC_AOV_CHECKPOINT_CONTENT_DOMAIN),
            emitted: 0,
        })
    }

    fn bytes(&mut self, bytes: &[u8]) -> Result<(), CinematicAovCheckpointWriteError<E>> {
        if self.buffer.len() + bytes.len() > STREAM_CHUNK_BYTES {
            self.flush()?;
        }
        self.hasher.update(bytes);
        self.buffer.extend_from_slice(bytes);
        Ok(())
    }

    fn flush(&mut self) -> Result<(), CinematicAovCheckpointWriteError<E>> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.cx.checkpoint().map_err(|_| {
            CinematicAovCheckpointWriteError::Checkpoint(CinematicAovCheckpointError::Cancelled)
        })?;
        (self.emit)(&self.buffer).map_err(CinematicAovCheckpointWriteError::Sink)?;
        self.emitted = self.emitted.checked_add(self.buffer.len() as u64).ok_or(
            CinematicAovCheckpointWriteError::Checkpoint(
                CinematicAovCheckpointError::LengthOverflow,
            ),
        )?;
        self.buffer.clear();
        Ok(())
    }

    fn finish(mut self) -> Result<(ContentHash, u64), CinematicAovCheckpointWriteError<E>> {
        self.flush()?;
        let content_hash = self.hasher.finalize();
        self.cx.checkpoint().map_err(|_| {
            CinematicAovCheckpointWriteError::Checkpoint(CinematicAovCheckpointError::Cancelled)
        })?;
        (self.emit)(content_hash.as_bytes()).map_err(CinematicAovCheckpointWriteError::Sink)?;
        let byte_len = self.emitted.checked_add(SEAL_BYTES).ok_or(
            CinematicAovCheckpointWriteError::Checkpoint(
                CinematicAovCheckpointError::LengthOverflow,
            ),
        )?;
        // A transactional sink may request cancellation while durably writing
        // the seal. Observe it before returning the publication receipt so the
        // caller can discard that still-uncommitted artifact.
        self.cx.checkpoint().map_err(|_| {
            CinematicAovCheckpointWriteError::Checkpoint(CinematicAovCheckpointError::Cancelled)
        })?;
        Ok((content_hash, byte_len))
    }

    fn u8(&mut self, value: u8) -> Result<(), CinematicAovCheckpointWriteError<E>> {
        self.bytes(&[value])
    }

    fn u16(&mut self, value: u16) -> Result<(), CinematicAovCheckpointWriteError<E>> {
        self.bytes(&value.to_le_bytes())
    }

    fn u32(&mut self, value: u32) -> Result<(), CinematicAovCheckpointWriteError<E>> {
        self.bytes(&value.to_le_bytes())
    }

    fn u64(&mut self, value: u64) -> Result<(), CinematicAovCheckpointWriteError<E>> {
        self.bytes(&value.to_le_bytes())
    }

    fn f64(&mut self, value: f64) -> Result<(), CinematicAovCheckpointWriteError<E>> {
        self.u64(value.to_bits())
    }

    fn hash(&mut self, value: ContentHash) -> Result<(), CinematicAovCheckpointWriteError<E>> {
        self.bytes(value.as_bytes())
    }
}

fn encode_film<E, F>(
    writer: &mut StreamWriter<'_, '_, '_, E, F>,
    film: &CinematicAovFilm,
) -> Result<(), CinematicAovCheckpointWriteError<E>>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    writer.bytes(MAGIC)?;
    writer.u16(CINEMATIC_AOV_CHECKPOINT_SCHEMA_VERSION)?;
    writer.u8(film.config.profile as u8)?;
    writer.u32(CINEMATIC_AOV_SEMANTICS_VERSION)?;
    writer.u32(CINEMATIC_AOV_CATEGORY_SEMANTICS_VERSION)?;
    writer.u32(CINEMATIC_AOV_ALBEDO_SEMANTICS_VERSION)?;
    writer.u32(TRACER_BIT_SEMANTICS_VERSION)?;
    writer.u32(MOTION_TRACER_BIT_SEMANTICS_VERSION)?;
    writer.u32(CINEMATIC_CAMERA_TRACER_BIT_SEMANTICS_VERSION)?;
    writer.u32(DIELECTRIC_TRACER_BIT_SEMANTICS_VERSION)?;
    writer.u32(LIGHTING_TRACER_BIT_SEMANTICS_VERSION)?;
    writer.u32(MOTION_VECTOR_SEMANTICS_VERSION)?;
    writer.u32(crate::charts::CHART_BACKEND_BIT_SEMANTICS_VERSION)?;
    writer.u32(film.beauty.width)?;
    writer.u32(film.beauty.height)?;
    writer.u32(film.beauty.spp_done)?;
    let provenance = film.config.provenance;
    writer.u64(provenance.frame_index)?;
    writer.f64(provenance.frame_time_s)?;
    writer.f64(provenance.previous_frame_time_s)?;
    writer.f64(provenance.next_frame_time_s)?;
    writer.hash(provenance.source_trajectory_identity)?;
    writer.hash(provenance.scene_identity)?;
    writer.hash(provenance.composition_identity)?;
    let limits = film.config.limits;
    writer.u64(limits.max_pixels)?;
    writer.u64(limits.max_retained_bytes)?;
    writer.u64(limits.max_export_plane_bytes)?;
    writer.u64(limits.max_export_metadata_bytes)?;
    writer.u64(limits.max_exr_encoder_scratch_bytes)?;
    writer.u64(limits.max_encoded_exr_bytes)?;
    writer.u32(limits.max_palette_entries)?;
    writer.hash(film.config.identity)?;

    if let Some(binding) = &film.binding {
        writer.u8(1)?;
        let settings = binding.settings;
        writer.u32(settings.width)?;
        writer.u32(settings.height)?;
        writer.u32(settings.spp)?;
        writer.u32(settings.max_depth)?;
        writer.u8(sampler_tag(settings.sampler))?;
        writer.u8(strategy_tag(settings.strategy))?;
        writer.u64(settings.seed)?;
        encode_shutter(writer, binding.shutter)?;
        writer.u64(binding.shot_id)?;
        writer.u8(match binding.cut_side {
            CutSide::Before => 0,
            CutSide::After => 1,
        })?;
        writer.hash(binding.continuity_fingerprint)?;
        let (object_count, material_count) = palette_counts(&binding.palette)
            .map_err(CinematicAovCheckpointWriteError::Checkpoint)?;
        writer.u32(object_count)?;
        writer.u32(material_count)?;
        for object_id in &binding.palette.object_ids {
            writer.u64(*object_id)?;
        }
        for material_identity in &binding.palette.material_identities {
            writer.hash(*material_identity)?;
        }
    } else {
        writer.u8(0)?;
    }

    for (pixel, beauty) in film.beauty.xyz.iter().enumerate() {
        if pixel.is_multiple_of(RESTORE_POLL_PIXELS) {
            writer.cx.checkpoint().map_err(|_| {
                CinematicAovCheckpointWriteError::Checkpoint(CinematicAovCheckpointError::Cancelled)
            })?;
        }
        for value in beauty {
            writer.f64(*value)?;
        }
        if let Some(common) = &film.common {
            encode_common(writer, common[pixel])?;
        }
        if let Some(final_diagnostic) = &film.final_diagnostic {
            encode_final(writer, final_diagnostic[pixel])?;
        }
    }
    Ok(())
}

fn encode_common<E, F>(
    writer: &mut StreamWriter<'_, '_, '_, E, F>,
    value: CommonPixel,
) -> Result<(), CinematicAovCheckpointWriteError<E>>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    for component in value.albedo_sum {
        writer.f64(component)?;
    }
    for component in value.shading_normal_sum {
        writer.f64(component)?;
    }
    writer.f64(value.depth_sum_m)?;
    for component in value.previous_motion_sum_pixels {
        writer.f64(component)?;
    }
    writer.f64(value.mean_y)?;
    writer.f64(value.m2_y)?;
    writer.u32(value.accepted_count)?;
    writer.u32(value.primary_count)?;
    writer.u32(value.albedo_count)?;
    writer.u32(value.authored_shading_normal_count)?;
    writer.u32(value.previous_motion_count)
}

fn encode_final<E, F>(
    writer: &mut StreamWriter<'_, '_, '_, E, F>,
    value: FinalPixel,
) -> Result<(), CinematicAovCheckpointWriteError<E>>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    for values in [
        value.geometric_normal_sum,
        value.direct_xyz_sum,
        value.indirect_xyz_sum,
        value.emission_xyz_sum,
    ] {
        for component in values {
            writer.f64(component)?;
        }
    }
    let nearest = value.nearest_primary;
    writer.u8(u8::from(nearest.present))?;
    writer.f64(nearest.distance_squared)?;
    writer.u32(nearest.absolute_sample)?;
    writer.u64(nearest.primitive_index)?;
    writer.u32(nearest.object_palette_index)?;
    writer.u32(nearest.material_palette_index)
}

fn encode_shutter<E, F>(
    writer: &mut StreamWriter<'_, '_, '_, E, F>,
    shutter: ShutterInterval,
) -> Result<(), CinematicAovCheckpointWriteError<E>>
where
    F: FnMut(&[u8]) -> Result<(), E>,
{
    writer.f64(shutter.open_s())?;
    writer.f64(shutter.close_s())?;
    writer.u8(match shutter.convention() {
        ShutterConvention::Centered => 0,
        ShutterConvention::FrontLoaded => 1,
        ShutterConvention::BackLoaded => 2,
    })?;
    match shutter.distribution() {
        ShutterDistribution::UniformCounterV1 => {
            writer.u8(0)?;
            writer.u32(0)
        }
        ShutterDistribution::StratifiedCounterV1 { strata } => {
            writer.u8(1)?;
            writer.u32(strata)
        }
    }
}

fn checkpoint_byte_len(film: &CinematicAovFilm) -> Result<u64, CinematicAovCheckpointError> {
    let (objects, materials, settings) = if let Some(binding) = &film.binding {
        let (objects, materials) = palette_counts(&binding.palette)?;
        (objects, materials, Some(binding.settings))
    } else {
        (0, 0, None)
    };
    checkpoint_byte_len_fields(
        film.beauty.width,
        film.beauty.height,
        film.config.profile,
        settings,
        objects,
        materials,
    )
}

fn checkpoint_byte_len_fields(
    width: u32,
    height: u32,
    profile: CinematicAovProfile,
    binding: Option<Settings>,
    object_count: u32,
    material_count: u32,
) -> Result<u64, CinematicAovCheckpointError> {
    let pixels = u64::from(width)
        .checked_mul(u64::from(height))
        .ok_or(CinematicAovCheckpointError::LengthOverflow)?;
    let pixel_bytes = BEAUTY_PIXEL_BYTES
        .checked_add(if profile.has_common() {
            COMMON_PIXEL_BYTES
        } else {
            0
        })
        .and_then(|bytes| {
            bytes.checked_add(if profile.has_final() {
                FINAL_PIXEL_BYTES
            } else {
                0
            })
        })
        .ok_or(CinematicAovCheckpointError::LengthOverflow)?;
    let palette_bytes = if binding.is_some() {
        u64::from(object_count)
            .checked_mul(8)
            .and_then(|bytes| bytes.checked_add(u64::from(material_count).checked_mul(32)?))
            .ok_or(CinematicAovCheckpointError::LengthOverflow)?
    } else if object_count == 0 && material_count == 0 {
        0
    } else {
        return Err(CinematicAovCheckpointError::InvalidEnvelope);
    };
    FIXED_HEADER_BYTES
        .checked_add(if binding.is_some() {
            BOUND_HEADER_BYTES
        } else {
            0
        })
        .and_then(|bytes| bytes.checked_add(palette_bytes))
        .and_then(|bytes| bytes.checked_add(pixels.checked_mul(pixel_bytes)?))
        .and_then(|bytes| bytes.checked_add(SEAL_BYTES))
        .ok_or(CinematicAovCheckpointError::LengthOverflow)
}

fn require_exact_length(actual: u64, expected: u64) -> Result<(), CinematicAovCheckpointError> {
    match actual.cmp(&expected) {
        core::cmp::Ordering::Less => Err(CinematicAovCheckpointError::Truncated),
        core::cmp::Ordering::Greater => Err(CinematicAovCheckpointError::TrailingBytes),
        core::cmp::Ordering::Equal => Ok(()),
    }
}

fn require_expected_binding(
    actual: Option<&CinematicAovRenderBinding>,
    actual_samples_per_pixel: u32,
    expected: Option<ExpectedRenderBinding>,
) -> Result<(), CinematicAovCheckpointError> {
    let (Some(actual), Some(expected)) = (actual, expected) else {
        return if actual.is_none() && expected.is_none() {
            Ok(())
        } else {
            Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
                field: "bound render state",
            })
        };
    };
    if actual.settings != expected.settings {
        return Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "render settings",
        });
    }
    if actual.shutter != expected.shutter {
        return Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "shutter interval",
        });
    }
    if actual.shot_id != expected.shot_id {
        return Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "shot identity",
        });
    }
    if actual.cut_side != expected.cut_side {
        return Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "cut-side convention",
        });
    }
    if actual_samples_per_pixel != expected.committed_samples_per_pixel {
        return Err(CinematicAovCheckpointError::ExpectedBindingMismatch {
            field: "committed sample prefix",
        });
    }
    Ok(())
}

fn validate_checkpoint_structure(
    film: &CinematicAovFilm,
) -> Result<(), CinematicAovCheckpointError> {
    let pixel_count = checked_pixel_count(film.beauty.width, film.beauty.height)?;
    if film.beauty.xyz.len() != pixel_count
        || film
            .common
            .as_ref()
            .is_some_and(|plane| plane.len() != pixel_count)
        || film
            .final_diagnostic
            .as_ref()
            .is_some_and(|plane| plane.len() != pixel_count)
        || film.common.is_some() != film.config.profile.has_common()
        || film.final_diagnostic.is_some() != film.config.profile.has_final()
        || film.retained_bytes != retained_bytes(pixel_count, film.config.profile)?
    {
        return Err(CinematicAovCheckpointError::InvalidState {
            field: "film shape or profile",
        });
    }
    match (film.beauty.spp_done, &film.binding, film.beauty.time_mode) {
        (0, None, FilmTimeMode::Uninitialized) => {}
        (
            samples,
            Some(binding),
            FilmTimeMode::Cinematic {
                shutter,
                stream_identity,
                shot_id,
            },
        ) if samples <= binding.settings.spp
            && samples != 0
            && binding.settings.spp != 0
            && binding.settings.width == film.beauty.width
            && binding.settings.height == film.beauty.height
            && binding.shutter == shutter
            && binding.settings.seed == stream_identity
            && binding.shot_id == shot_id
            && binding.shot_id != 0
            && binding
                .continuity_fingerprint
                .as_bytes()
                .iter()
                .any(|byte| *byte != 0)
            && binding.adaptive_policy.is_none() =>
        {
            validate_reference_times(film.config, binding.shutter)?;
            validate_palette_shape(&binding.palette, film.config)?;
        }
        _ => {
            return Err(CinematicAovCheckpointError::InvalidState {
                field: "uniform render binding",
            });
        }
    }

    if film.config.profile.has_final() && film.beauty.spp_done > MAX_EXACT_F32_INTEGER {
        return Err(CinematicAovCheckpointError::Aov(
            CinematicAovError::InexactSampleCount {
                samples: film.beauty.spp_done,
            },
        ));
    }

    Ok(())
}

#[allow(clippy::float_cmp)] // canonical checkpoint zero is an exact state invariant
fn validate_checkpoint_state(
    film: &CinematicAovFilm,
    cx: &Cx<'_>,
) -> Result<(), CinematicAovCheckpointError> {
    validate_checkpoint_structure(film)?;
    if let Some(binding) = &film.binding {
        validate_palette(&binding.palette, film.config, cx)?;
    }
    let palette = film.binding.as_ref().map(|binding| &binding.palette);
    for pixel in 0..film.beauty.xyz.len() {
        if pixel.is_multiple_of(RESTORE_POLL_PIXELS) {
            cx.checkpoint()
                .map_err(|_| CinematicAovCheckpointError::Cancelled)?;
        }
        let beauty = film.beauty.xyz[pixel];
        require_canonical_finite(beauty, "beauty")?;
        if film.beauty.spp_done == 0 && beauty != [0.0; 3] {
            return Err(CinematicAovCheckpointError::InvalidState {
                field: "unbound beauty",
            });
        }
        if let Some(common) = &film.common {
            validate_common(common[pixel], film.beauty.spp_done)?;
            if let Some(final_diagnostic) = &film.final_diagnostic {
                validate_final(
                    final_diagnostic[pixel],
                    common[pixel],
                    film.beauty.spp_done,
                    palette,
                )?;
            }
        }
    }
    Ok(())
}

#[allow(clippy::float_cmp)] // canonical sums/count-zero relationships are bit-exact
fn validate_common(
    value: CommonPixel,
    expected_samples: u32,
) -> Result<(), CinematicAovCheckpointError> {
    require_canonical_finite(value.albedo_sum, "common albedo")?;
    require_canonical_finite(value.shading_normal_sum, "common shading normal")?;
    require_canonical_finite([value.depth_sum_m], "common depth")?;
    require_canonical_finite(value.previous_motion_sum_pixels, "common motion")?;
    require_canonical_finite([value.mean_y, value.m2_y], "common moments")?;
    if value.depth_sum_m < 0.0
        || (value.primary_count != 0 && value.depth_sum_m <= 0.0)
        || value.m2_y < 0.0
        || value.accepted_count != expected_samples
        || value.primary_count > value.accepted_count
        || value.albedo_count > value.primary_count
        || value.authored_shading_normal_count != 0
        || value.previous_motion_count > value.primary_count
        || (value.albedo_count == 0 && value.albedo_sum != [0.0; 3])
        || (value.primary_count == 0
            && (value.shading_normal_sum != [0.0; 3] || value.depth_sum_m != 0.0))
        || (value.previous_motion_count == 0 && value.previous_motion_sum_pixels != [0.0; 2])
        || (expected_samples < 2 && value.m2_y != 0.0)
        || (expected_samples == 0 && value != CommonPixel::EMPTY)
    {
        return Err(CinematicAovCheckpointError::InvalidState {
            field: "common pixel",
        });
    }
    Ok(())
}

#[allow(clippy::float_cmp)] // canonical absent-state vectors are exact zeros
fn validate_final(
    value: FinalPixel,
    common: CommonPixel,
    expected_samples: u32,
    palette: Option<&CinematicAovPalette>,
) -> Result<(), CinematicAovCheckpointError> {
    require_canonical_finite(value.geometric_normal_sum, "geometric normal")?;
    require_canonical_finite(value.direct_xyz_sum, "direct contribution")?;
    require_canonical_finite(value.indirect_xyz_sum, "indirect contribution")?;
    require_canonical_finite(value.emission_xyz_sum, "emission contribution")?;
    let nearest = value.nearest_primary;
    if nearest.present {
        let palette = palette.ok_or(CinematicAovCheckpointError::InvalidState {
            field: "categorical palette",
        })?;
        require_canonical_finite([nearest.distance_squared], "categorical distance")?;
        if common.primary_count == 0
            || nearest.distance_squared < 0.0
            || nearest.distance_squared > 0.5
            || nearest.absolute_sample >= expected_samples
            || usize::try_from(nearest.primitive_index).is_err()
            || nearest.material_palette_index == 0
            || nearest.material_palette_index as usize > palette.material_identities.len()
            || nearest.object_palette_index as usize > palette.object_ids.len()
        {
            return Err(CinematicAovCheckpointError::InvalidState {
                field: "categorical primary",
            });
        }
    } else if nearest != CategoricalPrimary::NONE || common.primary_count != 0 {
        return Err(CinematicAovCheckpointError::InvalidState {
            field: "absent categorical primary",
        });
    }
    if common.primary_count == 0 && value.geometric_normal_sum != [0.0; 3] {
        return Err(CinematicAovCheckpointError::InvalidState {
            field: "geometric normal without primary",
        });
    }
    if expected_samples == 0 && value != FinalPixel::EMPTY {
        return Err(CinematicAovCheckpointError::InvalidState {
            field: "empty final pixel",
        });
    }
    Ok(())
}

fn validate_palette(
    palette: &CinematicAovPalette,
    config: CinematicAovConfig,
    cx: &Cx<'_>,
) -> Result<(), CinematicAovCheckpointError> {
    validate_palette_shape(palette, config)?;
    for (index, object_id) in palette.object_ids.iter().copied().enumerate() {
        if index.is_multiple_of(RESTORE_POLL_PIXELS) {
            cx.checkpoint()
                .map_err(|_| CinematicAovCheckpointError::Cancelled)?;
        }
        if object_id == 0 || index > 0 && palette.object_ids[index - 1] >= object_id {
            return Err(CinematicAovCheckpointError::InvalidState {
                field: "identity palette",
            });
        }
    }
    for (index, identity) in palette.material_identities.iter().copied().enumerate() {
        if index.is_multiple_of(RESTORE_POLL_PIXELS) {
            cx.checkpoint()
                .map_err(|_| CinematicAovCheckpointError::Cancelled)?;
        }
        if identity.as_bytes().iter().all(|byte| *byte == 0)
            || index > 0 && palette.material_identities[index - 1] >= identity
        {
            return Err(CinematicAovCheckpointError::InvalidState {
                field: "identity palette",
            });
        }
    }
    Ok(())
}

fn validate_palette_shape(
    palette: &CinematicAovPalette,
    config: CinematicAovConfig,
) -> Result<(), CinematicAovCheckpointError> {
    let (object_count, material_count) = palette_counts(palette)?;
    validate_palette_counts(object_count, material_count, config.limits)?;
    if !config.profile.has_final()
        && (!palette.object_ids.is_empty() || !palette.material_identities.is_empty())
    {
        return Err(CinematicAovCheckpointError::InvalidState {
            field: "identity palette",
        });
    }
    Ok(())
}

fn palette_counts(
    palette: &CinematicAovPalette,
) -> Result<(u32, u32), CinematicAovCheckpointError> {
    let object_count = u32::try_from(palette.object_ids.len())
        .map_err(|_| CinematicAovCheckpointError::LengthOverflow)?;
    let material_count = u32::try_from(palette.material_identities.len())
        .map_err(|_| CinematicAovCheckpointError::LengthOverflow)?;
    Ok((object_count, material_count))
}

fn validate_palette_counts(
    objects: u32,
    materials: u32,
    limits: CinematicAovLimits,
) -> Result<(), CinematicAovCheckpointError> {
    if objects > limits.max_palette_entries
        || materials > limits.max_palette_entries
        || objects >= MAX_EXACT_F32_INTEGER
        || materials >= MAX_EXACT_F32_INTEGER
    {
        return Err(CinematicAovCheckpointError::InvalidEnvelope);
    }
    Ok(())
}

fn require_canonical_finite(
    values: impl IntoIterator<Item = f64>,
    field: &'static str,
) -> Result<(), CinematicAovCheckpointError> {
    if values
        .into_iter()
        .any(|value| !value.is_finite() || value.to_bits() == (-0.0_f64).to_bits())
    {
        Err(CinematicAovCheckpointError::InvalidState { field })
    } else {
        Ok(())
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Reader<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn take(&mut self, count: usize) -> Result<&'a [u8], CinematicAovCheckpointError> {
        let end = self
            .offset
            .checked_add(count)
            .ok_or(CinematicAovCheckpointError::LengthOverflow)?;
        let value = self
            .bytes
            .get(self.offset..end)
            .ok_or(CinematicAovCheckpointError::Truncated)?;
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8, CinematicAovCheckpointError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, CinematicAovCheckpointError> {
        Ok(u16::from_le_bytes(
            self.take(2)?
                .try_into()
                .map_err(|_| CinematicAovCheckpointError::Truncated)?,
        ))
    }

    fn u32(&mut self) -> Result<u32, CinematicAovCheckpointError> {
        Ok(u32::from_le_bytes(
            self.take(4)?
                .try_into()
                .map_err(|_| CinematicAovCheckpointError::Truncated)?,
        ))
    }

    fn u64(&mut self) -> Result<u64, CinematicAovCheckpointError> {
        Ok(u64::from_le_bytes(
            self.take(8)?
                .try_into()
                .map_err(|_| CinematicAovCheckpointError::Truncated)?,
        ))
    }

    fn canonical_f64(&mut self) -> Result<f64, CinematicAovCheckpointError> {
        let value = f64::from_bits(self.u64()?);
        if !value.is_finite() || value.to_bits() == (-0.0_f64).to_bits() {
            return Err(CinematicAovCheckpointError::InvalidEnvelope);
        }
        Ok(value)
    }

    fn hash(&mut self) -> Result<ContentHash, CinematicAovCheckpointError> {
        ContentHash::from_slice(self.take(32)?).ok_or(CinematicAovCheckpointError::Truncated)
    }

    fn f64_array<const N: usize>(&mut self) -> Result<[f64; N], CinematicAovCheckpointError> {
        let mut values = [0.0; N];
        for value in &mut values {
            *value = self.canonical_f64()?;
        }
        Ok(values)
    }

    fn ascending_nonzero_u64s(
        &mut self,
        count: u32,
        cx: &Cx<'_>,
    ) -> Result<Vec<u64>, CinematicAovCheckpointError> {
        let len =
            usize::try_from(count).map_err(|_| CinematicAovCheckpointError::LengthOverflow)?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(len)
            .map_err(|_| CinematicAovCheckpointError::AllocationRefused)?;
        for index in 0..count {
            if index.is_multiple_of(RESTORE_POLL_PIXELS as u32) {
                cx.checkpoint()
                    .map_err(|_| CinematicAovCheckpointError::Cancelled)?;
            }
            let value = self.u64()?;
            if value == 0 || values.last().is_some_and(|prior| *prior >= value) {
                return Err(CinematicAovCheckpointError::InvalidEnvelope);
            }
            values.push(value);
        }
        Ok(values)
    }

    fn ascending_nonzero_hashes(
        &mut self,
        count: u32,
        cx: &Cx<'_>,
    ) -> Result<Vec<ContentHash>, CinematicAovCheckpointError> {
        let len =
            usize::try_from(count).map_err(|_| CinematicAovCheckpointError::LengthOverflow)?;
        let mut values = Vec::new();
        values
            .try_reserve_exact(len)
            .map_err(|_| CinematicAovCheckpointError::AllocationRefused)?;
        for index in 0..count {
            if index.is_multiple_of(RESTORE_POLL_PIXELS as u32) {
                cx.checkpoint()
                    .map_err(|_| CinematicAovCheckpointError::Cancelled)?;
            }
            let value = self.hash()?;
            if value.as_bytes().iter().all(|byte| *byte == 0)
                || values.last().is_some_and(|prior| *prior >= value)
            {
                return Err(CinematicAovCheckpointError::InvalidEnvelope);
            }
            values.push(value);
        }
        Ok(values)
    }

    fn is_done(&self) -> bool {
        self.offset == self.bytes.len()
    }
}

fn decode_common(reader: &mut Reader<'_>) -> Result<CommonPixel, CinematicAovCheckpointError> {
    Ok(CommonPixel {
        albedo_sum: reader.f64_array()?,
        shading_normal_sum: reader.f64_array()?,
        depth_sum_m: reader.canonical_f64()?,
        previous_motion_sum_pixels: reader.f64_array()?,
        mean_y: reader.canonical_f64()?,
        m2_y: reader.canonical_f64()?,
        accepted_count: reader.u32()?,
        primary_count: reader.u32()?,
        albedo_count: reader.u32()?,
        authored_shading_normal_count: reader.u32()?,
        previous_motion_count: reader.u32()?,
    })
}

fn decode_final(reader: &mut Reader<'_>) -> Result<FinalPixel, CinematicAovCheckpointError> {
    let geometric_normal_sum = reader.f64_array()?;
    let direct_xyz_sum = reader.f64_array()?;
    let indirect_xyz_sum = reader.f64_array()?;
    let emission_xyz_sum = reader.f64_array()?;
    let present = match reader.u8()? {
        0 => false,
        1 => true,
        _ => return Err(CinematicAovCheckpointError::InvalidEnvelope),
    };
    Ok(FinalPixel {
        geometric_normal_sum,
        direct_xyz_sum,
        indirect_xyz_sum,
        emission_xyz_sum,
        nearest_primary: CategoricalPrimary {
            present,
            distance_squared: reader.canonical_f64()?,
            absolute_sample: reader.u32()?,
            primitive_index: reader.u64()?,
            object_palette_index: reader.u32()?,
            material_palette_index: reader.u32()?,
        },
    })
}

fn decode_profile(tag: u8) -> Result<CinematicAovProfile, CinematicAovCheckpointError> {
    match tag {
        0 => Ok(CinematicAovProfile::BeautyOnly),
        1 => Ok(CinematicAovProfile::DailyCore),
        2 => Ok(CinematicAovProfile::FinalDiagnostic),
        _ => Err(CinematicAovCheckpointError::InvalidEnvelope),
    }
}

const fn sampler_tag(sampler: Sampler) -> u8 {
    match sampler {
        Sampler::Iid => 0,
        Sampler::OwenSobol => 1,
        Sampler::OwenSobolFullPath => 2,
    }
}

fn decode_sampler(tag: u8) -> Result<Sampler, CinematicAovCheckpointError> {
    match tag {
        0 => Ok(Sampler::Iid),
        1 => Ok(Sampler::OwenSobol),
        2 => Ok(Sampler::OwenSobolFullPath),
        _ => Err(CinematicAovCheckpointError::InvalidEnvelope),
    }
}

const fn strategy_tag(strategy: DirectStrategy) -> u8 {
    match strategy {
        DirectStrategy::NeeOnly => 0,
        DirectStrategy::BsdfOnly => 1,
        DirectStrategy::Mis => 2,
        DirectStrategy::PowerMis => 3,
    }
}

fn decode_strategy(tag: u8) -> Result<DirectStrategy, CinematicAovCheckpointError> {
    match tag {
        0 => Ok(DirectStrategy::NeeOnly),
        1 => Ok(DirectStrategy::BsdfOnly),
        2 => Ok(DirectStrategy::Mis),
        3 => Ok(DirectStrategy::PowerMis),
        _ => Err(CinematicAovCheckpointError::InvalidEnvelope),
    }
}

fn decode_shutter(reader: &mut Reader<'_>) -> Result<ShutterInterval, CinematicAovCheckpointError> {
    let open_s = reader.canonical_f64()?;
    let close_s = reader.canonical_f64()?;
    let convention = match reader.u8()? {
        0 => ShutterConvention::Centered,
        1 => ShutterConvention::FrontLoaded,
        2 => ShutterConvention::BackLoaded,
        _ => return Err(CinematicAovCheckpointError::InvalidEnvelope),
    };
    let distribution_tag = reader.u8()?;
    let strata = reader.u32()?;
    let distribution = match (distribution_tag, strata) {
        (0, 0) => ShutterDistribution::UniformCounterV1,
        (1, strata) if strata != 0 => ShutterDistribution::StratifiedCounterV1 { strata },
        _ => return Err(CinematicAovCheckpointError::InvalidEnvelope),
    };
    ShutterInterval::try_from_canonical_parts(open_s, close_s, convention, distribution)
        .map_err(|_| CinematicAovCheckpointError::InvalidEnvelope)
}
