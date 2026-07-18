//! fs-rand — counter-based Philox streams keyed by LOGICAL work identity
//! (plan §6.7; Decalogue P2's seed pillar).
//!
//! The design that makes e-raced tournaments bit-reproducible and MC results
//! scheduling-independent: a draw is a pure function of
//! `(seed, kernel, tile, index)` — never of which thread ran when. Streams
//! support RANDOM ACCESS by index ([`Stream::at`]), so replay, forking, and
//! out-of-order tile execution cannot perturb randomness.
//!
//! Strict distributions are built on fs-math's deterministic functions:
//! Box–Muller normals via `det::{ln, cos}` and exponentials via `det::ln`.
//! A ziggurat normal is available as an explicit fast-mode path; strict
//! callers stay on Box–Muller until the cross-ISA admission proof lands.
//!
//! Field widths (documented contract): seed 64 bits (Philox key), tile id
//! 32 bits, kernel id 32 bits (together counter words 2–3), draw index 64
//! bits (counter words 0–1). 2⁶⁴ draws per (seed, kernel, tile) stream.

pub mod cbc;
pub mod cbc_cert;
pub mod cbc_exec;
pub mod dist;
pub mod philox;
pub mod qmc;
pub mod ziggurat;

use fs_math::det;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// STREAM-SEMANTICS VERSION (bead y4pt): bump on ANY change that can
/// move the bits a downstream consumer draws from a given
/// (seed, kernel, tile, index) — counter advancement, key mapping,
/// Philox rounds, or distribution transforms. Downstream goldens
/// declare the version they were frozen against in
/// golden-couplings.json; `cargo run -p xtask -- check-goldens` fails
/// on drift until every dependent golden is deliberately re-frozen.
pub const STREAM_SEMANTICS_VERSION: u32 = 1;

/// Domain name for the complete seed/kernel/tile/index Philox mapping.
pub const STREAM_POSITION_IDENTITY_DOMAIN: &str = "org.frankensim.fs-rand.stream-position.v1";

/// Owner-local stream-position declaration consumed by
/// `xtask check-identities`.
pub const STREAM_POSITION_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-rand:stream-position",
    "version_const=STREAM_SEMANTICS_VERSION",
    "version=1",
    "domain=org.frankensim.fs-rand.stream-position.v1",
    "domain_const=STREAM_POSITION_IDENTITY_DOMAIN",
    "encoder=encode_stream_position",
    "encoder_helpers=none",
    "schema_constants=STREAM_SEMANTICS_VERSION,STREAM_POSITION_IDENTITY_DOMAIN,crates/fs-rand/src/philox.rs#M0,crates/fs-rand/src/philox.rs#M1,crates/fs-rand/src/philox.rs#W0,crates/fs-rand/src/philox.rs#W1",
    "schema_functions=crates/fs-rand/src/philox.rs#mulhilo,crates/fs-rand/src/philox.rs#round,crates/fs-rand/src/philox.rs#philox4x32_10,crates/fs-rand/src/philox.rs#philox4x32_10_batch",
    "schema_dependencies=none",
    "digest=philox4x32_10",
    "encoding=fixed-width-key",
    "sources=StreamPositionIdentityInput,StreamKey",
    "source_fields=StreamPositionIdentityInput.key:derived:nested-key-fields-encoded-separately,StreamPositionIdentityInput.index:semantic,StreamKey.seed:semantic,StreamKey.kernel:semantic,StreamKey.tile:semantic",
    "source_bindings=StreamPositionIdentityInput.index>index-low+index-high,StreamKey.seed>seed-low+seed-high,StreamKey.kernel>kernel,StreamKey.tile>tile",
    "external_semantic_fields=none",
    "semantic_fields=seed-low,seed-high,kernel,tile,index-low,index-high",
    "excluded_fields=worker:execution-order-only,thread:execution-order-only,schedule:execution-order-only",
    "consumers=Stream::at,Stream::next_u64,Stream::next_f64,Stream::next_below,Stream::next_normal,Stream::next_normal_ziggurat,Stream::next_exponential,Stream::fill_u64,Stream::fill_f64,fs-rand::dist,fs-rand::qmc,fs-rand:distribution-stream-golden,fs-rand:ziggurat-stream-golden",
    "mutations=seed-low:crates/fs-rand/src/lib.rs#stream_identity_mutation_battery,seed-high:crates/fs-rand/src/lib.rs#stream_identity_mutation_battery,kernel:crates/fs-rand/src/lib.rs#stream_identity_mutation_battery,tile:crates/fs-rand/src/lib.rs#stream_identity_mutation_battery,index-low:crates/fs-rand/src/lib.rs#stream_identity_mutation_battery,index-high:crates/fs-rand/src/lib.rs#stream_identity_mutation_battery",
    "nonsemantic_mutations=worker:crates/fs-rand/src/lib.rs#worker_shuffle_invariance,thread:crates/fs-rand/src/lib.rs#worker_shuffle_invariance,schedule:crates/fs-rand/src/lib.rs#worker_shuffle_invariance",
    "field_guard=classify_stream_key_identity_fields",
    "transport_guard=philox_position_words",
    "version_guard=crates/fs-rand/src/lib.rs#stale_checkpoint_versions_fail_closed",
    "coupling_surface=fs-rand:stream-semantics",
];

/// Version of the retained stream-checkpoint transport. A checkpoint is replay
/// authority only after both this transport version and
/// [`STREAM_SEMANTICS_VERSION`] are accepted exactly.
pub const STREAM_CHECKPOINT_VERSION: u32 = 1;

/// Domain name for retained stream checkpoint transport and replay admission.
pub const STREAM_CHECKPOINT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-rand.stream-checkpoint.v1";

/// Type tag at the start of every canonical retained stream checkpoint.
///
/// This refuses accidental cross-type decoding before any retained values are
/// admitted as replay authority.
pub const STREAM_CHECKPOINT_MAGIC: [u8; 8] = *b"FSRCKPT\0";

/// Exact width of [`StreamCheckpoint::to_canonical_le_bytes`].
///
/// Layout: magic (8), domain (43), checkpoint version (4), stream-semantics
/// version (4), seed (8), kernel (4), tile (4), and next index (8).
pub const STREAM_CHECKPOINT_CANONICAL_LEN: usize = STREAM_CHECKPOINT_MAGIC.len()
    + STREAM_CHECKPOINT_IDENTITY_DOMAIN.len()
    + core::mem::size_of::<u32>() * 4
    + core::mem::size_of::<u64>() * 2;

/// Owner-local checkpoint declaration consumed by `xtask check-identities`.
pub const STREAM_CHECKPOINT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-rand:stream-checkpoint",
    "version_const=STREAM_CHECKPOINT_VERSION",
    "version=1",
    "domain=org.frankensim.fs-rand.stream-checkpoint.v1",
    "domain_const=STREAM_CHECKPOINT_IDENTITY_DOMAIN",
    "encoder=StreamCheckpoint::to_canonical_le_bytes",
    "encoder_helpers=none",
    "schema_constants=STREAM_CHECKPOINT_VERSION,STREAM_CHECKPOINT_IDENTITY_DOMAIN,STREAM_CHECKPOINT_MAGIC,STREAM_CHECKPOINT_CANONICAL_LEN,STREAM_SEMANTICS_VERSION",
    "schema_functions=none",
    "schema_dependencies=fs-rand:stream-position",
    "digest=none-exact-canonical-transport",
    "encoding=canonical-transport-exact-bits",
    "sources=StreamCheckpoint,StreamKey",
    "source_fields=StreamCheckpoint.checkpoint_version:semantic,StreamCheckpoint.stream_semantics_version:semantic,StreamCheckpoint.key:derived:nested-key-fields-encoded-separately,StreamCheckpoint.index:semantic,StreamKey.seed:semantic,StreamKey.kernel:semantic,StreamKey.tile:semantic",
    "source_bindings=StreamCheckpoint.checkpoint_version>checkpoint-version,StreamCheckpoint.stream_semantics_version>stream-semantics-version,StreamCheckpoint.index>index,StreamKey.seed>seed,StreamKey.kernel>kernel,StreamKey.tile>tile",
    "external_semantic_fields=magic,domain",
    "semantic_fields=magic,domain,checkpoint-version,stream-semantics-version,seed,kernel,tile,index",
    "excluded_fields=none",
    "consumers=Stream::checkpoint,StreamCheckpoint::from_canonical_le_bytes,Stream::resume,Stream::resume_retained,resumable-samplers",
    "mutations=magic:crates/fs-rand/src/lib.rs#checkpoint_transport_mutation_battery,domain:crates/fs-rand/src/lib.rs#checkpoint_transport_mutation_battery,checkpoint-version:crates/fs-rand/src/lib.rs#checkpoint_transport_mutation_battery,stream-semantics-version:crates/fs-rand/src/lib.rs#checkpoint_transport_mutation_battery,seed:crates/fs-rand/src/lib.rs#checkpoint_transport_mutation_battery,kernel:crates/fs-rand/src/lib.rs#checkpoint_transport_mutation_battery,tile:crates/fs-rand/src/lib.rs#checkpoint_transport_mutation_battery,index:crates/fs-rand/src/lib.rs#checkpoint_transport_mutation_battery",
    "nonsemantic_mutations=none",
    "field_guard=classify_stream_checkpoint_identity_fields",
    "transport_guard=StreamCheckpoint::from_canonical_le_bytes",
    "version_guard=crates/fs-rand/src/lib.rs#stale_checkpoint_versions_fail_closed",
    "coupling_surface=fs-rand:stream-checkpoint",
];

/// The logical identity of a stream: the Cx-carried seed/kernel/tile key.
/// fs-exec's separate iteration/generation axis is not the within-stream draw
/// index and is refused by [`StreamKey::from_exec_parts`] when nonzero.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StreamKey {
    /// Study seed (one of the Five Explicits).
    pub seed: u64,
    /// Kernel identity (registry-assigned; stable across runs).
    pub kernel: u32,
    /// Logical tile identity (NOT the worker/thread id — the whole point).
    pub tile: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct StreamPositionIdentityInput {
    key: StreamKey,
    index: u64,
}

#[allow(dead_code)]
fn classify_stream_key_identity_fields(input: &StreamPositionIdentityInput, key: &StreamKey) {
    let StreamPositionIdentityInput { key: _, index: _ } = input;
    let StreamKey {
        seed: _,
        kernel: _,
        tile: _,
    } = key;
}

/// Versioned replay state for one logical stream.
///
/// The fields are public because checkpoints are transport data, not trusted
/// authority. [`Stream::resume`] validates both version fields before using the
/// key or index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StreamCheckpoint {
    /// Checkpoint transport schema version.
    pub checkpoint_version: u32,
    /// Bit-semantics version under which subsequent draws were defined.
    pub stream_semantics_version: u32,
    /// Complete logical stream identity.
    pub key: StreamKey,
    /// Next counter position to draw.
    pub index: u64,
}

#[allow(dead_code)]
fn classify_stream_checkpoint_identity_fields(checkpoint: &StreamCheckpoint) {
    let StreamCheckpoint {
        checkpoint_version: _,
        stream_semantics_version: _,
        key,
        index: _,
    } = checkpoint;
    let StreamKey {
        seed: _,
        kernel: _,
        tile: _,
    } = key;
}

impl StreamCheckpoint {
    /// Construct a checkpoint under this build's exact transport and stream
    /// semantics. Untrusted decoders may construct the public fields directly;
    /// resume still validates them.
    #[must_use]
    pub const fn current(key: StreamKey, index: u64) -> Self {
        Self {
            checkpoint_version: STREAM_CHECKPOINT_VERSION,
            stream_semantics_version: STREAM_SEMANTICS_VERSION,
            key,
            index,
        }
    }

    /// Encode the complete replay identity as one fixed-width canonical frame.
    ///
    /// Every integer is little-endian. The magic and domain are part of the
    /// bytes so another retained type or identity domain cannot be decoded as
    /// an fs-rand checkpoint accidentally.
    #[must_use]
    pub fn to_canonical_le_bytes(self) -> [u8; STREAM_CHECKPOINT_CANONICAL_LEN] {
        const MAGIC_END: usize = STREAM_CHECKPOINT_MAGIC.len();
        const DOMAIN_END: usize = MAGIC_END + STREAM_CHECKPOINT_IDENTITY_DOMAIN.len();
        const CHECKPOINT_VERSION_END: usize = DOMAIN_END + core::mem::size_of::<u32>();
        const STREAM_VERSION_END: usize = CHECKPOINT_VERSION_END + core::mem::size_of::<u32>();
        const SEED_END: usize = STREAM_VERSION_END + core::mem::size_of::<u64>();
        const KERNEL_END: usize = SEED_END + core::mem::size_of::<u32>();
        const TILE_END: usize = KERNEL_END + core::mem::size_of::<u32>();

        let mut bytes = [0; STREAM_CHECKPOINT_CANONICAL_LEN];
        bytes[..MAGIC_END].copy_from_slice(&STREAM_CHECKPOINT_MAGIC);
        bytes[MAGIC_END..DOMAIN_END].copy_from_slice(STREAM_CHECKPOINT_IDENTITY_DOMAIN.as_bytes());
        bytes[DOMAIN_END..CHECKPOINT_VERSION_END]
            .copy_from_slice(&self.checkpoint_version.to_le_bytes());
        bytes[CHECKPOINT_VERSION_END..STREAM_VERSION_END]
            .copy_from_slice(&self.stream_semantics_version.to_le_bytes());
        bytes[STREAM_VERSION_END..SEED_END].copy_from_slice(&self.key.seed.to_le_bytes());
        bytes[SEED_END..KERNEL_END].copy_from_slice(&self.key.kernel.to_le_bytes());
        bytes[KERNEL_END..TILE_END].copy_from_slice(&self.key.tile.to_le_bytes());
        bytes[TILE_END..].copy_from_slice(&self.index.to_le_bytes());
        bytes
    }

    /// Decode and admit one exact canonical retained checkpoint.
    ///
    /// This is deliberately not a permissive historical parser: the input
    /// must have the exact fixed width, magic, domain, checkpoint version, and
    /// stream-semantics version supported by this build. In particular,
    /// trailing bytes are refused rather than ignored.
    ///
    /// # Errors
    /// [`StreamReplayError`] naming the first transport or version admission
    /// rule that failed.
    pub fn from_canonical_le_bytes(bytes: &[u8]) -> Result<Self, StreamReplayError> {
        const MAGIC_END: usize = STREAM_CHECKPOINT_MAGIC.len();
        const DOMAIN_END: usize = MAGIC_END + STREAM_CHECKPOINT_IDENTITY_DOMAIN.len();
        const CHECKPOINT_VERSION_END: usize = DOMAIN_END + core::mem::size_of::<u32>();
        const STREAM_VERSION_END: usize = CHECKPOINT_VERSION_END + core::mem::size_of::<u32>();
        const SEED_END: usize = STREAM_VERSION_END + core::mem::size_of::<u64>();
        const KERNEL_END: usize = SEED_END + core::mem::size_of::<u32>();
        const TILE_END: usize = KERNEL_END + core::mem::size_of::<u32>();

        if bytes.len() != STREAM_CHECKPOINT_CANONICAL_LEN {
            return Err(StreamReplayError::InvalidCheckpointLength {
                actual: bytes.len(),
                expected: STREAM_CHECKPOINT_CANONICAL_LEN,
            });
        }
        if &bytes[..MAGIC_END] != STREAM_CHECKPOINT_MAGIC.as_slice() {
            return Err(StreamReplayError::InvalidCheckpointMagic);
        }
        if &bytes[MAGIC_END..DOMAIN_END] != STREAM_CHECKPOINT_IDENTITY_DOMAIN.as_bytes() {
            return Err(StreamReplayError::InvalidCheckpointDomain);
        }

        let checkpoint = Self {
            checkpoint_version: u32::from_le_bytes(
                bytes[DOMAIN_END..CHECKPOINT_VERSION_END]
                    .try_into()
                    .expect("checkpoint version has a fixed four-byte slice"),
            ),
            stream_semantics_version: u32::from_le_bytes(
                bytes[CHECKPOINT_VERSION_END..STREAM_VERSION_END]
                    .try_into()
                    .expect("stream version has a fixed four-byte slice"),
            ),
            key: StreamKey {
                seed: u64::from_le_bytes(
                    bytes[STREAM_VERSION_END..SEED_END]
                        .try_into()
                        .expect("seed has a fixed eight-byte slice"),
                ),
                kernel: u32::from_le_bytes(
                    bytes[SEED_END..KERNEL_END]
                        .try_into()
                        .expect("kernel has a fixed four-byte slice"),
                ),
                tile: u32::from_le_bytes(
                    bytes[KERNEL_END..TILE_END]
                        .try_into()
                        .expect("tile has a fixed four-byte slice"),
                ),
            },
            index: u64::from_le_bytes(
                bytes[TILE_END..]
                    .try_into()
                    .expect("index has a fixed eight-byte slice"),
            ),
        };
        checkpoint.validate_versions()?;
        Ok(checkpoint)
    }

    fn validate_versions(self) -> Result<(), StreamReplayError> {
        if self.checkpoint_version != STREAM_CHECKPOINT_VERSION {
            return Err(StreamReplayError::UnknownCheckpointVersion {
                declared: self.checkpoint_version,
                supported: STREAM_CHECKPOINT_VERSION,
            });
        }
        if self.stream_semantics_version != STREAM_SEMANTICS_VERSION {
            return Err(StreamReplayError::UnknownStreamSemanticsVersion {
                declared: self.stream_semantics_version,
                supported: STREAM_SEMANTICS_VERSION,
            });
        }
        Ok(())
    }
}

/// Why retained stream state cannot be replayed by this build.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamReplayError {
    /// The retained frame was truncated or carried trailing bytes.
    InvalidCheckpointLength {
        /// Width supplied by the caller.
        actual: usize,
        /// One exact width accepted by this build.
        expected: usize,
    },
    /// The retained frame is not tagged as an fs-rand stream checkpoint.
    InvalidCheckpointMagic,
    /// The retained frame belongs to another identity domain.
    InvalidCheckpointDomain,
    /// The checkpoint transport schema is unknown.
    UnknownCheckpointVersion {
        /// Version declared by the retained checkpoint.
        declared: u32,
        /// Exact version supported by this build.
        supported: u32,
    },
    /// The checkpoint was produced under different draw semantics.
    UnknownStreamSemanticsVersion {
        /// Version declared by the retained checkpoint.
        declared: u32,
        /// Exact version supported by this build.
        supported: u32,
    },
}

impl core::fmt::Display for StreamReplayError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidCheckpointLength { actual, expected } => write!(
                f,
                "stream checkpoint has {actual} bytes; canonical transport requires exactly {expected} and refuses truncation or trailing data"
            ),
            Self::InvalidCheckpointMagic => {
                f.write_str("stream checkpoint magic is invalid; retained type refused")
            }
            Self::InvalidCheckpointDomain => f.write_str(
                "stream checkpoint identity domain is invalid; cross-domain replay refused",
            ),
            Self::UnknownCheckpointVersion {
                declared,
                supported,
            } => write!(
                f,
                "stream checkpoint schema v{declared} is unsupported; this build accepts exactly v{supported}"
            ),
            Self::UnknownStreamSemanticsVersion {
                declared,
                supported,
            } => write!(
                f,
                "stream semantics v{declared} are unsupported; this build accepts exactly v{supported} and refuses to guess replay bits"
            ),
        }
    }
}

impl core::error::Error for StreamReplayError {}

/// Version of the fs-exec → fs-rand key bridge contract (field widths
/// and refusal rules below). Bump ONLY with a recorded justification —
/// replayability of ledgered keys depends on it.
pub const EXEC_KEY_BRIDGE_VERSION: u32 = 1;

/// Why an fs-exec logical key cannot become an fs-rand [`StreamKey`]
/// (bead wf9.7.1): the bridge REFUSES rather than truncates, because a
/// silent truncation would let two distinct logical streams collide.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecKeyBridgeError {
    /// `kernel_id` exceeds this key's u32 kernel slot.
    KernelOverflow {
        /// The offending value.
        kernel_id: u64,
    },
    /// `tile` exceeds this key's u32 tile slot.
    TileOverflow {
        /// The offending value.
        tile: u64,
    },
    /// fs-exec's iteration/generation axis has NO slot here: fs-rand's
    /// draw index is the WITHIN-stream counter, not an identity axis.
    /// Callers with generation-diverging streams must ledger the
    /// generation into the seed (e.g. fs-exec's `key128` path) rather
    /// than silently folding it.
    IterationUnrepresentable {
        /// The offending value.
        iteration: u64,
    },
}

impl core::fmt::Display for ExecKeyBridgeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ExecKeyBridgeError::KernelOverflow { kernel_id } => write!(
                f,
                "exec kernel_id {kernel_id} exceeds the u32 kernel slot (bridge v{EXEC_KEY_BRIDGE_VERSION}); refused rather than truncated"
            ),
            ExecKeyBridgeError::TileOverflow { tile } => write!(
                f,
                "exec tile {tile} exceeds the u32 tile slot (bridge v{EXEC_KEY_BRIDGE_VERSION}); refused rather than truncated"
            ),
            ExecKeyBridgeError::IterationUnrepresentable { iteration } => write!(
                f,
                "exec iteration {iteration} has no slot in fs-rand's key (bridge v{EXEC_KEY_BRIDGE_VERSION}): the draw index is a counter, not identity — ledger the generation into the seed instead"
            ),
        }
    }
}

impl core::error::Error for ExecKeyBridgeError {}

impl StreamKey {
    /// Open the stream at index 0.
    #[must_use]
    pub fn stream(self) -> Stream {
        Stream {
            key: self,
            index: 0,
        }
    }

    /// CHECKED bridge from fs-exec's four-u64 logical key fields
    /// (`seed`, `kernel_id`, `tile`, `iteration` — bead wf9.7.1,
    /// bridge v[`EXEC_KEY_BRIDGE_VERSION`]). Field-width contract:
    /// seed is lossless (u64 → u64); kernel and tile must fit their
    /// u32 slots; iteration must be 0 (no identity slot exists —
    /// see [`ExecKeyBridgeError::IterationUnrepresentable`]).
    /// Refusal, never truncation: replay must reconstruct the SAME
    /// stream from ledgered fields, so a lossy mapping is a collision
    /// generator, not a bridge.
    ///
    /// # Errors
    /// [`ExecKeyBridgeError`] naming the unrepresentable field.
    pub fn from_exec_parts(
        seed: u64,
        kernel_id: u64,
        tile: u64,
        iteration: u64,
    ) -> Result<StreamKey, ExecKeyBridgeError> {
        let kernel = u32::try_from(kernel_id)
            .map_err(|_| ExecKeyBridgeError::KernelOverflow { kernel_id })?;
        let tile32 = u32::try_from(tile).map_err(|_| ExecKeyBridgeError::TileOverflow { tile })?;
        if iteration != 0 {
            return Err(ExecKeyBridgeError::IterationUnrepresentable { iteration });
        }
        Ok(StreamKey {
            seed,
            kernel,
            tile: tile32,
        })
    }
}

/// A sequential view over the counter-based generator. `Copy` is deliberate:
/// forking a stream is just copying it (forks that must diverge should use
/// distinct tile/kernel ids instead — divergence by IDENTITY, not by state).
#[derive(Debug, Clone, Copy)]
pub struct Stream {
    key: StreamKey,
    index: u64,
}

const fn philox_position_words(key: StreamKey, index: u64) -> ([u32; 4], [u32; 2]) {
    encode_stream_position(StreamPositionIdentityInput { key, index })
}

const fn encode_stream_position(input: StreamPositionIdentityInput) -> ([u32; 4], [u32; 2]) {
    let StreamPositionIdentityInput { key, index } = input;
    let StreamKey { seed, kernel, tile } = key;
    (
        [index as u32, (index >> 32) as u32, tile, kernel],
        [seed as u32, (seed >> 32) as u32],
    )
}

impl Stream {
    /// RANDOM ACCESS: the 128 output bits at `index`, independent of any
    /// sequential position. The foundation of replay and shuffle-invariance.
    #[must_use]
    pub fn at(key: StreamKey, index: u64) -> [u32; 4] {
        let (counter, key_words) = philox_position_words(key, index);
        philox::philox4x32_10(counter, key_words)
    }

    /// Current index (for provenance records / resumable checkpoints).
    #[must_use]
    pub fn index(&self) -> u64 {
        self.index
    }

    /// Capture versioned state for provenance, persistence, or migration.
    #[must_use]
    pub const fn checkpoint(&self) -> StreamCheckpoint {
        StreamCheckpoint::current(self.key, self.index)
    }

    /// Resume only an exactly supported checkpoint.
    ///
    /// # Errors
    /// [`StreamReplayError`] if either the checkpoint transport or draw
    /// semantics version is stale or unknown. Replay never guesses how an old
    /// key/index pair should be interpreted.
    pub fn resume(checkpoint: StreamCheckpoint) -> Result<Stream, StreamReplayError> {
        checkpoint.validate_versions()?;
        Ok(Stream {
            key: checkpoint.key,
            index: checkpoint.index,
        })
    }

    /// Decode and resume one canonical retained checkpoint.
    ///
    /// # Errors
    /// [`StreamReplayError`] if framing, domain, length, or either version is
    /// not accepted exactly.
    pub fn resume_retained(bytes: &[u8]) -> Result<Stream, StreamReplayError> {
        Self::resume(StreamCheckpoint::from_canonical_le_bytes(bytes)?)
    }

    /// Next 64 uniform bits.
    #[must_use]
    pub fn next_u64(&mut self) -> u64 {
        let block = Self::at(self.key, self.index);
        self.index = self.index.wrapping_add(1);
        (u64::from(block[1]) << 32) | u64::from(block[0])
    }

    /// Uniform in [0, 1) with 53 random bits (the standard exact ladder).
    #[must_use]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0) // 2⁻⁵³
    }

    /// Uniform integer in [0, n) via Lemire's widening-multiply method with
    /// the DETERMINISTIC rejection contract: rejected draws advance the
    /// index like any other draw, so the consumed-count is a pure function
    /// of the stream content (replay-safe).
    #[must_use]
    pub fn next_below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "next_below(0) is meaningless");
        loop {
            let x = self.next_u64();
            let m = u128::from(x) * u128::from(n);
            let lo = m as u64;
            if lo >= n.wrapping_neg() % n {
                return (m >> 64) as u64;
            }
            // rejected: index already advanced; try the next block.
        }
    }

    /// Standard normal via Box–Muller on fs-math strict functions —
    /// cross-ISA deterministic sampled values. Consumes exactly 2 draws.
    #[must_use]
    pub fn next_normal(&mut self) -> f64 {
        // u ∈ (0,1]: guard the log; v ∈ [0,1).
        let u = 1.0 - self.next_f64();
        let v = self.next_f64();
        det::sqrt(-2.0 * det::ln(u)) * det::cos(2.0 * std::f64::consts::PI * v)
    }

    /// Standard normal via the ZIGGURAT (bead 1za9) — the FAST-MODE-ONLY perf
    /// path. Deterministic table + deterministic rejection consumption, but not
    /// admitted to strict mode until a cross-ISA bitwise proof lands; strict
    /// callers use [`Stream::next_normal`] (Box–Muller). See [`ziggurat`].
    #[must_use]
    pub fn next_normal_ziggurat(&mut self) -> f64 {
        ziggurat::normal(self)
    }

    /// Exponential(1) via inversion (consumes exactly 1 draw).
    #[must_use]
    pub fn next_exponential(&mut self) -> f64 {
        -det::ln(1.0 - self.next_f64())
    }

    /// The `L` counter blocks for indices `[base, base+L)` under this stream's
    /// key (all blocks share the key; consecutive counters differ in the low
    /// words) — the bulk-generation primitive.
    fn blocks_from<const L: usize>(&self, base: u64) -> [[u32; 4]; L] {
        let ctr: [[u32; 4]; L] = core::array::from_fn(|l| {
            let idx = base.wrapping_add(l as u64);
            philox_position_words(self.key, idx).0
        });
        let key_words = philox_position_words(self.key, base).1;
        philox::philox4x32_10_batch::<L>(&ctr, key_words)
    }

    /// BULK-fill a slice with uniform `[0,1)` values via 8-lane batched
    /// generation (auto-vectorizable), then a scalar tail. BITWISE-IDENTICAL to
    /// `out.len()` sequential [`Stream::next_f64`] calls, and the index advances
    /// by exactly `out.len()` (replay-safe).
    pub fn fill_f64(&mut self, out: &mut [f64]) {
        const L: usize = 8;
        let (chunks, tail) = out.as_chunks_mut::<L>();
        for chunk in chunks {
            let blocks = self.blocks_from::<L>(self.index);
            for (o, b) in chunk.iter_mut().zip(&blocks) {
                let u = (u64::from(b[1]) << 32) | u64::from(b[0]);
                *o = (u >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0); // 2⁻⁵³
            }
            self.index = self.index.wrapping_add(L as u64);
        }
        for o in tail {
            *o = self.next_f64();
        }
    }

    /// BULK-fill a slice with uniform 64-bit words (same batching + bitwise
    /// equivalence to sequential [`Stream::next_u64`]).
    pub fn fill_u64(&mut self, out: &mut [u64]) {
        const L: usize = 8;
        let (chunks, tail) = out.as_chunks_mut::<L>();
        for chunk in chunks {
            let blocks = self.blocks_from::<L>(self.index);
            for (o, b) in chunk.iter_mut().zip(&blocks) {
                *o = (u64::from(b[1]) << 32) | u64::from(b[0]);
            }
            self.index = self.index.wrapping_add(L as u64);
        }
        for o in tail {
            *o = self.next_u64();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const KEY: StreamKey = StreamKey {
        seed: 0x5EED_0001_DEAD_BEEF,
        kernel: 7,
        tile: 42,
    };

    const CHECKPOINT_CANONICAL_LE_KAT: [u8; STREAM_CHECKPOINT_CANONICAL_LEN] =
        *b"FSRCKPT\0org.frankensim.fs-rand.stream-checkpoint.v1\
          \x01\x00\x00\x00\x01\x00\x00\x00\xef\xcd\xab\x89\x67\x45\x23\x01\
          \xdf\x9b\x57\x13\xe0\xac\x68\x24\x80\x70\x60\x50\x40\x30\x20\x10";

    const CHECKPOINT_KAT: StreamCheckpoint = StreamCheckpoint::current(
        StreamKey {
            seed: 0x0123_4567_89ab_cdef,
            kernel: 0x1357_9bdf,
            tile: 0x2468_ace0,
        },
        0x1020_3040_5060_7080,
    );

    fn assert_only_transport_range_changed(
        base: &[u8; STREAM_CHECKPOINT_CANONICAL_LEN],
        changed: &[u8; STREAM_CHECKPOINT_CANONICAL_LEN],
        start: usize,
        end: usize,
        field: &str,
    ) {
        assert_eq!(
            &base[..start],
            &changed[..start],
            "{field} mutation moved earlier transport bytes"
        );
        assert_ne!(
            &base[start..end],
            &changed[start..end],
            "{field} mutation did not move its transport bytes"
        );
        assert_eq!(
            &base[end..],
            &changed[end..],
            "{field} mutation moved later transport bytes"
        );
    }

    #[test]
    fn random_access_matches_sequential() {
        let mut s = KEY.stream();
        let seq: Vec<u64> = (0..64).map(|_| s.next_u64()).collect();
        for (i, &want) in seq.iter().enumerate() {
            let block = Stream::at(KEY, i as u64);
            let got = (u64::from(block[1]) << 32) | u64::from(block[0]);
            assert_eq!(got, want, "random access diverged at index {i}");
        }
    }

    #[test]
    fn stream_at_position_mapping_known_answer() {
        // The nonzero Random123 vector pins fs-rand's public word mapping:
        // index -> counter[0..=1], tile -> counter[2], kernel -> counter[3],
        // seed -> key[0..=1], all low word first.
        let key = StreamKey {
            seed: 0x299f_31d0_a409_3822,
            kernel: 0x0370_7344,
            tile: 0x1319_8a2e,
        };
        let index = 0x85a3_08d3_243f_6a88;
        assert_eq!(
            Stream::at(key, index),
            [0xd16c_fe09, 0x94fd_cceb, 0x5001_e420, 0x2412_6ea1]
        );
    }

    #[test]
    fn worker_shuffle_invariance() {
        // Simulate tiles executed in three different worker orders; each
        // tile's draws must be identical regardless (the P2 property that
        // makes results independent of scheduling).
        let tiles: Vec<u32> = (0..16).collect();
        let draw_tile = |tile: u32| -> Vec<f64> {
            let mut s = StreamKey {
                seed: 1234,
                kernel: 3,
                tile,
            }
            .stream();
            (0..32).map(|_| s.next_f64()).collect()
        };
        let baseline: Vec<Vec<f64>> = tiles.iter().map(|&t| draw_tile(t)).collect();
        for order in [
            tiles.iter().rev().copied().collect::<Vec<_>>(),
            tiles
                .iter()
                .step_by(2)
                .chain(tiles.iter().skip(1).step_by(2))
                .copied()
                .collect(),
        ] {
            for &t in &order {
                let redo = draw_tile(t);
                assert!(
                    redo.iter()
                        .zip(&baseline[t as usize])
                        .all(|(a, b)| a.to_bits() == b.to_bits()),
                    "tile {t} draws depended on execution order"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-rand\",\"case\":\"shuffle-invariance\",\"verdict\":\"pass\",\"detail\":\"16 tiles x 3 orders bitwise identical\"}}"
        );
    }

    #[test]
    fn streams_with_different_identities_are_uncorrelated() {
        // Crude cross-correlation check between adjacent tiles/kernels/seeds.
        let corr = |a: StreamKey, b: StreamKey| -> f64 {
            let (mut sa, mut sb) = (a.stream(), b.stream());
            let n = 4096;
            let (mut ma, mut mb, mut cov, mut va, mut vb) = (0.0, 0.0, 0.0, 0.0, 0.0);
            let xs: Vec<(f64, f64)> = (0..n).map(|_| (sa.next_f64(), sb.next_f64())).collect();
            for &(x, y) in &xs {
                ma += x;
                mb += y;
            }
            ma /= f64::from(n);
            mb /= f64::from(n);
            for &(x, y) in &xs {
                cov += (x - ma) * (y - mb);
                va += (x - ma) * (x - ma);
                vb += (y - mb) * (y - mb);
            }
            cov / (va.sqrt() * vb.sqrt())
        };
        for (a, b) in [
            (KEY, StreamKey { tile: 43, ..KEY }),
            (KEY, StreamKey { kernel: 8, ..KEY }),
            (
                KEY,
                StreamKey {
                    seed: KEY.seed ^ 1,
                    ..KEY
                },
            ),
        ] {
            let c = corr(a, b);
            assert!(c.abs() < 0.06, "adjacent identities correlate: {c}");
        }
    }

    #[test]
    fn uniform_chi_square_and_moments() {
        const BINS: usize = 64;
        const N: usize = 64 * 1024;
        let mut s = KEY.stream();
        let mut counts = [0u32; BINS];
        let mut mean = 0.0;
        for _ in 0..N {
            let x = s.next_f64();
            assert!((0.0..1.0).contains(&x));
            counts[(x * BINS as f64) as usize] += 1;
            mean += x;
        }
        mean /= N as f64;
        let expect = (N / BINS) as f64;
        let chi2: f64 = counts
            .iter()
            .map(|&c| (f64::from(c) - expect).powi(2) / expect)
            .sum();
        // 63 dof: mean 63, sd ~11.2; accept within ±5 sd.
        assert!((10.0..=120.0).contains(&chi2), "chi2 {chi2} out of band");
        assert!((mean - 0.5).abs() < 0.005, "mean {mean}");
    }

    #[test]
    fn normal_and_exponential_moments() {
        const N: usize = 200_000;
        let mut s = KEY.stream();
        let (mut m1, mut m2, mut m4) = (0.0, 0.0, 0.0);
        for _ in 0..N {
            let z = s.next_normal();
            m1 += z;
            m2 += z * z;
            m4 += z * z * z * z;
        }
        let n = N as f64;
        assert!((m1 / n).abs() < 0.01, "normal mean {}", m1 / n);
        assert!((m2 / n - 1.0).abs() < 0.02, "normal var {}", m2 / n);
        assert!((m4 / n - 3.0).abs() < 0.12, "normal kurtosis {}", m4 / n);
        let (mut e1, mut e2) = (0.0, 0.0);
        for _ in 0..N {
            let x = s.next_exponential();
            assert!(x >= 0.0);
            e1 += x;
            e2 += x * x;
        }
        assert!((e1 / n - 1.0).abs() < 0.01, "exp mean {}", e1 / n);
        assert!((e2 / n - 2.0).abs() < 0.05, "exp 2nd moment {}", e2 / n);
    }

    #[test]
    fn next_below_is_unbiased_and_replayable() {
        let mut s = KEY.stream();
        let mut counts = [0u32; 7];
        for _ in 0..70_000 {
            counts[s.next_below(7) as usize] += 1;
        }
        for (i, &c) in counts.iter().enumerate() {
            assert!(
                (9_500..=10_500).contains(&c),
                "biased bucket {i}: {c} (expect ~10000)"
            );
        }
        // Replay: same key + index range → same values, even through
        // rejection loops (the consumed-count is content-determined).
        let checkpoint = StreamCheckpoint::current(KEY, 12345);
        let mut a = Stream::resume(checkpoint).expect("current checkpoint is supported");
        let mut b = Stream::resume(checkpoint).expect("current checkpoint is supported");
        for _ in 0..1000 {
            assert_eq!(a.next_below(1000), b.next_below(1000));
        }
        assert_eq!(
            a.index(),
            b.index(),
            "rejection consumption must be deterministic"
        );
    }

    #[test]
    fn checkpoint_resume_equality() {
        let mut s = KEY.stream();
        for _ in 0..100 {
            let _ = s.next_normal();
        }
        let retained = s.checkpoint().to_canonical_le_bytes();
        let tail_a: Vec<u64> = (0..50).map(|_| s.next_u64()).collect();
        let mut resumed =
            Stream::resume_retained(&retained).expect("current retained checkpoint is supported");
        let tail_b: Vec<u64> = (0..50).map(|_| resumed.next_u64()).collect();
        assert_eq!(
            tail_a, tail_b,
            "resume from checkpoint must continue identically"
        );
    }

    #[test]
    fn checkpoint_canonical_le_known_answer_and_length_refusal() {
        let bytes = CHECKPOINT_KAT.to_canonical_le_bytes();
        assert_eq!(bytes, CHECKPOINT_CANONICAL_LE_KAT);
        assert_eq!(
            StreamCheckpoint::from_canonical_le_bytes(&bytes),
            Ok(CHECKPOINT_KAT)
        );

        let truncated = &bytes[..STREAM_CHECKPOINT_CANONICAL_LEN - 1];
        assert_eq!(
            StreamCheckpoint::from_canonical_le_bytes(truncated),
            Err(StreamReplayError::InvalidCheckpointLength {
                actual: STREAM_CHECKPOINT_CANONICAL_LEN - 1,
                expected: STREAM_CHECKPOINT_CANONICAL_LEN,
            })
        );

        let mut trailing = [0; STREAM_CHECKPOINT_CANONICAL_LEN + 1];
        trailing[..STREAM_CHECKPOINT_CANONICAL_LEN].copy_from_slice(&bytes);
        assert_eq!(
            StreamCheckpoint::from_canonical_le_bytes(&trailing),
            Err(StreamReplayError::InvalidCheckpointLength {
                actual: STREAM_CHECKPOINT_CANONICAL_LEN + 1,
                expected: STREAM_CHECKPOINT_CANONICAL_LEN,
            })
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one auditable field-to-offset matrix for the complete frame
    fn checkpoint_transport_mutation_battery() {
        const MAGIC_END: usize = STREAM_CHECKPOINT_MAGIC.len();
        const DOMAIN_END: usize = MAGIC_END + STREAM_CHECKPOINT_IDENTITY_DOMAIN.len();
        const CHECKPOINT_VERSION_END: usize = DOMAIN_END + 4;
        const STREAM_VERSION_END: usize = CHECKPOINT_VERSION_END + 4;
        const SEED_END: usize = STREAM_VERSION_END + 8;
        const KERNEL_END: usize = SEED_END + 4;
        const TILE_END: usize = KERNEL_END + 4;

        let base = CHECKPOINT_KAT.to_canonical_le_bytes();

        let mut invalid_magic = base;
        invalid_magic[0] ^= 1;
        assert_eq!(
            StreamCheckpoint::from_canonical_le_bytes(&invalid_magic),
            Err(StreamReplayError::InvalidCheckpointMagic)
        );

        let mut invalid_domain = base;
        invalid_domain[MAGIC_END] ^= 1;
        assert_eq!(
            StreamCheckpoint::from_canonical_le_bytes(&invalid_domain),
            Err(StreamReplayError::InvalidCheckpointDomain)
        );

        let mutations = [
            (
                "checkpoint-version",
                StreamCheckpoint {
                    checkpoint_version: STREAM_CHECKPOINT_VERSION + 1,
                    ..CHECKPOINT_KAT
                }
                .to_canonical_le_bytes(),
                DOMAIN_END,
                CHECKPOINT_VERSION_END,
            ),
            (
                "stream-semantics-version",
                StreamCheckpoint {
                    stream_semantics_version: STREAM_SEMANTICS_VERSION + 1,
                    ..CHECKPOINT_KAT
                }
                .to_canonical_le_bytes(),
                CHECKPOINT_VERSION_END,
                STREAM_VERSION_END,
            ),
            (
                "seed",
                StreamCheckpoint {
                    key: StreamKey {
                        seed: CHECKPOINT_KAT.key.seed ^ 1,
                        ..CHECKPOINT_KAT.key
                    },
                    ..CHECKPOINT_KAT
                }
                .to_canonical_le_bytes(),
                STREAM_VERSION_END,
                SEED_END,
            ),
            (
                "kernel",
                StreamCheckpoint {
                    key: StreamKey {
                        kernel: CHECKPOINT_KAT.key.kernel ^ 1,
                        ..CHECKPOINT_KAT.key
                    },
                    ..CHECKPOINT_KAT
                }
                .to_canonical_le_bytes(),
                SEED_END,
                KERNEL_END,
            ),
            (
                "tile",
                StreamCheckpoint {
                    key: StreamKey {
                        tile: CHECKPOINT_KAT.key.tile ^ 1,
                        ..CHECKPOINT_KAT.key
                    },
                    ..CHECKPOINT_KAT
                }
                .to_canonical_le_bytes(),
                KERNEL_END,
                TILE_END,
            ),
            (
                "index",
                StreamCheckpoint {
                    index: CHECKPOINT_KAT.index ^ 1,
                    ..CHECKPOINT_KAT
                }
                .to_canonical_le_bytes(),
                TILE_END,
                STREAM_CHECKPOINT_CANONICAL_LEN,
            ),
        ];
        for (field, changed, start, end) in mutations {
            assert_only_transport_range_changed(&base, &changed, start, end, field);
        }

        for checkpoint in [
            StreamCheckpoint {
                key: StreamKey {
                    seed: CHECKPOINT_KAT.key.seed ^ 1,
                    ..CHECKPOINT_KAT.key
                },
                ..CHECKPOINT_KAT
            },
            StreamCheckpoint {
                key: StreamKey {
                    kernel: CHECKPOINT_KAT.key.kernel ^ 1,
                    ..CHECKPOINT_KAT.key
                },
                ..CHECKPOINT_KAT
            },
            StreamCheckpoint {
                key: StreamKey {
                    tile: CHECKPOINT_KAT.key.tile ^ 1,
                    ..CHECKPOINT_KAT.key
                },
                ..CHECKPOINT_KAT
            },
            StreamCheckpoint {
                index: CHECKPOINT_KAT.index ^ 1,
                ..CHECKPOINT_KAT
            },
        ] {
            assert_eq!(
                StreamCheckpoint::from_canonical_le_bytes(&checkpoint.to_canonical_le_bytes()),
                Ok(checkpoint),
                "current semantic payload must round-trip exactly"
            );
        }
    }

    #[test]
    fn stream_identity_mutation_battery() {
        let key = StreamKey {
            seed: 0x0123_4567_89ab_cdef,
            kernel: 0x1357_9bdf,
            tile: 0x2468_ace0,
        };
        let index = 0x1020_3040_5060_7080;
        let base = Stream::at(key, index);
        let mutations = [
            (
                "seed-low",
                Stream::at(
                    StreamKey {
                        seed: key.seed ^ 1,
                        ..key
                    },
                    index,
                ),
            ),
            (
                "seed-high",
                Stream::at(
                    StreamKey {
                        seed: key.seed ^ (1 << 32),
                        ..key
                    },
                    index,
                ),
            ),
            (
                "kernel",
                Stream::at(
                    StreamKey {
                        kernel: key.kernel ^ 1,
                        ..key
                    },
                    index,
                ),
            ),
            (
                "tile",
                Stream::at(
                    StreamKey {
                        tile: key.tile ^ 1,
                        ..key
                    },
                    index,
                ),
            ),
            ("index-low", Stream::at(key, index ^ 1)),
            ("index-high", Stream::at(key, index ^ (1 << 32))),
        ];
        for (field, mutated) in mutations {
            assert_ne!(
                base, mutated,
                "semantic field {field} did not move replay bits"
            );
        }

        let checkpoint = StreamCheckpoint::current(key, index);
        let mut first = Stream::resume(checkpoint).expect("current checkpoint is supported");
        let mut replay = Stream::resume(checkpoint).expect("current checkpoint is supported");
        assert_eq!(first.next_u64(), replay.next_u64());
        assert_eq!(first.next_f64().to_bits(), replay.next_f64().to_bits());
        assert_eq!(first.checkpoint(), replay.checkpoint());
    }

    #[test]
    fn stale_checkpoint_versions_fail_closed() {
        let current = StreamCheckpoint::current(KEY, 17);
        assert!(Stream::resume(current).is_ok());

        for declared in [0, STREAM_CHECKPOINT_VERSION + 1] {
            let unsupported = StreamCheckpoint {
                checkpoint_version: declared,
                ..current
            };
            let expected = StreamReplayError::UnknownCheckpointVersion {
                declared,
                supported: STREAM_CHECKPOINT_VERSION,
            };
            assert_eq!(
                Stream::resume(unsupported).expect_err("unsupported checkpoint version"),
                expected
            );
            assert_eq!(
                Stream::resume_retained(&unsupported.to_canonical_le_bytes())
                    .expect_err("unsupported retained checkpoint version"),
                expected
            );
        }

        for declared in [0, STREAM_SEMANTICS_VERSION + 1] {
            let unsupported = StreamCheckpoint {
                stream_semantics_version: declared,
                ..current
            };
            let expected = StreamReplayError::UnknownStreamSemanticsVersion {
                declared,
                supported: STREAM_SEMANTICS_VERSION,
            };
            assert_eq!(
                Stream::resume(unsupported).expect_err("unsupported stream semantics version"),
                expected
            );
            assert_eq!(
                Stream::resume_retained(&unsupported.to_canonical_le_bytes())
                    .expect_err("unsupported retained stream semantics version"),
                expected
            );
        }
    }

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }
}
