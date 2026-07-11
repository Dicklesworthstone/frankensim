//! Resumable solvers (plan §5.2 behavior 2): iterative solvers as EXPLICIT
//! state machines whose snapshots serialize, migrate, resume, and FORK —
//! the forkable-worlds enabler and the resource governor's
//! pause-serialize-resume primitive ("pause the LES, run the urgent trim
//! study, resume" must be routine, not heroic).
//!
//! Distribution-readiness (plan §17): the serialized representation is
//! self-contained bytes — no pointers, no shared-memory assumptions, large
//! artifacts referenced by content hash — so "migrate" can someday mean
//! "to another machine" without an API change.
//!
//! Determinism invariant (G4): pause → serialize → deserialize → resume
//! reproduces the uninterrupted trajectory BIT-EXACTLY. The conformance
//! suite asserts it on a reference iterative solver.

use crate::cx::Cx;

/// In-house, deterministic, little-endian state codec (P1: no serde).
/// Floats travel as raw bits (`to_bits`), so round-trips are bit-exact
/// including NaN payloads and signed zeros.
pub mod codec {
    use core::fmt;

    /// Structured decode failure (Decalogue P10).
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct CodecError {
        /// Byte offset where decoding failed.
        pub at: usize,
        /// What the decoder was reading.
        pub what: &'static str,
        /// Bytes it needed.
        pub needed: usize,
        /// Bytes that remained.
        pub remaining: usize,
    }

    impl fmt::Display for CodecError {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(
                f,
                "solver-state decode failed at byte {}: reading {} needs {} bytes but {} \
                 remain; the snapshot is truncated or from an incompatible encoder version",
                self.at, self.what, self.needed, self.remaining
            )
        }
    }

    impl core::error::Error for CodecError {}

    /// Append-only encoder.
    #[derive(Debug, Default)]
    pub struct Enc {
        buf: Vec<u8>,
    }

    impl Enc {
        /// Fresh encoder.
        #[must_use]
        pub fn new() -> Self {
            Enc::default()
        }

        /// Append a u32 (little-endian).
        pub fn put_u32(&mut self, v: u32) {
            self.buf.extend_from_slice(&v.to_le_bytes());
        }

        /// Append a u64 (little-endian).
        pub fn put_u64(&mut self, v: u64) {
            self.buf.extend_from_slice(&v.to_le_bytes());
        }

        /// Append an f64 as raw bits (bit-exact round-trip).
        pub fn put_f64(&mut self, v: f64) {
            self.put_u64(v.to_bits());
        }

        /// Append a length-prefixed f64 slice.
        pub fn put_f64_slice(&mut self, xs: &[f64]) {
            self.put_u64(xs.len() as u64);
            for &x in xs {
                self.put_f64(x);
            }
        }

        /// Finish, yielding the snapshot bytes.
        #[must_use]
        pub fn into_bytes(self) -> Vec<u8> {
            self.buf
        }
    }

    /// Cursor decoder over snapshot bytes.
    #[derive(Debug)]
    pub struct Dec<'a> {
        bytes: &'a [u8],
        at: usize,
    }

    impl<'a> Dec<'a> {
        /// Decode from `bytes`.
        #[must_use]
        pub fn new(bytes: &'a [u8]) -> Self {
            Dec { bytes, at: 0 }
        }

        fn take(&mut self, n: usize, what: &'static str) -> Result<&'a [u8], CodecError> {
            let remaining = self.bytes.len() - self.at;
            if remaining < n {
                return Err(CodecError {
                    at: self.at,
                    what,
                    needed: n,
                    remaining,
                });
            }
            let s = &self.bytes[self.at..self.at + n];
            self.at += n;
            Ok(s)
        }

        /// Read a u32.
        ///
        /// # Errors
        /// [`CodecError`] on truncation.
        pub fn get_u32(&mut self) -> Result<u32, CodecError> {
            Ok(u32::from_le_bytes(
                self.take(4, "u32")?.try_into().expect("length checked"),
            ))
        }

        /// Read a u64.
        ///
        /// # Errors
        /// [`CodecError`] on truncation.
        pub fn get_u64(&mut self) -> Result<u64, CodecError> {
            Ok(u64::from_le_bytes(
                self.take(8, "u64")?.try_into().expect("length checked"),
            ))
        }

        /// Read an f64 (from raw bits).
        ///
        /// # Errors
        /// [`CodecError`] on truncation.
        pub fn get_f64(&mut self) -> Result<f64, CodecError> {
            Ok(f64::from_bits(self.get_u64()?))
        }

        /// Read a length-prefixed f64 slice.
        ///
        /// # Errors
        /// [`CodecError`] on truncation (including an implausible length).
        pub fn get_f64_vec(&mut self) -> Result<Vec<f64>, CodecError> {
            let encoded_len = self.get_u64()?;
            let remaining = self.bytes.len() - self.at;
            let len = usize::try_from(encoded_len).map_err(|_| CodecError {
                at: self.at,
                what: "f64 slice length exceeds platform usize",
                needed: usize::MAX,
                remaining,
            })?;
            let needed = len.checked_mul(8).ok_or(CodecError {
                at: self.at,
                what: "f64 slice byte length overflow",
                needed: usize::MAX,
                remaining,
            })?;
            if remaining < needed {
                return Err(CodecError {
                    at: self.at,
                    what: "f64 slice body",
                    needed,
                    remaining,
                });
            }
            (0..len).map(|_| self.get_f64()).collect()
        }

        /// True when every byte was consumed (decoders should check this to
        /// reject trailing garbage).
        #[must_use]
        pub fn is_empty(&self) -> bool {
            self.at == self.bytes.len()
        }
    }
}

/// The snapshot ENVELOPE (bead wf9.8.2): magic, versions, type
/// identity, length, checksum, and provenance — all validated BEFORE
/// the payload decoder runs, so same-length bytes from another solver,
/// another schema version, a bit flip, a truncation, or an append can
/// never decode into plausible-but-wrong state.
pub mod envelope {
    use core::fmt;

    /// Envelope magic (8 bytes).
    pub const MAGIC: [u8; 8] = *b"FSEXSNAP";
    /// Envelope layout version. Bump only with a recorded migration.
    pub const ENVELOPE_VERSION: u32 = 1;
    /// Header size: magic + env version + type id + schema version +
    /// provenance + payload len + payload hash.
    pub const HEADER_LEN: usize = 8 + 4 + 8 + 4 + 8 + 8 + 8;

    /// Structured envelope refusal — never a wrong-state decode.
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub enum EnvelopeError {
        /// Not a snapshot envelope at all.
        BadMagic,
        /// Shorter than a header (or than the declared payload).
        Truncated {
            /// Bytes needed.
            needed: usize,
            /// Bytes present.
            have: usize,
        },
        /// Envelope layout from a different (unsupported) version.
        UnknownEnvelopeVersion {
            /// The version found.
            found: u32,
        },
        /// The snapshot belongs to a DIFFERENT state type.
        WrongTypeId {
            /// The expected stable type id.
            expected: u64,
            /// The id in the envelope.
            found: u64,
        },
        /// Same type, incompatible schema version: explicit refusal
        /// (the structured alternative to a silent wrong decode; write
        /// a migration when a version must remain readable).
        IncompatibleSchema {
            /// The reader's schema version.
            expected: u32,
            /// The snapshot's schema version.
            found: u32,
        },
        /// Declared payload length disagrees with the actual bytes
        /// (truncation past the header, or appended bytes).
        LengthMismatch {
            /// Length declared in the header.
            declared: u64,
            /// Bytes actually present after the header.
            actual: u64,
        },
        /// Payload bytes do not hash to the declared checksum.
        ChecksumMismatch {
            /// The declared hash.
            declared: u64,
            /// The computed hash.
            computed: u64,
        },
    }

    impl fmt::Display for EnvelopeError {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                EnvelopeError::BadMagic => write!(f, "not a solver snapshot (bad magic)"),
                EnvelopeError::Truncated { needed, have } => write!(
                    f,
                    "snapshot truncated: needs {needed} bytes, {have} present"
                ),
                EnvelopeError::UnknownEnvelopeVersion { found } => write!(
                    f,
                    "unknown snapshot envelope version {found} (this reader supports {ENVELOPE_VERSION})"
                ),
                EnvelopeError::WrongTypeId { expected, found } => write!(
                    f,
                    "snapshot is for state type {found:#018x}, not {expected:#018x} — refusing a cross-type decode"
                ),
                EnvelopeError::IncompatibleSchema { expected, found } => write!(
                    f,
                    "snapshot schema v{found} is incompatible with this reader (v{expected}); \
                     write an explicit migration or regenerate the snapshot"
                ),
                EnvelopeError::LengthMismatch { declared, actual } => write!(
                    f,
                    "snapshot payload length mismatch: header declares {declared}, {actual} bytes present"
                ),
                EnvelopeError::ChecksumMismatch { declared, computed } => write!(
                    f,
                    "snapshot payload checksum mismatch (declared {declared:#018x}, computed {computed:#018x}): corrupted bytes"
                ),
            }
        }
    }

    impl core::error::Error for EnvelopeError {}

    /// Seal a payload: canonical header + payload bytes.
    #[must_use]
    pub fn seal(type_id: u64, schema_version: u32, provenance: u64, payload: &[u8]) -> Vec<u8> {
        let mut out = Vec::with_capacity(HEADER_LEN + payload.len());
        out.extend_from_slice(&MAGIC);
        out.extend_from_slice(&ENVELOPE_VERSION.to_le_bytes());
        out.extend_from_slice(&type_id.to_le_bytes());
        out.extend_from_slice(&schema_version.to_le_bytes());
        out.extend_from_slice(&provenance.to_le_bytes());
        out.extend_from_slice(&(payload.len() as u64).to_le_bytes());
        out.extend_from_slice(&fs_obs::fnv1a64(payload).to_le_bytes());
        out.extend_from_slice(payload);
        out
    }

    /// Validate an envelope for (`type_id`, `schema_version`) and return
    /// `(payload, provenance)`. Every header field is checked before a
    /// single payload byte is interpreted.
    ///
    /// # Errors
    /// [`EnvelopeError`], each naming the exact refusal.
    pub fn open(
        bytes: &[u8],
        type_id: u64,
        schema_version: u32,
    ) -> Result<(&[u8], u64), EnvelopeError> {
        if bytes.len() < HEADER_LEN {
            if bytes.len() >= 8 && bytes[..8] != MAGIC {
                return Err(EnvelopeError::BadMagic);
            }
            return Err(EnvelopeError::Truncated {
                needed: HEADER_LEN,
                have: bytes.len(),
            });
        }
        if bytes[..8] != MAGIC {
            return Err(EnvelopeError::BadMagic);
        }
        let u32_at = |o: usize| u32::from_le_bytes(bytes[o..o + 4].try_into().expect("len"));
        let u64_at = |o: usize| u64::from_le_bytes(bytes[o..o + 8].try_into().expect("len"));
        let env_version = u32_at(8);
        if env_version != ENVELOPE_VERSION {
            return Err(EnvelopeError::UnknownEnvelopeVersion { found: env_version });
        }
        let found_type = u64_at(12);
        if found_type != type_id {
            return Err(EnvelopeError::WrongTypeId {
                expected: type_id,
                found: found_type,
            });
        }
        let found_schema = u32_at(20);
        if found_schema != schema_version {
            return Err(EnvelopeError::IncompatibleSchema {
                expected: schema_version,
                found: found_schema,
            });
        }
        let provenance = u64_at(24);
        let declared_len = u64_at(32);
        let payload = &bytes[HEADER_LEN..];
        if declared_len != payload.len() as u64 {
            return Err(EnvelopeError::LengthMismatch {
                declared: declared_len,
                actual: payload.len() as u64,
            });
        }
        let declared_hash = u64_at(40);
        let computed = fs_obs::fnv1a64(payload);
        if computed != declared_hash {
            return Err(EnvelopeError::ChecksumMismatch {
                declared: declared_hash,
                computed,
            });
        }
        Ok((payload, provenance))
    }
}

/// A snapshot failure: envelope refusal or payload decode error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SnapshotError {
    /// The envelope refused (wrong type/version/corruption) — the
    /// payload decoder never ran.
    Envelope(envelope::EnvelopeError),
    /// The envelope validated but the payload decoder failed (an
    /// encode/decode bug within one schema version).
    Payload(codec::CodecError),
}

impl core::fmt::Display for SnapshotError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SnapshotError::Envelope(e) => write!(f, "{e}"),
            SnapshotError::Payload(e) => write!(f, "{e}"),
        }
    }
}

impl core::error::Error for SnapshotError {}

/// A serializable solver snapshot. Implementations must be self-contained
/// (no pointers; artifact references by content hash) — see module docs.
/// Every state declares a STABLE type id and a schema version (bead
/// wf9.8.2): snapshots travel inside the [`envelope`], which is validated
/// before the payload decoder ever runs.
pub trait SolverState: Sized {
    /// Stable type identity. Never reuse across state types; never
    /// change for a type (that is what [`Self::SCHEMA_VERSION`] is for).
    const TYPE_ID: u64;
    /// Payload schema version. Bump on ANY layout change; readers
    /// refuse other versions structurally.
    const SCHEMA_VERSION: u32;

    /// Write the snapshot payload.
    fn encode(&self, enc: &mut codec::Enc);

    /// Read a snapshot payload.
    ///
    /// # Errors
    /// [`codec::CodecError`] on truncated/incompatible bytes.
    fn decode(dec: &mut codec::Dec<'_>) -> Result<Self, codec::CodecError>;

    /// Seal the snapshot with an explicit caller-ledgered provenance
    /// (run/ledger identity — e.g. a `RunId` value or ledger row id).
    fn seal(&self, provenance: u64) -> Vec<u8> {
        let mut enc = codec::Enc::new();
        self.encode(&mut enc);
        envelope::seal(
            Self::TYPE_ID,
            Self::SCHEMA_VERSION,
            provenance,
            &enc.into_bytes(),
        )
    }

    /// Validate the envelope, decode the payload, and return the state
    /// with its provenance. Rejects trailing payload garbage.
    ///
    /// # Errors
    /// [`SnapshotError`] — envelope refusals never reach the decoder.
    fn unseal(bytes: &[u8]) -> Result<(Self, u64), SnapshotError> {
        let (payload, provenance) = envelope::open(bytes, Self::TYPE_ID, Self::SCHEMA_VERSION)
            .map_err(SnapshotError::Envelope)?;
        let mut dec = codec::Dec::new(payload);
        let state = Self::decode(&mut dec).map_err(SnapshotError::Payload)?;
        if dec.is_empty() {
            Ok((state, provenance))
        } else {
            Err(SnapshotError::Payload(codec::CodecError {
                at: payload.len(),
                what: "end of snapshot payload",
                needed: 0,
                remaining: 1,
            }))
        }
    }

    /// The ENVELOPED snapshot bytes (ledger checkpoint payload) with
    /// unattributed provenance; ledger paths should prefer
    /// [`SolverState::seal`] with a real run identity.
    fn to_bytes(&self) -> Vec<u8> {
        self.seal(0)
    }

    /// Rebuild from enveloped snapshot bytes.
    ///
    /// # Errors
    /// [`SnapshotError`] on any envelope refusal or payload error.
    fn from_bytes(bytes: &[u8]) -> Result<Self, SnapshotError> {
        Self::unseal(bytes).map(|(state, _)| state)
    }

    /// Deterministic content hash of the ENVELOPED snapshot (FNV-1a
    /// until the BLAKE3-class ledger hash supersedes it — same upgrade
    /// path as fs-obs).
    fn content_hash(&self) -> u64 {
        fs_obs::fnv1a64(&self.to_bytes())
    }
}

/// One bounded step's verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepVerdict<T> {
    /// More steps remain.
    Continue,
    /// Converged/finished with a result.
    Done(T),
}

/// An iterative solver as an explicit state machine: `step` advances one
/// BOUNDED unit of work (an iteration, a sweep) — the pause granularity.
pub trait ResumableSolver {
    /// The serializable snapshot type.
    type State: SolverState;
    /// The final result type.
    type Out;

    /// Advance one bounded step. Implementations may poll `cx` internally
    /// for finer-grained cancellation inside expensive steps.
    fn step(&self, state: &mut Self::State, cx: &Cx<'_>) -> StepVerdict<Self::Out>;
}

/// The outcome of [`drive`]: finished, or paused holding the resumable
/// snapshot (the caller serializes it to the ledger and later resumes or
/// forks).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverProgress<S, T> {
    /// Ran to completion.
    Done(T),
    /// Cancellation/pause was requested; `state` resumes bit-exactly.
    Paused(S),
}

/// Drive a solver until completion or until the context's cancel gate is
/// requested — pause IS the cancellation path, which is what makes
/// "pause, run something urgent, resume" routine (graceful-degradation
/// hook for the session governor).
pub fn drive<R: ResumableSolver>(
    solver: &R,
    mut state: R::State,
    cx: &Cx<'_>,
) -> SolverProgress<R::State, R::Out> {
    loop {
        if cx.is_cancel_requested() {
            return SolverProgress::Paused(state);
        }
        match solver.step(&mut state, cx) {
            StepVerdict::Continue => {}
            StepVerdict::Done(out) => return SolverProgress::Done(out),
        }
    }
}

/// Fork a solver state by round-tripping it through its serialized form —
/// proving at fork time that the snapshot really is self-contained (a fork
/// that only works in-memory is a distribution bug waiting to happen).
///
/// # Errors
/// [`SnapshotError`] when the state's encode/decode disagree — a
/// serialization bug surfaced early.
pub fn fork<S: SolverState>(state: &S) -> Result<S, SnapshotError> {
    S::from_bytes(&state.to_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cx::{CancelGate, ExecMode, StreamKey};
    use asupersync::types::Budget;

    /// wf9.8.2 acceptance: the envelope refuses every corruption and
    /// misbinding class BEFORE the payload decoder runs.
    #[test]
    fn envelope_refuses_every_misbinding_class() {
        // A twin state with the IDENTICAL payload layout but its own
        // type id: same-length bytes must not cross-decode.
        #[derive(Debug, PartialEq)]
        struct TwinState {
            x: Vec<f64>,
            iter: u64,
        }
        impl SolverState for TwinState {
            const TYPE_ID: u64 = 0x5457_494e_0000_0001;
            const SCHEMA_VERSION: u32 = 1;
            fn encode(&self, enc: &mut codec::Enc) {
                enc.put_u64(self.iter);
                enc.put_f64_slice(&self.x);
            }
            fn decode(dec: &mut codec::Dec<'_>) -> Result<Self, codec::CodecError> {
                Ok(TwinState {
                    iter: dec.get_u64()?,
                    x: dec.get_f64_vec()?,
                })
            }
        }
        let state = JacobiState {
            x: vec![1.5, -2.25, 0.0],
            iter: 42,
        };
        let sealed = state.seal(0xABCD);
        // The happy path round-trips bit-exactly WITH provenance.
        let (back, prov) = JacobiState::unseal(&sealed).expect("valid seal");
        assert_eq!(back, state);
        assert_eq!(prov, 0xABCD);
        // Cross-type: identical payload layout, refused by TYPE ID.
        assert!(matches!(
            TwinState::unseal(&sealed),
            Err(SnapshotError::Envelope(
                envelope::EnvelopeError::WrongTypeId { .. }
            ))
        ));
        // Bit flip in the payload: checksum refuses.
        let mut flipped = sealed.clone();
        let last = flipped.len() - 1;
        flipped[last] ^= 0x40;
        assert!(matches!(
            JacobiState::unseal(&flipped),
            Err(SnapshotError::Envelope(
                envelope::EnvelopeError::ChecksumMismatch { .. }
            ))
        ));
        // Bit flip in the magic: not a snapshot.
        let mut bad_magic = sealed.clone();
        bad_magic[0] ^= 0x01;
        assert!(matches!(
            JacobiState::unseal(&bad_magic),
            Err(SnapshotError::Envelope(envelope::EnvelopeError::BadMagic))
        ));
        // Truncation: header-level and payload-level both refuse.
        assert!(matches!(
            JacobiState::unseal(&sealed[..10]),
            Err(SnapshotError::Envelope(
                envelope::EnvelopeError::Truncated { .. }
            ))
        ));
        assert!(matches!(
            JacobiState::unseal(&sealed[..sealed.len() - 3]),
            Err(SnapshotError::Envelope(
                envelope::EnvelopeError::LengthMismatch { .. }
            ))
        ));
        // Appended bytes: refused by the declared length.
        let mut appended = sealed.clone();
        appended.extend_from_slice(&[0u8; 5]);
        assert!(matches!(
            JacobiState::unseal(&appended),
            Err(SnapshotError::Envelope(
                envelope::EnvelopeError::LengthMismatch { .. }
            ))
        ));
        // Unknown envelope version.
        let mut future = sealed.clone();
        future[8..12].copy_from_slice(&9u32.to_le_bytes());
        assert!(matches!(
            JacobiState::unseal(&future),
            Err(SnapshotError::Envelope(
                envelope::EnvelopeError::UnknownEnvelopeVersion { found: 9 }
            ))
        ));
        // Stale schema version: structured refusal, not a wrong decode.
        let mut stale = sealed;
        stale[20..24].copy_from_slice(&7u32.to_le_bytes());
        assert!(matches!(
            JacobiState::unseal(&stale),
            Err(SnapshotError::Envelope(
                envelope::EnvelopeError::IncompatibleSchema {
                    expected: 1,
                    found: 7
                }
            ))
        ));
    }

    /// Reference solver: damped Jacobi on a fixed diagonally-dominant
    /// system (deterministic, non-trivial float trajectory).
    struct Jacobi {
        rhs: Vec<f64>,
        tol: f64,
    }

    #[derive(Debug, Clone, PartialEq)]
    struct JacobiState {
        x: Vec<f64>,
        iter: u64,
    }

    impl SolverState for JacobiState {
        const TYPE_ID: u64 = 0x4a41_434f_4249_0001;
        const SCHEMA_VERSION: u32 = 1;

        fn encode(&self, enc: &mut codec::Enc) {
            enc.put_u64(self.iter);
            enc.put_f64_slice(&self.x);
        }

        fn decode(dec: &mut codec::Dec<'_>) -> Result<Self, codec::CodecError> {
            Ok(JacobiState {
                iter: dec.get_u64()?,
                x: dec.get_f64_vec()?,
            })
        }
    }

    impl ResumableSolver for Jacobi {
        type State = JacobiState;
        type Out = (Vec<f64>, u64);

        fn step(&self, state: &mut JacobiState, _cx: &Cx<'_>) -> StepVerdict<(Vec<f64>, u64)> {
            let n = state.x.len();
            let mut next = vec![0.0f64; n];
            let mut residual = 0.0f64;
            for (i, slot) in next.iter_mut().enumerate() {
                let left = if i > 0 { state.x[i - 1] } else { 0.0 };
                let right = if i + 1 < n { state.x[i + 1] } else { 0.0 };
                *slot = state.x[i] + 0.6 * ((self.rhs[i] - left - right) / 4.0 - state.x[i]);
                residual = residual.max((*slot - state.x[i]).abs());
            }
            state.x = next;
            state.iter += 1;
            if residual < self.tol {
                StepVerdict::Done((state.x.clone(), state.iter))
            } else {
                StepVerdict::Continue
            }
        }
    }

    fn jacobi() -> (Jacobi, JacobiState) {
        let rhs: Vec<f64> = (0..32).map(|i| 1.0 + 0.25 * f64::from(i % 5)).collect();
        (
            Jacobi { rhs, tol: 1e-12 },
            JacobiState {
                x: vec![0.0; 32],
                iter: 0,
            },
        )
    }

    fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    #[test]
    fn codec_round_trips_are_bit_exact_and_reject_garbage() {
        let mut enc = codec::Enc::new();
        enc.put_u64(42);
        enc.put_f64(f64::NAN);
        enc.put_f64(-0.0);
        enc.put_f64_slice(&[1.5, f64::INFINITY, f64::MIN_POSITIVE]);
        let bytes = enc.into_bytes();
        let mut dec = codec::Dec::new(&bytes);
        assert_eq!(dec.get_u64().expect("u64"), 42);
        assert_eq!(
            dec.get_f64().expect("nan").to_bits(),
            f64::NAN.to_bits(),
            "NaN payload preserved"
        );
        assert_eq!(
            dec.get_f64().expect("neg zero").to_bits(),
            (-0.0f64).to_bits()
        );
        let v = dec.get_f64_vec().expect("slice");
        assert_eq!(v.len(), 3);
        assert!(dec.is_empty());
        // Truncation is a structured, teaching error.
        let err = codec::Dec::new(&bytes[..5])
            .get_u64()
            .expect_err("truncated");
        assert!(err.to_string().contains("truncated"), "{err}");
        let impossible_len = u64::MAX.to_le_bytes();
        let err = codec::Dec::new(&impossible_len)
            .get_f64_vec()
            .expect_err("wire lengths must fit usize and their byte extent");
        assert!(err.what.starts_with("f64 slice"), "{err}");
        #[cfg(target_pointer_width = "32")]
        assert_eq!(
            err.what, "f64 slice length exceeds platform usize",
            "32-bit readers must not truncate the u64 wire length"
        );
        // Trailing garbage is rejected by from_bytes.
        let (_, s0) = jacobi();
        let mut noisy = s0.to_bytes();
        noisy.push(0xFF);
        assert!(JacobiState::from_bytes(&noisy).is_err());
    }

    #[test]
    fn pause_serialize_resume_is_bit_exact_versus_uninterrupted() {
        let (solver, s0) = jacobi();
        // Uninterrupted reference.
        let gate = CancelGate::new();
        let SolverProgress::Done((x_ref, iters_ref)) =
            with_cx(&gate, |cx| drive(&solver, s0.clone(), cx))
        else {
            panic!("uninterrupted run must finish");
        };
        // Interrupted every step: advance ONE bounded step, then pause,
        // serialize, deserialize, resume — the maximally hostile schedule.
        let mut state = s0;
        let mut resumes = 0u64;
        let finished = loop {
            let g2 = CancelGate::new();
            let (st, verdict) = with_cx(&g2, |cx| {
                let mut st = state.clone();
                let verdict = solver.step(&mut st, cx);
                (st, verdict)
            });
            match verdict {
                StepVerdict::Done(out) => break out,
                StepVerdict::Continue => {
                    let bytes = st.to_bytes();
                    state = JacobiState::from_bytes(&bytes).expect("round trip");
                    resumes += 1;
                }
            }
        };
        assert_eq!(finished.1, iters_ref, "same iteration count");
        assert!(resumes > 10, "the trajectory must actually be interrupted");
        let bits_ref: Vec<u64> = x_ref.iter().map(|v| v.to_bits()).collect();
        let bits_paused: Vec<u64> = finished.0.iter().map(|v| v.to_bits()).collect();
        assert_eq!(bits_ref, bits_paused, "bit-exact continuation (G4 law)");
    }

    #[test]
    fn drive_pauses_on_cancel_and_resumes_to_the_same_answer() {
        let (solver, s0) = jacobi();
        let gate = CancelGate::new();
        let SolverProgress::Done((x_ref, _)) = with_cx(&gate, |cx| drive(&solver, s0.clone(), cx))
        else {
            panic!("reference finishes");
        };
        // Cancel mid-flight: drive must return Paused with usable state.
        let paused_state = {
            let gate = CancelGate::new();
            gate.request();
            match with_cx(&gate, |cx| drive(&solver, s0, cx)) {
                SolverProgress::Paused(s) => s,
                SolverProgress::Done(_) => panic!("pre-requested gate must pause"),
            }
        };
        let gate = CancelGate::new();
        let SolverProgress::Done((x_resumed, _)) =
            with_cx(&gate, |cx| drive(&solver, paused_state, cx))
        else {
            panic!("resume finishes");
        };
        assert_eq!(
            x_ref.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            x_resumed.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn forks_are_independent_and_serialization_proven() {
        let (solver, s0) = jacobi();
        // Advance 10 steps.
        let gate = CancelGate::new();
        let mut warm = s0;
        with_cx(&gate, |cx| {
            for _ in 0..10 {
                let _ = solver.step(&mut warm, cx);
            }
        });
        let fork_a = fork(&warm).expect("fork proves serializability");
        let fork_b = fork(&warm).expect("second fork");
        assert_eq!(fork_a.content_hash(), fork_b.content_hash());
        // Diverge: different subsequent inputs (different rhs) per fork.
        let solver_b = {
            let mut j = jacobi().0;
            j.rhs.iter_mut().for_each(|r| *r += 0.5);
            j
        };
        let SolverProgress::Done((xa, _)) = with_cx(&gate, |cx| drive(&solver, fork_a, cx)) else {
            panic!("fork A finishes");
        };
        let SolverProgress::Done((xb, _)) = with_cx(&gate, |cx| drive(&solver_b, fork_b, cx))
        else {
            panic!("fork B finishes");
        };
        assert_ne!(
            xa[0].to_bits(),
            xb[0].to_bits(),
            "forks with different inputs stay independent"
        );
    }
}
