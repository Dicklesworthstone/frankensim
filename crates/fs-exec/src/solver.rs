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
            let len = self.get_u64()? as usize;
            let remaining = self.bytes.len() - self.at;
            if remaining < len.saturating_mul(8) {
                return Err(CodecError {
                    at: self.at,
                    what: "f64 slice body",
                    needed: len.saturating_mul(8),
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

/// A serializable solver snapshot. Implementations must be self-contained
/// (no pointers; artifact references by content hash) — see module docs.
pub trait SolverState: Sized {
    /// Write the snapshot.
    fn encode(&self, enc: &mut codec::Enc);

    /// Read a snapshot.
    ///
    /// # Errors
    /// [`codec::CodecError`] on truncated/incompatible bytes.
    fn decode(dec: &mut codec::Dec<'_>) -> Result<Self, codec::CodecError>;

    /// The snapshot bytes (ledger checkpoint payload).
    fn to_bytes(&self) -> Vec<u8> {
        let mut enc = codec::Enc::new();
        self.encode(&mut enc);
        enc.into_bytes()
    }

    /// Rebuild from snapshot bytes, rejecting trailing garbage.
    ///
    /// # Errors
    /// [`codec::CodecError`] on truncation or trailing bytes.
    fn from_bytes(bytes: &[u8]) -> Result<Self, codec::CodecError> {
        let mut dec = codec::Dec::new(bytes);
        let state = Self::decode(&mut dec)?;
        if dec.is_empty() {
            Ok(state)
        } else {
            Err(codec::CodecError {
                at: bytes.len(),
                what: "end of snapshot",
                needed: 0,
                remaining: 1,
            })
        }
    }

    /// Deterministic content hash of the snapshot (FNV-1a until the
    /// BLAKE3-class ledger hash supersedes it — same upgrade path as
    /// fs-obs).
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
/// [`codec::CodecError`] when the state's encode/decode disagree — a
/// serialization bug surfaced early.
pub fn fork<S: SolverState>(state: &S) -> Result<S, codec::CodecError> {
    S::from_bytes(&state.to_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cx::{CancelGate, ExecMode, StreamKey};
    use asupersync::types::Budget;

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
            for i in 0..n {
                let left = if i > 0 { state.x[i - 1] } else { 0.0 };
                let right = if i + 1 < n { state.x[i + 1] } else { 0.0 };
                next[i] = state.x[i] + 0.6 * ((self.rhs[i] - left - right) / 4.0 - state.x[i]);
                residual = residual.max((next[i] - state.x[i]).abs());
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
        let (solver_b, _) = {
            let mut j = jacobi().0;
            j.rhs.iter_mut().for_each(|r| *r += 0.5);
            (j, ())
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
