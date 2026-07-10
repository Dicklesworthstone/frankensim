//! fs-session (plan §11.3): sessions, capability tokens, and the resource
//! GOVERNOR — budgets are ENFORCED, not advisory — plus the agent-proofing
//! trio: idempotency keys (a retry cannot double-spend), `estimate()` dry
//! runs (plan before you spend), and errors as GUIDANCE ("a refusal that
//! teaches is worth ten silent successes").
//!
//! Layer: L6 (HELM). Threading contract: the governor's hot paths are
//! `Send + Sync` (in-memory, mutex-guarded) so enforcement and idempotency
//! survive concurrent submission storms; ledger persistence is an explicit
//! single-threaded `flush_to_ledger` step because fsqlite connections are
//! `!Send` by design.

pub mod estimate;
pub mod gemm_tune;
pub mod governor;
pub mod guidance;
pub mod token;

pub use estimate::{CalibrationReport, Estimate, estimate};
pub use gemm_tune::{
    GemmDispatch, GemmTuneError, gemm_f64_session, gemm_kernel_key, gemm_shape_class,
};
pub use governor::{
    Charge, DegradationEvent, DegradationStep, Enforcement, Governor, StepPhase, SubmitOutcome,
};
pub use guidance::Guidance;
pub use token::{CapabilityToken, SessionId};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured session failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum SessionError {
    /// The session id is unknown to the governor.
    UnknownSession {
        /// The id.
        id: u64,
    },
    /// A resource grant, charge, or accumulated meter is outside its valid
    /// finite, non-negative domain.
    InvalidResource {
        /// The resource field.
        resource: &'static str,
        /// The rejected value.
        value: f64,
        /// The required domain.
        requirement: &'static str,
    },
    /// A submission failed structurally (parse/admission).
    Submission {
        /// Diagnosis.
        what: String,
    },
    /// Ledger persistence failed.
    Persistence {
        /// Diagnosis.
        what: String,
    },
    /// A memory-pressure level outside the declared ladder 1..=3
    /// (bead gp3.13: out-of-ladder levels are refused, never clamped).
    InvalidPressureLevel {
        /// The rejected level.
        level: u8,
    },
    /// Level-3 pressure targeted a session opened without a bound
    /// cancellation gate — a pause that cannot reach the computation
    /// is refused, not ledgered (bead gp3.13).
    UngatedSession {
        /// The id.
        id: u64,
    },
    /// A pause acknowledgement arrived with no outstanding pause
    /// request for the session (bead gp3.13).
    NoPendingPause {
        /// The id.
        id: u64,
    },
}

impl fmt::Display for SessionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SessionError::UnknownSession { id } => write!(f, "unknown session {id}"),
            SessionError::InvalidResource {
                resource,
                value,
                requirement,
            } => write!(
                f,
                "invalid {resource} value {value}: {requirement}; the governor did not mutate \
                 session state"
            ),
            SessionError::Submission { what } => write!(f, "submission failed: {what}"),
            SessionError::Persistence { what } => write!(f, "persistence failed: {what}"),
            SessionError::InvalidPressureLevel { level } => write!(
                f,
                "memory-pressure level {level} is outside the declared ladder 1..=3; \
                 out-of-ladder levels are refused, never clamped"
            ),
            SessionError::UngatedSession { id } => write!(
                f,
                "session {id} was opened without a cancellation gate; level-3 pressure \
                 (pause-serialize-resume) is refused — open with open_session_gated to \
                 bind the session's own gate"
            ),
            SessionError::NoPendingPause { id } => write!(
                f,
                "session {id} has no outstanding pause request to acknowledge"
            ),
        }
    }
}

impl std::error::Error for SessionError {}
