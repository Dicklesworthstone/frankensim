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
pub mod governor;
pub mod guidance;
pub mod token;

pub use estimate::{CalibrationReport, Estimate, estimate};
pub use governor::{
    Charge, DegradationEvent, DegradationStep, Enforcement, Governor, SubmitOutcome,
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
}

impl fmt::Display for SessionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SessionError::UnknownSession { id } => write!(f, "unknown session {id}"),
            SessionError::Submission { what } => write!(f, "submission failed: {what}"),
            SessionError::Persistence { what } => write!(f, "persistence failed: {what}"),
        }
    }
}

impl std::error::Error for SessionError {}
