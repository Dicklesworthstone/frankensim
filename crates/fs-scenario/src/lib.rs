//! fs-scenario: the boundary-condition and load-case ALGEBRA (plan
//! patch Rev D). A Region answers "what is the domain?"; a PhysicsModel
//! answers "which equations?"; a `Scenario` answers "WHAT IS BEING DONE
//! TO IT?" — as a typed value with dimensional analysis (fs-qty),
//! provenance (seed + canonical IR), and validity checks. Boundary
//! conditions are where simulations quietly become invalid; this crate
//! makes that class of mistake a structured, fixable refusal instead.
//!
//! Layer: L3 (FLUX support). Runtime deps: `std`, fs-blake3, fs-qty,
//! fs-rand, fs-cheb, fs-exec, fs-ga, fs-math. The Design Ledger stores
//! scenarios as canonical-IR artifacts — that integration lives ABOVE this layer
//! (exercised here via a dev-dependency in conformance tests).

pub mod bc;
pub mod ensemble;
pub mod frame;
pub mod ir;
pub mod scenario;
pub mod signal;

pub use bc::{BcKind, BcValue, BoundaryCondition, Compat, Physics};
pub use ensemble::{
    DEFAULT_REALIZATION_BUDGET, Realization, RealizationBudget, SpectrumModel, StochasticEnsemble,
};
pub use frame::{Frame, FrameId, FrameMotion, FrameTree, WORLD};
pub use scenario::{
    Combination, ContactLaw, ContactModel, DEFAULT_VALIDATION_BUDGET, Environment, LoadCase,
    Scenario, ValidationBudget, ValidationError, ValidationPlan, Violation,
};
pub use signal::{ChebProfile, Interp, TimeSignal};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured scenario failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum ScenarioError {
    /// A value's dimensions disagree with what the context demands.
    Dimensions {
        /// What was being evaluated.
        context: String,
        /// Expected SI exponents (m, kg, s, K, A, mol).
        expected: [i8; 6],
        /// Supplied exponents.
        got: [i8; 6],
    },
    /// A frame reference could not be resolved.
    Frame {
        /// Diagnosis.
        what: String,
    },
    /// An evaluation was structurally impossible (empty table, bad time).
    Evaluate {
        /// Diagnosis.
        what: String,
    },
    /// Canonical-IR text failed to parse.
    Parse {
        /// Byte offset of the failure.
        at: usize,
        /// Diagnosis.
        what: String,
    },
}

impl fmt::Display for ScenarioError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ScenarioError::Dimensions {
                context,
                expected,
                got,
            } => write!(
                f,
                "{context}: dimension mismatch (expected {expected:?}, got {got:?})"
            ),
            ScenarioError::Frame { what } => write!(f, "frame error: {what}"),
            ScenarioError::Evaluate { what } => write!(f, "evaluation error: {what}"),
            ScenarioError::Parse { at, what } => write!(f, "IR parse error at byte {at}: {what}"),
        }
    }
}

impl std::error::Error for ScenarioError {}
