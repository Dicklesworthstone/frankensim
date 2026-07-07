//! fs-regime (plan patch Rev A): the physics-regime and
//! nondimensionalization kernel — the formal layer answering "WHICH
//! SOLVER IS EVEN VALID for this physical situation?" before FLUX runs.
//! An agent can otherwise run a beautiful, fast, WRONG simulation;
//! regime checking turns that class of mistake into a structured,
//! alternatives-ranked refusal.
//!
//! Layer: L3. Runtime deps: `std`, fs-qty (dimension vectors), fs-evidence
//! (`ModelCard`/`ValidityDomain`/`Evidence`), fs-math.

pub mod cards;
pub mod groups;
pub mod pi;
pub mod report;
pub mod scaling;

pub use cards::{Admission, flux_model_cards};
pub use groups::{NamedGroup, Role, RoleInput, standard_groups};
pub use pi::{Input, PiBasis, PiGroup, pi_groups};
pub use report::{BenchmarkMatch, RegimeReport, assess};
pub use scaling::{ScalingMap, condition_number};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured regime failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum RegimeError {
    /// A quantity combination that must be dimensionless is not.
    NotDimensionless {
        /// Which group/product.
        context: String,
        /// The residual SI exponents.
        residual: [i8; 5],
    },
    /// The Pi machinery was given a degenerate input set.
    Degenerate {
        /// Diagnosis.
        what: String,
    },
    /// A required physical role is missing from the inputs.
    MissingRole {
        /// The role that was needed.
        role: &'static str,
        /// What needed it.
        context: String,
    },
    /// An unknown solver/model name was queried.
    UnknownModel {
        /// The name.
        name: String,
    },
    /// A value was non-finite or non-positive where positivity is required.
    BadValue {
        /// Diagnosis.
        what: String,
    },
}

impl fmt::Display for RegimeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RegimeError::NotDimensionless { context, residual } => {
                write!(f, "{context}: not dimensionless (residual {residual:?})")
            }
            RegimeError::Degenerate { what } => write!(f, "degenerate input set: {what}"),
            RegimeError::MissingRole { role, context } => {
                write!(f, "{context}: missing required role {role}")
            }
            RegimeError::UnknownModel { name } => write!(f, "unknown model {name:?}"),
            RegimeError::BadValue { what } => write!(f, "bad value: {what}"),
        }
    }
}

impl std::error::Error for RegimeError {}
