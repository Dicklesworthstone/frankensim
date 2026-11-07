//! Original sharp fusion on the shared molecular and native thermal timeline.
//!
//! The ordinary entry point retains a solid interface throughout. The explicit
//! [`boundary`] entry point stops the same physical owners at complete melting
//! and returns the time left for a caller-selected liquid interface.

pub use super::{GasTimelineEquilibriumSurface, advance_deposition_with_equilibrium};

/// Conservative complete melting on the original heat, reflection and gas clock.
#[path = "equilibrium/boundary.rs"]
pub mod boundary;
