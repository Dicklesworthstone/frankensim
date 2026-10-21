//! Layered thermal coupling.
#[path = "layered_feedback_api.rs"]
mod api;
#[path = "layered_feedback_error.rs"]
mod error;
#[path = "layered_feedback_result.rs"]
mod result;
mod types {
    pub(super) use super::api::*;
    pub(super) use super::error::*;
    pub(super) use super::result::*;
}
#[path = "layered_feedback_impl.rs"]
mod solve;

pub use api::{LayerThermalBinding,ThermalCoefficients,StripEnthalpyConfig};
pub use error::StripEnthalpyError;
pub use result::{StripEnthalpyReport,StripEnthalpyCandidate};
pub use solve::step_candidate;

/// Atomic coupled intervals and actual scheduled layer births.
pub mod growth;

/// Pinned material conductivity in the original self-heating residual.
pub mod transport;

/// Spatial cell/layer self-heating and gravity sag on the original bridge owner.
pub mod bridge;
