//! Layered thermal coupling.
mod api { include!("layered_feedback_api.rs"); }
mod state {
    include!("layered_feedback_error.rs");
    include!("layered_feedback_result.rs");
}
mod types {
    pub(super) use super::api::*;
    pub(super) use super::state::*;
}
mod solve { include!("layered_feedback_impl.rs"); }

pub use api::{LayerThermalBinding,ThermalCoefficients,StripEnthalpyConfig};
pub use state::{StripEnthalpyError,StripEnthalpyReport,StripEnthalpyCandidate};
pub use solve::step_candidate;
