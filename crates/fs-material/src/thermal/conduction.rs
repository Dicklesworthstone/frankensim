//! Quasi-steady planar conduction and its material/timeline adapters.
//!
//! The shared layer kernel owns conductivity integration and the interface solve.
//! Source adapters bind it to their original caloric laws without duplicating fluxes.

mod layer;
use fs_thermal::enthalpy::DepositionError;
pub use layer::{
    ConductivityCurve, LayerConductionConfig, LayerConductionResponse, MAX_CONDUCTIVITY_KNOTS,
    TwoLayerConduction,
};

/// Endpoint coefficients for the original conservative thermal network.
pub mod network;
/// Source-qualified preparation on an original pure-phase caloric branch.
pub mod source;
