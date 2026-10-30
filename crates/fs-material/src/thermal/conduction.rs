//! Quasi-steady conduction and its material/timeline adapters.
//!
//! The shared layer kernel owns conductivity integration and the interface solve.
//! Source adapters bind it to their original caloric laws without duplicating fluxes.

mod layer;
use fs_thermal::enthalpy::DepositionError;
pub use layer::{
    ConductivityCurve, LayerConductionConfig, LayerConductionResponse, MAX_CONDUCTIVITY_KNOTS,
    TwoLayerConduction,
};

mod stack;
pub use stack::{LayerStackConduction, LayerStackResponse, MAX_CONDUCTION_LAYERS};

/// Physical geometry lowered to the same conductivity-integral flux owner.
mod geometry;
pub use geometry::LayerGeometry;

/// Endpoint coefficients for the original conservative thermal network.
pub mod network;
/// Source-qualified preparation on an original pure-phase caloric branch.
pub mod source;
