//! Shared conservative thermal transport, independent of any spatial solver.
//!
//! The enthalpy timeline owns delivered mass, thermal storage, contact exchange,
//! source work, event order and cancellation. Material crates implement its
//! [`enthalpy::ThermalLaw`] using their existing constitutive charts. Electrical,
//! mechanical and FEM consumers consequently share the SAME history and types,
//! without depending on each other's discretizations or constitutive models.
//!
//! This is the existing reduced, fixed-geometry control-volume model, not a
//! continuum FEM solve, an equation of state, or a melt-flow model. Its original
//! cell-uniformity, timestep and supplied-conductance assumptions still apply.

pub mod enthalpy;
pub mod observers;
/// Mass-owned layer/cylinder/sphere fronts and thermostatted extinction.
pub mod phase_geometry;
/// Conservative solid/vapor exchange and implicit latent-energy feedback.
pub mod sublimation;
/// Continuum gas-film evaporation, condensation and original-chart latent cooling.
pub mod evaporation;

/// Compatibility path for the original conduction-owned timeline. This is a
/// re-export, not a second implementation or a conversion between state types.
pub mod transient {
    pub use crate::enthalpy as deposition;
}
