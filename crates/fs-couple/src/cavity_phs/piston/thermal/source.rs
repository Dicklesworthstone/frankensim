//! Compatibility path for the shared material-owned thermal condition binding.
//!
//! Electrical, conductive and radiative consumers use these identical types.
//! No source query, admission rule or physical state is duplicated here.
pub use fs_material::thermal::source::{ThermalMaterialError, ThermalMaterialPoint};
