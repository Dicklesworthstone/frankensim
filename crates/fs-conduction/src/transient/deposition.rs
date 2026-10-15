//! Compatibility facade for the shared conservative enthalpy timeline.
//!
//! The implementation and every public state/trait now have one owner in
//! `fs-thermal`. Re-exporting preserves Rust type identity and the original
//! conduction API while allowing electrical and mechanical adapters to use the
//! timeline without pulling in the FEM/FEEC/adjoint dependency graph.
//!
//! Constitutive adapters remain with the material type that owns the law;
//! executor, implicit FEM, bonding and fracture adapters remain in
//! [`super::enthalpy::deposition`]. No physical or numerical law is replaced.

pub use fs_thermal::enthalpy::*;
