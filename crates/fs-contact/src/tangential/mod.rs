//! Solver-facing tangential contact adapters.

mod linearization;

pub mod smooth;
/// Coupled, unregularized partial-slip islands with independent physical gates.
pub mod nonsmooth;

/// Exact ideal-Coulomb impulses for decoupled isotropic tangent mobility.
pub mod impact;
