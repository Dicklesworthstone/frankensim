//! Safe-Rust, solver-independent finite-patch normal-response laws.

/// Source-bound, displacement-controlled JKR adhesion with explicit jump events.
pub mod adhesive;
mod binding;
mod embed;
mod finite_gap;
mod law;
mod response_curve;
/// Source-bound point contact shared by generic potential and port consumers.
pub mod source;

pub use binding::*;
pub use embed::*;
pub use finite_gap::*;
pub use law::*;
pub use response_curve::*;
