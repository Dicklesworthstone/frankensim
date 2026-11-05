//! Safe-Rust, solver-independent finite-patch normal-response laws.

/// Source-bound, displacement-controlled JKR adhesion with explicit jump events.
pub mod adhesive;
/// Source-bound low-Tabor DMT contact with a conservative finite-range tail.
pub mod dmt;
mod binding;
mod embed;
mod finite_gap;
mod law;
mod response_curve;
/// Source-bound point contact shared by generic potential and port consumers.
pub mod source;
/// Original material hardness driving irreversible finite-gap indentation.
pub mod plastic_source;

pub use binding::*;
pub use embed::*;
pub use finite_gap::*;
pub use law::*;
pub use response_curve::*;

/// Finite-zone Maugis--Dugdale adhesive sphere contact.
pub mod maugis;

/// Source-preserving spatial tractions and conservative annular patch transfer.
pub mod traction;

/// Signed-gap source-law ports with conservative work and atomic receiver admission.
pub mod source_port;
