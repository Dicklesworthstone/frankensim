//! Safe-Rust, solver-independent finite-patch normal-response laws.

mod embed;
mod finite_gap;
mod law;
mod response_curve;

pub use embed::*;
pub use finite_gap::*;
pub use law::*;
pub use response_curve::*;
