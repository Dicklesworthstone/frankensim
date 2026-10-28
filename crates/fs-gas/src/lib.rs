#![forbid(unsafe_code)]

//! L3 evolving gas systems over the original L1 thermochemical source owners.
//!
//! The `experimental-euler` feature admits the estimated, inviscid one-dimensional
//! frozen-composition Euler lane. It is a conservative shock-capturing solver,
//! not an entropy/invariant-domain certificate or a complete reactive network.

/// Crate version retained by callers in their simulation records.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(feature = "experimental-euler")]
pub mod dynamics;
