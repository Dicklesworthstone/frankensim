#![forbid(unsafe_code)]

//! L3 evolving gas systems over the original L1 thermochemical source owners.
//!
//! The `experimental-euler` feature admits estimated, inviscid one-dimensional
//! frozen and multicomponent Euler flow, including finite-rate reaction coupling.
//! These are conservative shock-capturing engineering solvers, not fully discrete
//! entropy/invariant-domain certificates or complete reactive networks.

/// Crate version retained by callers in their simulation records.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub mod chemistry;

#[cfg(feature = "experimental-euler")]
pub mod dynamics;

/// Conservative transport of evolving species with their original caloric laws.
#[cfg(feature = "experimental-euler")]
pub mod multicomponent;

/// Finite-rate updates of the existing source-bearing material gas parcel.
#[cfg(feature = "material")]
pub mod material;
