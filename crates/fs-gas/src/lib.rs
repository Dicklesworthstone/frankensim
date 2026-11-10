#![forbid(unsafe_code)]

//! L3 evolving gas systems over the original L1 thermochemical source owners.
//!
//! The `fusion` owner evolves a bounded rigid D-T plasma parcel in the common
//! ionized-gas state, with explicit alpha thermalization and escaping neutrons.
//!
//! The `experimental-euler` feature admits estimated, inviscid one-dimensional
//! frozen and multicomponent Euler flow, including finite-rate reaction coupling.
//! These are conservative shock-capturing engineering solvers, not fully discrete
//! entropy/invariant-domain certificates or complete reactive networks.

/// Crate version retained by callers in their simulation records.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub mod chemistry;

/// Closed source-backed gas storage with pressure-work and heat ports.
#[cfg(feature = "component")]
pub mod component;

/// Reacting adiabatic throttles, compressors and expanders over original sources.
pub mod equilibrium_flow;

// Shared source-owned converging-nozzle pressure/sonic search.
mod nozzle;

/// Finite D-T fuel, helium ash and two-temperature nuclear energy accounting.
pub mod fusion;

/// Conservative distributed mechanical forcing of the original 1-D gas fields.
#[cfg(feature = "experimental-euler")]
pub mod body_force;

#[cfg(feature = "experimental-euler")]
pub mod dynamics;

/// Conservative transport of evolving species with their original caloric laws.
#[cfg(feature = "experimental-euler")]
pub mod multicomponent;

/// Finite-rate updates of the existing source-bearing material gas parcel.
#[cfg(feature = "material")]
pub mod material;
