//! fs-selfknow-e2e — the STRUCTURE & SELF-KNOWLEDGE battery (bead
//! knh1.7): Layer 4 exercised as ONE runnable script. The six stages
//! (interface types, symmetry, spectral health, abstraction ladder,
//! explanation objects, value-of-information) live in `tests/suite.rs`
//! behind the `selfknow-e2e` feature; this crate body carries only the
//! shared identity so the workspace can name the lane.

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// The suite's stage names, in battery order (ledger keys).
pub const STAGES: [&str; 6] = [
    "interface-types",
    "symmetry-harvest",
    "spectral-health",
    "abstraction-ladder",
    "explanation-objects",
    "value-of-information",
];
