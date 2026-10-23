//! The offline `grand_piano` reference's OWN source modules, compiled into this
//! crate unmodified (the in-tree precedent is `crates/fs-couple/examples/
//! piano_exterior.rs`). Nothing here is a copy: every `#[path]` points at the
//! reference file, so the live engine and the offline renderer share
//!
//! * the 88-key RT-0425 string scale and per-key hammer cards (`steinway_scale`),
//! * the stiff-string partial law and strike projection (`geometry::Course`),
//! * the Model D soundboard reconstruction and its FE eigen-solve
//!   (`steinway_d`, `board_geometry`),
//! * the energy-discrete WoolFelt contact adapter (`felt`),
//! * the felt pack's per-key hammers, Stulov's 88 measured Abel hammers (`hammer_felt`),
//! * the Prony creep realization (`relaxation`) and the geometry-derived
//!   hammer shank (`shank`),
//! * the reference string/board bank (`linear`), used here for its
//!   `BoardMode` type and, natively, as the validation oracle.
//!
//! The reference `engine.rs` (the complete `Instrument`) is NOT part of the
//! live path: measured at 48 kHz x 1 substep it costs about 22 ns per string
//! mode per step, i.e. roughly 450 string modes at real time natively.
#![allow(dead_code, unused_imports, clippy::all)]

#[path = "../../../../crates/fs-couple/examples/grand_piano/geometry.rs"]
pub mod geometry;
#[path = "../../../../crates/fs-couple/examples/grand_piano/steinway_scale.rs"]
pub mod steinway_scale;
#[path = "../../../../crates/fs-couple/examples/grand_piano/stretch_tuning.rs"]
pub mod stretch_tuning;
#[path = "../../../../crates/fs-couple/examples/grand_piano/linear.rs"]
pub mod linear;
#[path = "../../../../crates/fs-couple/examples/grand_piano/board_geometry.rs"]
pub mod board_geometry;
#[path = "../../../../crates/fs-couple/examples/grand_piano/steinway_d.rs"]
pub mod steinway_d;
#[path = "../../../../crates/fs-couple/examples/grand_piano/felt.rs"]
pub mod felt;
#[path = "../../../../crates/fs-couple/examples/grand_piano/hammer_felt.rs"]
pub mod hammer_felt;
#[path = "../../../../crates/fs-couple/examples/grand_piano/felt_relaxation.rs"]
pub mod relaxation;
#[path = "../../../../crates/fs-couple/examples/grand_piano/hammer_shank.rs"]
pub mod shank;
/// The complete offline reference engine (`Instrument`). Native only: it is
/// the validation oracle of `examples/validate.rs`, never the live path.
#[cfg(not(target_arch = "wasm32"))]
#[path = "../../../../crates/fs-couple/examples/grand_piano/engine.rs"]
pub mod engine;
#[cfg(not(target_arch = "wasm32"))]
#[path = "../../../../crates/fs-couple/examples/grand_piano/hammer_materials.rs"]
pub mod hammer_materials;

/// Per-partial intrinsic string-loss coefficients of a prepared reference bank
/// (`pub(super)` in `linear.rs`, hence this accessor lives in its parent
/// module). The live kernel takes them from the bank instead of restating the
/// authored loss card.
pub fn bank_string_losses(bank: &linear::Bank) -> (&[linear::string_damping::ModeLoss], &[f64]) {
    (&bank.string_loss.modes, &bank.string_loss.endpoint_drag)
}
