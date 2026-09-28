//! Spatial DKT plate dynamics and physical-parameter adjoints.
//!
//! This code needs `fs-time` and `fs-solver`. `fs-solver` depends on
//! `fs-feec`, whose `terminal-relative` feature depends on `fs-couple`, which
//! depends on `fs-plate`. As an optional feature of `fs-plate` it therefore
//! closed a package cycle that made every workspace cargo command fail. As its
//! own crate it depends on `fs-plate`, and nothing below it depends back.
//! It reaches `fs-plate` only through the public [`PlateModel`].

pub use fs_plate::PlateModel;

mod transient;
pub use transient::*;
