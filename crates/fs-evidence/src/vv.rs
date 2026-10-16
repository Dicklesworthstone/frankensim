//! Operational verification-and-validation artifacts.
//!
//! The schemas in this module keep context of use, experimental lineage,
//! calibration/validation separation, numerical verification, prediction
//! uncertainty, and runtime assumptions machine-checkable. Admission is
//! intentionally structural: it does not authenticate a laboratory or turn a
//! simulation comparison into physical validation.

mod codec;
mod model;

pub use codec::{MAX_VV_CANONICAL_BYTES, VvCodecError};
pub use model::*;
