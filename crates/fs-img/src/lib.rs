//! fs-img: in-house image plumbing (plan §10.5) — PNG and OpenEXR
//! writers/readers, an à-trous denoiser whose outputs are PERMANENTLY
//! labeled biased, and deterministic film/display transforms. Pure Rust
//! (P1), byte-exact deterministic encodes (P2).
//!
//! Layer: L5 (LUMEN). Runtime deps: `std`, fs-math (deterministic
//! `pow`/`exp` for the display transforms).

pub mod denoise;
pub mod exr;
pub mod film;
pub mod png;

pub use denoise::{DenoiseParams, LabeledPlane, PixelProvenance, atrous_denoise, mse};
pub use exr::{
    Channel, DecodedExr, PixelType, f16_bits_to_f32, f32_to_f16_bits, read_exr, write_exr,
};
pub use png::{DecodedPng, PngColor, read_png, write_png8, write_png16};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured image-plumbing failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ImgError {
    /// A buffer length disagrees with the declared shape.
    Shape {
        /// Expected element count.
        expected: usize,
        /// Supplied count.
        got: usize,
        /// What was being shaped.
        context: &'static str,
    },
    /// Structurally invalid bytes (corruption caught, never decoded
    /// silently).
    Malformed {
        /// Diagnosis.
        what: String,
    },
    /// Valid-looking bytes outside our implemented subset.
    Unsupported {
        /// What feature was encountered.
        what: String,
    },
}

impl fmt::Display for ImgError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ImgError::Shape {
                expected,
                got,
                context,
            } => {
                write!(f, "{context}: expected {expected} elements, got {got}")
            }
            ImgError::Malformed { what } => write!(f, "malformed image data: {what}"),
            ImgError::Unsupported { what } => write!(
                f,
                "unsupported: {what} — fs-img readers cover fs-img writers' subset \
                 (CONTRACT.md no-claims)"
            ),
        }
    }
}

impl std::error::Error for ImgError {}
