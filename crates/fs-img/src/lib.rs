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

mod cinematic_color;

pub use cinematic_color::{
    CINEMATIC_COLOR_CONFIG_CANONICAL_BYTES, CINEMATIC_COLOR_PIPELINE_VERSION, CinematicColorConfig,
    CinematicColorError, CinematicColorLimits, CinematicColorMetadata, CinematicDisplayTarget,
    CinematicGamutMap, CinematicInputValueClass, CinematicNegativePolicy, CinematicPreview,
    CinematicPreviewAuthority, CinematicPreviewSamples, CinematicToneCurve, CinematicWorkingSpace,
    MAX_CINEMATIC_GLARE_RADIUS_PX, MAX_CINEMATIC_PREVIEW_PIXELS, PreviewBitDepth, PreviewDither,
    PreviewGlare, transform_cinematic_preview,
};

pub use denoise::{DenoiseParams, LabeledPlane, PixelProvenance, atrous_denoise, mse};
pub use exr::{
    Channel, DecodedExr, ExrAttribute, ExrWriteLimits, ExrWriteRequirements, PixelType,
    SOURCE_ARTIFACT_HASH_ATTRIBUTE, exr_write_requirements, exr_write_requirements_for_layout,
    f16_bits_to_f32, f32_to_f16_bits, read_exr, write_exr, write_exr_with_attributes,
    write_exr_with_attributes_budgeted,
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
    /// A caller-supplied byte ceiling was too small for a requested operation.
    ResourceLimit {
        /// Stable name of the bounded resource.
        resource: &'static str,
        /// Exact logical bytes required by the operation.
        requested: u64,
        /// Caller-supplied byte ceiling.
        limit: u64,
    },
    /// The allocator refused an already admitted, fallible reservation.
    AllocationRefused {
        /// Stable name of the allocation.
        resource: &'static str,
        /// Logical bytes requested from the allocator.
        requested: u64,
    },
    /// Checked size arithmetic could not represent the requested artifact.
    SizeOverflow {
        /// Stable description of the overflowing quantity.
        context: &'static str,
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
            ImgError::ResourceLimit {
                resource,
                requested,
                limit,
            } => write!(
                f,
                "{resource} needs {requested} bytes above caller limit {limit}"
            ),
            ImgError::AllocationRefused {
                resource,
                requested,
            } => write!(
                f,
                "allocator refused {requested} bytes for {resource} after admission"
            ),
            ImgError::SizeOverflow { context } => {
                write!(f, "image size arithmetic overflowed for {context}")
            }
        }
    }
}

impl std::error::Error for ImgError {}
