//! Shared lowering seam for a reduced physical reed specimen.
//!
//! Structural producers own their material, geometry, pressure projection and
//! approximation admission. Acoustic consumers own the bore, contacts, gas and
//! time integration. Both depend on this description-level interface, never on
//! each other's solver crates. Implementing this trait is not a certification.

use crate::BeatingReed;

/// Compile an immutable specimen into the existing inward-striking reed image.
///
/// The source's own admission type makes truncation/validity criteria explicit.
/// Only mouth pressure and attack are playing controls; mechanical coefficients
/// must come from this source. No instrument name or nominal pitch participates.
pub trait BeatingReedSource {
    /// Producer-owned approximation criteria, retained by the consumer.
    type Admission: Clone + core::fmt::Debug;
    /// Original producer refusal; consumers must not replace it with defaults.
    type Error: std::error::Error + 'static;

    /// Resolve the supported playing image without mutating the specimen.
    ///
    /// # Errors
    /// Refuse unsupported approximation, nonphysical controls or source data.
    fn compile_beating_reed(
        &self,
        blowing_pressure_pa: f64,
        attack_s: f64,
        admission: &Self::Admission,
    ) -> Result<BeatingReed, Self::Error>;
}
