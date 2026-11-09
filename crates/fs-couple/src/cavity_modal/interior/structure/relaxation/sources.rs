//! Reusable original force/flow spectra for complete-state resonance queries.

use super::{
    StructuralInteriorRelaxationHarmonicError, StructuralInteriorSourceParameterJet,
    StructuralInteriorSources,
};

/// Derivatives needed by one bounded physical response evaluation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StructuralInteriorRelaxationSourceOrder {
    /// Values and first angular-frequency derivatives (two state solves).
    Frequency,
    /// Values and first/second frequency derivatives (three state solves).
    Curvature,
    /// All five frequency/physical-parameter orders (five state solves).
    Parameter,
}

/// One coherent physical source law, evaluated at positive omega [rad/s].
///
/// A provider may retain and reuse its own arrays; the returned borrowed jet
/// lives only through the current solve. Derivatives use the original force
/// and inward-flow coordinates and peak exp(-i*omega*t) phasors. None declares
/// an exactly zero derivative, never an unknown one. Requested orders must
/// describe the SAME differentiable source law on the entire query interval.
/// In particular, the final curvature/parameter call must reproduce values
/// and slopes previously supplied by a frequency call.
///
/// The order permits avoiding unused source derivatives during refinement.
/// The provider owns the domain and compute bounds of its source law; the
/// query polls cancellation before and after each call. A refusal publishes
/// no partial extremum, linewidth or survey. No source spectrum is inferred
/// from material parameters or from a desired resonance frequency.
pub trait StructuralInteriorRelaxationSourceSpectrum {
    fn jet(
        &mut self,
        angular_frequency_rad_s: f64,
        order: StructuralInteriorRelaxationSourceOrder,
    ) -> Result<StructuralInteriorSourceParameterJet<'_>, StructuralInteriorRelaxationHarmonicError>;
}

/// Fixed coherent force/flow phasors require no spectrum allocation.
impl StructuralInteriorRelaxationSourceSpectrum for StructuralInteriorSources<'_> {
    fn jet(
        &mut self,
        _angular_frequency_rad_s: f64,
        _order: StructuralInteriorRelaxationSourceOrder,
    ) -> Result<StructuralInteriorSourceParameterJet<'_>, StructuralInteriorRelaxationHarmonicError>
    {
        Ok(StructuralInteriorSourceParameterJet::fixed(*self))
    }
}
