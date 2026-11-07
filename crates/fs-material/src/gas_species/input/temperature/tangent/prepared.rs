//! Retained source-cell geometry for repeated original gas tangent queries.

use super::{
    GasRecord, IsobaricGasLaw, IsobaricGasTransportTangent, MaterialError, Source, Temperature,
    invalid, thermochemical,
};

/// Reusable analytic transport view borrowing one immutable isobaric gas law.
///
/// Selected temperature/pressure surfaces retain their original prepared
/// geometry and affine gradients. Engineering gases and scalar/curve claims
/// need no additional cache allocation. Each evaluation still obtains the
/// complete original gas record, source receipts and phase/domain admission;
/// preparation does not prove interval coverage or smoothness at a source kink.
#[derive(Debug)]
pub struct PreparedIsobaricGasTransport<'g> {
    law: &'g IsobaricGasLaw,
    thermochemical: Option<thermochemical::PreparedTransport>,
}

impl IsobaricGasLaw {
    /// Prepare the original selected surface cells once for repeated tangents.
    /// The gas law remains borrowed and cannot switch source, pressure or
    /// composition. Actual temperature validity is checked by `evaluate`.
    pub fn prepare_transport_temperature_tangent(
        &self,
    ) -> Result<PreparedIsobaricGasTransport<'_>, MaterialError> {
        let thermochemical = if let Source::Thermochemical { model, .. } = &self.source {
            let GasRecord::Thermochemical(snapshot) = &self.reference else {
                return Err(invalid(
                    "thermochemical preparation lost its original snapshot",
                ));
            };
            thermochemical::PreparedTransport::try_new(model, snapshot)?
        } else {
            None
        };
        Ok(PreparedIsobaricGasTransport {
            law: self,
            thermochemical,
        })
    }
}

impl PreparedIsobaricGasTransport<'_> {
    /// Evaluate original values and their fixed-pressure/composition tangents
    /// without preparing or copying the selected surface geometry again.
    pub fn evaluate(
        &self,
        temperature: Temperature,
    ) -> Result<IsobaricGasTransportTangent, MaterialError> {
        self.law
            .transport_temperature_tangent_with(temperature, self.thermochemical.as_ref())
    }
}
