//! Compatibility adapter; all gray/spectral phase arithmetic has one native owner.
use super::*;

pub(super) fn enclose(phase: &PhaseFunction, mu: Interval, cx: &Cx<'_>) -> Result<Interval, ScatterEnclosureError> {
    phase.enclosure_cosine(cx, mu).map_err(|error| match error {
        crate::volumes::transport::VolumeTransportError::Cancelled =>
            ScatterEnclosureError::Transport(PhotometryError::Cancelled),
        crate::volumes::transport::VolumeTransportError::Wavelength => ScatterEnclosureError::UnsupportedPhase,
        _ => ScatterEnclosureError::NumericalRange,
    })
}
