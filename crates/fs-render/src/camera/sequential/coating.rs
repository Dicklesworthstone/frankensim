//! Owned original spectral films for an allocation-free sequential lens path.

use std::sync::Arc;

use super::CameraError;
use crate::conductor::{
    CoatingLayer, CoherentAbsorbingCoating, CoherentCoatingStack, CoherentDielectricCoating,
    MAX_COHERENT_STACK_LAYERS,
};
use crate::dielectric::DielectricCoatingStack;
use crate::polarization::dielectric::PolarizedFresnel;
use fs_blake3::ContentHash;

#[cfg(test)]
mod tests;

/// One original transparent or absorbing film, without fitting its spectrum
/// to a compact Cauchy coating. Thickness and source authority stay with the
/// native film owner. Layers are ordered from sensor to scene on every face.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LensCoatingLayer {
    /// Original Cauchy, Sellmeier or tabulated phase-index law.
    Dielectric(CoherentDielectricCoating),
    /// Original tabulated real index and extinction coefficient.
    Absorbing(CoherentAbsorbingCoating),
}

impl LensCoatingLayer {
    /// Borrow this exact film for the native general multilayer evaluator.
    #[must_use]
    pub const fn as_borrowed(&self) -> CoatingLayer<'_> {
        match self {
            Self::Dielectric(film) => CoatingLayer::Dielectric(film),
            Self::Absorbing(film) => CoatingLayer::Absorbing(film),
        }
    }

    /// Closed support of the original film law in vacuum nanometres.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.as_borrowed().wavelength_bounds_nm()
    }

    /// Physical thickness normal to the face, in metres.
    #[must_use]
    pub const fn thickness_m(&self) -> f64 {
        match self {
            Self::Dielectric(film) => film.thickness_m(),
            Self::Absorbing(film) => film.thickness_m(),
        }
    }
}

impl From<CoherentDielectricCoating> for LensCoatingLayer {
    fn from(film: CoherentDielectricCoating) -> Self {
        Self::Dielectric(film)
    }
}

impl From<CoherentAbsorbingCoating> for LensCoatingLayer {
    fn from(film: CoherentAbsorbingCoating) -> Self {
        Self::Absorbing(film)
    }
}

/// Immutable shared stack of original spectral films for a real lens face.
///
/// One through 64 ordered layers retain their complete laws, physical
/// thicknesses and source identities. Their common support may extend outside
/// the visible band; no extrapolation or resampling is performed. Even a
/// zero-thickness film retains its source domain. Cloning shares the admitted
/// layer data. Tracing uses bounded borrowed descriptors with no allocation
/// and delegates interference, absorption and signed s/p phase to the native
/// [`DielectricCoatingStack`] evaluator.
///
/// Films are locally parallel and collapsed onto the lens face; they change
/// the reciprocal polarization weight, not the bulk Snell ray geometry.
/// Nonzero films require steady tracing because coherent group residence
/// delay is not modeled by the sequential camera.
#[derive(Clone, Debug, PartialEq)]
pub struct LensCoating {
    layers: Arc<[LensCoatingLayer]>,
    wavelength_bounds_nm: (f64, f64),
    identity: ContentHash,
    active: bool,
}

impl LensCoating {
    /// Admit a nonempty stack over the exact intersection of its source bands.
    /// An uncoated face needs no stack; an authored zero-thickness stack still
    /// carries a source and a wavelength domain. Disjoint bands refuse.
    pub fn try_new(layers: &[LensCoatingLayer]) -> Result<Self, CameraError> {
        if layers.is_empty() || layers.len() > MAX_COHERENT_STACK_LAYERS {
            return Err(CameraError::InvalidSequentialLens);
        }
        let band = layers.iter().fold((0.0_f64, f64::INFINITY), |band, film| {
            let support = film.wavelength_bounds_nm();
            (band.0.max(support.0), band.1.min(support.1))
        });
        let descriptors = borrowed_layers(layers);
        let borrowed = &descriptors[..layers.len()];
        DielectricCoatingStack::try_mixed_in_band(borrowed, band)
            .map_err(|_| CameraError::InvalidSequentialLens)?;
        let identity = CoherentCoatingStack::try_new(borrowed)
            .map_err(|_| CameraError::InvalidSequentialLens)?
            .content_identity();
        Ok(Self {
            layers: Arc::from(layers),
            wavelength_bounds_nm: band,
            identity,
            active: layers.iter().any(|layer| layer.thickness_m() != 0.0),
        })
    }

    /// Original films in sensor-to-scene order, without copied spectral fits.
    #[must_use]
    pub fn layers(&self) -> &[LensCoatingLayer] {
        &self.layers
    }

    /// Exact common support of the original films in vacuum nanometres.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.wavelength_bounds_nm
    }

    /// Native ordered-stack identity, including complete optical laws,
    /// thicknesses, source identities and their declared authority.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    pub(super) const fn is_active(&self) -> bool {
        self.active
    }

    pub(super) fn fresnel_polarized(
        &self,
        wavelength_nm: f64,
        reverse: bool,
        cosine: f64,
        eta_i: f64,
        eta_t: f64,
    ) -> Result<PolarizedFresnel, CameraError> {
        let descriptors = borrowed_layers(&self.layers);
        DielectricCoatingStack::try_mixed_in_band(
            &descriptors[..self.layers.len()],
            self.wavelength_bounds_nm,
        )
        .and_then(|stack| stack.fresnel_polarized(wavelength_nm, reverse, cosine, eta_i, eta_t))
        .map_err(|_| CameraError::SequentialLensOutsideModel)
    }
}

// Admission guarantees a nonempty stack with at most 64 layers. Only small
// borrowed descriptors enter the ray stack; the original spectral tables
// remain in the shared immutable allocation.
fn borrowed_layers(layers: &[LensCoatingLayer]) -> [CoatingLayer<'_>; MAX_COHERENT_STACK_LAYERS] {
    let mut descriptors = [layers[0].as_borrowed(); MAX_COHERENT_STACK_LAYERS];
    for (descriptor, layer) in descriptors.iter_mut().zip(layers) {
        *descriptor = layer.as_borrowed();
    }
    descriptors
}
