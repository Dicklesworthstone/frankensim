//! Uniform Lambertian emission attached to the scene's original mesh surfaces.
//!
//! A binding supplies total emitted lumens and the original visible spectrum.
//! Admission computes the actual world-space triangle area; the existing
//! primitive material still reflects incoming light. No display RGB value,
//! fictitious temperature or independently reconstructed surface is involved.

use super::source::{PhotometricSource, SourceError, SourceSpectrum};

/// Which geometric sides emit the supplied total luminous flux.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmissionSidedness {
    /// Emit into the hemisphere indicated by each triangle's winding normal.
    Front,
    /// Split total flux equally between the two geometric hemispheres.
    Both,
}

/// Uniform spectral emitter bound to one native mesh primitive.
///
/// For actual mesh area `A`, front-side luminance is `flux/(pi*A)` for
/// [`EmissionSidedness::Front`] and `flux/(2*pi*A)` on either side for
/// [`EmissionSidedness::Both`]. Area is measured in square metres at scene
/// admission, including the original rigid-instance placement. The represented
/// spectral/observer model is the same one used by [`SourceSpectrum`].
///
/// The binding does not replace the primitive's material or geometry. Bind a
/// separate primitive for each spatially uniform source spectrum. Mesh disks
/// remain polygonal sources; their declared flux uses their actual mesh area.
#[derive(Debug, Clone)]
pub struct PhotometricAreaEmitter {
    primitive_index: usize,
    emitted_flux_lm: f64,
    spectrum: SourceSpectrum,
    sidedness: EmissionSidedness,
}

impl PhotometricAreaEmitter {
    /// Admit finite nonnegative full-output lumens and a linear `[0,1]` dimmer.
    /// Primitive existence, area and material compatibility are checked when
    /// the source is admitted together with its borrowed scene.
    pub fn try_new(
        primitive_index: usize,
        lumens: f64,
        dimmer: f64,
        spectrum: SourceSpectrum,
        sidedness: EmissionSidedness,
    ) -> Result<Self, SourceError> {
        if !lumens.is_finite()
            || lumens < 0.0
            || !dimmer.is_finite()
            || !(0.0..=1.0).contains(&dimmer)
        {
            return Err(SourceError::InvalidFlux);
        }
        let emitted_flux_lm = lumens * dimmer;
        if !emitted_flux_lm.is_finite() {
            return Err(SourceError::NumericalRange);
        }
        Ok(Self {
            primitive_index,
            emitted_flux_lm,
            spectrum,
            sidedness,
        })
    }

    /// Primitive in the same borrowed scene used by the path job.
    #[must_use]
    pub const fn primitive_index(&self) -> usize {
        self.primitive_index
    }

    /// Total flux after the linear dimmer, summed over the emitting sides.
    #[must_use]
    pub const fn flux_lm(&self) -> f64 {
        self.emitted_flux_lm
    }

    /// Full normalized 80-bin spectrum, without an RGB approximation.
    #[must_use]
    pub const fn spectrum(&self) -> &SourceSpectrum {
        &self.spectrum
    }

    /// Emitting side convention relative to original triangle winding.
    #[must_use]
    pub const fn sidedness(&self) -> EmissionSidedness {
        self.sidedness
    }
}

/// Borrowed calibrated point fixtures and finite mesh-emitter bindings.
#[derive(Debug, Clone, Copy)]
pub struct PhotometricEmitters<'a> {
    /// Delta-position fixtures with flux-normalized angular intensity laws.
    pub points: &'a [PhotometricSource],
    /// Finite Lambertian emitters attached to original scene primitives.
    pub surfaces: &'a [PhotometricAreaEmitter],
}

impl<'a> PhotometricEmitters<'a> {
    /// Preserve the existing point-fixture source model without finite emitters.
    #[must_use]
    pub const fn points_only(points: &'a [PhotometricSource]) -> Self {
        Self {
            points,
            surfaces: &[],
        }
    }
}
