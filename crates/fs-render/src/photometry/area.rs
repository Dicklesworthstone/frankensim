//! Calibrated emission attached to the scene's original mesh surfaces.
//!
//! A binding supplies total emitted lumens and the original visible spectrum.
//! Admission computes the actual world-space triangle area; the existing
//! primitive material still reflects incoming light. No display RGB value,
//! fictitious temperature or independently reconstructed surface is involved.

use super::source::{PhotometricSource, SPECTRAL_BINS, SourceError, SourceSpectrum};

mod projector;
pub use projector::{
    GoboTransmission, MAX_GOBO_TEXELS, MAX_PROJECTOR_BEAMS, ParaxialBeamSplitter, ProjectorBeam,
    ThinLensProjector,
};

/// Which geometric sides emit the supplied total luminous flux.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmissionSidedness {
    /// Emit into the hemisphere indicated by each triangle's winding normal.
    Front,
    /// Split total flux equally between the two geometric hemispheres.
    Both,
}

/// Spectral emitter bound to one native mesh primitive.
///
/// By default, for actual mesh area `A`, front-side luminance is `flux/(pi*A)` for
/// [`EmissionSidedness::Front`] and `flux/(2*pi*A)` on either side for
/// [`EmissionSidedness::Both`]. Area is measured in square metres at scene
/// admission, including the original rigid-instance placement. The represented
/// spectral/observer model is the same one used by [`SourceSpectrum`].
///
/// The binding does not replace the primitive's material or geometry. Bind a
/// separate primitive for each spatially uniform source spectrum. Mesh disks
/// remain polygonal sources; their declared flux uses their actual mesh area.
/// [`Self::with_projector`] replaces the Lambertian law with a calibrated
/// directional thin-lens light field on a planar front-emitting aperture.
#[derive(Debug, Clone)]
pub struct PhotometricAreaEmitter {
    primitive_index: usize,
    clear_flux_lm: f64,
    emitted_flux_lm: f64,
    spectrum: SourceSpectrum,
    sidedness: EmissionSidedness,
    projector: Option<ThinLensProjector>,
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
            clear_flux_lm: emitted_flux_lm,
            emitted_flux_lm,
            spectrum,
            sidedness,
            projector: None,
        })
    }

    /// Use an ideal thin-lens projector in place of Lambertian emission.
    ///
    /// The original post-dimmer, spectrally filtered flux is clear-output flux.
    /// The passive gobo and optional splitter reduce emitted flux; blocked
    /// power is not renormalized into surviving rays. Replacing a projector uses
    /// the same clear-output calibration. Only front emission is supported;
    /// scene admission checks the aperture's plane, winding and numeric range.
    pub fn with_projector(mut self, projector: ThinLensProjector) -> Result<Self, SourceError> {
        if self.sidedness != EmissionSidedness::Front {
            return Err(SourceError::InvalidDistribution);
        }
        let transmission = projector.transmitted_fraction();
        let emitted_flux_lm = self.clear_flux_lm * transmission;
        if !emitted_flux_lm.is_finite()
            || (self.clear_flux_lm > 0.0 && transmission > 0.0 && emitted_flux_lm == 0.0)
        {
            return Err(SourceError::NumericalRange);
        }
        self.emitted_flux_lm = emitted_flux_lm;
        self.projector = Some(projector);
        Ok(self)
    }

    /// Apply a spatially uniform passive spectral filter to emitted ray power.
    ///
    /// Both clear-output and currently emitted lumens decrease by the photopic
    /// transmission of the original spectrum. The surviving spectrum is
    /// normalized without replacing lost power. Gobo and splitter data are
    /// retained, and replacing optics later retains the filtered calibration.
    /// This represents achromatic optics with a uniform spectral filter, not
    /// spatially varying color gobos or wavelength-dependent beam deflection.
    pub fn filtered(&self, transmittance: [f64; SPECTRAL_BINS]) -> Result<Self, SourceError> {
        let (spectrum, fraction) = self.spectrum.filtered(transmittance)?;
        let clear_flux_lm = self.clear_flux_lm * fraction;
        let emitted_flux_lm = self.emitted_flux_lm * fraction;
        if !clear_flux_lm.is_finite()
            || !emitted_flux_lm.is_finite()
            || (fraction > 0.0
                && ((self.clear_flux_lm > 0.0 && clear_flux_lm == 0.0)
                    || (self.emitted_flux_lm > 0.0 && emitted_flux_lm == 0.0)))
        {
            return Err(SourceError::NumericalRange);
        }
        let mut result = self.clone();
        result.clear_flux_lm = clear_flux_lm;
        result.emitted_flux_lm = emitted_flux_lm;
        result.spectrum = spectrum;
        Ok(result)
    }

    /// Primitive in the same borrowed scene used by the path job.
    #[must_use]
    pub const fn primitive_index(&self) -> usize {
        self.primitive_index
    }

    /// Total flux after dimming, spectral filtering and passive optical losses,
    /// summed over emitting sides.
    #[must_use]
    pub const fn flux_lm(&self) -> f64 {
        self.emitted_flux_lm
    }

    /// Post-dimmer, spectrally filtered clear flux before gobo or splitter losses.
    #[must_use]
    pub const fn clear_flux_lm(&self) -> f64 {
        self.clear_flux_lm
    }

    /// Ideal projector law, or `None` for the original uniform Lambertian law.
    #[must_use]
    pub const fn projector(&self) -> Option<&ThinLensProjector> {
        self.projector.as_ref()
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
    /// Finite Lambertian or ideal-projector emitters on original scene primitives.
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
