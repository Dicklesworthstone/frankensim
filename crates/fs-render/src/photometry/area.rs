//! Calibrated emission attached to the scene's original mesh surfaces.
//!
//! A binding supplies total emitted lumens and the original visible spectrum.
//! Admission computes the actual world-space triangle area; the existing
//! primitive material still reflects incoming light. No display RGB value,
//! fictitious temperature or independently reconstructed surface is involved.

use super::source::{PhotometricSource, SPECTRAL_BINS, SourceError, SourceSpectrum};

mod projector;
pub mod reference;
pub use projector::{
    GoboTransmission, MAX_GOBO_TEXELS, MAX_PROJECTOR_BEAMS, ParaxialBeamSplitter,
    ParaxialBoxDiffuser, ProjectorBeam, ThinLensProjector,
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
    /// In a fixture definition, `primitive_index` names a local geometry slot.
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
        if !emitted_flux_lm.is_finite() || (lumens > 0.0 && dimmer > 0.0 && emitted_flux_lm == 0.0)
        {
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

    /// Apply an additional linear `[0,1]` dimmer to this immutable source.
    /// Both clear-output and emitted lumens are attenuated, preserving the
    /// spectrum, sidedness and any projector's existing passive optical loss.
    /// Positive flux that would underflow to zero is refused.
    pub fn dimmed(&self, dimmer: f64) -> Result<Self, SourceError> {
        if !dimmer.is_finite() || !(0.0..=1.0).contains(&dimmer) {
            return Err(SourceError::InvalidFlux);
        }
        let clear_flux_lm = self.clear_flux_lm * dimmer;
        let emitted_flux_lm = self.emitted_flux_lm * dimmer;
        if !clear_flux_lm.is_finite()
            || !emitted_flux_lm.is_finite()
            || (dimmer > 0.0
                && ((self.clear_flux_lm > 0.0 && clear_flux_lm == 0.0)
                    || (self.emitted_flux_lm > 0.0 && emitted_flux_lm == 0.0)))
        {
            return Err(SourceError::NumericalRange);
        }
        let mut result = self.clone();
        result.clear_flux_lm = clear_flux_lm;
        result.emitted_flux_lm = emitted_flux_lm;
        Ok(result)
    }

    /// Rebind an already calibrated source during fixture assembly. Geometry
    /// existence and compatibility remain native scene-admission checks.
    pub(crate) fn bound_to_primitive(&self, primitive_index: usize) -> Self {
        let mut result = self.clone();
        result.primitive_index = primitive_index;
        result
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

    /// Change the existing ideal projector's physical lens distances while
    /// retaining its original mask, splitter, diffuser, frame and calibration.
    /// Clear and emitted lumens are preserved bit-for-bit: these distances
    /// redistribute ray power without changing passive transmission. A source
    /// without projector optics refuses instead of inventing a pupil model.
    pub fn with_lens_distances(
        &self,
        gobo_distance_m: f64,
        focal_length_m: f64,
    ) -> Result<Self, SourceError> {
        let projector = self
            .projector
            .as_ref()
            .ok_or(SourceError::InvalidDistribution)?
            .with_lens_distances(gobo_distance_m, focal_length_m)?;
        let mut result = self.clone();
        result.projector = Some(projector);
        Ok(result)
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

    /// Primitive slot in the binding's declared domain: the native scene for
    /// standalone sources, fixture geometry for fixture inputs, or rig geometry
    /// for assembled sources. Native admission resolves rig slots to the same
    /// mounted geometry instances in the borrowed scene.
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::instances::RigidTransform;

    fn source(lumens: f64, dimmer: f64) -> Result<PhotometricAreaEmitter, SourceError> {
        PhotometricAreaEmitter::try_new(
            7,
            lumens,
            dimmer,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
    }

    fn projector(transmission: f64) -> ThinLensProjector {
        let gobo = GoboTransmission::try_new(0.2, 0.2, 1, 1, vec![transmission]).unwrap();
        ThinLensProjector::try_new(RigidTransform::identity(), 1.0, 0.5, gobo).unwrap()
    }

    #[test]
    fn g0_dimming_preserves_projector_loss_and_filtered_clear_calibration() {
        let original = source(128.0, 0.5)
            .unwrap()
            .with_projector(projector(0.25))
            .unwrap();
        let dimmed = original
            .filtered([0.5; SPECTRAL_BINS])
            .unwrap()
            .dimmed(0.25)
            .unwrap();
        assert_eq!(original.clear_flux_lm(), 64.0);
        assert_eq!(original.flux_lm(), 16.0);
        assert_eq!(dimmed.clear_flux_lm(), 8.0);
        assert_eq!(dimmed.flux_lm(), 2.0);
        assert_eq!(dimmed.projector().unwrap().transmitted_fraction(), 0.25);
        assert_eq!(dimmed.sidedness(), EmissionSidedness::Front);
        assert_eq!(
            dimmed.spectrum().bins_w_per_lumen_nm(),
            original.spectrum().bins_w_per_lumen_nm()
        );
        let rebound = dimmed.bound_to_primitive(19);
        assert_eq!(dimmed.primitive_index(), 7);
        assert_eq!(rebound.primitive_index(), 19);
        assert_eq!(rebound.flux_lm(), dimmed.flux_lm());
        let replaced = rebound.with_projector(projector(1.0)).unwrap();
        assert_eq!(replaced.clear_flux_lm(), 8.0);
        assert_eq!(replaced.flux_lm(), 8.0);
    }

    #[test]
    fn g0_positive_dimming_underflow_refuses_and_explicit_zero_remains_valid() {
        let tiny = f64::from_bits(1);
        assert!(matches!(
            source(tiny, 0.5),
            Err(SourceError::NumericalRange)
        ));
        let clear = source(tiny, 1.0).unwrap();
        assert!(matches!(
            clear.dimmed(0.5),
            Err(SourceError::NumericalRange)
        ));
        let projected = source(2.0 * tiny, 1.0)
            .unwrap()
            .with_projector(projector(0.5))
            .unwrap();
        // The clear flux remains representable here; the emitted flux does not.
        assert!(matches!(
            projected.dimmed(0.5),
            Err(SourceError::NumericalRange)
        ));
        assert_eq!(source(tiny, 0.0).unwrap().flux_lm(), 0.0);
        assert_eq!(
            source(0.0, 1.0).unwrap().dimmed(0.5).unwrap().flux_lm(),
            0.0
        );
        let dark = projected.dimmed(0.0).unwrap();
        assert_eq!(dark.clear_flux_lm(), 0.0);
        assert_eq!(dark.flux_lm(), 0.0);
        assert!(dark.projector().is_some());
        for invalid in [f64::NAN, -0.1, 1.1] {
            assert!(matches!(
                clear.dimmed(invalid),
                Err(SourceError::InvalidFlux)
            ));
        }
    }
}
