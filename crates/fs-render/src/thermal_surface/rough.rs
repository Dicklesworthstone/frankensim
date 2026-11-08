//! Directional Kirchhoff emission from the original EON diffuse facets.
//!
//! Facet albedo is supplied explicitly. A measured hemispherical reflectance
//! or emissivity is not interpreted as microscopic facet data. At equilibrium
//! the directional emissivity is the EON directional absorptance, so reflected
//! blackbody illumination plus emitted radiance equals that same blackbody in
//! every direction. The roughness and multiple-return closure are those of
//! [`RoughDiffuseSurface`], not a resolved microgeometry or specimen claim.

use super::{ContentHash, DomainHasher, ThermalSurfaceError, planck_radiance_w_m2_sr_nm};
use crate::rough_diffuse::RoughDiffuseSurface;
use crate::spectral::SpectralReflectance;

/// An isothermal opaque rough diffuse source in an explicit finite band.
///
/// Vacuum wavelengths are nanometres and temperature is kelvin. Each physical
/// face emits into its own hemisphere. This source uses the original vacuum
/// Planck law; refractive hosts require a separate thermodynamic source model.
/// Observing it never advances heat or supplies an unmeasured bolometric tail.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RoughDiffuseThermalSurface {
    temperature_k: f64,
    facet_albedo: SpectralReflectance,
    reflection: RoughDiffuseSurface,
    source_identity: ContentHash,
}

impl RoughDiffuseThermalSurface {
    /// Bind explicit diffuse-facet single-scattering albedo to the original
    /// source temperature and identity. No emissivity-to-facet fit is made.
    pub fn try_new(
        temperature_k: f64,
        facet_albedo: SpectralReflectance,
        reflection: RoughDiffuseSurface,
        source_identity: ContentHash,
    ) -> Result<Self, ThermalSurfaceError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(ThermalSurfaceError::Temperature);
        }
        if source_identity == ContentHash([0; 32]) {
            return Err(ThermalSurfaceError::SourceIdentity);
        }
        if facet_albedo.knots().iter().any(|&(_, rho)| rho < 1.0) {
            super::spectral::admit_source_band_temperature(
                temperature_k,
                facet_albedo.wavelength_bounds_nm().0,
            )?;
        }
        Ok(Self {
            temperature_k,
            facet_albedo,
            reflection,
            source_identity,
        })
    }

    /// Accepted thermodynamic temperature, K.
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Original spectral facet albedo, without a macroscopic-reflectance fit.
    #[must_use]
    pub const fn facet_albedo(&self) -> &SpectralReflectance {
        &self.facet_albedo
    }

    /// The same admitted EON law used for native reciprocal reflection.
    #[must_use]
    pub const fn reflection(&self) -> RoughDiffuseSurface {
        self.reflection
    }

    /// Inclusive support in vacuum nanometres.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.facet_albedo.wavelength_bounds_nm()
    }

    /// Kirchhoff directional emissivity equals directional absorptance.
    /// The stable EON formula retains weak absorption near a white reflector.
    pub fn directional_emissivity(
        &self,
        wavelength_nm: f64,
        outgoing_cosine: f64,
    ) -> Result<f64, ThermalSurfaceError> {
        let rho = self
            .facet_albedo
            .at(wavelength_nm)
            .map_err(|_| ThermalSurfaceError::Wavelength)?;
        self.reflection
            .directional_absorptance(outgoing_cosine, rho)
            .map_err(|_| ThermalSurfaceError::Direction)
    }

    /// Projected-hemisphere mean emissivity, without wavelength extrapolation.
    pub fn hemispherical_emissivity(&self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        let rho = self
            .facet_albedo
            .at(wavelength_nm)
            .map_err(|_| ThermalSurfaceError::Wavelength)?;
        self.reflection
            .hemispherical_absorptance(rho)
            .map_err(|_| ThermalSurfaceError::Emissivity)
    }

    /// Emitted radiance in W/(m² sr nm) at the actual outward direction.
    pub fn directional_radiance_w_m2_sr_nm(
        &self,
        wavelength_nm: f64,
        outgoing_cosine: f64,
    ) -> Result<f64, ThermalSurfaceError> {
        self.emitted_radiance(
            wavelength_nm,
            self.directional_emissivity(wavelength_nm, outgoing_cosine)?,
        )
    }

    /// Projected-hemisphere mean radiance. This is suitable for source
    /// selection and exitance, but must not replace directional path radiance.
    pub fn hemispherical_mean_radiance_w_m2_sr_nm(
        &self,
        wavelength_nm: f64,
    ) -> Result<f64, ThermalSurfaceError> {
        self.emitted_radiance(wavelength_nm, self.hemispherical_emissivity(wavelength_nm)?)
    }

    /// Spectral exitance of one emitting face, W/(m² nm).
    pub fn spectral_exitance_w_m2_nm(
        &self,
        wavelength_nm: f64,
    ) -> Result<f64, ThermalSurfaceError> {
        let value =
            core::f64::consts::PI * self.hemispherical_mean_radiance_w_m2_sr_nm(wavelength_nm)?;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(ThermalSurfaceError::Unrepresentable)
        }
    }

    fn emitted_radiance(
        &self,
        wavelength_nm: f64,
        emissivity: f64,
    ) -> Result<f64, ThermalSurfaceError> {
        if emissivity == 0.0 {
            return Ok(0.0);
        }
        let value = emissivity * planck_radiance_w_m2_sr_nm(self.temperature_k, wavelength_nm)?;
        if value.is_finite() && value >= 0.0 {
            Ok(value)
        } else {
            Err(ThermalSurfaceError::Unrepresentable)
        }
    }

    /// Whether any original source knot has nonzero absorption/emission.
    #[must_use]
    pub fn is_emitting(&self) -> bool {
        self.facet_albedo.knots().iter().any(|&(_, rho)| rho < 1.0)
    }

    /// Original source, temperature, model, roughness, and facet-spectrum identity.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.render.rough-diffuse-thermal.v1");
        hash.update(self.source_identity.as_bytes());
        hash.update(&self.temperature_k.to_bits().to_le_bytes());
        hash.update(&self.reflection.roughness().to_bits().to_le_bytes());
        hash.update(&[self.facet_albedo.knots().len() as u8]);
        for &(wavelength, rho) in self.facet_albedo.knots() {
            hash.update(&wavelength.to_bits().to_le_bytes());
            hash.update(&rho.to_bits().to_le_bytes());
        }
        hash.finalize()
    }
}

#[cfg(test)]
#[path = "rough/tests.rs"]
mod tests;
