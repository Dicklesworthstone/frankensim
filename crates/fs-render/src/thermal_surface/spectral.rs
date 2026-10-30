//! Visible-band Kirchhoff emission and reflection from one accepted source law.
//!
//! These opaque diffuse surfaces declare wavelength-dependent emissivity over
//! the renderer's visible band. Their reciprocal reflectance is exactly
//! `1 - emissivity(lambda)`, and their equilibrium radiance is
//! `emissivity(lambda) * Planck(lambda, temperature)`. The declared law does not
//! establish measured-material authority, infrared behavior or total cooling.

use super::{
    BOLTZMANN_J_K, ContentHash, DomainHasher, GrayDiffuseThermalSurface, LIGHT_M_S,
    ThermalSurfaceError, planck_radiance_w_m2_sr_nm,
};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};

const MAX_KNOTS: usize = 16;

/// Bounded piecewise-linear emissivity, including both visible-band endpoints.
///
/// The accepted law covers exactly 380–780 nm with 2–16 strictly increasing
/// wavelength knots and finite emissivity in `[0,1]`. It never extrapolates,
/// interpolates RGB values or implies a law outside the declared band.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralEmissivity {
    knots: [(f64, f64); MAX_KNOTS],
    len: u8,
    emitting: bool,
}

impl SpectralEmissivity {
    /// Declare constant visible-band emissivity, with no broadband inference.
    pub fn constant(emissivity: f64) -> Result<Self, ThermalSurfaceError> {
        Self::table(&[(LAMBDA_MIN, emissivity), (LAMBDA_MAX, emissivity)])
    }

    /// Admit `(wavelength_nm, emissivity)` pairs without sorting or clipping.
    ///
    /// Endpoints must be exactly 380 and 780 nm. At most 16 knots are copied
    /// into inline storage; no allocation or caller-owned table survives.
    pub fn table(knots: &[(f64, f64)]) -> Result<Self, ThermalSurfaceError> {
        if !(2..=MAX_KNOTS).contains(&knots.len())
            || knots[0].0 != LAMBDA_MIN
            || knots[knots.len() - 1].0 != LAMBDA_MAX
        {
            return Err(ThermalSurfaceError::SpectralEmissivity);
        }
        let mut accepted = Self {
            knots: [(0.0, 0.0); MAX_KNOTS],
            len: knots.len() as u8,
            emitting: false,
        };
        for (index, &(wavelength_nm, emissivity)) in knots.iter().enumerate() {
            if !wavelength_nm.is_finite()
                || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm)
                || (index > 0 && wavelength_nm <= knots[index - 1].0)
            {
                return Err(ThermalSurfaceError::SpectralEmissivity);
            }
            if !emissivity.is_finite() || !(0.0..=1.0).contains(&emissivity) {
                return Err(ThermalSurfaceError::Emissivity);
            }
            // One physical zero has one spectral identity.
            let emissivity = if emissivity == 0.0 { 0.0 } else { emissivity };
            accepted.knots[index] = (wavelength_nm, emissivity);
            accepted.emitting |= emissivity > 0.0;
        }
        Ok(accepted)
    }

    /// Evaluate the accepted linear law at a wavelength in the closed band.
    pub fn at(&self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
            return Err(ThermalSurfaceError::Wavelength);
        }
        let knots = self.knots();
        let index = knots.partition_point(|&(wavelength, _)| wavelength <= wavelength_nm) - 1;
        let (left_wavelength, left_value) = knots[index];
        if wavelength_nm == left_wavelength {
            return Ok(left_value);
        }
        let (right_wavelength, right_value) = knots[index + 1];
        if left_value == right_value {
            return Ok(left_value);
        }
        let fraction = (wavelength_nm - left_wavelength) / (right_wavelength - left_wavelength);
        // The convex form preserves nonnegative support near a zero endpoint
        // without subtracting nearly equal emissivities or clamping the law.
        Ok((1.0 - fraction) * left_value + fraction * right_value)
    }

    /// Declared wavelength and emissivity knots, without unused inline slots.
    #[must_use]
    pub fn knots(&self) -> &[(f64, f64)] {
        &self.knots[..usize::from(self.len)]
    }

    /// Identity of the declared band, interpolation model and numerical table.
    /// Redundant knots retain distinct declaration identities.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.render.spectral-emissivity.v1");
        hash.update(&[self.len]);
        for &(wavelength_nm, emissivity) in self.knots() {
            hash.update(&wavelength_nm.to_bits().to_le_bytes());
            hash.update(&emissivity.to_bits().to_le_bytes());
        }
        hash.finalize()
    }
}

/// Accepted visible-band opaque diffuse surface in local thermal equilibrium.
///
/// The source identity binds this declared optical law and temperature to the
/// caller's accepted physical source; supplying a hash does not create measured
/// provenance. A hemispherical-total emissivity alone is insufficient to infer
/// the spectral knots. This view emits into vacuum/the exterior-air model and
/// keeps the existing two-sided, per-face convention of gray thermal surfaces.
///
/// No bolometric exitance or cooling API is provided: an infrared continuation
/// and its owning thermal model must be supplied before asserting total power.
/// Observing this spectrum never advances or debits the thermal owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralDiffuseThermalSurface {
    temperature_k: f64,
    emissivity: SpectralEmissivity,
    source_identity: ContentHash,
}

impl SpectralDiffuseThermalSurface {
    /// Bind a finite positive temperature to the accepted visible source law.
    ///
    /// An emitting surface requires a representable visible-band blackbody
    /// envelope. This is a numerical restriction, without a broadband cooling
    /// assumption; an identically zero source needs no blackbody evaluation.
    pub fn try_new(
        temperature_k: f64,
        emissivity: SpectralEmissivity,
        source_identity: ContentHash,
    ) -> Result<Self, ThermalSurfaceError> {
        if !temperature_k.is_finite() || temperature_k <= 0.0 {
            return Err(ThermalSurfaceError::Temperature);
        }
        if source_identity == ContentHash([0; 32]) {
            return Err(ThermalSurfaceError::SourceIdentity);
        }
        if emissivity.emitting {
            // exp(x)-1 >= x gives the Rayleigh–Jeans upper bound B_lambda <=
            // 2*c*k*T/lambda^4. Its largest value in this band occurs at 380 nm.
            // This check avoids both a broadband T^4 restriction and a source
            // that overflows before the accepted emissivity can be applied.
            let lambda2 = LAMBDA_MIN * LAMBDA_MIN;
            let upper_per_kelvin = 2.0 * LIGHT_M_S * BOLTZMANN_J_K * 1.0e27 / (lambda2 * lambda2);
            if temperature_k > f64::MAX / upper_per_kelvin {
                return Err(ThermalSurfaceError::Unrepresentable);
            }
            planck_radiance_w_m2_sr_nm(temperature_k, LAMBDA_MIN)?;
        }
        Ok(Self {
            temperature_k,
            emissivity,
            source_identity,
        })
    }

    /// Accepted thermodynamic temperature, K.
    #[must_use]
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    /// Accepted visible-band spectral emissivity law.
    #[must_use]
    pub const fn emissivity(&self) -> SpectralEmissivity {
        self.emissivity
    }

    /// Reciprocal opaque diffuse reflectance `1 - emissivity(lambda)`.
    pub fn reflectance(&self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        Ok(1.0 - self.emissivity.at(wavelength_nm)?)
    }

    /// Emitted spectral radiance in W/(m² sr nm), within 380–780 nm inclusive.
    pub fn radiance_w_m2_sr_nm(&self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        let emissivity = self.emissivity.at(wavelength_nm)?;
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

    /// Whether the accepted emissivity has positive support in the visible band.
    /// A sufficiently cold spectrum may still underflow numerically to zero.
    #[must_use]
    pub const fn is_emitting(&self) -> bool {
        self.emissivity.emitting
    }

    /// Source, temperature and declared spectral-diffuse-model identity.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.render.spectral-diffuse-thermal.v1");
        hash.update(self.source_identity.as_bytes());
        hash.update(&self.temperature_k.to_bits().to_le_bytes());
        hash.update(self.emissivity.identity().as_bytes());
        hash.finalize()
    }
}

/// Compact source view used by camera hits and prepared direct-light proposals.
///
/// Gray values retain their existing spectrum and identity exactly. Spectral
/// laws stay borrowed, so light triangles and samples do not copy knot tables.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DiffuseThermalSurface<'a> {
    /// Existing gray source, with its original all-wavelength behavior.
    Gray(GrayDiffuseThermalSurface),
    /// Accepted visible spectral source, borrowed from its material owner.
    Spectral(&'a SpectralDiffuseThermalSurface),
}

impl<'a> From<GrayDiffuseThermalSurface> for DiffuseThermalSurface<'a> {
    fn from(surface: GrayDiffuseThermalSurface) -> Self {
        Self::Gray(surface)
    }
}

impl<'a> From<&'a SpectralDiffuseThermalSurface> for DiffuseThermalSurface<'a> {
    fn from(surface: &'a SpectralDiffuseThermalSurface) -> Self {
        Self::Spectral(surface)
    }
}

impl DiffuseThermalSurface<'_> {
    /// Opaque reciprocal reflectance within the source's wavelength domain.
    pub fn reflectance(self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        match self {
            Self::Gray(surface) => {
                if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
                    return Err(ThermalSurfaceError::Wavelength);
                }
                Ok(surface.reflectance())
            }
            Self::Spectral(surface) => surface.reflectance(wavelength_nm),
        }
    }

    /// Emitted radiance in W/(m² sr nm) from the same accepted source law.
    pub fn radiance_w_m2_sr_nm(self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        match self {
            Self::Gray(surface) => surface.radiance_w_m2_sr_nm(wavelength_nm),
            Self::Spectral(surface) => surface.radiance_w_m2_sr_nm(wavelength_nm),
        }
    }

    /// Retain the original constitutive source identity without RGB conversion.
    #[must_use]
    pub fn identity(self) -> ContentHash {
        match self {
            Self::Gray(surface) => surface.identity(),
            Self::Spectral(surface) => surface.identity(),
        }
    }

    /// Whether the accepted source has any positive emissivity support.
    #[must_use]
    pub fn is_emitting(self) -> bool {
        match self {
            Self::Gray(surface) => surface.emissivity() > 0.0,
            Self::Spectral(surface) => surface.is_emitting(),
        }
    }
}

#[cfg(test)]
mod tests;
