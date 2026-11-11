//! Source-band Kirchhoff emission and reflection from one accepted source law.
//!
//! These opaque diffuse surfaces declare wavelength-dependent emissivity over
//! their explicitly declared vacuum-wavelength band. Reciprocal reflectance is
//! `1 - emissivity(lambda)`, and their equilibrium radiance is
//! `emissivity(lambda) * Planck(lambda, temperature)`. The declared law does not
//! establish measured-material authority, out-of-band behavior or total cooling.

use super::{
    BOLTZMANN_J_K, ContentHash, DomainHasher, GrayDiffuseThermalSurface, LIGHT_M_S,
    ThermalSurfaceError, planck_radiance_w_m2_sr_nm,
};
use crate::spectral::{
    LAMBDA_MAX, LAMBDA_MIN,
    reflectance::{SpectralFractionTable, SpectralFractionTableError},
};

/// Bounded piecewise-linear emissivity on an explicit finite source band.
///
/// The accepted law retains 2–16 strictly increasing positive vacuum-wavelength
/// knots and finite emissivity in `[0,1]`. It never extrapolates,
/// interpolates RGB values or implies a law outside the declared band.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralEmissivity {
    table: SpectralFractionTable,
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
        Ok(Self {
            table: SpectralFractionTable::table(knots).map_err(thermal_table_error)?,
        })
    }

    /// Retain the original finite emissivity band, including UV or infrared.
    /// The first and last positive finite vacuum wavelengths define support.
    /// The same bounded interpolation owner serves absorption and reflection.
    pub fn table_in_band(knots: &[(f64, f64)]) -> Result<Self, ThermalSurfaceError> {
        Ok(Self {
            table: SpectralFractionTable::table_in_band(knots).map_err(thermal_table_error)?,
        })
    }

    /// Constant emissivity only within the explicitly declared vacuum-nm band.
    /// This is a finite-band law, with no implied gray continuation.
    pub fn constant_in_band(
        emissivity: f64,
        wavelength_bounds_nm: (f64, f64),
    ) -> Result<Self, ThermalSurfaceError> {
        Self::table_in_band(&[
            (wavelength_bounds_nm.0, emissivity),
            (wavelength_bounds_nm.1, emissivity),
        ])
    }

    /// Inclusive support in vacuum nanometres.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.table.wavelength_bounds_nm()
    }

    /// Evaluate the accepted linear law at a wavelength in the closed band.
    pub fn at(&self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        self.table.at(wavelength_nm).map_err(thermal_table_error)
    }

    /// Declared wavelength and emissivity knots, without unused inline slots.
    #[must_use]
    pub fn knots(&self) -> &[(f64, f64)] {
        self.table.knots()
    }

    /// Identity of the declared band, interpolation model and numerical table.
    /// Redundant knots retain distinct declaration identities.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("org.frankensim.render.spectral-emissivity.v1");
        hash.update(&[self.knots().len() as u8]);
        for &(wavelength_nm, emissivity) in self.knots() {
            hash.update(&wavelength_nm.to_bits().to_le_bytes());
            hash.update(&emissivity.to_bits().to_le_bytes());
        }
        hash.finalize()
    }
}

fn thermal_table_error(error: SpectralFractionTableError) -> ThermalSurfaceError {
    match error {
        SpectralFractionTableError::Knots => ThermalSurfaceError::SpectralEmissivity,
        SpectralFractionTableError::Fraction => ThermalSurfaceError::Emissivity,
        SpectralFractionTableError::Wavelength => ThermalSurfaceError::Wavelength,
    }
}

// Shared finite-band blackbody representability admission. This is a source
// envelope check, never a conversion between facet and macroscopic optical laws.
pub(super) fn admit_source_band_temperature(
    temperature_k: f64,
    minimum_nm: f64,
) -> Result<(), ThermalSurfaceError> {
    // exp(x)-1 >= x gives the Rayleigh–Jeans upper bound B_lambda <=
    // 2*c*k*T/lambda^4. Its largest value is at the source band's lower
    // endpoint, which need not be a visible wavelength.
    // This check avoids both a broadband T^4 restriction and a source
    // that overflows before the accepted emissivity can be applied.
    let lambda2 = minimum_nm * minimum_nm;
    let upper_per_kelvin = 2.0 * LIGHT_M_S * BOLTZMANN_J_K * 1.0e27 / (lambda2 * lambda2);
    let unrepresentable = if upper_per_kelvin.is_finite() && upper_per_kelvin > 0.0 {
        temperature_k > f64::MAX / upper_per_kelvin
    } else {
        // Retain the same bound when a declared wavelength makes an
        // intermediate fourth power overflow or underflow.
        fs_math::det::ln(2.0 * LIGHT_M_S * BOLTZMANN_J_K * 1.0e27) + fs_math::det::ln(temperature_k)
            - 4.0 * fs_math::det::ln(minimum_nm)
            > fs_math::det::ln(f64::MAX)
    };
    if unrepresentable {
        return Err(ThermalSurfaceError::Unrepresentable);
    }
    planck_radiance_w_m2_sr_nm(temperature_k, minimum_nm)?;
    Ok(())
}

/// Accepted source-band opaque diffuse surface in local thermal equilibrium.
///
/// The source identity binds this declared optical law and temperature to the
/// caller's accepted physical source; supplying a hash does not create measured
/// provenance. A hemispherical-total emissivity alone is insufficient to infer
/// the spectral knots. This view emits into vacuum/the exterior-air model and
/// keeps the existing two-sided, per-face convention of gray thermal surfaces.
///
/// No bolometric exitance or cooling API is provided: the unmeasured remainder
/// and its owning thermal model must be supplied before asserting total power.
/// Observing this spectrum never advances or debits the thermal owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralDiffuseThermalSurface {
    temperature_k: f64,
    emissivity: SpectralEmissivity,
    source_identity: ContentHash,
}

impl SpectralDiffuseThermalSurface {
    /// Bind a finite positive temperature to the accepted finite source law.
    ///
    /// An emitting surface requires a representable source-band blackbody
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
        if emissivity.table.has_positive_value() {
            admit_source_band_temperature(temperature_k, emissivity.wavelength_bounds_nm().0)?;
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

    /// Accepted source-band spectral emissivity law.
    #[must_use]
    pub const fn emissivity(&self) -> SpectralEmissivity {
        self.emissivity
    }

    /// Borrow the declared knots for source integration without copying them.
    pub(crate) fn emissivity_knots(&self) -> &[(f64, f64)] {
        self.emissivity.knots()
    }

    /// Reciprocal opaque diffuse reflectance `1 - emissivity(lambda)`.
    pub fn reflectance(&self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        Ok(1.0 - self.emissivity.at(wavelength_nm)?)
    }

    /// Emitted spectral radiance in W/(m² sr nm), inside the declared band.
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

    /// Whether the accepted emissivity has positive support in its source band.
    /// A sufficiently cold spectrum may still underflow numerically to zero.
    #[must_use]
    pub const fn is_emitting(&self) -> bool {
        self.emissivity.table.has_positive_value()
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
/// and rough laws stay borrowed, so light triangles and samples do not copy
/// knot tables. Rough source queries require the actual outgoing direction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DiffuseThermalSurface<'a> {
    /// Existing gray source, with its original all-wavelength behavior.
    Gray(GrayDiffuseThermalSurface),
    /// Accepted finite-band spectral source, borrowed from its material owner.
    Spectral(&'a SpectralDiffuseThermalSurface),
    /// Original EON facets; actual outgoing angle is required for radiance.
    Rough(&'a super::RoughDiffuseThermalSurface),
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

impl<'a> From<&'a super::RoughDiffuseThermalSurface> for DiffuseThermalSurface<'a> {
    fn from(surface: &'a super::RoughDiffuseThermalSurface) -> Self {
        Self::Rough(surface)
    }
}

impl DiffuseThermalSurface<'_> {
    /// Accepted thermodynamic temperature of the original source, in kelvin.
    #[must_use]
    pub const fn temperature_k(self) -> f64 {
        match self {
            Self::Gray(surface) => surface.temperature_k(),
            Self::Spectral(surface) => surface.temperature_k(),
            Self::Rough(surface) => surface.temperature_k(),
        }
    }

    /// Opaque reciprocal reflectance within the source's wavelength domain.
    /// Rough sources refuse: facet albedo is not macroscopic reflectance, and
    /// their reflected path response belongs to the original EON BRDF.
    pub fn reflectance(self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        match self {
            Self::Gray(surface) => {
                if !wavelength_nm.is_finite() || wavelength_nm <= 0.0 {
                    return Err(ThermalSurfaceError::Wavelength);
                }
                Ok(surface.reflectance())
            }
            Self::Spectral(surface) => surface.reflectance(wavelength_nm),
            Self::Rough(_) => Err(ThermalSurfaceError::DirectionRequired),
        }
    }

    /// Emitted radiance in W/(m² sr nm) from the same accepted source law.
    /// Rough sources require `directional_radiance_w_m2_sr_nm` instead.
    pub fn radiance_w_m2_sr_nm(self, wavelength_nm: f64) -> Result<f64, ThermalSurfaceError> {
        match self {
            Self::Gray(surface) => surface.radiance_w_m2_sr_nm(wavelength_nm),
            Self::Spectral(surface) => surface.radiance_w_m2_sr_nm(wavelength_nm),
            Self::Rough(_) => Err(ThermalSurfaceError::DirectionRequired),
        }
    }

    /// Evaluate the original emission law in its actual outgoing hemisphere.
    /// The cosine is relative to the physical face, not a sampled microfacet.
    pub fn directional_radiance_w_m2_sr_nm(
        self,
        wavelength_nm: f64,
        outgoing_cosine: f64,
    ) -> Result<f64, ThermalSurfaceError> {
        if !outgoing_cosine.is_finite() || !(0.0..=1.0).contains(&outgoing_cosine) {
            return Err(ThermalSurfaceError::Direction);
        }
        match self {
            Self::Rough(surface) => {
                surface.directional_radiance_w_m2_sr_nm(wavelength_nm, outgoing_cosine)
            }
            _ => self.radiance_w_m2_sr_nm(wavelength_nm),
        }
    }

    /// Projected-hemisphere mean used for source proposals, never as a
    /// direction-independent replacement for a rough source's path radiance.
    pub fn hemispherical_mean_radiance_w_m2_sr_nm(
        self,
        wavelength_nm: f64,
    ) -> Result<f64, ThermalSurfaceError> {
        match self {
            Self::Rough(surface) => surface.hemispherical_mean_radiance_w_m2_sr_nm(wavelength_nm),
            _ => self.radiance_w_m2_sr_nm(wavelength_nm),
        }
    }

    /// Retain the original constitutive source identity without RGB conversion.
    #[must_use]
    pub fn identity(self) -> ContentHash {
        match self {
            Self::Gray(surface) => surface.identity(),
            Self::Spectral(surface) => surface.identity(),
            Self::Rough(surface) => surface.identity(),
        }
    }

    /// Whether the accepted source has any positive emissivity support.
    #[must_use]
    pub fn is_emitting(self) -> bool {
        match self {
            Self::Gray(surface) => surface.emissivity() > 0.0,
            Self::Spectral(surface) => surface.is_emitting(),
            Self::Rough(surface) => surface.is_emitting(),
        }
    }
}

#[cfg(test)]
mod tests;
