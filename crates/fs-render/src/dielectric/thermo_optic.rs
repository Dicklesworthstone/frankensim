//! Frozen-temperature glass using the original material thermo-optic law.
//!
//! The absolute reference dispersion and SCHOTT TIE-19 increment remain
//! separate, immutable sources. This composition feeds the same phase/group
//! interface as ordinary lenses, nested dielectrics, films and curved rays.
//! It never fits the heated index to new Sellmeier or tabulated coefficients.

use super::{
    AirRelativeSellmeierIor, CauchyIor, DielectricError, DielectricIor, SellmeierIor, TabulatedIor,
    admitted_dielectric_index,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_material::thermo_optic::{SchottThermoOpticLaw, ThermoOpticError};

// A nonrecursive reference retains the native inline, allocation-free law.
// Atmospheric/water states and an already heated index are not glass fits at
// the declared reference temperature, so they cannot silently receive this
// glass correction. Air-relative glass is converted by its original owner.
#[derive(Clone, Copy, Debug, PartialEq)]
enum ReferenceIor {
    Cauchy(CauchyIor),
    Sellmeier(SellmeierIor),
    Tabulated(TabulatedIor),
    AirRelativeSellmeier(AirRelativeSellmeierIor),
}

impl ReferenceIor {
    fn try_new(ior: DielectricIor) -> Result<Self, DielectricError> {
        match ior {
            DielectricIor::Cauchy(law) => Ok(Self::Cauchy(law)),
            DielectricIor::Sellmeier(law) => Ok(Self::Sellmeier(law)),
            DielectricIor::Tabulated(law) => Ok(Self::Tabulated(law)),
            DielectricIor::AirRelativeSellmeier(law) => Ok(Self::AirRelativeSellmeier(law)),
            _ => Err(DielectricError::InvalidIor),
        }
    }

    const fn as_ior(self) -> DielectricIor {
        match self {
            Self::Cauchy(law) => DielectricIor::Cauchy(law),
            Self::Sellmeier(law) => DielectricIor::Sellmeier(law),
            Self::Tabulated(law) => DielectricIor::Tabulated(law),
            Self::AirRelativeSellmeier(law) => DielectricIor::AirRelativeSellmeier(law),
        }
    }

    fn phase(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        match self {
            Self::Cauchy(law) => law.eval(wavelength_nm),
            Self::Sellmeier(law) => law.eval(wavelength_nm),
            Self::Tabulated(law) => law.eval(wavelength_nm),
            Self::AirRelativeSellmeier(law) => law.eval(wavelength_nm),
        }
    }

    fn group(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        match self {
            Self::Cauchy(law) => law.group_index(wavelength_nm),
            Self::Sellmeier(law) => law.group_index(wavelength_nm),
            Self::Tabulated(law) => law.group_index(wavelength_nm),
            Self::AirRelativeSellmeier(law) => law.group_index(wavelength_nm),
        }
    }

    fn update_identity(&self, hash: &mut DomainHasher) {
        match self {
            Self::Cauchy(law) => {
                hash.update(&[0]);
                for value in law.coefficients() {
                    hash.update(&value.to_bits().to_le_bytes());
                }
            }
            Self::Sellmeier(law) => {
                hash.update(&[1]);
                let (strengths, resonances) = law.coefficients();
                for value in strengths.into_iter().chain(resonances) {
                    hash.update(&value.to_bits().to_le_bytes());
                }
                law.update_identity_band(hash);
            }
            Self::Tabulated(law) => {
                hash.update(&[2]);
                hash.update(&law.max_knot_span_nm().to_bits().to_le_bytes());
                for knot in law.samples() {
                    hash.update(&knot.wavelength_nm().to_bits().to_le_bytes());
                    hash.update(&knot.index().to_bits().to_le_bytes());
                }
            }
            Self::AirRelativeSellmeier(law) => {
                hash.update(&[3]);
                hash.update(law.content_identity().as_bytes());
            }
        }
    }
}

/// Absolute phase/group dispersion of one homogeneous glass at a frozen T.
///
/// The reference must be an original Cauchy, Sellmeier, tabulated, or explicitly
/// air-corrected Sellmeier glass at the material law's reference temperature.
/// Temperature and spectral support come from that same law. Queries retain
/// its analytic group derivative, including the wavelength dependence of the
/// reference-index prefactor. The reference atmosphere is never heated along
/// with the glass. A refresh reevaluates the original reference, not the last
/// accepted index, so temperature cycles cannot accumulate optical drift.
///
/// Absorption, expansion, stress and source authentication are separate. The
/// constructor checks endpoints; every query checks its actual response. This
/// is not an interval certificate for arbitrary supplied coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SchottThermoOpticIor {
    reference: ReferenceIor,
    law: SchottThermoOpticLaw,
    temperature_k: f64,
    wavelength_bounds_nm: (f64, f64),
    identity: ContentHash,
}

impl SchottThermoOpticIor {
    /// Bind an original absolute-reference glass and explicitly supported T.
    /// An already thermally shifted law refuses, preventing double correction.
    pub fn try_new(
        reference: impl Into<DielectricIor>,
        law: SchottThermoOpticLaw,
        temperature_k: f64,
    ) -> Result<Self, DielectricError> {
        let reference = ReferenceIor::try_new(reference.into())?;
        let [min_t, max_t] = law.temperature_range_k();
        if !temperature_k.is_finite() || !(min_t..=max_t).contains(&temperature_k) {
            return Err(DielectricError::InvalidIor);
        }
        let (min_r, max_r) = reference.as_ior().wavelength_bounds_nm();
        let [min_l, max_l] = law.wavelength_range_nm();
        let band = (min_r.max(min_l), max_r.min(max_l));
        if band.0 > band.1 {
            return Err(DielectricError::InvalidWavelength);
        }
        let mut hash = DomainHasher::new("fs-render:schott-thermo-optic-ior:v1");
        reference.update_identity(&mut hash);
        hash.update(law.identity().as_bytes());
        hash.update(&temperature_k.to_bits().to_le_bytes());
        let result = Self {
            reference,
            law,
            temperature_k,
            wavelength_bounds_nm: band,
            identity: hash.finalize(),
        };
        result.eval(band.0)?;
        result.eval(band.1)?;
        Ok(result)
    }

    /// Exact original reference; no fitted heated coefficients are substituted.
    pub const fn reference_ior(&self) -> DielectricIor {
        self.reference.as_ior()
    }
    pub const fn law(&self) -> &SchottThermoOpticLaw {
        &self.law
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.wavelength_bounds_nm
    }
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Rebind from the immutable reference, preserving coefficients and source.
    pub fn at_temperature(&self, temperature_k: f64) -> Result<Self, DielectricError> {
        Self::try_new(self.reference_ior(), self.law, temperature_k)
    }

    /// Phase-only transport does not require a derivative of a tabulated knot.
    pub fn eval(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        self.admit_wavelength(wavelength_nm)?;
        let n0 = self.reference.phase(wavelength_nm)?;
        let sample = self
            .law
            .prepare(wavelength_nm, n0, None)
            .and_then(|law| law.sample(self.temperature_k))
            .map_err(map_error)?;
        if !admitted_dielectric_index(sample.phase_index) {
            return Err(DielectricError::InvalidIor);
        }
        Ok(sample.phase_index)
    }

    /// Original analytic `n - lambda dn/dlambda`, without finite differences.
    pub fn group_index(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        self.admit_wavelength(wavelength_nm)?;
        let n0 = self.reference.phase(wavelength_nm)?;
        let ng0 = self.reference.group(wavelength_nm)?;
        // A valid phase interpolant can lack an admitted positive group law.
        // Preserve phase-only transport just as the original table owner does.
        if !ng0.is_finite() || ng0 <= 0.0 {
            self.eval(wavelength_nm)?;
            return Err(DielectricError::UndefinedGroupIndex);
        }
        let sample = self
            .law
            .prepare(wavelength_nm, n0, Some(ng0))
            .and_then(|law| law.sample(self.temperature_k))
            .or_else(|error| {
                // Distinguish a failed phase law from unavailable thermal
                // group propagation; do not suppress genuine phase errors.
                self.eval(wavelength_nm)?;
                match error {
                    ThermoOpticError::NonPositiveIndex | ThermoOpticError::NumericalRange => {
                        Err(DielectricError::UndefinedGroupIndex)
                    }
                    _ => Err(map_error(error)),
                }
            })?;
        if !admitted_dielectric_index(sample.phase_index) {
            return Err(DielectricError::InvalidIor);
        }
        sample
            .group_index
            .ok_or(DielectricError::UndefinedGroupIndex)
    }

    fn admit_wavelength(&self, wavelength_nm: f64) -> Result<(), DielectricError> {
        let (minimum, maximum) = self.wavelength_bounds_nm;
        if !wavelength_nm.is_finite() || !(minimum..=maximum).contains(&wavelength_nm) {
            return Err(DielectricError::InvalidWavelength);
        }
        Ok(())
    }
}

fn map_error(error: ThermoOpticError) -> DielectricError {
    match error {
        ThermoOpticError::WavelengthOutsideRange => DielectricError::InvalidWavelength,
        _ => DielectricError::InvalidIor,
    }
}

impl From<SchottThermoOpticIor> for DielectricIor {
    fn from(law: SchottThermoOpticIor) -> Self {
        Self::ThermoOptic(law)
    }
}

#[cfg(test)]
mod tests;
