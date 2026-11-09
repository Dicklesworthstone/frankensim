//! Absolute atmospheric and air-relative catalog indices on native optical paths.
//!
//! The original NIST Ciddor material owner supplies the reference atmosphere.
//! Relative glass coefficients are never substituted for an absolute index or
//! refitted to a Cauchy/table approximation. The caller declares both the
//! reference-air state and the source's wavelength coordinate. In particular,
//! SCHOTT TIE-29 specifies relative indices and the spectral lines of its
//! Table 1; a source's air-wavelength convention cannot silently become vacuum.
//! <https://media.schott.com/api/public/content/aaa572afd854434fb7b3faa4bc46103f?v=7f40e8ec>
//!
//! NIST's coordinate relation is `lambda_air=lambda_vac/n_air`:
//! <https://emtoolbox.nist.gov/Wavelength/Documentation.asp>.
//! For a vacuum-parameterized relative fit, differentiation of `n_abs=n_air*n_rel`
//! gives `ng_abs=n_air*ng_rel+n_rel*(ng_air-n_air)`. A fit in air wavelength
//! has a separately declared coordinate atmosphere: its standard-air spectral
//! lines need not use the measurement's index-reference temperature. Writing
//! `A=n_index_air`, `B=n_wavelength_air`, the chain rule gives
//! `ng_abs=ng_A*n_rel+(A*ng_B/B)*(ng_rel-n_rel)` at `lambda_vac/B`.
//! Equal atmospheres simplify to `ng_abs=ng_air*ng_rel`. These are narrowband
//! group indices, not a pulse-reshaping model or a glass temperature correction.

use super::{DielectricError, DielectricIor};
use fs_blake3::{ContentHash, DomainHasher};
use fs_material::air_refraction::{AirRefractiveSample, CiddorAir, WAVELENGTH_RANGE_NM};
use fs_material::catalog_refraction::{
    AirRelativeSellmeier, CatalogRefractionError, CatalogWavelength, CiddorAirState,
};

/// Homogeneous atmospheric phase/group law at an explicitly frozen state.
/// Pressure is absolute Pa, temperature K, and water content is mole fraction
/// in the borrowed composition model. Absorption is independently authored in
/// `DielectricGlass`; this adapter never infers atmospheric transmission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CiddorAirIor {
    composition: CiddorAir,
    temperature_k: f64,
    pressure_pa: f64,
}

impl CiddorAirIor {
    pub fn try_new(
        composition: CiddorAir,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<Self, DielectricError> {
        let candidate = Self {
            composition,
            temperature_k,
            pressure_pa,
        };
        candidate.sample(WAVELENGTH_RANGE_NM[0])?;
        candidate.sample(WAVELENGTH_RANGE_NM[1])?;
        Ok(candidate)
    }

    pub const fn composition(&self) -> CiddorAir {
        self.composition
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn pressure_pa(&self) -> f64 {
        self.pressure_pa
    }
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        (WAVELENGTH_RANGE_NM[0], WAVELENGTH_RANGE_NM[1])
    }
    pub fn eval(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        Ok(self.sample(wavelength_nm)?.phase_index)
    }
    pub fn group_index(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        Ok(self.sample(wavelength_nm)?.group_index)
    }
    /// Complete source/formulation, composition and frozen-state identity.
    pub fn content_identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("fs-render:ciddor-air-ior:v1");
        hash.update(self.composition.identity().as_bytes());
        hash.update(&self.temperature_k.to_bits().to_le_bytes());
        hash.update(&self.pressure_pa.to_bits().to_le_bytes());
        hash.finalize()
    }

    fn sample(&self, wavelength_nm: f64) -> Result<AirRefractiveSample, DielectricError> {
        self.composition
            .prepare(wavelength_nm)
            .map_err(|_| DielectricError::InvalidWavelength)?
            .sample(self.temperature_k, self.pressure_pa)
            .map_err(|_| DielectricError::InvalidIor)
    }
}

/// Coordinate used to fit the original relative-index oscillator coefficients.
/// Air wavelength uses its own explicit atmosphere: standard-air spectral
/// lines can have a different temperature/pressure from the index measurement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AirRelativeWavelength {
    Vacuum,
    Air(CiddorAirIor),
}

/// Absolute glass law obtained from original `n_glass/n_reference_air` data.
///
/// Reference temperature, pressure and composition belong to the measurement
/// convention. Changing the scene atmosphere must not change this reference;
/// instead use a separate `CiddorAirIor` for the actual exterior medium. Source
/// authentication and glass thermo-optic evolution remain separate concerns.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AirRelativeSellmeierIor {
    // The lower material owner is shared with native electromagnetic consumers.
    catalog: AirRelativeSellmeier,
    reference_air: CiddorAirIor,
    wavelength_coordinate: AirRelativeWavelength,
}

impl AirRelativeSellmeierIor {
    /// Retain the original strengths, squared resonance wavelengths in µm²,
    /// and inclusive source-coordinate support in nm. The native Sellmeier
    /// pole/ordinary-index admission applies to the relative function.
    /// Effective vacuum support is its exact coordinate-mapped intersection
    /// with the Ciddor band; endpoints are rounded inward without extrapolation.
    pub fn try_new(
        strengths: [f64; 3],
        resonances_um2: [f64; 3],
        source_wavelength_bounds_nm: (f64, f64),
        reference_air: CiddorAirIor,
        wavelength_coordinate: AirRelativeWavelength,
    ) -> Result<Self, DielectricError> {
        let coordinate = match wavelength_coordinate {
            AirRelativeWavelength::Vacuum => CatalogWavelength::Vacuum,
            AirRelativeWavelength::Air(air) => CatalogWavelength::Air(catalog_air(air)?),
        };
        let catalog = AirRelativeSellmeier::try_new(
            strengths,
            resonances_um2,
            source_wavelength_bounds_nm,
            catalog_air(reference_air)?,
            coordinate,
        )
        .map_err(catalog_error)?;
        Ok(Self {
            catalog,
            reference_air,
            wavelength_coordinate,
        })
    }

    pub const fn coefficients(&self) -> ([f64; 3], [f64; 3]) {
        self.catalog.coefficients()
    }
    pub const fn reference_air(&self) -> &CiddorAirIor {
        &self.reference_air
    }
    pub const fn wavelength_coordinate(&self) -> AirRelativeWavelength {
        self.wavelength_coordinate
    }
    /// Source support in the explicitly declared fit coordinate.
    pub const fn source_wavelength_bounds_nm(&self) -> (f64, f64) {
        self.catalog.source_wavelength_bounds_nm()
    }
    /// Effective absolute-vacuum support consumed by all native optical paths.
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.catalog.wavelength_bounds_nm()
    }
    pub fn eval(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        self.catalog.eval(wavelength_nm).map_err(catalog_error)
    }
    pub fn group_index(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        self.catalog
            .group_index(wavelength_nm)
            .map_err(catalog_error)
    }
    /// The exact shared material law; its derivatives retain both atmospheres.
    pub const fn material_law(&self) -> &AirRelativeSellmeier {
        &self.catalog
    }
    pub fn content_identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("fs-render:air-relative-sellmeier:v1");
        hash.update(self.reference_air.content_identity().as_bytes());
        match self.wavelength_coordinate {
            AirRelativeWavelength::Vacuum => hash.update(&[0]),
            AirRelativeWavelength::Air(law) => {
                hash.update(&[1]);
                hash.update(law.content_identity().as_bytes());
            }
        }
        let (strengths, resonances) = self.coefficients();
        let (lo, hi) = self.source_wavelength_bounds_nm();
        for value in strengths.into_iter().chain(resonances).chain([lo, hi]) {
            hash.update(&value.to_bits().to_le_bytes());
        }
        hash.finalize()
    }
}

fn catalog_air(air: CiddorAirIor) -> Result<CiddorAirState, DielectricError> {
    CiddorAirState::try_new(air.composition, air.temperature_k, air.pressure_pa)
        .map_err(|_| DielectricError::InvalidIor)
}

fn catalog_error(error: CatalogRefractionError) -> DielectricError {
    match error {
        CatalogRefractionError::InvalidWavelength => DielectricError::InvalidWavelength,
        _ => DielectricError::InvalidIor,
    }
}

impl From<CiddorAirIor> for DielectricIor {
    fn from(ior: CiddorAirIor) -> Self {
        Self::Air(ior)
    }
}
impl From<AirRelativeSellmeierIor> for DielectricIor {
    fn from(ior: AirRelativeSellmeierIor) -> Self {
        Self::AirRelativeSellmeier(ior)
    }
}

#[cfg(test)]
mod tests;
