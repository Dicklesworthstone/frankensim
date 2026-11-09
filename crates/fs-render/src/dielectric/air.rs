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

use super::{DielectricError, DielectricIor, SellmeierIor, admitted_dielectric_index};
use fs_blake3::{ContentHash, DomainHasher};
use fs_material::air_refraction::{AirRefractiveSample, CiddorAir, WAVELENGTH_RANGE_NM};

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
    // Reuse the validated oscillator arithmetic in the explicit source
    // coordinate; this private value is not exposed as an absolute glass law.
    relative: SellmeierIor,
    reference_air: CiddorAirIor,
    wavelength_coordinate: AirRelativeWavelength,
    vacuum_bounds_nm: (f64, f64),
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
        let relative =
            SellmeierIor::try_new_in_band(strengths, resonances_um2, source_wavelength_bounds_nm)?;
        let mut candidate = Self {
            relative,
            reference_air,
            wavelength_coordinate,
            vacuum_bounds_nm: reference_air.wavelength_bounds_nm(),
        };
        let (source_lo, source_hi) = source_wavelength_bounds_nm;
        let (air_lo, air_hi) = reference_air.wavelength_bounds_nm();
        let (lo, hi) = match wavelength_coordinate {
            AirRelativeWavelength::Vacuum => (source_lo.max(air_lo), source_hi.min(air_hi)),
            AirRelativeWavelength::Air(_) => {
                let mapped_lo = candidate.source_wavelength(air_lo)?.0;
                let mapped_hi = candidate.source_wavelength(air_hi)?.0;
                if source_hi < mapped_lo || source_lo > mapped_hi {
                    return Err(DielectricError::InvalidWavelength);
                }
                (
                    if source_lo <= mapped_lo {
                        air_lo
                    } else {
                        candidate.inverse_air_wavelength(source_lo, true)?
                    },
                    if source_hi >= mapped_hi {
                        air_hi
                    } else {
                        candidate.inverse_air_wavelength(source_hi, false)?
                    },
                )
            }
        };
        if lo > hi {
            return Err(DielectricError::InvalidWavelength);
        }
        candidate.vacuum_bounds_nm = (lo, hi);
        // Both factors are positive and nonincreasing on the admitted bands;
        // lambda_vac/n_air is increasing. Endpoint admission covers this product.
        candidate.eval(lo)?;
        candidate.eval(hi)?;
        Ok(candidate)
    }

    pub const fn coefficients(&self) -> ([f64; 3], [f64; 3]) {
        self.relative.coefficients()
    }
    pub const fn reference_air(&self) -> &CiddorAirIor {
        &self.reference_air
    }
    pub const fn wavelength_coordinate(&self) -> AirRelativeWavelength {
        self.wavelength_coordinate
    }
    /// Source support in the explicitly declared fit coordinate.
    pub const fn source_wavelength_bounds_nm(&self) -> (f64, f64) {
        self.relative.wavelength_bounds_nm()
    }
    /// Effective absolute-vacuum support consumed by all native optical paths.
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.vacuum_bounds_nm
    }
    pub fn eval(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        let (coordinate, air, _) = self.inputs(wavelength_nm)?;
        let index = self.relative.eval(coordinate)? * air.phase_index;
        if !admitted_dielectric_index(index) {
            return Err(DielectricError::InvalidIor);
        }
        Ok(index)
    }
    pub fn group_index(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        let (coordinate, air, coordinate_air) = self.inputs(wavelength_nm)?;
        let relative_group = self.relative.group_index(coordinate)?;
        let group = match self.wavelength_coordinate {
            AirRelativeWavelength::Vacuum => {
                let relative_phase = self.relative.eval(coordinate)?;
                relative_group * air.phase_index
                    + relative_phase * (air.group_index - air.phase_index)
            }
            AirRelativeWavelength::Air(law) if law == self.reference_air => {
                relative_group * air.group_index
            }
            AirRelativeWavelength::Air(_) => {
                let relative_phase = self.relative.eval(coordinate)?;
                let coordinate_air = coordinate_air.expect("air wavelength sample");
                air.group_index * relative_phase
                    + (air.phase_index / coordinate_air.phase_index)
                        * coordinate_air.group_index
                        * (relative_group - relative_phase)
            }
        };
        if !group.is_finite() || group <= 0.0 {
            return Err(DielectricError::InvalidIor);
        }
        Ok(group)
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

    fn admit_wavelength(&self, wavelength_nm: f64) -> Result<(), DielectricError> {
        let (lo, hi) = self.vacuum_bounds_nm;
        if !wavelength_nm.is_finite() || !(lo..=hi).contains(&wavelength_nm) {
            return Err(DielectricError::InvalidWavelength);
        }
        Ok(())
    }

    fn inputs(
        &self,
        vacuum_nm: f64,
    ) -> Result<(f64, AirRefractiveSample, Option<AirRefractiveSample>), DielectricError> {
        self.admit_wavelength(vacuum_nm)?;
        let (coordinate, coordinate_air) = self.source_wavelength(vacuum_nm)?;
        let index_air = match self.wavelength_coordinate {
            AirRelativeWavelength::Air(law) if law == self.reference_air => {
                coordinate_air.expect("air wavelength sample")
            }
            _ => self.reference_air.sample(vacuum_nm)?,
        };
        Ok((coordinate, index_air, coordinate_air))
    }

    fn source_wavelength(
        &self,
        vacuum_nm: f64,
    ) -> Result<(f64, Option<AirRefractiveSample>), DielectricError> {
        match self.wavelength_coordinate {
            AirRelativeWavelength::Vacuum => Ok((vacuum_nm, None)),
            AirRelativeWavelength::Air(law) => {
                let air = law.sample(vacuum_nm)?;
                Ok((vacuum_nm / air.phase_index, Some(air)))
            }
        }
    }

    // n_air is positive with normal dispersion in the source domain, hence
    // d(lambda_vac/n_air)/d(lambda_vac)=ng_air/n_air² > 0. A bounded bisection
    // finds the inward representable band edge once, never during ray queries.
    fn inverse_air_wavelength(&self, source_nm: f64, lower: bool) -> Result<f64, DielectricError> {
        let (mut lo, mut hi) = self.reference_air.wavelength_bounds_nm();
        for _ in 0..64 {
            if lo.next_up() >= hi {
                break;
            }
            let middle = lo + 0.5 * (hi - lo);
            let mapped = self.source_wavelength(middle)?.0;
            if if lower {
                mapped < source_nm
            } else {
                mapped <= source_nm
            } {
                lo = middle;
            } else {
                hi = middle;
            }
        }
        Ok(if lower { hi } else { lo })
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
