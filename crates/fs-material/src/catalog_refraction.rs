//! Absolute optical glass from an explicitly air-relative Sellmeier catalog law.
//!
//! The index-reference atmosphere and wavelength-coordinate atmosphere are
//! independent original measurement conventions. Neither is the scene medium
//! or a glass temperature correction. No coefficient fitting or absorption is
//! inferred. SCHOTT TIE-29 describes relative indices and spectral-line data:
//! <https://media.schott.com/api/public/content/aaa572afd854434fb7b3faa4bc46103f?v=7f40e8ec>.
//! The original Ciddor owner supplies air and `lambda_air=lambda_vac/n_air`:
//! <https://emtoolbox.nist.gov/Wavelength/Documentation.asp>.
//!
//! For `N(lambda)=A(lambda) F(lambda/B(lambda))`, the wavelength jet applies
//! the full product and coordinate chain rules, including ordinary curvature.
//! Preparation maps source-band edges inward once; additional catalog curvature
//! and product/coordinate second derivatives are evaluated only when requested.
//! All storage and query work are bounded, with
//! no allocation. Finite-band derivatives do not certify causal pulse transport
//! or measured source accuracy. No unsafe code or hidden material state is used.

use crate::air_refraction::{
    AirRefractionError, AirRefractiveSample, CiddorAir, WAVELENGTH_RANGE_NM,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_math::{dd::Dd, det};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CatalogRefractionError {
    InvalidWavelength,
    InvalidIndex,
    Air(AirRefractionError),
}
impl core::fmt::Display for CatalogRefractionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "air-relative catalog refraction: {self:?}")
    }
}
impl core::error::Error for CatalogRefractionError {}

/// One immutable original atmosphere, independent of the actual glass state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CiddorAirState {
    composition: CiddorAir,
    temperature_k: f64,
    pressure_pa: f64,
}
impl CiddorAirState {
    pub fn try_new(
        composition: CiddorAir,
        temperature_k: f64,
        pressure_pa: f64,
    ) -> Result<Self, AirRefractionError> {
        let value = Self {
            composition,
            temperature_k,
            pressure_pa,
        };
        value.sample(WAVELENGTH_RANGE_NM[0])?;
        value.sample(WAVELENGTH_RANGE_NM[1])?;
        Ok(value)
    }
    pub const fn composition(self) -> CiddorAir {
        self.composition
    }
    pub const fn temperature_k(self) -> f64 {
        self.temperature_k
    }
    pub const fn pressure_pa(self) -> f64 {
        self.pressure_pa
    }
    pub fn sample(self, wavelength_nm: f64) -> Result<AirRefractiveSample, AirRefractionError> {
        self.composition
            .prepare(wavelength_nm)?
            .sample(self.temperature_k, self.pressure_pa)
    }
    fn hash_into(self, hash: &mut DomainHasher) {
        hash.update(self.composition.identity().as_bytes());
        hash.update(&self.temperature_k.to_bits().to_le_bytes());
        hash.update(&self.pressure_pa.to_bits().to_le_bytes());
    }
}

/// Coordinate of the original fit; its atmosphere need not be the index reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CatalogWavelength {
    Vacuum,
    Air(CiddorAirState),
}

/// Absolute vacuum index and its derivatives at fixed original source conventions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CatalogRefractiveIndexJet {
    pub phase_index: f64,
    pub group_index: f64,
    /// `[lambda dN/dlambda, lambda² d²N/dlambda²]`, ordinary scaled derivatives.
    pub scaled_wavelength_derivatives: [f64; 2],
}

/// Original `n_glass/n_reference_air` oscillators with an explicit fit coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AirRelativeSellmeier {
    strengths: [f64; 3],
    resonances_um2: [f64; 3],
    source_bounds_nm: (f64, f64),
    reference_air: CiddorAirState,
    coordinate: CatalogWavelength,
    vacuum_bounds_nm: (f64, f64),
}

impl AirRelativeSellmeier {
    /// Strengths are dimensionless, resonances are squared micrometres, and
    /// the inclusive fit band is in nanometres of the declared coordinate.
    /// Admit the ENTIRE original fit before intersecting its mapped air band.
    pub fn try_new(
        strengths: [f64; 3],
        resonances_um2: [f64; 3],
        source_bounds_nm: (f64, f64),
        reference_air: CiddorAirState,
        coordinate: CatalogWavelength,
    ) -> Result<Self, CatalogRefractionError> {
        let (source_lo, source_hi) = source_bounds_nm;
        let min_square = (source_lo * 1e-3) * (source_lo * 1e-3);
        let max_square = (source_hi * 1e-3) * (source_hi * 1e-3);
        if !source_lo.is_finite()
            || !source_hi.is_finite()
            || source_lo <= 0.0
            || source_lo > source_hi
            || min_square <= 0.0
            || !max_square.is_finite()
        {
            return Err(CatalogRefractionError::InvalidWavelength);
        }
        let mut value = Self {
            strengths,
            resonances_um2,
            source_bounds_nm,
            reference_air,
            coordinate,
            vacuum_bounds_nm: (WAVELENGTH_RANGE_NM[0], WAVELENGTH_RANGE_NM[1]),
        };
        for i in 0..3 {
            let (b, c) = (strengths[i], resonances_um2[i]);
            if !b.is_finite()
                || b < 0.0
                || !c.is_finite()
                || c < 0.0
                || (b > 0.0 && (min_square..=max_square).contains(&c))
            {
                return Err(CatalogRefractionError::InvalidIndex);
            }
            if b == 0.0 {
                value.strengths[i] = 0.0;
                value.resonances_um2[i] = 0.0;
            } else if c == 0.0 {
                value.resonances_um2[i] = 0.0;
            }
        }
        value.relative_phase(source_lo)?;
        value.relative_phase(source_hi)?;
        let [air_lo, air_hi] = WAVELENGTH_RANGE_NM;
        let (lo, hi) = match coordinate {
            CatalogWavelength::Vacuum => (source_lo.max(air_lo), source_hi.min(air_hi)),
            CatalogWavelength::Air(_) => {
                let mapped_lo = value.source_wavelength(air_lo)?.0;
                let mapped_hi = value.source_wavelength(air_hi)?.0;
                if source_hi < mapped_lo || source_lo > mapped_hi {
                    return Err(CatalogRefractionError::InvalidWavelength);
                }
                (
                    if source_lo <= mapped_lo {
                        air_lo
                    } else {
                        value.inverse_air_wavelength(source_lo, true)?
                    },
                    if source_hi >= mapped_hi {
                        air_hi
                    } else {
                        value.inverse_air_wavelength(source_hi, false)?
                    },
                )
            }
        };
        if lo > hi {
            return Err(CatalogRefractionError::InvalidWavelength);
        }
        value.vacuum_bounds_nm = (lo, hi);
        value.eval(lo)?;
        value.eval(hi)?;
        Ok(value)
    }

    pub const fn coefficients(&self) -> ([f64; 3], [f64; 3]) {
        (self.strengths, self.resonances_um2)
    }
    pub const fn reference_air(&self) -> CiddorAirState {
        self.reference_air
    }
    pub const fn wavelength_coordinate(&self) -> CatalogWavelength {
        self.coordinate
    }
    pub const fn source_wavelength_bounds_nm(&self) -> (f64, f64) {
        self.source_bounds_nm
    }
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.vacuum_bounds_nm
    }

    pub fn eval(&self, wavelength_nm: f64) -> Result<f64, CatalogRefractionError> {
        let (coordinate, air, _) = self.inputs(wavelength_nm)?;
        let index = self.relative_phase(coordinate)? * air.phase_index;
        if !admitted_index(index) {
            return Err(CatalogRefractionError::InvalidIndex);
        }
        Ok(index)
    }

    pub fn group_index(&self, wavelength_nm: f64) -> Result<f64, CatalogRefractionError> {
        let (coordinate, air, coordinate_air) = self.inputs(wavelength_nm)?;
        let phase = self.relative_phase(coordinate)?;
        let group = self.relative_group(coordinate, phase)?;
        self.absolute_group(phase, group, air, coordinate_air)
    }

    /// Full ordinary wavelength curvature, including both independent air laws.
    /// It is evaluated only when requested, without differencing or resampling.
    pub fn refractive_index_jet(
        &self,
        wavelength_nm: f64,
    ) -> Result<CatalogRefractiveIndexJet, CatalogRefractionError> {
        let (coordinate, air, coordinate_air) = self.inputs(wavelength_nm)?;
        let relative = self.relative_phase(coordinate)?;
        let group = self.absolute_group(
            relative,
            self.relative_group(coordinate, relative)?,
            air,
            coordinate_air,
        )?;
        let phase = relative * air.phase_index;
        if !admitted_index(phase) {
            return Err(CatalogRefractionError::InvalidIndex);
        }
        let square = (coordinate * 1e-3) * (coordinate * 1e-3);
        let mut first = Dd::ZERO;
        let mut second = Dd::ZERO;
        for (b, c) in self.strengths.into_iter().zip(self.resonances_um2) {
            if b == 0.0 || c == 0.0 {
                continue;
            }
            let denominator = square - c;
            let ratio = c / denominator;
            let moment = b * (square / denominator) * ratio;
            first = first + Dd::from_f64(-2.0 * moment);
            second = second + Dd::from_f64(moment * (6.0 + 8.0 * ratio));
        }
        let first = first.to_f64() / (2.0 * relative);
        let second = second.to_f64() / (2.0 * relative) - first * (first / relative);
        let a1 = wavelength_nm * air.wavelength_derivative_per_nm;
        let a2 = wavelength_nm.powi(2) * air.wavelength_second_derivative_per_nm2;
        // x=lambda/B; u=lambda*x'/x, v=lambda²*x''/x. These are
        // ordinary derivatives, so v includes the -2*lambda*B'/B term.
        let (u, v) = coordinate_air.map_or((1.0, 0.0), |b| {
            let b1 = wavelength_nm * b.wavelength_derivative_per_nm / b.phase_index;
            let b2 = wavelength_nm.powi(2) * b.wavelength_second_derivative_per_nm2 / b.phase_index;
            (1.0 - b1, -2.0 * b1 - b2 + 2.0 * b1 * b1)
        });
        let d = Dd::from_f64;
        let scaled_wavelength_derivatives = [
            (d(a1) * d(relative) + d(air.phase_index) * d(first) * d(u)).to_f64(),
            (d(a2) * d(relative)
                + d(2.0) * d(a1) * d(first) * d(u)
                + d(air.phase_index) * (d(second) * d(u) * d(u) + d(first) * d(v)))
            .to_f64(),
        ];
        if scaled_wavelength_derivatives.iter().any(|v| !v.is_finite()) {
            return Err(CatalogRefractionError::InvalidIndex);
        }
        Ok(CatalogRefractiveIndexJet {
            phase_index: phase,
            group_index: group,
            scaled_wavelength_derivatives,
        })
    }

    pub fn content_identity(&self) -> ContentHash {
        let mut hash = DomainHasher::new("fs-material:air-relative-sellmeier:v1");
        self.reference_air.hash_into(&mut hash);
        match self.coordinate {
            CatalogWavelength::Vacuum => hash.update(&[0]),
            CatalogWavelength::Air(air) => {
                hash.update(&[1]);
                air.hash_into(&mut hash);
            }
        }
        for v in self
            .strengths
            .into_iter()
            .chain(self.resonances_um2)
            .chain([self.source_bounds_nm.0, self.source_bounds_nm.1])
        {
            hash.update(&v.to_bits().to_le_bytes());
        }
        hash.finalize()
    }

    fn relative_phase(&self, wavelength_nm: f64) -> Result<f64, CatalogRefractionError> {
        if !wavelength_nm.is_finite()
            || !(self.source_bounds_nm.0..=self.source_bounds_nm.1).contains(&wavelength_nm)
        {
            return Err(CatalogRefractionError::InvalidWavelength);
        }
        let wavelength_um = wavelength_nm * 1e-3;
        let square = wavelength_um * wavelength_um;
        let mut n2 = 1.0;
        for i in 0..3 {
            n2 += self.strengths[i] * (square / (square - self.resonances_um2[i]));
        }
        if !n2.is_finite() || !(1.0..=f64::MAX / 64.0).contains(&n2) {
            return Err(CatalogRefractionError::InvalidIndex);
        }
        Ok(det::sqrt(n2))
    }
    fn relative_group(
        &self,
        wavelength_nm: f64,
        phase: f64,
    ) -> Result<f64, CatalogRefractionError> {
        let square = (wavelength_nm * 1e-3) * (wavelength_nm * 1e-3);
        let mut correction = 0.0;
        for i in 0..3 {
            let denominator = square - self.resonances_um2[i];
            correction +=
                self.strengths[i] * (self.resonances_um2[i] / denominator) * (square / denominator);
        }
        let group = phase + correction / phase;
        if !group.is_finite() {
            return Err(CatalogRefractionError::InvalidIndex);
        }
        Ok(group)
    }
    fn absolute_group(
        &self,
        relative: f64,
        relative_group: f64,
        air: AirRefractiveSample,
        coordinate_air: Option<AirRefractiveSample>,
    ) -> Result<f64, CatalogRefractionError> {
        let group = match self.coordinate {
            CatalogWavelength::Vacuum => {
                relative_group * air.phase_index + relative * (air.group_index - air.phase_index)
            }
            CatalogWavelength::Air(law) if law == self.reference_air => {
                relative_group * air.group_index
            }
            CatalogWavelength::Air(_) => {
                let b = coordinate_air.expect("air wavelength sample");
                air.group_index * relative
                    + (air.phase_index / b.phase_index)
                        * b.group_index
                        * (relative_group - relative)
            }
        };
        if !group.is_finite() || group <= 0.0 {
            return Err(CatalogRefractionError::InvalidIndex);
        }
        Ok(group)
    }
    fn inputs(
        &self,
        vacuum_nm: f64,
    ) -> Result<(f64, AirRefractiveSample, Option<AirRefractiveSample>), CatalogRefractionError>
    {
        if !vacuum_nm.is_finite()
            || !(self.vacuum_bounds_nm.0..=self.vacuum_bounds_nm.1).contains(&vacuum_nm)
        {
            return Err(CatalogRefractionError::InvalidWavelength);
        }
        let (coordinate, coordinate_air) = self.source_wavelength(vacuum_nm)?;
        let air = match self.coordinate {
            CatalogWavelength::Air(law) if law == self.reference_air => {
                coordinate_air.expect("air wavelength sample")
            }
            _ => self
                .reference_air
                .sample(vacuum_nm)
                .map_err(CatalogRefractionError::Air)?,
        };
        Ok((coordinate, air, coordinate_air))
    }
    fn source_wavelength(
        &self,
        vacuum_nm: f64,
    ) -> Result<(f64, Option<AirRefractiveSample>), CatalogRefractionError> {
        match self.coordinate {
            CatalogWavelength::Vacuum => Ok((vacuum_nm, None)),
            CatalogWavelength::Air(air) => {
                let sample = air.sample(vacuum_nm).map_err(CatalogRefractionError::Air)?;
                Ok((vacuum_nm / sample.phase_index, Some(sample)))
            }
        }
    }
    fn inverse_air_wavelength(
        &self,
        source_nm: f64,
        lower: bool,
    ) -> Result<f64, CatalogRefractionError> {
        let [mut lo, mut hi] = WAVELENGTH_RANGE_NM;
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

fn admitted_index(index: f64) -> bool {
    index.is_finite() && index >= 1.0 && index * index <= f64::MAX / 64.0
}

#[cfg(test)]
mod tests;
