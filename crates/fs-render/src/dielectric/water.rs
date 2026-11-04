//! Exact native IAPWS phase/group law at one frozen ordinary-water state.
//!
//! The underlying material owner retains the original R9-97 coefficients and
//! domain. No Cauchy fit, sampled interpolation, equation of state, absorption
//! coefficient, or phase-mixture closure is substituted by this adapter.

use super::{DielectricError, DielectricIor};
use fs_blake3::ContentHash;
use fs_material::water_refraction::{IapwsWaterRefraction, WAVELENGTH_RANGE_NM};

/// Absolute ordinary-water phase index at the original temperature and density.
/// Supply pure-fluid state coordinates from the accepted physical snapshot.
/// Absorption remains an independent explicit law in `DielectricGlass`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IapwsWaterIor {
    temperature_k: f64,
    density_kg_m3: f64,
}

impl IapwsWaterIor {
    pub fn try_new(temperature_k: f64, density_kg_m3: f64) -> Result<Self, DielectricError> {
        IapwsWaterRefraction::prepare(589.0)
            .map_err(|_| DielectricError::InvalidWavelength)?
            .sample(temperature_k, density_kg_m3)
            .map_err(|_| DielectricError::InvalidIor)?;
        Ok(Self {
            temperature_k,
            density_kg_m3,
        })
    }

    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }

    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    pub fn source_identity(&self) -> ContentHash {
        IapwsWaterRefraction::source_identity()
    }

    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        (WAVELENGTH_RANGE_NM[0], WAVELENGTH_RANGE_NM[1])
    }

    pub fn eval(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        Ok(IapwsWaterRefraction::prepare(wavelength_nm)
            .map_err(|_| DielectricError::InvalidWavelength)?
            .sample(self.temperature_k, self.density_kg_m3)
            .map_err(|_| DielectricError::InvalidIor)?
            .phase_index)
    }

    pub fn group_index(&self, wavelength_nm: f64) -> Result<f64, DielectricError> {
        Ok(IapwsWaterRefraction::prepare(wavelength_nm)
            .map_err(|_| DielectricError::InvalidWavelength)?
            .sample(self.temperature_k, self.density_kg_m3)
            .map_err(|_| DielectricError::InvalidIor)?
            .group_index)
    }

    pub const fn is_dispersive(&self) -> bool {
        self.density_kg_m3 > 0.0
    }
}

impl From<IapwsWaterIor> for DielectricIor {
    fn from(ior: IapwsWaterIor) -> Self {
        Self::Water(ior)
    }
}

#[cfg(test)]
mod tests;
