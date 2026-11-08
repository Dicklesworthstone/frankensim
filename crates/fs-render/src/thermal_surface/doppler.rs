//! Moving observations of the original thermal source and physical detector.
//!
//! This view leaves thermodynamic temperature, emissivity and source identity
//! owned by `DiffuseThermalSurface`. Its ray's `VacuumDoppler` determines the
//! rest-frame wavelength queried in that original law. The observer's sensor
//! retains its own wavelength response, photon energy and importance proposal.
//! No shifted tables, RGB conversion, resampling grid or per-query allocation
//! are introduced. A finite source band must cover the complete admitted
//! detector band after shifting; unsupported source wavelengths never become
//! invented black samples or extrapolated data.
//!
//! This is directional emitted radiance, before visibility, apparent projected
//! area, detector etendue and exposure. It does not solve retarded geometry,
//! transform surface normals, move existing scene instances, debit thermal
//! storage, or treat a dielectric/moving-medium path as vacuum.
//! [`retarded`] selects the original thermal history and emitting hemisphere
//! at the source event supplied by a finite-worldline light-cone solution, and
//! its exposure owner integrates these observations on the detector's clock.

use super::{DiffuseThermalSurface, ThermalSurfaceError, planck_radiance_w_m2_sr_nm};
use crate::spectral::{
    detector::{QuantumEfficiencySensor, SpectralSensorError},
    doppler::{DopplerError, VacuumDoppler},
};

/// Original source histories selected on the retarded material proper clock.
pub mod retarded;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ThermalDopplerError {
    Doppler(DopplerError),
    Source(ThermalSurfaceError),
    Sensor(SpectralSensorError),
    /// The receiver wavelength or complete sensor band is not covered by the
    /// original source's shifted finite wavelength support.
    SourceBand,
    NumericalRange,
}

impl core::fmt::Display for ThermalDopplerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "thermal Doppler observation: {self:?}")
    }
}
impl core::error::Error for ThermalDopplerError {}
impl From<DopplerError> for ThermalDopplerError {
    fn from(value: DopplerError) -> Self {
        Self::Doppler(value)
    }
}
impl From<ThermalSurfaceError> for ThermalDopplerError {
    fn from(value: ThermalSurfaceError) -> Self {
        Self::Source(value)
    }
}
impl From<SpectralSensorError> for ThermalDopplerError {
    fn from(value: SpectralSensorError) -> Self {
        Self::Sensor(value)
    }
}

/// Source-preserving directional emission into vacuum between inertial frames.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ObservedThermalSurface<'source> {
    source: DiffuseThermalSurface<'source>,
    doppler: VacuumDoppler,
    receiver_band_nm: Option<(f64, f64)>,
    apparent_temperature_k: f64,
}

impl<'source> ObservedThermalSurface<'source> {
    /// The original gray and spectral sources are direction independent in
    /// their rest frame. A rough thermal surface additionally needs a physical
    /// rest-frame surface normal, which this Doppler-only view does not own.
    pub fn try_new(
        source: DiffuseThermalSurface<'source>,
        doppler: VacuumDoppler,
    ) -> Result<Self, ThermalDopplerError> {
        let receiver_band_nm = match source {
            DiffuseThermalSurface::Gray(_) => None,
            DiffuseThermalSurface::Rough(_) => {
                return Err(ThermalSurfaceError::DirectionRequired.into());
            }
            DiffuseThermalSurface::Spectral(surface) => {
                let (low, high) = surface.emissivity().wavelength_bounds_nm();
                let band = (
                    doppler.receiver_wavelength_nm(low)?,
                    doppler.receiver_wavelength_nm(high)?,
                );
                if band.0 >= band.1 {
                    return Err(ThermalDopplerError::NumericalRange);
                }
                Some(band)
            }
        };
        let apparent_temperature_k = source.temperature_k() * doppler.frequency_ratio();
        if !apparent_temperature_k.is_finite() || apparent_temperature_k <= 0.0 {
            return Err(ThermalDopplerError::NumericalRange);
        }
        Ok(Self {
            source,
            doppler,
            receiver_band_nm,
            apparent_temperature_k,
        })
    }

    /// Original thermal source, including its unchanged actual temperature.
    pub const fn source(&self) -> DiffuseThermalSurface<'source> {
        self.source
    }

    pub const fn doppler(&self) -> VacuumDoppler {
        self.doppler
    }

    /// Closed support in the receiver's vacuum nanometres. `None` denotes the
    /// existing gray law's positive-wavelength domain, not a measured spectrum.
    pub const fn receiver_wavelength_bounds_nm(&self) -> Option<(f64, f64)> {
        self.receiver_band_nm
    }

    fn source_wavelength_nm(&self, receiver_nm: f64) -> Result<f64, ThermalDopplerError> {
        let source_nm = self.doppler.source_wavelength_nm(receiver_nm)?;
        if let Some((low, high)) = self.receiver_band_nm {
            if !(low..=high).contains(&receiver_nm) {
                return Err(ThermalDopplerError::SourceBand);
            }
            if let DiffuseThermalSurface::Spectral(surface) = self.source {
                let (source_low, source_high) = surface.emissivity().wavelength_bounds_nm();
                // The query was admitted in the transformed closed source
                // band. Retain its endpoints through divide/multiply roundoff;
                // no wavelength outside that band reaches this clamp.
                return Ok(source_nm.clamp(source_low, source_high));
            }
        }
        Ok(source_nm)
    }

    /// Receiver-frame spectral radiance in W/(m² sr nm).
    ///
    /// Algebraically `D^5 epsilon(D lambda) B_lambda(D lambda,T)` equals
    /// `epsilon(D lambda) B_lambda(lambda,D T)`. The latter uses the original
    /// stable Planck owner once and avoids underflowing a cold rest-frame Wien
    /// tail before its Lorentz amplification is applied. `D T` is only a
    /// directional color temperature; it never mutates the thermal state.
    pub fn radiance_w_m2_sr_nm(&self, receiver_nm: f64) -> Result<f64, ThermalDopplerError> {
        let source_nm = self.source_wavelength_nm(receiver_nm)?;
        self.radiance_at_source_wavelength(receiver_nm, source_nm)
    }

    fn radiance_at_source_wavelength(
        &self,
        receiver_nm: f64,
        source_nm: f64,
    ) -> Result<f64, ThermalDopplerError> {
        if self.doppler.frequency_ratio() == 1.0 {
            return Ok(self.source.radiance_w_m2_sr_nm(source_nm)?);
        }
        let emissivity = match self.source {
            DiffuseThermalSurface::Gray(surface) => surface.emissivity(),
            DiffuseThermalSurface::Spectral(surface) => surface.emissivity().at(source_nm)?,
            DiffuseThermalSurface::Rough(_) => {
                return Err(ThermalSurfaceError::DirectionRequired.into());
            }
        };
        if emissivity == 0.0 {
            return Ok(0.0);
        }
        let radiance =
            emissivity * planck_radiance_w_m2_sr_nm(self.apparent_temperature_k, receiver_nm)?;
        if radiance.is_finite() && radiance >= 0.0 {
            Ok(radiance)
        } else {
            Err(ThermalDopplerError::NumericalRange)
        }
    }

    /// Prepare the existing physical sensor in its OWN rest frame. The entire
    /// sensor band is checked before a random draw or spectral evaluation.
    pub fn with_sensor<'sensor, const C: usize>(
        self,
        sensor: &'sensor QuantumEfficiencySensor<C>,
    ) -> Result<ObservedThermalDetector<'source, 'sensor, C>, ThermalDopplerError> {
        let (low, high) = sensor.wavelength_bounds_nm();
        self.source_wavelength_nm(low)?;
        self.source_wavelength_nm(high)?;
        Ok(ObservedThermalDetector {
            source: self,
            sensor,
        })
    }
}

/// One immutable prepared pairing; no new source or detector data is copied.
pub struct ObservedThermalDetector<'source, 'sensor, const CHANNELS: usize> {
    source: ObservedThermalSurface<'source>,
    sensor: &'sensor QuantumEfficiencySensor<CHANNELS>,
}

/// One detector-proposal draw. Average independent or stratified draws to
/// estimate channel integrals; these values precede etendue and exposure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalDetectorSample<const CHANNELS: usize> {
    pub receiver_wavelength_nm: f64,
    pub source_wavelength_nm: f64,
    pub spectral_radiance_w_m2_sr_nm: f64,
    pub qe_weighted_radiance_w_m2_sr: [f64; CHANNELS],
    pub electron_rate_per_m2_sr: [f64; CHANNELS],
}

impl<const C: usize> ObservedThermalDetector<'_, '_, C> {
    /// O(log K + C) work, one existing Planck/source query, no allocation.
    /// Electron energy is h*c/lambda_RECEIVER, never the emitting wavelength.
    pub fn sample(&self, u: f64) -> Result<ThermalDetectorSample<C>, ThermalDopplerError> {
        let draw = self.sensor.sample_wavelength(u)?;
        let receiver_nm = draw.wavelength_nm();
        let source_nm = self.source.source_wavelength_nm(receiver_nm)?;
        let radiance = self
            .source
            .radiance_at_source_wavelength(receiver_nm, source_nm)?;
        let mut power = [0.0; C];
        let mut electrons = [0.0; C];
        for i in 0..C {
            power[i] = weighted(radiance, draw.radiance_weights_nm()[i])?;
            electrons[i] = weighted(radiance, draw.electron_weights_nm_per_joule()[i])?;
        }
        Ok(ThermalDetectorSample {
            receiver_wavelength_nm: receiver_nm,
            source_wavelength_nm: source_nm,
            spectral_radiance_w_m2_sr_nm: radiance,
            qe_weighted_radiance_w_m2_sr: power,
            electron_rate_per_m2_sr: electrons,
        })
    }
}

fn weighted(radiance: f64, weight: f64) -> Result<f64, ThermalDopplerError> {
    let value = radiance * weight;
    if value.is_finite() && (value > 0.0 || radiance == 0.0 || weight == 0.0) {
        Ok(value)
    } else {
        Err(ThermalDopplerError::NumericalRange)
    }
}

#[cfg(test)]
mod tests;
