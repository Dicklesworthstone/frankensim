//! Retarded material histories connected to the original thermal/QE owners.
//!
//! The callback is called once, at the admitted emission PROPER time. It
//! supplies the original thermal surface and the outward normal in that
//! source's rotation-free rest frame. The already-aberrated source direction
//! decides whether this side emits toward the receiver. Diffuse radiance gets
//! no additional cosine or clock multiplier here: projected area, visibility,
//! etendue and exposure remain the caller's separate geometry responsibilities.

use super::{
    ObservedThermalDetector, ObservedThermalSurface, ThermalDetectorSample, ThermalDopplerError,
};
use crate::spectral::{
    detector::QuantumEfficiencySensor,
    doppler::retarded::{
        PiecewiseVacuumWorldline, RetardedVacuumError, RetardedVacuumRay, UniformVacuumWorldline,
        VacuumReception,
    },
};
use crate::thermal_surface::{DiffuseThermalSurface, ThermalSurfaceError};
use fs_math::{dd::Dd, det};

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RetardedThermalError {
    Geometry(RetardedVacuumError),
    Source(ThermalSurfaceError),
    Observation(ThermalDopplerError),
    InvalidNormal,
    NonEmittingSide,
}

impl core::fmt::Display for RetardedThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "retarded thermal observation: {self:?}")
    }
}
impl core::error::Error for RetardedThermalError {}

/// One actual source state at an emission event, without copying its spectrum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalEmissionState<'source> {
    pub surface: DiffuseThermalSurface<'source>,
    /// Material rest-frame outward normal; normalized on admission.
    pub outward_normal: [f64; 3],
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RetardedThermalSurface<'source> {
    ray: RetardedVacuumRay,
    observed: ObservedThermalSurface<'source>,
    emission_cosine: f64,
}

impl<'source> RetardedThermalSurface<'source> {
    /// Resolve geometry before accessing the source. Invalid geometry and
    /// expired history never call `at_proper_time`; an unsupported source band
    /// continues to refuse through the original observation/sensor admission.
    pub fn observe(
        motion: &UniformVacuumWorldline,
        receiver: VacuumReception,
        max_residual_m: f64,
        at_proper_time: impl FnOnce(Dd) -> Result<ThermalEmissionState<'source>, ThermalSurfaceError>,
    ) -> Result<Self, RetardedThermalError> {
        let ray = motion
            .observe(receiver, max_residual_m)
            .map_err(RetardedThermalError::Geometry)?;
        Self::observe_ray(ray, at_proper_time)
    }

    /// Observe a finite nonuniform source history on its accumulated proper
    /// clock. The causal segment is selected before the single source query;
    /// velocity jumps use the trajectory's documented right-hand convention.
    pub fn observe_piecewise(
        motion: &PiecewiseVacuumWorldline,
        receiver: VacuumReception,
        max_residual_m: f64,
        at_proper_time: impl FnOnce(Dd) -> Result<ThermalEmissionState<'source>, ThermalSurfaceError>,
    ) -> Result<Self, RetardedThermalError> {
        let ray = motion
            .observe(receiver, max_residual_m)
            .map_err(RetardedThermalError::Geometry)?;
        Self::observe_ray(ray, at_proper_time)
    }

    fn observe_ray(
        ray: RetardedVacuumRay,
        at_proper_time: impl FnOnce(Dd) -> Result<ThermalEmissionState<'source>, ThermalSurfaceError>,
    ) -> Result<Self, RetardedThermalError> {
        let source =
            at_proper_time(ray.emission_proper_time_s()).map_err(RetardedThermalError::Source)?;
        let emission_cosine = cosine(source.outward_normal, ray.doppler().source_direction())?;
        let observed = ObservedThermalSurface::try_new(source.surface, ray.doppler())
            .map_err(RetardedThermalError::Observation)?;
        Ok(Self {
            ray,
            observed,
            emission_cosine,
        })
    }

    pub const fn ray(&self) -> &RetardedVacuumRay {
        &self.ray
    }

    /// The original spectral-radiance API, with this event's Lorentz view.
    pub const fn observed(&self) -> ObservedThermalSurface<'source> {
        self.observed
    }

    /// Source-frame hemisphere cosine; not an apparent projected-area map.
    pub const fn emission_cosine(&self) -> f64 {
        self.emission_cosine
    }

    pub fn with_sensor<'sensor, const C: usize>(
        self,
        sensor: &'sensor QuantumEfficiencySensor<C>,
    ) -> Result<RetardedThermalDetector<'source, 'sensor, C>, ThermalDopplerError> {
        Ok(RetardedThermalDetector {
            ray: self.ray,
            detector: self.observed.with_sensor(sensor)?,
            emission_cosine: self.emission_cosine,
        })
    }
}

/// Existing receiver-frame proposal and electron-rate conversion at one event.
pub struct RetardedThermalDetector<'source, 'sensor, const C: usize> {
    ray: RetardedVacuumRay,
    detector: ObservedThermalDetector<'source, 'sensor, C>,
    emission_cosine: f64,
}

impl<const C: usize> RetardedThermalDetector<'_, '_, C> {
    pub const fn ray(&self) -> &RetardedVacuumRay {
        &self.ray
    }

    pub const fn emission_cosine(&self) -> f64 {
        self.emission_cosine
    }

    pub fn sample(&self, u: f64) -> Result<ThermalDetectorSample<C>, ThermalDopplerError> {
        self.detector.sample(u)
    }
}

fn cosine(normal: [f64; 3], direction: [f64; 3]) -> Result<f64, RetardedThermalError> {
    if normal.iter().any(|x| !x.is_finite()) {
        return Err(RetardedThermalError::InvalidNormal);
    }
    let scale = normal.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
    if scale == 0.0 {
        return Err(RetardedThermalError::InvalidNormal);
    }
    let scaled = normal.map(|x| x / scale);
    let norm = det::sqrt(scaled.iter().map(|x| x * x).sum());
    let dot =
        (scaled[0] * direction[0] + scaled[1] * direction[1] + scaled[2] * direction[2]) / norm;
    if dot <= 0.0 {
        Err(RetardedThermalError::NonEmittingSide)
    } else {
        Ok(dot.min(1.0))
    }
}

#[cfg(test)]
mod tests;
