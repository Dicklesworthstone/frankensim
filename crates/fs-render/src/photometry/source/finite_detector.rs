//! Calibrated received flux over a finite native detector, including caustics.
//!
//! Source directions follow the represented angular intensity. Every packet
//! carries physical spectral power divided by its direction PDF; the native
//! curved scene owns propagation, absorption, polarization and stochastic
//! interface branches. A detector hit selects the receiving area. No second
//! cosine, inverse-square factor, or inverse ray-bundle Jacobian is applied.
//! Thus a focus does not create a point-irradiance division by zero here.
//!
//! Each independent Philox sample traces a complete 80-bin spectrum, sharing
//! its source direction and interface random numbers across wavelengths.
//! Broadband sample errors are measured AFTER spectral integration, retaining
//! those correlations. Reported standard errors are empirical Monte Carlo
//! estimates, not confidence bounds or geometry/spectral/model certificates.
//! In particular, no observed hits gives zero empirical error but does not
//! prove zero receiving power. Small acceptance angles may require many rays.
//!
//! All launches, including escapes and blocked paths, enter the denominator.
//! Surviving paths stopped by a numerical/physical work limit refuse the whole
//! measurement. Only the native homogeneous-medium escape proof makes a miss
//! complete; a GRIN field may bend an apparently outward ray back to a surface.

use super::sampling::SourceSamplingError;
use super::{PhotometricSource, SPECTRAL_BINS, wavelength_nm};
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedSceneError, CurvedTraceSettings,
};
use crate::polarization::dielectric::scene::measurement::MeterError;
use crate::polarization::dielectric::transport::{BeamError, PolarizationFrame};
use crate::polarization::{Stokes, photometry::PolarizedSource};
use fs_exec::Cx;
use fs_geom::Vec3;
use fs_math::det;
use fs_rand::philox::philox4x32_10;

/// Explicit ceiling on complete spectral samples in one measurement.
pub const MAX_FINITE_DETECTOR_SAMPLES: u32 = 1_048_576;

/// One actual scene detector and a bounded independent-sample measurement.
#[derive(Clone, Copy, Debug)]
pub struct FiniteDetectorRequest {
    pub detector_object_id: u64,
    /// World positive-Q axis projected into each launch direction. An exactly
    /// parallel sample refuses because that source polarization has no frame.
    pub source_transverse_axis: Vec3,
    pub initial_medium: usize,
    pub trace: CurvedTraceSettings,
    /// At least two complete spectra, so a sample variance is defined.
    pub samples: u32,
    pub seed: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FiniteDetectorError {
    InvalidRequest,
    Sampling(SourceSamplingError),
    SourceFrame(BeamError),
    Transport {
        sample: u32,
        bin: usize,
        error: CurvedSceneError,
    },
    /// This surviving path is unresolved, not an absorbed or escaped sample.
    IncompletePath {
        sample: u32,
        bin: usize,
        termination: PathTermination,
    },
    Meter(MeterError),
    Cancelled,
    NumericalRange,
}
impl core::fmt::Display for FiniteDetectorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite spectral detector: {self:?}")
    }
}
impl core::error::Error for FiniteDetectorError {}
impl From<fs_exec::Cancelled> for FiniteDetectorError {
    fn from(_: fs_exec::Cancelled) -> Self {
        Self::Cancelled
    }
}
impl From<MeterError> for FiniteDetectorError {
    fn from(error: MeterError) -> Self {
        Self::Meter(error)
    }
}
impl From<SourceSamplingError> for FiniteDetectorError {
    fn from(error: SourceSamplingError) -> Self {
        Self::Sampling(error)
    }
}

/// Complete received flux and empirical standard errors in detector frames.
/// Fixed-size storage; path count does not increase retained memory.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteDetectorSpectrum {
    spectral: [Stokes; SPECTRAL_BINS],
    spectral_error: [[f64; 4]; SPECTRAL_BINS],
    radiant: Stokes,
    photopic: Stokes,
    radiant_error: [f64; 4],
    photopic_error: [f64; 4],
    samples: u32,
    detected_paths: u64,
    received_samples: u32,
    max_contact_resolution_m: f64,
}
impl FiniteDetectorSpectrum {
    pub const fn spectral_power_w_nm(&self) -> &[Stokes; SPECTRAL_BINS] {
        &self.spectral
    }
    /// Componentwise standard errors ordered I, Q, U, V. These are not Stokes
    /// states and cannot be propagated as though their components covary freely.
    pub const fn spectral_standard_error_w_nm(&self) -> &[[f64; 4]; SPECTRAL_BINS] {
        &self.spectral_error
    }
    pub const fn radiant_stokes_w(&self) -> Stokes {
        self.radiant
    }
    pub const fn photopic_stokes_lm(&self) -> Stokes {
        self.photopic
    }
    pub const fn radiant_standard_error_w(&self) -> [f64; 4] {
        self.radiant_error
    }
    pub const fn photopic_standard_error_lm(&self) -> [f64; 4] {
        self.photopic_error
    }
    pub const fn sample_count(&self) -> u32 {
        self.samples
    }
    pub const fn launched_paths(&self) -> u64 {
        self.samples as u64 * SPECTRAL_BINS as u64
    }
    pub const fn detected_paths(&self) -> u64 {
        self.detected_paths
    }
    /// Independent spectral samples with a hit on this detector in any bin.
    pub const fn received_samples(&self) -> u32 {
        self.received_samples
    }
    /// Numerical native contact collar, not an accumulated optical error bound.
    pub const fn max_contact_resolution_m(&self) -> f64 {
        self.max_contact_resolution_m
    }

    /// Area average for an explicitly supplied actual receiving area in m².
    /// The caller owns that area and its match to the detector geometry. This
    /// is an aperture average, not the local irradiance at a selected point.
    pub fn area_average(
        &self,
        area_m2: f64,
    ) -> Result<AreaAveragedDetectorReading, FiniteDetectorError> {
        if !area_m2.is_finite() || area_m2 <= 0.0 {
            return Err(FiniteDetectorError::InvalidRequest);
        }
        let divide = |values: [f64; 4]| -> Result<[f64; 4], FiniteDetectorError> {
            let out = values.map(|value| value / area_m2);
            if out.iter().all(|value| value.is_finite()) {
                Ok(out)
            } else {
                Err(FiniteDetectorError::NumericalRange)
            }
        };
        Ok(AreaAveragedDetectorReading {
            radiant_w_m2: physical(divide(components(self.radiant))?)?,
            photopic_lux: physical(divide(components(self.photopic))?)?,
            radiant_standard_error_w_m2: divide(self.radiant_error)?,
            photopic_standard_error_lux: divide(self.photopic_error)?,
        })
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AreaAveragedDetectorReading {
    pub radiant_w_m2: Stokes,
    pub photopic_lux: Stokes,
    pub radiant_standard_error_w_m2: [f64; 4],
    pub photopic_standard_error_lux: [f64; 4],
}

impl PhotometricSource {
    /// Measure finite-detector power with an initially unpolarized fixture.
    /// Native interfaces retain induced polarization and physical branch power.
    pub fn measure_detector(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: FiniteDetectorRequest,
    ) -> Result<FiniteDetectorSpectrum, FiniteDetectorError> {
        measure(self, cx, scene, request, |_| Stokes::unpolarized(1.0))
    }
}
impl PolarizedSource {
    /// Measure the surviving source spectrum and polarization. Upstream
    /// filtering is already in the intensity source and is not applied again.
    pub fn measure_detector(
        &self,
        cx: &Cx<'_>,
        scene: &CurvedPolarizedScene<'_>,
        request: FiniteDetectorRequest,
    ) -> Result<FiniteDetectorSpectrum, FiniteDetectorError> {
        measure(self.intensity_source(), cx, scene, request, |bin| {
            let s = self.relative_stokes()[bin];
            if s.i == 0.0 {
                Stokes::unpolarized(1.0)
            } else {
                Stokes {
                    i: 1.0,
                    q: s.q / s.i,
                    u: s.u / s.i,
                    v: s.v / s.i,
                }
            }
        })
    }
}

fn measure(
    source: &PhotometricSource,
    cx: &Cx<'_>,
    scene: &CurvedPolarizedScene<'_>,
    request: FiniteDetectorRequest,
    polarization: impl Fn(usize) -> Stokes,
) -> Result<FiniteDetectorSpectrum, FiniteDetectorError> {
    cx.checkpoint()?;
    if !(2..=MAX_FINITE_DETECTOR_SAMPLES).contains(&request.samples) {
        return Err(FiniteDetectorError::InvalidRequest);
    }
    let source_axis = super::unit(request.source_transverse_axis)
        .map_err(|_| FiniteDetectorError::InvalidRequest)?;
    let empty_meter = scene.power_meter(request.detector_object_id)?;
    let sampler = source.prepare_direction_sampler(cx)?;
    let powers = source
        .spectrum()
        .bins_w_per_lumen_nm()
        .map(|value| value * source.flux_lm());
    let states: [Stokes; SPECTRAL_BINS] = core::array::from_fn(polarization);
    if powers.iter().any(|value| !value.is_finite()) {
        return Err(FiniteDetectorError::NumericalRange);
    }
    let mut spectral = [Moments::default(); SPECTRAL_BINS];
    let mut radiant = Moments::default();
    let mut photopic = Moments::default();
    let mut detected_paths = 0;
    let mut received_samples = 0;
    let mut max_contact_resolution_m = 0.0_f64;
    for sample in 0..request.samples {
        cx.checkpoint()?;
        // The source domain is disjoint from the native interface domain
        // 0x706f_6c32, even with the same seed and sample/path identifier.
        let random = philox4x32_10(
            [sample, 0, 0, 0x736f_7572],
            [request.seed as u32, (request.seed >> 32) as u32],
        );
        let direction = sampler.sample(core::array::from_fn(|i| {
            (f64::from(random[i]) + 0.5) / 4_294_967_296.0
        }))?;
        let frame = PolarizationFrame::try_new(direction.direction, source_axis)
            .map_err(FiniteDetectorError::SourceFrame)?;
        let mut meter = empty_meter.clone();
        for bin in 0..SPECTRAL_BINS {
            cx.checkpoint()?;
            let intensity = powers[bin] * direction.angular_weight();
            let state = states[bin];
            let power = physical([
                intensity,
                intensity * state.q,
                intensity * state.u,
                intensity * state.v,
            ])?;
            let transport_error = |error| FiniteDetectorError::Transport { sample, bin, error };
            let beam = scene
                .launch_beam(
                    source.position_m(),
                    frame,
                    wavelength_nm(bin),
                    request.initial_medium,
                    power,
                )
                .map_err(transport_error)?;
            let path = scene
                .trace(
                    cx,
                    beam,
                    request.initial_medium,
                    request.trace,
                    request.seed,
                    u64::from(sample),
                )
                .map_err(transport_error)?;
            if path.beam().power().i > 0.0
                && matches!(
                    path.termination(),
                    PathTermination::DistanceLimit
                        | PathTermination::SurfaceLimit
                        | PathTermination::VolumeLimit
                )
            {
                return Err(FiniteDetectorError::IncompletePath {
                    sample,
                    bin,
                    termination: path.termination(),
                });
            }
            meter.record_curved(bin, &path)?;
        }
        let reading = meter.reading()?.require_complete()?;
        let count = sample + 1;
        for (bin, moments) in spectral.iter_mut().enumerate() {
            moments.add(meter.spectral_mean_w_nm(bin)?, count)?;
        }
        radiant.add(reading.radiant_stokes_w(), count)?;
        photopic.add(reading.photopic_stokes_lm(), count)?;
        detected_paths += reading.detected_paths();
        received_samples += u32::from(reading.detected_paths() > 0);
        max_contact_resolution_m = max_contact_resolution_m.max(reading.max_contact_resolution_m());
    }
    let mut bins = [Stokes::unpolarized(0.0); SPECTRAL_BINS];
    let mut errors = [[0.0; 4]; SPECTRAL_BINS];
    for (bin, moments) in spectral.into_iter().enumerate() {
        bins[bin] = moments.mean()?;
        errors[bin] = moments.standard_error(request.samples)?;
    }
    cx.checkpoint()?;
    Ok(FiniteDetectorSpectrum {
        spectral: bins,
        spectral_error: errors,
        radiant: radiant.mean()?,
        photopic: photopic.mean()?,
        radiant_error: radiant.standard_error(request.samples)?,
        photopic_error: photopic.standard_error(request.samples)?,
        samples: request.samples,
        detected_paths,
        received_samples,
        max_contact_resolution_m,
    })
}

/// Scaled Welford moments avoid squaring dimensionful powers (which can
/// overflow even when both the mean and its standard error are representable).
#[derive(Clone, Copy, Debug, Default)]
struct Moments {
    scale: [f64; 4],
    mean: [f64; 4],
    m2: [f64; 4],
}
impl Moments {
    fn add(&mut self, value: Stokes, count: u32) -> Result<(), FiniteDetectorError> {
        for (axis, value) in components(value).into_iter().enumerate() {
            let magnitude = value.abs();
            if magnitude > self.scale[axis] {
                let ratio = self.scale[axis] / magnitude;
                self.mean[axis] *= ratio;
                self.m2[axis] *= ratio * ratio;
                self.scale[axis] = magnitude;
            }
            if self.scale[axis] == 0.0 {
                continue;
            }
            let scaled = value / self.scale[axis];
            let delta = scaled - self.mean[axis];
            self.mean[axis] += delta / f64::from(count);
            self.m2[axis] += delta * (scaled - self.mean[axis]);
            if !self.m2[axis].is_finite() || self.m2[axis] < 0.0 {
                return Err(FiniteDetectorError::NumericalRange);
            }
        }
        Ok(())
    }
    fn mean(self) -> Result<Stokes, FiniteDetectorError> {
        physical(core::array::from_fn(|axis| {
            self.scale[axis] * self.mean[axis]
        }))
    }
    fn standard_error(self, count: u32) -> Result<[f64; 4], FiniteDetectorError> {
        let denominator = f64::from(count) * f64::from(count - 1);
        let errors =
            core::array::from_fn(|axis| self.scale[axis] * det::sqrt(self.m2[axis] / denominator));
        if errors.iter().all(|value| value.is_finite()) {
            Ok(errors)
        } else {
            Err(FiniteDetectorError::NumericalRange)
        }
    }
}
fn components(value: Stokes) -> [f64; 4] {
    [value.i, value.q, value.u, value.v]
}
fn physical(value: [f64; 4]) -> Result<Stokes, FiniteDetectorError> {
    Stokes::try_new(value[0], value[1], value[2], value[3])
        .ok_or(FiniteDetectorError::NumericalRange)
}

#[cfg(test)]
mod tests;
