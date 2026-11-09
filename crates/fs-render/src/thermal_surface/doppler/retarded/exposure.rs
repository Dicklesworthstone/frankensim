//! Expected photoelectrons over finite detector-proper-time exposures.
//!
//! The supplied constant etendue is in the receiver rest frame and must already
//! include projected area, solid angle and throughput. No inverse-square,
//! incidence-cosine, aberration Jacobian or second Doppler factor is invented.
//! Temporal and spectral Gauss quadratures are estimates. QE knots and shifted
//! source-emissivity knots split every spectral panel, including narrow lines;
//! caller-authored time panels must resolve material-history changes. The
//! opt-in knot-aligned path resolves retarded source and detector turns too.
//! No visibility, shutter profile, charge noise, gain or saturation is inferred.

use super::{RetardedThermalError, RetardedThermalSurface, ThermalEmissionState};
use crate::spectral::{
    detector::{MAX_SENSOR_CHANNELS, QuantumEfficiencySensor},
    doppler::retarded::{
        PiecewiseVacuumWorldline, RetardedVacuumRay, UniformVacuumWorldline,
        exposure::{
            self, GAUSS4, RetardedExposureError, RetardedExposureIntegral, RetardedExposureSettings,
        },
    },
};
use crate::thermal_surface::{
    DiffuseThermalSurface, ThermalSurfaceError,
    doppler::{ObservedThermalDetector, ThermalDopplerError},
};
use fs_math::dd::Dd;

const CHANNELS: usize = 2 * MAX_SENSOR_CHANNELS;
const PHOTONS_PER_JOULE_PER_NM: f64 = 1.0e-9 / (6.626_070_15e-34 * 299_792_458.0);

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalExposureSettings {
    pub time: RetardedExposureSettings,
    pub receiver_etendue_m2_sr: f64,
    /// Positive uniform subdivisions of each merged smooth wavelength interval.
    pub wavelength_subdivisions: usize,
    /// Total Planck/QE point evaluations across the complete exposure.
    pub max_spectral_evaluations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ThermalExposureError {
    Etendue,
    SpectralQuadrature,
    SpectralBudget,
    NumericalRange,
    Cancelled,
    Thermal(RetardedThermalError),
    Detector(ThermalDopplerError),
}

/// Expected collected charge and QE-weighted optical energy, before gain/noise.
/// QE-weighted joules are a response integral, not absorbed heat in the sensor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RetardedThermalExposure<const C: usize> {
    pub expected_electrons: [f64; C],
    pub qe_weighted_energy_j: [f64; C],
    pub receiver_proper_time_s: [Dd; 2],
    pub ray_queries: usize,
    pub spectral_evaluations: usize,
}

impl<const C: usize> RetardedThermalExposure<C> {
    /// Integrate uniform endpoint histories. Original source state is read once
    /// per temporal node, after complete endpoint-history admission. No spectrum
    /// is resampled or allocated; unsupported bands or cancellation return no
    /// partial exposure. The callback must not advance the source's state.
    pub fn uniform<'source>(
        source: &UniformVacuumWorldline,
        receiver: &UniformVacuumWorldline,
        sensor: &QuantumEfficiencySensor<C>,
        panel_edges_s: &[Dd],
        settings: ThermalExposureSettings,
        mut at_proper_time: impl FnMut(Dd) -> Result<ThermalEmissionState<'source>, ThermalSurfaceError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RetardedExposureError<ThermalExposureError>> {
        admit(settings)?;
        let mut work = 0;
        let integral = exposure::integrate_uniform(
            source,
            receiver,
            panel_edges_s,
            settings.time,
            |ray, checkpoint| {
                rate(
                    ray,
                    sensor,
                    settings,
                    &mut work,
                    &mut at_proper_time,
                    checkpoint,
                )
            },
            cancelled,
        )?;
        Ok(Self::from_integral(integral, work))
    }

    /// The same exposure over finite source and detector trajectories. Panel
    /// boundaries are expressed on the detector's accumulated proper clock.
    pub fn piecewise<'source>(
        source: &PiecewiseVacuumWorldline,
        receiver: &PiecewiseVacuumWorldline,
        sensor: &QuantumEfficiencySensor<C>,
        panel_edges_s: &[Dd],
        settings: ThermalExposureSettings,
        mut at_proper_time: impl FnMut(Dd) -> Result<ThermalEmissionState<'source>, ThermalSurfaceError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RetardedExposureError<ThermalExposureError>> {
        admit(settings)?;
        let mut work = 0;
        let integral = exposure::integrate_piecewise(
            source,
            receiver,
            panel_edges_s,
            settings.time,
            |ray, checkpoint| {
                rate(
                    ray,
                    sensor,
                    settings,
                    &mut work,
                    &mut at_proper_time,
                    checkpoint,
                )
            },
            cancelled,
        )?;
        Ok(Self::from_integral(integral, work))
    }

    /// Also split time panels at detector turns and received source turns.
    /// Arrival planning and complete motion-history admission finish before
    /// reading the thermal source. Other material discontinuities still need
    /// authored panel edges; both quadratures remain numerical estimates.
    pub fn piecewise_with_knots<'source>(
        source: &PiecewiseVacuumWorldline,
        receiver: &PiecewiseVacuumWorldline,
        sensor: &QuantumEfficiencySensor<C>,
        panel_edges_s: &[Dd],
        settings: ThermalExposureSettings,
        mut at_proper_time: impl FnMut(Dd) -> Result<ThermalEmissionState<'source>, ThermalSurfaceError>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RetardedExposureError<ThermalExposureError>> {
        admit(settings)?;
        let mut work = 0;
        let integral = exposure::integrate_piecewise_with_knots(
            source,
            receiver,
            panel_edges_s,
            settings.time,
            |ray, checkpoint| {
                rate(
                    ray,
                    sensor,
                    settings,
                    &mut work,
                    &mut at_proper_time,
                    checkpoint,
                )
            },
            cancelled,
        )?;
        Ok(Self::from_integral(integral, work))
    }

    fn from_integral(integral: RetardedExposureIntegral<CHANNELS>, work: usize) -> Self {
        Self {
            expected_electrons: core::array::from_fn(|i| {
                integral.integrated[MAX_SENSOR_CHANNELS + i]
            }),
            qe_weighted_energy_j: core::array::from_fn(|i| integral.integrated[i]),
            receiver_proper_time_s: integral.receiver_proper_time_s,
            ray_queries: integral.ray_queries,
            spectral_evaluations: work,
        }
    }
}

fn admit(
    settings: ThermalExposureSettings,
) -> Result<(), RetardedExposureError<ThermalExposureError>> {
    if !settings.receiver_etendue_m2_sr.is_finite() || settings.receiver_etendue_m2_sr < 0.0 {
        return Err(RetardedExposureError::Observable(
            ThermalExposureError::Etendue,
        ));
    }
    if settings.wavelength_subdivisions == 0 {
        return Err(RetardedExposureError::Observable(
            ThermalExposureError::SpectralQuadrature,
        ));
    }
    if settings
        .wavelength_subdivisions
        .checked_mul(GAUSS4.len())
        .is_none_or(|n| n > settings.max_spectral_evaluations)
    {
        return Err(RetardedExposureError::Observable(
            ThermalExposureError::SpectralBudget,
        ));
    }
    Ok(())
}

fn rate<'source, const C: usize>(
    ray: &RetardedVacuumRay,
    sensor: &QuantumEfficiencySensor<C>,
    settings: ThermalExposureSettings,
    work: &mut usize,
    at_proper_time: &mut impl FnMut(Dd) -> Result<ThermalEmissionState<'source>, ThermalSurfaceError>,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<[f64; CHANNELS], ThermalExposureError> {
    let surface = RetardedThermalSurface::observe_ray(*ray, at_proper_time)
        .map_err(ThermalExposureError::Thermal)?;
    let detector = surface
        .observed()
        .with_sensor(sensor)
        .map_err(ThermalExposureError::Detector)?;
    spectrum(detector, settings, work, cancelled)
}

fn spectrum<const C: usize>(
    detector: ObservedThermalDetector<'_, '_, C>,
    settings: ThermalExposureSettings,
    work: &mut usize,
    cancelled: &mut dyn FnMut() -> bool,
) -> Result<[f64; CHANNELS], ThermalExposureError> {
    let source_knots = match detector.source.source() {
        DiffuseThermalSurface::Gray(_) => &[][..],
        DiffuseThermalSurface::Spectral(source) => source.emissivity_knots(),
        DiffuseThermalSurface::Rough(source) => source.facet_albedo().knots(),
    };
    let shifted = |index: usize| {
        detector
            .source
            .doppler()
            .receiver_wavelength_nm(source_knots[index].0)
            .map_err(|e| ThermalExposureError::Detector(e.into()))
    };
    let mut previous = 0.0;
    for index in 0..source_knots.len() {
        let wavelength = shifted(index)?;
        if wavelength <= previous {
            return Err(ThermalExposureError::NumericalRange);
        }
        previous = wavelength;
    }
    let sensor_knots = detector.sensor.knots();
    let (mut lower, band_upper) = detector.sensor.wavelength_bounds_nm();
    let (mut sensor_index, mut source_index) = (1, 0);
    let mut totals = [Dd::ZERO; CHANNELS];
    while lower < band_upper {
        while sensor_index < sensor_knots.len() && sensor_knots[sensor_index].0 <= lower {
            sensor_index += 1;
        }
        while source_index < source_knots.len() && shifted(source_index)? <= lower {
            source_index += 1;
        }
        let mut upper = sensor_knots.get(sensor_index).map_or(band_upper, |k| k.0);
        if source_index < source_knots.len() {
            upper = upper.min(shifted(source_index)?);
        }
        let width = upper - lower;
        if !width.is_finite() || width <= 0.0 {
            return Err(ThermalExposureError::NumericalRange);
        }
        let subdivisions = settings.wavelength_subdivisions as f64;
        for panel in 0..settings.wavelength_subdivisions {
            let panel_lower = lower + width * (panel as f64 / subdivisions);
            let panel_upper = if panel + 1 == settings.wavelength_subdivisions {
                upper
            } else {
                lower + width * ((panel + 1) as f64 / subdivisions)
            };
            let panel_width = panel_upper - panel_lower;
            previous = panel_lower;
            for (fraction, weight) in GAUSS4 {
                if cancelled() {
                    return Err(ThermalExposureError::Cancelled);
                }
                if *work >= settings.max_spectral_evaluations {
                    return Err(ThermalExposureError::SpectralBudget);
                }
                *work += 1;
                let wavelength = panel_lower + panel_width * fraction;
                if wavelength <= previous || wavelength >= panel_upper {
                    return Err(ThermalExposureError::NumericalRange);
                }
                previous = wavelength;
                let radiance = detector
                    .source
                    .radiance_w_m2_sr_nm(wavelength)
                    .map_err(ThermalExposureError::Detector)?;
                let qe = detector
                    .sensor
                    .quantum_efficiency(wavelength)
                    .map_err(|e| ThermalExposureError::Detector(e.into()))?;
                let measure = Dd::from_f64(radiance)
                    * ((Dd::from_f64(panel_width) * Dd::from_f64(weight))
                        * Dd::from_f64(settings.receiver_etendue_m2_sr));
                for (channel, efficiency) in qe.into_iter().enumerate() {
                    let energy = measure * Dd::from_f64(efficiency);
                    if radiance > 0.0
                        && efficiency > 0.0
                        && settings.receiver_etendue_m2_sr > 0.0
                        && energy == Dd::ZERO
                    {
                        return Err(ThermalExposureError::NumericalRange);
                    }
                    totals[channel] = totals[channel] + energy;
                    totals[MAX_SENSOR_CHANNELS + channel] = totals[MAX_SENSOR_CHANNELS + channel]
                        + energy
                            * (Dd::from_f64(wavelength) * Dd::from_f64(PHOTONS_PER_JOULE_PER_NM));
                }
            }
        }
        lower = upper;
    }
    let rates = totals.map(Dd::to_f64);
    if rates
        .iter()
        .zip(totals)
        .any(|(&x, dd)| !x.is_finite() || x < 0.0 || (x == 0.0 && dd != Dd::ZERO))
    {
        return Err(ThermalExposureError::NumericalRange);
    }
    Ok(rates)
}

#[cfg(test)]
mod tests;
