//! Finite-time photodetection of an explicitly coherent optical line spectrum.
//!
//! The detector observes the squared analytic optical field, normalized in
//! watts as in the monochromatic query. This is the usual cycle-averaged
//! photodetector law, not a subcycle measurement of the real electric field:
//! optical sum-frequency terms are outside this scalar model. Difference
//! frequencies survive through the shutter profile's analytic kernels and
//! bounded positive modal factor, without a time raster.
//!
//! Monochromatic intensity can omit the sensor's common Fresnel phase; a line
//! spectrum cannot. For x=d*theta, each line contributes
//! `h_f(theta) exp(i k_f d (1 + |theta|²/2) - i 2 pi f t)`. Pupil fields and
//! their relative source phases retain the original entrance-pupil reference.
//! A common frequency-independent phase may be dropped, but the wavelength-
//! dependent longitudinal piston and quadratic sensor chirp are retained.
//! Source-to-pupil and any additional lens/material piston remain authored.

use super::{
    CoherentDiffractionError, CoherentDiffractionLimits, CoherentDiffractionPixel,
    CoherentDiffractionPixelLimits, CoherentDiffractionPixelPower, CoherentPupilIllumination, Dd,
    PhysicalCamera, PreparedCoherentDiffraction, admit_coordinate, admit_count, pixel, poll,
    rounded, square,
};
use fs_fft::harmonic_window::{
    HarmonicWindow, HarmonicWindowProfile, PreparedHarmonicBand, ProperHarmonicProduct,
};
use fs_math::det;

mod interferometer;
pub use interferometer::CoherentInterferometerPixelExposure;
pub use interferometer::{
    CoherentInterferometerExposureLimits, CoherentInterferometerExposureSample,
};

mod jones;
pub use jones::PolarizedCoherentInterferometerExposureSample;
pub use jones::PolarizedCoherentInterferometerPixelExposure;
pub use jones::{
    CoherentJonesDetector, CoherentJonesDiffractionTone, CoherentJonesPupilIllumination,
    CoherentJonesPupilLayer, CoherentJonesPupilLimits, CoherentJonesPupilResponse,
    CoherentJonesPupilWave, JonesPupilCoherenceFactor, MAX_COHERENT_JONES_PUPIL_MODES,
    PolarizedCoherentDiffractionExposureSample, PolarizedCoherentDiffractionPixelExposure,
    PreparedCoherentJonesDiffractionExposure,
};
pub use jones::{CoherentJonesInterferometerResponse, PreparedCoherentJonesInterferometer};

const C_M_S: f64 = 299_792_458.0;
const TAU: Dd = Dd {
    hi: core::f64::consts::TAU,
    lo: 2.449_293_598_294_706_4e-16,
};

/// Hard spectral-line bound; observation scratch is a fixed array of this size.
/// Caller budgets may impose a smaller bound. Source waves remain borrowed.
pub const MAX_COHERENT_DIFFRACTION_TONES: usize = 64;

/// One positive-frequency line on the common entrance-pupil clock.
#[derive(Clone, Copy, Debug)]
pub struct CoherentDiffractionTone<'a> {
    /// Vacuum frequency in Hz; the original pupil receives wavelength c/f.
    /// Use a shared clock and exp(-i 2*pi*f*t). Equal frequencies are merged
    /// as fields before moments, even when authored as separate tone entries.
    pub frequency_hz: f64,
    /// Mean plus covariance columns. Column r denotes the SAME random source
    /// coefficient at every frequency; separate columns are uncorrelated.
    /// All tones must have the same column count; an absent line uses an empty
    /// wave slice in that column. No coherence is inferred from spectral proximity.
    pub illumination: CoherentPupilIllumination<'a>,
}

/// Physical shutter with unit peak transmission and an explicit temporal profile.
///
/// Uniform transmission is one throughout the interval. Triangular transmission
/// rises linearly from zero to one at its midpoint and falls to zero. Power
/// observations are normalized transmission-weighted averages; collected energy
/// multiplies them by `effective_duration_s()`, which is T or T/2 respectively.
/// Start and duration remain separately represented on the original source clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDiffractionShutter(HarmonicWindow);

impl CoherentDiffractionShutter {
    /// The start/end must be finite and the duration finite and positive.
    pub fn try_new(start_s: f64, duration_s: f64) -> Result<Self, CoherentDiffractionError> {
        Self::try_with_profile(start_s, duration_s, HarmonicWindowProfile::Uniform)
    }

    /// Admit the normalized observation profile and its physical unit-peak
    /// transmission. A shaped shutter does not transmit for the full duration
    /// at unit gain; its integrated transmission remains explicit.
    pub fn try_with_profile(
        start_s: f64,
        duration_s: f64,
        profile: HarmonicWindowProfile,
    ) -> Result<Self, CoherentDiffractionError> {
        Self::from_window(HarmonicWindow::with_profile(start_s, duration_s, profile)?)
    }

    /// Retain the entrance-pupil source epoch below one binary64 clock step.
    /// The original optical beats and shutter duration keep their own phases.
    pub fn try_from_clock_dd(
        start_s: Dd,
        duration_s: f64,
    ) -> Result<Self, CoherentDiffractionError> {
        Self::try_with_profile_from_clock_dd(start_s, duration_s, HarmonicWindowProfile::Uniform)
    }

    /// Retain the precise source epoch through a shaped physical shutter.
    pub fn try_with_profile_from_clock_dd(
        start_s: Dd,
        duration_s: f64,
        profile: HarmonicWindowProfile,
    ) -> Result<Self, CoherentDiffractionError> {
        Self::from_window(HarmonicWindow::with_profile_from_clock_dd(
            start_s, duration_s, profile,
        )?)
    }

    fn from_window(window: HarmonicWindow) -> Result<Self, CoherentDiffractionError> {
        let shutter = Self(window);
        if shutter.effective_duration_s() == 0.0 {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
        Ok(shutter)
    }

    #[must_use]
    pub fn start_s(self) -> f64 {
        self.0.start_s()
    }

    #[must_use]
    pub fn start_clock_dd(self) -> Dd {
        self.0.start_clock_dd()
    }

    #[must_use]
    pub fn duration_s(self) -> f64 {
        self.0.duration_s()
    }

    #[must_use]
    pub fn profile(self) -> HarmonicWindowProfile {
        self.0.profile()
    }

    /// Integral of the physical transmission over time, in seconds. This is
    /// the exposure multiplier for normalized weighted powers and coherencies.
    #[must_use]
    pub fn effective_duration_s(self) -> f64 {
        self.duration_s()
            * match self.profile() {
                HarmonicWindowProfile::Uniform => 1.0,
                HarmonicWindowProfile::Triangular => 0.5,
            }
    }
}

/// Preparation bounds, including the cost of one spatial observation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoherentDiffractionExposureLimits {
    pub maximum_tones: usize,
    /// Global authored wave count and shared covariance-column count.
    pub sources: CoherentDiffractionLimits,
    /// Temporal projection/product work across mean and all factors. A narrow
    /// band with N frequencies and M modes costs M*(N+1) per source bank;
    /// the wider-band pair contraction costs N*(N+1)/2. Checked before the
    /// temporal plan or pupil transfers are allocated.
    pub maximum_pair_contractions: usize,
}

/// Spatial and temporal work bounds for one collected pixel exposure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDiffractionExposurePixelLimits {
    pub spatial: CoherentDiffractionPixelLimits,
    /// Total frequency contractions over every spatial quadrature sample.
    pub maximum_pair_contractions: usize,
}

/// Transmission-weighted mean power per tangent-coordinate area, normalized
/// by the shutter's integrated transmission rather than its elapsed duration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDiffractionExposureSample {
    pub mean_power_density_w: f64,
    pub fluctuation_power_density_w: f64,
    pub power_density_w: f64,
}

/// Physical shutter energy collected over the complete finite rectangular pixel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDiffractionPixelExposure {
    pub average_power: CoherentDiffractionPixelPower,
    pub mean_energy_j: f64,
    pub fluctuation_energy_j: f64,
    pub energy_j: f64,
    pub pair_contractions: usize,
}

struct FrequencyGroup {
    frequency_hz: f64,
    // Keep the authored frequency and the original phase index through
    // propagation. A rounded c/f wavelength cannot preserve long-path beats.
    optical_frequency_hz: Dd,
    relative_optical_frequency_hz: Dd,
}

/// Cached original pupil transfers and analytic shutter kernels.
///
/// Equal frequencies are reduced in DD before squaring. Different frequencies
/// use analytic difference-only kernels. Narrow bands use a positive modal
/// factor of the shutter profile, retaining weak multi-tone interference
/// without subtracting nearly equal powers. The bounded modal truncation is
/// documented by `fs_fft::harmonic_window::PreparedHarmonicBand`; it is
/// not a relative error certificate at an arbitrary null. Observation allocates
/// nothing and retains bounded cancellation points. Unrepresentable contracted
/// powers are refused, never clamped to a plausible positive measurement.
pub struct PreparedCoherentDiffractionExposure<'camera, 'source> {
    tones: Vec<PreparedCoherentDiffraction<'camera, 'source>>,
    tone_groups: Vec<usize>,
    groups: Vec<FrequencyGroup>,
    kernels: Vec<ProperHarmonicProduct>,
    band: Option<PreparedHarmonicBand>,
    shutter: CoherentDiffractionShutter,
    sensor_distance_m: f64,
    minimum_scale: f64,
    wavenumber_span: f64,
    factors: usize,
    waves: usize,
    contractions: usize,
}

impl PhysicalCamera {
    /// Prepare finite-exposure scalar photodetection of a coherent line spectrum.
    ///
    /// `sensor_distance_m` is the pupil-to-sensor axial distance. For a physical
    /// focal-length camera it must agree with the original lens conjugate; for
    /// an angular/FOV camera it explicitly fixes the otherwise undeclared
    /// sensor plane (and ideal lens power at the declared focus distance).
    /// Every line retains original wavelength/material/defocus admission.
    /// All source phases are at the pupil centre at t=0, not at the detector.
    pub fn prepare_coherent_diffraction_exposure<'camera, 'source>(
        &'camera self,
        tones: &[CoherentDiffractionTone<'source>],
        sensor_distance_m: f64,
        shutter: CoherentDiffractionShutter,
        limits: CoherentDiffractionExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCoherentDiffractionExposure<'camera, 'source>, CoherentDiffractionError>
    {
        poll(&mut cancelled)?;
        admit_count(
            "spectral tones",
            tones.len(),
            limits.maximum_tones.min(MAX_COHERENT_DIFFRACTION_TONES),
        )?;
        admit_sensor(self, sensor_distance_m)?;
        let factors = tones
            .first()
            .map_or(0, |tone| tone.illumination.factors.len());
        admit_count("source factors", factors, limits.sources.maximum_factors)?;
        let mut waves = 0usize;
        let mut frequencies = [0.0; MAX_COHERENT_DIFFRACTION_TONES];
        let mut group_count = 0;
        let mut tone_groups = Vec::with_capacity(tones.len());
        for tone in tones {
            poll(&mut cancelled)?;
            if !tone.frequency_hz.is_finite() || tone.frequency_hz <= 0.0 {
                return Err(CoherentDiffractionError::InvalidFrequency);
            }
            if tone.illumination.factors.len() != factors {
                return Err(CoherentDiffractionError::InconsistentCoherenceFactors);
            }
            for bank in core::iter::once(tone.illumination.mean)
                .chain(tone.illumination.factors.iter().map(|factor| factor.waves))
            {
                poll(&mut cancelled)?;
                waves = add_count(
                    waves,
                    bank.len(),
                    "pupil waves",
                    limits.sources.maximum_waves,
                )?;
            }
            #[allow(clippy::float_cmp)]
            // Exact authored frequency identity defines a coherent line.
            let group = frequencies[..group_count]
                .iter()
                .position(|&f| f == tone.frequency_hz);
            let group = group.unwrap_or_else(|| {
                frequencies[group_count] = tone.frequency_hz;
                group_count += 1;
                group_count - 1
            });
            tone_groups.push(group);
        }
        let band_modes = shutter.0.band_mode_count_hz(&frequencies[..group_count])?;
        let component_work = band_modes.map_or(group_count * (group_count + 1) / 2, |modes| {
            modes * (group_count + 1)
        });
        let contractions = component_work.checked_mul(factors + 1).ok_or_else(|| {
            limit(
                "exposure pair contractions",
                limits.maximum_pair_contractions,
            )
        })?;
        admit_count(
            "exposure pair contractions",
            contractions,
            limits.maximum_pair_contractions,
        )?;
        let band = shutter.0.prepare_band_hz(&frequencies[..group_count])?;
        let mut prepared = Vec::with_capacity(tones.len());
        let mut groups = Vec::with_capacity(group_count);
        let mut minimum_scale = f64::INFINITY;
        let mut reference_frequency = Dd::ZERO;
        let mut minimum_k = 0.0_f64;
        let mut maximum_k = 0.0_f64;
        for (index, tone) in tones.iter().enumerate() {
            poll(&mut cancelled)?;
            let wavelength_nm = rounded(
                Dd::from_f64(C_M_S) / Dd::from_f64(tone.frequency_hz) * Dd::from_f64(1.0e9),
            )?;
            let pupil = self.prepare_coherent_diffraction(
                wavelength_nm,
                tone.illumination,
                limits.sources,
                &mut cancelled,
            )?;
            minimum_scale = minimum_scale.min(pupil.pupil.scale);
            if tone_groups[index] == groups.len() {
                // The phase index belongs to the original admitted pupil medium.
                // Vacuum keeps n exactly one. Retain f and subtraction in DD so
                // close carrier frequencies do not acquire false propagation beats.
                #[cfg(feature = "tracer")]
                let phase_index = self
                    .dispersive_lens
                    .as_ref()
                    .map_or(Ok(1.0), |lens| lens.exterior_phase_index(wavelength_nm))?;
                #[cfg(not(feature = "tracer"))]
                let phase_index = 1.0;
                let optical_frequency_hz =
                    Dd::from_f64(tone.frequency_hz) * Dd::from_f64(phase_index);
                rounded(optical_frequency_hz)?;
                if groups.is_empty() {
                    reference_frequency = optical_frequency_hz;
                }
                let relative_optical_frequency_hz = optical_frequency_hz - reference_frequency;
                // Radian conversion follows the frequency subtraction. This
                // wavenumber is only the spatial quadrature bandwidth; physical
                // phases reduce frequency-length cycles before multiplying TAU.
                let relative_wavenumber = TAU * relative_optical_frequency_hz / Dd::from_f64(C_M_S);
                let relative = rounded(relative_wavenumber)?;
                minimum_k = minimum_k.min(relative);
                maximum_k = maximum_k.max(relative);
                groups.push(FrequencyGroup {
                    frequency_hz: tone.frequency_hz,
                    optical_frequency_hz,
                    relative_optical_frequency_hz,
                });
            }
            prepared.push(pupil);
        }
        let mut kernels = Vec::new();
        if band.is_none() {
            kernels.reserve_exact(group_count * group_count.saturating_sub(1) / 2);
            for (i, group) in groups.iter().enumerate() {
                for previous in &groups[..i] {
                    poll(&mut cancelled)?;
                    kernels.push(
                        shutter
                            .0
                            .proper_product_kernel_hz(group.frequency_hz, previous.frequency_hz)?,
                    );
                }
            }
        }
        let wavenumber_span = maximum_k - minimum_k;
        if !wavenumber_span.is_finite() {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
        poll(&mut cancelled)?;
        Ok(PreparedCoherentDiffractionExposure {
            tones: prepared,
            tone_groups,
            groups,
            kernels,
            band,
            shutter,
            sensor_distance_m,
            minimum_scale,
            wavenumber_span,
            factors,
            waves,
            contractions,
        })
    }
}

impl PreparedCoherentDiffractionExposure<'_, '_> {
    // Per-source-bank work from the same temporal representation used by all
    // scalar, polarized, and mutual-coherence consumers. N is hard-bounded.
    fn power_component_work(&self) -> usize {
        self.band.as_ref().map_or_else(
            || self.groups.len() * (self.groups.len() + 1) / 2,
            PreparedHarmonicBand::power_work,
        )
    }

    fn cross_component_work(&self) -> usize {
        self.band.as_ref().map_or_else(
            || self.groups.len() * (2 * self.groups.len()).saturating_sub(1),
            PreparedHarmonicBand::cross_work,
        )
    }

    fn coherency_component_work(&self) -> usize {
        self.band.as_ref().map_or_else(
            || self.groups.len() * (3 * self.groups.len() + 1),
            PreparedHarmonicBand::coherency_work,
        )
    }

    #[must_use]
    pub const fn shutter(&self) -> CoherentDiffractionShutter {
        self.shutter
    }

    #[must_use]
    pub fn frequency_count(&self) -> usize {
        self.groups.len()
    }

    /// Analytic shutter average at one admitted observation coordinate.
    pub fn observe(
        &self,
        at: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentDiffractionExposureSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        admit_coordinate(at)?;
        let [mean, fluctuations] = self.moments(at, &mut cancelled)?;
        let result = CoherentDiffractionExposureSample {
            mean_power_density_w: rounded(mean)?,
            fluctuation_power_density_w: rounded(fluctuations)?,
            power_density_w: rounded(mean + fluctuations)?,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Integrate the original spatial rule and analytic shutter-weighted beats.
    ///
    /// Spatial resolution includes both the shortest physical wavelength's
    /// pupil bandwidth and the mixed-frequency sensor-chirp gradient. Time is
    /// integrated analytically; there is no carrier/beat time raster. All work
    /// bounds are checked before any pupil evaluation. Refine the spatial phase
    /// setting to assess quadrature error; it is not an error certificate.
    pub fn integrate_pixel(
        &self,
        pixel: CoherentDiffractionPixel,
        limits: CoherentDiffractionExposurePixelLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentDiffractionPixelExposure, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let (subdivisions, samples, pupil_evaluations, pair_contractions) =
            self.pixel_work(pixel, limits, self.waves, self.contractions)?;
        let power =
            pixel::integrate_moments(pixel, subdivisions, &mut cancelled, |at, cancelled| {
                self.moments(at, cancelled)
            })?;
        let average_power = CoherentDiffractionPixelPower {
            mean_power_w: rounded(power[0])?,
            fluctuation_power_w: rounded(power[1])?,
            power_w: rounded(power[0] + power[1])?,
            subdivisions,
            quadrature_samples: samples,
            pupil_evaluations,
        };
        let duration = Dd::from_f64(self.shutter.effective_duration_s());
        let result = CoherentDiffractionPixelExposure {
            average_power,
            mean_energy_j: rounded(power[0] * duration)?,
            fluctuation_energy_j: rounded(power[1] * duration)?,
            energy_j: rounded((power[0] + power[1]) * duration)?,
            pair_contractions,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    // The Jones detector shares the same pupil/sensor bandwidth and positive
    // spatial rule; only its source count and number of analyzed components differ.
    fn pixel_work(
        &self,
        pixel: CoherentDiffractionPixel,
        limits: CoherentDiffractionExposurePixelLimits,
        waves: usize,
        contractions: usize,
    ) -> Result<([usize; 2], usize, usize, usize), CoherentDiffractionError> {
        let minimum = pixel.minimum_tangent();
        let maximum = pixel.maximum_tangent();
        let maximum_step = std::array::from_fn(|axis| {
            let radius = minimum[axis].abs().max(maximum[axis].abs());
            let gradient =
                2.0 / self.minimum_scale + self.wavenumber_span * self.sensor_distance_m * radius;
            limits.spatial.maximum_phase_span_rad / gradient
        });
        let (subdivisions, samples, pupil_evaluations) =
            pixel::pixel_work(pixel, limits.spatial, maximum_step, waves)?;
        let pair_contractions = samples.checked_mul(contractions).ok_or_else(|| {
            limit(
                "pixel exposure pair contractions",
                limits.maximum_pair_contractions,
            )
        })?;
        admit_count(
            "pixel exposure pair contractions",
            pair_contractions,
            limits.maximum_pair_contractions,
        )?;
        Ok((subdivisions, samples, pupil_evaluations, pair_contractions))
    }

    fn rotations(
        &self,
        at: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[[Dd; 2]; MAX_COHERENT_DIFFRACTION_TONES], CoherentDiffractionError> {
        let distance = Dd::from_f64(self.sensor_distance_m)
            * (Dd::ONE
                + (Dd::from_f64(at[0]) * Dd::from_f64(at[0])
                    + Dd::from_f64(at[1]) * Dd::from_f64(at[1]))
                    * Dd::from_f64(0.5));
        rounded(distance)?;
        let mut rotations = [[Dd::ZERO; 2]; MAX_COHERENT_DIFFRACTION_TONES];
        for (i, group) in self.groups.iter().enumerate() {
            poll(cancelled)?;
            let phase = interferometer::frequency_length_phase(
                group.relative_optical_frequency_hz,
                distance,
            )?;
            let (sine, cosine) = det::sin_cos_dd(phase);
            rotations[i] = [cosine, sine];
        }
        Ok(rotations)
    }

    // One scalar analytic-field component. Jones projection happens BEFORE
    // this contraction, so weak analyzer residuals never subtract rounded powers.
    fn contract_component(
        &self,
        fields: &mut [[Dd; 2]; MAX_COHERENT_DIFFRACTION_TONES],
        rotations: &[[Dd; 2]; MAX_COHERENT_DIFFRACTION_TONES],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Dd, CoherentDiffractionError> {
        if let Some(band) = &self.band {
            for i in 0..self.groups.len() {
                poll(cancelled)?;
                let [ar, ai] = fields[i];
                let [cr, ci] = rotations[i];
                fields[i] = [ar * cr - ai * ci, ar * ci + ai * cr];
            }
            let moment = band.mean_square_parts_dd(&fields[..self.groups.len()])?;
            if rounded(moment)? < 0.0 {
                return Err(CoherentDiffractionError::Unrepresentable);
            }
            poll(cancelled)?;
            return Ok(moment);
        }
        let mut moment = Dd::ZERO;
        for i in 0..self.groups.len() {
            poll(cancelled)?;
            let [ar, ai] = fields[i];
            let [cr, ci] = rotations[i];
            fields[i] = [ar * cr - ai * ci, ar * ci + ai * cr];
            moment = moment + square(fields[i])?;
            for j in 0..i {
                poll(cancelled)?;
                // The shared kernel includes 1/2 for real peak phasors.
                // Analytic-field power removes that half, and the two
                // Hermitian cross terms supply another factor of two.
                moment = moment
                    + Dd::from_f64(4.0)
                        * self.kernels[i * (i - 1) / 2 + j]
                            .mean_real_product_parts_dd(fields[i], fields[j])?;
            }
        }
        if rounded(moment)? < 0.0 {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
        Ok(moment)
    }

    fn moments(
        &self,
        at: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 2], CoherentDiffractionError> {
        if self.waves == 0 {
            return Ok([Dd::ZERO; 2]);
        }
        let rotations = self.rotations(at, cancelled)?;
        let mut result = [Dd::ZERO; 2];
        for bank in 0..=self.factors {
            poll(cancelled)?;
            let mut fields = [[Dd::ZERO; 2]; MAX_COHERENT_DIFFRACTION_TONES];
            for (tone, &group) in self.tones.iter().zip(&self.tone_groups) {
                poll(cancelled)?;
                let waves = if bank == 0 {
                    tone.illumination.mean
                } else {
                    tone.illumination.factors[bank - 1].waves
                };
                let field = tone.field(waves, at, cancelled)?;
                fields[group][0] = fields[group][0] + field[0];
                fields[group][1] = fields[group][1] + field[1];
            }
            let moment = self.contract_component(&mut fields, &rotations, cancelled)?;
            let component = usize::from(bank != 0);
            result[component] = result[component] + moment;
        }
        Ok(result)
    }
}

pub(super) fn admit_sensor(
    camera: &PhysicalCamera,
    distance: f64,
) -> Result<(), CoherentDiffractionError> {
    if !distance.is_finite() || distance <= 0.0 {
        return Err(CoherentDiffractionError::InvalidSensorPlane);
    }
    let physical = camera
        .projection()
        .focal_length_m()
        .map(|f| f / (1.0 - f / camera.focus_distance_m()));
    #[cfg(feature = "tracer")]
    let physical = camera.dispersive_sensor_distance_m().or(physical);
    if physical.is_some_and(|expected| {
        !expected.is_finite()
            || expected <= 0.0
            || (distance - expected).abs() > 8.0 * f64::EPSILON * expected.max(distance)
    }) {
        return Err(CoherentDiffractionError::InvalidSensorPlane);
    }
    Ok(())
}

fn limit(quantity: &'static str, maximum: usize) -> CoherentDiffractionError {
    CoherentDiffractionError::Limit {
        quantity,
        requested: usize::MAX,
        maximum,
    }
}

fn add_count(
    value: usize,
    additional: usize,
    quantity: &'static str,
    maximum: usize,
) -> Result<usize, CoherentDiffractionError> {
    let total = value
        .checked_add(additional)
        .ok_or_else(|| limit(quantity, maximum))?;
    admit_count(quantity, total, maximum)?;
    Ok(total)
}
