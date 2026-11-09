//! Polarization-resolved interference through the original paraxial camera pupil.
//!
//! Jones components use one common orthonormal pupil x/y basis. A detector
//! insensitive to polarization measures |E_x|² + |E_y|², never |E_x + E_y|².
//! An ideal complex analyzer measures |a† E|² for a normalized Jones mode a.
//! The projection precedes temporal and spatial power integration, retaining
//! cross-polarization phase without subtracting already rounded intensities.
//!
//! This is two-component PARAXIAL diffraction through the original scalar
//! pupil, optionally carrying passive patterned Jones layers. It does not model
//! longitudinal fields, high-NA vector diffraction, thick anisotropic lenses or
//! transport between nonparallel ray frames.
//! See TU Delft, *Polarisation*, for the Jones-vector and analyzer conventions:
//! <https://interactivetextbooks.tudelft.nl/interactive-optics/content/Chap4_Polarisation/Polarization_2022_01Clean.html>.

use super::{
    CoherentDiffractionError, CoherentDiffractionExposureLimits,
    CoherentDiffractionExposurePixelLimits, CoherentDiffractionExposureSample,
    CoherentDiffractionLimits, CoherentDiffractionPixel, CoherentDiffractionPixelExposure,
    CoherentDiffractionPixelPower, CoherentDiffractionShutter, CoherentDiffractionTone,
    CoherentPupilIllumination, Dd, MAX_COHERENT_DIFFRACTION_TONES, PhysicalCamera,
    PreparedCoherentDiffractionExposure, add_count, admit_coordinate, admit_count, limit, pixel,
    poll, rounded,
};
use fs_math::c64::C64;

mod pupil;
use pupil::PreparedJonesPupil;
pub use pupil::{
    CoherentJonesPupilLayer, CoherentJonesPupilLimits, CoherentJonesPupilResponse,
    MAX_COHERENT_JONES_PUPIL_MODES,
};

mod coherency;
pub use coherency::{
    PolarizedCoherentDiffractionExposureSample, PolarizedCoherentDiffractionPixelExposure,
};


mod interferometer;
pub use interferometer::{
    CoherentJonesInterferometerResponse, PreparedCoherentJonesInterferometer,
    PolarizedCoherentInterferometerExposureSample, PolarizedCoherentInterferometerPixelExposure,
};

type JonesFields = [[[Dd; 2]; 2]; MAX_COHERENT_DIFFRACTION_TONES];

/// One tilted Jones field in the camera's common paraxial pupil x/y basis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentJonesPupilWave {
    /// Same dimensionless tangent coordinates and domain as the scalar query.
    pub tangent_direction: [f64; 2],
    /// Complex x/y amplitudes in sqrt(W), with the common pupil-centre phase
    /// and clock. |a_x|² + |a_y|² is this wave's clear-pupil power in isolation,
    /// before any installed Jones-pupil layers; passive loss is not renormalized.
    /// Relative polarization phase is physical and must not be discarded.
    pub amplitude_sqrt_w: [C64; 2],
}

/// One shared random coefficient multiplying all Jones waves in this column.
#[derive(Clone, Copy, Debug)]
pub struct JonesPupilCoherenceFactor<'a> {
    /// Both polarization components share this coefficient; different columns
    /// are uncorrelated. Independent x/y fluctuations need different columns.
    pub waves: &'a [CoherentJonesPupilWave],
}

/// Deterministic Jones mean and an explicitly authored covariance factorization.
#[derive(Clone, Copy, Debug)]
pub struct CoherentJonesPupilIllumination<'a> {
    pub mean: &'a [CoherentJonesPupilWave],
    /// Column identities are shared across frequencies AND polarizations.
    pub factors: &'a [JonesPupilCoherenceFactor<'a>],
}

/// One positive-frequency Jones line on the original exp(-i 2*pi*f*t) clock.
#[derive(Clone, Copy, Debug)]
pub struct CoherentJonesDiffractionTone<'a> {
    pub frequency_hz: f64,
    pub illumination: CoherentJonesPupilIllumination<'a>,
}

/// Ideal, frequency-independent polarization sensitivity on the sensor plane.
///
/// An analyzer selects one unit Jones mode with amplitude a† E. Its orthogonal
/// output is rejected; it is not another coherent input to this detector.
/// Changing the analyzer does not change source coherence or the pupil cache.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentJonesDetector {
    analyzer: Option<[[Dd; 2]; 2]>,
}

impl CoherentJonesDetector {
    /// Collect both orthogonal polarization powers without mixing their fields.
    pub const POLARIZATION_INSENSITIVE: Self = Self { analyzer: None };

    /// Select the supplied nonzero Jones direction, normalized in DD.
    ///
    /// Magnitude is discarded deliberately: this argument defines a MODE,
    /// not an attenuator or gain. [1, 1] selects +45-degree linear polarization;
    /// [1, i] selects that explicitly phased circular mode without a handedness
    /// convention inferred from viewing direction. Components must be finite.
    pub fn analyzer(mode: [C64; 2]) -> Result<Self, CoherentDiffractionError> {
        let values = [[mode[0].re, mode[0].im], [mode[1].re, mode[1].im]];
        let mut scale = 0.0_f64;
        for value in values.into_iter().flatten() {
            if !value.is_finite() {
                return Err(CoherentDiffractionError::InvalidPolarizationAnalyzer);
            }
            scale = scale.max(value.abs());
        }
        if scale == 0.0 {
            return Err(CoherentDiffractionError::InvalidPolarizationAnalyzer);
        }
        let mut axis = values.map(|z| z.map(|v| Dd::from_f64(v) / Dd::from_f64(scale)));
        let norm = axis
            .iter()
            .flatten()
            .fold(Dd::ZERO, |sum, &v| sum + v * v)
            .sqrt();
        for value in axis.iter_mut().flatten() {
            *value = *value / norm;
            rounded(*value)?;
        }
        Ok(Self {
            analyzer: Some(axis),
        })
    }

    fn components(self) -> usize {
        if self.analyzer.is_some() { 1 } else { 2 }
    }
}

/// Reusable Jones source spectrum through the original scalar pupil optics.
///
/// Each unmasked wave evaluates its original chromatic pupil once per spatial
/// sample, updating both DD components. Patterned pupils use one evaluation per
/// compiled Fourier mode, with explicit total work bounds. Observation allocates nothing.
/// The common transverse basis and scalar pupil are a paraxial model boundary,
/// not a claim of vector Maxwell transport or a relative-error certificate.
pub struct PreparedCoherentJonesDiffractionExposure<'camera, 'source> {
    // Empty scalar source banks retain the original optical cache, frequency
    // grouping, sensor phases and analytic shutter kernels without duplicating
    // their physical implementation or borrowing temporary projected sources.
    optics: PreparedCoherentDiffractionExposure<'camera, 'static>,
    sources: Vec<CoherentJonesPupilIllumination<'source>>,
    factors: usize,
    pupil_evaluations: usize,
    pupil_radius_m: f64,
    pupils: Option<Vec<PreparedJonesPupil>>,
    contractions_per_component: usize,
    maximum_pair_contractions: usize,
}

impl PhysicalCamera {
    /// Prepare coherent Jones illumination for finite-time photodetection.
    ///
    /// One tone is the monochromatic case. Each line retains the original
    /// wavelength, material, defocus, sensor-conjugate and paraxial admission.
    /// `sources.maximum_waves` counts vector waves, not their two components.
    /// Preparation's pair budget covers BOTH detector components so the same
    /// immutable preparation supports either detector. Pixel work charges the
    /// actual one or two analyzed components. No source phase or coherence is
    /// inferred from polarization, frequency proximity or scene radiance.
    /// Full coherency queries additionally preflight their four matrix moments
    /// against this same pair budget before evaluating any pupil.
    pub fn prepare_coherent_jones_diffraction_exposure<'camera, 'source>(
        &'camera self,
        tones: &[CoherentJonesDiffractionTone<'source>],
        sensor_distance_m: f64,
        shutter: CoherentDiffractionShutter,
        limits: CoherentDiffractionExposureLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCoherentJonesDiffractionExposure<'camera, 'source>, CoherentDiffractionError>
    {
        poll(&mut cancelled)?;
        admit_count(
            "spectral tones",
            tones.len(),
            limits.maximum_tones.min(MAX_COHERENT_DIFFRACTION_TONES),
        )?;
        let factors = tones
            .first()
            .map_or(0, |tone| tone.illumination.factors.len());
        admit_count("source factors", factors, limits.sources.maximum_factors)?;
        let mut waves = 0;
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
        }
        let component_banks = factors
            .checked_add(1)
            .and_then(|banks| banks.checked_mul(2))
            .ok_or_else(|| {
                limit(
                    "Jones exposure pair contractions",
                    limits.maximum_pair_contractions,
                )
            })?;
        for tone in tones {
            for bank in core::iter::once(tone.illumination.mean)
                .chain(tone.illumination.factors.iter().map(|factor| factor.waves))
            {
                poll(&mut cancelled)?;
                for wave in bank {
                    poll(&mut cancelled)?;
                    admit_coordinate(wave.tangent_direction)?;
                    if wave
                        .amplitude_sqrt_w
                        .iter()
                        .any(|a| !a.re.is_finite() || !a.im.is_finite())
                    {
                        return Err(CoherentDiffractionError::NonFiniteAmplitude);
                    }
                }
            }
        }
        let empty: Vec<CoherentDiffractionTone<'static>> = tones
            .iter()
            .map(|tone| CoherentDiffractionTone {
                frequency_hz: tone.frequency_hz,
                illumination: CoherentPupilIllumination {
                    mean: &[],
                    factors: &[],
                },
            })
            .collect();
        // Dividing the budget lets the original frequency grouping preflight
        // the full Jones contraction count BEFORE preparing any pupil/kernel.
        let optics = self
            .prepare_coherent_diffraction_exposure(
                &empty,
                sensor_distance_m,
                shutter,
                CoherentDiffractionExposureLimits {
                    maximum_tones: limits.maximum_tones,
                    sources: CoherentDiffractionLimits {
                        maximum_waves: 0,
                        maximum_factors: 0,
                    },
                    maximum_pair_contractions: limits.maximum_pair_contractions / component_banks,
                },
                &mut cancelled,
            )
            .map_err(|error| match error {
                CoherentDiffractionError::Limit {
                    quantity: "exposure pair contractions",
                    requested,
                    ..
                } => CoherentDiffractionError::Limit {
                    quantity: "Jones exposure pair contractions",
                    requested: requested.saturating_mul(component_banks),
                    maximum: limits.maximum_pair_contractions,
                },
                other => other,
            })?;
        let contractions_per_component =
            optics
                .contractions
                .checked_mul(factors + 1)
                .ok_or_else(|| {
                    limit(
                        "Jones exposure pair contractions",
                        limits.maximum_pair_contractions,
                    )
                })?;
        poll(&mut cancelled)?;
        Ok(PreparedCoherentJonesDiffractionExposure {
            optics,
            sources: tones.iter().map(|tone| tone.illumination).collect(),
            factors,
            pupil_evaluations: waves,
            pupil_radius_m: self.aperture().radius_m(),
            pupils: None,
            contractions_per_component,
            maximum_pair_contractions: limits.maximum_pair_contractions,
        })
    }
}

impl PreparedCoherentJonesDiffractionExposure<'_, '_> {
    #[must_use]
    pub const fn shutter(&self) -> CoherentDiffractionShutter {
        self.optics.shutter()
    }

    #[must_use]
    pub fn frequency_count(&self) -> usize {
        self.optics.frequency_count()
    }

    /// Shutter-average power density in the original tangent coordinates.
    /// Orthogonal components add power; an analyzer projects the coherent
    /// Jones fields before the shared analytic temporal contraction.
    pub fn observe(
        &self,
        at: [f64; 2],
        detector: CoherentJonesDetector,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentDiffractionExposureSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        admit_coordinate(at)?;
        let [mean, fluctuations] = self.moments(at, detector, &mut cancelled)?;
        let result = CoherentDiffractionExposureSample {
            mean_power_density_w: rounded(mean)?,
            fluctuation_power_density_w: rounded(fluctuations)?,
            power_density_w: rounded(mean + fluctuations)?,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Collect polarized average watts and shutter-integrated joules over a
    /// finite pixel. Uses the ORIGINAL positive Gauss rule, physical pupil
    /// bandwidth and cross-frequency sensor-chirp gradient. Resolution is not
    /// a quadrature certificate. No field or polarization state is assigned to
    /// the integrated power. All work is preflighted before pupil evaluation.
    pub fn integrate_pixel(
        &self,
        pixel: CoherentDiffractionPixel,
        detector: CoherentJonesDetector,
        limits: CoherentDiffractionExposurePixelLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentDiffractionPixelExposure, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let contractions = self.contractions_per_component * detector.components();
        let (subdivisions, samples, pupil_evaluations, pair_contractions) = self
            .optics
            .pixel_work(pixel, limits, self.pupil_evaluations, contractions)?;
        let power =
            pixel::integrate_moments(pixel, subdivisions, &mut cancelled, |at, cancelled| {
                self.moments(at, detector, cancelled)
            })?;
        let duration = Dd::from_f64(self.shutter().effective_duration_s());
        let result = CoherentDiffractionPixelExposure {
            average_power: CoherentDiffractionPixelPower {
                mean_power_w: rounded(power[0])?,
                fluctuation_power_w: rounded(power[1])?,
                power_w: rounded(power[0] + power[1])?,
                subdivisions,
                quadrature_samples: samples,
                pupil_evaluations,
            },
            mean_energy_j: rounded(power[0] * duration)?,
            fluctuation_energy_j: rounded(power[1] * duration)?,
            energy_j: rounded((power[0] + power[1]) * duration)?,
            pair_contractions,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn moments(
        &self,
        at: [f64; 2],
        detector: CoherentJonesDetector,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 2], CoherentDiffractionError> {
        if self.pupil_evaluations == 0 {
            return Ok([Dd::ZERO; 2]);
        }
        let rotations = self.optics.rotations(at, cancelled)?;
        let mut result = [Dd::ZERO; 2];
        for bank in 0..=self.factors {
            poll(cancelled)?;
            let fields = self.group_fields(bank, at, cancelled)?;
            let mut moment = Dd::ZERO;
            for component in 0..detector.components() {
                poll(cancelled)?;
                let mut analyzed = [[Dd::ZERO; 2]; MAX_COHERENT_DIFFRACTION_TONES];
                for (target, field) in analyzed
                    .iter_mut()
                    .zip(&fields)
                    .take(self.frequency_count())
                {
                    for value in field.iter().flatten() {
                        rounded(*value)?;
                    }
                    if let Some(axis) = detector.analyzer {
                        for (a, e) in axis.iter().zip(field) {
                            target[0] = target[0] + (a[0] * e[0] + a[1] * e[1]);
                            target[1] = target[1] + (a[0] * e[1] - a[1] * e[0]);
                        }
                    } else {
                        *target = field[component];
                    }
                }
                moment = moment
                    + self
                        .optics
                        .contract_component(&mut analyzed, &rotations, cancelled)?;
            }
            let component = usize::from(bank != 0);
            result[component] = result[component] + moment;
        }
        Ok(result)
    }

    // Shared original source grouping for detector power and full coherency.
    // Equal-frequency fields interfere here, before any temporal moments.
    fn group_fields(
        &self,
        bank: usize,
        at: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<JonesFields, CoherentDiffractionError> {
        let mut fields = [[[Dd::ZERO; 2]; 2]; MAX_COHERENT_DIFFRACTION_TONES];
        for ((source, tone), &group) in self
            .sources
            .iter()
            .zip(&self.optics.tones)
            .zip(&self.optics.tone_groups)
        {
            poll(cancelled)?;
            let waves = if bank == 0 {
                source.mean
            } else {
                source.factors[bank - 1].waves
            };
            for wave in waves {
                poll(cancelled)?;
                if let Some(pupils) = &self.pupils {
                    let value = pupils[group].field(&tone.pupil, wave, at, cancelled)?;
                    for (field, added) in fields[group].iter_mut().zip(value) {
                        field[0] = field[0] + added[0];
                        field[1] = field[1] + added[1];
                    }
                    continue;
                }
                let h = tone.pupil.amplitude([
                    at[0] - wave.tangent_direction[0],
                    at[1] - wave.tangent_direction[1],
                ])?;
                let hr = Dd::from_f64(h.re);
                let hi = Dd::from_f64(h.im);
                for (field, a) in fields[group].iter_mut().zip(wave.amplitude_sqrt_w) {
                    let ar = Dd::from_f64(a.re);
                    let ai = Dd::from_f64(a.im);
                    let term = [ar * hr - ai * hi, ar * hi + ai * hr];
                    if term.iter().all(|v| v.hi == 0.0 && v.lo == 0.0)
                        && (a.re != 0.0 || a.im != 0.0)
                        && (h.re != 0.0 || h.im != 0.0)
                    {
                        return Err(CoherentDiffractionError::Unrepresentable);
                    }
                    field[0] = field[0] + term[0];
                    field[1] = field[1] + term[1];
                }
            }
        }
        Ok(fields)
    }
}
