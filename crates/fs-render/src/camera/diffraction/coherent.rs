//! Coherent scalar illumination of the original physical camera pupil.
//!
//! Each authored wave is referenced to the SAME entrance-pupil plane and
//! polarization. Its pupil dependence is `a * exp(+i k theta_source.u)`;
//! the camera supplies the common signed quadratic phase. Fourier shifting
//! therefore gives `a * h(theta - theta_source)`. No separately shifted
//! detector chirp or unprovided object-to-pupil propagation phase is added.
//!
//! The deterministic mean and each covariance factor propagate as fields.
//! Only their squared magnitudes add: `|H mean|² + sum_r |H factor_r|²`.
//! This is a positive factorization of partial spatial coherence, not an
//! inference that equal wavelengths imply coherent sources. Weights belong
//! in factor amplitudes as square roots. The result is not renormalized to
//! the sum of isolated-wave powers, since their pupil fields can interfere.
//!
//! All waves share one vacuum wavelength and the declared source-plane
//! curvature. This clear-pupil, scalar, paraxial query does not make the
//! radiance path tracer coherent or infer phases from arbitrary scene points.
//! The companion Jones query retains both common-frame transverse components
//! and their coherency through this same scalar transfer. Finite rectangular pixels
//! integrate the same field with bandwidth-guided positive quadrature. The
//! separate spectral exposure query retains cross-frequency beats over an
//! actual shutter and the physical sensor plane's common Fresnel phase. Its
//! Jones counterpart preserves two paraxial polarization components and projects
//! analyzer fields before photodetection. See TU Delft,
//! *Interference and Coherence*, for field versus intensity superposition:
//! <https://interactivetextbooks.tudelft.nl/tn2421/content/Chap5_Interference/InterferenceCoherence_2022_01Clean.html>.

use super::{CameraError, PhysicalCamera, PreparedDiffractionAmplitude};
use fs_math::{c64::C64, dd::Dd};

mod exposure;
mod interferometer;
mod pixel;
mod polarization;
pub use interferometer::{CoherentDiffractionInterferometer, CoherentInterferometerSample};
pub use polarization::{
    PolarizedCoherentDiffractionPixelPower, PolarizedCoherentDiffractionSample,
    PolarizedDiffractionCoherency, PreparedPolarizedCoherentDiffraction,
};
pub use exposure::{
    CoherentDiffractionExposureLimits, CoherentDiffractionExposurePixelLimits,
    CoherentDiffractionExposureSample, CoherentDiffractionPixelExposure,
    CoherentDiffractionShutter, CoherentDiffractionTone, CoherentJonesDetector,
    CoherentJonesDiffractionTone, CoherentJonesPupilIllumination, CoherentJonesPupilLayer,
    CoherentJonesPupilLimits, CoherentJonesPupilResponse, CoherentJonesPupilWave,
    JonesPupilCoherenceFactor, MAX_COHERENT_DIFFRACTION_TONES, MAX_COHERENT_JONES_PUPIL_MODES,
    PolarizedCoherentDiffractionExposureSample, PolarizedCoherentDiffractionPixelExposure,
    PreparedCoherentDiffractionExposure, PreparedCoherentJonesDiffractionExposure,
};
pub use pixel::{
    CoherentDiffractionPixel, CoherentDiffractionPixelLimits, CoherentDiffractionPixelPower,
};

/// Largest admitted magnitude of a source or observation tangent coordinate.
/// This paraxial model boundary is not an optical error certificate.
pub const MAX_COHERENT_PUPIL_TANGENT: f64 = 0.1;

/// One tilted pupil field, with a common pupil-centre phase reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentPupilWave {
    /// Two dimensionless tangent coordinates along the camera pupil axes.
    pub tangent_direction: [f64; 2],
    /// Complex amplitude in sqrt(W). Its squared magnitude is the power this
    /// wave would transmit through this clear pupil IN ISOLATION. The common
    /// optical carrier is omitted; no peak/RMS factor is added to this power
    /// normalization. Relative source and propagation phases must be authored.
    pub amplitude_sqrt_w: C64,
}

/// One column of an explicitly authored source covariance factorization.
/// Waves inside this column interfere; distinct columns are uncorrelated.
#[derive(Clone, Copy, Debug)]
pub struct PupilCoherenceFactor<'a> {
    /// Amplitudes already include the square root of this factor's weight.
    pub waves: &'a [CoherentPupilWave],
}

/// Deterministic mean plus zero-mean, unit-covariance source factors.
/// Empty banks represent zero fields. Coherence is never created across banks.
#[derive(Clone, Copy, Debug)]
pub struct CoherentPupilIllumination<'a> {
    /// Coherent mean field on the common entrance-pupil reference plane.
    pub mean: &'a [CoherentPupilWave],
    /// Independent covariance columns; each column can contain many waves.
    pub factors: &'a [PupilCoherenceFactor<'a>],
}

/// Bounds on preparation and each observation. No source-pair matrix is built.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CoherentDiffractionLimits {
    /// Total authored waves in the mean and all factors, including zero waves.
    pub maximum_waves: usize,
    /// Number of covariance columns, including empty columns.
    pub maximum_factors: usize,
}

/// A point observation per `d(tan(theta_x)) d(tan(theta_y))`.
/// These densities are NOT irradiance per square metre of detector. A physical
/// image-plane conversion needs its coordinate Jacobian and pixel integration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDiffractionSample {
    /// Propagated deterministic mean in sqrt(W) per tangent-coordinate length.
    /// A common output-plane phase is omitted consistently for every source.
    pub mean_field: C64,
    /// Squared mean field, retained through DD field cancellation before rounding.
    pub mean_power_density_w: f64,
    /// Sum of squared propagated covariance columns, excluding the mean.
    pub fluctuation_power_density_w: f64,
    /// Ensemble-mean power density: mean plus fluctuation contributions.
    pub power_density_w: f64,
}

/// Invalid illumination, exceeded work bounds, or an unavailable pupil query.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoherentDiffractionError {
    /// Cancellation was observed; no sample was returned.
    Cancelled,
    /// A field coefficient is not finite.
    NonFiniteAmplitude,
    /// A source or observation coordinate is outside the declared paraxial domain.
    OutsideParaxialDomain,
    /// Pixel bounds are not finite and strictly increasing on both axes.
    InvalidPixel,
    /// The maximum subcell intensity phase span is outside `(0, 1]` radians.
    InvalidPixelResolution,
    /// A spectral line needs a finite, strictly positive vacuum frequency.
    InvalidFrequency,
    /// Sensor distance must be positive and agree with a declared lens conjugate.
    InvalidSensorPlane,
    /// Arm weights must be passive and optical path/combiner phase finite.
    InvalidInterferometer,
    /// Every tone must declare the same ordered covariance-factor identities.
    InconsistentCoherenceFactors,
    /// An ideal Jones analyzer needs a finite, nonzero polarization direction.
    InvalidPolarizationAnalyzer,
    /// Jones-pupil powers must lie in [0,1] and all layer parameters be finite.
    InvalidJonesPupil,
    /// Every distinct source frequency needs exactly one pupil response.
    InconsistentJonesPupilSpectrum,
    /// The original analytic shutter kernel refused its clock or contraction.
    Shutter(fs_fft::harmonic_window::WindowError),
    /// A preparation count exceeds the caller's explicit bound.
    Limit {
        /// Quantity whose admitted count was exceeded.
        quantity: &'static str,
        /// Requested count, or usize::MAX if addition overflowed.
        requested: usize,
        /// Caller-supplied maximum.
        maximum: usize,
    },
    /// A field contraction or its physical power exceeds floating-point range.
    Unrepresentable,
    /// The original camera/pupil/wavelength admission refused the query.
    Camera(CameraError),
}

impl core::fmt::Display for CoherentDiffractionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coherent pupil diffraction: {self:?}")
    }
}
impl std::error::Error for CoherentDiffractionError {}
impl From<CameraError> for CoherentDiffractionError {
    fn from(value: CameraError) -> Self {
        Self::Camera(value)
    }
}
impl From<fs_fft::harmonic_window::WindowError> for CoherentDiffractionError {
    fn from(value: fs_fft::harmonic_window::WindowError) -> Self {
        Self::Shutter(value)
    }
}

/// Reusable coherent image query with original pupil geometry and source banks.
///
/// Preparation caches the wavelength, material admission and pupil coefficients.
/// Every observation costs O(total authored waves), uses constant scratch and
/// allocates nothing. There is no dense source-pair expansion. DD accumulation
/// retains weak fields after strong destructive interference; it is not a
/// relative-error guarantee at a null. Source and camera data are borrowed
/// immutably, so changing either requires another preparation.
pub struct PreparedCoherentDiffraction<'camera, 'source> {
    camera: &'camera PhysicalCamera,
    pupil: PreparedDiffractionAmplitude<'camera>,
    illumination: CoherentPupilIllumination<'source>,
    wavelength_nm: f64,
    waves: usize,
}

impl PhysicalCamera {
    /// Prepare explicit coherent illumination for the selected focused or
    /// declared-plane defocused pupil. Layered scene-depth policies are refused:
    /// this query requires one shared pupil curvature and phase reference.
    pub fn prepare_coherent_diffraction<'camera, 'source>(
        &'camera self,
        wavelength_nm: f64,
        illumination: CoherentPupilIllumination<'source>,
        limits: CoherentDiffractionLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedCoherentDiffraction<'camera, 'source>, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        admit_count(
            "source factors",
            illumination.factors.len(),
            limits.maximum_factors,
        )?;
        let mut waves = illumination.mean.len();
        admit_count("pupil waves", waves, limits.maximum_waves)?;
        for factor in illumination.factors {
            poll(&mut cancelled)?;
            waves =
                waves
                    .checked_add(factor.waves.len())
                    .ok_or(CoherentDiffractionError::Limit {
                        quantity: "pupil waves",
                        requested: usize::MAX,
                        maximum: limits.maximum_waves,
                    })?;
            admit_count("pupil waves", waves, limits.maximum_waves)?;
        }
        for bank in core::iter::once(illumination.mean)
            .chain(illumination.factors.iter().map(|factor| factor.waves))
        {
            poll(&mut cancelled)?;
            for wave in bank {
                poll(&mut cancelled)?;
                admit_coordinate(wave.tangent_direction)?;
                if !wave.amplitude_sqrt_w.re.is_finite() || !wave.amplitude_sqrt_w.im.is_finite() {
                    return Err(CoherentDiffractionError::NonFiniteAmplitude);
                }
            }
        }
        let pupil = PreparedDiffractionAmplitude::new(self, wavelength_nm)?;
        poll(&mut cancelled)?;
        Ok(PreparedCoherentDiffraction {
            camera: self,
            pupil,
            illumination,
            wavelength_nm,
            waves,
        })
    }
}

impl PreparedCoherentDiffraction<'_, '_> {
    /// Original vacuum wavelength, in nanometres.
    #[must_use]
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }

    /// Number of original authored wave terms, including zero amplitudes.
    #[must_use]
    pub const fn wave_count(&self) -> usize {
        self.waves
    }

    /// Number of independent covariance columns, excluding the mean.
    #[must_use]
    pub fn factor_count(&self) -> usize {
        self.illumination.factors.len()
    }

    /// Observe one tangent-plane coordinate without allocation or state changes.
    /// The original source order defines the deterministic reduction order.
    /// A cancellation or numerical refusal publishes no partial sample and
    /// leaves the same preparation reusable.
    pub fn observe(
        &self,
        tangent_coordinate: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentDiffractionSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        admit_coordinate(tangent_coordinate)?;
        let (mean, mean_power, fluctuation_power) =
            self.field_moments(tangent_coordinate, &mut cancelled)?;
        let result = CoherentDiffractionSample {
            mean_field: C64::new(rounded(mean[0])?, rounded(mean[1])?),
            mean_power_density_w: rounded(mean_power)?,
            fluctuation_power_density_w: rounded(fluctuation_power)?,
            power_density_w: rounded(mean_power + fluctuation_power)?,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn field_moments(
        &self,
        tangent_coordinate: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([Dd; 2], Dd, Dd), CoherentDiffractionError> {
        let mean = self.field(self.illumination.mean, tangent_coordinate, cancelled)?;
        let mean_power = square(mean)?;
        let mut fluctuation_power = Dd::ZERO;
        for factor in self.illumination.factors {
            poll(cancelled)?;
            fluctuation_power = fluctuation_power
                + square(self.field(factor.waves, tangent_coordinate, cancelled)?)?;
        }
        rounded(fluctuation_power)?;
        Ok((mean, mean_power, fluctuation_power))
    }

    fn field(
        &self,
        waves: &[CoherentPupilWave],
        observation: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Dd; 2], CoherentDiffractionError> {
        let mut sum = [Dd::ZERO; 2];
        for wave in waves {
            poll(cancelled)?;
            let h = self.pupil.amplitude([
                observation[0] - wave.tangent_direction[0],
                observation[1] - wave.tangent_direction[1],
            ])?;
            let ar = Dd::from_f64(wave.amplitude_sqrt_w.re);
            let ai = Dd::from_f64(wave.amplitude_sqrt_w.im);
            let hr = Dd::from_f64(h.re);
            let hi = Dd::from_f64(h.im);
            sum[0] = sum[0] + (ar * hr - ai * hi);
            sum[1] = sum[1] + (ar * hi + ai * hr);
        }
        // Validate before products, including intermediate DD overflow.
        rounded(sum[0])?;
        rounded(sum[1])?;
        Ok(sum)
    }
}

fn admit_count(
    quantity: &'static str,
    requested: usize,
    maximum: usize,
) -> Result<(), CoherentDiffractionError> {
    if requested > maximum {
        Err(CoherentDiffractionError::Limit {
            quantity,
            requested,
            maximum,
        })
    } else {
        Ok(())
    }
}

fn admit_coordinate(value: [f64; 2]) -> Result<(), CoherentDiffractionError> {
    if value.iter().any(|x| !x.is_finite())
        || value[0] * value[0] + value[1] * value[1]
            > MAX_COHERENT_PUPIL_TANGENT * MAX_COHERENT_PUPIL_TANGENT
    {
        Err(CoherentDiffractionError::OutsideParaxialDomain)
    } else {
        Ok(())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), CoherentDiffractionError> {
    if cancelled() {
        Err(CoherentDiffractionError::Cancelled)
    } else {
        Ok(())
    }
}

fn rounded(value: Dd) -> Result<f64, CoherentDiffractionError> {
    let result = value.hi + value.lo;
    if value.hi.is_finite() && value.lo.is_finite() && result.is_finite() {
        Ok(result)
    } else {
        Err(CoherentDiffractionError::Unrepresentable)
    }
}

fn square(value: [Dd; 2]) -> Result<Dd, CoherentDiffractionError> {
    let result = value[0] * value[0] + value[1] * value[1];
    rounded(result)?;
    Ok(result)
}
