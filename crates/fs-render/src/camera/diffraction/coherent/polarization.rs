//! Jones fields through the original scalar paraxial pupil transfer.
//!
//! The Jones axes are camera right/up, with common incoming propagation along
//! -camera.forward: (right, up, -forward) is the native right-handed frame.
//! Both components use this common transverse frame and the same scalar h;
//! no direction-dependent vector diffraction is inferred within the paraxial bound.
//! Fields interfere within each declared source bank;
//! independent covariance banks add positive semidefinite coherency matrices.
//! Equal wavelength alone never creates mutual coherence. This is not vector
//! high-NA diffraction, a polarization-dependent pupil, or coherent ray tracing.
//! Inputs are the existing Jones spectral-exposure wave/illumination/factor types,
//! with the same shared source coefficients and pupil-centre phase reference.
//!
//! Cxy = Ex conj(Ey), I = Cxx+Cyy, Q = Cxx-Cyy, U = 2 Re(Cxy),
//! V = 2 Im(Cxy), matching the native Stokes/diagonal-Jones convention.
//! Finite pixels integrate these moments after interference, not Jones fields.
//! This can reduce polarization degree when polarization varies across a pixel.
//! See TU Delft, *Polarisation* and *Interference and Coherence*:
//! <https://interactivetextbooks.tudelft.nl/tn2421/content/Chap4_Polarisation/Polarization_2022_01Clean.html>
//! <https://interactivetextbooks.tudelft.nl/tn2421/content/Chap5_Interference/InterferenceCoherence_2022_01Clean.html>.

use super::{
    CoherentDiffractionError, CoherentDiffractionLimits, CoherentDiffractionPixel,
    CoherentDiffractionPixelLimits, CoherentJonesPupilIllumination, CoherentJonesPupilWave,
    PhysicalCamera, PreparedDiffractionAmplitude,
    admit_coordinate, admit_count, pixel, poll, rounded, square,
};
use crate::polarization::Stokes;
use fs_math::{c64::C64, dd::Dd};

type JonesField = [[Dd; 2]; 2];
type Moments = [Dd; 4];

/// A physical coherency result in the declared common transverse pupil frame.
///
/// Units follow the containing result: tangent-coordinate power density for
/// point observations, or watts for finite-pixel integrals. Matrix entries and
/// Stokes are rounded separately from the same DD moments. In particular, a
/// weak diagonal remains available even when rounded I and Q become equal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedDiffractionCoherency {
    axis_powers: [f64; 2],
    cross: C64,
    stokes: Stokes,
}

impl PolarizedDiffractionCoherency {
    /// Independent nonnegative powers on the first and second pupil axes.
    #[must_use]
    pub const fn axis_powers(self) -> [f64; 2] {
        self.axis_powers
    }

    /// Complex cross moment Ex conj(Ey); the other matrix entry is its conjugate.
    #[must_use]
    pub const fn cross(self) -> C64 {
        self.cross
    }

    /// Native Stokes result, formed before rounding the separate matrix entries.
    #[must_use]
    pub const fn stokes(self) -> Stokes {
        self.stokes
    }
}

/// Monochromatic point observation per d(tan(theta_x)) d(tan(theta_y)).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedCoherentDiffractionSample {
    /// Deterministic mean Jones field, in sqrt(W) per tangent-coordinate length.
    pub mean_field: [C64; 2],
    /// Coherency of the propagated deterministic mean field.
    pub mean_coherency_density_w: PolarizedDiffractionCoherency,
    /// Sum of the independent propagated covariance coherencies.
    pub fluctuation_coherency_density_w: PolarizedDiffractionCoherency,
    /// Total physical coherency and Stokes density, without power renormalization.
    pub coherency_density_w: PolarizedDiffractionCoherency,
}

/// Integrated monochromatic Stokes/coherency power over a physical pixel.
///
/// No mean Jones field is assigned to an intensity-integrating detector.
/// The original pixel tangent-coordinate mapping supplies the physical measure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PolarizedCoherentDiffractionPixelPower {
    /// Collected coherency of the propagated deterministic mean, in watts.
    pub mean_coherency_w: PolarizedDiffractionCoherency,
    /// Collected coherency of the independent source fluctuations, in watts.
    pub fluctuation_coherency_w: PolarizedDiffractionCoherency,
    /// Total collected coherency and native Stokes power, in watts.
    pub coherency_w: PolarizedDiffractionCoherency,
    /// Subcells in original pupil-axis order; zero for an empty illumination.
    pub subdivisions: [usize; 2],
    /// Number of positive spatial quadrature samples evaluated.
    pub quadrature_samples: usize,
    /// One original scalar transfer evaluation per authored wave and sample,
    /// shared by both Jones components.
    pub pupil_evaluations: usize,
}

/// Reusable Jones mean/covariance propagation through the original pupil.
///
/// Observations allocate nothing and cost O(authored waves + factors), with fixed
/// scratch independent of source count. Each wave evaluates the scalar h once.
/// DD field sums retain cancellation before powers/cross moments are formed;
/// this is not a relative-error guarantee at a dark fringe. The same source,
/// factor, paraxial, physical-pupil and pixel-work limits as the scalar query
/// remain authoritative. Source data are borrowed immutably.
pub struct PreparedPolarizedCoherentDiffraction<'camera, 'source> {
    pupil: PreparedDiffractionAmplitude<'camera>,
    illumination: CoherentJonesPupilIllumination<'source>,
    wavelength_nm: f64,
    waves: usize,
}

impl PhysicalCamera {
    /// Prepare explicit Jones illumination in the camera's common pupil frame.
    ///
    /// Focused and declared-plane defocused pupils share the existing transfer,
    /// including its original physical wavelength and signed quadratic phase.
    /// No coherence, polarization basis or object-to-pupil phase is inferred.
    pub fn prepare_polarized_coherent_diffraction<'camera, 'source>(
        &'camera self,
        wavelength_nm: f64,
        illumination: CoherentJonesPupilIllumination<'source>,
        limits: CoherentDiffractionLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPolarizedCoherentDiffraction<'camera, 'source>, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        admit_count("source factors", illumination.factors.len(), limits.maximum_factors)?;
        let mut waves = illumination.mean.len();
        admit_count("pupil waves", waves, limits.maximum_waves)?;
        for factor in illumination.factors {
            poll(&mut cancelled)?;
            waves = waves.checked_add(factor.waves.len()).ok_or(CoherentDiffractionError::Limit {
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
                if wave.amplitude_sqrt_w.iter().any(|a| !a.re.is_finite() || !a.im.is_finite()) {
                    return Err(CoherentDiffractionError::NonFiniteAmplitude);
                }
            }
        }
        let pupil = PreparedDiffractionAmplitude::new(self, wavelength_nm)?;
        poll(&mut cancelled)?;
        Ok(PreparedPolarizedCoherentDiffraction {
            pupil, illumination, wavelength_nm, waves,
        })
    }
}

impl PreparedPolarizedCoherentDiffraction<'_, '_> {
    #[must_use]
    pub const fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }

    /// Total original wave terms, including all zero Jones amplitudes.
    #[must_use]
    pub const fn wave_count(&self) -> usize {
        self.waves
    }

    #[must_use]
    pub fn factor_count(&self) -> usize {
        self.illumination.factors.len()
    }

    /// Observe full point Stokes/coherency without state changes or allocation.
    pub fn observe(
        &self,
        tangent_coordinate: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentDiffractionSample, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        admit_coordinate(tangent_coordinate)?;
        let (mean, coherent, fluctuations) = self.field_moments(tangent_coordinate, &mut cancelled)?;
        let result = PolarizedCoherentDiffractionSample {
            mean_field: [
                C64::new(rounded(mean[0][0])?, rounded(mean[0][1])?),
                C64::new(rounded(mean[1][0])?, rounded(mean[1][1])?),
            ],
            mean_coherency_density_w: coherency(coherent)?,
            fluctuation_coherency_density_w: coherency(fluctuations)?,
            coherency_density_w: coherency(std::array::from_fn(|i| coherent[i] + fluctuations[i]))?,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Integrate the positive coherency measure over the original finite pixel.
    ///
    /// Each quadrature sample sums fields inside each coherent bank, then adds
    /// independent coherencies. The original positive spatial rule integrates
    /// all four matrix moments with one shared work count. Refine the existing
    /// phase-span limit to assess quadrature error. Cancellation returns no
    /// partial result; the same preparation remains reusable.
    pub fn integrate_pixel(
        &self,
        pixel: CoherentDiffractionPixel,
        limits: CoherentDiffractionPixelLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PolarizedCoherentDiffractionPixelPower, CoherentDiffractionError> {
        poll(&mut cancelled)?;
        let maximum_step = self.pupil.scale * (0.5 * limits.maximum_phase_span_rad);
        let (subdivisions, quadrature_samples, pupil_evaluations) =
            pixel::pixel_work(pixel, limits, [maximum_step; 2], self.waves)?;
        let mut nonzero = [false; 4];
        let power: [Dd; 8] = pixel::integrate_moments(
            pixel, subdivisions, &mut cancelled, |at, cancelled| {
                let (_, mean, fluctuations) = self.field_moments(at, cancelled)?;
                for (flag, value) in nonzero.iter_mut().zip([mean[0], mean[1], fluctuations[0], fluctuations[1]]) {
                    *flag |= !is_zero(value);
                }
                Ok(std::array::from_fn(|i| if i < 4 { mean[i] } else { fluctuations[i - 4] }))
            },
        )?;
        for (active, index) in nonzero.into_iter().zip([0, 1, 4, 5]) {
            if active && is_zero(power[index]) {
                return Err(CoherentDiffractionError::Unrepresentable);
            }
        }
        let mean = std::array::from_fn(|i| power[i]);
        let fluctuations = std::array::from_fn(|i| power[i + 4]);
        let result = PolarizedCoherentDiffractionPixelPower {
            mean_coherency_w: coherency(mean)?,
            fluctuation_coherency_w: coherency(fluctuations)?,
            coherency_w: coherency(std::array::from_fn(|i| mean[i] + fluctuations[i]))?,
            subdivisions, quadrature_samples, pupil_evaluations,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn field_moments(
        &self,
        at: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(JonesField, Moments, Moments), CoherentDiffractionError> {
        let mean = self.field(self.illumination.mean, at, cancelled)?;
        let coherent = outer(mean)?;
        let mut fluctuations = [Dd::ZERO; 4];
        for factor in self.illumination.factors {
            poll(cancelled)?;
            let moment = outer(self.field(factor.waves, at, cancelled)?)?;
            for i in 0..4 {
                fluctuations[i] = fluctuations[i] + moment[i];
            }
        }
        for value in fluctuations {
            rounded(value)?;
        }
        Ok((mean, coherent, fluctuations))
    }

    fn field(
        &self,
        waves: &[CoherentJonesPupilWave],
        at: [f64; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<JonesField, CoherentDiffractionError> {
        let mut sum = [[Dd::ZERO; 2]; 2];
        for wave in waves {
            poll(cancelled)?;
            let h = self.pupil.amplitude([
                at[0] - wave.tangent_direction[0],
                at[1] - wave.tangent_direction[1],
            ])?;
            let hr = Dd::from_f64(h.re);
            let hi = Dd::from_f64(h.im);
            for (component, a) in sum.iter_mut().zip(wave.amplitude_sqrt_w) {
                let ar = Dd::from_f64(a.re);
                let ai = Dd::from_f64(a.im);
                let term = [ar * hr - ai * hi, ar * hi + ai * hr];
                // A product of two nonzero complex values cannot be zero.
                // Do not turn an underflowed source term into exact extinction.
                if term.into_iter().all(is_zero)
                    && (a.re != 0.0 || a.im != 0.0)
                    && (h.re != 0.0 || h.im != 0.0)
                {
                    return Err(CoherentDiffractionError::Unrepresentable);
                }
                component[0] = component[0] + term[0];
                component[1] = component[1] + term[1];
            }
        }
        for component in sum {
            for value in component {
                rounded(value)?;
            }
        }
        Ok(sum)
    }
}

fn is_zero(value: Dd) -> bool {
    value.hi == 0.0 && value.lo == 0.0
}

fn outer(field: JonesField) -> Result<Moments, CoherentDiffractionError> {
    let x = square(field[0])?;
    let y = square(field[1])?;
    // A represented nonzero field cannot physically have zero axis power.
    // Refuse arithmetic underflow rather than report false exact extinction.
    for (power, component) in [x, y].into_iter().zip(field) {
        if is_zero(power) && component.into_iter().any(|value| !is_zero(value)) {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
    }
    let cross_re = field[0][0] * field[1][0] + field[0][1] * field[1][1];
    let cross_im = field[0][1] * field[1][0] - field[0][0] * field[1][1];
    rounded(cross_re)?;
    rounded(cross_im)?;
    Ok([x, y, cross_re, cross_im])
}

fn coherency(moment: Moments) -> Result<PolarizedDiffractionCoherency, CoherentDiffractionError> {
    let axis_powers = [rounded(moment[0])?, rounded(moment[1])?];
    if axis_powers.iter().any(|x| *x < 0.0) {
        return Err(CoherentDiffractionError::Unrepresentable);
    }
    let cross = C64::new(rounded(moment[2])?, rounded(moment[3])?);
    let two = Dd::from_f64(2.0);
    let stokes = Stokes::try_new(
        rounded(moment[0] + moment[1])?,
        rounded(moment[0] - moment[1])?,
        rounded(two * moment[2])?,
        rounded(two * moment[3])?,
    ).ok_or(CoherentDiffractionError::Unrepresentable)?;
    Ok(PolarizedDiffractionCoherency { axis_powers, cross, stokes })
}
