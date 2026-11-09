//! Stationary, separable temporal coherence on the original narrowband receiver.
//!
//! For a normalized power spectrum with FWHM D in Hz, Fourier transformation
//! of a Lorentzian gives g(tau)=exp(-pi D |tau|); a Gaussian
//! exp(-4 ln(2) (nu/D)^2) gives g(tau)=exp(-(pi D tau)^2/(4 ln(2))).
//! These are field correlations, not power transmissions: arm energy stays
//! unchanged while the cross-arm coherency is multiplied by g.
//!
//! References: NISTIR 7117, Eq. 2.4.10 (round-trip delay 2R/c), and TU Delft,
//! Interference and Coherence, on spectral correlation and fringe visibility.
//! <https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir7117.pdf>
//! <https://interactivetextbooks.tudelft.nl/tn2421/content/Chap5_Interference/InterferenceCoherence_2022_01Clean.html>

use super::{CoherentDiffractionError, Dd, rounded};
use fs_math::det;

#[derive(Clone, Copy, Debug, PartialEq)]
enum LineShape {
    Lorentzian,
    Gaussian,
}

/// Explicit stationary linewidth and complete relative envelope delay.
///
/// The source spectrum is normalized to unit power, centered on the original
/// query's vacuum wavelength. Each original spatial source bank has the same
/// line profile; independent banks remain independent. The pupil transfer and
/// uniform Jones optics are frozen at that center wavelength. Only the declared
/// relative envelope delay varies across the line: higher-order dispersion,
/// chromatic pupil changes, finite-shutter noise and source turn-on are outside
/// this narrowband ensemble-mean model. It does not alter the sharp-tone query.
///
/// The FWHM must be at most 1% of the optical center frequency when observed.
/// This is an admission boundary, not an accuracy certificate. The ideal line
/// shapes extend over all real frequency DETUNINGS; they do not supply a broad
/// positive-frequency material spectrum. Use the explicit spectral exposure
/// owner when frequency-dependent optics or a finite realization is required.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InterferometerTemporalCoherence {
    shape: LineShape,
    fwhm_hz: f64,
    delay_s: f64,
    relay_delay_gradient_s: [f64; 2],
}

impl InterferometerTemporalCoherence {
    /// Lorentzian POWER-spectrum full width at half maximum in Hz.
    ///
    /// `relative_group_delay_s` is the complete arm-1 minus arm-0 envelope
    /// delay at the receiver centers, including all propagation and optics.
    /// It is deliberately separate from the receiver's optical PHASE length:
    /// dividing phase OPD by c is valid only for nondispersive paths. For a
    /// vacuum path mismatch L, supply L/c; a Michelson mirror motion d gives
    /// 2d/c. Source power is unchanged. Zero width recovers the original query.
    pub fn try_lorentzian(
        fwhm_hz: f64,
        relative_group_delay_s: f64,
    ) -> Result<Self, CoherentDiffractionError> {
        Self::new(LineShape::Lorentzian, fwhm_hz, relative_group_delay_s)
    }

    /// Gaussian POWER-spectrum FWHM in Hz, with the same explicit delay law
    /// and narrowband assumptions as [`Self::try_lorentzian`].
    pub fn try_gaussian(
        fwhm_hz: f64,
        relative_group_delay_s: f64,
    ) -> Result<Self, CoherentDiffractionError> {
        Self::new(LineShape::Gaussian, fwhm_hz, relative_group_delay_s)
    }

    fn new(
        shape: LineShape,
        fwhm_hz: f64,
        delay_s: f64,
    ) -> Result<Self, CoherentDiffractionError> {
        if !fwhm_hz.is_finite() || fwhm_hz < 0.0 || !delay_s.is_finite() {
            return Err(CoherentDiffractionError::InvalidTemporalCoherence);
        }
        Ok(Self {
            shape,
            fwhm_hz,
            delay_s,
            relay_delay_gradient_s: [0.0; 2],
        })
    }

    /// Author the complete linear envelope-delay variation over a relay pixel.
    ///
    /// Delay is `delay_at_centers + dot(gradient, relay_tangent_offset)` in
    /// seconds. Include any sensor-chirp group delay here; it is not inferred
    /// from phase index. A zero gradient explicitly selects a spatially uniform
    /// delay approximation. The original pixel sampler resolves this envelope
    /// as well as the optical fringe. No higher spatial delay terms are inferred.
    pub fn try_with_relay_delay_gradient(
        mut self,
        seconds_per_tangent: [f64; 2],
    ) -> Result<Self, CoherentDiffractionError> {
        if seconds_per_tangent.iter().any(|value| !value.is_finite()) {
            return Err(CoherentDiffractionError::InvalidTemporalCoherence);
        }
        self.relay_delay_gradient_s = seconds_per_tangent;
        Ok(self)
    }

    pub(in super::super) fn admit(
        self,
        wavelength_nm: f64,
    ) -> Result<(), CoherentDiffractionError> {
        if self.fwhm_hz > 0.01 * (299_792_458.0e9 / wavelength_nm) {
            Err(CoherentDiffractionError::InvalidTemporalCoherence)
        } else {
            Ok(())
        }
    }

    pub(in super::super) fn uniform_delay(self) -> bool {
        self.fwhm_hz == 0.0 || self.relay_delay_gradient_s == [0.0; 2]
    }

    pub(in super::super) fn weights(
        self,
        offset: [f64; 2],
    ) -> Result<TemporalWeights, CoherentDiffractionError> {
        if self.fwhm_hz == 0.0 {
            return Ok(TemporalWeights::COHERENT);
        }
        let delay = Dd::from_f64(self.delay_s)
            + Dd::from_f64(self.relay_delay_gradient_s[0]) * Dd::from_f64(offset[0])
            + Dd::from_f64(self.relay_delay_gradient_s[1]) * Dd::from_f64(offset[1]);
        let q = self.fwhm_hz * rounded(delay.abs())?;
        // Finite inputs may overflow this product; both analytic correlations
        // then round to zero. Saturation is the incoherent limit, not a failure.
        if q > 1000.0 {
            return Ok(TemporalWeights::INCOHERENT);
        }
        let x = core::f64::consts::PI * q;
        let exponent = match self.shape {
            LineShape::Lorentzian => x,
            LineShape::Gaussian => x * (x / (4.0 * core::f64::consts::LN_2)),
        };
        if exponent == 0.0 && delay != Dd::ZERO {
            return Err(CoherentDiffractionError::Unrepresentable);
        }
        // Retain whichever component is weak. Subtracting an f64 exp from one
        // would erase linewidth-induced dark-port leakage for tiny D*tau.
        let (coherent, incoherent) = if exponent < core::f64::consts::LN_2 {
            let loss = Dd::from_f64(-det::expm1(-exponent));
            (Dd::ONE - loss, loss)
        } else {
            let retained = Dd::from_f64(det::exp(-exponent));
            (retained, Dd::ONE - retained)
        };
        Ok(TemporalWeights {
            coherent,
            incoherent,
            field_scale: coherent.sqrt(),
        })
    }

    pub(in super::super) fn envelope_rate(
        self,
        axis: usize,
    ) -> Result<f64, CoherentDiffractionError> {
        // |dg/dtau| <= pi*D for BOTH line shapes. This bounds absolute
        // coherence variation through the Lorentzian cusp as well as a Gaussian.
        rounded(
            Dd::from_f64(self.fwhm_hz)
                * Dd::from_f64(self.relay_delay_gradient_s[axis]).abs()
                * Dd::from_f64(core::f64::consts::PI),
        )
    }
}

#[derive(Clone, Copy)]
pub(in super::super) struct TemporalWeights {
    pub coherent: Dd,
    pub incoherent: Dd,
    pub field_scale: Dd,
}

impl TemporalWeights {
    pub const COHERENT: Self = Self {
        coherent: Dd::ONE,
        incoherent: Dd::ZERO,
        field_scale: Dd::ONE,
    };
    const INCOHERENT: Self = Self {
        coherent: Dd::ZERO,
        incoherent: Dd::ONE,
        field_scale: Dd::ZERO,
    };
}
