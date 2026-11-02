//! Spectral metal wave response from the original visible complex-index source.
//!
//! The source uses absolute vacuum-referenced n and extinction k. With this
//! owner's exp(+i omega t) convention N=n-i k and epsilon=N^2/(mu c0^2).
//! Permeability is explicitly sourced in the SAME complete material snapshot;
//! no nonmagnetic default, DC-conductivity fit or second loss term is inserted.
//! Linear interpolation of the original n/k knots is confined to their band.
//!
//! Wave and absorption conventions: S. J. Byrnes, "Multilayer optical
//! calculations", sections 2--4, <https://arxiv.org/abs/1603.02720> (opposite
//! time convention). The general dispersion relation retains sourced mu.

use super::{ReflectionError, finite, positive};
use fs_material::{
    field_storage::linear_magnetic_permeability,
    state_point::{ResolvedMaterialStatePoint, VisibleConductorStatePoint},
};
use fs_math::{c64::C64, det};

/// Exact SI vacuum light speed, used only to interpret vacuum wavelengths.
const VACUUM_SPEED_M_S: f64 = 299_792_458.0;
const FREQUENCY_WAVELENGTH_NM: f64 = core::f64::consts::TAU * VACUUM_SPEED_M_S * 1.0e9;

/// Original spectral complex-index law with its explicitly sourced permeability.
///
/// Adopt a complete material bundle with
/// [`VisibleConductorStatePoint::try_from_material`]; that bundle must contain
/// the optical knots and absolute scalar magnetic permeability. The original
/// spectrum supplies all electric loss. Ohmic loss must not be added again.
/// This homogeneous isotropic monochromatic law does not supply a causal
/// time-domain realization or extrapolate the source band.
#[derive(Clone, Copy, Debug)]
pub struct OpticalConductorResponse<'a> {
    source: &'a VisibleConductorStatePoint,
    permeability_h_m: f64,
    index_to_sqrt_permittivity: f64,
}

impl<'a> OpticalConductorResponse<'a> {
    /// Validate the original permeability and retain the complete source owner.
    pub fn try_new(source: &'a VisibleConductorStatePoint) -> Result<Self, ReflectionError> {
        let permeability_h_m = linear_magnetic_permeability(source.resolved())?;
        let index_to_sqrt_permittivity =
            positive((VACUUM_SPEED_M_S * det::sqrt(permeability_h_m)).recip())?;
        Ok(Self {
            source,
            permeability_h_m,
            index_to_sqrt_permittivity,
        })
    }

    pub const fn source(&self) -> &'a VisibleConductorStatePoint {
        self.source
    }

    pub const fn material(&self) -> &'a ResolvedMaterialStatePoint {
        self.source.resolved()
    }

    pub const fn permeability_h_m(&self) -> f64 {
        self.permeability_h_m
    }

    /// Original n-i k at one angular frequency [rad/s], with no extrapolation.
    /// At most eight source intervals are visited and no allocation is made.
    pub fn complex_index(&self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        let omega = positive(omega_rad_s)?;
        let samples = self.source.samples();
        let first = samples[0].wavelength_nm;
        let last = samples[samples.len() - 1].wavelength_nm;
        if omega < FREQUENCY_WAVELENGTH_NM / last || omega > FREQUENCY_WAVELENGTH_NM / first {
            return Err(ReflectionError::Invalid(
                "optical wave frequency is outside the original n/k source band",
            ));
        }
        // Frequency admission precedes this rounding guard: it cannot extend
        // the source band, even by one representable frequency step.
        let wavelength = (FREQUENCY_WAVELENGTH_NM / omega).clamp(first, last);
        for pair in samples.windows(2) {
            let [left, right] = [pair[0], pair[1]];
            if wavelength <= right.wavelength_nm {
                let t =
                    (wavelength - left.wavelength_nm) / (right.wavelength_nm - left.wavelength_nm);
                let n = (1.0 - t) * left.eta + t * right.eta;
                let k = (1.0 - t) * left.k + t * right.k;
                return Ok(C64::new(n, -k));
            }
        }
        Err(ReflectionError::Invalid(
            "optical source interval is missing",
        ))
    }

    /// Absolute complex epsilon [F/m], including the original source loss.
    /// Negative real epsilon is valid for a passive metal; passivity requires
    /// the imaginary part to be nonpositive under exp(+i omega t).
    pub fn permittivity_f_m(&self, omega_rad_s: f64) -> Result<C64, ReflectionError> {
        let index = self.complex_index(omega_rad_s)?;
        // Scale before squaring, and factor the difference of squares so a
        // nearly zero real part does not lose the difference between n and k.
        let n = index.re * self.index_to_sqrt_permittivity;
        let k = -index.im * self.index_to_sqrt_permittivity;
        let epsilon = C64::new((n - k) * (n + k), -2.0 * n * k);
        if !finite(epsilon)
            || epsilon == C64::ZERO
            || epsilon.im > 0.0
            || (index.im < 0.0 && epsilon.im == 0.0)
        {
            return Err(ReflectionError::Invalid(
                "optical permittivity or its original loss is unrepresentable",
            ));
        }
        Ok(epsilon)
    }
}
