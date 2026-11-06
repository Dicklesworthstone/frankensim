//! Finite glass thickness with incoherent internal reflections.
//!
//! Each polarization sums its own geometric series before unpolarized powers
//! are averaged. The two smooth, uncoated faces have the same exterior medium.
//! This is the native real-index/Beer law, not coherent thin-film interference,
//! a rough-surface model, or an inference of absorption from scattering loss.
//! See SCHOTT TIE-35, section 1, for internal versus external transmittance;
//! here the internal-reflection sum retains absorption at every traversal.

use super::{
    AbsorptionModel, BeerLambertAbsorption, DielectricError, DielectricGlass, DielectricIor,
    SpectralAbsorption, dielectric_power_coefficients, fresnel_dielectric,
};
use fs_math::det;

/// Fractions of incident normal power, including every internal reflection.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabPowerBalance {
    /// Power returning to the incident exterior medium.
    pub reflectance: f64,
    /// Power leaving the rear face into the identical exterior medium.
    pub transmittance: f64,
    /// Power absorbed inside the glass, excluding interface reflection.
    /// Evaluated directly to retain optically thin nonzero heating.
    pub absorptance: f64,
}

/// Incoherent power response of one homogeneous parallel-sided glass plate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncoherentSlabResponse {
    /// Electric field perpendicular to the incidence plane (s polarization).
    pub perpendicular: SlabPowerBalance,
    /// Electric field parallel to the incidence plane (p polarization).
    pub parallel: SlabPowerBalance,
    /// Cosine of the refracted angle, or zero at total internal reflection.
    pub transmitted_cosine: f64,
    /// Absorption optical depth of a single refracted traversal.
    /// Infinite when no propagating internal traversal exists.
    pub optical_depth: f64,
}

impl IncoherentSlabResponse {
    /// Unpolarized power; averaging interface Fresnel values before summing
    /// internal bounces would give a different, incorrect oblique response.
    #[must_use]
    pub fn unpolarized(self) -> SlabPowerBalance {
        SlabPowerBalance {
            reflectance: 0.5 * self.perpendicular.reflectance + 0.5 * self.parallel.reflectance,
            transmittance: 0.5 * self.perpendicular.transmittance
                + 0.5 * self.parallel.transmittance,
            absorptance: 0.5 * self.perpendicular.absorptance + 0.5 * self.parallel.absorptance,
        }
    }
}

impl DielectricGlass {
    /// Evaluate wavelength-, angle- and thickness-dependent slab opacity.
    ///
    /// `thickness_m` is positive normal thickness, not distance along the ray.
    /// `incident_cosine` is in `(0,1]`; `exterior_index` is the same lossless
    /// ordinary real index on both sides. Original phase/absorption wavelength
    /// support is enforced even when the incident wave is totally reflected.
    /// No coherence, coating, scattering, lateral exit position, or traversal
    /// delay is modeled. The bulk uses the existing weak-absorption dielectric
    /// approximation; use complex-index optics for strongly absorbing faces.
    pub fn incoherent_slab_response(
        &self,
        wavelength_nm: f64,
        thickness_m: f64,
        incident_cosine: f64,
        exterior_index: f64,
    ) -> Result<IncoherentSlabResponse, DielectricError> {
        validate_thickness(thickness_m)?;
        if incident_cosine <= 0.0 {
            return Err(DielectricError::InvalidInterface);
        }
        let index = self.ior().eval(wavelength_nm)?;
        let alpha = self.absorption_ref().extinction_per_m(wavelength_nm)?;
        let interface = fresnel_dielectric(incident_cosine, exterior_index, index)?;
        if interface.total_internal_reflection {
            let reflected = SlabPowerBalance {
                reflectance: 1.0,
                transmittance: 0.0,
                absorptance: 0.0,
            };
            return Ok(IncoherentSlabResponse {
                perpendicular: reflected,
                parallel: reflected,
                transmitted_cosine: 0.0,
                optical_depth: f64::INFINITY,
            });
        }
        let cosine = interface.transmitted_cosine;
        let mut depth = alpha * thickness_m / cosine;
        if alpha > 0.0 && (depth == 0.0 || !depth.is_finite()) {
            // Recover a representable product when the intermediate multiply
            // underflows/overflows. A genuinely opaque depth may be infinite.
            depth = det::exp(det::ln(alpha) + det::ln(thickness_m) - det::ln(cosine));
        }
        let (reflection, transmission) = if exterior_index == index {
            // Preserve the original Fresnel null-interface fast path. At an
            // arbitrarily grazing matched ray, squared normal wave numbers
            // may underflow even though exact zero reflection is known.
            ([0.0; 2], [1.0; 2])
        } else {
            dielectric_power_coefficients(incident_cosine, exterior_index, index, cosine)
        };
        Ok(IncoherentSlabResponse {
            perpendicular: power_balance(reflection[0] * reflection[0], transmission[0], depth),
            parallel: power_balance(reflection[1] * reflection[1], transmission[1], depth),
            transmitted_cosine: cosine,
            optical_depth: depth,
        })
    }
}

fn power_balance(r: f64, t: f64, depth: f64) -> SlabPowerBalance {
    if depth == 0.0 {
        return SlabPowerBalance {
            reflectance: 2.0 * r / (1.0 + r),
            transmittance: t / (1.0 + r),
            absorptance: 0.0,
        };
    }
    let a = det::exp(-depth);
    let loss = -det::expm1(-depth);
    // D = 1 - r² a², written without subtracting nearly equal unit powers.
    // Native Poynting transmission t retains the complement when r rounds to 1.
    let denominator = t * (1.0 + r) + r * r * -det::expm1(-2.0 * depth);
    let transmittance = t * (t / denominator) * a;
    SlabPowerBalance {
        reflectance: r + r * a * transmittance,
        transmittance,
        // Cancel (1+r*a) analytically from t*(1-a)*(1+r*a)/D.
        absorptance: t * (loss / (t + r * loss)),
    }
}

impl BeerLambertAbsorption {
    /// Convert original spectral INTERNAL transmittance at a positive path
    /// length into the existing absorption law, without RGB lifting.
    ///
    /// Knots retain their vacuum-nanometre support. Values must be in `(0,1]`;
    /// zero is a measurement limit, not a finite absorption coefficient.
    /// Between knots alpha interpolates linearly, so a different thickness
    /// obeys Beer-Lambert scaling. Reflection/scattering corrections must
    /// already have been applied by the measurement's stated physical model.
    pub fn try_from_spectral_internal_transmittance_in_band(
        knots: &[(f64, f64)],
        distance_m: f64,
    ) -> Result<Self, DielectricError> {
        validate_thickness(distance_m)?;
        let mut spectrum = SpectralAbsorption::try_new_in_band(knots)?;
        for (_, transmission) in &mut spectrum.knots[..spectrum.count] {
            validate_transmission(*transmission)?;
            *transmission = absorption_coefficient(-det::ln(*transmission), distance_m)?;
        }
        Ok(Self {
            model: AbsorptionModel::Spectral(spectrum),
        })
    }

    /// Recover bulk absorption from normal-incidence EXTERNAL transmission of
    /// a smooth, uncoated, incoherent parallel plate in one exterior medium.
    ///
    /// The detector must collect all specularly transmitted internal-bounce
    /// orders. Reflection is removed using the supplied original phase-index
    /// law at each wavelength; this is not the approximation `-ln(T)/d`.
    /// Above-clear transmission, zero transmission, unsupported wavelengths,
    /// and unrepresentable finite absorption refuse rather than being clipped.
    /// The resulting law feeds existing glass and thermal transport directly.
    /// This numerical conversion adds no measurement/uncertainty authority.
    pub fn try_from_normal_slab_transmittance_in_band(
        knots: &[(f64, f64)],
        thickness_m: f64,
        phase_index: &DielectricIor,
        exterior_index: f64,
    ) -> Result<Self, DielectricError> {
        validate_thickness(thickness_m)?;
        let mut spectrum = SpectralAbsorption::try_new_in_band(knots)?;
        for (wavelength, transmission) in &mut spectrum.knots[..spectrum.count] {
            validate_transmission(*transmission)?;
            let index = phase_index.eval(*wavelength)?;
            let interface = fresnel_dielectric(1.0, exterior_index, index)?;
            let r = interface.reflectance;
            let t = interface.transmittance;
            let clear = t / (1.0 + r);
            if *transmission > clear {
                return Err(DielectricError::InvalidAbsorption);
            }
            let y = *transmission / clear;
            let deficit = (clear - *transmission) / clear;
            let b = t * (1.0 + r);
            let root = det::hypot(b, 2.0 * r * y);
            // Solve y = (1-r²)*a/(1-r²*a²), a = exp(-alpha*d).
            // Rationalizing 1-a preserves weak loss near the clear ceiling.
            let loss = (2.0 * b * deficit) / (b + 2.0 * r * r * y + root);
            let depth = if loss <= 0.25 {
                negative_log_complement(loss)
            } else {
                -det::ln((2.0 * y) / (b + root))
            };
            *transmission = absorption_coefficient(depth, thickness_m)?;
        }
        Ok(Self {
            model: AbsorptionModel::Spectral(spectrum),
        })
    }
}

fn validate_thickness(value: f64) -> Result<(), DielectricError> {
    if !value.is_finite() || value <= 0.0 {
        Err(DielectricError::InvalidAbsorption)
    } else {
        Ok(())
    }
}

fn validate_transmission(value: f64) -> Result<(), DielectricError> {
    if !value.is_finite() || value <= 0.0 || value > 1.0 {
        Err(DielectricError::InvalidAbsorption)
    } else {
        Ok(())
    }
}

fn absorption_coefficient(depth: f64, distance_m: f64) -> Result<f64, DielectricError> {
    let coefficient = depth / distance_m;
    if !coefficient.is_finite() || coefficient < 0.0 || (depth > 0.0 && coefficient == 0.0) {
        Err(DielectricError::InvalidAbsorption)
    } else {
        Ok(if coefficient == 0.0 { 0.0 } else { coefficient })
    }
}

fn negative_log_complement(loss: f64) -> f64 {
    // -ln(1-loss) = 2*atanh(loss/(2-loss)); loss <= 1/4 means x <= 1/7.
    // Twelve fixed terms leave truncation below f64 rounding, including a
    // representable positive loss whose complement rounds exactly to one.
    let x = loss / (2.0 - loss);
    let squared = x * x;
    let mut term = x;
    let mut sum = x;
    for k in 1..12 {
        term *= squared;
        sum += term / f64::from(2 * k + 1);
    }
    2.0 * sum
}

#[cfg(test)]
mod tests;
