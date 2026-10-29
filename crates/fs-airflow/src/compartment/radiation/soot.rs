//! Soot absorption from an explicit particle population and optical source.
//!
//! For dilute Rayleigh spheres with a wavelength-independent complex refractive
//! index `m = n - i k`, `k >= 0`, the absorption coefficient is
//! `kappa_lambda = C0 f_v / lambda`, where
//! `C0 = 36 pi n k / ((n²-k²+2)² + 4 n² k²)`.
//! The same small-particle expression appears in `POD` in RADCAL:
//! <https://github.com/firemodels/radcal/blob/abe2a8f1e84ef0b3afda73d0c08d5c2e2964e401/Source/rcal.f90>.
//! RADCAL additionally resolves wavelength-dependent refractive indices; this
//! reduced model deliberately requires a caller-sourced constant pair instead.
//!
//! Integrating against a Planck spectrum at the supplied spectral temperature
//! gives `kappa_P = C0 f_v T [4 zeta(5)/zeta(4)] / c2`, `c2 = h c / k_B`.
//! This follows from `integral_0^infinity x^p/(exp(x)-1) dx = Gamma(p+1) zeta(p+1)`
//! after setting `x = c2/(lambda T)`. The weighting and wavelength conventions
//! follow FDS 6.11.1's gas absorption appendix:
//! <https://github.com/firemodels/fds/blob/ff928dbf68a2623ced3a007be75391c6dc6c9f90/Manuals/FDS_Technical_Reference_Guide/Appendices.tex>.
//!
//! The spectral temperature describes the radiation being averaged. Use the
//! local particle temperature for its emission mean; use a source's blackbody
//! temperature for that source's incident mean. Those are generally different.
//! Converting one mean to `LayerAbsorption` is an **estimated gray closure**:
//! it does not exactly represent a nonequilibrium spectrum, spectral filtering
//! along a finite path, or a nongray mixed enclosure.
//!
//! The caller must justify small primary particles, dilution, negligible
//! scattering, constant optical properties over the relevant spectrum, and
//! thermal equilibration if using the gas layer's temperature. No particle-size
//! cutoff, soot yield, refractive index, density, or molecular gas absorption is
//! guessed. Soot concentration is condensed particulate mass per mixture volume;
//! these methods never add it to ideal-gas species, moles, pressure or energy.

use super::{
    AbsorptionCoefficient, LayerAbsorption, RadiationError, finite, nonnegative, positive,
};
use fs_qty::{Density, Length, Temperature};

// Exact defining SI values; no rounded c2 or fitted emission coefficient.
const SECOND_RADIATION_CONSTANT_M_K: f64 = 6.626_070_15e-34 * 299_792_458.0 / 1.380_649e-23;
const ZETA_5: f64 = 1.036_927_755_143_37;
const PLANCK_WAVENUMBER_FACTOR: f64 = 360.0 * ZETA_5
    / (core::f64::consts::PI
        * core::f64::consts::PI
        * core::f64::consts::PI
        * core::f64::consts::PI);

/// Explicit source domain. Admission checks these numerical bounds; it does
/// not certify the source's Rayleigh, dilution or spectral assumptions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SootOpticalDomain {
    /// Inclusive range for the Planck weighting temperature [K]. This is not
    /// automatically the gas or particle temperature.
    pub spectral_temperature: [Temperature; 2],
    /// Largest admitted particulate volume fraction, strictly between 0 and 1.
    /// The source must justify that its value is dilute enough for this model.
    pub maximum_volume_fraction: f64,
}

/// Prepared, source-bearing constant-index Rayleigh absorption model.
/// Construction allocates its source once; coefficient evaluations allocate
/// nothing and use constant storage and arithmetic per population.
#[derive(Clone, Debug, PartialEq)]
pub struct RayleighSoot {
    refractive_index: [f64; 2],
    particle_density: Density,
    domain: SootOpticalDomain,
    source: String,
    spectral_factor: f64,
}

/// An evaluated Planck mean retaining its exact population and spectral source.
/// This borrows the immutable optical model and allocates no source string.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanckSootAbsorption<'a> {
    /// Soot-only Planck absorption [1/m], excluding scattering and molecular gas.
    pub coefficient: AbsorptionCoefficient,
    /// Particulate mass per total mixture volume [kg/m³].
    pub mass_concentration: Density,
    /// Particulate volume per total mixture volume, `concentration / density`.
    pub volume_fraction: f64,
    /// Actual blackbody weighting temperature, not an inferred gas temperature.
    pub spectral_temperature: Temperature,
    model: &'a RayleighSoot,
}

impl RayleighSoot {
    /// Bind `n > 0`, absorptive `k >= 0`, material density and declared domain
    /// to an explicit optical source. Zero `k` is a nonabsorbing population.
    /// A material density is needed even when callers later use volume fraction,
    /// so every report has an unambiguous particulate mass concentration.
    pub fn new(
        real_index: f64,
        absorptive_index: f64,
        particle_density: Density,
        domain: SootOpticalDomain,
        source: impl Into<String>,
    ) -> Result<Self, RadiationError> {
        let source = source.into();
        let [minimum, maximum] = domain.spectral_temperature.map(Temperature::value);
        if !positive(real_index)
            || !nonnegative(absorptive_index)
            || !positive(particle_density.value())
            || !positive(minimum)
            || !positive(maximum)
            || minimum > maximum
            || !positive(domain.maximum_volume_fraction)
            || domain.maximum_volume_fraction >= 1.0
            || source.trim().is_empty()
        {
            return Err(RadiationError::Invalid(
                "sourced Rayleigh indices, density and optical domain required",
            ));
        }
        // Divide by the norm twice rather than squaring it, retaining the
        // small absorptive contribution without overflowing a squared norm.
        let real = finite(real_index * real_index - absorptive_index * absorptive_index + 2.0)?;
        let imaginary = finite(2.0 * real_index * absorptive_index)?;
        let norm = finite(real.hypot(imaginary))?;
        if !positive(norm) {
            return Err(RadiationError::Invalid(
                "unresolvable Rayleigh optical denominator",
            ));
        }
        let spectral_factor = finite(18.0 * core::f64::consts::PI * (imaginary / norm) / norm)?;
        if absorptive_index > 0.0 && spectral_factor == 0.0 {
            return Err(RadiationError::Invalid(
                "positive soot absorption underflow",
            ));
        }
        Ok(Self {
            refractive_index: [real_index, absorptive_index],
            particle_density,
            domain,
            source,
            spectral_factor,
        })
    }

    /// Original `[real n, nonnegative absorptive k]` in `m = n - i k`.
    pub const fn refractive_index(&self) -> [f64; 2] {
        self.refractive_index
    }
    /// Original condensed-particle material density [kg/m³].
    pub const fn particle_density(&self) -> Density {
        self.particle_density
    }
    /// Original optical-source domain, without promotion to physical validation.
    pub const fn domain(&self) -> SootOpticalDomain {
        self.domain
    }
    /// Original optical/density source and its stated Rayleigh assumptions.
    pub fn source(&self) -> &str {
        &self.source
    }

    fn volume_fraction(&self, volume_fraction: f64) -> Result<(), RadiationError> {
        if !nonnegative(volume_fraction) || volume_fraction > self.domain.maximum_volume_fraction {
            return Err(RadiationError::Invalid(
                "soot volume fraction outside optical source domain",
            ));
        }
        Ok(())
    }

    fn temperature(&self, temperature: Temperature) -> Result<(), RadiationError> {
        let t = temperature.value();
        let [minimum, maximum] = self.domain.spectral_temperature.map(Temperature::value);
        if !positive(t) || t < minimum || t > maximum {
            return Err(RadiationError::Invalid(
                "spectral temperature outside soot optical source domain",
            ));
        }
        Ok(())
    }

    fn coefficient(
        &self,
        volume_fraction: f64,
        inverse_wavelength: f64,
    ) -> Result<AbsorptionCoefficient, RadiationError> {
        if volume_fraction == 0.0 || self.spectral_factor == 0.0 {
            return Ok(AbsorptionCoefficient::new(0.0));
        }
        let coefficient = finite(self.spectral_factor * volume_fraction * inverse_wavelength)?;
        if !positive(coefficient) {
            return Err(RadiationError::Invalid(
                "positive soot absorption underflow",
            ));
        }
        Ok(AbsorptionCoefficient::new(coefficient))
    }

    /// Spectral absorption at an explicitly supplied vacuum wavelength [m].
    /// This applies the same declared constant-index Rayleigh model. The caller
    /// must justify its use at this wavelength; this is absorption, not total
    /// extinction or a visibility prediction. No Planck temperature is used.
    pub fn spectral_coefficient(
        &self,
        volume_fraction: f64,
        wavelength: Length,
    ) -> Result<AbsorptionCoefficient, RadiationError> {
        self.volume_fraction(volume_fraction)?;
        if !positive(wavelength.value()) {
            return Err(RadiationError::Invalid(
                "positive finite soot spectral wavelength required",
            ));
        }
        let inverse = finite(1.0 / wavelength.value())?;
        self.coefficient(volume_fraction, inverse)
    }

    /// Allocation-free Planck mean from the actual condensed particulate volume
    /// fraction and a supplied blackbody spectral temperature.
    pub fn planck_from_volume_fraction(
        &self,
        volume_fraction: f64,
        spectral_temperature: Temperature,
    ) -> Result<PlanckSootAbsorption<'_>, RadiationError> {
        self.volume_fraction(volume_fraction)?;
        self.temperature(spectral_temperature)?;
        let concentration = finite(self.particle_density.value() * volume_fraction)?;
        if volume_fraction > 0.0 && concentration == 0.0 {
            return Err(RadiationError::Invalid(
                "positive soot concentration underflow",
            ));
        }
        let inverse = finite(
            PLANCK_WAVENUMBER_FACTOR * spectral_temperature.value() / SECOND_RADIATION_CONSTANT_M_K,
        )?;
        Ok(PlanckSootAbsorption {
            coefficient: self.coefficient(volume_fraction, inverse)?,
            mass_concentration: Density::new(concentration),
            volume_fraction,
            spectral_temperature,
            model: self,
        })
    }

    /// Planck mean from measured or transported soot mass per mixture volume.
    /// No soot yield or ideal-gas species/pressure contribution is inferred.
    pub fn planck_from_mass_concentration(
        &self,
        concentration: Density,
        spectral_temperature: Temperature,
    ) -> Result<PlanckSootAbsorption<'_>, RadiationError> {
        if !nonnegative(concentration.value()) {
            return Err(RadiationError::Invalid(
                "finite nonnegative soot mass concentration required",
            ));
        }
        let fraction = finite(concentration.value() / self.particle_density.value())?;
        if concentration.value() > 0.0 && fraction == 0.0 {
            return Err(RadiationError::Invalid(
                "positive soot volume fraction underflow",
            ));
        }
        let mut report = self.planck_from_volume_fraction(fraction, spectral_temperature)?;
        report.mass_concentration = concentration;
        Ok(report)
    }
}

impl PlanckSootAbsorption<'_> {
    /// Immutable model/source used to construct this mean.
    pub fn model(&self) -> &RayleighSoot {
        self.model
    }

    /// Explicitly choose this Planck mean as the soot-only gray layer opacity.
    /// Allocates one provenance string recording the model, weighting state and
    /// population. Molecular absorption must be supplied separately by its owner.
    /// Gas and soot must share a temperature if the existing one-temperature
    /// enclosure's gas energy account is used for this particulate emission.
    pub fn layer_absorption(&self) -> Result<LayerAbsorption, RadiationError> {
        LayerAbsorption::new(
            self.coefficient,
            format!(
                "constant-index Rayleigh soot Planck mean; {}; n={}; k={}; particle_density_kg_m3={}; spectral_temperature_K={}; soot_volume_fraction={}; soot_mass_concentration_kg_m3={}",
                self.model.source,
                self.model.refractive_index[0],
                self.model.refractive_index[1],
                self.model.particle_density.value(),
                self.spectral_temperature.value(),
                self.volume_fraction,
                self.mass_concentration.value(),
            ),
        )
    }
}

#[cfg(test)]
mod tests;
