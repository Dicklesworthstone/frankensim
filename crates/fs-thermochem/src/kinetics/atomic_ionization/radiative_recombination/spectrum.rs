//! Differential Milne emission and conditional tails on the original source.
//!
//! The event coefficient is per photon energy, not per frequency or solid
//! angle. Multiplication by electron and product-ion number densities gives
//! the angle-integrated spontaneous photon emissivity. Angular transport and
//! populations remain the consumer's responsibility.

use super::{
    AtomicRadiativeRecombinationV1, Error, LIGHT_SPEED_M_S,
    MAX_RADIATIVE_RECOMBINATION_THERMAL_RATIO_V1,
};
use crate::electron_collision::{AVOGADRO_PER_MOL, BOLTZMANN_J_PER_K, ELECTRON_MASS_KG};
use fs_math::det::{exp, ln};
use fs_qty::{Area, Energy, Temperature};

/// Original reciprocal source prepared at one actual electron temperature.
#[derive(Clone, Copy, Debug)]
pub struct RecombinationSpectrumV1<'a> {
    source: &'a AtomicRadiativeRecombinationV1,
    thermal: f64,
    gap: f64,
    log_normalization: f64,
}

/// Differential event and energy coefficients at one original photon energy.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RecombinationSpectralDensityV1 {
    /// Differential capture coefficient [m³ s⁻¹ J⁻¹].
    pub coefficient_m3_s_per_j: f64,
    /// Captured electron kinetic-energy density [J m³ s⁻¹ J⁻¹ = m³/s].
    pub captured_kinetic_energy_m3_s: f64,
    /// Original photoionization cross section at this same energy [m²].
    pub photoionization_cross_section: Area,
}

/// Absolute omitted-tail bounds, conditional on the retained global envelope.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct RecombinationSpectrumTailV1 {
    /// Event coefficient above the requested energy [m³/s].
    pub coefficient_m3_s: f64,
    /// Captured electron energy above that energy [J m³/s].
    pub captured_kinetic_energy_j_m3_s: f64,
    /// Emitted photon energy above that energy [J m³/s].
    pub emitted_photon_energy_j_m3_s: f64,
}

impl AtomicRadiativeRecombinationV1 {
    /// Prepare the differential Milne law without extrapolating the original
    /// cross section. No quadrature, angular distribution, or escape is implied.
    pub fn spectrum(
        &self,
        electron_temperature: Temperature,
    ) -> Result<RecombinationSpectrumV1<'_>, Error> {
        let channel = self.source().channel();
        channel
            .admit_source_temperature(channel.electron_index(), electron_temperature)
            .map_err(Error::Atomic)?;
        let thermal = BOLTZMANN_J_PER_K * electron_temperature.value();
        if !thermal.is_finite()
            || thermal <= 0.0
            || thermal / (ELECTRON_MASS_KG * LIGHT_SPEED_M_S * LIGHT_SPEED_M_S)
                > MAX_RADIATIVE_RECOMBINATION_THERMAL_RATIO_V1
        {
            return Err(Error::Invalid(
                "nonrelativistic Maxwellian thermal-energy limit",
            ));
        }
        let weights = channel.statistical_weights();
        let log_normalization = ln(f64::from(weights.neutral) / f64::from(weights.ion))
            + 0.5 * ln(2.0 / core::f64::consts::PI)
            - 1.5 * ln(ELECTRON_MASS_KG)
            - 2.0 * ln(LIGHT_SPEED_M_S)
            - 1.5 * ln(thermal);
        if !log_normalization.is_finite() {
            return Err(Error::Arithmetic("spectral Milne normalization"));
        }
        Ok(RecombinationSpectrumV1 {
            source: self,
            thermal,
            gap: channel.ionization_energy_j_mol() / AVOGADRO_PER_MOL,
            log_normalization,
        })
    }
}

impl<'a> RecombinationSpectrumV1<'a> {
    /// Original source including its finite cross-section interval and envelope.
    pub const fn source(&self) -> &'a AtomicRadiativeRecombinationV1 {
        self.source
    }

    /// k_B T_e [J] used by this immutable prepared spectrum.
    pub const fn thermal_energy(&self) -> Energy {
        Energy::new(self.thermal)
    }

    /// Evaluate C E² sigma(E) exp[-(E-chi)/(kT)]/(kT)^(3/2).
    /// The original onset, including a rounded onset above chi, is retained.
    /// Positive underflow refuses; consumers should first bound and omit an
    /// insignificant far tail, rather than sample it as artificial zero data.
    pub fn evaluate(&self, photon_energy: Energy) -> Result<RecombinationSpectralDensityV1, Error> {
        let sigma = self
            .source
            .source()
            .cross_section(photon_energy)
            .map_err(Error::Photoionization)?;
        if sigma.value() == 0.0 {
            return Ok(RecombinationSpectralDensityV1 {
                photoionization_cross_section: sigma,
                ..RecombinationSpectralDensityV1::default()
            });
        }
        let energy = photon_energy.value();
        let kinetic = (energy - self.gap).max(0.0);
        let log_density =
            self.log_normalization + 2.0 * ln(energy) + ln(sigma.value()) - kinetic / self.thermal;
        let density = exp(log_density);
        let captured = kinetic * density;
        if !density.is_finite()
            || density <= 0.0
            || !captured.is_finite()
            || (kinetic > 0.0 && captured == 0.0)
        {
            return Err(Error::Arithmetic("positive spectral Milne density"));
        }
        Ok(RecombinationSpectralDensityV1 {
            coefficient_m3_s_per_j: density,
            captured_kinetic_energy_m3_s: captured,
            photoionization_cross_section: sigma,
        })
    }

    /// Integrate the constant global cross-section envelope from `lower` to
    /// infinity analytically. This bounds omitted photons even ABOVE original
    /// source coverage; it never evaluates or extends that cross section.
    ///
    /// Below the original onset the lower limit is raised to that onset.
    /// Positive bounds include an outward roundoff allowance and retain the
    /// least positive f64 on underflow; they are conditional numerical bounds,
    /// not an interval or physical-source uncertainty certificate.
    pub fn tail_bound(&self, lower: Energy) -> Result<RecombinationSpectrumTailV1, Error> {
        if !lower.value().is_finite() || lower.value() < 0.0 {
            return Err(Error::Invalid("spectral tail lower energy"));
        }
        let bound = self.source.global_cross_section_upper_bound().value();
        if bound == 0.0 {
            return Ok(RecombinationSpectrumTailV1::default());
        }
        let energy = lower.value().max(self.source.source().threshold().value());
        let kinetic = energy - self.gap;
        let photon_scale = energy + self.thermal;
        let kinetic_scale = kinetic + self.thermal;
        if !photon_scale.is_finite() || !kinetic_scale.is_finite() || kinetic < 0.0 {
            return Err(Error::Arithmetic("spectral tail scales"));
        }
        let p = energy / photon_scale;
        let q = self.thermal / photon_scale;
        let k = kinetic / kinetic_scale;
        let r = self.thermal / kinetic_scale;
        let event = p * p + 2.0 * p * q + 2.0 * q * q;
        let captured = k * event + r * (p * p + 4.0 * p * q + 6.0 * q * q);
        let terms = [
            self.log_normalization,
            ln(bound),
            ln(self.thermal),
            2.0 * ln(photon_scale),
            -kinetic / self.thermal,
        ];
        let log_scale: f64 = terms.iter().sum();
        let magnitude: f64 = terms.iter().map(|value| value.abs()).sum();
        let outward = |extra: f64| {
            let log_value = log_scale + extra;
            if log_value == f64::NEG_INFINITY {
                return f64::from_bits(1);
            }
            // Scale the allowance to the log arithmetic, not only to the
            // final exponential, whose small result can hide cancellation.
            let padding = 16.0 * f64::EPSILON * (1.0 + magnitude + extra.abs());
            exp(log_value + padding).max(f64::from_bits(1))
        };
        let coefficient_m3_s = outward(ln(event));
        let captured_kinetic_energy_j_m3_s = outward(ln(kinetic_scale) + ln(captured));
        let emitted_photon_energy_j_m3_s = (self.gap * coefficient_m3_s
            + captured_kinetic_energy_j_m3_s)
            * (1.0 + 8.0 * f64::EPSILON);
        if [
            coefficient_m3_s,
            captured_kinetic_energy_j_m3_s,
            emitted_photon_energy_j_m3_s,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(Error::Arithmetic("spectral tail bounds"));
        }
        Ok(RecombinationSpectrumTailV1 {
            coefficient_m3_s,
            captured_kinetic_energy_j_m3_s,
            emitted_photon_energy_j_m3_s,
        })
    }
}
