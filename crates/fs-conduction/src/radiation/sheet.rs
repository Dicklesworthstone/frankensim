//! Passive, reciprocal two-sided diffuse sheets with explicit optical sources.
//!
//! Both faces share one temperature and the same external refractive-index /
//! radiometric convention. A sheet may absorb and reflect differently on its
//! two sides, but has one reciprocal hemispherical transmittance. Kirchhoff's
//! law sets each side's emissivity to that side's absorptance.
//!
//! For irradiation H and blackbody emission B, the outgoing radiosity is
//! J_i = a_i B + r_i H_i + t H_other. Only a_i (H_i - B) heats the sheet;
//! H_i - J_i also contains through-sheet transport and is not a per-face heat
//! source. Equal-area consumers must sum the two thermal face powers once.
//!
//! These are explicit gray or band-constant diffuse data, not an inference
//! from visible alpha, total emissivity or a refractive index. Specular paths,
//! coherent films, angular optics and non-isothermal through-thickness transfer
//! need their respective transport owners.

use std::sync::Arc;

use super::{BlackbodyBand, STEFAN_BOLTZMANN_W_M2_K4, radiation_error};
use crate::ConductionError;

/// Source-bearing reciprocal optics, constant on the declared temperature
/// interval and, for spectral use, the explicitly supplied wavelength band.
#[derive(Clone, Debug, PartialEq)]
pub struct ReciprocalSheetOptics {
    absorptance: [f64; 2],
    transmittance: f64,
    temperature_support_k: [f64; 2],
    source: Arc<str>,
}

/// Constitutive energy accounts per unit area of each physical face.
/// Irradiation and every reported flux use the same gray or finite-band scope.
#[derive(Clone, Debug, PartialEq)]
pub struct SheetExchange {
    pub optics: ReciprocalSheetOptics,
    pub temperature_k: f64,
    /// None means explicitly gray; Some retains the actual finite band.
    pub wavelengths_m: Option<[f64; 2]>,
    pub irradiation_w_m2: [f64; 2],
    pub blackbody_w_m2: f64,
    /// No finite-window renormalization is performed. A spectral consumer must
    /// enforce its own omitted-spectrum budget over every emitting source.
    pub omitted_blackbody_fraction: f64,
    /// Outgoing emitted + reflected + radiation transmitted from the OTHER side.
    pub outgoing_w_m2: [f64; 2],
    pub absorbed_w_m2: [f64; 2],
    pub emitted_w_m2: [f64; 2],
    pub reflected_w_m2: [f64; 2],
    /// Indexed by the incident face: this power exits the opposite face.
    pub transmitted_w_m2: [f64; 2],
    /// Signed absorbed-minus-emitted heat, positive into the material.
    pub net_heat_w_m2: [f64; 2],
    /// Analytic derivative at fixed irradiation and fixed band endpoints.
    pub net_heat_temperature_derivative_w_m2_k: [f64; 2],
}

impl ReciprocalSheetOptics {
    /// Admit a passive partition without silently repairing source values.
    ///
    /// Each absorptance must fit in the representable remainder 1 - t, so
    /// deriving reflectance cannot produce a negative coefficient. Temperature
    /// support is a finite positive interval, inclusive at evaluation.
    pub fn new(
        absorptance: [f64; 2],
        transmittance: f64,
        temperature_support_k: [f64; 2],
        source: impl Into<String>,
    ) -> Result<Self, ConductionError> {
        let source = source.into();
        let [low, high] = temperature_support_k;
        if !transmittance.is_finite()
            || !(0.0..=1.0).contains(&transmittance)
            || absorptance
                .iter()
                .any(|&a| !a.is_finite() || a < 0.0 || a > 1.0 - transmittance)
            || !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high <= low
            || source.trim().is_empty()
        {
            return Err(invalid(
                "need a passive [absorptance, reflectance, transmittance] partition, a positive temperature interval and an explicit source",
            ));
        }
        Ok(Self {
            absorptance,
            transmittance,
            temperature_support_k,
            source: source.into(),
        })
    }

    /// Kirchhoff absorptance and emissivity in face order.
    pub const fn absorptance(&self) -> [f64; 2] {
        self.absorptance
    }

    pub const fn transmittance(&self) -> f64 {
        self.transmittance
    }

    pub fn reflectance(&self) -> [f64; 2] {
        let reflected_or_absorbed = 1.0 - self.transmittance;
        self.absorptance.map(|a| reflected_or_absorbed - a)
    }

    pub const fn temperature_support_k(&self) -> [f64; 2] {
        self.temperature_support_k
    }

    pub fn source(&self) -> &str {
        &self.source
    }

    pub fn validate_temperature(&self, temperature_k: f64) -> Result<(), ConductionError> {
        let [low, high] = self.temperature_support_k;
        if !temperature_k.is_finite() || !(low..=high).contains(&temperature_k) {
            return Err(invalid(
                "sheet temperature is outside its optical source support",
            ));
        }
        Ok(())
    }

    /// Explicit full-spectrum gray evaluation with two prescribed irradiations.
    pub fn evaluate_gray(
        &self,
        temperature_k: f64,
        irradiation_w_m2: [f64; 2],
    ) -> Result<SheetExchange, ConductionError> {
        self.validate_temperature(temperature_k)?;
        // Multiply the small SI constant first: T^4 can overflow before the
        // physically representable product sigma*T^4 does.
        let b = (((STEFAN_BOLTZMANN_W_M2_K4 * temperature_k) * temperature_k) * temperature_k)
            * temperature_k;
        let db = ((4.0 * STEFAN_BOLTZMANN_W_M2_K4 * temperature_k) * temperature_k) * temperature_k;
        self.evaluate(temperature_k, irradiation_w_m2, b, db, None, 0.0)
    }

    /// Evaluate one explicit finite wavelength band using the shared Planck
    /// primitive, including its analytic temperature derivative and coverage.
    /// The source must actually support constant optics on this band.
    pub fn evaluate_band(
        &self,
        band: BlackbodyBand,
        temperature_k: f64,
        irradiation_w_m2: [f64; 2],
    ) -> Result<SheetExchange, ConductionError> {
        self.validate_temperature(temperature_k)?;
        let b = band.evaluate(temperature_k)?;
        self.evaluate(
            temperature_k,
            irradiation_w_m2,
            b.power_w_m2,
            b.temperature_derivative_w_m2_k,
            Some([band.lower_wavelength_m(), band.upper_wavelength_m()]),
            b.short_wave_fraction + b.long_wave_fraction,
        )
    }

    fn evaluate(
        &self,
        temperature_k: f64,
        irradiation: [f64; 2],
        blackbody: f64,
        blackbody_derivative: f64,
        wavelengths_m: Option<[f64; 2]>,
        omitted_blackbody_fraction: f64,
    ) -> Result<SheetExchange, ConductionError> {
        if irradiation
            .iter()
            .chain([&blackbody, &blackbody_derivative])
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "irradiation and blackbody emission must be finite and nonnegative",
            ));
        }
        let reflectance = self.reflectance();
        let absorbed = std::array::from_fn(|i| self.absorptance[i] * irradiation[i]);
        let emitted = self.absorptance.map(|a| a * blackbody);
        let reflected = std::array::from_fn(|i| reflectance[i] * irradiation[i]);
        let transmitted = irradiation.map(|h| self.transmittance * h);
        let outgoing = std::array::from_fn(|i| emitted[i] + reflected[i] + transmitted[1 - i]);
        // Subtract fields before scaling to retain near-equilibrium heat.
        // In particular, zero absorption yields exactly zero material heating.
        let net_heat = std::array::from_fn(|i| self.absorptance[i] * (irradiation[i] - blackbody));
        let tangent = self.absorptance.map(|a| -a * blackbody_derivative);
        if outgoing
            .iter()
            .chain(net_heat.iter())
            .chain(tangent.iter())
            .any(|v| !v.is_finite())
        {
            return Err(invalid("unrepresentable sheet energy accounts"));
        }
        Ok(SheetExchange {
            optics: self.clone(),
            temperature_k,
            wavelengths_m,
            irradiation_w_m2: irradiation,
            blackbody_w_m2: blackbody,
            omitted_blackbody_fraction,
            outgoing_w_m2: outgoing,
            absorbed_w_m2: absorbed,
            emitted_w_m2: emitted,
            reflected_w_m2: reflected,
            transmitted_w_m2: transmitted,
            net_heat_w_m2: net_heat,
            net_heat_temperature_derivative_w_m2_k: tangent,
        })
    }
}

fn invalid(reason: &str) -> ConductionError {
    radiation_error(
        "reciprocal-sheet-optics",
        reason,
        "supply sourced passive two-sided diffuse optics with one supported sheet temperature",
    )
}

#[cfg(test)]
mod tests;
