//! Original liquid calorics matched to the saturated-vapor interface.
//!
//! The nonlinear vapor owner transfers the original saturated-liquid enthalpy
//! `h_l,sat(T)` with phase mass. A resolved surrounding-liquid approximation
//! must use that SAME enthalpy reference on its thermal and advective ports.
//! This view retains the original IF97 liquid branch, its industrial transport,
//! and the analytic saturation-path derivative
//! `dh_l,sat/dT = cp + (v - T dv/dT|p) dp_sat/dT`.
//!
//! Using this curve as an incompressible liquid caloric law is an explicit
//! near-saturation approximation. It neither assigns the actual liquid pressure
//! to `p_sat(T)` nor supplies a metastable EOS. Its capacity differs from the
//! physical fixed-pressure `cp`: the relative difference remains exposed for
//! consumer admission. Thermal expansion, pressure-motion heating and the
//! pressure departure from this reference curve require separate treatment or
//! an explicitly bounded approximation. No heat equation or integrator lives
//! in this immutable material view.
//!
//! Sources are the original IAPWS-IF97 Regions 1/2/4 and the existing original
//! water viscosity/conductivity implementation; no coefficient is fitted here.

use crate::MaterialError;
use crate::cavitation::invalid;
use fs_qty::Temperature;
use fs_thermochem::water::if97::{
    If97SaturationPair, regions::If97FluidState, saturation_pressure_derivative,
};
use fs_thermochem::water::transport::If97TransportState;

/// One original saturation pair and matching liquid transport/caloric tangent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturatedWaterLiquidThermal {
    saturation: If97SaturationPair,
    transport: If97TransportState,
    enthalpy_temperature_derivative_j_kg_k: f64,
    pressure_path_capacity_fraction: f64,
}

impl SaturatedWaterLiquidThermal {
    /// Sample the original reference curve without a root solve or allocation
    /// on the successful path. The shared thermodynamic/transport domain is
    /// 273.16..=623.15 K; the transport source refuses IF97's colder extension.
    pub fn at_temperature_k(temperature_k: f64) -> Result<Self, MaterialError> {
        let pair = If97SaturationPair::at_temperature(Temperature::new(temperature_k))
            .map_err(|error| invalid(&format!("water liquid thermal reference: {error}")))?;
        Self::try_from_saturation(pair)
    }

    /// Retain an already evaluated original phase pair and its caloric datum.
    pub fn try_from_saturation(saturation: If97SaturationPair) -> Result<Self, MaterialError> {
        let liquid = saturation.liquid();
        let transport = If97TransportState::from_state(If97FluidState::Gibbs(*liquid))
            .map_err(|error| invalid(&format!("water liquid thermal transport: {error}")))?;
        let pressure_tangent = saturation_pressure_derivative(saturation.temperature())
            .map_err(|error| invalid(&format!("water liquid saturation tangent: {error}")))?
            .value();
        let correction = liquid.dh_dp_at_temperature().value() * pressure_tangent;
        let capacity = liquid.cp().value() + correction;
        let fraction = correction.abs() / liquid.cp().value();
        if !capacity.is_finite() || capacity <= 0.0 || !fraction.is_finite() {
            return Err(invalid(
                "unrepresentable liquid saturation-path heat capacity",
            ));
        }
        Ok(Self {
            saturation,
            transport,
            enthalpy_temperature_derivative_j_kg_k: capacity,
            pressure_path_capacity_fraction: fraction,
        })
    }

    /// Original phase pair, preserving released endpoints and source mismatch.
    pub const fn saturation(&self) -> &If97SaturationPair {
        &self.saturation
    }
    /// Original transport evaluated at this exact saturated-liquid state.
    pub const fn transport(&self) -> &If97TransportState {
        &self.transport
    }
    pub fn temperature_k(self) -> f64 {
        self.saturation.temperature().value()
    }
    pub fn saturation_pressure_pa(self) -> f64 {
        self.saturation.pressure().value()
    }
    /// Original liquid enthalpy in the SAME datum used by vapor mass exchange.
    pub fn enthalpy_j_kg(self) -> f64 {
        self.saturation.liquid().specific_enthalpy().value()
    }
    /// Caloric slope along the reference curve, not fixed-pressure cp [J/(kg K)].
    pub const fn enthalpy_temperature_derivative_j_kg_k(self) -> f64 {
        self.enthalpy_temperature_derivative_j_kg_k
    }
    /// Original Gibbs heat capacity at fixed pressure [J/(kg K)].
    pub fn cp_j_kg_k(self) -> f64 {
        self.saturation.liquid().cp().value()
    }
    /// `|(dh_sat/dT)-cp|/cp`, a visible local approximation scale.
    pub const fn pressure_path_capacity_fraction(self) -> f64 {
        self.pressure_path_capacity_fraction
    }
    pub fn conductivity_w_m_k(self) -> f64 {
        self.transport.thermal_conductivity().value()
    }
    pub fn density_kg_m3(self) -> f64 {
        self.saturation.liquid().density().value()
    }
    pub fn dynamic_viscosity_pa_s(self) -> f64 {
        self.transport.dynamic_viscosity().value()
    }
    pub fn sound_speed_m_s(self) -> f64 {
        self.saturation.liquid().sound_speed().value()
    }
    /// Signed original volumetric thermal expansion [1/K].
    pub fn thermal_expansion_per_k(self) -> f64 {
        self.saturation.liquid().thermal_expansion().value()
    }
    /// Original isothermal compressibility [1/Pa].
    pub fn isothermal_compressibility_per_pa(self) -> f64 {
        self.saturation
            .liquid()
            .isothermal_compressibility()
            .value()
    }
}

#[cfg(test)]
#[path = "liquid/tests.rs"]
mod tests;
