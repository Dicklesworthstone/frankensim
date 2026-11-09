//! Saturated-water vapor calorics from the original IF97 source.
//!
//! A phase-changing cavity needs more than a frozen vapor pressure. Its mass
//! `m = rho_v(T) V` and energy `U = m u_v(T)` change with both temperature and
//! volume. This view retains the original liquid/vapor pair and differentiates
//! it along the ORIGINAL Region 4 saturation curve. It does not replace that
//! derivative with the Clapeyron relation: IF97's independently released
//! equations have a small, visible consistency defect.
//!
//! With incoming saturated-liquid enthalpy `h_l`, the open vapor first law is
//! `dU = dQ + h_l dm - p dV`. Its temperature coefficient per vapor volume is
//! `C_phase = rho_v du_v/dT + (u_v - h_l) d rho_v/dT`, and the volume
//! coefficient is `rho_v (h_v - h_l)`. These are analytic source coefficients,
//! not a heat-transfer law, evaporation-rate law or finite liquid inventory.
//! In particular, `C_phase` is NOT the isochoric heat capacity of a closed gas.
//!
//! Source: IAPWS R7-97(2012), Regions 1, 2 and 4:
//! <https://iapws.org/technical-guidance/release/IF97-Rev>.
//! This view retains the original 273.15..=623.15 K two-Gibbs-branch domain.
//! Consumers must separately admit density ratio, interface equilibrium,
//! mechanical pressure differences and their actual transport approximation.

use super::invalid;
use crate::MaterialError;
use fs_qty::Temperature;
use fs_thermochem::water::if97::{If97SaturationPair, saturation_pressure_derivative};

/// Matched liquid enthalpy and transport for resolved cavitation heat fields.
pub mod liquid;

/// One immutable original saturation pair and its allocation-free first jet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaturatedWaterVapor {
    saturation: If97SaturationPair,
    pressure_temperature_derivative_pa_k: f64,
    vapor_density_temperature_derivative_kg_m3_k: f64,
    vapor_internal_energy_temperature_derivative_j_kg_k: f64,
    phase_temperature_capacity_j_m3_k: f64,
}

impl SaturatedWaterVapor {
    /// Evaluate the original source at an absolute interface temperature [K].
    /// No transport parameters, gas EOS or phase-rate coefficient are inferred.
    pub fn at_temperature_k(temperature_k: f64) -> Result<Self, MaterialError> {
        let pair = If97SaturationPair::at_temperature(Temperature::new(temperature_k))
            .map_err(|error| invalid(&format!("water vapor saturation: {error}")))?;
        Self::try_from_saturation(pair)
    }

    /// Retain the supplied original pair without changing its caloric datum,
    /// pressure anchor or raw phase-equilibrium inconsistency.
    pub fn try_from_saturation(saturation: If97SaturationPair) -> Result<Self, MaterialError> {
        let vapor = saturation.vapor();
        let t = saturation.temperature().value();
        let p = saturation.pressure().value();
        let rho = vapor.density().value();
        let dp_dt = saturation_pressure_derivative(saturation.temperature())
            .map_err(|error| invalid(&format!("water saturation tangent: {error}")))?
            .value();
        let dv_dt = vapor.dv_dt_at_pressure().value();
        let dv_dp = vapor.dv_dp_at_temperature().value();
        let drho_dt = -rho * rho * (dv_dt + dv_dp * dp_dt);
        // u = h-pv, h_T|p = cp, h_p|T = v-T v_T|p.
        let du_dt = vapor.cp().value() - p * dv_dt - (t * dv_dt + p * dv_dp) * dp_dt;
        let capacity = rho * du_dt
            + (vapor.specific_internal_energy().value()
                - saturation.liquid().specific_enthalpy().value())
                * drho_dt;
        if ![dp_dt, drho_dt, du_dt, capacity]
            .into_iter()
            .all(f64::is_finite)
            || dp_dt <= 0.0
            || capacity <= 0.0
            || saturation.latent_enthalpy().value() <= 0.0
        {
            return Err(invalid(
                "unrepresentable or nonpositive water phase temperature capacity",
            ));
        }
        Ok(Self {
            saturation,
            pressure_temperature_derivative_pa_k: dp_dt,
            vapor_density_temperature_derivative_kg_m3_k: drho_dt,
            vapor_internal_energy_temperature_derivative_j_kg_k: du_dt,
            phase_temperature_capacity_j_m3_k: capacity,
        })
    }

    /// Original pair, including its raw latent enthalpy and Gibbs discrepancy.
    pub const fn saturation(&self) -> &If97SaturationPair {
        &self.saturation
    }
    /// Common saturation temperature [K].
    pub fn temperature_k(self) -> f64 {
        self.saturation.temperature().value()
    }
    /// Original absolute saturation pressure [Pa].
    pub fn pressure_pa(self) -> f64 {
        self.saturation.pressure().value()
    }
    /// Saturated liquid density [kg/m3]; not a metastable bulk-liquid EOS.
    pub fn liquid_density_kg_m3(self) -> f64 {
        self.saturation.liquid().density().value()
    }
    /// Saturated vapor density [kg/m3].
    pub fn vapor_density_kg_m3(self) -> f64 {
        self.saturation.vapor().density().value()
    }
    /// Vapor specific internal energy in the original IF97 datum [J/kg].
    pub fn vapor_internal_energy_j_kg(self) -> f64 {
        self.saturation.vapor().specific_internal_energy().value()
    }
    /// Saturated-liquid mass-transfer enthalpy in the SAME datum [J/kg].
    pub fn liquid_enthalpy_j_kg(self) -> f64 {
        self.saturation.liquid().specific_enthalpy().value()
    }
    /// Original vapor-minus-liquid specific enthalpy [J/kg].
    pub fn latent_enthalpy_j_kg(self) -> f64 {
        self.saturation.latent_enthalpy().value()
    }
    /// Original saturated-vapor sound speed [m/s].
    pub fn vapor_sound_speed_m_s(self) -> f64 {
        self.saturation.vapor().sound_speed().value()
    }
    /// Original Region 4 `dp_sat/dT` [Pa/K].
    pub const fn pressure_temperature_derivative_pa_k(self) -> f64 {
        self.pressure_temperature_derivative_pa_k
    }
    /// Total vapor-density derivative ALONG saturation [kg/(m3 K)].
    pub const fn vapor_density_temperature_derivative_kg_m3_k(self) -> f64 {
        self.vapor_density_temperature_derivative_kg_m3_k
    }
    /// Total vapor-specific-energy derivative ALONG saturation [J/(kg K)].
    pub const fn vapor_internal_energy_temperature_derivative_j_kg_k(self) -> f64 {
        self.vapor_internal_energy_temperature_derivative_j_kg_k
    }
    /// Open-system temperature coefficient `rho u' + (u-h_l) rho'` [J/(m3 K)].
    /// Multiply by vapor volume, not the total volume of a two-phase mixture.
    pub const fn phase_temperature_capacity_j_m3_k(self) -> f64 {
        self.phase_temperature_capacity_j_m3_k
    }
}

#[cfg(test)]
#[path = "vapor/tests.rs"]
mod tests;
