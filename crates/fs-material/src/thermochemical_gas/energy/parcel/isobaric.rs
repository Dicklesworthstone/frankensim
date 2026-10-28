//! Closed-parcel heating against a constant, quasistatic pressure load.
//!
//! With fixed mass/composition and only pressure-volume work, Q = Delta H
//! and Delta U + p Delta V = Q. Both balances are checked independently from
//! the original NASA-9 laws and the actual accepted volume. See NASA Glenn,
//! https://www.grc.nasa.gov/WWW/K-12/BGP/enthalpy.html.
//! The existing parcel remains the sole state owner; this is not piston
//! dynamics, an open-system mixing step or a chemical heat-release model.

use super::ClosedGasParcel;
use super::super::{GasEnergyError, GasEnergySolve, poll, positive, updated_point};
use crate::thermochemical_gas::failure;
use fs_qty::{Pressure, Temperature};

/// Accepted constant-pressure heat and boundary-work balance in coherent SI.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasIsobaricStep {
    /// Constant gas and opposing quasistatic load pressure [Pa].
    pub pressure_pa: f64,
    /// Conserved parcel mass [kg].
    pub mass_kg: f64,
    pub initial_temperature_k: f64,
    pub final_temperature_k: f64,
    pub initial_volume_m3: f64,
    pub final_volume_m3: f64,
    /// Signed heat entering the gas [J].
    pub supplied_heat_j: f64,
    /// Original-law m Delta h, with common formation references canceled [J].
    pub enthalpy_change_j: f64,
    /// Original-law m Delta u, not subtraction of rounded absolute energies [J].
    pub stored_energy_change_j: f64,
    /// Actual constant pressure times accepted volume change [J].
    /// Positive work is delivered to the external load.
    pub work_by_gas_j: f64,
    /// m Delta h - Q [J].
    pub enthalpy_residual_j: f64,
    /// m Delta u + p Delta V - Q [J].
    pub energy_residual_j: f64,
    /// Change of GAS entropy [J/K], not total entropy production. The caller
    /// must account for its thermal reservoirs before making a passivity claim.
    pub entropy_change_j_k: f64,
    /// Interior iterations in the shared nonlinear enthalpy inversion.
    pub iterations: usize,
}

impl ClosedGasParcel {
    /// Add heat while a quasistatic external load holds the CURRENT gas pressure.
    ///
    /// Frozen mass/composition, ideal-gas phase, no kinetic energy or shaft work.
    /// The original caloric law solves Delta h = Q/m; the same EOS determines
    /// the volume. Original-law Delta u and actual p Delta V independently close
    /// the first law within the caller's absolute Joule tolerance. No constant
    /// Cp or gamma is substituted. Full source references and extra query axes
    /// survive. Cooling contracts the gas and gives negative delivered work.
    ///
    /// A nonzero operation requires a stable continuous bracket containing the
    /// initial temperature. Source/transport refusal, missing representable
    /// temperature or volume change, failed balances and cancellation leave the
    /// original parcel unchanged. Zero heat is an exact no-op after control and
    /// cancellation checks. This finite thermodynamic operation does not supply
    /// a pressure-controller trajectory, piston acceleration or heat-transfer rate.
    pub fn add_heat_at_constant_pressure(
        &mut self, heat_j: f64, controls: GasEnergySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasIsobaricStep, GasEnergyError> {
        controls.validate()?;
        poll(&mut cancelled)?;
        if !heat_j.is_finite() {
            return Err(GasEnergyError::Invalid("finite signed isobaric heat required"));
        }
        let old = &self.parcel;
        let (mass, volume0) = (old.mass_kg, old.volume_m3);
        let (t0, pressure) = (old.state.state().temperature, old.state.state().pressure);
        let mut report = GasIsobaricStep {
            pressure_pa: pressure, mass_kg: mass,
            initial_temperature_k: t0, final_temperature_k: t0,
            initial_volume_m3: volume0, final_volume_m3: volume0,
            supplied_heat_j: heat_j, enthalpy_change_j: 0.0, stored_energy_change_j: 0.0,
            work_by_gas_j: 0.0, enthalpy_residual_j: 0.0, energy_residual_j: 0.0,
            entropy_change_j_k: 0.0, iterations: 0,
        };
        if heat_j == 0.0 { return Ok(report); }
        let change = heat_j / mass;
        let tolerance = mass * controls.energy_tolerance_j_kg;
        if !change.is_finite() || change == 0.0 || !positive(tolerance) {
            return Err(GasEnergyError::Invalid("isobaric heat or Joule tolerance is unrepresentable"));
        }
        let solved = old.model.reconstruct_pressure_enthalpy_increment(
            &old.point, change, controls, &mut cancelled,
        )?;
        let endpoint = solved.state().state();
        let t1 = endpoint.temperature;
        let volume1 = mass / endpoint.density;
        let work = pressure * (volume1 - volume0);
        let enthalpy = mass * solved.enthalpy_change_j_kg();
        let stored = mass * old.model.thermodynamics().specific_internal_energy_change(
            Temperature::new(t0), Temperature::new(t1),
        ).map_err(|error| failure("isobaric internal-energy change", error))?.value();
        let entropy = mass * old.model.thermodynamics().specific_entropy_change(
            Temperature::new(t0), Pressure::new(pressure),
            Temperature::new(t1), Pressure::new(endpoint.pressure),
        ).map_err(|error| failure("isobaric entropy change", error))?.value();
        let enthalpy_residual = enthalpy - heat_j;
        let energy_residual = (stored + work) - heat_j;
        if !positive(volume1) || volume1 == volume0 || t1 == t0
            || [work, enthalpy, stored, entropy, enthalpy_residual, energy_residual]
                .iter().any(|value| !value.is_finite())
            || [work, enthalpy, stored, entropy].iter().any(|&value|
                value == 0.0 || value.is_sign_positive() != heat_j.is_sign_positive())
            || enthalpy_residual.abs() > tolerance || energy_residual.abs() > tolerance
        {
            return Err(GasEnergyError::Invalid(
                "isobaric endpoint did not close represented enthalpy, storage and volume work",
            ));
        }
        let point = updated_point(&old.point, t1, endpoint.pressure)?;
        report.final_temperature_k = t1;
        report.final_volume_m3 = volume1;
        report.enthalpy_change_j = enthalpy;
        report.stored_energy_change_j = stored;
        report.work_by_gas_j = work;
        report.enthalpy_residual_j = enthalpy_residual;
        report.energy_residual_j = energy_residual;
        report.entropy_change_j_k = entropy;
        report.iterations = solved.iterations();
        poll(&mut cancelled)?;
        self.parcel.state = solved.into_state();
        self.parcel.point = point;
        self.parcel.volume_m3 = volume1;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
