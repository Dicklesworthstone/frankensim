//! A finite constant-capacity heat store on the original adaptive timeline.
//!
//! The reservoir keeps a compensated energy CHANGE from its original datum.
//! Its rounded published temperature is never converted back into stored heat.
//! Each implicit exchange eliminates the reservoir endpoint analytically;
//! only the original bath/phase scalar root and time controller are required.

use super::heat_transfer::{HeatReceipt, WindowHeat};
use super::{
    SolubilityCaloricErrorV1, SolubilityThermalKineticControlsV1,
    SolubilityThermalKineticOutcomeV1, SolubilityThermalKineticStateV1,
    SolubilityThermalKineticsV1, poll, positive,
};
use fs_qty::{Temperature, Time};

/// Immutable sensible reservoir with explicit constant heat capacity [J/K].
///
/// This is a caller-qualified lumped heat store, not inferred material Cp,
/// solvent inventory, a phase model or an infinite thermostat. Positive finite
/// temperature and capacity are required. Absolute C*T is never formed, so a
/// representable exchange does not require its arbitrary absolute datum to fit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SensibleHeatReservoirV1 {
    reference_temperature: Temperature,
    temperature: Temperature,
    heat_capacity_j_k: f64,
    energy: HeatReceipt,
    energy_change_j: f64,
}
impl SensibleHeatReservoirV1 {
    pub fn new(
        temperature: Temperature,
        heat_capacity_j_k: f64,
    ) -> Result<Self, SolubilityCaloricErrorV1> {
        if !positive(temperature.value()) || !positive(heat_capacity_j_k) {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "positive finite sensible reservoir temperature and heat capacity",
            ));
        }
        Ok(Self {
            reference_temperature: temperature,
            temperature,
            heat_capacity_j_k,
            energy: HeatReceipt::default(),
            energy_change_j: 0.0,
        })
    }

    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub const fn reference_temperature(&self) -> Temperature {
        self.reference_temperature
    }
    pub const fn heat_capacity_j_k(&self) -> f64 {
        self.heat_capacity_j_k
    }
    /// Stored sensible-energy change [J] from the original reservoir datum.
    /// Continuation retains its compensation even when this scalar getter or
    /// the published temperature rounds a small subsequent exchange away.
    pub const fn energy_change_j(&self) -> f64 {
        self.energy_change_j
    }

    /// Direct energy increment between states of the SAME capacity and datum.
    /// This uses retained energy, not C times rounded endpoint temperatures.
    pub fn enthalpy_change_j(&self, final_state: &Self) -> Result<f64, SolubilityCaloricErrorV1> {
        if self.reference_temperature.value().to_bits()
            != final_state.reference_temperature.value().to_bits()
            || self.heat_capacity_j_k.to_bits() != final_state.heat_capacity_j_k.to_bits()
        {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "reservoir energy increment requires the same capacity and original datum",
            ));
        }
        final_state.energy.subtract(self.energy)?.value()
    }

    fn energy_after(
        self,
        heat_to_bath: HeatReceipt,
    ) -> Result<HeatReceipt, SolubilityCaloricErrorV1> {
        self.energy.subtract(heat_to_bath)
    }

    fn temperature_after(
        self,
        heat_to_bath: HeatReceipt,
    ) -> Result<Temperature, SolubilityCaloricErrorV1> {
        let temperature = self
            .energy_after(heat_to_bath)?
            .shifted_quotient(self.heat_capacity_j_k, self.reference_temperature.value())?;
        if !positive(temperature) {
            return Err(SolubilityCaloricErrorV1::Arithmetic(
                "positive finite reservoir endpoint temperature is not representable",
            ));
        }
        Ok(Temperature::new(temperature))
    }

    fn after_heat(self, heat_to_bath: HeatReceipt) -> Result<Self, SolubilityCaloricErrorV1> {
        if heat_to_bath == HeatReceipt::default() {
            return Ok(self);
        }
        let energy = self.energy_after(heat_to_bath)?;
        if energy == self.energy && heat_to_bath.value()? != 0.0 {
            return Err(SolubilityCaloricErrorV1::Arithmetic(
                "reservoir energy increment is below retained precision",
            ));
        }
        Ok(Self {
            reference_temperature: self.reference_temperature,
            temperature: self.temperature_after(heat_to_bath)?,
            heat_capacity_j_k: self.heat_capacity_j_k,
            energy,
            energy_change_j: energy.value()?,
        })
    }
}

/// One accepted closed bath/reservoir window, with separate energy receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct SolubilityThermalReservoirOutcomeV1 {
    balance: SolubilityThermalKineticOutcomeV1,
    reservoir: SensibleHeatReservoirV1,
    reservoir_enthalpy_change_j: f64,
    closed_system_energy_residual_j: f64,
}
impl SolubilityThermalReservoirOutcomeV1 {
    pub const fn balance(&self) -> &SolubilityThermalKineticOutcomeV1 {
        &self.balance
    }
    pub const fn reservoir(&self) -> &SensibleHeatReservoirV1 {
        &self.reservoir
    }
    pub const fn reservoir_enthalpy_change_j(&self) -> f64 {
        self.reservoir_enthalpy_change_j
    }
    /// Direct bath enthalpy plus independently retained reservoir energy [J].
    pub const fn closed_system_energy_residual_j(&self) -> f64 {
        self.closed_system_energy_residual_j
    }
    pub fn into_parts(self) -> (SolubilityThermalKineticOutcomeV1, SensibleHeatReservoirV1) {
        (self.balance, self.reservoir)
    }
}

impl SolubilityThermalKineticsV1 {
    /// Exchange heat with a finite sensible store through conductance [W/K].
    ///
    /// The SAME native stages eliminate the reservoir endpoint analytically:
    /// `Q = G*dt/(1+G*dt/C) * (T_res,start-T_bath,end)`. The harmonic coefficient
    /// remains bounded by C, even when G*dt overflows. All accepted half stages
    /// carry independent heat receipts; the original full/two-half controller
    /// includes BOTH bath and reservoir temperature disagreement.
    ///
    /// Source support still applies to the bath; the separate explicit
    /// constant-capacity reservoir may start outside it. No external heat is
    /// supplied to the closed pair. Failure or cancellation changes neither
    /// input. Pass the returned reservoir into later windows to retain its
    /// original datum and sub-temperature-resolution energy changes.
    pub fn advance_with_finite_heat_reservoir(
        &self,
        initial: &SolubilityThermalKineticStateV1,
        duration: Time,
        reservoir: &SensibleHeatReservoirV1,
        conductance_w_k: f64,
        controls: SolubilityThermalKineticControlsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SolubilityThermalReservoirOutcomeV1, SolubilityCaloricErrorV1> {
        poll(&mut cancelled)?;
        if !conductance_w_k.is_finite() || conductance_w_k < 0.0 {
            return Err(SolubilityCaloricErrorV1::Invalid(
                "nonnegative finite reservoir heat-transfer conductance",
            ));
        }
        let (balance, heat) = if conductance_w_k == 0.0 {
            (
                self.advance(initial, duration, 0.0, controls, &mut cancelled)?,
                HeatReceipt::default(),
            )
        } else {
            self.advance_with_source(
                initial,
                duration,
                WindowHeat::Finite(FiniteReservoir {
                    reservoir: *reservoir,
                    conductance_w_k,
                }),
                controls,
                &mut cancelled,
            )?
        };
        let final_reservoir = reservoir.after_heat(heat)?;
        let reservoir_change = final_reservoir.energy.subtract(reservoir.energy)?;
        let reservoir_enthalpy_change_j = reservoir_change.value()?;
        let closed_system_energy_residual_j =
            reservoir_change.add(balance.enthalpy_change_j())?.value()?;
        if closed_system_energy_residual_j.abs() > controls.root.energy_tolerance_j {
            return Err(SolubilityCaloricErrorV1::EnergyResidual {
                residual_j: closed_system_energy_residual_j,
            });
        }
        poll(&mut cancelled)?;
        Ok(SolubilityThermalReservoirOutcomeV1 {
            balance,
            reservoir: final_reservoir,
            reservoir_enthalpy_change_j,
            closed_system_energy_residual_j,
        })
    }
}

#[derive(Clone, Copy)]
pub(super) struct FiniteReservoir {
    reservoir: SensibleHeatReservoirV1,
    conductance_w_k: f64,
}
impl FiniteReservoir {
    pub(super) fn heat(
        self,
        prefix: HeatReceipt,
        duration: Time,
        endpoint: Temperature,
    ) -> Result<f64, SolubilityCaloricErrorV1> {
        self.reservoir.temperature_after(prefix)?;
        let capacity = self.reservoir.heat_capacity_j_k;
        let difference = self.reservoir.energy_after(prefix)?.shifted_quotient(
            capacity,
            self.reservoir.reference_temperature.value() - endpoint.value(),
        )?;
        if difference == 0.0 || duration.value() == 0.0 {
            return Ok(0.0);
        }
        let g = self.conductance_w_k;
        let dt = duration.value();
        let exposure = g * dt;
        let value = if exposure <= capacity {
            let divisor = 1.0 + exposure / capacity;
            // This order also retains physical heat when G*dt underflows but
            // its product with the temperature difference is representable.
            let reduced_g = g / divisor;
            let candidates = [
                (reduced_g * dt) * difference,
                (reduced_g * difference) * dt,
                (dt * difference) * reduced_g,
            ];
            let direct = candidates.into_iter().find(|x| x.is_finite() && *x != 0.0);
            if reduced_g.is_normal() && direct.is_some() {
                direct
            } else {
                // Reducing a subnormal G first can erase or coarsely round
                // a perfectly representable complete physical heat product.
                [
                    (g * (dt / divisor)) * difference,
                    (g * difference) * (dt / divisor),
                    exposure * (difference / divisor),
                ]
                .into_iter()
                .find(|x| x.is_finite() && *x != 0.0)
                .or(direct)
            }
        } else {
            let inverse = [
                (capacity / g) / dt,
                (capacity / dt) / g,
                capacity / exposure,
            ]
            .into_iter()
            .find(|x| x.is_finite() && *x > 0.0)
            .unwrap_or(0.0);
            let effective = capacity / (1.0 + inverse);
            let value = effective * difference;
            if effective.is_normal() && value.is_finite() && value != 0.0 {
                Some(value)
            } else {
                [
                    (capacity * difference) / (1.0 + inverse),
                    capacity * (difference / (1.0 + inverse)),
                    value,
                ]
                .into_iter()
                .find(|x| x.is_finite() && *x != 0.0)
            }
        };
        value.ok_or(SolubilityCaloricErrorV1::Arithmetic(
            "nonzero finite-reservoir stage heat is not representable",
        ))
    }

    pub(super) fn temporal_error(
        self,
        full: HeatReceipt,
        fine: HeatReceipt,
        controls: SolubilityThermalKineticControlsV1,
    ) -> Result<f64, SolubilityCaloricErrorV1> {
        let full_temperature = self.reservoir.temperature_after(full)?.value();
        let fine_temperature = self.reservoir.temperature_after(fine)?.value();
        // Independent receipts preserve this difference even when both
        // published reservoir temperatures round to the same f64 value.
        let delta = full
            .subtract(fine)?
            .divided_by(self.reservoir.heat_capacity_j_k)?
            .abs();
        if delta == 0.0 {
            return Ok(0.0);
        }
        let absolute = controls.temperature_absolute_tolerance_k;
        let magnitude = full_temperature.max(fine_temperature);
        let scale = magnitude.max(absolute);
        let error = (delta / scale)
            / (absolute / scale + controls.relative_tolerance * (magnitude / scale));
        if error.is_nan() {
            Err(SolubilityCaloricErrorV1::Arithmetic(
                "finite reservoir temporal error",
            ))
        } else {
            // An infinite disagreement is a valid request to refine under
            // the original controller, not an arithmetic stage failure.
            Ok(error)
        }
    }
}

#[cfg(test)]
mod tests;
