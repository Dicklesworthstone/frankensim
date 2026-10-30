//! Freezing/thawing on the existing conservative enthalpy timeline.
//!
//! Preparation proves positive Cp bounds over both complete sensible branches
//! of the immutable pressure chart. The physical fusion plateau contributes
//! zero dT/dh, so the maximum reciprocal Cp is a global thermal CFL bound.
//! The same original enthalpy and phase owner supplies T(h) and integral dh/T.
//!
//! Pressure is externally maintained. Timeline geometry is reference geometry
//! for heat exchange; the material state separately reports actual ice/liquid
//! volumes and pressure work. This adapter adds neither a rigid vessel model
//! nor an ice-front, contact-conductance, supercooling or nucleation law.

use super::FreezingWater;
use core::fmt;
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};
use fs_thermochem::water::ice::bounds::{
    IceIhHeatCapacityBound, IceIhHeatCapacityBoundError, IceIhHeatCapacityBoundSolve,
    isobaric_heat_capacity_bound as ice_heat_capacity_bound,
};
use fs_thermochem::water::if97::If97Region;
use fs_thermochem::water::if97::bounds::{
    If97HeatCapacityBound, If97HeatCapacityBoundError, If97HeatCapacityBoundSolve,
    isobaric_heat_capacity_bound as liquid_heat_capacity_bound,
};

/// Independent bounded preparation work for the two original source equations.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FreezingThermalControls {
    pub ice_heat_capacity: IceIhHeatCapacityBoundSolve,
    pub liquid_heat_capacity: If97HeatCapacityBoundSolve,
}

/// A failed source bound or cancellation yields no prepared thermal law.
#[derive(Clone, Debug, PartialEq)]
pub enum FreezingThermalError {
    IceHeatCapacity(IceIhHeatCapacityBoundError),
    LiquidHeatCapacity(If97HeatCapacityBoundError),
    Cancelled,
}
impl fmt::Display for FreezingThermalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "water freezing thermal preparation: {self:?}")
    }
}
impl std::error::Error for FreezingThermalError {}
impl From<IceIhHeatCapacityBoundError> for FreezingThermalError {
    fn from(error: IceIhHeatCapacityBoundError) -> Self {
        if error == IceIhHeatCapacityBoundError::Cancelled {
            Self::Cancelled
        } else {
            Self::IceHeatCapacity(error)
        }
    }
}
impl From<If97HeatCapacityBoundError> for FreezingThermalError {
    fn from(error: If97HeatCapacityBoundError) -> Self {
        if error == If97HeatCapacityBoundError::Cancelled {
            Self::Cancelled
        } else {
            Self::LiquidHeatCapacity(error)
        }
    }
}

/// Original freezing chart and reusable bounds for the shared thermal solver.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingWaterThermalLaw {
    chart: FreezingWater,
    ice: Option<IceIhHeatCapacityBound>,
    liquid: Option<If97HeatCapacityBound>,
}
impl FreezingWaterThermalLaw {
    /// Prepare each sensible interval once; a latent-only chart needs neither.
    pub fn new(
        chart: FreezingWater,
        controls: FreezingThermalControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FreezingThermalError> {
        if cancelled() {
            return Err(FreezingThermalError::Cancelled);
        }
        let [minimum, maximum] = chart.temperature_bounds_k();
        let melting = chart.melting().temperature();
        let pressure = Pressure::new(chart.pressure_pa());
        let ice = if minimum < melting.value() {
            Some(ice_heat_capacity_bound(
                pressure,
                Temperature::new(minimum),
                melting,
                controls.ice_heat_capacity,
                &mut cancelled,
            )?)
        } else {
            None
        };
        let liquid = if maximum > melting.value() {
            Some(liquid_heat_capacity_bound(
                pressure,
                If97Region::Region1,
                melting,
                Temperature::new(maximum),
                controls.liquid_heat_capacity,
                &mut cancelled,
            )?)
        } else {
            None
        };
        if cancelled() {
            return Err(FreezingThermalError::Cancelled);
        }
        Ok(Self { chart, ice, liquid })
    }

    pub const fn chart(&self) -> &FreezingWater {
        &self.chart
    }
    pub const fn ice_heat_capacity_bound(&self) -> Option<&IceIhHeatCapacityBound> {
        self.ice.as_ref()
    }
    pub const fn liquid_heat_capacity_bound(&self) -> Option<&If97HeatCapacityBound> {
        self.liquid.as_ref()
    }
    /// Conservative global upper bound on dT/dh, reused without source queries.
    pub fn maximum_temperature_slope(&self) -> f64 {
        self.ice
            .as_ref()
            .map(IceIhHeatCapacityBound::maximum_temperature_slope)
            .into_iter()
            .chain(
                self.liquid
                    .as_ref()
                    .map(If97HeatCapacityBound::maximum_temperature_slope),
            )
            .fold(0.0, f64::max)
    }
}

impl ThermalLaw for FreezingWaterThermalLaw {
    fn temperature_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        self.chart
            .state_at_specific_enthalpy(specific_h_j_kg, || false)
            .map(|solution| solution.state().temperature_k())
            .map_err(|_| DepositionError::MaterialDomain)
    }

    fn maximum_temperature_slope(&self) -> f64 {
        Self::maximum_temperature_slope(self)
    }

    fn specific_entropy_change_j_kg_k(
        &self,
        initial_h_j_kg: f64,
        final_h_j_kg: f64,
    ) -> Result<Option<f64>, DepositionError> {
        self.chart
            .specific_entropy_change_j_kg_k(initial_h_j_kg, final_h_j_kg, || false)
            .map(Some)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

#[cfg(test)]
mod tests;
