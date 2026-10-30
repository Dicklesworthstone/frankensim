//! Source-backed heat conduction through ice, water and the fusion plateau.
//!
//! Caloric state remains the original `FreezingWater` state. Conductivity is an
//! independent physical property: liquid uses IAPWS SR6-08(2011), while ordinary
//! polycrystalline freshwater ice uses the empirical Yen (1981) correlation.
//! This composition is deliberately restricted to standard atmospheric pressure
//! and the intersection of the source temperature windows. The ice correlation
//! supplies neither a pressure correction nor an experimental uncertainty bound.
//!
//! A partially frozen cell needs an explicit microstructure assumption. Series
//! and parallel layers use the actual liquid VOLUME fraction derived from the
//! distinct source densities. They are different directional closures, not a
//! claim that conductivity follows uniquely from phase fraction. No apparent
//! heat capacity, additional latent heat or shifted caloric datum is introduced.

use super::{FreezingWater, FreezingWaterError, FreezingWaterState};
use fs_qty::{Pressure, Temperature};
use fs_thermal::enthalpy::DepositionError;
use fs_thermal::enthalpy::conduction::ConductivityLaw;
use fs_thermochem::water::conductivity::{
    WaterConductivityError, WaterConductivityState, ice_at_temperature,
    liquid_at_pressure_temperature,
};
use fs_thermochem::water::ice::NORMAL_PRESSURE_PA;

/// Directional closure for the material inside a mixed-phase cell.
/// The caller must choose; neither is an inferred isotropic microstructure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhaseConductivityMixture {
    /// Ice and liquid layers normal to heat flow: reciprocal conductivities add
    /// with their volume fractions. Appropriate as a declared unresolved planar
    /// front closure, without locating either phase inside a half-cell.
    Series,
    /// Ice and liquid layers parallel to heat flow: conductivities add with
    /// their volume fractions. This is the distinct parallel-layer assumption.
    Parallel,
}

/// Source, caloric, model-binding or cancellation refusal; no state is changed.
#[derive(Clone, Debug, PartialEq)]
pub enum FreezingConductivityError {
    Caloric(FreezingWaterError),
    Source(WaterConductivityError),
    /// No ice pressure dependence is invented outside the declared atmosphere.
    PressureOutsideModel {
        pressure_pa: f64,
    },
    ChartMismatch,
    Unrepresentable,
    Cancelled,
}
impl From<FreezingWaterError> for FreezingConductivityError {
    fn from(error: FreezingWaterError) -> Self {
        if error == FreezingWaterError::Cancelled {
            Self::Cancelled
        } else {
            Self::Caloric(error)
        }
    }
}
impl From<WaterConductivityError> for FreezingConductivityError {
    fn from(error: WaterConductivityError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for FreezingConductivityError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "water freezing conductivity: {self:?}")
    }
}
impl std::error::Error for FreezingConductivityError {}

/// Immutable original caloric chart and independently sourced transport law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingWaterConductivity {
    chart: FreezingWater,
    mixture: PhaseConductivityMixture,
    ice_at_melting: WaterConductivityState,
    liquid_at_melting: WaterConductivityState,
}

/// One material state and the precise phase data used for its conductivity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreezingConductivityState {
    material: FreezingWaterState,
    ice: Option<WaterConductivityState>,
    liquid: Option<WaterConductivityState>,
    conductivity_w_m_k: f64,
}
impl FreezingConductivityState {
    pub const fn material_state(&self) -> &FreezingWaterState {
        &self.material
    }
    /// Effective directional conductivity in W/(m K), with no added heat term.
    pub const fn conductivity_w_m_k(&self) -> f64 {
        self.conductivity_w_m_k
    }
    /// Absent phases do not acquire an extrapolated transport state.
    pub const fn ice_source(&self) -> Option<&WaterConductivityState> {
        self.ice.as_ref()
    }
    pub const fn liquid_source(&self) -> Option<&WaterConductivityState> {
        self.liquid.as_ref()
    }
}

impl FreezingWaterConductivity {
    /// Admit the complete caloric chart against both transport source windows.
    /// This preparation is constant work; melting-point phase values are cached
    /// for subsequent plateau queries. Source domain checks remain active on
    /// every sensible-phase query. No sampled conductivity bound is needed by
    /// the shared explicit owner, which refreshes each positive conductance at
    /// every actual thermal substep before computing its timestep restriction.
    pub fn new(
        chart: FreezingWater,
        mixture: PhaseConductivityMixture,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FreezingConductivityError> {
        if cancelled() {
            return Err(FreezingConductivityError::Cancelled);
        }
        if chart.pressure_pa() != NORMAL_PRESSURE_PA {
            return Err(FreezingConductivityError::PressureOutsideModel {
                pressure_pa: chart.pressure_pa(),
            });
        }
        let [minimum, maximum] = chart.temperature_bounds_k();
        let pressure = Pressure::new(chart.pressure_pa());
        let melting = chart.melting().temperature();
        ice_at_temperature(Temperature::new(minimum))?;
        let ice_at_melting = ice_at_temperature(melting)?;
        let liquid_at_melting = liquid_at_pressure_temperature(pressure, melting)?;
        liquid_at_pressure_temperature(pressure, Temperature::new(maximum))?;
        if cancelled() {
            return Err(FreezingConductivityError::Cancelled);
        }
        Ok(Self {
            chart,
            mixture,
            ice_at_melting,
            liquid_at_melting,
        })
    }

    pub const fn chart(&self) -> &FreezingWater {
        &self.chart
    }
    pub const fn mixture(&self) -> PhaseConductivityMixture {
        self.mixture
    }

    /// Resolve enthalpy with the same bounded original source inverse used by
    /// thermal storage. A topology adapter should evaluate this once per cell
    /// per substep, then share it across all incident conductive faces.
    pub fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FreezingConductivityState, FreezingConductivityError> {
        let material = self
            .chart
            .state_at_specific_enthalpy(specific_h_j_kg, &mut cancelled)?
            .into_state();
        self.state_at_material_state(material, cancelled)
    }

    /// Reuse an already resolved original caloric state without another inverse.
    /// Chart identity includes pressure, source matching and numerical domain.
    pub fn state_at_material_state(
        &self,
        material: FreezingWaterState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FreezingConductivityState, FreezingConductivityError> {
        if cancelled() {
            return Err(FreezingConductivityError::Cancelled);
        }
        if material.law() != &self.chart {
            return Err(FreezingConductivityError::ChartMismatch);
        }
        let x = material.liquid_mass_fraction();
        let temperature = Temperature::new(material.temperature_k());
        let (ice, liquid, conductivity_w_m_k) = if x == 0.0 {
            let source = if temperature == self.chart.melting().temperature() {
                self.ice_at_melting
            } else {
                ice_at_temperature(temperature)?
            };
            (Some(source), None, source.conductivity().value())
        } else if x == 1.0 {
            let source = if temperature == self.chart.melting().temperature() {
                self.liquid_at_melting
            } else {
                liquid_at_pressure_temperature(
                    Pressure::new(self.chart.pressure_pa()),
                    temperature,
                )?
            };
            (None, Some(source), source.conductivity().value())
        } else {
            let phi = material.liquid_volume_fraction();
            if !phi.is_finite() || !(0.0..=1.0).contains(&phi) {
                return Err(FreezingConductivityError::Unrepresentable);
            }
            let ki = self.ice_at_melting.conductivity().value();
            let kl = self.liquid_at_melting.conductivity().value();
            let k = match self.mixture {
                PhaseConductivityMixture::Series => 1.0 / ((1.0 - phi) / ki + phi / kl),
                PhaseConductivityMixture::Parallel => ki + phi * (kl - ki),
            };
            (Some(self.ice_at_melting), Some(self.liquid_at_melting), k)
        };
        if !conductivity_w_m_k.is_finite() || conductivity_w_m_k <= 0.0 {
            return Err(FreezingConductivityError::Unrepresentable);
        }
        if cancelled() {
            return Err(FreezingConductivityError::Cancelled);
        }
        Ok(FreezingConductivityState {
            material,
            ice,
            liquid,
            conductivity_w_m_k,
        })
    }
}

impl ConductivityLaw for FreezingWaterConductivity {
    fn conductivity_w_m_k(&self, specific_h_j_kg: f64) -> Result<f64, DepositionError> {
        self.state_at_specific_enthalpy(specific_h_j_kg, || false)
            .map(|state| state.conductivity_w_m_k())
            .map_err(|_| DepositionError::MaterialDomain)
    }
}

#[cfg(test)]
mod tests;
