//! Finite conductance to an explicitly prescribed, constant-temperature bath.
//!
//! The homogeneous parcel exchanges `Qdot = G (Tb - T)` with this reservoir
//! while its original reaction laws evolve. The bath has no finite heat
//! capacity or evolving state; its entropy change is `-Q/Tb`. Conductance is
//! supplied by the caller, not inferred from a material name or geometry.

use super::*;

/// Original-source finite caloric bodies, including sharp latent plateaus.
#[cfg(feature = "material")]
pub mod material;

/// Total parcel-to-reservoir thermal conductance, W/K (not conductivity).
pub type ThermalConductance = fs_qty::Qty<2, 1, -3, -1, 0, 0>;

/// Passive, finite thermal exchange with one held physical reservoir.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatBath {
    conductance: ThermalConductance,
    temperature: Temperature,
}

/// Shared thermal constraint; both forms enter the same implicit gas stages.
#[derive(Clone, Copy)]
pub(super) enum HeatExchange<'a> {
    Bath(&'a HeatBath),
    #[cfg(feature = "material")]
    Body(material::FiniteHeatBody<'a>),
}

// Finite-body state parameters disappear from the default bath-only match.
#[cfg_attr(not(feature = "material"), allow(unused_variables))]
impl HeatExchange<'_> {
    fn conductance(self) -> f64 {
        match self {
            Self::Bath(bath) => bath.conductance.value(),
            #[cfg(feature = "material")]
            Self::Body(body) => body.conductance().value(),
        }
    }

    pub fn increment_density(
        self,
        step: f64,
        temperature: f64,
        heat: Sum,
        volume: f64,
    ) -> Result<f64, IsochoricKineticsError> {
        match self {
            Self::Bath(bath) => Ok(bath.increment_density(step, temperature, volume)),
            #[cfg(feature = "material")]
            Self::Body(body) => Ok((step * body.conductance().value() / volume)
                * (body.sample(heat.value() * volume)?.temperature_k - temperature)),
        }
    }

    pub fn entropy_transfer(
        self,
        initial: Sum,
        final_heat: Sum,
        volume: f64,
    ) -> Result<f64, IsochoricKineticsError> {
        match self {
            Self::Bath(bath) => {
                Ok(final_heat.difference(initial) * volume / bath.temperature.value())
            }
            #[cfg(feature = "material")]
            Self::Body(body) => body.entropy_transfer(initial, final_heat, volume),
        }
    }

    pub fn additional_error(
        self,
        base: Sum,
        full: Sum,
        fine: Sum,
        volume: f64,
        controls: IsochoricKineticsSolve,
    ) -> Result<f64, IsochoricKineticsError> {
        match self {
            Self::Bath(_) => Ok(0.0),
            #[cfg(feature = "material")]
            Self::Body(body) => body.additional_error(base, full, fine, volume, controls),
        }
    }

    pub fn admit_heat(
        self,
        heat: Sum,
        volume: f64,
        tolerance: f64,
    ) -> Result<(), IsochoricKineticsError> {
        match self {
            Self::Bath(_) => Ok(()),
            #[cfg(feature = "material")]
            Self::Body(body) => body.admit_heat(heat.value() * volume, tolerance),
        }
    }
}

impl HeatBath {
    pub fn new(
        conductance: ThermalConductance,
        temperature: Temperature,
    ) -> Result<Self, IsochoricKineticsError> {
        if !conductance.value().is_finite()
            || conductance.value() < 0.0
            || !positive(temperature.value())
        {
            return Err(IsochoricKineticsError::Invalid(
                "invalid thermal bath conductance or temperature",
            ));
        }
        Ok(Self {
            conductance,
            temperature,
        })
    }

    pub const fn conductance(self) -> ThermalConductance {
        self.conductance
    }
    pub const fn temperature(self) -> Temperature {
        self.temperature
    }

    pub(super) fn increment_density(self, step: f64, temperature: f64, volume: f64) -> f64 {
        if self.conductance.value() == 0.0 {
            0.0
        } else {
            (step * self.conductance.value() / volume) * (self.temperature.value() - temperature)
        }
    }
}

/// Eliminate the linear heat equation inside the original caloric root.
/// Every trial owns its fixed forcing; rejected stages leave no heat history.
#[derive(Clone, Copy, Default)]
pub(super) struct HeatStage<'a> {
    fixed: Sum,
    capacity_density_j_m3_k: f64,
    bath_temperature: f64,
    #[cfg_attr(not(feature = "material"), allow(dead_code))]
    exchange: Option<HeatExchange<'a>>,
    #[cfg_attr(not(feature = "material"), allow(dead_code))]
    volume: f64,
}

impl<'a> HeatStage<'a> {
    pub fn held(heat: Sum) -> Self {
        Self {
            fixed: heat,
            ..Self::default()
        }
    }

    pub fn implicit(
        mut heat: Sum,
        exchange: Option<HeatExchange<'a>>,
        volume: f64,
        diagonal_step: f64,
        forcing: f64,
    ) -> Result<Self, IsochoricKineticsError> {
        heat.add(forcing);
        let (capacity, temperature) = exchange.map_or((0.0, 0.0), |exchange| {
            let temperature = match exchange {
                HeatExchange::Bath(bath) => bath.temperature.value(),
                #[cfg(feature = "material")]
                HeatExchange::Body(_) => 0.0,
            };
            (diagonal_step * exchange.conductance() / volume, temperature)
        });
        if !heat.value().is_finite() || !capacity.is_finite() || capacity < 0.0 {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable thermal stage forcing or conductance",
            ));
        }
        Ok(Self {
            fixed: heat,
            capacity_density_j_m3_k: capacity,
            bath_temperature: temperature,
            exchange,
            volume,
        })
    }

    pub fn heat(self, temperature: f64) -> Result<(Sum, f64), IsochoricKineticsError> {
        #[cfg(feature = "material")]
        if let Some(HeatExchange::Body(body)) = self.exchange {
            return body.stage_heat(
                self.fixed,
                self.volume,
                self.capacity_density_j_m3_k,
                temperature,
            );
        }
        let mut heat = self.fixed;
        if self.capacity_density_j_m3_k != 0.0 {
            heat.add(self.capacity_density_j_m3_k * (self.bath_temperature - temperature));
        }
        Ok((heat, self.capacity_density_j_m3_k))
    }

    pub fn temperature_interval(
        self,
        lower: f64,
        upper: f64,
    ) -> Result<(f64, f64), IsochoricKineticsError> {
        #[cfg(feature = "material")]
        if let Some(HeatExchange::Body(body)) = self.exchange {
            return body.stage_temperature_interval(
                self.fixed,
                self.volume,
                self.capacity_density_j_m3_k,
                lower,
                upper,
            );
        }
        Ok((lower, upper))
    }
}
