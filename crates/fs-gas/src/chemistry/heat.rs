//! Finite conductance to an explicitly prescribed, constant-temperature bath.
//!
//! The homogeneous parcel exchanges `Qdot = G (Tb - T)` with this reservoir
//! while its original reaction laws evolve. The bath has no finite heat
//! capacity or evolving state; its entropy change is `-Q/Tb`. Conductance is
//! supplied by the caller, not inferred from a material name or geometry.

use super::*;

/// Total parcel-to-reservoir thermal conductance, W/K (not conductivity).
pub type ThermalConductance = fs_qty::Qty<2, 1, -3, -1, 0, 0>;

/// Passive, finite thermal exchange with one held physical reservoir.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeatBath {
    conductance: ThermalConductance,
    temperature: Temperature,
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
pub(super) struct HeatStage {
    fixed: Sum,
    pub capacity_density_j_m3_k: f64,
    bath_temperature: f64,
}

impl HeatStage {
    pub fn held(heat: Sum) -> Self {
        Self {
            fixed: heat,
            ..Self::default()
        }
    }

    pub fn implicit(
        mut heat: Sum,
        bath: Option<HeatBath>,
        volume: f64,
        diagonal_step: f64,
        forcing: f64,
    ) -> Result<Self, IsochoricKineticsError> {
        heat.add(forcing);
        let (capacity, temperature) = bath.map_or((0.0, 0.0), |bath| {
            (
                diagonal_step * bath.conductance.value() / volume,
                bath.temperature.value(),
            )
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
        })
    }

    pub fn heat(self, temperature: f64) -> Sum {
        let mut heat = self.fixed;
        if self.capacity_density_j_m3_k != 0.0 {
            heat.add(self.capacity_density_j_m3_k * (self.bath_temperature - temperature));
        }
        heat
    }
}
