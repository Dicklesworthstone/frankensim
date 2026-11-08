//! Conductive and gray radiative exchange with one constant-temperature bath.
//!
//! The homogeneous parcel exchanges `Qdot = G (Tb - T) + sigma epsilon A
//! (Tb^4 - T^4)` with this reservoir while its original reaction laws evolve.
//! The bath has no finite heat capacity or evolving state; its entropy change
//! is `-Q/Tb`. Conductance and effective radiating area are caller-supplied
//! closures, not inferred from a material name, gas composition or geometry.

use super::*;
use fs_qty::Area;

// Same SI convention as the native gray surface/gap owners in fs-thermal.
// W m^-2 K^-4, rounded from the post-2019 defining SI constants.
const STEFAN_BOLTZMANN: f64 = 5.670_374_419e-8;

/// Finite original NASA-mixture storage with independent pressure/volume constraint.
pub mod gas;

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
    effective_radiating_area: Area,
}

/// Shared thermal partners entering the same implicit gas stages.
#[derive(Clone, Copy)]
pub(super) enum HeatExchange<'a> {
    Bath(&'a HeatBath),
    Gas(&'a gas::FiniteGasHeatBody<'a>),
    #[cfg(feature = "material")]
    Body(material::FiniteHeatBody<'a>),
}

// Finite-body state parameters disappear from the default bath-only match.
#[cfg_attr(not(feature = "material"), allow(unused_variables))]
impl HeatExchange<'_> {
    fn conductance(self) -> f64 {
        match self {
            Self::Bath(bath) => bath.conductance.value(),
            Self::Gas(body) => body.conductance().value(),
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
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, IsochoricKineticsError> {
        match self {
            Self::Bath(bath) => bath.increment_density(step, temperature, volume),
            Self::Gas(body) => {
                body.increment_density(step, temperature, heat, volume, controls, cancelled)
            }
            #[cfg(feature = "material")]
            Self::Body(body) => body.increment_density(step, temperature, heat, volume),
        }
    }

    pub fn entropy_transfer(
        self,
        initial: Sum,
        final_heat: Sum,
        volume: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, IsochoricKineticsError> {
        match self {
            Self::Bath(bath) => {
                Ok(final_heat.difference(initial) * volume / bath.temperature.value())
            }
            Self::Gas(body) => {
                body.entropy_transfer(initial, final_heat, volume, controls, cancelled)
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
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, IsochoricKineticsError> {
        match self {
            Self::Bath(_) => Ok(0.0),
            Self::Gas(body) => body.additional_error(base, full, fine, volume, controls, cancelled),
            #[cfg(feature = "material")]
            Self::Body(body) => body.additional_error(base, full, fine, volume, controls),
        }
    }

    pub fn admit_heat(
        self,
        heat: Sum,
        volume: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), IsochoricKineticsError> {
        match self {
            Self::Bath(_) => Ok(()),
            Self::Gas(body) => body.admit_heat(heat.value() * volume, controls, cancelled),
            #[cfg(feature = "material")]
            Self::Body(body) => body.admit_heat(heat.value() * volume, controls.energy_tolerance_j),
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
            effective_radiating_area: Area::new(0.0),
        })
    }

    /// Add or replace gray thermal radiation to this same held reservoir.
    ///
    /// `area` is nonnegative and `effective_emissivity` is in [0,1]. The caller
    /// supplies the complete gray geometry/emissivity closure: only their
    /// product, the effective radiating area, is retained. Setting either to
    /// zero restores the original linear-bath arithmetic. In particular, this
    /// is not a species-derived participating-gas opacity or spectral model.
    /// The declared area is held during an advance, including pressure-driven
    /// gas expansion; moving geometry requires an explicit bath schedule.
    pub fn with_gray_radiation(
        mut self,
        area: Area,
        effective_emissivity: f64,
    ) -> Result<Self, IsochoricKineticsError> {
        if !area.value().is_finite()
            || area.value() < 0.0
            || !effective_emissivity.is_finite()
            || !(0.0..=1.0).contains(&effective_emissivity)
        {
            return Err(IsochoricKineticsError::Invalid(
                "invalid gray bath area or effective emissivity",
            ));
        }
        self.effective_radiating_area = Area::new(area.value() * effective_emissivity);
        Ok(self)
    }

    /// Nonradiative total conductance. Radiation has its own nonlinear slope.
    pub const fn conductance(self) -> ThermalConductance {
        self.conductance
    }
    pub const fn temperature(self) -> Temperature {
        self.temperature
    }

    /// Caller-declared `epsilon * A` [m²]; zero disables gray radiation.
    pub const fn effective_radiating_area(self) -> Area {
        self.effective_radiating_area
    }

    fn radiation_capacity_density(self, step: f64, volume: f64) -> f64 {
        (step * (STEFAN_BOLTZMANN * self.effective_radiating_area.value())) / volume
    }

    pub(super) fn increment_density(
        self,
        step: f64,
        temperature: f64,
        volume: f64,
    ) -> Result<f64, IsochoricKineticsError> {
        let linear = if self.conductance.value() == 0.0 {
            0.0
        } else {
            (step * self.conductance.value() / volume) * (self.temperature.value() - temperature)
        };
        if self.effective_radiating_area.value() == 0.0 {
            return Ok(linear);
        }
        let (radiation, _) = radiation_heat(
            self.radiation_capacity_density(step, volume),
            self.temperature.value(),
            temperature,
        )?;
        let increment = linear + radiation;
        if !increment.is_finite() {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable radiative bath heat increment",
            ));
        }
        Ok(increment)
    }
}

/// Return stage radiation and its positive caloric feedback `-dQ/dT`.
fn radiation_heat(
    coefficient: f64,
    bath_temperature: f64,
    temperature: f64,
) -> Result<(f64, f64), IsochoricKineticsError> {
    // Factoring the fourth-power difference preserves the represented Kelvin
    // separation near equilibrium; two separately rounded T^4 values do not.
    let secant = coefficient
        * (bath_temperature + temperature)
        * bath_temperature.mul_add(bath_temperature, temperature * temperature);
    let heat = secant * (bath_temperature - temperature);
    let feedback = ((coefficient * temperature) * temperature) * (4.0 * temperature);
    if !heat.is_finite() || !feedback.is_finite() || feedback < 0.0 {
        return Err(IsochoricKineticsError::Invalid(
            "unrepresentable radiative bath heat or tangent",
        ));
    }
    Ok((heat, feedback))
}

/// Eliminate the complete bath heat equation inside the original caloric root.
/// Every trial owns its fixed forcing; rejected stages leave no heat history.
#[derive(Clone, Copy, Default)]
pub(super) struct HeatStage<'a> {
    fixed: Sum,
    capacity_density_j_m3_k: f64,
    radiation_capacity_density_j_m3_k4: f64,
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
                HeatExchange::Gas(_) => 0.0,
                #[cfg(feature = "material")]
                HeatExchange::Body(_) => 0.0,
            };
            (diagonal_step * exchange.conductance() / volume, temperature)
        });
        let radiation_capacity = match exchange {
            Some(HeatExchange::Bath(bath)) if bath.effective_radiating_area.value() != 0.0 => {
                bath.radiation_capacity_density(diagonal_step, volume)
            }
            Some(HeatExchange::Gas(body)) if body.effective_radiating_area().value() != 0.0 => {
                body.radiation_capacity_density(diagonal_step, volume)
            }
            #[cfg(feature = "material")]
            Some(HeatExchange::Body(body)) if body.effective_radiating_area().value() != 0.0 => {
                body.radiation_capacity_density(diagonal_step, volume)
            }
            _ => 0.0,
        };
        if !heat.value().is_finite()
            || !capacity.is_finite()
            || capacity < 0.0
            || !radiation_capacity.is_finite()
            || radiation_capacity < 0.0
        {
            return Err(IsochoricKineticsError::Invalid(
                "unrepresentable thermal stage forcing or conductance",
            ));
        }
        Ok(Self {
            fixed: heat,
            capacity_density_j_m3_k: capacity,
            radiation_capacity_density_j_m3_k4: radiation_capacity,
            bath_temperature: temperature,
            exchange,
            volume,
        })
    }

    #[cfg_attr(not(feature = "material"), allow(unused_variables))]
    pub fn heat(
        self,
        temperature: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Sum, f64), IsochoricKineticsError> {
        if let Some(HeatExchange::Gas(body)) = self.exchange {
            return body.stage_heat(
                self.fixed,
                self.volume,
                self.capacity_density_j_m3_k,
                self.radiation_capacity_density_j_m3_k4,
                temperature,
                controls,
                cancelled,
            );
        }
        #[cfg(feature = "material")]
        if let Some(HeatExchange::Body(body)) = self.exchange {
            if self.radiation_capacity_density_j_m3_k4 != 0.0 {
                return body.radiative_stage_heat(
                    self.fixed,
                    self.volume,
                    self.capacity_density_j_m3_k,
                    self.radiation_capacity_density_j_m3_k4,
                    temperature,
                    controls,
                    cancelled,
                );
            }
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
        let mut feedback = self.capacity_density_j_m3_k;
        if self.radiation_capacity_density_j_m3_k4 != 0.0 {
            let (radiation, radiation_feedback) = radiation_heat(
                self.radiation_capacity_density_j_m3_k4,
                self.bath_temperature,
                temperature,
            )?;
            heat.add(radiation);
            feedback += radiation_feedback;
            if !heat.value().is_finite() || !feedback.is_finite() {
                return Err(IsochoricKineticsError::Invalid(
                    "unrepresentable combined bath heat or tangent",
                ));
            }
        }
        Ok((heat, feedback))
    }

    #[cfg_attr(not(feature = "material"), allow(unused_variables))]
    pub fn temperature_interval(
        self,
        lower: f64,
        upper: f64,
        controls: IsochoricKineticsSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), IsochoricKineticsError> {
        if let Some(HeatExchange::Gas(body)) = self.exchange {
            return body.stage_temperature_interval(
                self.fixed,
                self.volume,
                self.capacity_density_j_m3_k,
                self.radiation_capacity_density_j_m3_k4,
                [lower, upper],
                controls,
                cancelled,
            );
        }
        #[cfg(feature = "material")]
        if let Some(HeatExchange::Body(body)) = self.exchange {
            if self.radiation_capacity_density_j_m3_k4 != 0.0 {
                return body.radiative_stage_temperature_interval(
                    self.fixed,
                    self.volume,
                    self.capacity_density_j_m3_k,
                    self.radiation_capacity_density_j_m3_k4,
                    [lower, upper],
                    controls,
                    cancelled,
                );
            }
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
