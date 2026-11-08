//! Continuous co-current water/steam heat exchange through sensible and latent zones.
//!
//! With normalized cumulative conductance x in [0,1], the steady equations are
//! m_dot_water dh_water/dx = -UA*(T_water-T_sensible), with the opposite heat
//! on the sensible stream. The existing conservative thermal owner solves
//! exactly these equations using auxiliary masses m_dot*tau and auxiliary time
//! t=x*tau. Here tau=1 second; it is a coordinate scale, not fluid holdup or
//! physical transit time. No exchanger storage, pressure loss or mixing is added.
//!
//! Both original caloric laws supply T(h), including the true water plateau.
//! Native adaptive SSPRK2 owns all stages, enthalpy/temperature error estimates,
//! cancellation and work limits. Its local estimate is not a global error bound;
//! phase corners require refinement. A pairwise inverse-slope limit also keeps
//! every stage on the passive side of co-current equilibrium. Native material
//! calorimetry then independently audits h, u+p*v and both entropy conventions.

use super::phase::{WaterPhaseExchangerConfig, WaterPhaseExchangerError};
use super::{ExchangerError, ExchangerStream, endpoint, nonnegative, poll, positive};
use fs_exec::Cx;
use fs_material::phase::water::flow::across::WaterFlowState;
use fs_material::phase::water::{WaterError, WaterState};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::advection::Advected;
use fs_thermal::enthalpy::advection::explicit::{
    AdaptiveSsprk2Controls, AdaptiveSsprk2Error, AdaptiveSsprk2Report, Ssprk2Controls,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

/// Native accuracy/work controls expressed over normalized cumulative UA.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoCurrentWaterExchangerConfig {
    pub endpoint: WaterPhaseExchangerConfig,
    pub enthalpy_absolute_tolerance_j_kg: f64,
    pub temperature_absolute_tolerance_k: Option<f64>,
    /// Relative to the local enthalpy change, independent of enthalpy datum.
    pub relative_tolerance: f64,
    pub maximum_conductance_fraction: f64,
    pub initial_conductance_fraction: f64,
    pub minimum_conductance_fraction: f64,
    pub maximum_trials: usize,
    pub maximum_stage_evaluations: usize,
}

impl Default for CoCurrentWaterExchangerConfig {
    fn default() -> Self {
        Self {
            endpoint: WaterPhaseExchangerConfig {
                absolute_heat_tolerance_w: 1e-5,
                ..WaterPhaseExchangerConfig::default()
            },
            enthalpy_absolute_tolerance_j_kg: 1e-3,
            temperature_absolute_tolerance_k: Some(1e-6),
            relative_tolerance: 1e-8,
            maximum_conductance_fraction: 0.1,
            initial_conductance_fraction: 0.05,
            minimum_conductance_fraction: 1e-10,
            maximum_trials: 16_000,
            maximum_stage_evaluations: 80_000,
        }
    }
}

/// Original source endpoint calorimetry shared by continuous water exchangers.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterExchangerState {
    pub sensible_outlet_temperature_k: f64,
    pub water_outlet: WaterState,
    /// One shared Q, positive from water to the sensible stream.
    pub heat_rate_water_to_sensible_w: f64,
    pub sensible_specific_enthalpy_change_j_kg: f64,
    pub water_specific_enthalpy_change_j_kg: f64,
    pub water_internal_energy_change_rate_w: f64,
    pub water_pressure_flow_work_change_rate_w: f64,
    /// Native relative source calorimetry minus retained heat, [h, u+p*v].
    pub water_source_energy_residuals_w: [f64; 2],
    /// Independently evaluated represented endpoints, [h, u+p*v].
    pub water_endpoint_energy_residuals_w: [f64; 2],
    pub sensible_endpoint_energy_residual_w: f64,
    pub sensible_endpoint_roundoff_allowance_w: f64,
    pub water_pressure_closure_defects_j_kg: [f64; 2],
    pub sensible_entropy_change_rate_w_per_k: f64,
    pub water_source_entropy_change_rate_w_per_k: f64,
    pub source_entropy_generation_w_per_k: f64,
    pub water_thermal_entropy_correction_rate_w_per_k: f64,
    pub thermal_entropy_generation_w_per_k: f64,
    /// Explicit endpoint-inverse/rounding allowance, not source uncertainty.
    pub entropy_numerical_allowance_w_per_k: f64,
    /// Native report on the declared tau=1 s auxiliary coordinate. Energy
    /// defects divided by tau are watt defects. This is not exchanger storage.
    /// None when no native SSPRK2 trajectory is needed, including duty rating
    /// and countercurrent quadrature, zero UA or equal inlet temperatures.
    pub integration: Option<AdaptiveSsprk2Report>,
}

pub use WaterExchangerState as CoCurrentWaterExchangerState;

#[derive(Clone, Debug, PartialEq)]
pub enum CoCurrentWaterExchangerError {
    Caloric(ExchangerError),
    Water(WaterError),
    Thermal(AdaptiveSsprk2Error),
    Endpoint(WaterPhaseExchangerError),
    Unresolved(&'static str),
    Interrupted,
}
impl From<ExchangerError> for CoCurrentWaterExchangerError {
    fn from(error: ExchangerError) -> Self {
        match error {
            ExchangerError::Interrupted => Self::Interrupted,
            other => Self::Caloric(other),
        }
    }
}
impl From<WaterError> for CoCurrentWaterExchangerError {
    fn from(error: WaterError) -> Self {
        match error {
            WaterError::Cancelled => Self::Interrupted,
            other => Self::Water(other),
        }
    }
}
impl From<AdaptiveSsprk2Error> for CoCurrentWaterExchangerError {
    fn from(error: AdaptiveSsprk2Error) -> Self {
        match error {
            AdaptiveSsprk2Error::Deposition(DepositionError::Cancelled) => Self::Interrupted,
            other => Self::Thermal(other),
        }
    }
}
impl From<DepositionError> for CoCurrentWaterExchangerError {
    fn from(error: DepositionError) -> Self {
        AdaptiveSsprk2Error::from(error).into()
    }
}
impl From<WaterPhaseExchangerError> for CoCurrentWaterExchangerError {
    fn from(error: WaterPhaseExchangerError) -> Self {
        match error {
            WaterPhaseExchangerError::Interrupted => Self::Interrupted,
            other => Self::Endpoint(other),
        }
    }
}
impl core::fmt::Display for CoCurrentWaterExchangerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original co-current water exchanger: {self:?}")
    }
}
impl std::error::Error for CoCurrentWaterExchangerError {}

/// Steady co-current plug flow with original Gibbs water and one Cp source.
/// Pressure and composition remain fixed independently on each side.
#[derive(Clone, Debug, PartialEq)]
pub struct CoCurrentWaterHeatExchanger {
    sensible: ExchangerStream,
    water_inlet: WaterState,
    water_mass_flow_kg_per_s: f64,
    conductance_w_per_k: f64,
}

impl CoCurrentWaterHeatExchanger {
    pub fn new(
        sensible: ExchangerStream,
        water_inlet: WaterState,
        water_mass_flow_kg_per_s: f64,
        conductance_w_per_k: f64,
    ) -> Result<Self, CoCurrentWaterExchangerError> {
        positive(water_mass_flow_kg_per_s, "water mass flow")?;
        nonnegative(conductance_w_per_k, "conductance")?;
        Ok(Self {
            sensible,
            water_inlet,
            water_mass_flow_kg_per_s,
            conductance_w_per_k,
        })
    }

    pub const fn sensible_stream(&self) -> &ExchangerStream {
        &self.sensible
    }
    pub const fn water_inlet(&self) -> WaterState {
        self.water_inlet
    }
    pub const fn water_mass_flow_kg_per_s(&self) -> f64 {
        self.water_mass_flow_kg_per_s
    }
    pub const fn conductance_w_per_k(&self) -> f64 {
        self.conductance_w_per_k
    }

    pub fn solve(
        &self,
        cx: &Cx<'_>,
        config: CoCurrentWaterExchangerConfig,
    ) -> Result<CoCurrentWaterExchangerState, CoCurrentWaterExchangerError> {
        poll(cx)?;
        config.endpoint.validate()?;
        let law = self.water_inlet.law();
        let sensible_h = self
            .sensible
            .caloric
            .specific_enthalpy_j_kg(self.sensible.inlet_temperature_k)
            .map_err(ExchangerError::from)?;
        let pair_rate = self.conductance_w_per_k
            * (law.maximum_temperature_slope() / self.water_mass_flow_kg_per_s
                + self.sensible.caloric.maximum_temperature_slope()
                    / self.sensible.mass_flow_kg_per_s);
        nonnegative(pair_rate, "pair inverse-slope rate")?;
        // In auxiliary seconds, tau=1 s. A complete Euler stage cannot cross
        // equilibrium when dt * UA * (L_water/m_water + L_sensible/m_sensible) < 1.
        let maximum_fraction = if pair_rate == 0.0 {
            1.0
        } else {
            (0.8 / pair_rate).min(1.0)
        };
        let controls = config.controls(maximum_fraction)?;
        if self.conductance_w_per_k == 0.0
            || self.water_inlet.temperature_k() == self.sensible.inlet_temperature_k
        {
            return self.finish(
                cx,
                0.0,
                self.sensible.inlet_temperature_k,
                None,
                config.endpoint,
            );
        }
        let checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
        let births = [
            cell(
                0,
                self.sensible.mass_flow_kg_per_s,
                sensible_h,
                self.sensible.inlet_temperature_k,
            ),
            cell(
                1,
                self.water_mass_flow_kg_per_s,
                self.water_inlet.specific_enthalpy_j_kg(),
                self.water_inlet.temperature_k(),
            ),
        ];
        let schedule = vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(births[0].clone()),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(births[1].clone()),
            },
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 0,
                    a: CellId(0),
                    b: CellId(1),
                    area_m2: 1.0,
                    conductance_w_k: self.conductance_w_per_k,
                }),
            },
        ];
        let thermal = DepositionSimulation::new(
            schedule,
            &[
                (MaterialId(0), self.sensible.caloric.as_ref()),
                (MaterialId(1), law),
            ],
            DepositionBudget {
                max_cells: 2,
                max_links: 1,
                max_events: 3,
                max_substeps: config.maximum_stage_evaluations,
            },
        )?;
        let (thermal, _) = thermal.advance_to(0.0, controls.step.step, checkpoint)?;
        let owner = Advected::new((), vec![], 0)?;
        let (outlet, _, report) =
            owner.advance_ssprk2_adaptive_to(&thermal, 1.0, controls, checkpoint)?;
        let water_h = outlet.cells()[&CellId(1)].specific_h_j_kg;
        let heat = -(Dd::from_f64(self.water_mass_flow_kg_per_s)
            * (Dd::from_f64(water_h) - Dd::from_f64(self.water_inlet.specific_enthalpy_j_kg())))
        .to_f64();
        if !heat.is_finite() || heat == 0.0 {
            return Err(CoCurrentWaterExchangerError::Unresolved(
                "nonzero retained heat",
            ));
        }
        let sensible_outlet = outlet.temperature_k(CellId(0))?;
        let water_outlet = outlet.temperature_k(CellId(1))?;
        let initial_gap = self.water_inlet.temperature_k() - self.sensible.inlet_temperature_k;
        if (water_outlet - sensible_outlet) * initial_gap < 0.0 {
            return Err(CoCurrentWaterExchangerError::Unresolved(
                "co-current equilibrium crossing",
            ));
        }
        self.finish(cx, heat, sensible_outlet, Some(report), config.endpoint)
    }

    pub(super) fn finish(
        &self,
        cx: &Cx<'_>,
        heat: f64,
        sensible_outlet: f64,
        integration: Option<AdaptiveSsprk2Report>,
        config: WaterPhaseExchangerConfig,
    ) -> Result<CoCurrentWaterExchangerState, CoCurrentWaterExchangerError> {
        let tolerance = positive(
            config.absolute_heat_tolerance_w + config.relative_heat_tolerance * heat.abs(),
            "combined heat tolerance",
        )?;
        let sensible = endpoint(&self.sensible, sensible_outlet, heat, tolerance)?;
        // The original material performs its own stable sensible/latent source
        // increments and independent h and u+p*v audits on the same tau=1 s map.
        let water = self.water_inlet.law().advance_heat(
            self.water_inlet,
            self.water_mass_flow_kg_per_s,
            -heat,
            tolerance,
            || cx.checkpoint().is_err(),
        )?;
        let defects = [self.water_inlet, water.state]
            .map(|s| WaterFlowState::from(s).pressure_closure_defect_j_kg());
        if defects
            .iter()
            .any(|d| !d.is_finite() || d.abs() > config.maximum_pressure_closure_defect_j_kg)
        {
            return Err(WaterPhaseExchangerError::PressureClosure {
                defects_j_kg: defects,
            }
            .into());
        }
        let source_entropy =
            (sensible.entropy + Dd::from_f64(water.source_entropy_change_j_k)).to_f64();
        let thermal_entropy =
            (sensible.entropy + Dd::from_f64(water.thermal_entropy_change_j_k)).to_f64();
        let minimum_water_temperature = self
            .water_inlet
            .temperature_k()
            .min(water.state.temperature_k());
        let entropy_allowance = sensible.entropy_allowance
            + water.enthalpy_residual_j.abs() / minimum_water_temperature
            + 128.0
                * f64::EPSILON
                * (water.source_entropy_change_j_k.abs()
                    + water.thermal_entropy_change_j_k.abs()
                    + heat.abs() / minimum_water_temperature);
        if ![source_entropy, thermal_entropy, entropy_allowance]
            .into_iter()
            .all(f64::is_finite)
            || thermal_entropy < -entropy_allowance
        {
            return Err(CoCurrentWaterExchangerError::Unresolved(
                "original thermal entropy audit",
            ));
        }
        if source_entropy < -config.source_entropy_tolerance_w_per_k - entropy_allowance {
            return Err(WaterPhaseExchangerError::SourceEntropyDecrease {
                generation_w_per_k: source_entropy,
            }
            .into());
        }
        poll(cx)?;
        Ok(CoCurrentWaterExchangerState {
            sensible_outlet_temperature_k: sensible_outlet,
            water_outlet: water.state,
            heat_rate_water_to_sensible_w: heat,
            sensible_specific_enthalpy_change_j_kg: sensible.heat_per_mass,
            water_specific_enthalpy_change_j_kg: -heat / self.water_mass_flow_kg_per_s,
            water_internal_energy_change_rate_w: water.internal_energy_change_j,
            water_pressure_flow_work_change_rate_w: water.boundary_work_j,
            water_source_energy_residuals_w: [water.enthalpy_residual_j, water.energy_residual_j],
            water_endpoint_energy_residuals_w: [
                water.represented_enthalpy_residual_j,
                water.represented_energy_residual_j,
            ],
            sensible_endpoint_energy_residual_w: sensible.energy_residual,
            sensible_endpoint_roundoff_allowance_w: sensible.energy_roundoff,
            water_pressure_closure_defects_j_kg: defects,
            sensible_entropy_change_rate_w_per_k: sensible.entropy.to_f64(),
            water_source_entropy_change_rate_w_per_k: water.source_entropy_change_j_k,
            source_entropy_generation_w_per_k: source_entropy,
            water_thermal_entropy_correction_rate_w_per_k: water.entropy_seam_correction_change_j_k,
            thermal_entropy_generation_w_per_k: thermal_entropy,
            entropy_numerical_allowance_w_per_k: entropy_allowance,
            integration,
        })
    }
}

impl CoCurrentWaterExchangerConfig {
    fn controls(
        self,
        passive_maximum: f64,
    ) -> Result<AdaptiveSsprk2Controls, CoCurrentWaterExchangerError> {
        for value in [
            self.maximum_conductance_fraction,
            self.initial_conductance_fraction,
            self.minimum_conductance_fraction,
        ] {
            positive(value, "conductance progress fraction")?;
            if value > 1.0 {
                return Err(
                    ExchangerError::Invalid("conductance progress fraction above one").into(),
                );
            }
        }
        if self.minimum_conductance_fraction > self.initial_conductance_fraction
            || self.initial_conductance_fraction > self.maximum_conductance_fraction
        {
            return Err(ExchangerError::Invalid("conductance progress fraction order").into());
        }
        positive(self.enthalpy_absolute_tolerance_j_kg, "enthalpy accuracy")?;
        if let Some(t) = self.temperature_absolute_tolerance_k {
            positive(t, "temperature accuracy")?;
        }
        nonnegative(self.relative_tolerance, "relative accuracy")?;
        if self.relative_tolerance >= 1.0
            || self.maximum_trials == 0
            || self.maximum_stage_evaluations == 0
        {
            return Err(ExchangerError::Invalid("accuracy or work limits").into());
        }
        let maximum = positive(
            self.maximum_conductance_fraction.min(passive_maximum),
            "passive conductance progress",
        )?;
        Ok(AdaptiveSsprk2Controls {
            step: Ssprk2Controls {
                step: StepConfig {
                    max_dt_s: maximum,
                    energy_tolerance_j: self.endpoint.absolute_heat_tolerance_w,
                },
                maximum_stage_evaluations: self.maximum_stage_evaluations,
            },
            initial_dt_s: self.initial_conductance_fraction.min(maximum),
            minimum_dt_s: self.minimum_conductance_fraction,
            maximum_trials: self.maximum_trials,
            enthalpy_absolute_tolerance_j_kg: self.enthalpy_absolute_tolerance_j_kg,
            temperature_absolute_tolerance_k: self.temperature_absolute_tolerance_k,
            relative_tolerance: self.relative_tolerance,
        })
    }
}

// Only mass and the prescribed total link conductance enter these equations.
// Unit auxiliary geometry satisfies native topology admission; it supplies no
// inferred channel dimensions, density, heat coefficient or external heat.
fn cell(id: u64, mass_kg: f64, h: f64, temperature_k: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(id),
        mass_kg,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [1.0, 0.0, 0.0],
            volume_m3: 1.0,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: temperature_k,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}
