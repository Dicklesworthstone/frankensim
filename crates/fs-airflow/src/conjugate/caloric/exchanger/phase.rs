//! Finite-UA condensation and boiling on an original IF97 saturation plateau.
//!
//! One finite water/steam mass flow exchanges heat with a finite sensible
//! stream. The original Cp owner integrates the latter against the original
//! saturation temperature; the original water p,h inverse resolves quality.
//! Both water endpoints must remain on that same plateau. Exhausting latent
//! throughput refuses the entire result, rather than treating the phase stream
//! as an unlimited reservoir or extending it into an unmodelled sensible zone.
//!
//! Prescribed UA is independent of quality in this model. Co-current and
//! counterflow arrangements coincide while the entire water path is isothermal.
//! There is no wall storage, pressure loss, mass mixing, or boiling correlation.

use super::{ExchangerError, ExchangerStream, endpoint, nonnegative, poll, positive};
use fs_exec::Cx;
use fs_material::phase::water::flow::WaterFlowError;
use fs_material::phase::water::flow::across::{AcrossWaterFlowError, WaterFlowState};
use fs_math::dd::Dd;

/// Endpoint budgets, separate from the original source solvers' own controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterPhaseExchangerConfig {
    pub absolute_heat_tolerance_w: f64,
    pub relative_heat_tolerance: f64,
    /// Bounds each original state's h-u-p_requested*v defect. Dense IF97
    /// pressure inverses retain their actual source pressure, not a projection.
    pub maximum_pressure_closure_defect_j_kg: f64,
    /// Permitted negative RAW original-source entropy production. The local
    /// thermal seam adjustment is reported separately and never hides this gate.
    pub source_entropy_tolerance_w_per_k: f64,
}

impl Default for WaterPhaseExchangerConfig {
    fn default() -> Self {
        Self {
            absolute_heat_tolerance_w: 1e-7,
            relative_heat_tolerance: 1e-10,
            maximum_pressure_closure_defect_j_kg: 1e-4,
            source_entropy_tolerance_w_per_k: 0.0,
        }
    }
}

/// One reciprocal heat transfer and independently audited original endpoints.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterPhaseExchangerState {
    pub sensible_outlet_temperature_k: f64,
    pub water_outlet: WaterFlowState,
    /// Positive for condensation that heats the sensible stream; negative for
    /// boiling driven by that stream. The sole exchanged heat rate.
    pub heat_rate_water_to_sensible_w: f64,
    /// Relative heat coordinates survive rounding of absolute T and h.
    pub sensible_specific_enthalpy_change_j_kg: f64,
    pub water_specific_enthalpy_change_j_kg: f64,
    pub water_internal_energy_change_rate_w: f64,
    pub water_pressure_flow_work_change_rate_w: f64,
    pub sensible_endpoint_energy_residual_w: f64,
    pub water_endpoint_energy_residual_w: f64,
    /// m_water*(delta_u + p*delta_v) + Q, with flow work included once.
    pub water_first_law_residual_w: f64,
    pub sensible_endpoint_roundoff_allowance_w: f64,
    pub water_endpoint_roundoff_allowance_w: f64,
    pub water_pressure_closure_defects_j_kg: [f64; 2],
    /// Additional allowance for the independent first-law gate, exactly the
    /// mass rate times the sum of the two admitted absolute closure defects.
    pub water_pressure_closure_allowance_w: f64,
    /// Original Cp/T integral with its stated endpoint rounding correction.
    pub sensible_entropy_change_rate_w_per_k: f64,
    /// Retained heat times the ORIGINAL plateau's delta_s/delta_h. This is the
    /// source-linear entropy change before absolute endpoint rounding.
    pub water_source_entropy_change_rate_w_per_k: f64,
    pub source_entropy_generation_w_per_k: f64,
    /// Local fixed-pressure thermal primitive adjustment: delta_x times the
    /// original saturation seam. It is not a globally shifted source entropy.
    pub water_thermal_entropy_correction_rate_w_per_k: f64,
    /// Sensible entropy minus Q/T_sat, audited independently of raw source s.
    pub thermal_entropy_generation_w_per_k: f64,
    pub entropy_roundoff_allowance_w_per_k: f64,
    pub enthalpy_iterations: usize,
    pub dense_source_evaluations: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub enum WaterPhaseExchangerError {
    Caloric(ExchangerError),
    Water(AcrossWaterFlowError),
    OutsideWaterPlateau,
    LatentThroughputExceeded,
    PressureClosure { defects_j_kg: [f64; 2] },
    SourceEntropyDecrease { generation_w_per_k: f64 },
    Unresolved(&'static str),
    Interrupted,
}

impl From<ExchangerError> for WaterPhaseExchangerError {
    fn from(error: ExchangerError) -> Self {
        match error {
            ExchangerError::Interrupted => Self::Interrupted,
            other => Self::Caloric(other),
        }
    }
}
impl From<AcrossWaterFlowError> for WaterPhaseExchangerError {
    fn from(error: AcrossWaterFlowError) -> Self {
        match error {
            AcrossWaterFlowError::Flow(WaterFlowError::Cancelled) => Self::Interrupted,
            other => Self::Water(other),
        }
    }
}
impl core::fmt::Display for WaterPhaseExchangerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original water plateau exchanger: {self:?}")
    }
}
impl std::error::Error for WaterPhaseExchangerError {}

/// Immutable original-source condenser/boiler with prescribed finite UA.
#[derive(Clone, Debug, PartialEq)]
pub struct WaterPhaseHeatExchanger {
    sensible: ExchangerStream,
    water_inlet: WaterFlowState,
    water_mass_flow_kg_per_s: f64,
    conductance_w_per_k: f64,
}

impl WaterPhaseHeatExchanger {
    pub fn new(
        sensible: ExchangerStream,
        water_inlet: WaterFlowState,
        water_mass_flow_kg_per_s: f64,
        conductance_w_per_k: f64,
    ) -> Result<Self, WaterPhaseExchangerError> {
        positive(water_mass_flow_kg_per_s, "water mass flow")?;
        nonnegative(conductance_w_per_k, "conductance")?;
        if water_inlet.vapor_quality().is_none() {
            return Err(WaterPhaseExchangerError::OutsideWaterPlateau);
        }
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
    pub const fn water_inlet(&self) -> WaterFlowState {
        self.water_inlet
    }
    pub const fn water_mass_flow_kg_per_s(&self) -> f64 {
        self.water_mass_flow_kg_per_s
    }
    pub const fn conductance_w_per_k(&self) -> f64 {
        self.conductance_w_per_k
    }

    /// Integrate the original nonlinear Cp exactly as its native exchange law,
    /// resolve original p,h quality, then audit h and u+p*v independently.
    /// All admission, source-domain and entropy gates precede publication.
    pub fn solve(
        &self,
        cx: &Cx<'_>,
        config: WaterPhaseExchangerConfig,
    ) -> Result<WaterPhaseExchangerState, WaterPhaseExchangerError> {
        poll(cx)?;
        config.validate()?;
        let inlet = self.water_inlet;
        let chart = inlet.law();
        let liquid = chart.state_at_quality(0.0, || cx.checkpoint().is_err())?;
        let vapor = chart.state_at_quality(1.0, || cx.checkpoint().is_err())?;
        let latent = Dd::from_f64(vapor.specific_enthalpy_j_kg())
            - Dd::from_f64(liquid.specific_enthalpy_j_kg());
        positive(latent.to_f64(), "original latent enthalpy")?;
        let saturation_k = inlet.temperature_k();
        let exposure = self.conductance_w_per_k / self.sensible.mass_flow_kg_per_s;
        nonnegative(exposure, "conductance per sensible mass rate")?;
        if self.conductance_w_per_k != 0.0 && exposure == 0.0 {
            return Err(WaterPhaseExchangerError::Unresolved("nonzero exposure"));
        }
        let exchanged = self
            .sensible
            .caloric
            .linear_heat_exchange(self.sensible.inlet_temperature_k, saturation_k, exposure)
            .map_err(ExchangerError::from)?;
        poll(cx)?;
        let heat = self.sensible.mass_flow_kg_per_s * exchanged.specific_enthalpy_change_j_kg;
        let water_heat = -heat / self.water_mass_flow_kg_per_s;
        let expected_nonzero =
            self.conductance_w_per_k > 0.0 && self.sensible.inlet_temperature_k != saturation_k;
        if !heat.is_finite()
            || !water_heat.is_finite()
            || (expected_nonzero && (heat == 0.0 || water_heat == 0.0))
        {
            return Err(WaterPhaseExchangerError::Unresolved("retained heat"));
        }
        // Compare before rounding absolute h: even a tiny nonzero request
        // beyond saturated liquid/vapor cannot borrow nonexistent phase mass.
        let target = Dd::from_f64(inlet.specific_enthalpy_j_kg()) + Dd::from_f64(water_heat);
        if target.lt(Dd::from_f64(liquid.specific_enthalpy_j_kg()))
            || Dd::from_f64(vapor.specific_enthalpy_j_kg()).lt(target)
        {
            return Err(WaterPhaseExchangerError::LatentThroughputExceeded);
        }
        let (outlet, enthalpy_iterations, dense_source_evaluations) = if heat == 0.0 {
            (inlet, 0, 0)
        } else {
            chart.state_at_specific_enthalpy(target.to_f64(), || cx.checkpoint().is_err())?
        };
        if outlet.vapor_quality().is_none() || outlet.temperature_k() != saturation_k {
            return Err(WaterPhaseExchangerError::OutsideWaterPlateau);
        }
        let tolerance = positive(
            config.absolute_heat_tolerance_w + config.relative_heat_tolerance * heat.abs(),
            "combined heat tolerance",
        )?;
        let sensible = endpoint(&self.sensible, exchanged.temperature_k, heat, tolerance)?;
        let defects = [
            inlet.pressure_closure_defect_j_kg(),
            outlet.pressure_closure_defect_j_kg(),
        ];
        if defects
            .iter()
            .any(|d| !d.is_finite() || d.abs() > config.maximum_pressure_closure_defect_j_kg)
        {
            return Err(WaterPhaseExchangerError::PressureClosure {
                defects_j_kg: defects,
            });
        }
        let mass = Dd::from_f64(self.water_mass_flow_kg_per_s);
        let pressure = Dd::from_f64(inlet.pressure_pa());
        let enthalpy = mass
            * (Dd::from_f64(outlet.specific_enthalpy_j_kg())
                - Dd::from_f64(inlet.specific_enthalpy_j_kg()));
        let internal = mass
            * (Dd::from_f64(outlet.specific_internal_energy_j_kg())
                - Dd::from_f64(inlet.specific_internal_energy_j_kg()));
        let flow_work = mass
            * pressure
            * (Dd::from_f64(outlet.specific_volume_m3_kg())
                - Dd::from_f64(inlet.specific_volume_m3_kg()));
        let h_residual = (enthalpy + Dd::from_f64(heat)).to_f64();
        let first_law_residual = (internal + flow_work + Dd::from_f64(heat)).to_f64();
        let roundoff = 128.0
            * f64::EPSILON
            * (heat.abs()
                + self.water_mass_flow_kg_per_s
                    * (inlet.specific_enthalpy_j_kg().abs()
                        + outlet.specific_enthalpy_j_kg().abs()
                        + inlet.specific_internal_energy_j_kg().abs()
                        + outlet.specific_internal_energy_j_kg().abs()
                        + inlet.pressure_pa()
                            * (inlet.specific_volume_m3_kg() + outlet.specific_volume_m3_kg())));
        let closure_allowance =
            self.water_mass_flow_kg_per_s * (defects[0].abs() + defects[1].abs());
        if ![h_residual, first_law_residual, roundoff, closure_allowance]
            .into_iter()
            .all(f64::is_finite)
            || h_residual.abs() > tolerance + roundoff
            || first_law_residual.abs() > tolerance + roundoff + closure_allowance
        {
            return Err(WaterPhaseExchangerError::Unresolved(
                "water endpoint energy audit",
            ));
        }
        let source_entropy_span = Dd::from_f64(vapor.specific_entropy_j_kg_k())
            - Dd::from_f64(liquid.specific_entropy_j_kg_k());
        let water_entropy = -Dd::from_f64(heat) * source_entropy_span / latent;
        let thermal_water_entropy = -Dd::from_f64(heat) / Dd::from_f64(saturation_k);
        let seam = thermal_water_entropy - water_entropy;
        let source_entropy = (sensible.entropy + water_entropy).to_f64();
        let thermal_entropy = (sensible.entropy + thermal_water_entropy).to_f64();
        let entropy_allowance = sensible.entropy_allowance
            + 128.0
                * f64::EPSILON
                * heat.abs()
                * ((vapor.specific_entropy_j_kg_k().abs()
                    + liquid.specific_entropy_j_kg_k().abs())
                    / latent.to_f64()
                    + 1.0 / saturation_k);
        if ![
            source_entropy,
            thermal_entropy,
            entropy_allowance,
            seam.to_f64(),
        ]
        .into_iter()
        .all(f64::is_finite)
            || thermal_entropy < -entropy_allowance
        {
            return Err(WaterPhaseExchangerError::Unresolved(
                "local thermal entropy audit",
            ));
        }
        if source_entropy < -config.source_entropy_tolerance_w_per_k - entropy_allowance {
            return Err(WaterPhaseExchangerError::SourceEntropyDecrease {
                generation_w_per_k: source_entropy,
            });
        }
        poll(cx)?;
        Ok(WaterPhaseExchangerState {
            sensible_outlet_temperature_k: exchanged.temperature_k,
            water_outlet: outlet,
            heat_rate_water_to_sensible_w: heat,
            sensible_specific_enthalpy_change_j_kg: sensible.heat_per_mass,
            water_specific_enthalpy_change_j_kg: water_heat,
            water_internal_energy_change_rate_w: internal.to_f64(),
            water_pressure_flow_work_change_rate_w: flow_work.to_f64(),
            sensible_endpoint_energy_residual_w: sensible.energy_residual,
            water_endpoint_energy_residual_w: h_residual,
            water_first_law_residual_w: first_law_residual,
            sensible_endpoint_roundoff_allowance_w: sensible.energy_roundoff,
            water_endpoint_roundoff_allowance_w: roundoff,
            water_pressure_closure_defects_j_kg: defects,
            water_pressure_closure_allowance_w: closure_allowance,
            sensible_entropy_change_rate_w_per_k: sensible.entropy.to_f64(),
            water_source_entropy_change_rate_w_per_k: water_entropy.to_f64(),
            source_entropy_generation_w_per_k: source_entropy,
            water_thermal_entropy_correction_rate_w_per_k: seam.to_f64(),
            thermal_entropy_generation_w_per_k: thermal_entropy,
            entropy_roundoff_allowance_w_per_k: entropy_allowance,
            enthalpy_iterations,
            dense_source_evaluations,
        })
    }
}

impl WaterPhaseExchangerConfig {
    fn validate(self) -> Result<(), WaterPhaseExchangerError> {
        positive(self.absolute_heat_tolerance_w, "absolute heat tolerance")?;
        nonnegative(self.relative_heat_tolerance, "relative heat tolerance")?;
        nonnegative(
            self.maximum_pressure_closure_defect_j_kg,
            "pressure closure budget",
        )?;
        nonnegative(
            self.source_entropy_tolerance_w_per_k,
            "source entropy tolerance",
        )?;
        if self.relative_heat_tolerance >= 1.0 {
            return Err(ExchangerError::Invalid("relative heat tolerance").into());
        }
        Ok(())
    }
}
