//! Temperature-dependent transport and sensible enthalpy on the same flow path.
//!
//! Each isothermal-wall region uses one explicitly selected representative
//! property temperature. A bounded scalar solve closes its outlet against the
//! coefficient evaluated from that outlet. This is a regional approximation,
//! not a resolved variable-property boundary layer; refinement remains physical
//! model refinement. The existing caloric interval marcher owns energy transport.
//! Mass flow and pressure are prescribed, not solved or silently held as volume.
//!
//! Bulk/film property conventions are explicit, as in the heat-transfer models
//! described at <https://doc.comsol.com/6.4/doc/com.comsol.help.pipe/pipe_ug_heattransfer.06.17.html>.
//! The arithmetic endpoint mean here is a declared regional closure, not an
//! exact spatial mean or a new correlation-accuracy claim.

pub mod solid;
pub mod pressure;
/// Finite, source-bound coolant storage and recirculating pump feedback.
pub mod recirculation;
/// Joint source-enthalpy storage in the solid and recirculating coolant.
pub mod stored_solid;

use fs_convection::{CorrelationId, CorrelationRobinBoundary, NusseltEvaluation, ThermalDirection};
use fs_exec::Cx;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{EquilibriumPhaseState, ResolvedHeatCapacityCurve};
use fs_material::state_point::MaterialStatePointError;
use fs_qty::{MassFlowRate, Temperature, VolumetricFlowRate};

use crate::graph::thermal::medium::{FlowHeatGeometry, FrozenFlowMedium, MediumFlowError};
use crate::graph::thermal::medium::series::{FlowHeatSection, MAX_MATERIAL_FLOW_SECTIONS};
use super::super::{LiquidEnthalpyError, checkpoint, march_region, positive, state_at_temperature};

/// Explicit regional property-evaluation approximation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LiquidPropertyTemperature {
    /// Arithmetic mean of entering and leaving bulk-fluid temperatures.
    BulkMean,
    /// Arithmetic mean of wall temperature and the above bulk mean.
    FilmMean,
}
impl LiquidPropertyTemperature {
    fn evaluate(self, inlet: f64, outlet: f64, wall: f64) -> f64 {
        let bulk = 0.5 * inlet + 0.5 * outlet;
        match self { Self::BulkMean => bulk, Self::FilmMean => 0.5 * wall + 0.5 * bulk }
    }
}

/// Region data; heating/cooling direction is derived from the actual wall state.
#[derive(Clone, Copy, Debug)]
pub struct LiquidHeatSection<'a> {
    /// Unique existing solid trace name.
    pub region: &'a str,
    /// Flow-normal and wetted areas are independent, as in the existing compiler.
    pub geometry: FlowHeatGeometry,
    /// Explicit existing correlation, never chosen from the material name.
    pub correlation: CorrelationId,
}

/// Bounded local coefficient/outlet solve. These are numerical, not source errors.
#[derive(Clone, Copy, Debug)]
pub struct LiquidFeedbackBudget {
    /// Maximum outlet fixed-point residual in kelvin.
    pub temperature_tolerance_k: f64,
    /// Maximum relative coefficient change on a final property re-evaluation.
    pub coefficient_relative_tolerance: f64,
    /// Bisection or heated-outlet prediction updates per region, in 1..=128
    /// (plus endpoint/final checks).
    pub max_region_iterations: usize,
}
impl LiquidFeedbackBudget {
    fn admit(self) -> Result<(), LiquidFeedbackError> {
        if !positive(self.temperature_tolerance_k)
            || !positive(self.coefficient_relative_tolerance)
            || self.coefficient_relative_tolerance >= 1.0
            || !(1..=128).contains(&self.max_region_iterations)
        {
            return Err(LiquidFeedbackError::Invalid("invalid local liquid feedback budget"));
        }
        Ok(())
    }
}

/// Original physics refusals or an explicitly unclosed local solve.
#[derive(Debug)]
pub enum LiquidFeedbackError {
    /// Invalid geometry, source binding or numerical control.
    Invalid(&'static str),
    /// Original source-resolved material refusal, without fallback to constants.
    Material(MaterialStatePointError),
    /// Original correlation/channel refusal.
    Flow(MediumFlowError),
    /// Original caloric/cancellation refusal.
    Enthalpy(LiquidEnthalpyError),
    /// No locally self-consistent property state in the admitted iteration budget.
    NotConverged { region: String },
    /// Joint material/solid solve exhausted its explicit work budget.
    CouplingNotConverged {
        /// Solid solves performed.
        iterations: usize,
        /// Unrelaxed wall residual, K.
        maximum_wall_residual_k: f64,
        /// Effective coolant-reference residual, K.
        maximum_reference_residual_k: f64,
        /// Relative h residual.
        maximum_coefficient_residual_relative: f64,
        /// Largest regional power mismatch, W.
        maximum_heat_imbalance_w: f64,
    },
}
impl core::fmt::Display for LiquidFeedbackError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "liquid material feedback: {self:?}")
    }
}
impl std::error::Error for LiquidFeedbackError {}
impl From<MaterialStatePointError> for LiquidFeedbackError {
    fn from(e: MaterialStatePointError) -> Self { Self::Material(e) }
}
impl From<MediumFlowError> for LiquidFeedbackError {
    fn from(e: MediumFlowError) -> Self { Self::Flow(e) }
}
impl From<LiquidEnthalpyError> for LiquidFeedbackError {
    fn from(e: LiquidEnthalpyError) -> Self { Self::Enthalpy(e) }
}

#[derive(Clone, Debug)]
struct Section {
    name: String,
    geometry: FlowHeatGeometry,
    correlation: CorrelationId,
}

/// One conserved mass stream and one fixed-context material law.
///
/// The source chart and liquid law must select the same Cp/density claims and
/// nonthermal coordinates. Wall and fluid temperatures must be in the caloric
/// span; every property trial must also be supported by its pinned claims. This
/// deliberately conservative admission never extrapolates a failed trial.
#[derive(Debug)]
pub struct LiquidHeatPath<'a, 'card> {
    law: &'a LiquidTemperatureLaw<'card>,
    source: &'a ResolvedHeatCapacityCurve,
    inlet: EquilibriumPhaseState,
    mass_flow: f64,
    sections: Vec<Section>,
    property_temperature: LiquidPropertyTemperature,
}

impl<'a, 'card> LiquidHeatPath<'a, 'card> {
    /// Compile using the existing source/path admission, without another fluid law.
    /// The mass rate stays fixed; local volume flow becomes `mdot/rho(T)`.
    ///
    /// # Errors
    /// Refuses source mismatches, invalid sections, unsupported correlations or inlet.
    pub fn try_new(
        cx: &Cx<'_>, law: &'a LiquidTemperatureLaw<'card>, source: &'a ResolvedHeatCapacityCurve,
        mass_flow: MassFlowRate, inlet: Temperature, sections: &[LiquidHeatSection<'_>],
        property_temperature: LiquidPropertyTemperature,
    ) -> Result<Self, LiquidFeedbackError> {
        checkpoint(cx)?;
        if !positive(mass_flow.value()) || sections.is_empty() || sections.len() > MAX_MATERIAL_FLOW_SECTIONS {
            return Err(LiquidFeedbackError::Invalid("positive mass flow and bounded nonempty sections are required"));
        }
        let descriptions: Vec<_> = sections.iter().map(|s| FlowHeatSection {
            region: s.region, geometry: s.geometry, correlation: s.correlation,
            direction: ThermalDirection::HeatingFluid,
        }).collect();
        let reference = law.reference_state();
        let frozen = FrozenFlowMedium::Liquid(reference).channel_series(
            VolumetricFlowRate::new(mass_flow.value() / reference.density_kg_m3()), inlet, &descriptions,
        )?;
        frozen.admit_caloric_source(cx, source, law.temperature_axis())?;
        let inlet = state_at_temperature(source.curve(), inlet.value())?;
        law.at_temperature(Temperature::new(inlet.temperature_k()))?;
        let sections = sections.iter().map(|s| Section {
            name: s.region.to_owned(), geometry: s.geometry, correlation: s.correlation,
        }).collect();
        checkpoint(cx)?;
        Ok(Self { law, source, inlet, mass_flow: mass_flow.value(), sections, property_temperature })
    }

    /// One conserved mass-flow rate [kg/s].
    #[must_use]
    pub const fn mass_flow_kg_per_s(&self) -> f64 { self.mass_flow }
    /// Entering source-caloric state.
    #[must_use]
    pub const fn inlet(&self) -> EquilibriumPhaseState { self.inlet }

    /// Recompute local density, viscosity, conductivity, Cp, Re, Pr and h from
    /// temperature, then transport the same material's enthalpy downstream.
    /// Every returned coefficient has been checked again at its predicted outlet.
    /// No candidate mutates the input or returns a partial path on refusal.
    ///
    /// # Errors
    /// Refuses unsupported trial states, local nonconvergence, or cancellation.
    pub fn march(
        &self, cx: &Cx<'_>, walls: &[f64], budget: LiquidFeedbackBudget,
    ) -> Result<LiquidHeatMarch<'a>, LiquidFeedbackError> {
        checkpoint(cx)?;
        budget.admit()?;
        if walls.len() != self.sections.len() {
            return Err(LiquidFeedbackError::Invalid("one wall temperature is required per section"));
        }
        // Validate all walls before doing trial work on the first region.
        for &wall in walls { state_at_temperature(self.source.curve(), wall)?; }
        let mut inlet = self.inlet;
        let mut regions = Vec::with_capacity(walls.len());
        let mut total = 0.0;
        for (section, &wall) in self.sections.iter().zip(walls) {
            let region = self.solve_region(cx, section, inlet, wall, budget)?;
            total += region.heat_rate_w;
            if !total.is_finite() { return Err(LiquidFeedbackError::Invalid("total liquid heat is unrepresentable")); }
            inlet = region.outlet;
            regions.push(region);
        }
        checkpoint(cx)?;
        Ok(LiquidHeatMarch { source: self.source, regions, mass_flow_kg_per_s: self.mass_flow, total_heat_rate_w: total })
    }

    fn solve_region(
        &self, cx: &Cx<'_>, section: &Section, inlet: EquilibriumPhaseState,
        wall: f64, budget: LiquidFeedbackBudget,
    ) -> Result<LiquidHeatRegion, LiquidFeedbackError> {
        let mut low = inlet.temperature_k().min(wall);
        let mut high = inlet.temperature_k().max(wall);
        for guess in [low, high] {
            let candidate = self.evaluate_region(cx, section, inlet, wall, guess)?;
            if let Some(accepted) = self.verify_region(cx, section, inlet, wall, candidate, budget)? {
                return Ok(accepted);
            }
        }
        for _ in 0..budget.max_region_iterations {
            checkpoint(cx)?;
            let guess = 0.5 * low + 0.5 * high;
            let candidate = self.evaluate_region(cx, section, inlet, wall, guess)?;
            let residual = candidate.temperature_residual_k;
            if let Some(accepted) = self.verify_region(cx, section, inlet, wall, candidate, budget)? {
                return Ok(accepted);
            }
            if guess == low || guess == high { break; }
            // At the physical lower endpoint prediction >= guess, at the upper
            // endpoint prediction <= guess. Keep this bracket, not a monotonicity
            // claim for the material law. Multiple nonlinear roots are not excluded.
            if residual > 0.0 { low = guess; } else { high = guess; }
        }
        Err(LiquidFeedbackError::NotConverged { region: section.name.clone() })
    }

    fn verify_region(
        &self, cx: &Cx<'_>, section: &Section, inlet: EquilibriumPhaseState, wall: f64,
        candidate: LiquidHeatRegion, budget: LiquidFeedbackBudget,
    ) -> Result<Option<LiquidHeatRegion>, LiquidFeedbackError> {
        if candidate.temperature_residual_k.abs() > budget.temperature_tolerance_k { return Ok(None); }
        let mut checked = self.evaluate_region(cx, section, inlet, wall, candidate.outlet.temperature_k())?;
        let before = candidate.boundary.coefficient().value.value();
        let after = checked.boundary.coefficient().value.value();
        checked.coefficient_residual_relative = (after - before).abs() / before.max(after);
        if checked.temperature_residual_k.abs() <= budget.temperature_tolerance_k
            && checked.coefficient_residual_relative <= budget.coefficient_relative_tolerance
        { Ok(Some(checked)) } else { Ok(None) }
    }

    fn evaluate_region(
        &self, cx: &Cx<'_>, section: &Section, inlet: EquilibriumPhaseState, wall: f64, guess: f64,
    ) -> Result<LiquidHeatRegion, LiquidFeedbackError> {
        checkpoint(cx)?;
        let property_temperature_k = self.property_temperature.evaluate(inlet.temperature_k(), guess, wall);
        let material = self.law.at_temperature(Temperature::new(property_temperature_k))?;
        let volume_flow = self.mass_flow / material.density_kg_m3();
        let channel = FrozenFlowMedium::Liquid(&material).channel(
            &section.name, VolumetricFlowRate::new(volume_flow), Temperature::new(inlet.temperature_k()),
            section.geometry, section.correlation,
            if wall >= inlet.temperature_k() { ThermalDirection::HeatingFluid } else { ThermalDirection::CoolingFluid },
        )?;
        let segment = &channel.path().segments()[0];
        let outlet = march_region(cx, self.source.curve(), inlet, self.mass_flow, segment, wall)?;
        let heat_rate_w = self.mass_flow * (outlet.specific_enthalpy_j_kg() - inlet.specific_enthalpy_j_kg());
        let reference = wall - heat_rate_w / (segment.htc_w_per_m2_k() * segment.area_m2());
        let velocity_m_s = volume_flow / section.geometry.flow_area.value();
        if !heat_rate_w.is_finite() || !positive(reference) || !positive(velocity_m_s) {
            return Err(LiquidFeedbackError::Invalid("regional heat, reference or velocity is unrepresentable"));
        }
        let boundary = channel.robin_boundary(Temperature::new(reference))?;
        let correlation = channel.correlation().clone();
        let temperature_residual_k = outlet.temperature_k() - guess;
        checkpoint(cx)?;
        Ok(LiquidHeatRegion {
            region: section.name.clone(), material, property_temperature_k,
            volume_flow_m3_s: volume_flow, velocity_m_s, inlet, outlet, heat_rate_w,
            boundary, correlation, temperature_residual_k, coefficient_residual_relative: 0.0,
        })
    }
}

/// Accepted regional snapshot. Source and correlation evidence remain distinct.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidHeatRegion {
    /// Original solid trace.
    pub region: String,
    /// Actual newly resolved liquid and its four source receipts.
    pub material: ResolvedLiquidState,
    /// Temperature used to evaluate these material properties, K.
    pub property_temperature_k: f64,
    /// Volume flow adjusted to preserve the prescribed mass flow, m3/s.
    pub volume_flow_m3_s: f64,
    /// Local speed from the same density and flow area, m/s.
    pub velocity_m_s: f64,
    /// Entering caloric state.
    pub inlet: EquilibriumPhaseState,
    /// Leaving caloric state.
    pub outlet: EquilibriumPhaseState,
    /// Signed wall heat received by the fluid, W. In a dissipative march the
    /// separate hydraulic report supplies the additional pressure-loss heat.
    pub heat_rate_w: f64,
    /// Actual energy-equivalent Robin boundary, not a separately authored h.
    pub boundary: CorrelationRobinBoundary,
    /// Original admitted Re/Pr/Nu point and correlation evidence.
    pub correlation: NusseltEvaluation,
    /// Signed local outlet fixed-point residual, K.
    pub temperature_residual_k: f64,
    /// Relative h change from the final outlet-based property re-evaluation.
    pub coefficient_residual_relative: f64,
}

/// Complete source-backed coolant trajectory after local material feedback.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidHeatMarch<'a> {
    /// Same Cp/density chart as the transport source; never substituted.
    pub source: &'a ResolvedHeatCapacityCurve,
    /// Regions in original inlet-to-outlet order.
    pub regions: Vec<LiquidHeatRegion>,
    /// One conserved mass-flow rate, kg/s.
    pub mass_flow_kg_per_s: f64,
    /// Sum of regional signed wall heat transfers, W.
    pub total_heat_rate_w: f64,
}
