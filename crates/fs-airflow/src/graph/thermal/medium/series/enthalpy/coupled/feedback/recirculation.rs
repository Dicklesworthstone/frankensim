//! Finite coolant inventory closing a recirculating, quasi-steady flow circuit.
//!
//! At fixed thermodynamic source pressure, a well-mixed, constant-mass reservoir
//! obeys M dh/dt = sum(mdot * (h_return - h_supply)) + P_external. This owner
//! uses backward Euler on that balance and the EXISTING source enthalpy chart.
//! Each inlet trial re-solves the actual pump/branch problem at that inlet.
//! No frozen Cp, temperature averaging, replacement pump or fluid law is used.
//!
//! Only the reservoir stores energy here. Pipes/manifolds have zero residence
//! time/storage; volume may change at fixed pressure. This is not a rigid sealed
//! vessel, boiling model or full mechanical-energy balance. Branch hydraulic
//! dissipation remains omitted and is bounded in both watts and step joules.

use fs_exec::Cx;
use fs_material::phase::{EquilibriumPhaseState, ResolvedHeatCapacityCurve};
use fs_qty::{MassFlowRate, Temperature};

use super::{LiquidFeedbackBudget, LiquidFeedbackError, LiquidHeatPath, LiquidHeatSection,
    LiquidEnthalpyError, checkpoint, positive, state_at_temperature};
use super::pressure::porous::LiquidResistanceModel;
use super::pressure::parallel::{ParallelLiquidBranch, ParallelLiquidNetwork,
    ParallelLiquidSolution, ParallelPressureBudget, ParallelPressureError};
use super::pressure::parallel::pump::{PumpHeadCurve, PumpOperatingControls,
    PumpOperatingError, PumpOperatingPoint};

/// Committed inventory on one immutable caloric/source chart.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidReservoir<'a> {
    source: &'a ResolvedHeatCapacityCurve,
    state: EquilibriumPhaseState,
    mass_kg: f64,
    time_s: f64,
}
impl<'a> LiquidReservoir<'a> {
    /// Original Cp/density chart and its source receipts.
    pub const fn source(&self) -> &'a ResolvedHeatCapacityCurve { self.source }
    /// Committed source-bound state; not an unaccepted inlet trial.
    pub const fn state(&self) -> EquilibriumPhaseState { self.state }
    /// Constant reservoir mass, kg; loop circulation is not a mass source.
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Time of the last accepted state, seconds.
    pub const fn time_s(&self) -> f64 { self.time_s }
}

/// One implicit reservoir step, with explicit physical and numerical limits.
#[derive(Clone, Copy, Debug)]
pub struct ReservoirStepControls {
    /// Supported inlet-temperature bracket, K; not expanded or extrapolated.
    pub temperature_bracket_k: (f64, f64),
    /// Absolute balance residual, J. A small temperature bracket is insufficient.
    pub energy_absolute_tolerance_j: f64,
    /// Maximum complete pump/circuit evaluations, including endpoints, 2..=128.
    pub max_inlet_evaluations: usize,
    /// Maximum ignored hydraulic dissipation over this step, J. Together with
    /// the network's watt cap this also limits every nested pressure trial.
    pub maximum_neglected_hydraulic_energy_j: f64,
}

/// No failure publishes either a new reservoir state or a candidate circuit.
#[derive(Debug)]
pub enum RecirculationError {
    Invalid(&'static str),
    Physics(PumpOperatingError),
    NotBracketed { lower_residual_j: f64, upper_residual_j: f64 },
    NotConverged { evaluations: usize, last_residual_j: f64 },
    WorkBudget { flow_evaluations: usize, solid_solves: usize },
}
impl core::fmt::Display for RecirculationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "recirculating liquid: {self:?}")
    }
}
impl std::error::Error for RecirculationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Physics(e) => Some(e), _ => None }
    }
}
impl From<PumpOperatingError> for RecirculationError {
    fn from(e: PumpOperatingError) -> Self { Self::Physics(e) }
}
impl From<ParallelPressureError> for RecirculationError {
    fn from(e: ParallelPressureError) -> Self { Self::Physics(e.into()) }
}
impl From<LiquidFeedbackError> for RecirculationError {
    fn from(e: LiquidFeedbackError) -> Self { Self::from(ParallelPressureError::from(e)) }
}
impl From<LiquidEnthalpyError> for RecirculationError {
    fn from(e: LiquidEnthalpyError) -> Self { Self::from(ParallelPressureError::from(e)) }
}

/// Complete accepted step, retaining the original circuit/solid payload.
#[derive(Debug)]
pub struct RecirculatingLiquidStep<T> {
    pub circuit: T,
    pub initial: EquilibriumPhaseState,
    pub final_state: EquilibriumPhaseState,
    pub start_time_s: f64,
    pub end_time_s: f64,
    /// M * (h_new - h_old), J; the enthalpy origin cancels before multiplying.
    pub stored_enthalpy_change_j: f64,
    /// Integrated NET return-minus-supply enthalpy flux, J.
    pub returned_heat_j: f64,
    /// Caller-supplied reservoir heat over this step, counted once, J.
    pub external_heat_j: f64,
    /// Stored change minus return heat minus external heat, J.
    pub energy_residual_j: f64,
    /// Reported hydraulic energy omitted from the caloric balance, J.
    pub neglected_hydraulic_energy_j: f64,
    pub inlet_evaluations: usize,
    pub total_flow_evaluations: usize,
    pub total_solid_solves: usize,
}

impl<'p, 'a, 'card, H: LiquidResistanceModel> ParallelLiquidNetwork<'p, 'a, 'card, H> {
    /// Bind an initial finite inventory to this network's actual liquid source.
    /// The reservoir is well mixed, fixed mass and at the same constant source
    /// pressure. Its initial temperature need not equal the template inlet.
    pub fn reservoir(
        &self, cx: &Cx<'_>, mass_kg: f64, temperature: Temperature, time_s: f64,
    ) -> Result<LiquidReservoir<'a>, RecirculationError> {
        checkpoint(cx)?;
        if !positive(mass_kg) || !time_s.is_finite() {
            return Err(RecirculationError::Invalid("reservoir needs positive finite mass and finite time"));
        }
        let path = self.branches()[0].path;
        let state = state_at_temperature(path.source.curve(), temperature.value())?;
        path.law.at_temperature(temperature).map_err(LiquidFeedbackError::from)?;
        checkpoint(cx)?;
        Ok(LiquidReservoir { source: path.source, state, mass_kg, time_s })
    }

    /// Advance a recirculating reservoir through the existing pump/fixed-wall
    /// circuit. Each trial changes the COMMON supply temperature, recomputing
    /// inlet density, pump applicability, flow split, transport and outlet mix.
    /// `external_heat_j` is a separate signed step input to the reservoir only.
    /// Work is shared across every inlet, pressure and branch trial in the step.
    /// The reservoir changes only on success; the returned circuit was evaluated
    /// at exactly the committed reservoir temperature, not an unevaluated guess.
    pub fn advance_reservoir_with_pump<'s>(
        &self, cx: &Cx<'_>, reservoir: &mut LiquidReservoir<'a>, end_time_s: f64,
        external_heat_j: f64, walls: &[f64], pump: &'s PumpHeadCurve,
        local: LiquidFeedbackBudget, work: ParallelPressureBudget,
        pump_controls: PumpOperatingControls, step: ReservoirStepControls,
    ) -> Result<RecirculatingLiquidStep<PumpOperatingPoint<'s, ParallelLiquidSolution<'a>>>, RecirculationError> {
        let duration = self.admit_reservoir_step(cx, reservoir, end_time_s, external_heat_j, step)?;
        let mut work = work;
        work.maximum_neglected_hydraulic_power_w = limited_power(work.maximum_neglected_hydraulic_power_w, step, duration)?;
        let result = integrate(cx, reservoir, end_time_s, external_heat_j, step, |temperature, flows, _| {
            let remaining = work.max_flow_evaluations.saturating_sub(flows);
            if remaining < 2 { return Err(RecirculationError::WorkBudget { flow_evaluations: flows, solid_solves: 0 }); }
            self.with_reservoir_inlet(cx, temperature, |network| {
                network.solve_with_pump(cx, walls, pump, local,
                    ParallelPressureBudget { max_flow_evaluations: remaining, ..work }, pump_controls)
            })
        })?;
        checkpoint(cx)?;
        reservoir.state = result.final_state;
        reservoir.time_s = end_time_s;
        Ok(result)
    }

    fn admit_reservoir_step(
        &self, cx: &Cx<'_>, reservoir: &LiquidReservoir<'a>, end: f64,
        external: f64, controls: ReservoirStepControls,
    ) -> Result<f64, RecirculationError> {
        checkpoint(cx)?;
        let (lo, hi) = controls.temperature_bracket_k;
        let duration = end - reservoir.time_s;
        let source = self.branches()[0].path.source;
        if source.curve().identity() != reservoir.source.curve().identity()
            || !end.is_finite() || !positive(duration) || !external.is_finite()
            || !positive(lo) || !positive(hi) || lo >= hi
            || !positive(controls.energy_absolute_tolerance_j)
            || !(2..=128).contains(&controls.max_inlet_evaluations)
            || !controls.maximum_neglected_hydraulic_energy_j.is_finite()
            || controls.maximum_neglected_hydraulic_energy_j < 0.0
        { return Err(RecirculationError::Invalid("incompatible reservoir source, time interval or step controls")); }
        // Complete source preflight precedes the first flow or solid evaluation.
        for t in [lo, hi] {
            state_at_temperature(source.curve(), t)?;
            self.branches()[0].path.law.at_temperature(Temperature::new(t))
                .map_err(LiquidFeedbackError::from)?;
        }
        Ok(duration)
    }

    // The original laws, chart, region descriptions, and resistance assignments
    // survive inlet rebinding. Local paths are only descriptions, not new sources.
    fn with_reservoir_inlet<T>(
        &self, cx: &Cx<'_>, temperature: Temperature,
        evaluate: impl FnOnce(&ParallelLiquidNetwork<'_, 'a, 'card, H>) -> Result<T, PumpOperatingError>,
    ) -> Result<T, RecirculationError> {
        let mut paths = Vec::with_capacity(self.branches().len());
        for branch in self.branches() {
            checkpoint(cx)?;
            let path = branch.path;
            let sections: Vec<_> = path.sections.iter().map(|s| LiquidHeatSection {
                region: &s.name, geometry: s.geometry, correlation: s.correlation,
            }).collect();
            paths.push(LiquidHeatPath::try_new(cx, path.law, path.source,
                MassFlowRate::new(path.mass_flow), temperature, &sections, path.property_temperature)?);
        }
        let branches: Vec<_> = self.branches().iter().zip(&paths).map(|(old, path)| ParallelLiquidBranch {
            name: old.name, path, resistances: old.resistances, pressure_budget: old.pressure_budget,
        }).collect();
        let network = ParallelLiquidNetwork::try_new(cx, &branches)?;
        Ok(evaluate(&network)?)
    }
}

fn limited_power(watts: f64, controls: ReservoirStepControls, duration: f64) -> Result<f64, RecirculationError> {
    if !watts.is_finite() || watts < 0.0 {
        return Err(RecirculationError::Invalid("invalid original hydraulic-power allowance"));
    }
    Ok(watts.min(controls.maximum_neglected_hydraulic_energy_j / duration))
}

// Private adapters expose actual retained results, not caller-authored heat or
// work counters. The common-solid adapter uses the same time/balance algorithm.
trait ReservoirCircuit {
    fn fluid(&self) -> &ParallelLiquidSolution<'_>;
    fn flow_evaluations(&self) -> usize;
    fn solid_solves(&self) -> usize;
}
impl ReservoirCircuit for PumpOperatingPoint<'_, ParallelLiquidSolution<'_>> {
    fn fluid(&self) -> &ParallelLiquidSolution<'_> { &self.state }
    fn flow_evaluations(&self) -> usize { self.total_flow_evaluations }
    fn solid_solves(&self) -> usize { self.total_solid_solves }
}

fn integrate<T: ReservoirCircuit>(
    cx: &Cx<'_>, reservoir: &LiquidReservoir<'_>, end: f64, external: f64,
    controls: ReservoirStepControls,
    mut evaluate: impl FnMut(Temperature, usize, usize) -> Result<T, RecirculationError>,
) -> Result<RecirculatingLiquidStep<T>, RecirculationError> {
    let duration = end - reservoir.time_s;
    let mut total_flow = 0usize;
    let mut total_solid = 0usize;
    let mut probe = |temperature: f64, calls: usize| -> Result<RecirculatingLiquidStep<T>, RecirculationError> {
        checkpoint(cx)?;
        let final_state = state_at_temperature(reservoir.source.curve(), temperature)?;
        let circuit = evaluate(Temperature::new(temperature), total_flow, total_solid)?;
        let fluid = circuit.fluid();
        // Net circuit uptake already subtracts the supply enthalpy. Counting the
        // full returning mdot*h here would heat the tank merely by circulating it.
        let stored = reservoir.mass_kg * (final_state.specific_enthalpy_j_kg() - reservoir.state.specific_enthalpy_j_kg());
        let returned = duration * fluid.total_heat_rate_w;
        let neglected = duration * fluid.total_dissipated_power_w;
        let residual = (stored - returned) - external;
        if [stored, returned, neglected, residual].iter().any(|v| !v.is_finite())
            || neglected > controls.maximum_neglected_hydraulic_energy_j
            || fluid.branches.iter().any(|b| b.solution.thermal.regions.first()
                .is_none_or(|r| r.inlet != final_state))
        { return Err(RecirculationError::Invalid("inconsistent inlet or unrepresentable/excessive step energy")); }
        total_flow = total_flow.checked_add(circuit.flow_evaluations())
            .ok_or(RecirculationError::Invalid("step flow-evaluation count overflow"))?;
        total_solid = total_solid.checked_add(circuit.solid_solves())
            .ok_or(RecirculationError::Invalid("step solid-evaluation count overflow"))?;
        checkpoint(cx)?;
        Ok(RecirculatingLiquidStep { circuit, initial: reservoir.state, final_state,
            start_time_s: reservoir.time_s, end_time_s: end,
            stored_enthalpy_change_j: stored, returned_heat_j: returned, external_heat_j: external,
            energy_residual_j: residual, neglected_hydraulic_energy_j: neglected,
            inlet_evaluations: calls, total_flow_evaluations: total_flow, total_solid_solves: total_solid })
    };
    let closed = |s: &RecirculatingLiquidStep<T>| s.energy_residual_j.abs() <= controls.energy_absolute_tolerance_j;
    let (mut lo, mut hi) = controls.temperature_bracket_k;
    let lower = probe(lo, 1)?;
    if closed(&lower) { return Ok(lower); }
    let mut lower_residual = lower.energy_residual_j;
    drop(lower);
    let upper = probe(hi, 2)?;
    if closed(&upper) { return Ok(upper); }
    let mut last = upper.energy_residual_j;
    let mut upper_residual = last;
    if lower_residual.is_sign_positive() == last.is_sign_positive() {
        return Err(RecirculationError::NotBracketed { lower_residual_j: lower_residual, upper_residual_j: last });
    }
    drop(upper);
    let mut evaluations = 2;
    for calls in 3..=controls.max_inlet_evaluations {
        // Near-linear storage should not require dozens of complete solid/pump
        // solves. Use an interior secant, with a mandatory midpoint every fifth
        // evaluation to prevent endpoint stagnation for strongly curved maps.
        let scale = lower_residual.abs().max(upper_residual.abs());
        let left = lower_residual.abs()/scale;
        let right = upper_residual.abs()/scale;
        let fraction = left/(left+right);
        let secant = (1.0-fraction)*lo + fraction*hi;
        let middle = if calls % 5 != 0 && secant > lo && secant < hi {
            secant
        } else { 0.5*lo + 0.5*hi };
        if middle == lo || middle == hi { break; }
        let result = probe(middle, calls)?;
        evaluations = calls;
        if closed(&result) { return Ok(result); }
        last = result.energy_residual_j;
        if last.is_sign_positive() == lower_residual.is_sign_positive() {
            lo = middle; lower_residual = last;
        } else { hi = middle; upper_residual = last; }
    }
    Err(RecirculationError::NotConverged { evaluations, last_residual_j: last })
}

use fs_convection::CorrelationRobinBoundary;
use crate::AirflowError;
use crate::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use super::pressure::parallel::solid::CoupledParallelLiquidSolution;
use super::pressure::parallel::pump::solid::PumpCouplingControls;

impl<R> ReservoirCircuit for PumpOperatingPoint<'_, CoupledParallelLiquidSolution<'_, R>> {
    fn fluid(&self) -> &ParallelLiquidSolution<'_> { &self.state.fluid }
    fn flow_evaluations(&self) -> usize { self.total_flow_evaluations }
    fn solid_solves(&self) -> usize { self.total_solid_solves }
}

impl<'p, 'a, 'card, H: LiquidResistanceModel> ParallelLiquidNetwork<'p, 'a, 'card, H> {
    /// Advance finite coolant storage while resolving ONE common solid field,
    /// the pump operating point and all material-driven branches at each inlet.
    /// Reuses the same enthalpy time step as the fixed-wall path; every original
    /// regional heat, coefficient, wall, reference, mass and pressure gate stays
    /// active. The exact accepted solid payload and its sent boundaries survive.
    ///
    /// Every callback must evaluate from the SAME committed solid history for
    /// this interval, without externally committing a trial. On success callers
    /// may adopt the returned solid payload; on failure the reservoir is unchanged.
    /// A steady solid callback makes only the reservoir transient. A transient
    /// solid callback must itself implement the same backward-Euler interval;
    /// this adapter does not invent or certify its storage discretization.
    /// Flow work and actual solid calls are capped across ALL inlet trials, not
    /// reset for each pressure solve. Pump heat remains explicitly omitted.
    pub fn advance_reservoir_with_pump_coupled_iqn<'s, R, F>(
        &self, cx: &Cx<'_>, reservoir: &mut LiquidReservoir<'a>, end_time_s: f64,
        external_heat_j: f64, initial_walls_k: &[f64], pump: &'s PumpHeadCurve,
        controls: &PumpCouplingControls, step: ReservoirStepControls, mut solid: F,
    ) -> Result<RecirculatingLiquidStep<PumpOperatingPoint<'s, CoupledParallelLiquidSolution<'a, R>>>, RecirculationError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        let duration = self.admit_reservoir_step(cx, reservoir, end_time_s, external_heat_j, step)?;
        if controls.max_total_solid_solves == 0 {
            return Err(RecirculationError::Invalid("reservoir solid callback budget must be positive"));
        }
        let mut controls = *controls;
        controls.network.work.maximum_neglected_hydraulic_power_w = limited_power(
            controls.network.work.maximum_neglected_hydraulic_power_w, step, duration)?;
        let result = integrate(cx, reservoir, end_time_s, external_heat_j, step, |temperature, flows, solids| {
            let remaining_flow = controls.network.work.max_flow_evaluations.saturating_sub(flows);
            let remaining_solid = controls.max_total_solid_solves.saturating_sub(solids);
            if remaining_flow < 2 || remaining_solid == 0 {
                return Err(RecirculationError::WorkBudget { flow_evaluations: flows, solid_solves: solids });
            }
            let mut inner = controls;
            inner.network.work.max_flow_evaluations = remaining_flow;
            inner.max_total_solid_solves = remaining_solid;
            self.with_reservoir_inlet(cx, temperature, |network| {
                network.solve_with_pump_coupled_iqn(cx, initial_walls_k, pump, &inner, &mut solid)
            })
        })?;
        checkpoint(cx)?;
        reservoir.state = result.final_state;
        reservoir.time_s = end_time_s;
        Ok(result)
    }
}
