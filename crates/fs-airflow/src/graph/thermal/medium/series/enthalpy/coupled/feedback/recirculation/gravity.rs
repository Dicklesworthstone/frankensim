//! Finite-inventory thermosiphons on the original liquid pressure/enthalpy path.
//!
//! The reservoir owns storage; the loop remains quasi-steady with no pipe
//! residence time or inertial state. At every implicit inlet trial, source rho(T)
//! selects gravity head and source mu(T) selects hydraulic loss. This module
//! connects those existing owners to the existing reservoir energy integrator.

use fs_exec::Cx;
use fs_material::phase::EquilibriumPhaseState;
use fs_qty::Pressure;

use crate::graph::thermal::medium::series::hydraulic::pressure::SeriesGravity;
use super::super::{LiquidFeedbackBudget, LiquidHeatMarch, LiquidHeatPath, checkpoint};
use super::super::pressure::{
    LiquidResistanceSection, PressureDrivenBudget, PressureDrivenError, PressureDrivenSolution,
};
use super::super::pressure::direction::{DirectedLiquidSolution, LiquidFlowSelection};
use super::{
    LiquidReservoir, RecirculatingLiquidStep, RecirculationError, ReservoirCircuit,
    ReservoirExchange, ReservoirStepControls, integrate, limited_power,
};

/// Existing local/pressure/storage controls plus one shared flow-work limit.
/// The pressure bracket bounds a positive flow magnitude and is not expanded.
/// The directional entry point selects its orientation separately.
#[derive(Clone, Copy, Debug)]
pub struct GravityReservoirControls {
    pub local: LiquidFeedbackBudget,
    pub pressure: PressureDrivenBudget,
    pub reservoir: ReservoirStepControls,
    /// Actual thermal-path evaluations across ALL implicit inlet trials.
    pub max_total_flow_evaluations: usize,
}

impl<'a, 'card> LiquidHeatPath<'a, 'card> {
    /// Advance finite coolant storage through a closed buoyant circulation path.
    ///
    /// Every reservoir-temperature trial recomputes the original material,
    /// enthalpy trajectory, gravity head and pressure-selected flow. Zero imposed
    /// pressure supports passive circulation; a signed pressure may assist or
    /// oppose it. The accepted circuit uses exactly the committed reservoir inlet.
    ///
    /// The first and last geometry vertices must coincide. Conservative gravity
    /// then transfers no net mechanical energy to the reservoir; hydraulic heat
    /// remains omitted only within BOTH the watt and step-joule allowances.
    /// Reservoir mass and source pressure are fixed; volume may vary. Walls are
    /// prescribed thermal reservoirs, and only this coolant inventory stores heat.
    ///
    /// This is a positive-orientation, quasi-steady single-phase model. A missing
    /// positive root refuses atomically; no zero-flow conduction, flow reversal,
    /// inertial stability, pipe storage or complete mechanical balance is implied.
    /// Every failure/cancellation leaves reservoir state AND time unchanged.
    pub fn advance_reservoir_with_gravity(
        &self, cx: &Cx<'_>, reservoir: &mut LiquidReservoir<'a>, end_time_s: f64,
        external_heat_j: f64, walls: &[f64], resistances: &[LiquidResistanceSection<'_>],
        drive: Pressure, gravity: &SeriesGravity, controls: GravityReservoirControls,
    ) -> Result<RecirculatingLiquidStep<PressureDrivenSolution<LiquidHeatMarch<'a>>>, RecirculationError> {
        self.advance_gravity_reservoir(cx, reservoir, end_time_s, external_heat_j, gravity,
            controls, |path, pressure| path.solve_pressure_gravity_driven(cx, walls,
                resistances, drive, gravity, controls.local, pressure))
    }

    /// Advance the same finite inventory with an explicitly selected circulation
    /// direction, or the sole admitted orientation bracket. Reverse flow revisits
    /// the ORIGINAL regions in reverse order and recomputes every material state.
    /// Heat remains the actual return-minus-supply enthalpy, irrespective of sign.
    ///
    /// Both orientation searches share the step's flow budget. Every implicit
    /// inlet trial must satisfy the requested selection; ambiguity, a missing
    /// bracket or a source failure leaves inventory and time unchanged. Switching
    /// branches between accepted steps is a quasi-steady change, not a resolved
    /// passage through zero flow or an inertial/pipe-storage trajectory.
    pub fn advance_reservoir_circulation(
        &self, cx: &Cx<'_>, reservoir: &mut LiquidReservoir<'a>, end_time_s: f64,
        external_heat_j: f64, walls: &[f64], resistances: &[LiquidResistanceSection<'_>],
        drive: Pressure, gravity: &SeriesGravity, selection: LiquidFlowSelection,
        controls: GravityReservoirControls,
    ) -> Result<RecirculatingLiquidStep<DirectedLiquidSolution<'a>>, RecirculationError> {
        self.advance_gravity_reservoir(cx, reservoir, end_time_s, external_heat_j, gravity,
            controls, |path, pressure| path.solve_gravity_circulation(cx, walls,
                resistances, drive, gravity, selection, controls.local, pressure))
    }

    fn advance_gravity_reservoir<T: ReservoirCircuit>(
        &self, cx: &Cx<'_>, reservoir: &mut LiquidReservoir<'a>, end_time_s: f64,
        external_heat_j: f64, gravity: &SeriesGravity, controls: GravityReservoirControls,
        mut evaluate: impl FnMut(&Self, PressureDrivenBudget) -> Result<T, PressureDrivenError>,
    ) -> Result<RecirculatingLiquidStep<T>, RecirculationError> {
        let duration = self.admit_reservoir_step(cx, reservoir, end_time_s, external_heat_j, controls.reservoir)?;
        if gravity.section_count() != self.sections.len()
            || gravity.vertices_m().first() != gravity.vertices_m().last()
        {
            return Err(RecirculationError::Invalid("recirculation needs one gravity segment per region and coincident endpoints"));
        }
        if controls.max_total_flow_evaluations < 2
            || !(2..=128).contains(&controls.pressure.max_flow_evaluations)
        {
            return Err(RecirculationError::Invalid("invalid local or shared circulation work budget"));
        }
        let mut pressure = controls.pressure;
        pressure.maximum_neglected_hydraulic_power_w = limited_power(
            pressure.maximum_neglected_hydraulic_power_w, controls.reservoir, duration,
        )?;
        let result = integrate(cx, reservoir, end_time_s, external_heat_j, controls.reservoir,
            |temperature, performed, _| {
                let remaining = controls.max_total_flow_evaluations.saturating_sub(performed);
                if remaining < 2 {
                    return Err(RecirculationError::WorkBudget { flow_evaluations: performed, solid_solves: 0 });
                }
                let path = self.at_reservoir_inlet(cx, temperature)?;
                Ok(evaluate(&path, PressureDrivenBudget {
                    max_flow_evaluations: pressure.max_flow_evaluations.min(remaining), ..pressure
                })?)
            })?;
        checkpoint(cx)?;
        reservoir.state = result.final_state;
        reservoir.time_s = end_time_s;
        Ok(result)
    }
}

impl ReservoirCircuit for PressureDrivenSolution<LiquidHeatMarch<'_>> {
    fn exchange(&self, inlet: EquilibriumPhaseState) -> ReservoirExchange {
        ReservoirExchange {
            net_heat_rate_w: self.thermal.total_heat_rate_w,
            dissipated_power_w: self.hydraulics.total_dissipated_power_w,
            inlet_matches: self.thermal.regions.first().is_some_and(|r| r.inlet == inlet),
        }
    }
    fn flow_evaluations(&self) -> usize { self.flow_evaluations }
    fn solid_solves(&self) -> usize { 0 }
}

impl ReservoirCircuit for DirectedLiquidSolution<'_> {
    fn exchange(&self, inlet: EquilibriumPhaseState) -> ReservoirExchange {
        self.operating_point().exchange(inlet)
    }
    fn flow_evaluations(&self) -> usize { self.total_flow_evaluations() }
    fn solid_solves(&self) -> usize { 0 }
}
