//! Original pump head data select the native source-bound series flow.
//! The existing pump/network matcher owns bracketing and pressure acceptance.
//! Source material resolution occurs before matching; trials borrow checked
//! snapshots and final publication moves the ORIGINAL liquid owners unchanged.

use super::*;
use fs_exec::Cx;
use fs_material::liquid::ResolvedLiquidState;
use crate::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{
    PumpHeadCurve, PumpState, PumpOperatingControls, PumpOperatingError,
    PumpOperatingPoint, match_pump,
};

struct Trial<'a> {
    flow: SeriesPressureFlow<&'a FluidTransportState>,
    inlet: &'a ResolvedLiquidState,
}
impl PumpState for Trial<'_> {
    fn inlet_material(&self) -> &ResolvedLiquidState { self.inlet }
    fn inlet_volume_flow_m3_s(&self) -> f64 { self.flow.evaluations()[0].volume_flow_m3_s }
    fn driving_pressure_pa(&self) -> f64 { self.flow.driving_pressure_pa() }
    fn flow_evaluations(&self) -> usize { 1 }
}

impl SeriesPressureFlow<ResolvedLiquidState> {
    /// Intersect the ORIGINAL fixed-speed pump curve with these liquid sections.
    /// Section zero is the declared pump inlet; its actual density converts head
    /// to pressure and its viscosity gates the original equipment applicability.
    /// All sections carry the same mass rate, with their own local volume rates.
    ///
    /// Reuses the same bounded matcher as steady parallel/common-solid cooling.
    /// Every pressure candidate calls the original series law and Reynolds gates.
    /// Each original source supplies one checked snapshot; candidates borrow
    /// those snapshots, not cloned material cards. The returned series object
    /// retains the original owners and the EXACT matched candidate evaluations,
    /// with no post-acceptance replay. Its existing convection methods apply.
    ///
    /// At most controls.max_network_evaluations series solves are admitted.
    /// Each solve is O(sections); no independent pressure matrix or hydraulic
    /// time integrator is added. This is frozen-state nominal matching, not a
    /// transient pump, speed scaling, efficiency, compressibility or cavitation
    /// model. Hydraulic power is reported, NOT automatically deposited as heat.
    /// Negative/zero-pressure branches require other equipment data/owners;
    /// no extrapolation or last-good-flow fallback occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_pump<'p>(
        sections: Vec<PressureSection<ResolvedLiquidState>>,
        pump: &'p PumpHeadCurve,
        series_budget: SeriesPressureBudget,
        controls: PumpOperatingControls,
        cx: &Cx<'_>,
        gate: &CancelGate,
    ) -> core::result::Result<PumpOperatingPoint<'p, Self>, PumpOperatingError> {
        let checkpoint = || {
            if cx.checkpoint().is_err() { gate.request(); }
            poll(gate).map_err(PumpOperatingError::Series)
        };
        checkpoint()?;
        if sections.is_empty() || sections.len() > series_budget.maximum_sections {
            return Err(PumpOperatingError::Series(invalid("nonempty bounded pump series required")));
        }
        let mut snapshots = Vec::with_capacity(sections.len());
        for (index, section) in sections.iter().enumerate() {
            checkpoint()?;
            snapshots.push(section.medium.transport_state().map_err(|source|
                PumpOperatingError::Series(SeriesPressureError::Section {
                    index, source: ResistanceError::Material(source),
                }))?);
        }
        let matched = match_pump(cx, pump, controls, |drive| {
            checkpoint()?;
            let borrowed = sections.iter().zip(&snapshots).map(|(section, snapshot)|
                PressureSection { resistance: section.resistance, medium: snapshot }).collect();
            let flow = SeriesPressureFlow::solve(borrowed, drive.value(), series_budget, gate)
                .map_err(PumpOperatingError::Series)?;
            checkpoint()?;
            Ok(Trial { flow, inlet: &sections[0].medium })
        })?;
        let PumpOperatingPoint { curve, state, supply, network_pressure_pa,
            pressure_residual_pa, network_evaluations, total_flow_evaluations,
            total_solid_solves } = matched;
        let Trial { flow, inlet: _ } = state;
        let SeriesPressureFlow { sections: borrowed, evaluations, driving_pressure_pa,
            mass_flow_kg_s, pressure_residual_pa: series_residual, dissipated_power_w,
            dpressure_dmass_flow, gravity, hydrostatic, boussinesq } = flow;
        drop(borrowed);
        checkpoint()?;
        let state = Self { sections, evaluations, driving_pressure_pa, mass_flow_kg_s,
            pressure_residual_pa: series_residual, dissipated_power_w,
            dpressure_dmass_flow, gravity, hydrostatic, boussinesq };
        Ok(PumpOperatingPoint { curve, state, supply, network_pressure_pa,
            pressure_residual_pa, network_evaluations, total_flow_evaluations,
            total_solid_solves })
    }
}
