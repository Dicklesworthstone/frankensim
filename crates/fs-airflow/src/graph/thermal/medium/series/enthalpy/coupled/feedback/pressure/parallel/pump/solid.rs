//! Pump selection around the SAME common-solid/material/branch coupling solve.
//!
//! Every source-pressure trial re-solves the complete common solid. Callbacks
//! must start from the same accepted physical history; failed trials must not
//! publish external state. No independent solid model or acceleration is added.

use fs_convection::CorrelationRobinBoundary;
use fs_exec::Cx;
use crate::AirflowError;
use crate::graph::thermal::medium::series::enthalpy::coupled::LiquidSolidResponse;
use super::super::LiquidResistanceModel;
use super::super::solid::{CoupledParallelLiquidSolution, ParallelLiquidCouplingControls};
use super::{PumpHeadCurve, PumpLoad, PumpOperatingControls, PumpOperatingError,
    PumpOperatingPoint, ParallelLiquidNetwork, ParallelLiquidSolution,
    ParallelPressureError, match_pump};

/// All inner gates remain in force; resource caps span ALL source-pressure trials.
#[derive(Clone, Copy, Debug)]
pub struct PumpCouplingControls {
    /// Existing common-solid, branch-flow, material, power and work controls.
    /// Its work allowance is shared across the entire pump solve, not reset.
    pub network: ParallelLiquidCouplingControls,
    /// Pump-pressure root bracket and residual controls.
    pub operating_point: PumpOperatingControls,
    /// Maximum actual solid callbacks over every pressure candidate; nonzero.
    pub max_total_solid_solves: usize,
}

impl<R> PumpLoad for CoupledParallelLiquidSolution<'_, R> {
    fn fluid(&self) -> &ParallelLiquidSolution<'_> { &self.fluid }
    fn flow_evaluations(&self) -> usize { self.total_flow_evaluations }
    fn solid_solves(&self) -> usize { self.iterations }
}

impl<'p, 'a, 'card, H: LiquidResistanceModel> ParallelLiquidNetwork<'p, 'a, 'card, H> {
    /// Intersect a pump curve with fully coupled solid/material/branch states.
    /// The result retains the exact solid field, sent boundary states, recomputed
    /// coolant states and pump evaluation. Pump matching cannot bypass any inner
    /// regional power, wall, reference, coefficient or mass-flow check.
    ///
    /// Supply head/flow applicability and inlet viscosity are enforced without
    /// fallback. The pump hydraulic-power output is reported, not added to heat;
    /// all existing local and network omitted-dissipation allowances still apply.
    /// Callbacks MUST be disposable evaluations from the same accepted history.
    pub fn solve_with_pump_coupled_iqn<'s, R, F>(
        &self, cx: &Cx<'_>, initial_walls_k: &[f64], pump: &'s PumpHeadCurve,
        controls: &PumpCouplingControls, mut solid: F,
    ) -> Result<PumpOperatingPoint<'s, CoupledParallelLiquidSolution<'a, R>>, PumpOperatingError>
    where F: FnMut(&Cx<'_>, &[CorrelationRobinBoundary]) -> Result<LiquidSolidResponse<R>, AirflowError>,
    {
        controls.operating_point.admit()?;
        controls.network.work.admit()?;
        if controls.max_total_solid_solves == 0 {
            return Err(PumpOperatingError::Invalid("pump solid callback budget must be positive"));
        }
        self.admit_pump_inlet(cx, pump)?;
        let mut flows = 0usize;
        let mut solids = 0usize;
        match_pump(cx, pump, controls.operating_point, |drive| {
            let remaining = controls.network.work.max_flow_evaluations.saturating_sub(flows);
            if remaining < 2 {
                return Err(ParallelPressureError::FlowBudget {
                    performed: flows, maximum: controls.network.work.max_flow_evaluations,
                }.into());
            }
            if solids == controls.max_total_solid_solves {
                return Err(PumpOperatingError::SolidBudget { performed: solids, maximum: controls.max_total_solid_solves });
            }
            let mut inner = controls.network;
            inner.work.max_flow_evaluations = remaining;
            let mut exhausted = false;
            let result = self.solve_coupled_iqn(cx, initial_walls_k, drive, &inner, |cx, rows| {
                if solids == controls.max_total_solid_solves {
                    exhausted = true;
                    return Err(AirflowError::InvalidConjugateInput {
                        field: "pump total solid callback budget exhausted", value_bits: 0,
                    });
                }
                solids += 1;
                solid(cx, rows)
            });
            if exhausted {
                return Err(PumpOperatingError::SolidBudget { performed: solids, maximum: controls.max_total_solid_solves });
            }
            let state = result?;
            flows = flows.checked_add(state.total_flow_evaluations)
                .ok_or(PumpOperatingError::Invalid("pump flow-evaluation count overflow"))?;
            Ok(state)
        })
    }
}
