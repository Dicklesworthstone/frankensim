//! Hold declared .fsim boundary intervals through the original scenario owner.
//! Geometry, row kinds and shared Dirichlet precedence remain those of the
//! already-lowered native partition. No second boundary or time integrator.

use fs_conduction::bc::history::{ThermalBoundaryHistory, ThermalRegionHistory};
use fs_conduction::{HTC_DIMS, ScalarField, ThermalBc, ThermalBoundary};
use fs_qty::QtyAny;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use fs_project::TransientBoundaryTarget;

use super::*;

pub(super) struct BoundaryHistory<'m> {
    native: ThermalBoundaryHistory<'m>,
    pub(super) receipt: String,
}

impl<'m> BoundaryHistory<'m> {
    pub(super) fn bind(
        cx: &Cx<'_>, spec: &ProjectSpec, policy: &ConductionTransient,
        problem: ConductionProblem<'m>, deadline: Option<(Instant, f64)>,
    ) -> Result<Option<Self>, SolveRefusal> {
        if policy.boundary_schedules.is_empty() { return Ok(None); }
        // The original project time-grid owner validates dimensions, finite
        // ordered intervals and exact horizon coverage before allocating tables.
        policy.coarse_step_ends_s().map_err(|error| bad(error.what))?;
        let envelope = spec.envelope.as_ref()
            .ok_or_else(|| bad("thermal boundary histories require the declared operating envelope"))?;
        // All Robin references below are explicit signals, so no ambient or
        // gravity fallback participates in this thermal-only lowering.
        let mut environment = Environment::earth_lab();
        environment.ambient_temperature = envelope.ambient_lo;
        environment.ambient_pressure = envelope.pressure;
        let mut histories = Vec::with_capacity(policy.boundary_schedules.len());
        let mut receipts = Vec::with_capacity(policy.boundary_schedules.len());
        for schedule in &policy.boundary_schedules {
            poll(cx, deadline)?;
            let mut starts = Vec::with_capacity(schedule.steps.len());
            let mut values = Vec::with_capacity(schedule.steps.len());
            let mut intervals = Vec::with_capacity(schedule.steps.len());
            let mut start = 0.0;
            for step in &schedule.steps {
                starts.push(start);
                values.push(step.value.value);
                intervals.push(format!("{{\"until_s\":{},\"value_si\":{}}}",
                    number(step.until.value)?, number(step.value.value)?));
                start = step.until.value;
            }
            // Hold is right-continuous at each authored switch. Both grids
            // land there and call at(interval_start), never at the future end.
            let signal = TimeSignal::Table {
                times: starts, values, dims: schedule.kind.dims(), interp: Interp::Hold,
            };
            let (kind, value, reference_temperature) = match schedule.kind {
                TransientBoundaryTarget::FixedTemperature =>
                    (BcKind::Dirichlet, BcValue::Signal(signal), None),
                TransientBoundaryTarget::HeatFlux =>
                    (BcKind::Neumann, BcValue::Signal(signal), None),
                TransientBoundaryTarget::ConvectionTemperature => {
                    let index = problem.boundary.region_names().iter()
                        .position(|name| name == &schedule.target)
                        .ok_or_else(|| bad("thermal history target is absent from the retained boundary"))?;
                    let ThermalBc::Robin { htc: ScalarField::Uniform(h), .. } =
                        &problem.boundary.conditions()[index]
                        else { return Err(bad("reference-temperature history requires the original uniform Robin coefficient")); };
                    (BcKind::Robin, BcValue::Uniform(QtyAny::new(*h, HTC_DIMS)), Some(signal))
                }
            };
            histories.push(ThermalRegionHistory {
                row: BoundaryCondition {
                    region: schedule.target.clone(), physics: Physics::Thermal,
                    kind, value: Some(value), compatibility: None, frame: 0,
                },
                reference_temperature,
            });
            receipts.push(format!(
                "{{\"target\":{},\"kind\":{},\"source\":{},\"steps\":[{}]}}",
                json_string(&schedule.target), json_string(schedule.kind.name()),
                json_string(&schedule.source), intervals.join(","),
            ));
        }
        let native = ThermalBoundaryHistory::from_boundary(
            cx, problem.mesh, problem.boundary, &environment, histories, 0.0,
        ).map_err(lower)?;
        poll(cx, deadline)?;
        Ok(Some(Self { native, receipt: format!(
            "{{\"sampling\":\"held-interval-start\",\"partition\":\"original-native-boundary\",\"initial_state\":\"unchanged-declared-physical-history\",\"schedules\":[{}]}}",
            receipts.join(","),
        ) }))
    }

    pub(super) fn at(&self, cx: &Cx<'_>, time_s: f64) -> Result<ThermalBoundary, SolveRefusal> {
        self.native.at(cx, time_s).map_err(lower)
    }
}
