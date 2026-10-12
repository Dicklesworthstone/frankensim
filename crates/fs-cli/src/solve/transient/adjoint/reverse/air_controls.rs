//! Physical inlet and independent fan-bank controls of a capacity trajectory.
//!
//! The accepted complete temperature dual uses watt residual coordinates.
//! Neither the original Robin bars nor its nodal load multiplier receive dt.
//! Natural and radiative reservoirs retain their separate physical owners.

use std::collections::BTreeSet;
use fs_airflow::conjugate::{
    AirPath, goal::{CoupledGoalError, maximum::pullback_inlet_temperatures},
};
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc, adjoint::RobinGradient};
use fs_project::ProjectSpec;
use crate::solve::{conjugate::ConjugatePath, nominal_adjoint::fan_speed};
use super::{Cx, Instant, SolveRefusal, bad, finite, poll};
use super::super::super::{json_string, number};

struct Inlet {
    branch: String,
    temperature_k: f64,
    derivative: f64,
}

struct Fan {
    bank: String,
    speed_ratio: f64,
    derivative: Result<f64, String>,
}

/// Constant-through-trajectory controls, accumulated over the selected prefix.
/// An initial sampled maximum has zero future-control influence and never
/// needs to query a post-peak hydraulic differential.
pub(super) struct Controls {
    inlets: Vec<Inlet>,
    fans: Vec<Fan>,
    paths: Vec<AirPath>,
    steps: usize,
}

/// Bounded port/control and hydraulic-law scratch. The existing hydraulic
/// owner admits at most 8192 member-curve points; no nodal dense matrix or
/// per-time-step control arrays are retained. This is an allowance, not RSS.
pub(super) fn storage_bytes(
    spec: &ProjectSpec, air_ports: usize,
) -> Result<usize, SolveRefusal> {
    if air_ports == 0 { return Ok(0); }
    let banks = spec.cooling.as_ref().and_then(|c| c.fan_system.as_ref())
        .map_or(0, |system| system.banks.len());
    if air_ports > 64 || banks > 64 {
        return Err(bad("capacity airflow controls admit at most 64 ports and fan banks"));
    }
    let fan_scratch = if banks == 0 { 0 } else { 8192usize * 512 };
    // Member lowering clones declared provenance as well as curve points.
    // Citations have no fixed length cap, so a per-bank constant cannot admit
    // those live input/output/temporary copies. Single-bank affinity does not
    // reconstruct member curves and needs no such copies.
    let mut provenance_bytes = 0usize;
    if let Some(system) = spec.cooling.as_ref().and_then(|c| c.fan_system.as_ref())
        .filter(|system| system.banks.len() > 1)
    {
        for bank in &system.banks {
            provenance_bytes = bank.curve.source.len().checked_add(bank.curve.source_id.len())
                .and_then(|bytes| bytes.checked_mul(8))
                .and_then(|bytes| provenance_bytes.checked_add(bytes))
                .ok_or_else(|| bad("capacity fan provenance workspace overflows"))?;
        }
    }
    air_ports.checked_add(banks).and_then(|v| v.checked_mul(8192))
        .and_then(|v| v.checked_add(fan_scratch))
        .and_then(|v| v.checked_add(provenance_bytes))
        .filter(|&bytes| isize::try_from(bytes).is_ok())
        .ok_or_else(|| bad("capacity airflow-control workspace overflows"))
}

impl Controls {
    pub(super) fn new(
        cx: &Cx<'_>, spec: &ProjectSpec, path: &ConjugatePath,
        deadline: Option<(Instant, f64)>,
    ) -> Result<Self, SolveRefusal> {
        poll(cx, deadline)?;
        storage_bytes(spec, path.segments.len())?;
        let paths = path.air_paths();
        let names = path.branch_names();
        if paths.is_empty() || paths.len() != names.len() {
            return Err(bad("capacity inlet controls lost their physical branch identities"));
        }
        let mut seen = BTreeSet::new();
        let mut segments = 0usize;
        let mut inlets = Vec::with_capacity(paths.len());
        for (path, branch) in paths.iter().zip(names) {
            poll(cx, deadline)?;
            if branch.is_empty() || !seen.insert(branch) || path.segments().is_empty() {
                return Err(bad("capacity inlet controls need distinct nonempty branches"));
            }
            segments = segments.checked_add(path.segments().len())
                .ok_or_else(|| bad("capacity inlet control count overflows"))?;
            inlets.push(Inlet {
                branch: branch.to_owned(), temperature_k: finite(path.inlet_temperature_k())?,
                derivative: 0.0,
            });
        }
        if segments != path.segments.len() || segments > 64 {
            return Err(bad("capacity inlet controls omit retained streamwise segments"));
        }
        let fans = spec.cooling.as_ref().and_then(|c| c.fan_system.as_ref())
            .map(|system| system.banks.iter().map(|bank| Fan {
                bank: bank.bank_id.clone(), speed_ratio: bank.speed_ratio, derivative: Ok(0.0),
            }).collect()).unwrap_or_default();
        poll(cx, deadline)?;
        Ok(Self { inlets, fans, paths, steps: 0 })
    }

    /// Add one endpoint after the complete physical transpose passes its true
    /// residual gate. The boundary and reference/log(h) bars are ORIGINAL
    /// convection, after radiation remapping; the field and watt-coordinate
    /// multiplier remain unchanged. The shared fan owner includes transported
    /// air, local h, and the actual independent hydraulic speed partials.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn add(
        &mut self, cx: &Cx<'_>, spec: &ProjectSpec, problem: ConductionProblem<'_>,
        temperature: &[f64], ports: &[&str], gradient: &RobinGradient,
        deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        let nodes = problem.mesh.vertex_count();
        if temperature.len() != nodes || gradient.nodal_load.len() != nodes
            || ports.is_empty() || ports.len() > 64
            || gradient.references.len() != ports.len() || gradient.log_htc.len() != ports.len()
            || !gradient.relative_residual.is_finite() || gradient.relative_residual < 0.0
        {
            return Err(bad("capacity airflow controls lost their complete endpoint gradient"));
        }
        for (vertex, (&temperature, &weight)) in temperature.iter().zip(&gradient.nodal_load).enumerate() {
            if vertex % 256 == 0 { poll(cx, deadline)?; }
            finite(temperature)?; finite(weight)?;
        }
        let mut names = Vec::new();
        let mut reference_bars = Vec::new();
        for path in &self.paths {
            for segment in path.segments() {
                poll(cx, deadline)?;
                let name = segment.region();
                let port = ports.iter().position(|&found| found == name)
                    .ok_or_else(|| bad("capacity inlet control has no selected air trace"))?;
                if ports.iter().filter(|&&found| found == name).count() != 1 {
                    return Err(bad("capacity inlet control repeats a selected air trace"));
                }
                let region = problem.boundary.region_names().iter().position(|found| found == name)
                    .ok_or_else(|| bad("capacity inlet control has no original convection trace"))?;
                let Some(ThermalBc::Robin {
                    htc: ScalarField::Uniform(h), t_ref: ScalarField::Uniform(reference),
                }) = problem.boundary.conditions().get(region)
                else { return Err(bad("capacity inlet controls require original uniform convection")); };
                if h.to_bits() != segment.htc_w_per_m2_k().to_bits()
                    || !reference.is_finite() || *reference <= 0.0
                {
                    return Err(bad("capacity inlet control cannot substitute a combined radiative trace"));
                }
                names.push(name);
                reference_bars.push(finite(gradient.references[port])?);
            }
        }
        let inlets = pullback_inlet_temperatures(cx, &self.paths, &names, &reference_bars, 64)
            .map_err(|error| match error {
                CoupledGoalError::Interrupted => crate::solve::conduction_error(
                    "cli-solve-cancelled", "capacity air-inlet differentiation was cancelled",
                    "resume the accepted pipeline prefix",
                ),
                error => bad(format!("capacity air-inlet derivative refused: {error}")),
            })?;
        if inlets.len() != self.inlets.len() {
            return Err(bad("capacity inlet pullback changed its branch count"));
        }
        for (total, step) in self.inlets.iter_mut().zip(inlets) {
            total.derivative = finite(total.derivative + step)?;
        }
        // NO dt: StepLinearization already exposes dJ/d(nodal load in W).
        let fans = fan_speed::coupled_gradient(
            cx, spec, problem.mesh, problem.boundary, temperature, &self.paths, &gradient.nodal_load,
        )?;
        if fans.len() != self.fans.len() {
            return Err(bad("capacity fan-speed pullback lost a declared bank"));
        }
        for (ordinal, (total, step)) in self.fans.iter_mut().zip(fans).enumerate() {
            if step.ordinal != ordinal || step.bank_id != total.bank {
                return Err(bad("capacity fan-speed pullback changed a bank identity"));
            }
            match (&mut total.derivative, step.derivative) {
                (Ok(sum), Ok(value)) => *sum = finite(*sum + value)?,
                (current @ Ok(_), Err(reason)) => *current = Err(reason),
                (Err(_), _) => {}
            }
        }
        self.steps = self.steps.checked_add(1)
            .ok_or_else(|| bad("capacity airflow-control step count overflows"))?;
        poll(cx, deadline)
    }

    /// Fields for the enclosing convection object; natural controls and
    /// independent radiation surfaces are rendered by their own adapters.
    pub(super) fn fragment(&self, expected_steps: usize) -> Result<String, SolveRefusal> {
        if self.steps != expected_steps {
            return Err(bad("capacity airflow controls do not cover the selected reverse prefix"));
        }
        let inlets = self.inlets.iter().map(|control| Ok(format!(
            "{{\"branch\":{},\"inlet_temperature_k\":{},\"dtemperature_dinlet_temperature_k_per_k\":{}}}",
            json_string(&control.branch), number(control.temperature_k)?, number(control.derivative)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        let fans = self.fans.iter().map(|control| {
            let state = match &control.derivative {
                Ok(value) => format!("\"status\":\"available\",\"dtemperature_dspeed_ratio_k\":{}", number(*value)?),
                Err(reason) => format!("\"status\":\"unsupported\",\"reason\":{}", json_string(reason)),
            };
            Ok(format!(
                "{{\"bank\":{},\"speed_ratio\":{},{state}}}",
                json_string(&control.bank), number(control.speed_ratio)?,
            ))
        }).collect::<Result<Vec<_>, SolveRefusal>>()?;
        Ok(format!(
            ",\"air_control_steps\":{},\"air_parameter_scope\":\"constant-through-selected-trajectory-physical-branch-inlets-and-independent-absolute-fan-bank-speeds\",\"air_inlets\":[{}],\"fan_banks\":[{}]",
            self.steps, inlets.join(","), fans.join(","),
        ))
    }
}
