//! Physical reservoir and fan parameters of the complete trajectory transpose.
//! All contractions use the SAME accepted nodal-load pullback, with prescribed
//! residual rows masked and dt included exactly once by the shared core.

use fs_airflow::conjugate::{AirPath, goal::maximum::pullback_inlet_temperatures};
use fs_conduction::{ConductionProblem, ScalarField, ThermalBc};
use fs_project::ProjectSpec;

use super::super::{json_string, number};
use super::{
    Bound, Cx, Instant, SolveRefusal, air_error, airflow, bad, buoyancy, finite, poll, zeros,
};
use crate::solve::nominal_adjoint::fan_speed;

struct Reservoir {
    target: String,
    ambient_k: f64,
    derivative: f64,
}

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

/// Constant-through-trajectory parameters, accumulated only over the selected
/// prefix. Unsupported fan coordinates retain their identity and reason.
pub(in crate::solve::transient::adjoint::enthalpy) struct Controls {
    natural: Vec<Reservoir>,
    inlets: Vec<Inlet>,
    fans: Vec<Fan>,
    paths: Vec<AirPath>,
    steps: usize,
}

/// Bounded fan-law and aggregate-control scratch in addition to the already
/// admitted nodal factors/trace arrays. The shared hydraulic owner admits at
/// most 8192 curve points; this allowance is not a measured process RSS.
pub(in crate::solve::transient) fn storage_bytes(
    spec: &ProjectSpec,
    ports: usize,
) -> Result<usize, SolveRefusal> {
    if ports == 0 {
        return Ok(0);
    }
    let banks = spec
        .cooling
        .as_ref()
        .and_then(|c| c.fan_system.as_ref())
        .map_or(0, |system| system.banks.len());
    if ports > 64 || banks > 64 {
        return Err(bad(
            "native enthalpy boundary controls admit at most 64 ports and fan banks",
        ));
    }
    let fan_scratch = if banks == 0 { 0 } else { 8192usize * 512 };
    ports
        .checked_add(banks)
        .and_then(|v| v.checked_mul(8192))
        .and_then(|v| v.checked_add(fan_scratch))
        .ok_or_else(|| bad("enthalpy boundary-control workspace allowance overflows"))
}

impl Controls {
    pub(in crate::solve::transient::adjoint::enthalpy) fn new(
        cx: &Cx<'_>,
        spec: &ProjectSpec,
        natural: &[buoyancy::NaturalLaw],
        air: Option<&airflow::ConjugatePath>,
        deadline: Option<(Instant, f64)>,
    ) -> Result<Option<Self>, SolveRefusal> {
        poll(cx, deadline)?;
        let count = super::ports(natural, air)?;
        if count == 0 {
            return Ok(None);
        }
        storage_bytes(spec, count)?;
        let paths = air.map(|path| path.air_paths()).unwrap_or_default();
        let branch_names = air.map(|path| path.branch_names()).unwrap_or_default();
        if paths.len() != branch_names.len() {
            return Err(bad(
                "enthalpy inlet controls lost their physical branch identities",
            ));
        }
        let inlets = paths
            .iter()
            .zip(branch_names)
            .map(|(path, branch)| Inlet {
                branch: branch.to_string(),
                temperature_k: path.inlet_temperature_k(),
                derivative: 0.0,
            })
            .collect();
        let natural = natural
            .iter()
            .map(|law| Reservoir {
                target: law.target.clone(),
                ambient_k: law.ambient_k,
                derivative: 0.0,
            })
            .collect();
        let fans = spec
            .cooling
            .as_ref()
            .and_then(|c| c.fan_system.as_ref())
            .map(|system| {
                system
                    .banks
                    .iter()
                    .map(|bank| Fan {
                        bank: bank.bank_id.clone(),
                        speed_ratio: bank.speed_ratio,
                        derivative: Ok(0.0),
                    })
                    .collect()
            })
            .unwrap_or_default();
        Ok(Some(Self {
            natural,
            inlets,
            fans,
            paths,
            steps: 0,
        }))
    }

    /// Add one accepted endpoint's parameter contractions after its complete
    /// physical transpose has passed the true-residual gate. No solve occurs.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(in crate::solve::transient::adjoint::enthalpy) fn add(
        &mut self,
        cx: &Cx<'_>,
        spec: &ProjectSpec,
        problem: ConductionProblem<'_>,
        temperature: &[f64],
        bound: &Bound,
        load_weights: &[f64],
        deadline: Option<(Instant, f64)>,
    ) -> Result<(), SolveRefusal> {
        poll(cx, deadline)?;
        let n = problem.mesh.vertex_count();
        let m = bound.names.len();
        let air_count: usize = self.paths.iter().map(|path| path.segments().len()).sum();
        if temperature.len() != n
            || load_weights.len() != n
            || m != self.natural.len() + air_count
            || bound.ambient_slopes.len() != m
        {
            return Err(bad(
                "enthalpy boundary controls lost their complete endpoint/port layout",
            ));
        }
        for (vertex, &weight) in load_weights.iter().enumerate() {
            if vertex % 256 == 0 {
                poll(cx, deadline)?;
            }
            finite(weight)?;
            finite(temperature[vertex])?;
        }
        let mut htc = zeros(m)?;
        let mut reference_k = zeros(m)?;
        for (index, name) in bound.names.iter().enumerate() {
            poll(cx, deadline)?;
            let region = problem
                .boundary
                .region_names()
                .iter()
                .position(|found| found == name)
                .ok_or_else(|| bad("enthalpy boundary control has no original physical trace"))?;
            let ThermalBc::Robin {
                htc: ScalarField::Uniform(h),
                t_ref: ScalarField::Uniform(r),
            } = &problem.boundary.conditions()[region]
            else {
                return Err(bad(
                    "enthalpy parameter controls require original uniform convection",
                ));
            };
            htc[index] = finite(*h)?;
            reference_k[index] = finite(*r)?;
            if *h <= 0.0 {
                return Err(bad("enthalpy parameter control has nonpositive convection"));
            }
        }
        let (references, log_htc) = super::super::robin_controls::partials(
            cx, problem, temperature, &bound.names, load_weights, deadline,
        )?;
        for (index, control) in self.natural.iter_mut().enumerate() {
            if bound.names[index] != control.target || reference_k[index] != control.ambient_k {
                return Err(bad(
                    "enthalpy natural reservoir control changed its declaration",
                ));
            }
            // Explicit reference plus the coefficient's ambient partial at
            // fixed Tw. Its indirect Tw dependence is already in the dual.
            let derivative = finite(
                references[index]
                    + log_htc[index] * finite(bound.ambient_slopes[index] / htc[index])?,
            )?;
            control.derivative = finite(control.derivative + derivative)?;
        }
        if !self.paths.is_empty() {
            let offset = self.natural.len();
            let names: Vec<&str> = bound.names[offset..].iter().map(String::as_str).collect();
            let bars =
                pullback_inlet_temperatures(cx, &self.paths, &names, &references[offset..], 64)
                    .map_err(air_error)?;
            if bars.len() != self.inlets.len() {
                return Err(bad("enthalpy inlet pullback changed its branch count"));
            }
            for (control, bar) in self.inlets.iter_mut().zip(bars) {
                control.derivative = finite(control.derivative + bar)?;
            }
        }
        let fans = fan_speed::coupled_gradient(
            cx,
            spec,
            problem.mesh,
            problem.boundary,
            temperature,
            &self.paths,
            load_weights,
        )?;
        if fans.len() != self.fans.len() {
            return Err(bad("enthalpy fan-speed pullback lost a bank"));
        }
        for (ordinal, (total, step)) in self.fans.iter_mut().zip(fans).enumerate() {
            if step.ordinal != ordinal || step.bank_id != total.bank {
                return Err(bad("enthalpy fan-speed control changed its bank identity"));
            }
            match (&mut total.derivative, step.derivative) {
                (Ok(sum), Ok(value)) => *sum = finite(*sum + value)?,
                (current @ Ok(_), Err(reason)) => *current = Err(reason),
                (Err(_), _) => {}
            }
        }
        self.steps = self
            .steps
            .checked_add(1)
            .ok_or_else(|| bad("enthalpy boundary-control step count overflows"))?;
        poll(cx, deadline)
    }

    pub(in crate::solve::transient::adjoint::enthalpy) fn fragment(
        &self,
        expected_steps: usize,
    ) -> Result<String, SolveRefusal> {
        if self.steps != expected_steps {
            return Err(bad(
                "enthalpy boundary controls do not cover the selected reverse prefix",
            ));
        }
        let natural = self.natural.iter().map(|control| Ok(format!(
            "{{\"target\":{},\"ambient_temperature_k\":{},\"dtemperature_dambient_temperature_k_per_k\":{}}}",
            json_string(&control.target), number(control.ambient_k)?, number(control.derivative)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        let inlets = self.inlets.iter().map(|control| Ok(format!(
            "{{\"branch\":{},\"inlet_temperature_k\":{},\"dtemperature_dinlet_temperature_k_per_k\":{}}}",
            json_string(&control.branch), number(control.temperature_k)?, number(control.derivative)?,
        ))).collect::<Result<Vec<_>, SolveRefusal>>()?;
        let fans = self
            .fans
            .iter()
            .map(|control| {
                let state = match &control.derivative {
                    Ok(value) => format!(
                        "\"status\":\"available\",\"dtemperature_dspeed_ratio_k\":{}",
                        number(*value)?
                    ),
                    Err(reason) => format!(
                        "\"status\":\"unsupported\",\"reason\":{}",
                        json_string(reason)
                    ),
                };
                Ok(format!(
                    "{{\"bank\":{},\"speed_ratio\":{},{state}}}",
                    json_string(&control.bank),
                    number(control.speed_ratio)?
                ))
            })
            .collect::<Result<Vec<_>, SolveRefusal>>()?;
        Ok(format!(
            ",\"parameter_scope\":\"constant-through-selected-trajectory-physical-reservoirs-and-independent-absolute-fan-bank-speeds\",\"control_steps\":{},\"natural_reservoirs\":[{}],\"air_inlets\":[{}],\"fan_banks\":[{}]",
            self.steps,
            natural.join(","),
            inlets.join(","),
            fans.join(","),
        ))
    }
}
