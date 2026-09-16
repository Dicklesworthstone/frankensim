//! Size cooling or workload power against a complete transient, not steady state.
//! Only actually evaluated sampled peaks decide feasibility. Adaptive samples
//! may differ between candidates; neither local estimates nor a tight parameter
//! bracket certify a continuous-time peak or global optimality.

use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Control { FanSpeed, WorkloadPower }

#[derive(Debug)]
pub(super) struct Config {
    control: Control,
    minimum: f64,
    maximum: f64,
    multiplier_tolerance: f64,
    temperature_tolerance_k: f64,
    max_evaluations: usize,
}

impl Config {
    pub(super) fn parse(value: &J) -> Result<Self> {
        object(value, &["min_speed_multiplier", "max_speed_multiplier", "speed_multiplier_tolerance",
            "temperature_tolerance_k", "max_evaluations"], "transient.fan_speed_design")?;
        let result = Self {
            control: Control::FanSpeed,
            minimum: positive(get(value,"min_speed_multiplier")?,"min_speed_multiplier")?,
            maximum: positive(get(value,"max_speed_multiplier")?,"max_speed_multiplier")?,
            multiplier_tolerance: positive(get(value,"speed_multiplier_tolerance")?,"speed_multiplier_tolerance")?,
            temperature_tolerance_k: positive(get(value,"temperature_tolerance_k")?,"fan design temperature_tolerance_k")?,
            max_evaluations: count(get(value,"max_evaluations")?,"transient fan max_evaluations",256)?,
        };
        if result.minimum >= result.maximum { return Err(bad("transient fan multiplier bounds must be strictly ordered")); }
        Ok(result)
    }

    pub(super) fn parse_power(value: &J) -> Result<Self> {
        object(value, &["min_power_multiplier","max_power_multiplier","power_multiplier_tolerance",
            "temperature_tolerance_k","max_evaluations"], "transient.power_design")?;
        let result = Self {
            control: Control::WorkloadPower,
            minimum: number(get(value,"min_power_multiplier")?,"min_power_multiplier")?,
            maximum: positive(get(value,"max_power_multiplier")?,"max_power_multiplier")?,
            multiplier_tolerance: positive(get(value,"power_multiplier_tolerance")?,"power_multiplier_tolerance")?,
            temperature_tolerance_k: positive(get(value,"temperature_tolerance_k")?,"power design temperature_tolerance_k")?,
            max_evaluations: count(get(value,"max_evaluations")?,"transient power max_evaluations",256)?,
        };
        if result.minimum < 0.0 || result.minimum >= result.maximum {
            return Err(bad("transient power multiplier bounds require 0 <= minimum < maximum"));
        }
        Ok(result)
    }

    /// Validate the entire schedule's endpoint speeds before running any trial.
    /// A valid speed domain does not waive flow/correlation admission at runtime.
    pub(super) fn validate(&self, schedule: &Schedule, fan: Option<&fan_drive::FanDrive>) -> Result<()> {
        if schedule.adjoint.is_some() { return Err(bad("transient adjoints cannot be nested inside sizing")); }
        if schedule.limit.is_none() { return Err(bad("transient sizing requires transient.temperature_limit_k")); }
        if self.control == Control::WorkloadPower {
            for interval in &schedule.intervals { scaled_workload(&interval.workload,self.maximum)?; }
            return Ok(());
        }
        let fan = fan.ok_or_else(||bad("transient fan sizing requires hydraulics.fan"))?;
        for interval in &schedule.intervals {
            let base = interval.speed.ok_or_else(||bad("transient fan sizing requires every interval's base speed"))?;
            for multiplier in [self.minimum,self.maximum] {
                fan.bank(finite(base*multiplier)?)?;
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
struct Trial {
    multiplier: f64,
    peak_k: f64,
    peak_time_s: f64,
    steps: usize,
    solid_solves: usize,
}

struct Selected {
    passing: Trajectory,
    multiplier: f64,
    failed: Option<Trial>,
    width: f64,
    history: Vec<Trial>,
}

fn evaluate(cx: &Cx<'_>, config: &Config, multiplier: f64, history: &mut Vec<Trial>,
    run: &mut impl FnMut(f64) -> Result<Trajectory>) -> Result<Trajectory>
{
    poll(cx)?;
    if history.len() >= config.max_evaluations {
        return Err(budget("transient design evaluation budget exhausted; no partial design published"));
    }
    let evaluated = run(multiplier)?;
    if !evaluated.peak_k.is_finite() || !evaluated.peak_time_s.is_finite() {
        return Err(producer("nonfinite transient design objective"));
    }
    history.push(Trial {multiplier,peak_k:evaluated.peak_k,peak_time_s:evaluated.peak_time_s,
        steps:evaluated.steps,solid_solves:evaluated.solid_solves});
    poll(cx)?;
    Ok(evaluated)
}

/// The callable is the actual complete simulation in production. A scalar-only
/// test seam exercises brackets/refusals without inventing a second simulator.
fn search(cx: &Cx<'_>, config: &Config, limit: f64,
    mut run: impl FnMut(f64) -> Result<Trajectory>) -> Result<Selected>
{
    let mut history = Vec::new();
    let power = config.control == Control::WorkloadPower;
    let preferred = if power {config.maximum}else{config.minimum};
    let fallback = if power {config.minimum}else{config.maximum};
    let first = evaluate(cx,config,preferred,&mut history,&mut run)?;
    if first.peak_k <= limit {
        return Ok(Selected {passing:first,multiplier:preferred,failed:None,width:0.0,history});
    }
    let mut failed = history[0];
    drop(first);
    let mut passing = evaluate(cx,config,fallback,&mut history,&mut run)?;
    if passing.peak_k > limit {
        return Err(Failure {code:if power {"cooling-network-transient-power-bracket"}else{"cooling-network-transient-fan-bracket"},message:format!(
            "neither schedule endpoint passes the sampled peak limit {limit} K: preferred {} K, fallback {} K; no global infeasibility proof",
            failed.peak_k,passing.peak_k)});
    }
    let mut chosen = fallback;
    loop {
        poll(cx)?;
        let width = finite(chosen-failed.multiplier)?.abs();
        let slack = finite(limit-passing.peak_k)?;
        if width <= config.multiplier_tolerance && slack <= config.temperature_tolerance_k {
            return Ok(Selected {passing,multiplier:chosen,failed:Some(failed),width,history});
        }
        let candidate = 0.5*failed.multiplier + 0.5*chosen;
        if !(candidate > chosen.min(failed.multiplier) && candidate < chosen.max(failed.multiplier)) {
            return Err(Failure {code:if power {"cooling-network-transient-power-resolution"}else{"cooling-network-transient-fan-resolution"},
                message:"no representable multiplier can meet both design tolerances".into()});
        }
        let evaluated = evaluate(cx,config,candidate,&mut history,&mut run)?;
        if evaluated.peak_k <= limit { chosen=candidate; passing=evaluated; }
        else { failed=*history.last().ok_or_else(||bad("missing transient design trial"))?; }
    }
}

pub(super) fn solve(request: &Request, cx: &Cx<'_>, schedule: &Schedule, config: &Config) -> Result<String> {
    poll(cx)?;
    config.validate(schedule,request.fan.as_ref())?;
    let limit = schedule.limit.ok_or_else(||bad("missing transient design temperature limit"))?;
    let power = config.control == Control::WorkloadPower;
    let selected = search(cx,config,limit,|multiplier| {
        if power {
            let scaled = power_schedule(cx,schedule,multiplier)?;
            simulate(request,cx,&scaled,1.0)
        } else { simulate(request,cx,schedule,multiplier) }
    })?;
    let multiplier_key = if power {"power_multiplier"}else{"speed_multiplier"};
    let failed = selected.failed.map(|trial|render_trial(&trial,multiplier_key)).transpose()?.unwrap_or_else(||"null".into());
    let trials = selected.history.iter().map(|t|render_trial(t,multiplier_key)).collect::<Result<Vec<_>>>()?.join(",");
    let mut total_steps=0_usize;
    let mut total_solves=0_usize;
    for trial in &selected.history {
        total_steps=total_steps.checked_add(trial.steps).ok_or_else(||budget("design step count overflow"))?;
        total_solves=total_solves.checked_add(trial.solid_solves).ok_or_else(||budget("design solid work count overflow"))?;
    }
    let applied = schedule.intervals.iter().enumerate().map(|(index,interval)| {
        poll(cx)?;
        if power {
            let workload=scaled_workload(&interval.workload,selected.multiplier)?;
            Ok(format!("{{\"interval\":{index},{},\"fan_speed_ratio\":{}}}",workload.render()?,optional(interval.speed)?))
        } else {
            let base=interval.speed.ok_or_else(||bad("missing schedule speed"))?;
            Ok(format!("{{\"interval\":{index},\"base_speed_ratio\":{},\"selected_speed_ratio\":{}}}",
                num(base)?,num(finite(base*selected.multiplier)?)?))
        }
    }).collect::<Result<Vec<_>>>()?.join(",");
    let name=if power {"transient_power_design"}else{"transient_fan_speed_design"};
    let failed_key=if power {"failed_upper"}else{"failed_lower"};
    let status=if selected.failed.is_some(){"target-bracketed"}
        else if power {"maximum-feasible"}else{"minimum-feasible"};
    let scope=if power {
        "passing evaluated sampled workload at fixed fan schedule; every source in every interval is multiplied, with fixed footprints and durations; repeated-cycle warm-up peaks remain part of feasibility; no continuous-time compliance, globally maximal workload, electrical-power or transient-adjoint claim"
    } else {
        "passing evaluated sampled trajectory; each candidate restarts the same initial field and rescales the whole base fan schedule; repeated cycles carry heat and all their samples enter feasibility; no continuous-time compliance, global minimum-speed, electrical-power or fan-speed-adjoint claim"
    };
    let prefix=selected.passing.output.strip_suffix("}\n").ok_or_else(||bad("internal trajectory framing"))?;
    let output=format!("{prefix},\"{name}\":{{\"selected_{multiplier_key}\":{},\"status\":{},\"temperature_limit_k\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"{failed_key}\":{},\"multiplier_bracket_width\":{},\"evaluations\":{},\"total_accepted_steps\":{},\"total_solid_solves\":{},\"schedule\":[{}],\"history\":[{}],\"scope\":{},\"time_comparison\":\"adaptive samples may differ between candidates; local estimates and parameter brackets are not trajectory error bounds\"}}}}\n",
        num(selected.multiplier)?,quote(status),num(limit)?,num(selected.passing.peak_k)?,num(selected.passing.peak_time_s)?,
        failed,num(selected.width)?,selected.history.len(),total_steps,total_solves,applied,trials,quote(scope));
    poll(cx)?;
    Ok(output)
}

fn scaled_workload(workload: &Workload, multiplier: f64) -> Result<Workload> {
    if !multiplier.is_finite() || multiplier < 0.0 { return Err(bad("power multiplier must be finite and nonnegative")); }
    match workload {
        Workload::Scale(base) => Ok(Workload::Scale(finite(base*multiplier)?)),
        Workload::Components(powers) => powers.iter().map(|(name,watts)|
            Ok((name.clone(),finite(watts*multiplier)?))).collect::<Result<BTreeMap<_,_>>>().map(Workload::Components),
    }
}

/// Preserve footprints and let the original PowerMap construct the actual
/// source. Initial temperatures and fan speeds are never power-scaled.
fn power_schedule(cx: &Cx<'_>, schedule: &Schedule, multiplier: f64) -> Result<Schedule> {
    let intervals=schedule.intervals.iter().map(|interval| {
        poll(cx)?;
        Ok(Interval {duration:interval.duration,steps:interval.steps,speed:interval.speed,
            workload:scaled_workload(&interval.workload,multiplier)?})
    }).collect::<Result<Vec<_>>>()?;
    Ok(Schedule {initial:schedule.initial.clone(),capacities:schedule.capacities.clone(),intervals,
        limit:schedule.limit,total_steps:schedule.total_steps,max_step_s:schedule.max_step_s,max_steps:schedule.max_steps,
        adaptive:schedule.adaptive,nonlinear:schedule.nonlinear,adjoint:schedule.adjoint,
        fan_speed_design:None,power_design:None,repeat:schedule.repeat})
}

fn render_trial(trial: &Trial, multiplier_key: &str) -> Result<String> {
    Ok(format!("{{\"{multiplier_key}\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"accepted_steps\":{},\"solid_solves\":{}}}",
        num(trial.multiplier)?,num(trial.peak_k)?,num(trial.peak_time_s)?,trial.steps,trial.solid_solves))
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod power_tests;
