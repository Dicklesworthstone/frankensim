//! Size a fan schedule against a complete workload trajectory, not steady state.
//! Only actually evaluated sampled peaks decide feasibility. Adaptive samples
//! may differ between candidates; neither local estimates nor a tight parameter
//! bracket certify a continuous-time peak or global minimum fan speed.

use super::*;

#[derive(Debug)]
pub(super) struct Config {
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
            minimum: positive(get(value,"min_speed_multiplier")?,"min_speed_multiplier")?,
            maximum: positive(get(value,"max_speed_multiplier")?,"max_speed_multiplier")?,
            multiplier_tolerance: positive(get(value,"speed_multiplier_tolerance")?,"speed_multiplier_tolerance")?,
            temperature_tolerance_k: positive(get(value,"temperature_tolerance_k")?,"fan design temperature_tolerance_k")?,
            max_evaluations: count(get(value,"max_evaluations")?,"transient fan max_evaluations",256)?,
        };
        if result.minimum >= result.maximum { return Err(bad("transient fan multiplier bounds must be strictly ordered")); }
        Ok(result)
    }

    /// Validate the entire schedule's endpoint speeds before running any trial.
    /// A valid speed domain does not waive flow/correlation admission at runtime.
    pub(super) fn validate(&self, schedule: &Schedule, fan: Option<&fan_drive::FanDrive>) -> Result<()> {
        let fan = fan.ok_or_else(||bad("transient fan sizing requires hydraulics.fan"))?;
        if schedule.limit.is_none() { return Err(bad("transient fan sizing requires transient.temperature_limit_k")); }
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
    failed_lower: Option<Trial>,
    width: f64,
    history: Vec<Trial>,
}

fn evaluate(cx: &Cx<'_>, config: &Config, multiplier: f64, history: &mut Vec<Trial>,
    run: &mut impl FnMut(f64) -> Result<Trajectory>) -> Result<Trajectory>
{
    poll(cx)?;
    if history.len() >= config.max_evaluations {
        return Err(budget("transient fan design evaluation budget exhausted; no partial design published"));
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
    let lower = evaluate(cx,config,config.minimum,&mut history,&mut run)?;
    if lower.peak_k <= limit {
        return Ok(Selected {passing:lower,multiplier:config.minimum,failed_lower:None,width:0.0,history});
    }
    let mut failed = history[0];
    drop(lower);
    let mut passing = evaluate(cx,config,config.maximum,&mut history,&mut run)?;
    if passing.peak_k > limit {
        return Err(Failure {code:"cooling-network-transient-fan-bracket",message:format!(
            "neither schedule endpoint passes the sampled peak limit {limit} K: lower {} K, upper {} K; no global infeasibility proof",
            failed.peak_k,passing.peak_k)});
    }
    let mut high = config.maximum;
    loop {
        poll(cx)?;
        let width = finite(high-failed.multiplier)?;
        let slack = finite(limit-passing.peak_k)?;
        if width <= config.multiplier_tolerance && slack <= config.temperature_tolerance_k {
            return Ok(Selected {passing,multiplier:high,failed_lower:Some(failed),width,history});
        }
        let candidate = 0.5*failed.multiplier + 0.5*high;
        if !(candidate > failed.multiplier && candidate < high) {
            return Err(Failure {code:"cooling-network-transient-fan-resolution",
                message:"no representable speed multiplier can meet both design tolerances".into()});
        }
        let evaluated = evaluate(cx,config,candidate,&mut history,&mut run)?;
        if evaluated.peak_k <= limit { high=candidate; passing=evaluated; }
        else { failed=*history.last().ok_or_else(||bad("missing transient design trial"))?; }
    }
}

pub(super) fn solve(request: &Request, cx: &Cx<'_>, schedule: &Schedule, config: &Config) -> Result<String> {
    poll(cx)?;
    config.validate(schedule,request.fan.as_ref())?;
    let limit = schedule.limit.ok_or_else(||bad("missing transient design temperature limit"))?;
    let selected = search(cx,config,limit,|multiplier|simulate(request,cx,schedule,multiplier))?;
    let failed = selected.failed_lower.map(|trial|render_trial(&trial)).transpose()?.unwrap_or_else(||"null".into());
    let trials = selected.history.iter().map(render_trial).collect::<Result<Vec<_>>>()?.join(",");
    let mut total_steps=0_usize;
    let mut total_solves=0_usize;
    for trial in &selected.history {
        total_steps=total_steps.checked_add(trial.steps).ok_or_else(||budget("design step count overflow"))?;
        total_solves=total_solves.checked_add(trial.solid_solves).ok_or_else(||budget("design solid work count overflow"))?;
    }
    let speeds = schedule.intervals.iter().enumerate().map(|(index,interval)| {
        let base=interval.speed.ok_or_else(||bad("missing schedule speed"))?;
        Ok(format!("{{\"interval\":{index},\"base_speed_ratio\":{},\"selected_speed_ratio\":{}}}",
            num(base)?,num(finite(base*selected.multiplier)?)?))
    }).collect::<Result<Vec<_>>>()?.join(",");
    let prefix=selected.passing.output.strip_suffix("}\n").ok_or_else(||bad("internal trajectory framing"))?;
    let output=format!("{prefix},\"transient_fan_speed_design\":{{\"selected_speed_multiplier\":{},\"status\":{},\"temperature_limit_k\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"failed_lower\":{},\"multiplier_bracket_width\":{},\"evaluations\":{},\"total_accepted_steps\":{},\"total_solid_solves\":{},\"schedule\":[{}],\"history\":[{}],\"scope\":\"passing evaluated sampled trajectory; each candidate restarts the same initial field and rescales the whole base fan schedule; adaptive meshes may differ in time; no continuous-time compliance, global minimum-speed, electrical-power or fan-speed-adjoint claim\"}}}}\n",
        num(selected.multiplier)?,quote(if selected.failed_lower.is_none(){"minimum-feasible"}else{"target-bracketed"}),
        num(limit)?,num(selected.passing.peak_k)?,num(selected.passing.peak_time_s)?,failed,num(selected.width)?,
        selected.history.len(),total_steps,total_solves,speeds,trials);
    poll(cx)?;
    Ok(output)
}

fn render_trial(trial: &Trial) -> Result<String> {
    Ok(format!("{{\"speed_multiplier\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"accepted_steps\":{},\"solid_solves\":{}}}",
        num(trial.multiplier)?,num(trial.peak_k)?,num(trial.peak_time_s)?,trial.steps,trial.solid_solves))
}

#[cfg(test)]
mod tests;
