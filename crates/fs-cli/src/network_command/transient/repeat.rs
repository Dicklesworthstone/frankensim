//! Repeated duty cycles use the actual transient map, not repeated cold starts.
//! Only the last cycle keeps full output; bounded scalar summaries retain the
//! warm-up peaks and energy of every completed cycle for design decisions.

use super::*;

#[derive(Debug, Clone, Copy)]
pub(super) struct Config {
    cycles: usize,
    max_total_steps: usize,
}

impl Config {
    pub(super) fn parse(value: &J, planned_steps: usize, adaptive: bool) -> Result<Self> {
        object(value, &["cycles", "max_total_steps"], "transient.repeat")?;
        let config = Self {
            cycles: count(get(value,"cycles")?,"repeat.cycles",4096)?,
            max_total_steps: count(get(value,"max_total_steps")?,"repeat.max_total_steps",1_000_000)?,
        };
        let minimum = planned_steps.checked_mul(if adaptive {2} else {1})
            .and_then(|steps|steps.checked_mul(config.cycles))
            .ok_or_else(||budget("repeated-cycle step count overflow"))?;
        if minimum > config.max_total_steps {
            return Err(budget("planned repeated cycles exceed the total accepted-step budget"));
        }
        Ok(config)
    }
}

pub(super) fn simulate(request: &Request, cx: &Cx<'_>, schedule: &Schedule,
    speed_multiplier: f64, config: Config) -> Result<Trajectory>
{
    poll(cx)?;
    let mut field = schedule.initial.clone();
    let mut elapsed = 0.0;
    let mut peak = f64::NEG_INFINITY;
    let mut peak_time = 0.0;
    let mut first_violation = None;
    let mut steps = 0_usize;
    let mut work = 0_usize;
    let mut input = 0.0;
    let mut stored = 0.0;
    let mut exhaust = 0.0;
    let mut summaries = Vec::new();
    for cycle_index in 0..config.cycles {
        poll(cx)?;
        let remaining = config.max_total_steps.checked_sub(steps)
            .filter(|&left|left>0).ok_or_else(||budget("repeated-cycle accepted-step budget exhausted"))?;
        let cycle = simulate_cycle(request,cx,schedule,speed_multiplier,&field,remaining)?;
        let residual = field_residual(cx,&field,&cycle.final_temperature)?;
        let end = finite(elapsed + cycle.duration_s)?;
        if end <= elapsed { return Err(bad("elapsed cycle time is no longer representable")); }
        let global_peak_time = finite(elapsed + cycle.trajectory.peak_time_s)?;
        if cycle.trajectory.peak_k > peak {
            peak = cycle.trajectory.peak_k;
            peak_time = global_peak_time;
        }
        if first_violation.is_none() {
            first_violation = cycle.first_violation_s.map(|time|finite(elapsed+time)).transpose()?;
        }
        steps = steps.checked_add(cycle.trajectory.steps).ok_or_else(||budget("repeated step count overflow"))?;
        work = work.checked_add(cycle.trajectory.solid_solves).ok_or_else(||budget("repeated solid-work count overflow"))?;
        input = finite(input+cycle.input_j)?;
        stored = finite(stored+cycle.stored_j)?;
        exhaust = finite(exhaust+cycle.exhaust_j)?;
        let energy_residual = finite(stored-input+exhaust)?;
        if energy_residual.abs()>finite(request.limits.heat*end)? {
            return Err(producer("repeated-cycle cumulative energy gate failed"));
        }
        summaries.push(format!("{{\"cycle\":{},\"start_time_s\":{},\"end_time_s\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"start_to_end_field_residual_k\":{},\"stored_energy_change_j\":{},\"input_energy_j\":{},\"air_energy_gain_j\":{},\"accepted_steps\":{},\"solid_solves\":{}}}",
            cycle_index+1,num(elapsed)?,num(end)?,num(cycle.trajectory.peak_k)?,num(global_peak_time)?,
            num(residual)?,num(cycle.stored_j)?,num(cycle.input_j)?,num(cycle.exhaust_j)?,
            cycle.trajectory.steps,cycle.trajectory.solid_solves));
        if cycle_index+1==config.cycles {
            let prefix = cycle.trajectory.output.strip_suffix("}\n").ok_or_else(||bad("internal cycle result framing"))?;
            let output = format!("{prefix},\"repeated_cycles\":{{\"status\":\"fixed-count-complete\",\"cycles_completed\":{},\"cycle_duration_s\":{},\"elapsed_time_s\":{},\"last_cycle_start_time_s\":{},\"total_accepted_steps\":{},\"total_solid_solves\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"temperature_limit_k\":{},\"first_sampled_violation_s\":{},\"input_energy_j\":{},\"stored_energy_change_j\":{},\"air_energy_gain_j\":{},\"energy_residual_j\":{},\"cycles\":[{}],\"scope\":\"all cycles inherit the prior accepted nodal field; peaks include initial state and every accepted sample; transient contains only the final cycle in local time; no settled-periodic, future-cycle or continuous-time bound\"}}}}\n",
                config.cycles,num(cycle.duration_s)?,num(end)?,num(elapsed)?,steps,work,num(peak)?,num(peak_time)?,
                optional(schedule.limit)?,optional(first_violation)?,num(input)?,num(stored)?,num(exhaust)?,num(energy_residual)?,summaries.join(","));
            poll(cx)?;
            return Ok(Trajectory {output,peak_k:peak,peak_time_s:peak_time,solid_solves:work,steps});
        }
        field = cycle.final_temperature;
        elapsed = end;
    }
    Err(bad("repeat.cycles must be positive"))
}

fn field_residual(cx: &Cx<'_>, start: &[f64], end: &[f64]) -> Result<f64> {
    if start.is_empty() || start.len()!=end.len() { return Err(bad("cycle boundary fields have incompatible lengths")); }
    let mut residual = 0.0_f64;
    for (index,(&a,&b)) in start.iter().zip(end).enumerate() {
        if index%512==0 { poll(cx)?; }
        finite(a)?; finite(b)?;
        residual = residual.max(finite(b-a)?.abs());
    }
    Ok(residual)
}

#[cfg(test)]
mod tests;
