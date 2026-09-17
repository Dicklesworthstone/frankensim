//! Successive fixed-grid trajectory solves, not local step-error estimates.
//! All old endpoints survive exact integer subdivision. Compare the complete
//! nodal field there AND the all-trajectory sampled peak, which includes new
//! endpoints. Neither a cooled final value nor a steady solve can stand in for
//! the trajectory. Two consecutive passing comparisons are required.
use super::*;

#[derive(Debug, Clone, Copy)]
pub(super) struct Config {
    max_refinements: usize,
    consecutive: usize,
    tolerance_k: f64,
    max_total_steps: usize,
    max_trace_bytes: usize,
    cycles: usize,
    max_trajectory_steps: usize,
}

impl Config {
    pub(super) fn parse(value: &J, schedule: &J, max_steps: usize) -> Result<Self> {
        object(value, &["max_refinements", "consecutive_passes", "temperature_tolerance_k",
            "max_total_steps", "max_trace_bytes"], "transient.time_convergence")?;
        for key in ["adaptive", "adjoint", "fan_speed_design", "power_design"] {
            if schedule.get(key).is_some() {
                return Err(bad("time_convergence requires fixed timesteps without adjoints or nested sizing"));
            }
        }
        let (cycles, max_trajectory_steps) = match schedule.get("repeat") {
            None => (1, max_steps),
            Some(repeat) => {
                if repeat.get("until_periodic").is_some() || repeat.get("fan_controller").is_some() {
                    return Err(bad("time_convergence requires a fixed cycle count and fixed controls"));
                }
                (count(get(repeat,"cycles")?,"repeat.cycles",4096)?,
                    count(get(repeat,"max_total_steps")?,"repeat.max_total_steps",1_000_000)?)
            }
        };
        let max_refinements = count(get(value,"max_refinements")?,"time max_refinements",16)?;
        let consecutive = count(get(value,"consecutive_passes")?,"time consecutive_passes",16)?;
        if consecutive < 2 || consecutive > max_refinements {
            return Err(bad("time convergence requires 2 <= consecutive_passes <= max_refinements"));
        }
        Ok(Self {
            max_refinements, consecutive,
            tolerance_k: positive(get(value,"temperature_tolerance_k")?,"time temperature_tolerance_k")?,
            max_total_steps: count(get(value,"max_total_steps")?,"time max_total_steps",1_000_000)?,
            max_trace_bytes: count(get(value,"max_trace_bytes")?,"time max_trace_bytes",1_073_741_824)?,
            cycles, max_trajectory_steps,
        })
    }
}

/// One flat allocation: time followed by every nodal temperature per row.
/// No matrices, nonlinear trial fields or discarded adaptive endpoints enter it.
struct Trace {
    width: usize,
    rows: usize,
    data: Vec<f64>,
}
impl Trace {
    fn values(vertices: usize, steps: usize) -> Result<usize> {
        vertices.checked_add(1).and_then(|width| steps.checked_add(1)
            .and_then(|rows| rows.checked_mul(width)))
            .ok_or_else(|| budget("time-study trace size overflow"))
    }
    fn bytes(values: usize) -> Result<usize> {
        values.checked_mul(std::mem::size_of::<f64>())
            .and_then(|n| n.checked_add(std::mem::size_of::<Self>()))
            .ok_or_else(|| budget("time-study trace byte count overflow"))
    }
    fn new(vertices: usize, steps: usize) -> Result<Self> {
        let values = Self::values(vertices,steps)?;
        let mut data = Vec::new();
        data.try_reserve_exact(values).map_err(|_| budget("cannot allocate the admitted time-study trace"))?;
        Ok(Self { width: vertices+1, rows: steps+1, data })
    }
    fn record(&mut self, time: f64, field: &[f64]) -> Result<()> {
        if field.len()+1 != self.width || self.data.len()/self.width >= self.rows
            || !time.is_finite() || field.iter().any(|t| !t.is_finite() || *t <= 0.0) {
            return Err(producer("time-study observer received an invalid accepted field"));
        }
        if !self.data.is_empty() && time <= self.data[self.data.len()-self.width] {
            return Err(producer("time-study observation times are not strictly increasing"));
        }
        self.data.push(time);
        self.data.extend_from_slice(field);
        Ok(())
    }
    fn complete(&self) -> Result<()> {
        if self.data.len() != self.rows*self.width {
            return Err(producer("time-study trace lacks a complete accepted trajectory"));
        }
        Ok(())
    }
    fn compare(&self, coarse: &Self, cx: &Cx<'_>) -> Result<f64> {
        self.complete()?; coarse.complete()?;
        if self.width != coarse.width || self.rows-1 != 2*(coarse.rows-1) {
            return Err(producer("time-study grids are not nested by complete interval subdivision"));
        }
        let mut largest = 0.0_f64;
        for (row, old) in coarse.data.chunks_exact(coarse.width).enumerate() {
            if row%64 == 0 { poll(cx)?; }
            let start = 2*row*self.width;
            let new = &self.data[start..start+self.width];
            // No interpolation across source/fan switches, rounded matching,
            // or comparing different phases of a repeated trajectory.
            if old[0].to_bits() != new[0].to_bits() {
                return Err(producer("timestep refinement changed a common physical observation time"));
            }
            for (&a,&b) in old[1..].iter().zip(&new[1..]) {
                largest = largest.max(finite(a-b)?.abs());
            }
        }
        Ok(largest)
    }
}

fn refined(schedule: &Schedule) -> Result<Schedule> {
    let mut intervals = Vec::with_capacity(schedule.intervals.len());
    let mut total_steps = 0_usize;
    for interval in &schedule.intervals {
        let steps = interval.steps.checked_mul(2).ok_or_else(|| budget("refined step count overflow"))?;
        total_steps = total_steps.checked_add(steps).ok_or_else(|| budget("refined step count overflow"))?;
        if total_steps > schedule.max_steps { return Err(budget("refined trajectory exceeds transient.max_steps")); }
        intervals.push(Interval { duration: interval.duration, workload: interval.workload.clone(),
            speed: interval.speed, steps });
    }
    Ok(Schedule { initial: schedule.initial.clone(), capacities: schedule.capacities.clone(), intervals,
        limit: schedule.limit, total_steps, max_step_s: schedule.max_step_s, max_steps: schedule.max_steps,
        adaptive: None, nonlinear: schedule.nonlinear, adjoint: None, time_convergence: None,
        fan_speed_design: None, power_design: None, repeat: schedule.repeat })
}

/// Admit representable endpoint times BEFORE any solve on the proposed grid.
fn admit_times(schedule: &Schedule) -> Result<()> {
    let mut start = 0.0;
    for interval in &schedule.intervals {
        let end = finite(start+interval.duration)?;
        let mut previous = start;
        for step in 1..=interval.steps {
            let time = if step == interval.steps { end }
                else { start+interval.duration*(step as f64/interval.steps as f64) };
            if !(time.is_finite() && time > previous) {
                return Err(budget("refined physical timesteps are not representable"));
            }
            previous = time;
        }
        start = end;
    }
    Ok(())
}

pub(super) fn solve(request: &Request, cx: &Cx<'_>, base: &Schedule, config: Config) -> Result<String> {
    let mut owned: Option<Schedule> = None;
    let mut previous: Option<(Trace,f64)> = None;
    let mut history = Vec::new();
    let mut steps_spent = 0_usize;
    let mut solves_spent = 0_usize;
    let mut streak = 0;
    for level in 0..=config.max_refinements {
        poll(cx)?;
        let schedule = owned.as_ref().unwrap_or(base);
        admit_times(schedule)?;
        let steps = schedule.total_steps.checked_mul(config.cycles)
            .ok_or_else(|| budget("time-study trajectory step count overflow"))?;
        if steps > config.max_trajectory_steps {
            return Err(budget("refined trajectory exceeds its original repeated-cycle step cap"));
        }
        let after = steps_spent.checked_add(steps).ok_or_else(|| budget("time-study step count overflow"))?;
        if after > config.max_total_steps {
            return Err(budget("time-study cumulative accepted-step budget exhausted; no convergence result published"));
        }
        let values = Trace::values(request.mesh.vertex_count(),steps)?;
        let old_bytes = previous.as_ref().map(|(t,_)| Trace::bytes(t.data.capacity()))
            .transpose()?.unwrap_or(0);
        if Trace::bytes(values)?.checked_add(old_bytes).is_none_or(|n| n > config.max_trace_bytes) {
            return Err(budget("time-study retained field traces exceed max_trace_bytes"));
        }
        let mut trace = Trace::new(request.mesh.vertex_count(),steps)?;
        if Trace::bytes(trace.data.capacity())?.checked_add(old_bytes).is_none_or(|n| n > config.max_trace_bytes) {
            return Err(budget("time-study allocator capacity exceeds max_trace_bytes"));
        }
        trace.record(0.0,&schedule.initial)?;
        let run = simulate_observed(request,cx,schedule,1.0,
            Some(&mut |time,field| trace.record(time,field)))?;
        trace.complete()?;
        if run.steps != steps { return Err(producer("time-study producer returned the wrong complete step count")); }
        steps_spent = after;
        solves_spent = solves_spent.checked_add(run.solid_solves).ok_or_else(|| budget("time-study solve count overflow"))?;
        let (field_change,peak_change) = match &previous {
            None => (None,None),
            Some((coarse,peak)) => (Some(trace.compare(coarse,cx)?),Some(finite(run.peak_k-peak)?.abs())),
        };
        let passed = field_change.is_some_and(|d| d <= config.tolerance_k)
            && peak_change.is_some_and(|d| d <= config.tolerance_k);
        streak = if passed {streak+1} else {0};
        let counts = schedule.intervals.iter().map(|i| i.steps.to_string()).collect::<Vec<_>>().join(",");
        let maximum_dt = schedule.intervals.iter().map(|i| i.duration/i.steps as f64).fold(0.0_f64,f64::max);
        history.push(format!("{{\"level\":{level},\"steps_per_interval\":[{counts}],\"accepted_steps\":{steps},\"maximum_nominal_step_s\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"common_endpoint_field_change_k\":{},\"sampled_peak_change_k\":{},\"solid_solves\":{},\"consecutive_passes\":{streak}}}",
            num(maximum_dt)?,num(run.peak_k)?,num(run.peak_time_s)?,optional(field_change)?,optional(peak_change)?,run.solid_solves));
        if streak >= config.consecutive {
            poll(cx)?;
            let prefix = run.output.strip_suffix("}\n").ok_or_else(|| bad("time-study result framing"))?;
            return Ok(format!("{prefix},\"time_convergence\":{{\"status\":\"successive-time-grid-tolerance-met\",\"method\":\"interval-doubling-backward-euler\",\"refinements\":{level},\"trajectories_solved\":{},\"temperature_tolerance_k\":{},\"required_consecutive_passes\":{},\"total_accepted_steps\":{steps_spent},\"total_solid_solves\":{solves_spent},\"final_steps_per_interval\":[{counts}],\"history\":[{}],\"scope\":\"observed full-field agreement at common endpoints plus all-trajectory sampled-peak agreement; every level restarts the original initial state and carries history only within its fixed cycles; new intermediate endpoints enter the peak; no interpolation over workload/fan switches, Richardson certificate, continuous-time maximum bound, mesh-error claim or physical validation; max_trace_bytes covers the two retained field buffers, not total solver memory; to replay remove transient.time_convergence and set each interval.steps to final_steps_per_interval without changing other inputs\"}}}}\n",
                history.len(),num(config.tolerance_k)?,config.consecutive,history.join(",")));
        }
        previous = Some((trace,run.peak_k));
        if level == config.max_refinements { break; }
        owned = Some(refined(schedule)?);
    }
    Err(budget("time refinement limit exhausted before consecutive full-field and peak agreement; no convergence result published"))
}

#[cfg(test)]
mod tests;
