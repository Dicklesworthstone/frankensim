//! Repeated duty cycles use the actual transient map, not repeated cold starts.
//! Only the last cycle keeps full output; bounded scalar summaries retain the
//! warm-up peaks and energy of every completed cycle for design decisions.
//! Fixed-count adjoints additionally retain one admitted global state tape.
//!
//! An optional cycle-rate thermostat samples one explicit solid vertex at the
//! start of each complete cycle and selects a declared low/high fan multiplier.
//! Hysteresis state is carried between cycles; it never changes inside a cycle.

use super::*;

#[derive(Debug, Clone, Copy)]
pub(super) struct Config {
    cycles: usize,
    max_total_steps: usize,
    periodic: Option<Periodic>,
    controller: Option<FanController>,
}

#[derive(Debug, Clone, Copy)]
struct Periodic { tolerance_k: f64, consecutive: usize }

impl Periodic {
    fn observe(self, residual: f64, streak: &mut usize) -> bool {
        *streak = if residual <= self.tolerance_k { *streak + 1 } else { 0 };
        *streak >= self.consecutive
    }
}

#[derive(Debug, Clone, Copy)]
struct FanController {
    sensor_vertex: usize,
    low_temperature_k: f64,
    high_temperature_k: f64,
    low_speed_multiplier: f64,
    high_speed_multiplier: f64,
    initial_speed_multiplier: f64,
}

#[derive(Debug, Clone, Copy)]
struct FanControlState {
    multiplier: f64,
    switches: usize,
}

impl FanController {
    fn parse(value: &J) -> Result<Self> {
        object(value, &["sensor_vertex", "low_temperature_k", "high_temperature_k",
            "low_speed_multiplier", "high_speed_multiplier", "initial_speed_multiplier"],
            "repeat.fan_controller")?;
        let controller = Self {
            sensor_vertex: integer_raw(get(value, "sensor_vertex")?, "fan_controller.sensor_vertex")?,
            low_temperature_k: positive(get(value, "low_temperature_k")?, "fan_controller.low_temperature_k")?,
            high_temperature_k: positive(get(value, "high_temperature_k")?, "fan_controller.high_temperature_k")?,
            low_speed_multiplier: positive(get(value, "low_speed_multiplier")?, "fan_controller.low_speed_multiplier")?,
            high_speed_multiplier: positive(get(value, "high_speed_multiplier")?, "fan_controller.high_speed_multiplier")?,
            initial_speed_multiplier: positive(get(value, "initial_speed_multiplier")?, "fan_controller.initial_speed_multiplier")?,
        };
        if controller.low_temperature_k >= controller.high_temperature_k {
            return Err(bad("fan-controller temperature thresholds must be strictly ordered"));
        }
        if controller.low_speed_multiplier >= controller.high_speed_multiplier {
            return Err(bad("fan-controller low/high speed multipliers must be strictly ordered"));
        }
        if controller.initial_speed_multiplier < controller.low_speed_multiplier
            || controller.initial_speed_multiplier > controller.high_speed_multiplier {
            return Err(bad("fan-controller initial speed multiplier must lie within the declared low/high range"));
        }
        Ok(controller)
    }

    fn validate(self, request: &Request, schedule: &Schedule, outer_multiplier: f64) -> Result<()> {
        if self.sensor_vertex >= request.mesh.vertex_count() {
            return Err(bad("fan-controller sensor_vertex is outside the solid mesh"));
        }
        let fan = request.fan.as_ref().ok_or_else(||bad("repeat.fan_controller requires hydraulics.fan"))?;
        for interval in &schedule.intervals {
            let base = interval.speed.ok_or_else(||bad("fan-controlled repetition requires every interval's declared base fan speed"))?;
            for factor in [self.low_speed_multiplier, self.high_speed_multiplier, self.initial_speed_multiplier] {
                fan.bank(finite(base * outer_multiplier * factor)?)?;
            }
        }
        Ok(())
    }

    fn next_multiplier(self, sensor_k: f64, current: f64) -> Result<f64> {
        finite(sensor_k)?;
        finite(current)?;
        if sensor_k >= self.high_temperature_k {
            Ok(self.high_speed_multiplier)
        } else if sensor_k <= self.low_temperature_k {
            Ok(self.low_speed_multiplier)
        } else {
            Ok(current)
        }
    }

    fn update(self, sensor_k: f64, state: &mut FanControlState) -> Result<bool> {
        let next = self.next_multiplier(sensor_k, state.multiplier)?;
        let switched = next.to_bits() != state.multiplier.to_bits();
        if switched {
            state.switches = state.switches.checked_add(1)
                .ok_or_else(||budget("fan-controller switch count overflow"))?;
            state.multiplier = next;
        }
        Ok(switched)
    }
}

impl Config {
    pub(super) fn parse(value: &J, planned_steps: usize, adaptive: bool) -> Result<Self> {
        object(value, &["cycles", "until_periodic", "max_total_steps", "fan_controller"], "transient.repeat")?;
        let (cycles, periodic) = match (value.get("cycles"), value.get("until_periodic")) {
            (Some(cycles), None) => (count(cycles,"repeat.cycles",4096)?, None),
            (None, Some(value)) => {
                object(value, &["max_cycles", "temperature_tolerance_k", "consecutive_cycles"], "repeat.until_periodic")?;
                let cycles = count(get(value,"max_cycles")?,"periodic.max_cycles",4096)?;
                let periodic = Periodic {
                    tolerance_k: positive(get(value,"temperature_tolerance_k")?,"periodic.temperature_tolerance_k")?,
                    consecutive: count(get(value,"consecutive_cycles")?,"periodic.consecutive_cycles",4096)?,
                };
                if periodic.consecutive < 2 || periodic.consecutive > cycles {
                    return Err(bad("periodic convergence requires 2 <= consecutive_cycles <= max_cycles"));
                }
                (cycles, Some(periodic))
            }
            _ => return Err(bad("repeat requires either cycles or until_periodic, not both")),
        };
        let config = Self { cycles, periodic,
            max_total_steps: count(get(value,"max_total_steps")?,"repeat.max_total_steps",1_000_000)?,
            controller: value.get("fan_controller").map(FanController::parse).transpose()?,
        };
        let minimum = planned_steps.checked_mul(if adaptive {2} else {1})
            .and_then(|steps|steps.checked_mul(config.periodic.map_or(config.cycles, |p|p.consecutive)))
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
    let mut tape = if schedule.adjoint.is_some() {
        if config.periodic.is_some() || config.controller.is_some() || speed_multiplier != 1.0 {
            return Err(bad("repeated adjoints require fixed cycles and declared speeds without periodic stopping or a controller"));
        }
        let (initial,vertex) = initial_objective(request,cx,&schedule.initial)?;
        adjoint::Tape::for_cycles(request,schedule,initial,vertex,config.cycles)?
    } else { None };
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
    let mut streak = 0_usize;
    let mut last_residual = None;
    let mut controller_state = config.controller.map(|controller| FanControlState {
        multiplier: controller.initial_speed_multiplier,
        switches: 0,
    });
    if let Some(controller) = config.controller {
        controller.validate(request, schedule, speed_multiplier)?;
    }
    for cycle_index in 0..config.cycles {
        poll(cx)?;
        let sensor_start = config.controller.map(|controller| field[controller.sensor_vertex]);
        let mut switched = false;
        let controller_multiplier = match (config.controller, controller_state.as_mut()) {
            (Some(controller), Some(state)) => {
                switched = controller.update(sensor_start.expect("controller has sensor"), state)?;
                state.multiplier
            }
            (None, None) => 1.0,
            _ => return Err(bad("internal fan-controller state mismatch")),
        };
        let applied_multiplier = finite(speed_multiplier * controller_multiplier)?;
        if let Some(controller) = config.controller {
            controller.validate(request, schedule, speed_multiplier)?;
            let fan = request.fan.as_ref().ok_or_else(||bad("repeat.fan_controller requires hydraulics.fan"))?;
            for interval in &schedule.intervals {
                let base = interval.speed.ok_or_else(||bad("fan-controlled repetition requires base fan speeds"))?;
                fan.bank(finite(base * applied_multiplier)?)?;
            }
        }
        let remaining = config.max_total_steps.checked_sub(steps)
            .filter(|&left|left>0).ok_or_else(||budget("repeated-cycle accepted-step budget exhausted"))?;
        let cycle = simulate_cycle_recorded(request,cx,schedule,applied_multiplier,&field,remaining,
            tape.as_mut().map(|tape| (tape,elapsed)))?;
        let residual = field_residual(cx,&field,&cycle.final_temperature)?;
        last_residual = Some(residual);
        let sensor_end = config.controller.map(|controller| cycle.final_temperature[controller.sensor_vertex]);
        let controller_stable = match (config.controller, controller_state) {
            (Some(controller), Some(state)) => controller.next_multiplier(
                sensor_end.expect("controller has end sensor"), state.multiplier)?.to_bits() == state.multiplier.to_bits(),
            (None, None) => true,
            _ => return Err(bad("internal fan-controller state mismatch")),
        };
        let complete = match config.periodic {
            Some(periodic) if controller_stable => periodic.observe(residual, &mut streak),
            Some(_) => { streak = 0; false },
            None => cycle_index+1==config.cycles,
        };
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
        summaries.push(format!("{{\"cycle\":{},\"start_time_s\":{},\"end_time_s\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"start_to_end_field_residual_k\":{},\"controller_sensor_start_k\":{},\"controller_sensor_end_k\":{},\"controller_speed_multiplier\":{},\"applied_speed_multiplier\":{},\"controller_switched\":{},\"controller_state_stable\":{},\"stored_energy_change_j\":{},\"input_energy_j\":{},\"air_energy_gain_j\":{},\"accepted_steps\":{},\"solid_solves\":{}}}",
            cycle_index+1,num(elapsed)?,num(end)?,num(cycle.trajectory.peak_k)?,num(global_peak_time)?,
            num(residual)?,optional(sensor_start)?,optional(sensor_end)?,
            if config.controller.is_some(){num(controller_multiplier)?}else{"null".into()},
            num(applied_multiplier)?, if config.controller.is_some(){if switched{"true"}else{"false"}}else{"null"},
            if config.controller.is_some(){if controller_stable{"true"}else{"false"}}else{"null"},
            num(cycle.stored_j)?,num(cycle.input_j)?,num(cycle.exhaust_j)?,
            cycle.trajectory.steps,cycle.trajectory.solid_solves));
        if complete {
            let status = if config.periodic.is_some() {"periodic-field-tolerance-met"} else {"fixed-count-complete"};
            let periodic = match config.periodic {
                None => "null".to_string(),
                Some(p) => format!("{{\"temperature_tolerance_k\":{},\"full_field_residual_k\":{},\"consecutive_cycles_met\":{},\"required_consecutive_cycles\":{},\"max_cycles\":{},\"controller_state_stable\":{},\"criterion\":\"same-phase maximum absolute nodal start/end difference plus unchanged hysteretic controller state; a cycle-map residual, not distance to the infinite-cycle solution\"}}",
                    num(p.tolerance_k)?, num(residual)?, streak, p.consecutive, config.cycles,
                    if controller_stable{"true"}else{"false"}),
            };
            let fan_controller = match (config.controller, controller_state) {
                (None, None) => "null".to_string(),
                (Some(controller), Some(state)) => format!(
                    "{{\"sensor_vertex\":{},\"low_temperature_k\":{},\"high_temperature_k\":{},\"low_speed_multiplier\":{},\"high_speed_multiplier\":{},\"initial_speed_multiplier\":{},\"final_speed_multiplier\":{},\"switches\":{},\"update_phase\":\"cycle-start from prior accepted cycle-end sensor; multiplier held for the complete next cycle\",\"scope\":\"ideal sampled thermostat; no sensor lag/noise, actuator dynamics, PWM, electrical-power or sub-cycle feedback model\"}}",
                    controller.sensor_vertex,num(controller.low_temperature_k)?,num(controller.high_temperature_k)?,
                    num(controller.low_speed_multiplier)?,num(controller.high_speed_multiplier)?,
                    num(controller.initial_speed_multiplier)?,num(state.multiplier)?,state.switches),
                _ => return Err(bad("internal fan-controller state mismatch")),
            };
            // Every forward cycle and the cumulative energy gate have passed.
            // Reverse reconstructs endpoints but never advances physical history.
            let forward_work = work;
            let (adjoint,reverse_work) = match tape.take() {
                Some(tape) => {
                    let engine = BackwardEuler::per_element(cx,&request.mesh,&schedule.capacities).map_err(producer)?;
                    tape.reverse(request,cx,schedule,&engine)?
                }
                None => ("null".into(),0),
            };
            work = work.checked_add(reverse_work).ok_or_else(||budget("repeated adjoint work count overflow"))?;
            let prefix = cycle.trajectory.output.strip_suffix("}\n").ok_or_else(||bad("internal cycle result framing"))?;
            let output = format!("{prefix},\"repeated_cycles\":{{\"status\":{},\"periodic\":{},\"fan_controller\":{},\"cycles_completed\":{},\"cycle_duration_s\":{},\"elapsed_time_s\":{},\"last_cycle_start_time_s\":{},\"total_accepted_steps\":{},\"total_solid_solves\":{},\"forward_solid_solves\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"temperature_limit_k\":{},\"first_sampled_violation_s\":{},\"input_energy_j\":{},\"stored_energy_change_j\":{},\"air_energy_gain_j\":{},\"energy_residual_j\":{},\"cycles\":[{}],\"adjoint\":{},\"scope\":\"all cycles inherit the prior accepted nodal field; optional hysteretic fan control samples one declared vertex at cycle start and holds speed through that cycle; peaks include initial state and every accepted sample; transient contains only the final cycle in local time; periodic stopping checks field and controller state; an optional fixed-count adjoint uses the complete history and global time, not just the last cycle; no infinite-cycle, future-peak or continuous-time bound\"}}}}\n",
                quote(status),periodic,fan_controller,cycle_index+1,num(cycle.duration_s)?,num(end)?,num(elapsed)?,steps,work,forward_work,num(peak)?,num(peak_time)?,
                optional(schedule.limit)?,optional(first_violation)?,num(input)?,num(stored)?,num(exhaust)?,num(energy_residual)?,summaries.join(","),adjoint);
            poll(cx)?;
            return Ok(Trajectory {output,peak_k:peak,peak_time_s:peak_time,solid_solves:work,steps});
        }
        field = cycle.final_temperature;
        elapsed = end;
    }
    match (config.periodic, last_residual) {
        (Some(p), Some(residual)) => Err(budget(&format!(
            "periodic cycle budget exhausted after {} cycles: full-field residual {residual} K, tolerance {} K, consecutive passes {streak}/{}; no converged trajectory published",
            config.cycles,p.tolerance_k,p.consecutive))),
        _ => Err(bad("repeat.cycles must be positive")),
    }
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

#[cfg(test)]
mod periodic_tests;

#[cfg(test)]
mod control_tests;
