//! Repeated schedules carry accepted h through the same cycle-local producer.
//! Scalar summaries retain every warm-up peak; an admitted fixed-count tape
//! owns the complete chronological history and shared interval controls.

use super::*;

pub(super) fn simulate(
    request: &Request,
    cx: &Cx<'_>,
    schedule: &Schedule,
    config: &Config,
    engine: &Prepared<'_, '_>,
    policy: repeat::EnthalpyPolicy,
) -> Result<Trajectory> {
    poll(cx)?;
    let mut h = config.initial_h.clone();
    let mut temperature = config.initial_temperatures()?;
    let (initial, vertex) = initial_objective(request, cx, &temperature)?;
    let mut tape = adjoint::Tape::for_cycles(request, schedule, initial, vertex, policy.cycles)?;
    let mut elapsed = 0.0;
    let mut peak = f64::NEG_INFINITY;
    let mut peak_time = 0.0;
    let mut first_violation = None;
    let mut steps = 0_usize;
    let mut work = 0_usize;
    let (mut input, mut stored, mut exhaust, mut radiative) = (0.0, 0.0, 0.0, 0.0);
    let mut fresh_exhaust = 0.0;
    let mut summaries = Vec::new();
    let mut streak = 0_usize;
    let mut last_residuals = None;
    for cycle_index in 0..policy.cycles {
        poll(cx)?;
        let remaining = policy
            .max_total_steps
            .checked_sub(steps)
            .filter(|&left| left > 0)
            .ok_or_else(|| budget("repeated enthalpy accepted-step budget exhausted"))?;
        // The global offset never enters the physical dt calculation. The
        // accepted h field is the only storage history passed to the producer.
        let cycle = simulate_cycle(
            request,
            cx,
            schedule,
            config,
            engine,
            &h,
            remaining,
            tape.as_mut().map(|tape| (tape, elapsed)),
        )?;
        let temperature_residual = field_residual(cx, &temperature, &cycle.final_temperature)?;
        let enthalpy_residual = field_residual(cx, &h, &cycle.final_h)?;
        last_residuals = Some((temperature_residual, enthalpy_residual));
        let complete = match policy.periodic {
            Some(periodic) => {
                periodic.observe(temperature_residual, enthalpy_residual, &mut streak)?
            }
            None => cycle_index + 1 == policy.cycles,
        };
        let end = finite(elapsed + cycle.duration_s)?;
        if end <= elapsed {
            return Err(bad(
                "elapsed enthalpy cycle time is no longer representable",
            ));
        }
        let global_peak_time = finite(elapsed + cycle.trajectory.peak_time_s)?;
        if cycle.trajectory.peak_k > peak {
            peak = cycle.trajectory.peak_k;
            peak_time = global_peak_time;
        }
        if first_violation.is_none() {
            first_violation = cycle
                .first_violation_s
                .map(|time| finite(elapsed + time))
                .transpose()?;
        }
        steps = steps
            .checked_add(cycle.trajectory.steps)
            .ok_or_else(|| budget("repeated enthalpy step count overflow"))?;
        work = work
            .checked_add(cycle.trajectory.solid_solves)
            .ok_or_else(|| budget("repeated enthalpy solid-work count overflow"))?;
        input = finite(input + cycle.input_j)?;
        stored = finite(stored + cycle.stored_j)?;
        exhaust = finite(exhaust + cycle.exhaust_j)?;
        fresh_exhaust = finite(fresh_exhaust + cycle.fresh_exhaust_j)?;
        radiative = finite(radiative + cycle.radiative_j)?;
        let energy_residual = finite(stored - input + exhaust + radiative)?;
        if energy_residual.abs() > finite(request.limits.heat * end)? {
            return Err(producer("repeated enthalpy cumulative energy gate failed"));
        }
        let fresh_energy_residual = finite(stored - input + fresh_exhaust + radiative)?;
        if request.recirculation.is_some()
            && fresh_energy_residual.abs() > finite(request.limits.heat * end)?
        {
            return Err(producer(
                "repeated enthalpy cumulative fresh/exhaust energy gate failed",
            ));
        }
        let radiation_field = if request.radiation.is_some() {
            format!(",\"radiative_energy_loss_j\":{}", num(cycle.radiative_j)?)
        } else {
            String::new()
        };
        let recirculation_field = if request.recirculation.is_some() {
            format!(
                ",\"fresh_exhaust_energy_gain_j\":{},\"fresh_exhaust_energy_residual_j\":{}",
                num(cycle.fresh_exhaust_j)?,
                num(finite(
                    cycle.stored_j - cycle.input_j + cycle.fresh_exhaust_j + cycle.radiative_j,
                )?)?,
            )
        } else {
            String::new()
        };
        summaries.push(format!(
            "{{\"cycle\":{},\"start_time_s\":{},\"end_time_s\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"start_to_end_field_residual_k\":{},\"start_to_end_specific_enthalpy_residual_j_kg\":{},\"stored_energy_change_j\":{},\"input_energy_j\":{},\"air_energy_gain_j\":{},\"accepted_steps\":{},\"solid_solves\":{}{radiation_field}{recirculation_field}}}",
            cycle_index + 1,
            num(elapsed)?,
            num(end)?,
            num(cycle.trajectory.peak_k)?,
            num(global_peak_time)?,
            num(temperature_residual)?,
            num(enthalpy_residual)?,
            num(cycle.stored_j)?,
            num(cycle.input_j)?,
            num(cycle.exhaust_j)?,
            cycle.trajectory.steps,
            cycle.trajectory.solid_solves,
        ));
        if complete {
            let status = if policy.periodic.is_some() {
                "periodic-field-tolerance-met"
            } else {
                "fixed-count-complete"
            };
            let periodic = match policy.periodic {
                None => "null".into(),
                Some(p) => format!(
                    "{{\"temperature_tolerance_k\":{},\"specific_enthalpy_tolerance_j_kg\":{},\"full_field_residual_k\":{},\"full_specific_enthalpy_residual_j_kg\":{},\"consecutive_cycles_met\":{},\"required_consecutive_cycles\":{},\"max_cycles\":{},\"criterion\":\"same-phase maximum absolute nodal start/end differences in both temperature and specific enthalpy must pass on consecutive cycles; a cycle-map residual, not distance to the infinite-cycle solution\"}}",
                    num(p.temperature_tolerance_k)?,
                    num(p.specific_enthalpy_tolerance_j_kg)?,
                    num(temperature_residual)?,
                    num(enthalpy_residual)?,
                    streak,
                    p.consecutive_cycles,
                    policy.cycles,
                ),
            };
            // All forward cycles and their cumulative energy gate have passed.
            // Reverse replay never advances the accepted physical history.
            let forward_work = work;
            let (adjoint, reverse_work, design_gradient) = match tape.take() {
                Some(tape) => tape.reverse(request, cx, schedule, config, engine)?,
                None => ("null".into(), 0, None),
            };
            work = work
                .checked_add(reverse_work)
                .ok_or_else(|| budget("repeated enthalpy adjoint work count overflow"))?;
            let prefix = cycle
                .trajectory
                .output
                .strip_suffix("}\n")
                .ok_or_else(|| bad("internal enthalpy cycle result framing"))?;
            let radiation_field = if request.radiation.is_some() {
                format!(",\"radiative_energy_loss_j\":{}", num(radiative)?)
            } else {
                String::new()
            };
            let recirculation_field = if request.recirculation.is_some() {
                format!(
                    ",\"fresh_exhaust_energy_gain_j\":{},\"fresh_exhaust_energy_residual_j\":{}",
                    num(fresh_exhaust)?,
                    num(fresh_energy_residual)?,
                )
            } else {
                String::new()
            };
            let output = format!(
                "{prefix},\"repeated_cycles\":{{\"status\":{},\"periodic\":{},\"fan_controller\":null,\"cycles_completed\":{},\"cycle_duration_s\":{},\"elapsed_time_s\":{},\"last_cycle_start_time_s\":{},\"total_accepted_steps\":{},\"total_solid_solves\":{},\"forward_solid_solves\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"temperature_limit_k\":{},\"first_sampled_violation_s\":{},\"input_energy_j\":{},\"stored_energy_change_j\":{},\"air_energy_gain_j\":{},\"energy_residual_j\":{},\"cycles\":[{}],\"adjoint\":{}{radiation_field}{recirculation_field},\"scope\":\"all cycles inherit the prior accepted specific enthalpy; peaks include the original initial state, every cycle boundary and every accepted endpoint, including warm-up; transient contains only the final cycle in local time; periodic stopping requires both full nodal temperature and enthalpy residuals; an optional fixed-count adjoint uses the complete h history, shared interval controls and global time; no controller, variable-stopping derivative, infinite-cycle or continuous-time peak bound\"}}}}\n",
                quote(status),
                periodic,
                cycle_index + 1,
                num(cycle.duration_s)?,
                num(end)?,
                num(elapsed)?,
                steps,
                work,
                forward_work,
                num(peak)?,
                num(peak_time)?,
                optional(schedule.limit)?,
                optional(first_violation)?,
                num(input)?,
                num(stored)?,
                num(exhaust)?,
                num(energy_residual)?,
                summaries.join(","),
                adjoint,
            );
            poll(cx)?;
            return Ok(Trajectory {
                output,
                peak_k: peak,
                peak_time_s: peak_time,
                solid_solves: work,
                steps,
                design_gradient,
            });
        }
        h = cycle.final_h;
        temperature = cycle.final_temperature;
        elapsed = end;
    }
    match (policy.periodic, last_residuals) {
        (Some(p), Some((temperature, enthalpy))) => Err(budget(&format!(
            "periodic enthalpy cycle budget exhausted after {} cycles: full-field residual {temperature} K (tolerance {} K), specific-enthalpy residual {enthalpy} J/kg (tolerance {} J/kg), consecutive passes {streak}/{}; no converged trajectory published",
            policy.cycles,
            p.temperature_tolerance_k,
            p.specific_enthalpy_tolerance_j_kg,
            p.consecutive_cycles,
        ))),
        _ => Err(bad("repeat.cycles must be positive")),
    }
}

fn field_residual(cx: &Cx<'_>, start: &[f64], end: &[f64]) -> Result<f64> {
    if start.is_empty() || start.len() != end.len() {
        return Err(bad(
            "enthalpy cycle boundary fields have incompatible lengths",
        ));
    }
    let mut residual = 0.0_f64;
    for (index, (&a, &b)) in start.iter().zip(end).enumerate() {
        if index % 512 == 0 {
            poll(cx)?;
        }
        finite(a)?;
        finite(b)?;
        residual = residual.max(finite(b - a)?.abs());
    }
    Ok(residual)
}
