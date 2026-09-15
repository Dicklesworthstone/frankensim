//! Piecewise-constant workload/fan schedules with backward-Euler solid storage.
//! Air is quasi-steady at each endpoint. Coupling iterations always reuse the
//! previous accepted solid field. No transient derivative or continuous-time
//! peak bound is inferred from sampled temperatures.

mod workload;
mod adaptive;
mod sizing;
use workload::Workload;

use super::*;
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig, StepSolution};
use fs_conduction::transient::VolumetricHeatCapacity;

#[derive(Debug)]
pub(super) struct Interval {
    duration: f64,
    workload: Workload,
    speed: Option<f64>,
    steps: usize,
}
#[derive(Debug)]
pub(super) struct Schedule {
    initial: Vec<f64>,
    capacities: Vec<VolumetricHeatCapacity>,
    intervals: Vec<Interval>,
    limit: Option<f64>,
    total_steps: usize,
    max_step_s: f64,
    max_steps: usize,
    adaptive: Option<adaptive::Config>,
    fan_speed_design: Option<sizing::Config>,
    power_design: Option<sizing::Config>,
}

impl Schedule {
    pub(super) fn parse(value: &J, vertices: usize, elements: usize, fan: Option<&fan_drive::FanDrive>) -> Result<Self> {
        object(value, &["initial_temperature_k", "initial_temperatures_k", "volumetric_heat_capacity_j_m3_k",
            "element_heat_capacities_j_m3_k", "max_step_s", "max_steps", "intervals", "temperature_limit_k", "adaptive", "fan_speed_design", "power_design"], "transient")?;
        let initial = match (value.get("initial_temperature_k"), value.get("initial_temperatures_k")) {
            (Some(t),None) => vec![positive(t,"initial_temperature_k")?;vertices],
            (None,Some(ts)) => {
                let values = array(ts,"initial_temperatures_k",vertices)?;
                if values.len()!=vertices { return Err(bad("one initial temperature per solid vertex required")); }
                values.iter().map(|t| positive(t,"initial temperature")).collect::<Result<Vec<_>>>()?
            }
            _ => return Err(bad("choose one uniform or nodal initial temperature")),
        };
        let capacities = match (value.get("volumetric_heat_capacity_j_m3_k"), value.get("element_heat_capacities_j_m3_k")) {
            (Some(c),None) => vec![VolumetricHeatCapacity::declared(positive(c,"volumetric heat capacity")?).map_err(producer)?;elements],
            (None,Some(cs)) => {
                let values = array(cs,"element heat capacities",elements)?;
                if values.len()!=elements { return Err(bad("one heat capacity per tetrahedron required")); }
                values.iter().map(|c| VolumetricHeatCapacity::declared(positive(c,"element heat capacity")?).map_err(producer)).collect::<Result<Vec<_>>>()?
            }
            _ => return Err(bad("choose one uniform or per-element volumetric heat capacity")),
        };
        let max_dt = positive(get(value,"max_step_s")?,"max_step_s")?;
        let max_steps = count(get(value,"max_steps")?,"max_steps",10_000)?;
        let adaptive = value.get("adaptive").map(|v| adaptive::Config::parse(v, max_dt)).transpose()?;
        if value.get("fan_speed_design").is_some() && value.get("power_design").is_some() {
            return Err(bad("choose transient fan-speed sizing or workload-power sizing, not both"));
        }
        let power_design = value.get("power_design").map(sizing::Config::parse_power).transpose()?;
        let fan_speed_design = value.get("fan_speed_design").map(sizing::Config::parse).transpose()?;
        let limit = value.get("temperature_limit_k").map(|t| positive(t,"temperature_limit_k")).transpose()?;
        let mut intervals = Vec::new();
        let mut total_steps = 0_usize;
        let mut time = 0.0;
        for entry in array(get(value,"intervals")?,"intervals",4096)? {
            object(entry,&["duration_s","power_scale","component_powers_w","fan_speed_ratio"],"transient interval")?;
            let duration = positive(get(entry,"duration_s")?,"duration_s")?;
            let workload = Workload::parse(entry)?;
            let speed = match (fan,entry.get("fan_speed_ratio")) {
                (Some(fan),Some(s)) => {
                    let speed = positive(s,"fan_speed_ratio")?;
                    fan.bank(speed)?;
                    Some(speed)
                }
                (None,None) => None,
                _ => return Err(bad("each fan-driven interval requires a speed; pressure-driven intervals must omit it")),
            };
            let count_f = (duration/max_dt).ceil().max(1.0);
            if !count_f.is_finite() || count_f > max_steps as f64 {
                return Err(budget("planned transient steps exceed max_steps"));
            }
            let steps = count_f as usize;
            total_steps = total_steps.checked_add(steps).ok_or_else(|| budget("transient step count overflow"))?;
            if total_steps > max_steps { return Err(budget("planned transient steps exceed max_steps")); }
            let end = finite(time+duration)?;
            // Check the exact endpoint formula used during execution before
            // doing physics. Tiny intervals must not become zero-duration steps.
            let mut previous = time;
            for i in 1..=steps {
                let endpoint = if i==steps { end } else { time+duration*(i as f64/steps as f64) };
                if !(endpoint.is_finite() && endpoint > previous) {
                    return Err(bad("time resolution cannot represent every requested step"));
                }
                previous = endpoint;
            }
            time=end;
            intervals.push(Interval {duration,workload,speed,steps});
        }
        if intervals.is_empty() { return Err(bad("at least one transient interval required")); }
        if adaptive.is_some() && total_steps > max_steps / 2 {
            return Err(budget("adaptive half-step endpoints require at least twice the planned full-step count"));
        }
        let schedule = Self {initial,capacities,intervals,limit,total_steps,max_step_s:max_dt,max_steps,adaptive,fan_speed_design,power_design};
        for design in [&schedule.fan_speed_design,&schedule.power_design].into_iter().flatten() {
            design.validate(&schedule,fan)?;
        }
        Ok(schedule)
    }
}

fn finite(value:f64)->Result<f64> {
    if value.is_finite(){Ok(value)}else{Err(producer("nonfinite transient arithmetic"))}
}
fn budget(message: &str)->Failure { Failure {code:"cooling-network-transient-budget",message:message.into()} }

fn scaled_source(request:&Request, scale:f64)->Result<ScalarField> {
    match &request.solid_data.nodal_source {
        Some(ScalarField::Nodal(values)) => Ok(ScalarField::Nodal(values.iter().map(|v|finite(v*scale)).collect::<Result<_>>()?)),
        Some(ScalarField::Uniform(value)) => Ok(ScalarField::Uniform(finite(value*scale)?)),
        None => Ok(ScalarField::Uniform(finite(request.source*scale)?)),
    }
}

fn initial_objective(request:&Request,cx:&Cx<'_>,field:&[f64])->Result<(f64,Option<usize>)> {
    match &request.objective {
        objective::Objective::MeanWall(name) => {
            let surface=request.surfaces.iter().find(|s|s.name==*name).ok_or_else(||bad("unknown initial mean surface"))?;
            let mut mean=0.0;
            for face in request.mesh.boundary() {
                poll(cx)?;
                if surface.faces.contains(&face.vertices) {
                    for v in face.vertices { mean=finite(mean+(face.area/surface.area/3.0)*field[v as usize])?; }
                }
            }
            Ok((mean,None))
        }
        _ => {let state=request.objective.evaluate(cx,field,&[])?;Ok((state.value,state.vertex))}
    }
}

fn advance(request:&Request,cx:&Cx<'_>,engine:&BackwardEuler<'_>,network:&TransportNetwork<'_>,
    coefficients:&BTreeMap<String,f64>,old:&[f64],source:&ScalarField,dt:f64)
    ->Result<(CoupledTransportSolution,StepSolution)>
{
    let material=fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
    let names=network.regions();
    let gate=ConjugateConfig {max_iterations:request.limits.coupling,
        temperature_tolerance_k:request.limits.temperature,
        balance_tolerance_w:request.limits.heat/(names.len() as f64+1.0),balance_relative_tolerance:0.0,
        relaxation:Relaxation::Fixed {omega:request.limits.relaxation}};
    let config=StepConfig {linear:fs_conduction::LinearConfig {tolerance:request.limits.relative,
        max_iterations:request.limits.linear,restart:60},energy_tolerance_j:finite(request.limits.heat*dt)?};
    let mut last=None;
    let mut failure=None;
    let coupled=solve_coupled_transport(cx,network,&gate,|cx,references| {
        let evaluated=(||->Result<Vec<SolidRegionState>> {
            let boundary=request.boundary(&names,references,coefficients)?;
            let solution=engine.advance(cx,ConductionProblem {mesh:&request.mesh,boundary:&boundary,
                material:&material,element_materials:request.solid_data.element_materials.as_ref(),source},
                request.contacts.as_ref().map(|c|&c.interfaces),old,dt,config).map_err(producer)?;
            let states=names.iter().map(|name|solution.robin_fluxes.iter().find(|f|f.region==*name)
                .map(SolidRegionState::from_robin_flux).ok_or_else(||bad("transient solid lacks a cooling region")))
                .collect::<Result<Vec<_>>>()?;
            last=Some(solution);
            Ok(states)
        })();
        evaluated.map_err(|error|{failure=Some(error);AirflowError::Cancelled {iteration:0,references_k:references.to_vec()}})
    });
    if let Some(error)=failure {return Err(error);}
    let coupled=coupled.map_err(producer)?;
    let solid=last.ok_or_else(||bad("transient coupling returned without a solid response"))?;
    let residual=finite(solid.stored_energy_change_j-dt*(solid.source_w-coupled.transport.external_heat_gain_w))?;
    if residual.abs()>config.energy_tolerance_j {
        return Err(producer(format!("coupled transient energy residual {residual} J exceeds {} J",config.energy_tolerance_j)));
    }
    poll(cx)?;
    Ok((coupled,solid))
}

/// One fully admitted trajectory. Scalar design decisions do not parse JSON
/// receipts or substitute the final-state objective for the sampled peak.
struct Trajectory {
    output: String,
    peak_k: f64,
    peak_time_s: f64,
    solid_solves: usize,
    steps: usize,
}

pub(super) fn solve(request:&Request,cx:&Cx<'_>,schedule:&Schedule)->Result<String> {
    match (&schedule.fan_speed_design,&schedule.power_design) {
        (Some(design),None) | (None,Some(design)) => sizing::solve(request,cx,schedule,design),
        (None,None) => simulate(request,cx,schedule,1.0).map(|run|run.output),
        _ => Err(bad("transient design controls are mutually exclusive")),
    }
}

/// Always begins from the declared initial field, never from another candidate.
fn simulate(request:&Request,cx:&Cx<'_>,schedule:&Schedule,speed_multiplier:f64)->Result<Trajectory> {
    poll(cx)?;
    if request.gradient || request.design.is_some() || request.fan_speed_design.is_some() {
        return Err(bad("transient runs do not reuse steady adjoints or steady target searches"));
    }
    // Bind every named load before any hydraulic or thermal work, including
    // later intervals which might otherwise fail after a long simulation.
    for interval in &schedule.intervals { interval.workload.validate(request,cx)?; }
    let engine=BackwardEuler::per_element(cx,&request.mesh,&schedule.capacities).map_err(producer)?;
    let mut old=schedule.initial.clone();
    let (initial,initial_vertex)=initial_objective(request,cx,&old)?;
    let mut peak=initial;
    let mut peak_time=0.0;
    let mut first_violation=schedule.limit.filter(|&limit|initial>limit).map(|_|0.0);
    let mut history=vec![format!("{{\"time_s\":0,\"objective_temperature_k\":{},\"active_vertex\":{},\"initial_state\":true}}",
        num(initial)?,initial_vertex.map_or_else(||"null".into(),|v|v.to_string()))];
    let mut time=0.0;
    let mut stored=0.0;
    let mut input=0.0;
    let mut exhaust=0.0;
    let mut work=0_usize;
    let mut completed=0_usize;
    let mut adaptive_stats=adaptive::Stats::default();
    let mut final_result=None;
    for (ordinal,interval) in schedule.intervals.iter().enumerate() {
        poll(cx)?;
        let speed=interval.speed.map(|base|finite(base*speed_multiplier)).transpose()?;
        let flow=match (&request.fan,speed) {
            (Some(fan),Some(speed))=>fan.solve(cx,&request.graph,request.limits,speed)?,
            (None,None)=>request.flow(cx)?,
            _=>return Err(bad("transient drive/speed mismatch")),
        };
        let base=request.surfaces.iter().map(|s|(s.name.clone(),s.h)).collect();
        let (coefficients,mut derived)=convection::resolve(request,cx,&flow,&base)?;
        let network=request.transport(cx,&flow,&coefficients)?;
        let load=interval.workload.prepare(request,cx)?;
        let workload_json=interval.workload.render()?;
        let source=&load.source;
        let start=time;
        let end=finite(start+interval.duration)?;
        let mut i=0_usize;
        let mut suggested=interval.duration/interval.steps as f64;
        while time < end {
            poll(cx)?;
            let needed=if schedule.adaptive.is_some(){2}else{1};
            if completed > schedule.max_steps.saturating_sub(needed) || needed > schedule.max_steps {
                return Err(budget("accepted transient endpoint budget exhausted; no partial trajectory published"));
            }
            // All trials, including discarded coarse/rejected ones, undergo the
            // same physical gates. Count their work, never their heat as history.
            let mut trial=|old:&[f64],dt:f64| {
                let (coupled,solid)=advance(request,cx,&engine,&network,&coefficients,old,source,dt)?;
                work=work.checked_add(coupled.iterations).ok_or_else(||budget("transient work count overflow"))?;
                for c in &derived {c.check_direction(&coupled.solid,request.limits.heat)?;}
                if solid.temperature.iter().any(|&t|!t.is_finite()||t<=0.0) {
                    return Err(producer("transient FEM produced a nonpositive absolute temperature; refine the step/mesh"));
                }
                if let Some(expected)=load.expected_power_w {
                    if (solid.source_w-expected).abs()>request.limits.heat {
                        return Err(producer("transient source disagrees with the selected component workload"));
                    }
                }
                Ok((coupled,solid))
            };
            let (samples,estimate)=if let Some(config)=schedule.adaptive {
                let accepted=adaptive::step(cx,&old,time,end,suggested,schedule.max_step_s,
                    config,&mut adaptive_stats,&mut trial)?;
                suggested=accepted.next_trial_s;
                (Vec::from(accepted.samples),Some(accepted.error_ratio))
            } else {
                i+=1;
                let endpoint=if i==interval.steps {end}else{start+interval.duration*(i as f64/interval.steps as f64)};
                (vec![(endpoint,trial(&old,endpoint-time)?)],None)
            };
            // The estimate concerns the pair's final endpoint. Its midpoint is
            // retained for peak detection but has no separate local-error test.
            let last_sample=samples.len()-1;
            for (sample_index,(endpoint,(coupled,solid))) in samples.into_iter().enumerate() {
                let dt=endpoint-time;
                let state=request.objective.evaluate(cx,&solid.temperature,&coupled.solid)?;
                if state.value>peak {peak=state.value;peak_time=endpoint;}
                if first_violation.is_none() && schedule.limit.is_some_and(|limit|state.value>limit) {first_violation=Some(endpoint);}
                stored=finite(stored+solid.stored_energy_change_j)?;
                input=finite(input+dt*solid.source_w)?;
                exhaust=finite(exhaust+dt*coupled.transport.external_heat_gain_w)?;
                completed+=1;
                history.push(format!("{{\"time_s\":{},\"dt_s\":{},\"interval\":{},{},\"fan_speed_ratio\":{},\"objective_temperature_k\":{},\"active_vertex\":{},\"source_w\":{},\"air_heat_gain_w\":{},\"stored_energy_change_j\":{},\"solid_energy_residual_j\":{},\"coupled_energy_residual_j\":{},\"coupling_iterations\":{},\"estimated_local_error_ratio\":{}}}",
                    num(endpoint)?,num(dt)?,ordinal,workload_json,optional(speed)?,num(state.value)?,
                    state.vertex.map_or_else(||"null".into(),|v|v.to_string()),num(solid.source_w)?,
                    num(coupled.transport.external_heat_gain_w)?,num(solid.stored_energy_change_j)?,num(solid.energy_residual_j)?,
                    num(solid.stored_energy_change_j-dt*(solid.source_w-coupled.transport.external_heat_gain_w))?,coupled.iterations,
                    optional(estimate.filter(|_|sample_index==last_sample))?));
                // Atomic numerical boundary: no history changes before every solid,
                // air, objective and energy check for this step has passed.
                old.clone_from(&solid.temperature);
                time=endpoint;
                if ordinal+1==schedule.intervals.len() && endpoint==end {
                    let objective=state.value;
                    let evaluated=Evaluation {coupled,temperatures:solid.temperature,gradient:None,
                        objective,objective_state:state,robin_total_w:solid.robin_out_w,source_total_w:solid.source_w,
                        htc:network.regions().iter().map(|name|coefficients[*name]).collect(),
                        convection:std::mem::take(&mut derived),contact_fluxes:solid.contact_fluxes};
                    let result=render(request,&flow,&evaluated)?;
                    final_result=Some(match (&request.fan,speed) {
                        (Some(fan),Some(speed))=>fan.attach(result,&flow,speed)?,_=>result,
                    });
                }
            }
        }
    }
    if schedule.adaptive.is_none() && completed!=schedule.total_steps {
        return Err(bad("completed fixed timestep count differs from the admitted schedule"));
    }
    let residual=finite(stored-input+exhaust)?;
    if residual.abs()>finite(request.limits.heat*time)? {return Err(producer("whole-window transient energy gate failed"));}
    let result=final_result.ok_or_else(||bad("transient run has no completed final step"))?;
    let prefix=result.strip_suffix("}\n").ok_or_else(||bad("internal transient result framing"))?;
    let output=format!("{prefix},\"transient\":{{\"scheme\":\"backward-euler\",\"air_model\":\"quasi-steady endpoint mixing; no fluid storage or travel delay\",\"time_s\":{},\"steps\":{},\"total_solid_solves\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"temperature_limit_k\":{},\"first_sampled_violation_s\":{},\"stored_energy_change_j\":{},\"input_energy_j\":{},\"air_energy_gain_j\":{},\"energy_residual_j\":{},\"history\":[{}],\"adaptive\":{},\"scope\":\"initial state and accepted endpoints only; no inter-step peak/crossing certificate, transient adjoint, air inertia, ramp model, or time-discretization error bound; solid_inputs are base declarations; each interval selects a global power scale or absolute named component watts\"}}}}\n",
        num(time)?,completed,work,num(peak)?,num(peak_time)?,optional(schedule.limit)?,optional(first_violation)?,
        num(stored)?,num(input)?,num(exhaust)?,num(residual)?,history.join(","),adaptive_stats.render(schedule.adaptive)?);
    poll(cx)?;
    Ok(Trajectory {output,peak_k:peak,peak_time_s:peak_time,solid_solves:work,steps:completed})
}

#[cfg(test)]
mod tests;
