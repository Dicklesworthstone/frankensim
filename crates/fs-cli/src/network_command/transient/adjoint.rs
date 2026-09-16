//! Fixed-grid discrete trajectory adjoint. Only accepted fields/references
//! are retained; matrices and coupled adjoints are reconstructed one endpoint
//! at a time. No primal iteration or discarded adaptive trial is differentiated.
use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Observable { Final, SampledPeak }

#[derive(Debug, Clone, Copy)]
pub(super) struct Config {
    observable: Observable,
    max_checkpoint_bytes: usize,
}
impl Config {
    pub(super) fn parse(value: &J) -> Result<Self> {
        object(value, &["qoi", "max_checkpoint_bytes"], "transient.adjoint")?;
        let observable = match get(value, "qoi")?.as_str() {
            Some("final") => Observable::Final,
            Some("sampled-peak") => Observable::SampledPeak,
            _ => return Err(bad("transient.adjoint.qoi must be final or sampled-peak")),
        };
        Ok(Self { observable, max_checkpoint_bytes: count(get(value,"max_checkpoint_bytes")?,
            "adjoint.max_checkpoint_bytes", 512 * 1024 * 1024)? })
    }
}

struct Frame {
    temperature: Vec<f64>,
    references: Vec<f64>,
    time: f64,
    dt: f64,
    interval: usize,
    objective: f64,
    vertex: Option<usize>,
}

pub(super) struct Tape {
    config: Config,
    frames: Vec<Frame>,
    planned: usize,
    vertices: usize,
    regions: usize,
    charged_bytes: usize,
    peak: f64,
    selected_peak: usize,
    initial_vertex: Option<usize>,
}

impl Tape {
    pub(super) fn new(request: &Request, schedule: &Schedule, initial: f64,
        initial_vertex: Option<usize>) -> Result<Option<Self>> {
        let Some(config) = schedule.adjoint else { return Ok(None); };
        if schedule.adaptive.is_some() || schedule.repeat.is_some()
            || schedule.power_design.is_some() || schedule.fan_speed_design.is_some() {
            return Err(bad("transient adjoints require a fixed, nonrepeated schedule without nested design searches"));
        }
        let vertices = request.mesh.vertex_count();
        let regions = request.surfaces.len();
        let charged_bytes = vertices.checked_add(regions)
            .and_then(|n| n.checked_mul(std::mem::size_of::<f64>()))
            .and_then(|n| n.checked_add(std::mem::size_of::<Frame>()))
            .and_then(|n| n.checked_mul(schedule.total_steps))
            .ok_or_else(|| budget("transient adjoint checkpoint size overflow"))?;
        if charged_bytes > config.max_checkpoint_bytes {
            return Err(budget("accepted-field adjoint checkpoints exceed max_checkpoint_bytes"));
        }
        let mut frames = Vec::new();
        frames.try_reserve_exact(schedule.total_steps)
            .map_err(|_| budget("cannot allocate admitted transient adjoint checkpoints"))?;
        Ok(Some(Self { config, frames, planned: schedule.total_steps, vertices, regions,
            charged_bytes, peak: initial, selected_peak: 0, initial_vertex }))
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn record(&mut self, temperature: &[f64], references: &[f64], time: f64,
        dt: f64, interval: usize, objective: f64, vertex: Option<usize>) -> Result<()> {
        if self.frames.len() == self.planned || temperature.len() != self.vertices
            || references.len() != self.regions {
            return Err(bad("accepted endpoint does not match the admitted adjoint tape"));
        }
        if objective > self.peak { self.peak = objective; self.selected_peak = self.frames.len() + 1; }
        self.frames.push(Frame { temperature: temperature.to_vec(), references: references.to_vec(),
            time, dt, interval, objective, vertex });
        Ok(())
    }

    pub(super) fn reverse(self, request: &Request, cx: &Cx<'_>, schedule: &Schedule,
        engine: &BackwardEuler<'_>) -> Result<(String, usize)> {
        poll(cx)?;
        if self.frames.len() != self.planned { return Err(bad("incomplete trajectory has no adjoint")); }
        let selected = match self.config.observable {
            Observable::Final => self.frames.len(), Observable::SampledPeak => self.selected_peak,
        };
        let mut carry = vec![0.0; self.vertices];
        let mut powers = vec![0.0; schedule.intervals.len()];
        let mut fan_speeds = vec![request.fan.as_ref().map(|_| 0.0); schedule.intervals.len()];
        let mut inlets = vec![0.0; request.graph.node_count()];
        let mut capacity = 0.0;
        let mut reconstructed = 0_usize;
        let mut adjoint_sweeps = 0_usize;
        let mut worst_residual = 0.0_f64;
        if selected == 0 {
            seed_initial(request, cx, self.initial_vertex, &mut carry)?;
        } else {
            let material = fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
            for ordinal in (0..=self.frames[selected-1].interval).rev() {
                poll(cx)?;
                let interval = &schedule.intervals[ordinal];
                let flow = match (&request.fan, interval.speed) {
                    (Some(fan), Some(speed)) => fan.solve(cx, &request.graph, request.limits, speed)?,
                    (None, None) => request.flow(cx)?,
                    _ => return Err(bad("transient adjoint drive mismatch")),
                };
                let base = request.surfaces.iter().map(|s| (s.name.clone(),s.h)).collect();
                let (coefficients, derived) = convection::resolve(request,cx,&flow,&base)?;
                let network = request.transport(cx,&flow,&coefficients)?;
                let names = network.regions();
                let load = interval.workload.prepare(request,cx)?;
                let gate = ConjugateConfig { max_iterations: request.limits.coupling,
                    temperature_tolerance_k: request.limits.temperature,
                    balance_tolerance_w: request.limits.heat/(names.len() as f64+1.0),
                    balance_relative_tolerance: 0.0,
                    relaxation: Relaxation::Fixed { omega: request.limits.relaxation } };
                for index in (0..selected).rev().filter(|&i| self.frames[i].interval == ordinal) {
                    poll(cx)?;
                    let frame = &self.frames[index];
                    let old = if index == 0 { &schedule.initial } else { &self.frames[index-1].temperature };
                    let boundary = request.boundary(&names,&frame.references,&coefficients)?;
                    let config = StepConfig { linear: fs_conduction::LinearConfig {
                        tolerance: request.limits.relative, max_iterations: request.limits.linear, restart: 60 },
                        energy_tolerance_j: finite(request.limits.heat*frame.dt)? };
                    let step = engine.linearize_step(cx,ConductionProblem { mesh: &request.mesh,
                        boundary: &boundary, material: &material,
                        element_materials: request.solid_data.element_materials.as_ref(), source: &load.source },
                        request.contacts.as_ref().map(|c| &c.interfaces), old, frame.dt, config,
                        schedule.nonlinear.map(|c| c.policy), &names).map_err(producer)?;
                    // An operator at a nearby field is not the adjoint of the
                    // retained trajectory. Re-execute the same endpoint only,
                    // then insist on exact same-profile temperature replay.
                    if step.primal().temperature.iter().zip(&frame.temperature)
                        .any(|(a,b)| a.to_bits() != b.to_bits()) {
                        return Err(producer("transient adjoint endpoint reconstruction changed accepted temperature bits"));
                    }
                    reconstructed += 1;
                    let binding = CoupledLinearization::new(cx,&network,&step,&gate).map_err(producer)?;
                    let mut weights = binding.zero_objective();
                    weights.nodal_temperatures.clone_from(&carry);
                    if index + 1 == selected {
                        // The port vectors use NETWORK order, which can differ
                        // from the boundary declaration order in the report.
                        let states = names.iter().map(|name| step.primal().robin_fluxes.iter()
                            .find(|flux| flux.region == *name)
                            .map(SolidRegionState::from_robin_flux)
                            .ok_or_else(|| bad("reconstructed endpoint lacks an objective port")))
                            .collect::<Result<Vec<_>>>()?;
                        let objective = request.objective.evaluate(cx,&frame.temperature,&states)?;
                        if objective.vertex != frame.vertex || objective.value.to_bits() != frame.objective.to_bits() {
                            return Err(producer("transient adjoint objective branch changed during reconstruction"));
                        }
                        objective.seed(&mut weights);
                    }
                    // The existing fan chain rule also applies to an endpoint
                    // response containing storage, at fixed OLD temperature.
                    let gradient = fan_gradient::pullback(request,cx,&binding,&weights,&names,&derived)?;
                    powers[ordinal] = finite(powers[ordinal]
                        + step.source_multiplier_pullback(cx,&gradient.nodal_load).map_err(producer)?)?;
                    capacity = finite(capacity
                        + step.capacity_multiplier_pullback(cx,&gradient.nodal_load).map_err(producer)?)?;
                    fan_speeds[ordinal] = match (fan_speeds[ordinal],gradient.log_speed()) {
                        (Some(sum),Some(value)) => Some(finite(sum+value)?), _ => None,
                    };
                    if inlets.len() != gradient.inlets.len() { return Err(bad("adjoint inlet arity changed")); }
                    for (sum,value) in inlets.iter_mut().zip(&gradient.inlets) { *sum = finite(*sum+value)?; }
                    carry = step.previous_temperature_pullback(cx,&gradient.nodal_load).map_err(producer)?;
                    adjoint_sweeps = adjoint_sweeps.checked_add(gradient.iterations)
                        .ok_or_else(|| budget("transient adjoint sweep count overflow"))?;
                    worst_residual = worst_residual.max(gradient.interface_residual);
                }
            }
        }
        let uniform_initial = carry.iter().try_fold(0.0, |sum,value| finite(sum+value))?;
        let rows = powers.iter().zip(&fan_speeds).enumerate().map(|(i,(power,speed))| Ok(format!(
            "{{\"interval\":{i},\"dtemperature_dpower_multiplier_k\":{},\"dtemperature_dlog_fan_speed_ratio_k\":{}}}",
            num(*power)?, optional(*speed)?))).collect::<Result<Vec<_>>>()?.join(",");
        let (time,value,vertex) = if selected == 0 { (0.0,self.peak,self.initial_vertex) }
            else { let frame = &self.frames[selected-1]; (frame.time,frame.objective,frame.vertex) };
        poll(cx)?;
        let report = format!(
            "{{\"method\":\"discrete-backward-euler-coupled-adjoint\",\"qoi\":{},\"value_k\":{},\"time_s\":{},\"state_index\":{},\"active_vertex\":{},\"dtemperature_dinitial_temperatures\":{},\"dtemperature_duniform_initial_k\":{},\"dtemperature_dcapacity_multiplier_k\":{},\"dtemperature_dinlet_temperatures\":{},\"intervals\":[{}],\"checkpoint_bytes\":{},\"reconstructed_solid_endpoints\":{},\"adjoint_sweeps\":{},\"max_interface_residual\":{},\"scope\":\"fixed accepted time grid; selected final or earliest sampled-maximum branch, not a continuous-time maximum or unique derivative at ties; full storage, material K-prime, contact and mixed-air feedback; interval power multiplies its entire declared load at multiplier one, capacity multiplies the complete matrix at one; fan controls include single-bank affinity and supported convection response, null when unavailable; fixed geometry/material laws/contact resistance/fluid properties; derivative and linear iteration budgets apply per reverse endpoint under the original wall deadline; checkpoint bytes bound retained fields/references and frame storage, not total solver workspace; separate from the null steady-gradient fields\"}}",
            quote(match self.config.observable { Observable::Final => "final", Observable::SampledPeak => "sampled-peak" }),
            num(value)?,num(time)?,selected,vertex.map_or_else(||"null".into(),|v|v.to_string()),
            numbers(&carry)?,num(uniform_initial)?,num(capacity)?,numbers(&inlets)?,rows,
            self.charged_bytes,reconstructed,adjoint_sweeps,num(worst_residual)?,
        );
        Ok((report, reconstructed))
    }
}

fn seed_initial(request: &Request, cx: &Cx<'_>, vertex: Option<usize>, weights: &mut [f64]) -> Result<()> {
    if let objective::Objective::MeanWall(name) = &request.objective {
        let surface = request.surfaces.iter().find(|s| &s.name == name)
            .ok_or_else(|| bad("missing initial objective surface"))?;
        for face in request.mesh.boundary() {
            poll(cx)?;
            if surface.faces.contains(&face.vertices) {
                for v in face.vertices { weights[v as usize] = finite(weights[v as usize]+face.area/surface.area/3.0)?; }
            }
        }
    } else {
        let vertex = vertex.filter(|&v|v<weights.len()).ok_or_else(||bad("missing initial objective vertex"))?;
        weights[vertex] = 1.0;
    }
    Ok(())
}
