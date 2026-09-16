//! Fixed-grid discrete trajectory adjoint. Only accepted fields/references
//! are retained; matrices and coupled adjoints are reconstructed one endpoint
//! at a time. Repeated schedules share ONE chronological history: equal
//! interval labels in different cycles must never reorder the reverse sweep.
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

#[derive(Clone, Copy)]
struct Selection {
    state: usize,
    time: f64,
    value: f64,
    vertex: Option<usize>,
    // Cycle starts reuse the prior accepted field but evaluate their spatial
    // objective directly. Preserve that rounding/selection path at a boundary.
    boundary: bool,
}

pub(super) struct Tape {
    config: Config,
    frames: Vec<Frame>,
    planned: usize,
    cycles: usize,
    vertices: usize,
    regions: usize,
    charged_bytes: usize,
    peak: Selection,
}

impl Tape {
    pub(super) fn new(request: &Request, schedule: &Schedule, initial: f64,
        initial_vertex: Option<usize>) -> Result<Option<Self>> {
        if schedule.adjoint.is_some() && schedule.repeat.is_some() {
            return Err(bad("repeated adjoints require the complete duty-cycle tape"));
        }
        Self::for_cycles(request, schedule, initial, initial_vertex, 1)
    }

    /// The repeated driver admits a fixed count and rejects controllers or
    /// periodic stopping before calling. Charge ALL cycles before a PDE solve.
    pub(super) fn for_cycles(request: &Request, schedule: &Schedule, initial: f64,
        initial_vertex: Option<usize>, cycles: usize) -> Result<Option<Self>> {
        let Some(config) = schedule.adjoint else { return Ok(None); };
        if schedule.adaptive.is_some() || schedule.power_design.is_some()
            || schedule.fan_speed_design.is_some() || cycles == 0 {
            return Err(bad("transient adjoints require fixed timesteps and an admitted cycle count without nested design searches"));
        }
        let planned = schedule.total_steps.checked_mul(cycles)
            .ok_or_else(|| budget("transient adjoint endpoint count overflow"))?;
        let vertices = request.mesh.vertex_count();
        let regions = request.surfaces.len();
        let charged_bytes = vertices.checked_add(regions)
            .and_then(|n| n.checked_mul(std::mem::size_of::<f64>()))
            .and_then(|n| n.checked_add(std::mem::size_of::<Frame>()))
            .and_then(|n| n.checked_mul(planned))
            .ok_or_else(|| budget("transient adjoint checkpoint size overflow"))?;
        if charged_bytes > config.max_checkpoint_bytes {
            return Err(budget("accepted-field adjoint checkpoints exceed max_checkpoint_bytes"));
        }
        let mut frames = Vec::new();
        frames.try_reserve_exact(planned)
            .map_err(|_| budget("cannot allocate admitted transient adjoint checkpoints"))?;
        Ok(Some(Self { config, frames, planned, cycles, vertices, regions, charged_bytes,
            peak: Selection { state: 0, time: 0.0, value: finite(initial)?,
                vertex: initial_vertex, boundary: true } }))
    }

    /// Do not discard cycle-start samples used by the existing forward peak
    /// calculation. They have no new physical step and need no copied field.
    pub(super) fn begin_cycle(&mut self, initial: f64, vertex: Option<usize>, time: f64) -> Result<()> {
        finite(initial)?;
        finite(time)?;
        if initial > self.peak.value {
            self.peak = Selection { state: self.frames.len(), time, value: initial, vertex, boundary: true };
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn record(&mut self, temperature: &[f64], references: &[f64], time: f64,
        dt: f64, interval: usize, objective: f64, vertex: Option<usize>) -> Result<()> {
        if self.frames.len() == self.planned || temperature.len() != self.vertices
            || references.len() != self.regions || !dt.is_finite() || dt <= 0.0 {
            return Err(bad("accepted endpoint does not match the admitted adjoint tape"));
        }
        finite(time)?;
        finite(objective)?;
        if objective > self.peak.value {
            self.peak = Selection { state: self.frames.len() + 1, time,
                value: objective, vertex, boundary: false };
        }
        self.frames.push(Frame { temperature: temperature.to_vec(), references: references.to_vec(),
            time, dt, interval, objective, vertex });
        Ok(())
    }

    pub(super) fn reverse(self, request: &Request, cx: &Cx<'_>, schedule: &Schedule,
        engine: &BackwardEuler<'_>) -> Result<(String, usize)> {
        poll(cx)?;
        if self.frames.len() != self.planned { return Err(bad("incomplete trajectory has no adjoint")); }
        let selected = match self.config.observable {
            Observable::SampledPeak => self.peak,
            Observable::Final => {
                let frame = self.frames.last().ok_or_else(|| bad("empty trajectory has no final adjoint"))?;
                Selection { state: self.frames.len(), time: frame.time,
                    value: frame.objective, vertex: frame.vertex, boundary: false }
            }
        };
        let mut carry = vec![0.0; self.vertices];
        let mut powers = vec![0.0; schedule.intervals.len()];
        let mut fan_speeds = vec![request.fan.as_ref().map(|_| 0.0); schedule.intervals.len()];
        let mut inlets = vec![0.0; request.graph.node_count()];
        let mut capacity = 0.0;
        let mut reconstructed = 0_usize;
        let mut adjoint_sweeps = 0_usize;
        let mut worst_residual = 0.0_f64;
        if selected.state == 0 {
            seed_initial(request, cx, selected.vertex, &mut carry)?;
        } else {
            let material = fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
            let mut end = selected.state;
            while end > 0 {
                poll(cx)?;
                let ordinal = self.frames[end-1].interval;
                let interval = schedule.intervals.get(ordinal)
                    .ok_or_else(|| bad("retained endpoint names an unknown interval"))?;
                // Only CONTIGUOUS occurrences can share this reverse block.
                // Grouping all equal labels would reverse cycle 3/interval 1,
                // then cycle 2/interval 1 before cycle 3/interval 0: wrong history.
                let mut start = end - 1;
                while start > 0 && self.frames[start-1].interval == ordinal { start -= 1; }
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
                for index in (start..end).rev() {
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
                    if step.primal().temperature.len() != frame.temperature.len()
                        || step.primal().temperature.iter().zip(&frame.temperature)
                            .any(|(a,b)| a.to_bits() != b.to_bits()) {
                        return Err(producer("transient adjoint endpoint reconstruction changed accepted temperature bits"));
                    }
                    reconstructed += 1;
                    let binding = CoupledLinearization::new(cx,&network,&step,&gate).map_err(producer)?;
                    let mut weights = binding.zero_objective();
                    weights.nodal_temperatures.clone_from(&carry);
                    if index + 1 == selected.state {
                        if selected.boundary {
                            let (value, vertex) = initial_objective(request,cx,&frame.temperature)?;
                            if vertex != selected.vertex || value.to_bits() != selected.value.to_bits() {
                                return Err(producer("cycle-boundary objective changed during adjoint reconstruction"));
                            }
                            seed_initial(request,cx,vertex,&mut weights.nodal_temperatures)?;
                        } else {
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
                    }
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
                end = start;
            }
        }
        let uniform_initial = carry.iter().try_fold(0.0, |sum,value| finite(sum+value))?;
        let rows = powers.iter().zip(&fan_speeds).enumerate().map(|(i,(power,speed))| Ok(format!(
            "{{\"interval\":{i},\"dtemperature_dpower_multiplier_k\":{},\"dtemperature_dlog_fan_speed_ratio_k\":{}}}",
            num(*power)?, optional(*speed)?))).collect::<Result<Vec<_>>>()?.join(",");
        poll(cx)?;
        let report = format!(
            "{{\"method\":\"discrete-backward-euler-coupled-adjoint\",\"qoi\":{},\"value_k\":{},\"time_s\":{},\"state_index\":{},\"active_vertex\":{},\"cycles\":{},\"dtemperature_dinitial_temperatures\":{},\"dtemperature_duniform_initial_k\":{},\"dtemperature_dcapacity_multiplier_k\":{},\"dtemperature_dinlet_temperatures\":{},\"intervals\":[{}],\"checkpoint_bytes\":{},\"reconstructed_solid_endpoints\":{},\"adjoint_sweeps\":{},\"max_interface_residual\":{},\"scope\":\"fixed accepted time grid and cycle count; selected final or earliest all-cycle sampled-maximum branch, not a continuous-time maximum or unique derivative at ties; full storage, material K-prime, contact and mixed-air feedback across cycle boundaries; each interval control changes every occurrence of that base interval, including earlier warm-up cycles; initial controls change only the original initial field, capacity and inlets apply throughout; fan controls include single-bank affinity and supported convection response, null when unavailable; fixed geometry/material laws/contact resistance/fluid properties; derivative and linear iteration budgets apply per reverse endpoint under the original wall deadline; checkpoint bytes bound all retained fields/references and frame storage, not total solver workspace; separate from the null steady-gradient fields\"}}",
            quote(match self.config.observable { Observable::Final => "final", Observable::SampledPeak => "sampled-peak" }),
            num(selected.value)?,num(selected.time)?,selected.state,
            selected.vertex.map_or_else(||"null".into(),|v|v.to_string()),self.cycles,
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
        weights[vertex] = finite(weights[vertex] + 1.0)?;
    }
    Ok(())
}
