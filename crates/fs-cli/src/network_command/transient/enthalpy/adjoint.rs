//! Fixed-grid physical h-history transpose with complete endpoint air/radiation
//! feedback. Replay checks both h and T; no temperature-to-enthalpy inverse is
//! used for reverse carry, including inside latent plateaus.
use super::super::adjoint::{Config as AdjointConfig, Observable, controls, seed_initial};
use super::*;
use fs_airflow::graph::thermal::coupled_transport::sensitivity::enthalpy::CoupledEnthalpyLinearization;
use fs_conduction::transient::enthalpy::adjoint::EnthalpyStepLinearization;

struct Frame {
    h: Vec<f64>,
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
    // Cycle starts evaluate the field directly, while accepted endpoints use
    // the producer's convective quadrature for a wall-mean objective.
    boundary: bool,
}

pub(super) struct Tape {
    config: AdjointConfig,
    frames: Vec<Frame>,
    planned: usize,
    cycles: usize,
    vertices: usize,
    regions: usize,
    charged_bytes: usize,
    peak: Selection,
}

impl Tape {
    pub(super) fn new(
        request: &Request,
        schedule: &Schedule,
        initial: f64,
        vertex: Option<usize>,
    ) -> Result<Option<Self>> {
        if schedule.adjoint.is_some() && schedule.repeat.is_some() {
            return Err(bad(
                "repeated enthalpy adjoints require the complete duty-cycle tape",
            ));
        }
        Self::for_cycles(request, schedule, initial, vertex, 1)
    }

    pub(super) fn for_cycles(
        request: &Request,
        schedule: &Schedule,
        initial: f64,
        vertex: Option<usize>,
        cycles: usize,
    ) -> Result<Option<Self>> {
        let Some(config) = schedule.adjoint else {
            return Ok(None);
        };
        if schedule.adaptive.is_some()
            || schedule.power_design.is_some()
            || schedule.fan_speed_design.is_some()
            || cycles == 0
        {
            return Err(bad(
                "enthalpy adjoints require fixed timesteps and an admitted cycle count without nested design searches",
            ));
        }
        if let Some(repeat) = schedule.repeat {
            repeat.validate_adjoint()?;
        }
        let vertices = request.mesh.vertex_count();
        let regions = request.surfaces.len();
        let control_bytes = config.controls.storage_bytes(request, schedule)?;
        let planned = schedule
            .total_steps
            .checked_mul(cycles)
            .ok_or_else(|| budget("enthalpy adjoint endpoint count overflow"))?;
        let frame_bytes = vertices
            .checked_mul(2)
            .and_then(|n| n.checked_add(regions))
            .and_then(|n| n.checked_mul(std::mem::size_of::<f64>()))
            .and_then(|n| n.checked_add(std::mem::size_of::<Frame>()))
            .and_then(|n| n.checked_mul(planned));
        let accumulator_bytes = vertices
            .checked_add(schedule.intervals.len())
            .and_then(|n| n.checked_add(request.graph.node_count()))
            .and_then(|n| n.checked_mul(std::mem::size_of::<f64>()))
            .and_then(|n| n.checked_add(3 * std::mem::size_of::<Vec<f64>>()))
            .and_then(|n| {
                schedule
                    .intervals
                    .len()
                    .checked_mul(std::mem::size_of::<Option<f64>>())
                    .and_then(|fan| n.checked_add(fan))
            })
            .and_then(|n| n.checked_add(std::mem::size_of::<Vec<Option<f64>>>()));
        let charged_bytes = frame_bytes
            .zip(accumulator_bytes)
            .and_then(|(a, b)| a.checked_add(b))
            .and_then(|n| n.checked_add(control_bytes))
            .ok_or_else(|| budget("enthalpy adjoint checkpoint size overflow"))?;
        if charged_bytes > config.max_checkpoint_bytes {
            return Err(budget(
                "enthalpy adjoint h/temperature/reference checkpoints and requested controls exceed max_checkpoint_bytes",
            ));
        }
        let mut frames = Vec::new();
        frames
            .try_reserve_exact(planned)
            .map_err(|_| budget("cannot allocate admitted enthalpy adjoint checkpoints"))?;
        Ok(Some(Self {
            config,
            frames,
            planned,
            cycles,
            vertices,
            regions,
            charged_bytes,
            peak: Selection {
                state: 0,
                time: 0.0,
                value: finite(initial)?,
                vertex,
                boundary: true,
            },
        }))
    }

    pub(super) fn begin_cycle(
        &mut self,
        initial: f64,
        vertex: Option<usize>,
        time: f64,
    ) -> Result<()> {
        finite(initial)?;
        finite(time)?;
        if initial > self.peak.value {
            self.peak = Selection {
                state: self.frames.len(),
                time,
                value: initial,
                vertex,
                boundary: true,
            };
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn record(
        &mut self,
        h: &[f64],
        temperature: &[f64],
        references: &[f64],
        time: f64,
        dt: f64,
        interval: usize,
        objective: f64,
        vertex: Option<usize>,
    ) -> Result<()> {
        if self.frames.len() >= self.planned
            || h.len() != self.vertices
            || temperature.len() != self.vertices
            || references.len() != self.regions
            || !dt.is_finite()
            || dt <= 0.0
        {
            return Err(bad(
                "accepted endpoint does not match the admitted enthalpy adjoint tape",
            ));
        }
        finite(time)?;
        finite(objective)?;
        if objective > self.peak.value {
            self.peak = Selection {
                state: self.frames.len() + 1,
                time,
                value: objective,
                vertex,
                boundary: false,
            };
        }
        self.frames.push(Frame {
            h: h.to_vec(),
            temperature: temperature.to_vec(),
            references: references.to_vec(),
            time,
            dt,
            interval,
            objective,
            vertex,
        });
        Ok(())
    }

    pub(super) fn reverse(
        self,
        request: &Request,
        cx: &Cx<'_>,
        schedule: &Schedule,
        config: &Config,
        engine: &Prepared<'_, '_>,
    ) -> Result<(String, usize, Option<design_sensitivity::DesignSensitivity>)> {
        poll(cx)?;
        if self.frames.len() != self.planned {
            return Err(bad("incomplete enthalpy trajectory has no adjoint"));
        }
        let selected = match self.config.observable {
            Observable::SampledPeak => self.peak,
            Observable::Final => {
                let frame = self
                    .frames
                    .last()
                    .ok_or_else(|| bad("empty enthalpy trajectory"))?;
                Selection {
                    state: self.frames.len(),
                    time: frame.time,
                    value: frame.objective,
                    vertex: frame.vertex,
                    boundary: false,
                }
            }
        };
        let mut carry = vec![0.0; self.vertices];
        let mut powers = vec![0.0; schedule.intervals.len()];
        let mut fan_speeds = vec![request.fan.as_ref().map(|_| 0.0); schedule.intervals.len()];
        let mut inlets = vec![0.0; request.graph.node_count()];
        let mut controls = controls::Accumulation::new(self.config.controls, request, schedule)?;
        let mut reconstruction = Work::default();
        let mut reconstructed = 0_usize;
        let mut sweeps = 0_usize;
        let mut krylov = 0_usize;
        let mut worst_residual = 0.0_f64;
        if selected.state == 0 {
            seed_initial(request, cx, selected.vertex, &mut carry)?;
            for (vertex, weight) in carry.iter_mut().enumerate() {
                poll(cx)?;
                if *weight == 0.0 {
                    continue;
                }
                let curve = &config.materials[config.vertex_material_ids[vertex]].curve;
                let h = config.initial_h[vertex];
                check_initial_branch(curve, h)?;
                *weight = finite(
                    *weight
                        * curve
                            .temperature_derivative_at_specific_enthalpy(h)
                            .map_err(producer)?,
                )?;
            }
        } else {
            let material =
                fs_conduction::ConductivityModel::isotropic_declared(request.conductivity)
                    .map_err(producer)?;
            let radiation = request
                .radiation
                .as_ref()
                .map(|p| p.enthalpy_controls(request))
                .transpose()?;
            let mut end = selected.state;
            while end > 0 {
                poll(cx)?;
                let ordinal = self.frames[end - 1].interval;
                let interval = schedule
                    .intervals
                    .get(ordinal)
                    .ok_or_else(|| bad("enthalpy checkpoint names an unknown interval"))?;
                let mut start = end - 1;
                while start > 0 && self.frames[start - 1].interval == ordinal {
                    start -= 1;
                }
                let flow = match (&request.fan, interval.speed) {
                    (Some(fan), Some(speed)) => {
                        fan.solve(cx, &request.graph, request.limits, speed)?
                    }
                    (None, None) => request.flow(cx)?,
                    _ => return Err(bad("enthalpy adjoint drive mismatch")),
                };
                let base = request
                    .surfaces
                    .iter()
                    .map(|s| (s.name.clone(), s.h))
                    .collect();
                let (coefficients, derived) = convection::resolve(request, cx, &flow, &base)?;
                let network = request.transport(cx, &flow, &coefficients)?;
                let names = network.regions();
                let load = interval.workload.prepare(request, cx)?;
                let gate = ConjugateConfig {
                    max_iterations: request.limits.coupling,
                    temperature_tolerance_k: request.limits.temperature,
                    balance_tolerance_w: request.limits.heat / (names.len() as f64 + 1.0),
                    balance_relative_tolerance: 0.0,
                    relaxation: Relaxation::Fixed {
                        omega: request.limits.relaxation,
                    },
                };
                for index in (start..end).rev() {
                    poll(cx)?;
                    let frame = &self.frames[index];
                    let old = if index == 0 {
                        &config.initial_h
                    } else {
                        &self.frames[index - 1].h
                    };
                    let boundary = request.boundary(&names, &frame.references, &coefficients)?;
                    let problem = ConductionProblem {
                        mesh: &request.mesh,
                        boundary: &boundary,
                        material: &material,
                        element_materials: request.solid_data.element_materials.as_ref(),
                        source: &load.source,
                    };
                    let interfaces = request.contacts.as_ref().map(|c| &c.interfaces);
                    let step_config = EnthalpyStepConfig {
                        newton: config.newton,
                        max_newton_iterations: config.max_iterations,
                        energy_tolerance_j: finite(request.limits.heat * frame.dt)?,
                    };
                    let (accepted, objective_state) = if let Some((patches, policy)) = &radiation {
                        let solved = engine
                            .advance_with_ambient_radiation(
                                cx,
                                problem,
                                interfaces,
                                old,
                                frame.dt,
                                step_config,
                                patches,
                                *policy,
                            )
                            .map_err(producer)?;
                        reconstruction.add(
                            solved.radiation.iterations,
                            solved.radiation.solid_iterations,
                            solved.radiation.krylov_iterations,
                        )?;
                        let state = replay_objective(
                            request,
                            cx,
                            frame,
                            &names,
                            &solved.convective_robin_fluxes,
                            index + 1 == selected.state && !selected.boundary,
                        )?;
                        (solved.conduction, state)
                    } else {
                        let solved = engine
                            .advance(cx, problem, interfaces, old, frame.dt, step_config)
                            .map_err(producer)?;
                        let columns = solved
                            .newton
                            .history
                            .iter()
                            .try_fold(0_usize, |sum, row| sum.checked_add(row.linear_iterations))
                            .ok_or_else(|| budget("enthalpy replay work overflow"))?;
                        reconstruction.add(1, solved.newton.iterations, columns)?;
                        let state = replay_objective(
                            request,
                            cx,
                            frame,
                            &names,
                            &solved.robin_fluxes,
                            index + 1 == selected.state && !selected.boundary,
                        )?;
                        (solved, state)
                    };
                    check_replay(frame, &accepted)?;
                    // Both owners expose the same checked physical tangent. Keep
                    // the owning value alive while the convection response borrows it.
                    let (radiative, ordinary) = if let Some((patches, _)) = &radiation {
                        let entries = self
                            .vertices
                            .checked_mul(patches.len())
                            .and_then(|n| n.checked_mul(2))
                            .ok_or_else(|| budget("enthalpy radiation feedback size overflow"))?;
                        (
                            Some(
                                engine
                                    .linearize_accepted_with_ambient_radiation(
                                        cx,
                                        problem,
                                        interfaces,
                                        old,
                                        frame.dt,
                                        step_config,
                                        patches,
                                        accepted,
                                        entries,
                                    )
                                    .map_err(producer)?,
                            ),
                            None,
                        )
                    } else {
                        (
                            None,
                            Some(
                                engine
                                    .linearize_accepted(
                                        cx,
                                        problem,
                                        interfaces,
                                        old,
                                        frame.dt,
                                        step_config,
                                        accepted,
                                    )
                                    .map_err(producer)?,
                            ),
                        )
                    };
                    let step: &EnthalpyStepLinearization<'_> = match (&radiative, &ordinary) {
                        (Some(step), None) => step.transport(),
                        (None, Some(step)) => step.transport(),
                        _ => return Err(bad("missing enthalpy endpoint tangent")),
                    };
                    let response = step
                        .robin_response(
                            cx,
                            &names,
                            fs_conduction::LinearConfig {
                                tolerance: request.limits.relative,
                                max_iterations: request.limits.linear,
                                restart: 60,
                            },
                        )
                        .map_err(producer)?;
                    let binding = CoupledEnthalpyLinearization::new(cx, &network, &response, &gate)
                        .map_err(producer)?;
                    let mut objective = binding.zero_objective();
                    if index + 1 == selected.state && selected.boundary {
                        let (value, vertex) = initial_objective(request, cx, &frame.temperature)?;
                        if vertex != selected.vertex || value.to_bits() != selected.value.to_bits()
                        {
                            return Err(producer(
                                "enthalpy cycle-boundary objective changed during adjoint reconstruction",
                            ));
                        }
                        seed_initial(request, cx, vertex, &mut objective.nodal_temperatures)?;
                    } else if let Some(state) = objective_state {
                        state.seed(&mut objective);
                    }
                    let interface_config = InterfaceSolveConfig {
                        max_iterations: request.limits.derivative,
                        absolute_tolerance: request.limits.relative,
                        relative_tolerance: request.limits.relative,
                        relaxation: request.limits.relaxation,
                    };
                    let (gradient, speed) = if request.fan.is_some() {
                        let response = binding
                            .pullback_flow_scale_iqn(
                                cx,
                                &objective,
                                &carry,
                                interface_config,
                                acceleration::POLICY,
                            )
                            .map_err(producer)?;
                        let speed = fan_gradient::log_speed_from_flow_response(
                            cx,
                            &names,
                            &derived,
                            &response.thermal.log_htc,
                            response.log_flow_scale,
                        )?;
                        (response.thermal, speed)
                    } else {
                        (
                            binding
                                .pullback_iqn(
                                    cx,
                                    &objective,
                                    &carry,
                                    interface_config,
                                    acceleration::POLICY,
                                )
                                .map_err(producer)?,
                            None,
                        )
                    };
                    fan_speeds[ordinal] = match (fan_speeds[ordinal], speed) {
                        (Some(sum), Some(value)) => Some(finite(sum + value)?),
                        _ => None,
                    };
                    for (vertex, value) in gradient.solid.source_density.iter().enumerate() {
                        if vertex % 512 == 0 {
                            poll(cx)?;
                        }
                        powers[ordinal] = finite(powers[ordinal] + value * load.source.at(vertex))?;
                    }
                    if gradient.inlets.len() != inlets.len() {
                        return Err(bad("enthalpy inlet derivative arity changed"));
                    }
                    for (sum, value) in inlets.iter_mut().zip(&gradient.inlets) {
                        *sum = finite(*sum + value)?;
                    }
                    controls.record_pullback(
                        request,
                        cx,
                        &frame.temperature,
                        ordinal,
                        &gradient.solid.source_density,
                        &gradient.nodal_load,
                    )?;
                    carry = gradient.solid.previous_specific_enthalpy;
                    reconstructed += 1;
                    sweeps = sweeps
                        .checked_add(gradient.iterations)
                        .ok_or_else(|| budget("enthalpy adjoint sweep count overflow"))?;
                    krylov = krylov
                        .checked_add(gradient.solid_krylov_iterations)
                        .ok_or_else(|| budget("enthalpy adjoint Krylov count overflow"))?;
                    worst_residual = worst_residual.max(gradient.interface_residual);
                }
                end = start;
            }
        }
        let uniform = carry
            .iter()
            .try_fold(0.0, |sum, value| finite(sum + value))?;
        let rows = powers.iter().zip(&fan_speeds).enumerate().map(|(i,(power,speed))| Ok(format!(
            "{{\"interval\":{i},\"dtemperature_dpower_multiplier_k\":{},\"dtemperature_dlog_fan_speed_ratio_k\":{}}}", num(*power)?, optional(*speed)?)))
            .collect::<Result<Vec<_>>>()?.join(",");
        let design = if self.config.observable == Observable::SampledPeak {
            Some(design_sensitivity::DesignSensitivity::from_intervals(
                selected.value,
                &powers,
                &fan_speeds,
            )?)
        } else {
            None
        };
        let controls_fragment = controls.report_fragment(request, cx, schedule)?;
        let report = format!(
            "{{\"method\":\"discrete-backward-euler-coupled-enthalpy-adjoint\",\"qoi\":{},\"value_k\":{},\"time_s\":{},\"state_index\":{},\"active_vertex\":{},\"cycles\":{},\"dtemperature_dinitial_specific_enthalpies_k_kg_j\":{},\"dtemperature_duniform_initial_specific_enthalpy_k_kg_j\":{},\"dtemperature_dinlet_temperatures\":{},\"intervals\":[{}],\"checkpoint_bytes\":{},\"reconstructed_solid_endpoints\":{},\"reconstruction_solid_solves\":{},\"adjoint_sweeps\":{},\"adjoint_krylov_iterations\":{},\"max_interface_residual\":{}{controls_fragment},\"scope\":\"fixed accepted grid and cycle count; final or earliest sampled-maximum branch in global time, no continuous peak bound or unique derivative at ties; full specific-enthalpy history across every cycle, contact, mixed-air and declared ambient-radiation feedback; shared interval controls affect every repeated occurrence and multiply their actual source without division by baseline power; initial h derivatives describe the original initial field and have units K kg/J; inlet controls apply throughout; fan controls include single-bank affinity, capacity transport and supported convection Reynolds response, null when unavailable; separately requested component and contact controls report their own units and scope; fixed reference densities, charts, conductivity/contact laws, geometry, fluid properties and radiation controls; slope corners and validity endpoints refuse classical endpoint derivatives; replay verifies accepted h and T bits and rechecks physical residual/energy gates; checkpoint bytes bound retained h/T/references for all cycles and control accumulators, not total workspace; per-endpoint derivative budgets share the original wall deadline\"}}",
            quote(match self.config.observable {
                Observable::Final => "final",
                Observable::SampledPeak => "sampled-peak",
            }),
            num(selected.value)?,
            num(selected.time)?,
            selected.state,
            selected
                .vertex
                .map_or_else(|| "null".into(), |v| v.to_string()),
            self.cycles,
            numbers(&carry)?,
            num(uniform)?,
            numbers(&inlets)?,
            rows,
            self.charged_bytes,
            reconstructed,
            reconstruction.solves,
            sweeps,
            krylov,
            num(worst_residual)?
        );
        poll(cx)?;
        Ok((report, reconstruction.solves, design))
    }
}

fn replay_objective(
    request: &Request,
    cx: &Cx<'_>,
    frame: &Frame,
    names: &[&str],
    convection: &[fs_conduction::RobinFlux],
    selected: bool,
) -> Result<Option<objective::ObjectiveState>> {
    if !selected {
        return Ok(None);
    }
    // Preserve the producer's quadrature order when verifying a wall mean.
    // Normalized port sums can differ from that mean in their last bit.
    let states = names
        .iter()
        .map(|name| {
            convection
                .iter()
                .find(|row| row.region == *name)
                .map(SolidRegionState::from_robin_flux)
                .ok_or_else(|| bad("enthalpy replay lacks an objective convection port"))
        })
        .collect::<Result<Vec<_>>>()?;
    let state = request
        .objective
        .evaluate(cx, &frame.temperature, &states)?;
    if state.vertex != frame.vertex || state.value.to_bits() != frame.objective.to_bits() {
        return Err(producer(
            "enthalpy adjoint objective branch changed during reconstruction",
        ));
    }
    Ok(Some(state))
}

fn check_replay(frame: &Frame, accepted: &EnthalpyStepSolution) -> Result<()> {
    for (old, replayed) in [
        (&frame.h, &accepted.specific_enthalpy_j_kg),
        (&frame.temperature, &accepted.temperature),
    ] {
        if old.len() != replayed.len()
            || old
                .iter()
                .zip(replayed)
                .any(|(a, b)| a.to_bits() != b.to_bits())
        {
            return Err(producer(
                "enthalpy adjoint reconstruction changed accepted h or temperature bits",
            ));
        }
    }
    Ok(())
}

fn check_initial_branch(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> Result<()> {
    let knots = curve.knots();
    if let Some(index) = knots.iter().position(|k| k.specific_enthalpy_j_kg == h) {
        if index == 0 || index + 1 == knots.len() {
            return Err(producer(
                "initial enthalpy objective lies at a chart validity endpoint",
            ));
        }
        let left = curve
            .temperature_derivative_at_specific_enthalpy(knots[index - 1].specific_enthalpy_j_kg)
            .map_err(producer)?;
        let right = curve
            .temperature_derivative_at_specific_enthalpy(h)
            .map_err(producer)?;
        if left != right {
            return Err(producer(
                "initial enthalpy objective lies at a chart slope corner",
            ));
        }
    }
    Ok(())
}
