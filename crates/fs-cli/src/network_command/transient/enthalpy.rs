//! Fixed-grid total-enthalpy trajectories. Accepted h, never temperature, is
//! physical history. Air and radiation iterations borrow the same old h.
use super::*;
use fs_blake3::ContentHash;
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig, EnthalpyStepSolution,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_solver::{Globalization, LineSearchConfig, NewtonKrylovConfig};

#[derive(Debug)]
pub(super) struct Config {
    curve: EquilibriumEnthalpyPhaseCurve,
    source: String,
    density: f64,
    initial_h: Vec<f64>,
    newton: NewtonKrylovConfig,
    max_iterations: usize,
    max_backtracks: usize,
}

impl Config {
    pub(super) fn parse(value: &J, schedule: &J, vertices: usize) -> Result<Self> {
        for key in [
            "initial_temperature_k",
            "initial_temperatures_k",
            "volumetric_heat_capacity_j_m3_k",
            "element_heat_capacities_j_m3_k",
            "nonlinear",
            "adaptive",
            "repeat",
            "adjoint",
            "time_convergence",
            "fan_speed_design",
            "power_design",
        ] {
            if schedule.get(key).is_some() {
                return Err(bad(format!(
                    "transient.enthalpy does not admit transient.{key}; use explicit h and a fixed forward schedule"
                )));
            }
        }
        object(
            value,
            &[
                "material_card_identity",
                "source",
                "reference_density_kg_m3",
                "knots",
                "initial_specific_enthalpy_j_kg",
                "initial_specific_enthalpies_j_kg",
                "newton",
            ],
            "transient.enthalpy",
        )?;
        let identity = ContentHash::from_hex(&string(
            get(value, "material_card_identity")?,
            "material_card_identity",
        )?)
        .ok_or_else(|| bad("enthalpy material_card_identity requires 64 hex characters"))?;
        let source = string(get(value, "source")?, "enthalpy.source")?;
        let density = positive(
            get(value, "reference_density_kg_m3")?,
            "reference_density_kg_m3",
        )?;
        let mut knots = Vec::new();
        for row in array(get(value, "knots")?, "enthalpy.knots", 4096)? {
            object(
                row,
                &[
                    "specific_enthalpy_j_kg",
                    "temperature_k",
                    "liquid_mass_fraction",
                ],
                "enthalpy knot",
            )?;
            knots.push(EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: number(
                    get(row, "specific_enthalpy_j_kg")?,
                    "specific enthalpy",
                )?,
                temperature_k: positive(get(row, "temperature_k")?, "knot temperature")?,
                liquid_mass_fraction: number(
                    get(row, "liquid_mass_fraction")?,
                    "liquid mass fraction",
                )?,
                bulk_density_kg_m3: density,
            });
        }
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(identity, knots).map_err(producer)?;
        let initial_h = match (
            value.get("initial_specific_enthalpy_j_kg"),
            value.get("initial_specific_enthalpies_j_kg"),
        ) {
            (Some(h), None) => vec![number(h, "initial specific enthalpy")?; vertices],
            (None, Some(h)) => {
                let rows = array(h, "initial specific enthalpies", vertices)?;
                if rows.len() != vertices {
                    return Err(bad(
                        "one initial specific enthalpy per solid vertex required",
                    ));
                }
                rows.iter()
                    .map(|h| number(h, "initial specific enthalpy"))
                    .collect::<Result<Vec<_>>>()?
            }
            _ => {
                return Err(bad(
                    "choose exactly one uniform or nodal initial specific enthalpy",
                ));
            }
        };
        let policy = get(value, "newton")?;
        object(
            policy,
            &[
                "max_iterations",
                "residual_rtol",
                "residual_atol_j",
                "linear_restart",
                "max_linear_cycles",
                "armijo_c",
                "shrink",
                "max_backtracks",
            ],
            "enthalpy.newton",
        )?;
        let max_iterations = count(
            get(policy, "max_iterations")?,
            "enthalpy Newton attempts",
            10_000,
        )?;
        let max_backtracks = integer(get(policy, "max_backtracks")?, "enthalpy backtracks", 64)?;
        let armijo = positive(get(policy, "armijo_c")?, "enthalpy Armijo coefficient")?;
        let contraction = positive(get(policy, "shrink")?, "enthalpy line-search shrink")?;
        let absolute_tolerance =
            number(get(policy, "residual_atol_j")?, "enthalpy residual_atol_j")?;
        let relative_tolerance = positive(get(policy, "residual_rtol")?, "enthalpy residual_rtol")?;
        if armijo >= 1.0
            || contraction >= 1.0
            || absolute_tolerance < 0.0
            || relative_tolerance >= 1.0
        {
            return Err(bad(
                "enthalpy Newton requires Armijo, shrink and relative tolerance in (0,1), and nonnegative absolute tolerance",
            ));
        }
        let mut minimum_step = 1.0;
        for _ in 0..max_backtracks {
            minimum_step *= contraction;
        }
        if minimum_step <= 0.0 {
            return Err(bad(
                "enthalpy minimum line-search step is not representable",
            ));
        }
        let newton = NewtonKrylovConfig {
            absolute_tolerance,
            relative_tolerance,
            linear_restart: count(
                get(policy, "linear_restart")?,
                "enthalpy linear restart",
                vertices.min(256),
            )?,
            max_linear_cycles: count(
                get(policy, "max_linear_cycles")?,
                "enthalpy linear cycles",
                100_000,
            )?,
            globalization: Globalization::LineSearch(LineSearchConfig {
                armijo,
                contraction,
                minimum_step,
            }),
            ..NewtonKrylovConfig::default()
        };
        let config = Self {
            curve,
            source,
            density,
            initial_h,
            newton,
            max_iterations,
            max_backtracks,
        };
        config.initial_temperatures()?;
        Ok(config)
    }

    pub(super) fn initial_temperatures(&self) -> Result<Vec<f64>> {
        self.initial_h
            .iter()
            .map(|&h| {
                self.curve
                    .state_at_specific_enthalpy(h)
                    .map(|state| state.temperature_k())
                    .map_err(producer)
            })
            .collect()
    }

    pub(super) fn admit(&self, request: &Request, schedule: &Schedule) -> Result<()> {
        if request.gradient
            || request.design.is_some()
            || request.fan_speed_design.is_some()
            || request.mesh_convergence.is_some()
            || request.recirculation.is_some()
            || schedule.adaptive.is_some()
            || schedule.repeat.is_some()
            || schedule.adjoint.is_some()
            || schedule.time_convergence.is_some()
            || schedule.power_design.is_some()
            || schedule.fan_speed_design.is_some()
            || schedule.nonlinear.is_some()
        {
            return Err(bad(
                "enthalpy supports fixed forward schedules without design, adjoint, adaptive/repeated/study or recirculation modes",
            ));
        }
        let columns = self
            .max_iterations
            .checked_mul(self.newton.linear_restart)
            .and_then(|n| n.checked_mul(self.newton.max_linear_cycles))
            .ok_or_else(|| budget("enthalpy Krylov bound overflows"))?;
        if columns > request.limits.linear {
            return Err(budget(
                "enthalpy max_iterations * linear_restart * max_linear_cycles exceeds budgets.linear_iterations per inner solid solve",
            ));
        }
        if let Some(policy) = &request.radiation {
            policy.enthalpy_controls(request)?;
        }
        Ok(())
    }
}

#[derive(Default)]
struct Work {
    solves: usize,
    updates: usize,
    krylov: usize,
}
impl Work {
    fn add(&mut self, solves: usize, updates: usize, krylov: usize) -> Result<()> {
        for (slot, value) in [
            (&mut self.solves, solves),
            (&mut self.updates, updates),
            (&mut self.krylov, krylov),
        ] {
            *slot = slot
                .checked_add(value)
                .ok_or_else(|| budget("enthalpy work counter overflow"))?;
        }
        Ok(())
    }
}

struct Endpoint {
    coupled: CoupledTransportSolution,
    solid: EnthalpyStepSolution,
    radiation: Option<radiation::EndpointHeat>,
    physical_residual_j: f64,
    energy_residual_j: f64,
}

#[allow(clippy::too_many_arguments)]
fn advance(
    request: &Request,
    cx: &Cx<'_>,
    config: &Config,
    engine: &EnthalpyBackwardEuler<'_, '_>,
    network: &TransportNetwork<'_>,
    coefficients: &BTreeMap<String, f64>,
    old_h: &[f64],
    source: &ScalarField,
    dt: f64,
    work: &mut Work,
) -> Result<Endpoint> {
    let material = fs_conduction::ConductivityModel::isotropic_declared(request.conductivity)
        .map_err(producer)?;
    let names = network.regions();
    let coupling = ConjugateConfig {
        max_iterations: request.limits.coupling,
        temperature_tolerance_k: request.limits.temperature,
        balance_tolerance_w: request.limits.heat / (names.len() as f64 + 1.0),
        balance_relative_tolerance: 0.0,
        relaxation: Relaxation::Fixed {
            omega: request.limits.relaxation,
        },
    };
    let step_config = EnthalpyStepConfig {
        newton: config.newton,
        max_newton_iterations: config.max_iterations,
        energy_tolerance_j: finite(request.limits.heat * dt)?,
    };
    let controls = request
        .radiation
        .as_ref()
        .map(|policy| policy.enthalpy_controls(request))
        .transpose()?;
    let mut last = None;
    let mut failure = None;
    let coupled = solve_coupled_transport(cx, network, &coupling, |cx, references| {
        let result = (|| -> Result<Vec<SolidRegionState>> {
            let boundary = request.boundary(&names, references, coefficients)?;
            let problem = ConductionProblem {
                mesh: &request.mesh,
                boundary: &boundary,
                material: &material,
                element_materials: request.solid_data.element_materials.as_ref(),
                source,
            };
            let interfaces = request.contacts.as_ref().map(|c| &c.interfaces);
            // Neither a changed air reference nor a radiative trial advances h.
            let (solid, states, heat, residual, energy) = if let Some((patches, policy)) = &controls
            {
                let solved = engine
                    .advance_with_ambient_radiation(
                        cx,
                        problem,
                        interfaces,
                        old_h,
                        dt,
                        step_config,
                        patches,
                        *policy,
                    )
                    .map_err(producer)?;
                work.add(
                    solved.radiation.iterations,
                    solved.radiation.solid_iterations,
                    solved.radiation.krylov_iterations,
                )?;
                let heat = radiation::EndpointHeat::from_enthalpy(&solved.radiation)?;
                let states = names
                    .iter()
                    .map(|name| {
                        solved
                            .convective_robin_fluxes
                            .iter()
                            .find(|row| row.region == *name)
                            .map(SolidRegionState::from_robin_flux)
                            .ok_or_else(|| bad("enthalpy lacks a network convection region"))
                    })
                    .collect::<Result<Vec<_>>>()?;
                (
                    solved.conduction,
                    states,
                    Some(heat),
                    solved.physical_residual_norm_j,
                    solved.physical_energy_residual_j,
                )
            } else {
                let solid = engine
                    .advance(cx, problem, interfaces, old_h, dt, step_config)
                    .map_err(producer)?;
                let krylov = solid
                    .newton
                    .history
                    .iter()
                    .try_fold(0_usize, |sum, row| sum.checked_add(row.linear_iterations))
                    .ok_or_else(|| budget("enthalpy Krylov work overflow"))?;
                work.add(1, solid.newton.iterations, krylov)?;
                let states = names
                    .iter()
                    .map(|name| {
                        solid
                            .robin_fluxes
                            .iter()
                            .find(|row| row.region == *name)
                            .map(SolidRegionState::from_robin_flux)
                            .ok_or_else(|| bad("enthalpy lacks a network convection region"))
                    })
                    .collect::<Result<Vec<_>>>()?;
                let (residual, energy) = (solid.newton.residual_norm, solid.energy_residual_j);
                (solid, states, None, residual, energy)
            };
            last = Some((solid, heat, residual, energy));
            Ok(states)
        })();
        result.map_err(|error| {
            failure = Some(error);
            AirflowError::Cancelled {
                iteration: 0,
                references_k: references.to_vec(),
            }
        })
    });
    if let Some(error) = failure {
        return Err(error);
    }
    let coupled = coupled.map_err(producer)?;
    let (solid, radiation, physical_residual_j, energy_residual_j) =
        last.ok_or_else(|| bad("enthalpy coupling returned without a solid response"))?;
    let radiation_w = radiation.as_ref().map_or(0.0, |heat| heat.outward_w);
    let defect = finite(
        solid.stored_energy_change_j
            - dt * (solid.source_w - coupled.transport.external_heat_gain_w - radiation_w),
    )?;
    if defect.abs() > step_config.energy_tolerance_j {
        return Err(producer(format!(
            "coupled enthalpy energy residual {defect} J exceeds {} J",
            step_config.energy_tolerance_j
        )));
    }
    poll(cx)?;
    Ok(Endpoint {
        coupled,
        solid,
        radiation,
        physical_residual_j,
        energy_residual_j,
    })
}

fn summary(cx: &Cx<'_>, config: &Config, masses: &[f64], h: &[f64]) -> Result<(String, f64)> {
    let (mut minimum, mut maximum) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut mass, mut liquid, mut energy) = (0.0, 0.0, 0.0);
    for (i, (&h, &m)) in h.iter().zip(masses).enumerate() {
        if i % 512 == 0 {
            poll(cx)?;
        }
        minimum = minimum.min(h);
        maximum = maximum.max(h);
        let phase = config
            .curve
            .state_at_specific_enthalpy(h)
            .map_err(producer)?;
        mass = finite(mass + m)?;
        liquid = finite(liquid + m * phase.liquid_mass_fraction())?;
        energy = finite(energy + m * h)?;
    }
    Ok((
        format!(
            "\"minimum_specific_enthalpy_j_kg\":{},\"maximum_specific_enthalpy_j_kg\":{},\"mean_liquid_mass_fraction\":{}",
            num(minimum)?,
            num(maximum)?,
            num(finite(liquid / mass)?)?
        ),
        energy,
    ))
}

pub(super) fn solve(
    request: &Request,
    cx: &Cx<'_>,
    schedule: &Schedule,
    config: &Config,
) -> Result<String> {
    config.admit(request, schedule)?;
    poll(cx)?;
    let engine = EnthalpyBackwardEuler::uniform(
        cx,
        &request.mesh,
        &config.curve,
        config.density,
        EnthalpyBudget {
            max_vertices: 20_000,
            max_elements: 100_000,
        },
    )
    .map_err(producer)?;
    let masses = engine.reference_nodal_masses_kg();
    let mut h = config.initial_h.clone();
    let (initial, initial_vertex) = initial_objective(request, cx, &schedule.initial)?;
    let (phase, initial_total) = summary(cx, config, masses, &h)?;
    let mut history = vec![format!(
        "{{\"time_s\":0,\"objective_temperature_k\":{},\"active_vertex\":{},\"initial_state\":true,{phase}}}",
        num(initial)?,
        initial_vertex.map_or_else(|| "null".into(), |v| v.to_string())
    )];
    let (mut peak, mut peak_time) = (initial, 0.0);
    let mut first_violation = schedule.limit.filter(|&limit| initial > limit).map(|_| 0.0);
    let (mut time, mut stored, mut input, mut exhaust, mut radiative) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut completed = 0_usize;
    let mut work = Work::default();
    let mut final_result = None;
    let mut final_liquid = Vec::new();
    for (ordinal, interval) in schedule.intervals.iter().enumerate() {
        poll(cx)?;
        interval.workload.validate(request, cx)?;
        let flow = match (&request.fan, interval.speed) {
            (Some(fan), Some(speed)) => fan.solve(cx, &request.graph, request.limits, speed)?,
            (None, None) => request.flow(cx)?,
            _ => return Err(bad("enthalpy drive/speed mismatch")),
        };
        let base = request
            .surfaces
            .iter()
            .map(|s| (s.name.clone(), s.h))
            .collect();
        let (coefficients, mut derived) = convection::resolve(request, cx, &flow, &base)?;
        let network = request.transport(cx, &flow, &coefficients)?;
        let load = interval.workload.prepare(request, cx)?;
        let workload_json = interval.workload.render()?;
        let start = time;
        let end = finite(start + interval.duration)?;
        for step in 1..=interval.steps {
            poll(cx)?;
            if completed >= schedule.max_steps {
                return Err(budget("enthalpy endpoint budget exhausted"));
            }
            let endpoint = if step == interval.steps {
                end
            } else {
                start + interval.duration * (step as f64 / interval.steps as f64)
            };
            let dt = finite(endpoint - time)?;
            if dt <= 0.0 {
                return Err(bad("enthalpy endpoint time is not representable"));
            }
            let solved = advance(
                request,
                cx,
                config,
                &engine,
                &network,
                &coefficients,
                &h,
                &load.source,
                dt,
                &mut work,
            )?;
            for c in &derived {
                c.check_direction(&solved.coupled.solid, request.limits.heat)?;
            }
            if let Some(expected) = load.expected_power_w {
                if (solved.solid.source_w - expected).abs() > request.limits.heat {
                    return Err(producer(
                        "enthalpy source disagrees with component workload",
                    ));
                }
            }
            let radiation_w = solved.radiation.as_ref().map_or(0.0, |heat| heat.outward_w);
            let state =
                request
                    .objective
                    .evaluate(cx, &solved.solid.temperature, &solved.coupled.solid)?;
            if state.value > peak {
                peak = state.value;
                peak_time = endpoint;
            }
            if first_violation.is_none() && schedule.limit.is_some_and(|limit| state.value > limit)
            {
                first_violation = Some(endpoint);
            }
            stored = finite(stored + solved.solid.stored_energy_change_j)?;
            input = finite(input + dt * solved.solid.source_w)?;
            exhaust = finite(exhaust + dt * solved.coupled.transport.external_heat_gain_w)?;
            radiative = finite(radiative + dt * radiation_w)?;
            let (phase, _) = summary(cx, config, masses, &solved.solid.specific_enthalpy_j_kg)?;
            let radiation_field = if solved.radiation.is_some() {
                format!(",\"radiative_heat_w\":{}", num(radiation_w)?)
            } else {
                String::new()
            };
            history.push(format!("{{\"time_s\":{},\"dt_s\":{},\"interval\":{ordinal},{workload_json},\"fan_speed_ratio\":{},\"objective_temperature_k\":{},\"active_vertex\":{},\"source_w\":{},\"air_heat_gain_w\":{},\"stored_energy_change_j\":{},\"solid_energy_residual_j\":{},\"coupled_energy_residual_j\":{},\"physical_residual_norm_j\":{},\"coupling_iterations\":{},\"estimated_local_error_ratio\":null,{phase}{radiation_field}}}",
                num(endpoint)?,num(dt)?,optional(interval.speed)?,num(state.value)?,
                state.vertex.map_or_else(||"null".into(),|v|v.to_string()),num(solved.solid.source_w)?,
                num(solved.coupled.transport.external_heat_gain_w)?,num(solved.solid.stored_energy_change_j)?,
                num(solved.energy_residual_j)?,num(solved.solid.stored_energy_change_j-dt*(solved.solid.source_w-solved.coupled.transport.external_heat_gain_w-radiation_w))?,
                num(solved.physical_residual_j)?,solved.coupled.iterations));
            // Publication follows complete solid, air, radiation and energy acceptance.
            h.clone_from(&solved.solid.specific_enthalpy_j_kg);
            time = endpoint;
            completed += 1;
            if ordinal + 1 == schedule.intervals.len() && step == interval.steps {
                final_liquid = solved.solid.liquid_mass_fraction;
                let evaluated = Evaluation {
                    coupled: solved.coupled,
                    temperatures: solved.solid.temperature,
                    gradient: None,
                    objective: state.value,
                    objective_state: state,
                    robin_total_w: solved.solid.robin_out_w,
                    source_total_w: solved.solid.source_w,
                    htc: network
                        .regions()
                        .iter()
                        .map(|name| coefficients[*name])
                        .collect(),
                    convection: std::mem::take(&mut derived),
                    contact_fluxes: solved.solid.contact_fluxes,
                };
                let mut result = render(request, &flow, &evaluated)?;
                if let (Some(policy), Some(heat)) = (&request.radiation, &solved.radiation) {
                    let prefix = result
                        .strip_suffix("}\n")
                        .ok_or_else(|| bad("internal enthalpy result framing"))?;
                    result = format!(
                        "{prefix},\"radiation\":{}}}\n",
                        policy.endpoint_report(heat)?
                    );
                }
                final_result = Some(match (&request.fan, interval.speed) {
                    (Some(fan), Some(speed)) => fan.attach(result, &flow, speed)?,
                    _ => result,
                });
            }
        }
    }
    if completed != schedule.total_steps {
        return Err(bad(
            "enthalpy completed endpoint count differs from schedule",
        ));
    }
    let defect = finite(stored - input + exhaust + radiative)?;
    if defect.abs() > finite(request.limits.heat * time)? {
        return Err(producer("whole-window enthalpy energy gate failed"));
    }
    let (_, final_total) = summary(cx, config, masses, &h)?;
    let result = final_result.ok_or_else(|| bad("enthalpy trajectory has no final step"))?;
    let prefix = result
        .strip_suffix("}\n")
        .ok_or_else(|| bad("internal enthalpy result framing"))?;
    let radiation_field = if request.radiation.is_some() {
        format!(",\"radiative_energy_loss_j\":{}", num(radiative)?)
    } else {
        String::new()
    };
    let Globalization::LineSearch(line) = config.newton.globalization else {
        return Err(bad("enthalpy requires its admitted line search"));
    };
    let policy = format!(
        "{{\"material_card_identity\":{},\"chart_identity\":{},\"source\":{},\"reference_density_kg_m3\":{},\"initial_total_enthalpy_j\":{},\"final_total_enthalpy_j\":{},\"max_iterations\":{},\"residual_rtol\":{},\"residual_atol_j\":{},\"linear_restart\":{},\"max_linear_cycles\":{},\"armijo_c\":{},\"shrink\":{},\"max_backtracks\":{},\"solid_solves\":{},\"newton_updates\":{},\"krylov_iterations\":{},\"scope\":\"one caller-declared fixed-density equilibrium chart; reference masses and geometry fixed; source/identifier are declarations, not verified material evidence; sensible and latent energy counted once; final nodal h is restart history\"}}",
        quote(&config.curve.material_card_identity().to_hex()),
        quote(&config.curve.identity().to_hex()),
        quote(&config.source),
        num(config.density)?,
        num(initial_total)?,
        num(final_total)?,
        config.max_iterations,
        num(config.newton.relative_tolerance)?,
        num(config.newton.absolute_tolerance)?,
        config.newton.linear_restart,
        config.newton.max_linear_cycles,
        num(line.armijo)?,
        num(line.contraction)?,
        config.max_backtracks,
        work.solves,
        work.updates,
        work.krylov
    );
    let output = format!(
        "{prefix},\"solid_specific_enthalpies_j_kg\":{},\"solid_liquid_mass_fractions\":{},\"transient\":{{\"scheme\":\"backward-euler-total-enthalpy\",\"air_model\":\"quasi-steady endpoint mixing; no fluid storage or travel delay\",\"time_s\":{},\"steps\":{completed},\"total_solid_solves\":{},\"forward_solid_solves\":{},\"sampled_peak_objective_k\":{},\"sampled_peak_time_s\":{},\"temperature_limit_k\":{},\"first_sampled_violation_s\":{},\"stored_energy_change_j\":{},\"input_energy_j\":{},\"air_energy_gain_j\":{},\"energy_residual_j\":{},\"history\":[{}],\"adaptive\":null,\"nonlinear\":null,\"adjoint\":null,\"enthalpy\":{policy}{radiation_field},\"scope\":\"fixed forward workload/fan schedule; accepted enthalpy is physical history; temperatures and mass-weighted phase summaries observe that state; sampled endpoints do not bound inter-step peaks; no moving geometry, melt flow, adaptive/repeated/design/study/adjoint or enclosure mode\"}}}}\n",
        numbers(&h)?,
        numbers(&final_liquid)?,
        num(time)?,
        work.solves,
        work.solves,
        num(peak)?,
        num(peak_time)?,
        optional(schedule.limit)?,
        optional(first_violation)?,
        num(stored)?,
        num(input)?,
        num(exhaust)?,
        num(defect)?,
        history.join(",")
    );
    poll(cx)?;
    Ok(output)
}
