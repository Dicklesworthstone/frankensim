//! Multi-branch resistance design using the existing projected AL engine.
//!
//! Minimize a declared weighted squared log-change from the original losses,
//! subject to a mean or explicitly smoothed peak limit and multiplier bounds.
//! Every sample re-solves hydraulics, flow-derived convection and the coupled
//! solid. The base cost penalizes design changes; an optional term charges
//! nominal fan air power. Neither term is a manufacturing-cost model.
//! A feasible incumbent is independently replayed before publication; terminal
//! optimizer KKT diagnostics are never attributed to a different incumbent.

use super::super::*;
use fs_ascent::projected_al::{
    ProjectedAlError, ProjectedAlOptions, ProjectedAlSample, ProjectedAlState,
    ProjectedAlStop, ProjectedAlWork,
};
use fs_math::det;
use std::ops::ControlFlow;

mod fan_control;
mod power;
mod operating;
use fan_control::FanControl;

const SCHEMA: &str = "frankensim.cooling-resistance-design.v1";
const HELP: &str = "Usage: frankensim [--json] cooling-network resistance-design <model.json> <design.json>\n\nRequires --features cooling-design. Optimize bounded branch-resistance\nmultipliers and optional fan_speed_control against a mean or hotspot limit,\nusing complete hydraulic/convection/solid adjoints. Optional air_power_cost\nadds a declared weight on nominal fan air power. The model must set objective.gradient=true.\nMaximum selectors require explicit maximum_smoothing_k in the design policy;\nthe smooth upper envelope and actual selected peak must both pass on replay.\nOptional operating_envelope checks every declared source/inlet case, including nominal.\nSee examples/cooling-network/RESISTANCE_DESIGN.md and OPERATING_ENVELOPE.md.\n";

#[derive(Debug)]
struct Control {
    branch: usize,
    name: String,
    baseline: f64,
    minimum: f64,
    maximum: f64,
    weight: f64,
}

#[derive(Debug)]
struct Specification {
    controls: Vec<Control>,
    fan_control: Option<FanControl>,
    air_power_weight: f64,
    operating_envelope: Option<operating::Cases>,
    limit: f64,
    scale: f64,
    maximum_smoothing_k: Option<f64>,
    evaluations: usize,
    iterations: usize,
    options: ProjectedAlOptions,
}

impl Specification {
    fn parse(text: &str, request: &Request) -> Result<Self> {
        if !request.gradient {
            return Err(bad("resistance design requires objective.gradient=true"));
        }
        if request.design.is_some() || request.fan_speed_design.is_some()
            || request.transient.is_some() || request.mesh_convergence.is_some()
            || request.radiation.is_some() || request.recirculation.is_some()
        {
            return Err(bad("resistance design admits only a standalone steady once-through model; no nested search, radiation, return or transient model"));
        }
        if text.len() as u64 > MAX_INPUT_BYTES { return Err(bad("design request exceeds 16 MiB")); }
        let root = J::parse(text).map_err(bad_parse)?;
        object(&root, &["schema", "temperature_limit_k", "temperature_scale_k", "controls",
            "max_evaluations", "max_iterations", "max_multiplier_updates", "max_backtracks",
            "kkt_tolerance", "maximum_smoothing_k", "fan_speed_control", "air_power_cost", "operating_envelope"], "resistance design")?;
        if get(&root, "schema")?.as_str() != Some(SCHEMA) { return Err(bad("unsupported resistance-design schema")); }
        let maximum_smoothing_k = root.get("maximum_smoothing_k")
            .map(|value| positive(value, "maximum_smoothing_k")).transpose()?;
        match (request.objective.is_mean(), maximum_smoothing_k) {
            (true, Some(_)) => return Err(bad("maximum_smoothing_k cannot change a mean-wall requirement into a peak requirement")),
            (false, None) => return Err(bad("maximum-temperature resistance design requires explicit maximum_smoothing_k")),
            _ => {},
        }
        let mut controls = Vec::new();
        let mut selected = BTreeSet::new();
        for row in array(get(&root, "controls")?, "controls", 64)? {
            object(row, &["branch", "minimum_multiplier", "maximum_multiplier", "weight"], "control")?;
            let name = string(get(row, "branch")?, "control.branch")?;
            let branch = request.graph.branches().iter().position(|b| b.loss.name == name)
                .ok_or_else(|| bad(format!("unknown resistance-design branch {name}")))?;
            if !selected.insert(branch) { return Err(bad("resistance-design branch repeated")); }
            let baseline = request.graph.branches()[branch].loss.resistance.value();
            let minimum = positive(get(row, "minimum_multiplier")?, "minimum_multiplier")?;
            let maximum = positive(get(row, "maximum_multiplier")?, "maximum_multiplier")?;
            if minimum >= maximum || minimum > 1.0 || maximum < 1.0
                || baseline * minimum <= 0.0 || !(baseline * maximum).is_finite()
            { return Err(bad("finite positive resistance bounds must contain the original multiplier 1 and have nonzero width")); }
            controls.push(Control { branch, name, baseline, minimum, maximum,
                weight: positive(get(row, "weight")?, "control.weight")? });
        }
        let fan_control = root.get("fan_speed_control")
            .map(|value| FanControl::parse(value, request)).transpose()?;
        let air_power_weight = power::parse(root.get("air_power_cost"), request)?;
        let evaluations = count(get(&root, "max_evaluations")?, "max_evaluations", 4096)?;
        let limit = positive(get(&root, "temperature_limit_k")?, "temperature_limit_k")?;
        let operating_envelope = root.get("operating_envelope")
            .map(|value| operating::Cases::parse(value, request, limit, maximum_smoothing_k)).transpose()?;
        let cases = operating_envelope.as_ref().map_or(1, operating::Cases::len);
        if evaluations < 2 * cases {
            return Err(bad("max_evaluations must reserve every initialization case and every final replay case"));
        }
        let tolerance = positive(get(&root, "kkt_tolerance")?, "kkt_tolerance")?;
        let options = ProjectedAlOptions {
            max_dimension: 65,
            // The optimizer counts full case sweeps; the public budget counts
            // physical case evaluations, including a reserved final sweep.
            max_evaluations: evaluations / cases - 1,
            max_multiplier_updates: count(get(&root, "max_multiplier_updates")?, "max_multiplier_updates", 4096)?,
            max_backtracks: count(get(&root, "max_backtracks")?, "max_backtracks", 64)?,
            tolerance,
            // Versioned algorithm policy, not additional user work allowances.
            initial_penalty: 10.0, max_penalty: 1e10,
            inner_tolerance: tolerance.max(1e-2),
            min_spectral_step: 1e-8, max_spectral_step: 1e8,
        };
        options.validate::<Failure>(controls.len() + usize::from(fan_control.is_some())).map_err(optimizer_error)?;
        Ok(Self { controls, fan_control, air_power_weight, operating_envelope, maximum_smoothing_k,
            limit,
            scale: positive(get(&root, "temperature_scale_k")?, "temperature_scale_k")?,
            evaluations,
            iterations: count(get(&root, "max_iterations")?, "max_iterations", 4096)?, options })
    }

    fn case_count(&self) -> usize { self.operating_envelope.as_ref().map_or(1, operating::Cases::len) }

    fn dimension(&self) -> usize { self.controls.len() + usize::from(self.fan_control.is_some()) }

    fn bounds(&self) -> (Vec<f64>, Vec<f64>) {
        let mut lower: Vec<_> = self.controls.iter().map(|c| det::ln(c.minimum)).collect();
        let mut upper: Vec<_> = self.controls.iter().map(|c| det::ln(c.maximum)).collect();
        if let Some(fan) = &self.fan_control {
            let (lo, hi) = fan.bounds();
            lower.push(lo);
            upper.push(hi);
        }
        (lower, upper)
    }

    fn install(&self, request: &mut Request, original: &LossGraph, point: &[f64]) -> Result<()> {
        if point.len() != self.dimension() { return Err(bad("resistance design coordinate shape changed")); }
        let speed = self.fan_control.as_ref()
            .map(|control| control.candidate(request, point[self.controls.len()])).transpose()?;
        let mut branches = original.branches().to_vec();
        for (control, &coordinate) in self.controls.iter().zip(point) {
            if !coordinate.is_finite() || coordinate < det::ln(control.minimum)
                || coordinate > det::ln(control.maximum)
            { return Err(bad("resistance design point is outside its declared box")); }
            // Clamp only exp/log endpoint rounding; coordinates themselves must
            // already lie in the box. Never project an invalid trial silently.
            let multiplier = det::exp(coordinate).clamp(control.minimum, control.maximum);
            branches[control.branch].loss.resistance = LossResistance::new(control.baseline * multiplier);
        }
        let graph = LossGraph::new(original.node_count(), branches).map_err(producer)?;
        // Operating cases install their own admitted objective and smoothing.
        if self.operating_envelope.is_none() {
            if let Some(tau) = self.maximum_smoothing_k { request.objective.enable_upper_envelope(tau)?; }
        }
        if let Some(speed) = speed {
            request.fan.as_mut().ok_or_else(|| bad("joint design lost its fan drive"))?.speed_ratio = speed;
        }
        request.graph = graph;
        Ok(())
    }

    fn passes(&self, physical: &Evaluation) -> bool {
        // The optimizer uses the smooth upper envelope, but publication also
        // checks the exact maximum on the caller's original vertex selection.
        physical.objective <= self.limit && physical.objective_state.physical_value() <= self.limit
    }
}

fn finite(value: f64) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(producer("nonfinite resistance-design arithmetic")) }
}

fn optimizer_error(error: ProjectedAlError<Failure>) -> Failure {
    match error {
        ProjectedAlError::Evaluation(error) => error,
        ProjectedAlError::Cancelled => Failure { code: "cooling-network-cancelled", message: "resistance design cancelled before publication".into() },
        ProjectedAlError::Invalid(message) => bad(message),
    }
}

fn checkpoint(cx: &Cx<'_>, _: ProjectedAlWork) -> ControlFlow<()> {
    if cx.checkpoint().is_ok() { ControlFlow::Continue(()) } else { ControlFlow::Break(()) }
}

fn evaluate(
    request: &mut Request, original: &LossGraph, design: &Specification, point: &[f64], cx: &Cx<'_>,
) -> Result<(ProjectedAlSample, GraphSolution, Evaluation)> {
    poll(cx)?;
    design.install(request, original, point)?;
    let flow = request.flow(cx)?;
    let coefficients = request.surfaces.iter().map(|s| (s.name.clone(), s.h)).collect();
    let value = request.evaluate(cx, &flow, &coefficients, true)?;
    let derivative = value.gradient.as_ref().ok_or_else(|| bad("resistance design requires the complete hydraulic adjoint"))?;
    let mut objective = 0.0;
    let mut gradient = Vec::with_capacity(point.len());
    let mut constraint_gradient = Vec::with_capacity(point.len());
    for (control, &x) in design.controls.iter().zip(point) {
        poll(cx)?;
        let partial = derivative.log_resistance(control.branch)
            .ok_or_else(|| bad("hydraulic derivative unavailable at a resistance-design trial; no frozen-flow fallback"))?;
        objective = finite(objective + 0.5 * control.weight * x * x)?;
        gradient.push(finite(control.weight * x)?);
        constraint_gradient.push(finite(partial / design.scale)?);
    }
    if let Some(control) = &design.fan_control {
        let x = point[design.controls.len()];
        let partial = derivative.log_speed()
            .ok_or_else(|| bad("fan-speed derivative unavailable at a joint-design trial; no frozen-flow fallback"))?;
        objective = finite(objective + 0.5 * control.weight * x * x)?;
        gradient.push(finite(control.weight * x)?);
        constraint_gradient.push(finite(partial / design.scale)?);
    }
    let constraint = finite((value.objective - design.limit) / design.scale)?;
    let mut sample = ProjectedAlSample { objective, gradient, constraint, constraint_gradient };
    power::add_cost(request, &flow, design, cx, &mut sample)?;
    poll(cx)?;
    Ok((sample, flow, value))
}

#[derive(Debug)]
struct Incumbent {
    point: Vec<f64>, sample: ProjectedAlSample, cases: Vec<operating::Signature>,
}

/// Only this owned request is mutated. Caller files and the source graph are
/// unchanged; every new trial starts from the original, not the preceding trial.
fn solve(mut request: Request, design: Specification, cx: &Cx<'_>) -> Result<String> {
    poll(cx)?;
    let original = request.graph.clone();
    let (lower, upper) = design.bounds();
    let mut incumbent: Option<Incumbent> = None;
    let mut initial_temperature = None;
    let mut evaluations = 0_usize;
    let mut solid_solves = 0_usize;
    let (state, report) = {
        let mut sample = |point: &[f64]| -> Result<Option<ProjectedAlSample>> {
            evaluations += design.case_count();
            let candidate = operating::evaluate_candidate(&mut request, &original, &design, point, cx, false)?;
            initial_temperature.get_or_insert(candidate.actual_temperature);
            solid_solves = solid_solves.checked_add(candidate.solid_solves)
                .ok_or_else(|| bad("resistance-design solid work overflow"))?;
            // Rejected optimizer steps are still real solved candidates. Keep
            // the cheapest actually feasible sample, not a predicted crossing.
            if candidate.passes
                && incumbent.as_ref().is_none_or(|old| candidate.sample.objective < old.sample.objective)
            { incumbent = Some(Incumbent { point: point.to_vec(), sample: candidate.sample.clone(),
                cases: candidate.signatures }); }
            Ok(Some(candidate.sample))
        };
        let mut state = ProjectedAlState::try_new(&vec![0.0; design.dimension()], &lower, &upper,
            design.options, &mut sample, |work| checkpoint(cx, work)).map_err(optimizer_error)?;
        let report = state.try_run(design.iterations, &mut sample, |work| checkpoint(cx, work))
            .map_err(optimizer_error)?;
        (state, report)
    };
    let best = incumbent.ok_or_else(|| Failure { code: "cooling-network-design-no-feasible-sample",
        message: format!("resistance search stopped at {:?} after {evaluations} evaluations without an evaluated passing design; not a proof of infeasibility", report.stop) })?;
    let replay_evaluations = design.case_count();
    if evaluations + replay_evaluations > design.evaluations {
        return Err(bad("resistance-design replay budget was not reserved"));
    }
    evaluations += replay_evaluations;
    let selected = operating::evaluate_candidate(&mut request, &original, &design, &best.point, cx, true)?;
    solid_solves = solid_solves.checked_add(selected.solid_solves)
        .ok_or_else(|| bad("resistance-design solid work overflow"))?;
    if !selected.passes || selected.sample != best.sample || selected.signatures != best.cases {
        return Err(producer("selected resistance design did not reproduce every case's feasible value and gradient"));
    }
    let replayed = &selected.sample;
    let flow = &selected.flow;
    let physical = &selected.physical;
    let operating_envelope = selected.report.as_deref().unwrap_or("null");
    let converged = report.stop == ProjectedAlStop::Converged && state.point() == best.point.as_slice();
    let mut output = render(&request, flow, physical)?;
    if let Some(fan) = &request.fan { output = fan.attach(output, flow, fan.speed_ratio)?; }
    let prefix = output.strip_suffix("}\n").ok_or_else(|| bad("resistance-design result framing"))?;
    let controls = design.controls.iter().zip(&best.point).map(|(c, &x)| -> Result<String> {
        let resistance = request.graph.branches()[c.branch].loss.resistance.value();
        Ok(format!("{{\"branch\":{},\"baseline_resistance_pa_s2_m6\":{},\"selected_resistance_pa_s2_m6\":{},\"multiplier\":{},\"log_multiplier\":{},\"minimum_multiplier\":{},\"maximum_multiplier\":{},\"weight\":{}}}",
            quote(&c.name), num(c.baseline)?, num(resistance)?, num(resistance / c.baseline)?, num(x)?,
            num(c.minimum)?, num(c.maximum)?, num(c.weight)?))
    }).collect::<Result<Vec<_>>>()?.join(",");
    let fan_control = design.fan_control.as_ref()
        .map(|control| control.render(&request, best.point[design.controls.len()], design.controls.len()))
        .transpose()?.unwrap_or_else(|| "null".into());
    let air_power_cost = power::render(&request, flow, design.air_power_weight)?;
    let constraint_detail = format!("\"constraint_form\":{},\"selected_temperature_k\":{},\"selected_constraint_temperature_k\":{},\"maximum_smoothing_k\":{}",
        quote(if design.operating_envelope.is_some() { "finite-operating-case-upper-envelope" }
            else if design.maximum_smoothing_k.is_some() { "smooth-maximum-upper-envelope" } else { "area-mean" }),
        num(selected.actual_temperature)?, num(selected.constraint_temperature)?, optional(design.maximum_smoothing_k)?);
    let result = format!("{prefix},\"resistance_design\":{{\"schema\":{},\"status\":{},\"selected_cost\":{},\"temperature_limit_k\":{},\"temperature_scale_k\":{},\"initial_temperature_k\":{},\"controls\":[{}],\"fan_speed_control\":{fan_control},\"air_power_cost\":{air_power_cost},{constraint_detail},\"complete_evaluations\":{},\"replay_evaluations\":{replay_evaluations},\"operating_envelope\":{operating_envelope},\"optimizer_evaluation_unit\":\"complete-case-sweep\",\"total_solid_solves\":{},\"optimizer_terminal\":{{\"stop\":{},\"point_log_multipliers\":{},\"cost\":{},\"scaled_constraint\":{},\"stationarity\":{},\"feasibility\":{},\"complementarity\":{},\"iterations\":{},\"evaluations\":{},\"multiplier_updates\":{},\"rejected_trials\":{}}},\"policy\":\"projected-al-resistance-v1; weighted squared log-change plus optional declared fan air-power cost; initial penalty 10, maximum penalty 1e10, initial inner tolerance max(kkt_tolerance,0.01), spectral steps [1e-8,1e8]\",\"search_claim\":\"independently replayed nominally feasible selected temperature design over every declared case; peak designs require both the explicit upper envelope and actual selected discrete maximum to pass; feasible-incumbent is not optimizer convergence; terminal KKT diagnostics belong only to the terminal point; no global optimum, unsampled-operating-condition, continuum-temperature, hardware, geometry-realizability or uncertainty guarantee\"}}}}\n",
        quote(SCHEMA), quote(if converged { "locally-converged" } else { "feasible-incumbent" }),
        num(replayed.objective)?, num(design.limit)?, num(design.scale)?, optional(initial_temperature)?, controls,
        evaluations, solid_solves, quote(&format!("{:?}", report.stop)), numbers(state.point())?,
        num(report.objective)?, num(report.constraint)?, num(report.kkt.stationarity)?,
        num(report.kkt.feasibility)?, num(report.kkt.complementarity)?, report.work.iterations,
        report.work.evaluations, report.work.multiplier_updates, report.work.rejected_trials);
    poll(cx)?;
    Ok(result)
}

fn read(path: &OsString) -> Result<String> {
    let mut text = String::new();
    File::open(std::path::Path::new(path)).map_err(|e| bad(e.to_string()))?.take(MAX_INPUT_BYTES + 1)
        .read_to_string(&mut text).map_err(|e| bad(e.to_string()))?;
    if text.len() as u64 > MAX_INPUT_BYTES { return Err(bad("resistance-design input exceeds 16 MiB")); }
    Ok(text)
}

pub(in crate::network_command) fn run(args: &[OsString], json_mode: bool) -> CommandOutput {
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        return CommandOutput { exit_code: exit::SUCCESS, stdout: if json_mode {
            format!("{{\"schema\":{},\"help\":{}}}\n", quote(SCHEMA), quote(HELP))
        } else { HELP.into() }, stderr: String::new() };
    }
    if args.len() != 2 { return diagnostic(exit::USAGE, bad(HELP), json_mode); }
    let parsed = (|| {
        let request = Request::parse(&read(&args[0])?)?;
        let design = Specification::parse(&read(&args[1])?, &request)?;
        Ok::<_, Failure>((request, design))
    })();
    let (request, design) = match parsed { Ok(pair) => pair,
        Err(error) => return diagnostic(exit::REFUSED, error, json_mode) };
    let duration = Duration::from_secs_f64(request.limits.wall_seconds);
    let gate = CancelGate::new();
    let start = Instant::now();
    let result = std::thread::scope(|scope| {
        let (stop, stopped) = std::sync::mpsc::channel::<()>();
        let gate_ref = &gate;
        scope.spawn(move || {
            if matches!(stopped.recv_timeout(duration), Err(std::sync::mpsc::RecvTimeoutError::Timeout)) { gate_ref.request(); }
        });
        let result = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cx = Cx::new(&gate, arena, StreamKey { seed: request.seed, kernel_id: 718, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic);
            solve(request, design, &cx)
        });
        let _ = stop.send(());
        result
    });
    if start.elapsed() >= duration || gate.is_requested() {
        return diagnostic(exit::BUDGET, Failure { code: "cooling-network-time-budget",
            message: "resistance-design wall-time budget exhausted before publication".into() }, json_mode);
    }
    match result {
        Ok(stdout) => CommandOutput { exit_code: exit::SUCCESS, stdout, stderr: String::new() },
        Err(error) => diagnostic(exit::REFUSED, error, json_mode),
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_peak;

#[cfg(test)]
mod tests_joint;
