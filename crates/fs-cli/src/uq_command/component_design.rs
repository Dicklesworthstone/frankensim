//! Independent component/interval allocation through the existing cooling CLI.
//! The orchestration shares UQ's bounded same-executable runner, not its
//! probability model. Every feasibility decision uses complete trajectories.
use super::*;
use std::time::{Duration, Instant};
mod input;
mod search;
mod thermal_constraints;
mod constraint_execution;
use input::Plan;

const DESIGN_SCHEMA: &str = "frankensim.cooling-component-design.result.v1";
const DESIGN_HELP: &str = "Usage: frankensim [--json] cooling-component-design <base-request.json> <allocation.json>\n\nAllocate absolute component watts in declared priority order. Each control names\none base interval and affects every occurrence in its fixed repeated schedule.\nEvery candidate runs the real cooling-network command from the original initial\nfield. A requested sampled-peak component adjoint may guide safeguarded trials;\nwithout one, bisection performs no hidden reverse solves. Bounds, thermal slack\nand watt-bracket tolerances are explicit. Optional thermal_constraints add named\ncomponent-footprint or existing temperature objectives, each with its OWN limit.\nEvery limit must pass. The worst signed temperature excess guides the bracket,\nnot the largest absolute temperature. Each objective reruns the same full model\nand its own requested adjoint; every run consumes the cumulative step budget.\nEvaluation, cumulative trajectory-step and wall budgets cover ALL priorities.\nBudget exhaustion after a passing baseline returns that accepted allocation with\nstatus=budget-exhausted and exit 6, never a completed search. Model failures refuse.\nThe result retains cooling_result and resolved_request for ordinary exact replay.\nThis is not a global optimum or continuous-time temperature certificate.\nSee examples/cooling-network/COMPONENT_POWER_ALLOCATION.md and MULTI_LIMIT_ALLOCATION.md.\n";

pub(in crate) fn run(args: &[OsString], json_mode: bool) -> CommandOutput {
    if args.len() == 1 && (args[0] == "--help" || args[0] == "-h") {
        return CommandOutput { exit_code: exit::SUCCESS, stdout: if json_mode {
            format!("{{\"schema\":{},\"help\":{}}}\n", quote(DESIGN_SCHEMA), quote(DESIGN_HELP))
        } else { DESIGN_HELP.into() }, stderr: String::new() };
    }
    if args.len() != 2 { return failure(exit::USAGE, bad(DESIGN_HELP), json_mode); }
    if !cfg!(unix) { return failure(exit::REFUSED, bad("component allocation currently requires a Unix cooling child"), json_mode); }
    let result = (|| -> Result<(u8,String)> {
        let base = J::parse(&read(&args[0], MAX_BASE_BYTES)?).map_err(|e| bad(e.to_string()))?;
        let spec = J::parse(&read(&args[1], MAX_UQ_BYTES)?).map_err(|e| bad(e.to_string()))?;
        let (plan, constraints) = parse_plan(&base, &spec)?;
        let deadline = Instant::now().checked_add(Duration::from_secs_f64(plan.wall_seconds))
            .ok_or_else(|| bad("allocation wall deadline is not representable"))?;
        let criteria = constraints.len() + 1;
        let mut result = search::allocate_with_work(&plan, deadline, criteria, |values, work|
            constraint_execution::evaluate(&plan, &constraints, values, deadline, work))?;
        let resolved = plan.request(&result.values)?;
        let selected = plan.axes.iter().zip(&result.values).map(|(axis,&value)| Ok(J::Object(vec![
            ("component".into(),J::Str(axis.component.clone())), ("interval".into(),usize_value(axis.interval)),
            ("min_power_w".into(),number_value(axis.minimum)?), ("max_power_w".into(),number_value(axis.maximum)?),
            ("selected_power_w".into(),number_value(value)?),
        ]))).collect::<Result<Vec<_>>>()?;
        let thermal_constraints = result.passing.constraint_report()?;
        if result.reason.is_none() && Instant::now() >= deadline {
            result.reason = Some("allocation deadline exhausted during finalization".into());
        }
        let complete = result.reason.is_none();
        let document = J::Object(vec![
            ("schema".into(),J::Str(DESIGN_SCHEMA.into())), ("authority".into(),J::Str("nominal-estimate".into())),
            ("status".into(),J::Str(if complete { "priority-allocation-complete" } else { "budget-exhausted" }.into())),
            ("reason".into(),result.reason.map(J::Str).unwrap_or(J::Null)),
            ("temperature_limit_k".into(),number_value(plan.limit)?),
            ("power_tolerance_w".into(),number_value(plan.power_tolerance)?),
            ("temperature_tolerance_k".into(),number_value(plan.temperature_tolerance)?),
            ("baseline_sampled_peak_k".into(),number_value(result.baseline_peak)?),
            ("selected_sampled_peak_k".into(),number_value(result.passing.peak)?),
            ("thermal_constraints".into(),thermal_constraints),
            ("trajectories_per_candidate".into(),usize_value(criteria)),
            ("completed_priorities".into(),usize_value(result.completed)),
            ("selected".into(),J::Array(selected)), ("priority_decisions".into(),J::Array(result.decisions)),
            ("evaluations_attempted".into(),usize_value(result.attempted)),
            ("evaluations_completed".into(),usize_value(result.history.len())),
            ("trajectory_evaluations_attempted".into(),usize_value(result.work.attempted)),
            ("trajectory_evaluations_completed".into(),usize_value(result.work.completed)),
            ("total_completed_trajectory_steps".into(),usize_value(result.steps)),
            ("total_completed_trajectory_solid_solves".into(),usize_value(result.solves)),
            ("newton_trials".into(),usize_value(result.newton_trials)),
            ("search_method".into(),J::Str(if plan.adjoint { "priority-component-adjoint-brackets" } else { "priority-component-bisection" }.into())),
            ("history".into(),J::Array(result.history)), ("resolved_request".into(),resolved),
            ("cooling_result".into(),result.passing.document),
            ("scope".into(),J::Str("declared-order local coordinate allocation; every selected vector actually passes every declared initial/accepted-endpoint temperature limit on the fixed numerical model; selected_sampled_peak_k and cooling_result retain the primary objective, not an artificial combined temperature; each extra criterion reruns the same physical experiment and consumes its own requested adjoint; completed work includes finished criteria of an interrupted candidate set, but that incomplete set never replaces a passing allocation; later loads can change earlier conditional brackets; no global or lexicographic optimality, monotonicity, continuous-time peak, mesh/time error or physical-validation certificate; no work count is invented for a killed child; resolved_request is replayable but not a durable optimizer checkpoint".into())),
        ]);
        publish(document, complete, deadline)
    })();
    match result {
        Ok((exit_code,stdout)) => CommandOutput { exit_code, stdout, stderr: String::new() },
        Err(error) => failure(if error.code == "cooling-network-uq-budget" { exit::BUDGET } else { exit::REFUSED },error,json_mode),
    }
}

fn parse_plan(base: &J, spec: &J) -> Result<(Plan, Vec<thermal_constraints::Constraint>)> {
    let constraints = thermal_constraints::parse(base, spec.get("thermal_constraints"))?;
    let mut core = spec.clone();
    let J::Object(fields) = &mut core else { return Err(bad("component design must be an object")); };
    fields.retain(|(key, _)| key != "thermal_constraints");
    // Reuse the existing source/control/fixed-horizon admission unchanged.
    Ok((Plan::parse(base, &core)?, constraints))
}

fn publish(mut document: J, complete: bool, deadline: Instant) -> Result<(u8,String)> {
    let mut stdout = serialize(&document)?;
    if complete && Instant::now() >= deadline {
        input::put(&mut document,"status",J::Str("budget-exhausted".into()))?;
        input::put(&mut document,"reason",J::Str("allocation deadline exhausted during result serialization".into()))?;
        stdout = serialize(&document)?;
        return Ok((exit::BUDGET,stdout));
    }
    Ok((if complete {exit::SUCCESS} else {exit::BUDGET},stdout))
}

fn failure(exit_code: u8, error: Failure, json_mode: bool) -> CommandOutput {
    CommandOutput { exit_code, stdout: String::new(), stderr: if json_mode {
        format!("{{\"schema\":\"frankensim.cooling-component-design.diagnostic.v1\",\"cause\":{},\"message\":{}}}\n",quote(error.code),quote(&error.message))
    } else { format!("{error}\n") } }
}

pub(super) struct Evaluation {
    document: J,
    peak: f64,
    peak_time: f64,
    steps: usize,
    solves: usize,
    slopes: Vec<Option<f64>>,
    constraints: Option<thermal_constraints::Envelope>,
}
impl Evaluation {
    fn margin(&self, plan: &Plan) -> Result<f64> {
        self.constraints.as_ref().map_or_else(|| checked(self.peak-plan.limit), |e| Ok(e.margin))
    }
    fn active_constraint(&self) -> &str {
        self.constraints.as_ref().map_or("primary", |e| e.active_name())
    }
    fn constraint_report(&self) -> Result<J> {
        self.constraints.as_ref().map(|e| e.report()).transpose().map(|e| e.unwrap_or(J::Null))
    }
}

fn inspect(plan: &Plan, values: &[f64], document: J) -> Result<Evaluation> {
    if values.len() != plan.axes.len() { return Err(model_failure("component observation vector length mismatch")); }
    let peak = plan.qoi.extract(&document)?;
    let repeated = matches!(plan.qoi, model::Qoi::TransientPeak { repeated_cycles: Some(_) });
    let phase = field(&document, if repeated { "repeated_cycles" } else { "transient" })?;
    let steps = integer_raw(field(phase, if repeated { "total_accepted_steps" } else { "steps" })?, "returned trajectory steps")?;
    if steps != plan.planned_steps { return Err(model_failure("allocation candidate did not return the declared complete fixed grid")); }
    let solves = integer_raw(field(phase,"total_solid_solves")?,"returned solid solves")?;
    let peak_time = number(field(phase,"sampled_peak_time_s")?,"returned peak time")?;
    let mut slopes = vec![None;plan.axes.len()];
    if plan.adjoint {
        let adjoint = field(phase,"adjoint")?;
        if adjoint.str_field("method") != Some("discrete-backward-euler-coupled-adjoint")
            || adjoint.str_field("qoi") != Some("sampled-peak")
            || number(field(adjoint,"value_k")?,"adjoint peak")?.to_bits() != peak.to_bits()
            || number(field(adjoint,"time_s")?,"adjoint peak time")?.to_bits() != peak_time.to_bits()
            || integer_raw(field(adjoint,"cycles")?,"adjoint cycles")? != plan.cycles {
            return Err(model_failure("component adjoint does not describe this candidate's all-cycle sampled peak"));
        }
        let intervals = array(field(field(adjoint,"component_power_sensitivities")?,"intervals")?,"component derivative intervals",plan.intervals)?;
        if intervals.len() != plan.intervals { return Err(model_failure("incomplete component derivative interval set")); }
        for (index,axis) in plan.axes.iter().enumerate() {
            let matches: Vec<_> = intervals.iter().filter(|row| row.get("interval")
                .and_then(J::number_raw).and_then(|v|v.parse::<usize>().ok()) == Some(axis.interval)).collect();
            if matches.len()!=1 { return Err(model_failure("ambiguous component derivative interval")); }
            let rows = array(field(matches[0],"rows")?,"component derivative rows",4096)?;
            let matches: Vec<_> = rows.iter().filter(|row| row.str_field("component")==Some(axis.component.as_str())).collect();
            if matches.len()!=1 { return Err(model_failure("missing or ambiguous component derivative")); }
            if number(field(matches[0],"applied_power_w")?,"adjoint component watts")?.to_bits() != values[index].to_bits() {
                return Err(model_failure("component derivative is bound to different candidate watts"));
            }
            slopes[index] = Some(number(field(matches[0],"dtemperature_dpower_w_k_per_w")?,"absolute component derivative")?);
        }
    }
    Ok(Evaluation { document, peak, peak_time, steps, solves, slopes, constraints: None })
}

fn checked(value: f64) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(model_failure("nonfinite component allocation arithmetic")) }
}
fn number_value(value: f64) -> Result<J> { Ok(J::Number { value: checked(value)?, raw: value.to_string() }) }
fn usize_value(value: usize) -> J { J::Number { value: value as f64, raw: value.to_string() } }
fn number_array(values: &[f64]) -> Result<J> { values.iter().map(|&v| number_value(v)).collect::<Result<Vec<_>>>().map(J::Array) }
fn serialize(value: &J) -> Result<String> {
    fn write(value: &J, out: &mut String) -> Result<()> {
        match value {
            J::Null => out.push_str("null"), J::Bool(v) => out.push_str(if *v {"true"} else {"false"}),
            J::Number {value,raw} => { checked(*value)?; out.push_str(raw); }, J::Str(s) => out.push_str(&quote(s)),
            J::Array(rows) => { out.push('['); for (i,v) in rows.iter().enumerate() {if i>0 {out.push(',');} write(v,out)?;} out.push(']'); },
            J::Object(rows) => { out.push('{'); for (i,(k,v)) in rows.iter().enumerate() {if i>0 {out.push(',');} out.push_str(&quote(k));out.push(':');write(v,out)?;}out.push('}'); },
        }
        Ok(())
    }
    let mut out=String::new(); write(value,&mut out)?;out.push('\n');Ok(out)
}

#[cfg(test)]
mod tests;
