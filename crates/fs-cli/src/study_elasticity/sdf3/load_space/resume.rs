//! Restore the existing load-space AL state and re-solve physical endpoints.
use super::*;
use fs_ascent::projected_al::{ProjectedAlCheckpoint, ProjectedAlSample, ProjectedAlState};
use fs_topopt::{GradientDirection, GradientProbe, MultiLoadGradientCheck};

pub(super) const KIND: &str = "study-sdf3-load-space-checkpoint";
pub(super) const MODE: &str = "optimizer-state-restoration-v1";
const SCHEMA: &str = "sdf3-load-space-state.v1";

pub(super) struct Recovery {
    pub(super) checkpoint: LoadSpaceDesignCheckpoint3,
    pub(super) spent: Spent,
}

fn invalid(message: impl Into<String>) -> Failure {
    fail("cli-study-sdf3-load-space-checkpoint", message)
}
fn field<'a>(value: &'a JsonValue, key: &str) -> Result<&'a JsonValue> {
    value.get(key).ok_or_else(|| invalid(format!("missing checkpoint {key}")))
}
fn real(value: &JsonValue, key: &str) -> Result<f64> {
    field(value, key)?.as_f64().filter(|x| x.is_finite())
        .ok_or_else(|| invalid(format!("nonfinite checkpoint {key}")))
}
fn vector(value: &JsonValue, n: usize) -> Result<Vec<f64>> {
    let values = value.as_array().ok_or_else(|| invalid("expected a numerical vector"))?;
    if values.len() != n { return Err(invalid("checkpoint vector has the wrong cell count")); }
    values.iter().map(|v| v.as_f64().filter(|v| v.is_finite())
        .ok_or_else(|| invalid("nonfinite checkpoint vector"))).collect()
}
fn document(bytes: &[u8]) -> Result<JsonValue> {
    JsonValue::parse(std::str::from_utf8(bytes).map_err(|e| invalid(e.to_string()))?)
        .map_err(|e| invalid(e.to_string()))
}
fn work(value: &JsonValue) -> Result<ProjectedAlWork> {
    Ok(ProjectedAlWork { iterations: integer(value, "iterations")?,
        evaluations: integer(value, "evaluations")?,
        multiplier_updates: integer(value, "multiplier_updates")?,
        rejected_trials: integer(value, "rejected_trials")? })
}

pub(super) fn encode(spec: &Spec, state: &State, producer: ContentHash,
    design: &str, iterations: &str) -> String {
    let c = &state.checkpoint;
    let a = &c.optimizer;
    let previous = a.previous_violation.map_or_else(|| "null".into(), |v| format!("{v:.17e}"));
    let best = c.best_feasible_density.as_ref().map_or_else(|| "null".into(), |v| format!("{v:?}"));
    format!(
        "{{\"schema\":{SCHEMA:?},\"study_id\":{:?},\"producer\":{:?},\"design\":{:?},\"iterations\":{:?},\"restoration_evaluations\":{},\"point\":{:?},\"best_feasible_density\":{best},\"sample\":{{\"objective\":{:.17e},\"constraint\":{:.17e},\"gradient\":{:?},\"constraint_gradient\":{:?}}},\"work\":{{\"iterations\":{},\"evaluations\":{},\"multiplier_updates\":{},\"rejected_trials\":{}}},\"multiplier\":{:.17e},\"penalty\":{:.17e},\"spectral_step\":{:.17e},\"inner_tolerance\":{:.17e},\"previous_violation\":{previous}}}",
        spec.id.to_hex(), producer.to_hex(), hash_bytes(design.as_bytes()).to_hex(),
        hash_bytes(iterations.as_bytes()).to_hex(), c.restoration_evaluations, a.point,
        a.sample.objective, a.sample.constraint, a.sample.gradient, a.sample.constraint_gradient,
        a.work.iterations, a.work.evaluations, a.work.multiplier_updates, a.work.rejected_trials,
        a.multiplier, a.penalty, a.spectral_step, a.inner_tolerance,
    )
}

fn gradient(report: &JsonValue, spec: &Spec) -> Result<LoadSpaceGradientCheck3> {
    let config = spec.load_space.ok_or_else(|| invalid("retained source has no load-space model"))?;
    let data = field(report, "gradient_check")?;
    let rows = field(data, "probes")?.as_array().ok_or_else(|| invalid("missing initial gradient gate"))?;
    if rows.is_empty() || rows.len() > 2 || integer(data, "evaluations")? != 1 + 2 * rows.len()
        || data.get("passed") != Some(&JsonValue::Bool(true))
        || real(data, "step")? != config.options.gradient_check.step
        || real(data, "relative_tolerance")? != config.options.gradient_check.relative_tolerance {
        return Err(invalid("retained initial gradient experiment is incomplete or changed"));
    }
    let probes = rows.iter().map(|r| -> Result<GradientProbe> {
        Ok(GradientProbe {
            direction: match r.str_field("direction") {
                Some("Increasing") => GradientDirection::Increasing,
                Some("Decreasing") => GradientDirection::Decreasing,
                _ => return Err(invalid("unknown gradient direction")),
            },
            active_densities: integer(r, "active_densities")?,
            compliance_analytic: real(r, "spectral_analytic")?,
            compliance_difference: real(r, "spectral_difference")?,
            compliance_relative_error: real(r, "spectral_relative_error")?,
            volume_analytic: real(r, "material_analytic")?,
            volume_difference: real(r, "material_difference")?,
            volume_relative_error: real(r, "material_relative_error")?,
        })
    }).collect::<Result<Vec<_>>>()?;
    let w = field(data, "work")?;
    Ok(LoadSpaceGradientCheck3 {
        policy: config.cap(spec.loads.len())?.policy().clone(),
        robustness: config.options.robustness,
        numerical: MultiLoadGradientCheck {
            params: spec.schedule[0], options: config.options.gradient_check,
            baseline_compliance: real(data, "baseline_spectral")?,
            baseline_volume_fraction: real(data, "baseline_material")?, probes,
            work: SolveWork { linear_iterations: integer(w, "linear_iterations")?,
                linear_solves: integer(w, "linear_solves")?,
                preconditioner_operator_applications: integer(w, "preconditioner_operator_applications")?,
                preconditioner_galerkin_products: integer(w, "preconditioner_galerkin_products")? },
        },
    })
}

fn recover(spec: &Spec, ledger: &Ledger, old: &Loaded, producer: ContentHash) -> Result<Recovery> {
    let receipt = &old.value;
    if receipt.str_field("resume_mode") != Some(MODE)
        || receipt.get("resume_supported") != Some(&JsonValue::Bool(true)) {
        return Err(fail("cli-study-sdf3-resume-unsupported",
            "this historical receipt has no load-space optimizer checkpoint; report/package remain available"));
    }
    if receipt.str_field("driver") != Some(LOAD_SPACE3_DRIVER)
        || receipt.str_field("study_id") != Some(spec.id.to_hex().as_str())
        || receipt.str_field("producer") != Some(producer.to_hex().as_str())
        || integer(receipt, "target_iterations")? != spec.updates {
        return Err(invalid("load-space checkpoint source, executable or original target changed"));
    }
    let hash = receipt.str_field("checkpoint").and_then(ContentHash::from_hex)
        .ok_or_else(|| invalid("missing load-space optimizer artifact"))?;
    let op = ledger.artifact_output_seal(&old.hash)?.ok_or_else(|| invalid("unsealed load-space receipt"))?;
    if !ledger.edge_exists(op, &hash, EdgeRole::Out)? {
        return Err(invalid("load-space optimizer artifact is not a sealed output"));
    }
    let data = document(&artifact(ledger, hash, KIND)?)?;
    if data.str_field("schema") != Some(SCHEMA)
        || data.str_field("study_id") != receipt.str_field("study_id")
        || data.str_field("producer") != receipt.str_field("producer")
        || data.str_field("design") != receipt.str_field("design")
        || data.str_field("iterations") != receipt.str_field("iterations") {
        return Err(invalid("load-space checkpoint does not bind its source, fields and history"));
    }
    let point = field(&data, "point")?.as_array().ok_or_else(|| invalid("missing optimizer point"))?;
    let n = point.len();
    if n == 0 || n > spec.leaves { return Err(invalid("optimizer point exceeds its original cell allowance")); }
    let sample = field(&data, "sample")?;
    let optimizer = ProjectedAlCheckpoint {
        point: vector(field(&data, "point")?, n)?,
        sample: ProjectedAlSample { objective: real(sample, "objective")?,
            constraint: real(sample, "constraint")?, gradient: vector(field(sample, "gradient")?, n)?,
            constraint_gradient: vector(field(sample, "constraint_gradient")?, n)? },
        work: work(field(&data, "work")?)?, multiplier: real(&data, "multiplier")?,
        penalty: real(&data, "penalty")?, spectral_step: real(&data, "spectral_step")?,
        inner_tolerance: real(&data, "inner_tolerance")?,
        previous_violation: match field(&data, "previous_violation")? {
            JsonValue::Null => None, _ => Some(real(&data, "previous_violation")?),
        },
    };
    let counters = field(receipt, "optimizer_work")?;
    let restorations = integer(&data, "restoration_evaluations")?;
    if optimizer.work != work(counters)?
        || optimizer.work.iterations != integer(receipt, "iterations_completed")?
        || optimizer.work.iterations > spec.updates
        || restorations != integer(counters, "restoration_evaluations")? {
        return Err(invalid("retained load-space optimizer counters disagree"));
    }
    let options = spec.load_space.ok_or_else(|| invalid("retained source has no load-space optimizer"))?.options;
    ProjectedAlState::try_restore::<EvaluationStop>(optimizer.clone(), &vec![options.density_floor; n],
        &vec![1.0; n], options.optimizer).map_err(|e| invalid(e.to_string()))?;
    let report = document(&linked(ledger, receipt, "report_json", "study-report-json")?)?;
    let gradient_check = gradient(&report, spec)?;
    let rows = document(&linked(ledger, receipt, "iterations", "study-iterations")?)?;
    if rows.str_field("schema") != Some("sdf3-load-space-iterations.v1") {
        return Err(invalid("unsupported load-space history"));
    }
    let rows = field(&rows, "history")?.as_array().ok_or_else(|| invalid("missing load-space history"))?;
    if rows.len() != optimizer.work.iterations + 1 { return Err(invalid("incomplete load-space history")); }
    let history = rows.iter().map(|r| -> Result<LoadSpaceDesignIteration3> {
        Ok(LoadSpaceDesignIteration3 { iteration: integer(r, "iteration")?,
            smooth_maximum: real(r, "smooth_maximum_ratio")?, maximum: real(r, "maximum_ratio")?,
            dilated_volume: real(r, "dilated_volume_fraction")?,
            feasible: match field(r, "feasible")? {
                JsonValue::Bool(v) => *v, _ => return Err(invalid("missing feasibility state")),
            },
        })
    }).collect::<Result<Vec<_>>>()?;
    let checkpoint = LoadSpaceDesignCheckpoint3 { optimizer, history, gradient_check,
        restoration_evaluations: restorations,
        best_feasible_density: match field(&data, "best_feasible_density")? {
            JsonValue::Null => None, value => Some(vector(value, n)?),
        },
    };
    let spent = Spent::read(receipt)?;
    spent.validate(spec)?;
    Ok(Recovery { checkpoint, spent })
}

pub(super) fn drive(spec: &Spec, ledger: &Ledger, cap: Option<usize>, gate: &CancelGate,
    old: Option<&Loaded>) -> Result<Outcome> {
    let admission = Instant::now();
    let producer = checkpoint::producer_identity(gate)?;
    let recovery = old.map(|old| recover(spec, ledger, old, producer)).transpose()?;
    let limit = cap.unwrap_or(spec.updates).min(spec.updates);
    if let (Some(old), Some(recovery)) = (old, &recovery) {
        let status = match old.value.str_field("status") {
            Some("completed") => "completed", Some("no-feasible-descent") => "no-feasible-descent",
            Some("budget-exhausted") => "budget-exhausted", Some("cancelled") => "cancelled",
            Some("checkpointed") => "checkpointed", Some("numerical-failure") => "numerical-failure",
            _ => return Err(invalid("unknown retained load-space status")),
        };
        if matches!(status, "completed" | "no-feasible-descent")
            || recovery.checkpoint.optimizer.work.iterations == spec.updates {
            return Ok(Outcome { pointer: format!("study-{}", old.hash.to_hex()), receipt: old.bytes.clone(), status });
        }
        recovery.spent.require_remaining(spec)?;
        let maximum = spec.load_space.unwrap().options.optimizer.max_evaluations;
        if recovery.checkpoint.restoration_cost() > maximum - recovery.checkpoint.optimizer.work.evaluations {
            return Err(Failure { code: "cli-study-sdf3-resume-budget",
                message: "the original load-space evaluation allowance cannot rebuild the retained endpoints".into(),
                exit: exit::BUDGET });
        }
    }
    let mut previous = old.map(|old| old.hash);
    let run = compute_observed(spec, gate, limit, recovery.as_ref(), admission.elapsed().as_secs_f64(), |study, state| {
        let out = output::persist(spec, ledger, study, state, producer, previous, limit)?;
        previous = out.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        use std::io::Write as _;
        let _ = writeln!(std::io::stderr().lock(),
            "{{\"schema\":\"frankensim.cli.sdf3-load-space-progress.v1\",\"run_id\":{},\"iterations_completed\":{}}}",
            quoted(&out.pointer), state.optimizer.iterations);
        Ok(())
    }).map_err(|mut error| {
        if let Some(previous) = previous {
            error.message.push_str(&format!("; last durable result: study-{}", previous.to_hex()));
        }
        error
    })?;
    output::persist(spec, ledger, &run.study, &run.state, producer, previous, limit)
}

#[cfg(test)]
mod tests;
