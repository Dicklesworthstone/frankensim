//! Restore the affine AL state and independently re-solve its physical endpoints.
use super::*;
use fs_ascent::projected_al::{ProjectedAlCheckpoint, ProjectedAlSample, ProjectedAlState};
use std::io::Write as _;

pub(super) const KIND: &str = "study-sdf3-affine-load-space-checkpoint";
pub(super) const MODE: &str = "optimizer-state-restoration-v1";
const SCHEMA: &str = "sdf3-affine-load-space-state.v1";

pub(super) struct Recovery {
    pub(super) checkpoint: AffineLoadSpaceDesignCheckpoint3,
    pub(super) spent: Spent,
    pub(super) optimization_status: &'static str,
    pub(super) optimizer_stop: Option<ProjectedAlStop>,
    pub(super) assessment_checked: bool,
    pub(super) refinement: Option<grid_refinement::Context>,
}
impl Recovery {
    pub(super) fn optimizer_finished(&self, spec: &Spec) -> bool {
        matches!(self.optimization_status, "completed" | "no-feasible-descent")
            || self.checkpoint.optimizer.work.iterations == spec.updates
    }
}

fn invalid(message: impl Into<String>) -> Failure {
    fail("cli-study-sdf3-affine-load-space-checkpoint", message)
}
pub(super) fn status(value: Option<&str>) -> Result<&'static str> {
    match value {
        Some("completed") => Ok("completed"),
        Some("no-feasible-descent") => Ok("no-feasible-descent"),
        Some("budget-exhausted") => Ok("budget-exhausted"),
        Some("cancelled") => Ok("cancelled"),
        Some("checkpointed") => Ok("checkpointed"),
        Some("numerical-failure") => Ok("numerical-failure"),
        Some("grid-unresolved") => Ok("grid-unresolved"),
        _ => Err(invalid("unknown retained affine status")),
    }
}
fn field<'a>(value: &'a JsonValue, key: &str) -> Result<&'a JsonValue> {
    value
        .get(key)
        .ok_or_else(|| invalid(format!("missing checkpoint {key}")))
}
fn real(value: &JsonValue, key: &str) -> Result<f64> {
    field(value, key)?
        .as_f64()
        .filter(|x| x.is_finite())
        .ok_or_else(|| invalid(format!("nonfinite checkpoint {key}")))
}
fn vector(value: &JsonValue, n: usize) -> Result<Vec<f64>> {
    let values = value
        .as_array()
        .ok_or_else(|| invalid("expected a numerical vector"))?;
    if values.len() != n {
        return Err(invalid("checkpoint vector has the wrong cell count"));
    }
    values
        .iter()
        .map(|v| {
            v.as_f64()
                .filter(|v| v.is_finite())
                .ok_or_else(|| invalid("nonfinite checkpoint vector"))
        })
        .collect()
}
fn document(bytes: &[u8]) -> Result<JsonValue> {
    JsonValue::parse(std::str::from_utf8(bytes).map_err(|e| invalid(e.to_string()))?)
        .map_err(|e| invalid(e.to_string()))
}
fn work(value: &JsonValue) -> Result<ProjectedAlWork> {
    Ok(ProjectedAlWork {
        iterations: integer(value, "iterations")?,
        evaluations: integer(value, "evaluations")?,
        multiplier_updates: integer(value, "multiplier_updates")?,
        rejected_trials: integer(value, "rejected_trials")?,
    })
}

pub(super) fn encode(
    spec: &Spec,
    state: &State,
    producer: ContentHash,
    design: &str,
    iterations: &str,
) -> String {
    let c = &state.checkpoint;
    let a = &c.optimizer;
    let previous = a
        .previous_violation
        .map_or_else(|| "null".into(), |v| format!("{v:.17e}"));
    let best = c
        .best_feasible_density
        .as_ref()
        .map_or_else(|| "null".into(), |v| format!("{v:?}"));
    let refinement = state.refinement.as_ref().map_or_else(String::new,
        |context| format!(",\"grid_refinement\":{}", context.json()));
    format!(
        "{{\"schema\":{SCHEMA:?},\"study_id\":{:?},\"producer\":{:?},\"binding\":{:?},\"design\":{:?},\"iterations\":{:?},\"restoration_evaluations\":{},\"point\":{:?},\"best_feasible_density\":{best},\"sample\":{{\"objective\":{:.17e},\"constraint\":{:.17e},\"gradient\":{:?},\"constraint_gradient\":{:?}}},\"work\":{{\"iterations\":{},\"evaluations\":{},\"multiplier_updates\":{},\"rejected_trials\":{}}},\"multiplier\":{:.17e},\"penalty\":{:.17e},\"spectral_step\":{:.17e},\"inner_tolerance\":{:.17e},\"previous_violation\":{previous}{refinement}}}",
        spec.id.to_hex(),
        producer.to_hex(),
        c.binding.to_hex(),
        hash_bytes(design.as_bytes()).to_hex(),
        hash_bytes(iterations.as_bytes()).to_hex(),
        c.restoration_evaluations,
        a.point,
        a.sample.objective,
        a.sample.constraint,
        a.sample.gradient,
        a.sample.constraint_gradient,
        a.work.iterations,
        a.work.evaluations,
        a.work.multiplier_updates,
        a.work.rejected_trials,
        a.multiplier,
        a.penalty,
        a.spectral_step,
        a.inner_tolerance,
    )
}

fn gradient(report: &JsonValue, spec: &Spec) -> Result<AffineLoadSpaceGradientCheck3> {
    let config = spec
        .affine_load_space
        .ok_or_else(|| invalid("retained source has no affine load-space model"))?;
    let data = field(report, "gradient_check")?;
    let directions = integer(data, "directions")?;
    if !(1..=2).contains(&directions)
        || data.get("passed") != Some(&JsonValue::Bool(true))
        || real(data, "step")? != config.options.gradient_check.step
        || real(data, "relative_tolerance")? != config.options.gradient_check.relative_tolerance
    {
        return Err(invalid(
            "retained initial gradient experiment is incomplete or changed",
        ));
    }
    let check = AffineLoadSpaceGradientCheck3 {
        params: spec.schedule[0],
        goal: GOAL,
        options: config.options.gradient_check,
        constraints: config.constraint_policy(),
        baseline_objective: real(data, "baseline_objective")?,
        baseline_constraint: real(data, "baseline_constraint")?,
        objective_relative_error: real(data, "objective_relative_error")?,
        constraint_relative_error: real(data, "constraint_relative_error")?,
        directions,
        work: SolveWork {
            linear_iterations: integer(data, "linear_iterations")?,
            linear_solves: integer(data, "linear_solves")?,
            preconditioner_operator_applications: integer(
                data,
                "preconditioner_operator_applications",
            )?,
            preconditioner_galerkin_products: integer(data, "preconditioner_galerkin_products")?,
        },
    };
    if !check.passed() {
        return Err(invalid("retained affine gradient experiment did not pass"));
    }
    Ok(check)
}

/// Reconstruct only the immutable source topology, never its physics or AL
/// history. This prevents a retained mark list from silently changing which
/// cells the resumed candidate enriches.
fn verify_refinement_source(
    spec: &Spec, ledger: &Ledger, old: &Loaded, producer: ContentHash,
    context: &grid_refinement::Context,
) -> Result<()> {
    let source = load(ledger, &format!("study-{}", context.source_receipt.to_hex()))?;
    let receipt = &source.value;
    let op = ledger.artifact_output_seal(&old.hash)?
        .ok_or_else(|| invalid("unsealed refinement receipt"))?;
    if source.hash == old.hash || !ledger.edge_exists(op, &source.hash, EdgeRole::In)?
        || receipt.str_field("driver") != Some(AFFINE_LOAD_SPACE3_DRIVER)
        || receipt.str_field("study_id") != Some(spec.id.to_hex().as_str())
        || receipt.str_field("producer") != Some(producer.to_hex().as_str())
        || receipt.str_field("source") != Some(hash_bytes(spec.canonical.as_bytes()).to_hex().as_str())
        || receipt.get("grid_refinement").is_some()
        || integer(receipt, "target_iterations")? != spec.updates
        || receipt.get("selected_feasible") != Some(&JsonValue::Bool(true))
        || (context.uses_source_selection() && old.value.str_field("design") != receipt.str_field("design"))
        || (!context.uses_source_selection()
            && old.value.str_field("design") != old.value.str_field("working_design"))
    {
        return Err(invalid("refinement source, selection or original grid identity changed"));
    }
    let report = document(&linked(ledger, receipt, "report_json", "study-report-json")?)?;
    let assessment = field(&report, "grid_assessment")?;
    let finished = match field(&report, "stop")?.str_field("reason") {
        Some("Converged") => true,
        Some("IterationLimit") => integer(receipt, "iterations_completed")? == spec.updates,
        _ => false,
    };
    if !finished || report.get("selected_feasible") != Some(&JsonValue::Bool(true))
        || assessment.str_field("status") != Some("checked")
        || assessment.get("passed") != Some(&JsonValue::Bool(false))
    {
        return Err(invalid("refinement source is not a completed feasible unresolved grid endpoint"));
    }
    let design = document(&linked(ledger, receipt, "design", "study-design")?)?;
    let cells = field(field(&design, "selected")?, "cells")?.as_array()
        .ok_or_else(|| invalid("missing source active cells"))?;
    if cells.is_empty() || cells.len() != context.source_marks.len() {
        return Err(invalid("refinement marks do not cover the source active cells"));
    }
    let tree = Octree3::uniform(spec.level as u8, spec.max_level as u8, spec.leaves)
        .map_err(|e| invalid(e.to_string()))?;
    let background: Vec<_> = tree.leaves().iter().copied().collect();
    let mut marks = Vec::with_capacity(cells.len());
    for cell in cells {
        if integer(cell, "level")? != spec.level as usize {
            return Err(invalid("refinement source is not the original uniform grid"));
        }
        let index = field(cell, "index")?.as_array()
            .filter(|v| v.len() == 3).ok_or_else(|| invalid("invalid source cell index"))?;
        let mut coordinates = [0u32; 3];
        for (coordinate, value) in coordinates.iter_mut().zip(index) {
            *coordinate = value.as_f64().filter(|v| v.is_finite() && *v >= 0.0
                && *v < (1u32 << spec.level) as f64 && v.fract() == 0.0)
                .ok_or_else(|| invalid("invalid source cell coordinate"))? as u32;
        }
        marks.push(background.iter().position(|leaf| leaf.index() == coordinates)
            .ok_or_else(|| invalid("source cell is absent from original background"))?);
    }
    marks.sort_unstable();
    if marks != context.source_marks {
        return Err(invalid("retained refinement marks changed the source topology"));
    }
    Ok(())
}

pub(super) fn recover(original: &Spec, ledger: &Ledger, old: &Loaded, producer: ContentHash) -> Result<Recovery> {
    let receipt = &old.value;
    let refinement = grid_refinement::Context::read(original, receipt)?;
    let effective = grid_refinement::effective_spec(original, refinement.as_ref())?;
    let spec = &effective;
    let design_key = if refinement.is_some() { "working_design" } else { "design" };
    if receipt.str_field("resume_mode") != Some(MODE)
        || receipt.get("resume_supported") != Some(&JsonValue::Bool(true))
    {
        return Err(fail(
            "cli-study-sdf3-resume-unsupported",
            "this historical receipt has no affine optimizer checkpoint; report/package remain available",
        ));
    }
    if receipt.str_field("driver") != Some(AFFINE_LOAD_SPACE3_DRIVER)
        || receipt.str_field("study_id") != Some(spec.id.to_hex().as_str())
        || receipt.str_field("producer") != Some(producer.to_hex().as_str())
        || integer(receipt, "target_iterations")? != spec.updates
    {
        return Err(invalid(
            "affine checkpoint source, executable or original target changed",
        ));
    }
    let hash = receipt
        .str_field("checkpoint")
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| invalid("missing affine optimizer artifact"))?;
    let op = ledger
        .artifact_output_seal(&old.hash)?
        .ok_or_else(|| invalid("unsealed affine receipt"))?;
    if !ledger.edge_exists(op, &hash, EdgeRole::Out)? {
        return Err(invalid("affine optimizer artifact is not a sealed output"));
    }
    let data = document(&artifact(ledger, hash, KIND)?)?;
    if data.str_field("schema") != Some(SCHEMA)
        || data.str_field("study_id") != receipt.str_field("study_id")
        || data.str_field("producer") != receipt.str_field("producer")
        || data.str_field("design") != receipt.str_field(design_key)
        || data.get("grid_refinement") != receipt.get("grid_refinement")
        || data.str_field("iterations") != receipt.str_field("iterations")
    {
        return Err(invalid(
            "affine checkpoint does not bind its source, fields and history",
        ));
    }
    // Validate the same sealed working artifact that the native checkpoint
    // binds. The public source design may belong to a different grid.
    let design_hash = receipt.str_field(design_key).and_then(ContentHash::from_hex)
        .ok_or_else(|| invalid("missing retained working design"))?;
    if !ledger.edge_exists(op, &design_hash, EdgeRole::Out)? {
        return Err(invalid("working design is not an output of the checkpoint operation"));
    }
    linked(ledger, receipt, design_key, "study-design")?;
    if let Some(context) = &refinement {
        verify_refinement_source(original, ledger, old, producer, context)?;
    }
    let binding = data
        .str_field("binding")
        .and_then(ContentHash::from_hex)
        .ok_or_else(|| invalid("missing native affine problem binding"))?;
    let n = field(&data, "point")?
        .as_array()
        .ok_or_else(|| invalid("missing optimizer point"))?
        .len();
    if n == 0 || n > spec.leaves {
        return Err(invalid(
            "optimizer point exceeds its original cell allowance",
        ));
    }
    let sample = field(&data, "sample")?;
    let optimizer = ProjectedAlCheckpoint {
        point: vector(field(&data, "point")?, n)?,
        sample: ProjectedAlSample {
            objective: real(sample, "objective")?,
            constraint: real(sample, "constraint")?,
            gradient: vector(field(sample, "gradient")?, n)?,
            constraint_gradient: vector(field(sample, "constraint_gradient")?, n)?,
        },
        work: work(field(&data, "work")?)?,
        multiplier: real(&data, "multiplier")?,
        penalty: real(&data, "penalty")?,
        spectral_step: real(&data, "spectral_step")?,
        inner_tolerance: real(&data, "inner_tolerance")?,
        previous_violation: match field(&data, "previous_violation")? {
            JsonValue::Null => None,
            _ => Some(real(&data, "previous_violation")?),
        },
    };
    let counters = field(receipt, "optimizer_work")?;
    let restorations = integer(&data, "restoration_evaluations")?;
    if optimizer.work != work(counters)?
        || optimizer.work.iterations != integer(receipt, "iterations_completed")?
        || optimizer.work.iterations > spec.updates
        || restorations != integer(counters, "restoration_evaluations")?
    {
        return Err(invalid("retained affine optimizer counters disagree"));
    }
    let options = spec
        .affine_load_space
        .ok_or_else(|| invalid("retained source has no affine optimizer"))?
        .options;
    ProjectedAlState::try_restore::<EvaluationStop>(
        optimizer.clone(),
        &vec![options.density_floor; n],
        &vec![1.0; n],
        options.optimizer,
    )
    .map_err(|e| invalid(e.to_string()))?;
    let report = document(&linked(
        ledger,
        receipt,
        "report_json",
        "study-report-json",
    )?)?;
    if grid_refinement::Context::read(original, &report)? != refinement {
        return Err(invalid("retained report and checkpoint grid contexts disagree"));
    }
    let optimization_status = status(receipt.str_field("optimization_status")
        .or_else(|| receipt.str_field("status")))?;
    let optimizer_stop = match field(&report, "stop")?.str_field("reason") {
        Some("Converged") => Some(ProjectedAlStop::Converged),
        Some("IterationLimit") => Some(ProjectedAlStop::IterationLimit),
        Some("EvaluationLimit") => Some(ProjectedAlStop::EvaluationLimit),
        Some("MultiplierLimit") => Some(ProjectedAlStop::MultiplierLimit),
        Some("PenaltyLimit") => Some(ProjectedAlStop::PenaltyLimit),
        Some("Stalled") => Some(ProjectedAlStop::Stalled),
        Some("EvaluationStopped") => None,
        _ => return Err(invalid("unknown retained optimizer stop")),
    };
    if optimization_status == "grid-unresolved"
        || report.str_field("optimization_status").is_some_and(|s| s != optimization_status)
    {
        return Err(invalid("retained optimizer and assessment termination disagree"));
    }
    let assessment_checked = if spec.affine_grid_assessment.is_some() {
        match field(&report, "grid_assessment")?.str_field("status") {
            Some("checked") => true,
            Some("pending" | "budget-exhausted" | "cancelled" | "numerical-failure") => false,
            _ => return Err(invalid("missing or unknown grid-assessment status")),
        }
    } else { false };
    let gradient_check = gradient(&report, spec)?;
    let rows = document(&linked(ledger, receipt, "iterations", "study-iterations")?)?;
    if rows.str_field("schema") != Some("sdf3-affine-load-space-iterations.v1") {
        return Err(invalid("unsupported affine history"));
    }
    let rows = field(&rows, "history")?
        .as_array()
        .ok_or_else(|| invalid("missing affine history"))?;
    if rows.len() != optimizer.work.iterations + 1 {
        return Err(invalid("incomplete affine history"));
    }
    let history = rows
        .iter()
        .map(|r| -> Result<AffineLoadSpaceDesignIteration3> {
            Ok(AffineLoadSpaceDesignIteration3 {
                iteration: integer(r, "iteration")?,
                objective: real(r, "objective")?,
                constraint: real(r, "constraint")?,
                maximum: real(r, "maximum_ratio")?,
                worst_compliance: real(r, "worst_compliance_j")?,
                dilated_volume_fraction: real(r, "dilated_volume_fraction")?,
                dilated_mass_kg: None,
                feasible: match field(r, "feasible")? {
                    JsonValue::Bool(v) => *v,
                    _ => return Err(invalid("missing feasibility state")),
                },
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let checkpoint = AffineLoadSpaceDesignCheckpoint3 {
        optimizer,
        binding,
        history,
        gradient_check,
        restoration_evaluations: restorations,
        best_feasible_density: match field(&data, "best_feasible_density")? {
            JsonValue::Null => None,
            value => Some(vector(value, n)?),
        },
    };
    let spent = Spent::read(receipt)?;
    spent.validate(spec)?;
    Ok(Recovery { checkpoint, spent, optimization_status, optimizer_stop, assessment_checked, refinement })
}

pub(super) fn drive(
    spec: &Spec,
    ledger: &Ledger,
    cap: Option<usize>,
    gate: &CancelGate,
    old: Option<&Loaded>,
) -> Result<Outcome> {
    let admission = Instant::now();
    let producer = checkpoint::producer_identity(gate)?;
    let recovery = old
        .map(|old| recover(spec, ledger, old, producer))
        .transpose()?;
    let limit = cap.unwrap_or(spec.updates).min(spec.updates);
    if let (Some(old), Some(recovery)) = (old, &recovery) {
        let status = status(old.value.str_field("status"))?;
        let assessment_done = spec.affine_grid_assessment.is_none() || recovery.assessment_checked;
        if recovery.optimizer_finished(spec) && assessment_done {
            return Ok(Outcome {
                pointer: format!("study-{}", old.hash.to_hex()),
                receipt: old.bytes.clone(),
                status,
            });
        }
        recovery.spent.require_remaining(spec)?;
        let maximum = spec
            .affine_load_space
            .unwrap()
            .options
            .optimizer
            .max_evaluations;
        if recovery.checkpoint.restoration_cost()
            > maximum.saturating_sub(recovery.checkpoint.optimizer.work.evaluations)
        {
            return Err(Failure {
                code: "cli-study-sdf3-resume-budget",
                message:
                    "the original affine evaluation allowance cannot rebuild the retained endpoints"
                        .into(),
                exit: exit::BUDGET,
            });
        }
    }
    let mut previous = old.map(|old| old.hash);
    let run = compute_observed(
        spec, gate, limit, recovery.as_ref(), admission.elapsed().as_secs_f64(),
        |study, state| {
            let out = output::persist(spec, ledger, study, state, producer, previous, limit)?;
            previous = out.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
            let _ = writeln!(std::io::stderr().lock(),
                "{{\"schema\":\"frankensim.cli.sdf3-affine-load-space-progress.v1\",\"run_id\":{},\"iterations_completed\":{}}}",
                quoted(&out.pointer), state.optimizer.iterations);
            Ok(())
        },
    )
    .map_err(|mut error| {
        if let Some(previous) = previous {
            error.message.push_str(&format!("; last durable result: study-{}", previous.to_hex()));
        }
        error
    })?;
    output::persist(
        spec, ledger, &run.study, &run.state, producer, previous, limit,
    )
}

#[cfg(test)]
#[path = "resume/tests.rs"]
mod tests;
