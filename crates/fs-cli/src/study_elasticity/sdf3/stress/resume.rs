//! Restart the accepted stress optimizer; only its physical endpoints re-solve.
use super::*;
use fs_ascent::projected_al::{ProjectedAlCheckpoint, ProjectedAlSample, ProjectedAlState};

pub(super) const KIND: &str = "study-sdf3-stress-checkpoint";
pub(super) const MODE: &str = "optimizer-state-restoration-v1";
const SCHEMA: &str = "sdf3-stress-state.v1";

pub(super) struct Recovery {
    pub(super) checkpoint: StressDesignCheckpoint3,
    pub(super) audit: Audit,
    pub(super) spent: Spent,
}

fn invalid(message: impl Into<String>) -> Failure {
    fail("cli-study-sdf3-stress-checkpoint", message)
}
fn field<'a>(value: &'a JsonValue, key: &str) -> Result<&'a JsonValue> {
    value.get(key).ok_or_else(|| invalid(format!("missing checkpoint {key}")))
}
fn real(value: &JsonValue, key: &str) -> Result<f64> {
    field(value, key)?.as_f64().filter(|x| x.is_finite())
        .ok_or_else(|| invalid(format!("nonfinite checkpoint {key}")))
}
fn vector(value: &JsonValue, n: usize) -> Result<Vec<f64>> {
    let values = value.as_array().ok_or_else(|| invalid("expected a density or gradient array"))?;
    if values.len() != n { return Err(invalid("checkpoint vector has the wrong cell count")); }
    values.iter().map(|v| v.as_f64().filter(|v| v.is_finite())
        .ok_or_else(|| invalid("nonfinite checkpoint vector"))).collect()
}
fn document(bytes: &[u8]) -> Result<JsonValue> {
    JsonValue::parse(std::str::from_utf8(bytes).map_err(|e| invalid(e.to_string()))?)
        .map_err(|e| invalid(e.to_string()))
}
fn work(value: &JsonValue) -> Result<ProjectedAlWork> {
    Ok(ProjectedAlWork {
        iterations: integer(value, "iterations")?, evaluations: integer(value, "evaluations")?,
        multiplier_updates: integer(value, "multiplier_updates")?,
        rejected_trials: integer(value, "rejected_trials")?,
    })
}

pub(super) fn encode(spec: &Spec, state: &State, producer: ContentHash,
    design: &str, iterations: &str) -> Option<String> {
    let c = state.checkpoint.as_ref()?;
    let a = &c.optimizer;
    let previous = a.previous_violation.map_or_else(|| "null".into(), |v| format!("{v:.17e}"));
    let best = c.best_feasible_density.as_ref().map_or_else(|| "null".into(), |v| format!("{v:?}"));
    Some(format!(
        "{{\"schema\":{SCHEMA:?},\"study_id\":{:?},\"producer\":{:?},\"design\":{:?},\"iterations\":{:?},\"restoration_evaluations\":{},\"point\":{:?},\"best_feasible_density\":{best},\"sample\":{{\"objective\":{:.17e},\"constraint\":{:.17e},\"gradient\":{:?},\"constraint_gradient\":{:?}}},\"work\":{{\"iterations\":{},\"evaluations\":{},\"multiplier_updates\":{},\"rejected_trials\":{}}},\"multiplier\":{:.17e},\"penalty\":{:.17e},\"spectral_step\":{:.17e},\"inner_tolerance\":{:.17e},\"previous_violation\":{previous}}}",
        spec.id.to_hex(), producer.to_hex(), hash_bytes(design.as_bytes()).to_hex(),
        hash_bytes(iterations.as_bytes()).to_hex(), c.restoration_evaluations, a.point,
        a.sample.objective, a.sample.constraint, a.sample.gradient, a.sample.constraint_gradient,
        a.work.iterations, a.work.evaluations, a.work.multiplier_updates, a.work.rejected_trials,
        a.multiplier, a.penalty, a.spectral_step, a.inner_tolerance,
    ))
}

fn audit(report: &JsonValue, max_evaluations: usize) -> Result<Audit> {
    let data = field(report, "gradient_check")?;
    let rows = field(data, "probes")?.as_array().ok_or_else(|| invalid("missing initial gradient gate"))?;
    let evaluations = integer(data, "evaluations")?;
    if rows.is_empty() || rows.len() > 2 || evaluations != 1 + 2 * rows.len()
        || evaluations >= max_evaluations || data.get("passed") != Some(&JsonValue::Bool(true))
        || real(data, "step")? != 1e-4 || real(data, "relative_tolerance")? != 5e-4
    { return Err(invalid("initial gradient gate is incomplete or exceeds its original allowance")); }
    let probes = rows.iter().map(|r| -> Result<Probe> {
        Ok(Probe {
            direction: match r.str_field("direction") {
                Some("increasing") => "increasing", Some("decreasing") => "decreasing",
                _ => return Err(invalid("unknown initial gradient direction")),
            },
            active: integer(r, "active_densities")?,
            stress_analytic: real(r, "stress_analytic")?, stress_difference: real(r, "stress_difference")?,
            stress_relative_error: real(r, "stress_relative_error")?,
            volume_analytic: real(r, "volume_analytic")?, volume_difference: real(r, "volume_difference")?,
            volume_relative_error: real(r, "volume_relative_error")?,
        })
    }).collect::<Result<Vec<_>>>()?;
    if probes.iter().any(|p| p.active == 0 || !(0.0..=5e-4).contains(&p.stress_relative_error)
        || !(0.0..=5e-4).contains(&p.volume_relative_error))
    { return Err(invalid("retained initial gradient comparison did not pass")); }
    Ok(Audit { probes, evaluations, passed: true })
}

fn recover(spec: &Spec, ledger: &Ledger, old: &Loaded, producer: ContentHash) -> Result<Recovery> {
    let receipt = &old.value;
    if receipt.str_field("resume_mode") != Some(MODE)
        || receipt.get("resume_supported") != Some(&JsonValue::Bool(true))
    {
        return Err(fail("cli-study-sdf3-resume-unsupported",
            "this historical receipt has no complete stress optimizer checkpoint; report/package remain available"));
    }
    if receipt.str_field("driver") != Some(STRESS3_DRIVER)
        || receipt.str_field("study_id") != Some(spec.id.to_hex().as_str())
        || receipt.str_field("producer") != Some(producer.to_hex().as_str())
        || integer(receipt, "target_iterations")? != spec.updates
    { return Err(invalid("stress checkpoint source, executable or original target changed")); }
    let hash = receipt.str_field("checkpoint").and_then(ContentHash::from_hex)
        .ok_or_else(|| invalid("missing stress optimizer artifact"))?;
    let op = ledger.artifact_output_seal(&old.hash)?.ok_or_else(|| invalid("unsealed stress receipt"))?;
    if !ledger.edge_exists(op, &hash, EdgeRole::Out)? {
        return Err(invalid("stress optimizer artifact is not a sealed output"));
    }
    let data = document(&artifact(ledger, hash, KIND)?)?;
    if data.str_field("schema") != Some(SCHEMA)
        || data.str_field("study_id") != receipt.str_field("study_id")
        || data.str_field("producer") != receipt.str_field("producer")
        || data.str_field("design") != receipt.str_field("design")
        || data.str_field("iterations") != receipt.str_field("iterations")
    { return Err(invalid("stress state does not bind its source, fields and history")); }
    let point = field(&data, "point")?.as_array().ok_or_else(|| invalid("missing optimizer point"))?;
    let n = point.len();
    if n == 0 || n > spec.leaves { return Err(invalid("optimizer point exceeds the original cell allowance")); }
    let sample = field(&data, "sample")?;
    let optimizer = ProjectedAlCheckpoint {
        point: vector(field(&data, "point")?, n)?,
        sample: ProjectedAlSample {
            objective: real(sample, "objective")?, constraint: real(sample, "constraint")?,
            gradient: vector(field(sample, "gradient")?, n)?,
            constraint_gradient: vector(field(sample, "constraint_gradient")?, n)?,
        },
        work: work(field(&data, "work")?)?, multiplier: real(&data, "multiplier")?,
        penalty: real(&data, "penalty")?, spectral_step: real(&data, "spectral_step")?,
        inner_tolerance: real(&data, "inner_tolerance")?,
        previous_violation: match field(&data, "previous_violation")? {
            JsonValue::Null => None, _ => Some(real(&data, "previous_violation")?),
        },
    };
    if optimizer.work != work(field(receipt, "optimizer_work")?)?
        || optimizer.work.iterations != integer(receipt, "iterations_completed")?
        || optimizer.work.iterations > spec.updates
    { return Err(invalid("retained optimizer counters disagree")); }
    let options = spec.stress.ok_or_else(|| invalid("retained source is no longer a stress study"))?;
    let report = document(&linked(ledger, receipt, "report_json", "study-report-json")?)?;
    let audit = audit(&report, options.optimizer.max_evaluations)?;
    if audit.probes.iter().any(|p| p.active > n) { return Err(invalid("gradient gate exceeds the cell count")); }
    let mut policy = options.optimizer;
    policy.max_evaluations -= audit.evaluations;
    ProjectedAlState::try_restore::<StressError3>(optimizer.clone(), &vec![options.density_floor; n],
        &vec![1.0; n], policy).map_err(|e| invalid(e.to_string()))?;
    let rows = document(&linked(ledger, receipt, "iterations", "study-iterations")?)?;
    let history = field(&rows, "history")?.as_array().ok_or_else(|| invalid("missing stress history"))?;
    if history.len() != optimizer.work.iterations + 1 { return Err(invalid("stress history is incomplete")); }
    let history = history.iter().map(|r| -> Result<StressDesignIteration3> {
        Ok(StressDesignIteration3 {
            iteration: integer(r, "iteration")?, volume_fraction: real(r, "volume_fraction")?,
            stress_aggregate: real(r, "stress_aggregate_pa")?,
            sampled_relaxed_max: real(r, "sampled_relaxed_max_pa")?,
            sampled_physical_max: real(r, "sampled_physical_max_pa")?,
            constraint_violation: real(r, "constraint_violation")?,
            feasible: match field(r, "feasible")? {
                JsonValue::Bool(value) => *value, _ => return Err(invalid("missing feasibility state")),
            },
        })
    }).collect::<Result<Vec<_>>>()?;
    let checkpoint = StressDesignCheckpoint3 {
        optimizer, history, restoration_evaluations: integer(&data, "restoration_evaluations")?,
        best_feasible_density: match field(&data, "best_feasible_density")? {
            JsonValue::Null => None, value => Some(vector(value, n)?),
        },
    };
    let spent = Spent::read(receipt)?;
    spent.validate(spec)?;
    Ok(Recovery { checkpoint, audit, spent })
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
            _ => return Err(invalid("unknown retained stress status")),
        };
        if matches!(status, "completed" | "no-feasible-descent")
            || recovery.checkpoint.optimizer.work.iterations == spec.updates {
            return Ok(Outcome { pointer: format!("study-{}", old.hash.to_hex()), receipt: old.bytes.clone(), status });
        }
        recovery.spent.require_remaining(spec)?;
        let used = recovery.audit.evaluations + recovery.checkpoint.optimizer.work.evaluations;
        if recovery.checkpoint.restoration_cost() > spec.stress.unwrap().optimizer.max_evaluations - used {
            return Err(Failure { code: "cli-study-sdf3-resume-budget",
                message: "the original stress evaluation allowance cannot rebuild the retained endpoints".into(), exit: exit::BUDGET });
        }
    }
    let mut predecessor = old.map(|old| old.hash);
    let result = compute_observed(spec, gate, limit, recovery.as_ref(), admission.elapsed().as_secs_f64(), |study, state| {
        let out = output::persist(spec, ledger, study, state, producer, predecessor, limit)?;
        predecessor = out.pointer.strip_prefix("study-").and_then(ContentHash::from_hex);
        use std::io::Write as _;
        let _ = writeln!(std::io::stderr().lock(),
            "{{\"schema\":\"frankensim.cli.sdf3-stress-progress.v1\",\"run_id\":{},\"iterations_completed\":{}}}",
            quoted(&out.pointer), state.iterations());
        Ok(())
    }).map_err(|mut e| {
        if let Some(previous) = predecessor {
            e.message.push_str(&format!("; last durable result: study-{}", previous.to_hex()));
        }
        e
    })?;
    output::persist(spec, ledger, &result.study, &result.state, producer, predecessor, limit)
}

#[test]
fn g4_stress_checkpoint_binds_the_original_source_budget_and_executable() {
    let source = include_str!("../../../../../../examples/marquee/bracket-3d-stress.fsim");
    let spec = spec::parse(source).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let first = drive(&spec, &ledger, Some(1), &gate, None).unwrap();
    let old = load(&ledger, &first.pointer).unwrap();
    let producer = checkpoint::producer_identity(&gate).unwrap();
    recover(&spec, &ledger, &old, producer).unwrap();
    for changed in [
        source.replace(":stress-limit-pa 8.0", ":stress-limit-pa 9.0"),
        source.replace(":max-evaluations 2000", ":max-evaluations 3000"),
        source.replace(":linear-iterations 250000", ":linear-iterations 300000"),
    ] {
        assert!(recover(&spec::parse(&changed).unwrap(), &ledger, &old, producer).is_err());
    }
    assert!(recover(&spec, &ledger, &old, hash_bytes(b"different executable")).is_err());
    assert_eq!(load(&ledger, &first.pointer).unwrap().bytes, old.bytes);
}
