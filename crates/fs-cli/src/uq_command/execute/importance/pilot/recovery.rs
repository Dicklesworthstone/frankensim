//! Persist the entire learning study using the original atomic file transport.
use super::*;
use std::io::Read as _;

// At most 16 fits with 256 shifts, plus current QoI/weight observations.
// Other sampling formats retain their existing, smaller file-read allowances.
const MAX_BYTES: u64 = 65_536 + 16 * MAX_PRODUCT_SAMPLES as u64;

fn failure(message: impl Into<String>) -> Failure {
    Failure { code: "cooling-network-uq-checkpoint", message: message.into() }
}

pub(super) fn execute(
    base_text: &str, base: &J, config: &Config, options: &Options, settings: Settings,
) -> Result<ExecutionOutput> {
    let plan = config.plan();
    let policy = settings.library_config(plan.parameters.len());
    let mut execution = PilotImportanceExecution::new(&plan, policy.clone()).map_err(bad)?;
    let identity = checkpoint::model_identity(base_text, &config.render_parameters()?)?;
    if let Some(path) = &options.resume {
        let file = std::fs::File::open(path)
            .map_err(|e| failure(format!("{}: {e}", path.display())))?;
        let mut bytes = Vec::new();
        file.take(MAX_BYTES + 1).read_to_end(&mut bytes).map_err(|e| failure(e.to_string()))?;
        if bytes.len() as u64 > MAX_BYTES {
            return Err(failure("pilot importance checkpoint exceeds its observation/history cap"));
        }
        execution = PilotImportanceExecution::restore(&plan, policy, identity, &bytes)
            .map_err(|e| failure(format!("{}: {e}", path.display())))?;
    }
    // Validate the exact physical model, law AND learning policy before output.
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    let save = |execution: &PilotImportanceExecution| -> Result<()> {
        if let Some(output) = &output {
            let bytes = execution.checkpoint(identity).map_err(|e| failure(e.to_string()))?;
            output.publish(&bytes)?;
        }
        Ok(())
    };
    save(&execution)?;
    let initial = execution.evaluations_attempted();
    let allowance = options.max_new_samples.unwrap_or(config.samples).min(config.samples - initial);
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    let mut report = execution.report();
    loop {
        if report.status == UqStatus::Complete { break; }
        let paid = execution.evaluations_attempted();
        // A completed physical design can still await publication. Permit that
        // pure transition on resume, but explicit zero chunk remains a no-op.
        let finish_paid = paid == config.samples && options.max_new_samples != Some(0);
        if paid - initial >= allowance && !finish_paid { break; }
        report = execution.advance_interruptible(1, || Instant::now() >= deadline,
            |values| super::super::evaluate(base, config, values, deadline));
        if report.status == UqStatus::Refused {
            let refusal = model_failure(report.rejection_reason
                .unwrap_or_else(|| "importance pilot fit or physical solve refused".into()));
            if let Some(output) = &output {
                output.invalidate(&refusal.to_string()).map_err(|mut error| {
                    error.message.push_str(&format!("; original failure: {refusal}")); error
                })?;
            }
            return Err(refusal);
        }
        save(&execution)?;
        if matches!(report.status, UqStatus::Complete | UqStatus::Cancelled) { break; }
    }
    if report.status == UqStatus::Complete {
        // Reuse the existing production-only result; no report schema or
        // estimator change, including terminal restore under a tiny deadline.
        return Ok(ExecutionOutput {
            stdout: super::render(base, config, settings, &report)?, exit_code: exit::SUCCESS,
        });
    }
    let Some(path) = &options.checkpoint else {
        return Err(budget(format!(
            "importance pilot study stopped after {} paid calls; no partial probability published; use a fresh --checkpoint to retain new progress",
            report.model_calls_evaluated,
        )));
    };
    let (stage, round) = match report.stage {
        fs_uq::PilotImportanceStage::Pilot { round } => ("pilot", round.to_string()),
        fs_uq::PilotImportanceStage::Production => ("production", "null".into()),
        fs_uq::PilotImportanceStage::Complete => unreachable!("complete report returned above"),
    };
    let termination = if report.status == UqStatus::Cancelled { "wall-time-budget" } else { "sample-chunk" };
    let pending = execution.current().report().status == UqStatus::Complete;
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.importance-pilot.progress.v1\",\"status\":\"budget-truncated\",\"termination\":{},\"qoi\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"samples_evaluated_this_run\":{},\"remaining_physical_calls\":{},\"stage\":{},\"pilot_round\":{},\"completed_pilots\":{},\"current_seed\":{},\"current_proposal\":{},\"current_evaluations\":{},\"phase_work_pending\":{},\"checkpoint\":{},\"no_claim\":\"retained pilot history and current physical observations only; no shortened-design or pooled-pilot probability; resume the same model, laws, fitting policy and total budget\"}}\n",
        quote(termination), config.qoi.render(objective_kind(base)?),
        report.model_calls_planned, report.model_calls_evaluated,
        report.model_calls_evaluated - initial, config.samples - report.model_calls_evaluated,
        quote(stage), round, report.completed_pilots.len(), quote(&report.current_seed.to_string()),
        proposal_json(&report.current_proposal)?, report.current_evaluations, pending,
        quote(&path.to_string_lossy()),
    );
    Ok(ExecutionOutput { stdout, exit_code: exit::BUDGET })
}
