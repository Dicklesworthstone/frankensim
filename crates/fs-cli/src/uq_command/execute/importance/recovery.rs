//! Fixed-proposal cooling recovery using the original atomic file transport.
use super::*;
use std::io::Read as _;

fn failure(message: impl Into<String>) -> Failure {
    Failure { code: "cooling-network-uq-checkpoint", message: message.into() }
}

pub(super) fn execute(
    worker: &Path, base_text: &str, base: &J, config: &Config, options: &Options, proposal: &ImportanceConfig,
) -> Result<ExecutionOutput> {
    let plan = config.plan();
    if let Some(policy) = options.compliance { policy.validate_plan(&plan)?; }
    let mut execution = ImportanceExecution::new(&plan, proposal.clone()).map_err(bad)?;
    let identity = if options.checkpoint.is_some() || options.resume.is_some() {
        let parameters = config.render_parameters()?;
        let binding = options.compliance.map_or_else(|| parameters.clone(),
            |policy| decision::binding(policy, &parameters));
        Some(checkpoint::model_identity(worker, base_text, &binding)?)
    } else { None };
    if let Some(path) = &options.resume {
        // Importance records retain QoI AND scaled weight. Keep the existing
        // scalar transport's smaller read cap unchanged for its other formats.
        let cap = 4096 + 16 * MAX_PRODUCT_SAMPLES as u64;
        let file = std::fs::File::open(path).map_err(|e| failure(format!("{}: {e}", path.display())))?;
        let mut bytes = Vec::new();
        file.take(cap + 1).read_to_end(&mut bytes).map_err(|e| failure(e.to_string()))?;
        if bytes.len() as u64 > cap { return Err(failure("importance checkpoint exceeds its observation cap")); }
        execution = ImportanceExecution::restore(&plan, proposal.clone(), identity.expect("resume identity"), &bytes)
            .map_err(|e| failure(format!("{}: {e}", path.display())))?;
    }
    let mut assessment = options.compliance
        .map(|policy| decision::assess(policy, Some(&execution))).transpose()?;
    // Decode and check model/proposal/policy identity BEFORE reserving any output.
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    let save = |execution: &ImportanceExecution| -> Result<()> {
        if let Some(output) = &output {
            let bytes = execution.checkpoint(identity.expect("output identity")).map_err(|e| failure(e.to_string()))?;
            output.publish(&bytes)?;
        }
        Ok(())
    };
    save(&execution)?;
    let initial = execution.observations().len();
    let allowance = options.max_new_samples.unwrap_or(config.samples).min(config.samples - initial);
    // A fresh evaluation-time allowance does not refund completed model calls.
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    for _ in 0..allowance {
        if assessment.as_ref().is_some_and(decision::Assessment::reached) { break; }
        let report = execution.advance_interruptible(1, || Instant::now() >= deadline,
            |values| evaluate(worker, base, config, values, deadline));
        if report.status == UqStatus::Refused {
            let refusal = model_failure(report.rejection_reason.unwrap_or_else(|| "importance physical solve refused".into()));
            if let Some(output) = &output {
                output.invalidate(&refusal.to_string()).map_err(|mut error| {
                    error.message.push_str(&format!("; original failure: {refusal}")); error
                })?;
            }
            return Err(refusal);
        }
        save(&execution)?;
        assessment = options.compliance
            .map(|policy| decision::assess(policy, Some(&execution))).transpose()?;
        if matches!(report.status, UqStatus::Complete | UqStatus::Cancelled) { break; }
    }
    let report = execution.report();
    if let (Some(policy), Some(assessment)) = (options.compliance, assessment) {
        let termination = decision::termination(&assessment, report.status);
        let stdout = format!(
            "{{\"schema\":\"frankensim.cooling-network-uq.importance-compliance.v1\",\"sampling\":\"fixed-defensive-gaussian\",\"qoi\":{},\"seed\":{},\"proposal\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"production_samples_evaluated\":{},\"temperature_limit_k\":{},\"correlation\":{},\"parameters\":[{}],{}}}\n",
            config.qoi.render(objective_kind(base)?), quote(&config.seed.to_string()),
            proposal_json(proposal)?, config.samples, report.evaluations_attempted,
            report.evaluations_accepted, number_json(report.threshold)?,
            quote(config.correlation_label), config.render_parameters()?,
            decision::fields(config, policy, &assessment, termination, options.checkpoint.as_deref())?,
        );
        return Ok(ExecutionOutput { stdout, exit_code: decision::exit_code(&assessment) });
    }
    if report.status == UqStatus::Complete {
        // Completed recovery performs no model work, even with a tiny new wall
        // allowance. The final schema/bytes match the non-recovery command.
        return Ok(ExecutionOutput { stdout: render(base, config, proposal, &report)?, exit_code: exit::SUCCESS });
    }
    let Some(path) = &options.checkpoint else {
        return Err(budget(format!("importance stopped after {} accepted evaluations; no partial probability published; use a fresh --checkpoint to retain further progress", report.evaluations_accepted)));
    };
    let termination = if report.status == UqStatus::Cancelled { "wall-time-budget" } else { "sample-chunk" };
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.importance.progress.v1\",\"status\":\"budget-truncated\",\"termination\":{},\"qoi\":{},\"seed\":{},\"proposal\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"samples_evaluated_this_run\":{},\"next_sample_ordinal\":{},\"checkpoint\":{},\"no_claim\":\"retained fixed-proposal observations and likelihood weights only; no shortened-design probability or confidence estimate; resume the identical model, proposal, law and original budget\"}}\n",
        quote(termination), config.qoi.render(objective_kind(base)?), quote(&config.seed.to_string()),
        proposal_json(proposal)?, report.evaluations_planned, report.evaluations_attempted,
        report.evaluations_attempted - initial, report.evaluations_accepted, quote(&path.to_string_lossy()),
    );
    Ok(ExecutionOutput { stdout, exit_code: exit::BUDGET })
}
