//! Concurrent conditional pairs through the original cooling child owner.
use super::*;
use fs_uq::product_shapley::ShapleyBatchExecution;
use std::io::Read as _;

pub(super) fn execute(
    worker: &Path, base_text: &str, base: &J, config: &Config, options: &Options, workers: usize,
) -> Result<ExecutionOutput> {
    let plan = config.plan();
    let mut run = ShapleyBatchExecution::new(&plan).map_err(bad)?;
    let identity = if options.checkpoint.is_some() || options.resume.is_some() {
        Some(checkpoint::model_identity(worker, base_text, &config.render_parameters()?)?)
    } else { None };
    if let Some(path) = &options.resume {
        let cap = 4096 + 8 * MAX_PRODUCT_SAMPLES as u64;
        let file = std::fs::File::open(path).map_err(|e| recovery_error(format!("{}: {e}", path.display())))?;
        let mut bytes = Vec::new();
        file.take(cap + 1).read_to_end(&mut bytes).map_err(|e| recovery_error(e.to_string()))?;
        if bytes.len() as u64 > cap { return Err(recovery_error("Shapley batch checkpoint exceeds its observation cap")); }
        run = ShapleyBatchExecution::restore(&plan, identity.expect("resume identity"), &bytes)
            .map_err(|e| recovery_error(format!("{}: {e}", path.display())))?;
    }
    // Re-admit the original conditional law and model BEFORE reserving output.
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    let save = |run: &ShapleyBatchExecution| -> Result<()> {
        if let Some(output) = &output {
            let bytes = run.checkpoint(identity.expect("output identity"))
                .map_err(|e| recovery_error(e.to_string()))?;
            output.publish(&bytes)?;
        }
        Ok(())
    };
    save(&run)?;
    let initial = run.evaluations_attempted();
    let allowance = options.max_new_samples.unwrap_or(config.samples).min(config.samples - initial);
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    let mut report = run.report();
    while report.status != UqStatus::Complete && run.evaluations_attempted() - initial < allowance {
        let width = workers.min(allowance - (run.evaluations_attempted() - initial));
        report = run.advance_batch::<String>(width, || Instant::now() >= deadline, |jobs, slots| {
            drain_wave(jobs, slots, &|job| evaluate_job(worker, base, config, job, deadline));
        }).map_err(bad)?;
        if report.status == UqStatus::Refused {
            let refusal = model_failure(report.rejection_reason.clone()
                .unwrap_or_else(|| "parallel Shapley cooling evaluation refused".into()));
            if let Some(output) = &output {
                output.invalidate(&refusal.to_string()).map_err(|mut error| {
                    error.message.push_str(&format!("; original failure: {refusal}")); error
                })?;
            }
            return Err(refusal);
        }
        // The shared launcher has joined all workers and drained their children.
        save(&run)?;
        if report.status == UqStatus::Cancelled { break; }
    }
    if let Some(result) = &report.result {
        return Ok(ExecutionOutput {
            stdout: super::super::shapley::render(base, config, result)?, exit_code: exit::SUCCESS,
        });
    }
    let Some(path) = &options.checkpoint else {
        return Err(budget(format!("parallel Shapley stopped after {} accepted physical solves; no partial attribution; use --checkpoint to retain successful conditional siblings", report.evaluations_accepted)));
    };
    let next = report.contiguous_observations;
    let slot = next % report.evaluations_per_row;
    let termination = if report.status == UqStatus::Cancelled { "wall-time-budget" } else { "sample-chunk" };
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.shapley.parallel.progress.v1\",\"status\":\"budget-truncated\",\"termination\":{},\"qoi\":{},\"workers\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"samples_evaluated_this_run\":{},\"base_samples\":{},\"evaluations_per_row\":{},\"completed_rows\":{},\"contiguous_observations\":{},\"retained_out_of_order\":{},\"next_sample_ordinal\":{},\"next_row_ordinal\":{},\"next_subset_mask\":{},\"next_pair_half\":{},\"checkpoint\":{},\"no_claim\":\"retained original conditional pairs including later successful siblings; no selected-worker or shortened-design attribution; resume the same physical model, joint law and allocation, optionally changing worker count\"}}\n",
        quote(termination), config.qoi.render(objective_kind(base)?), workers,
        report.evaluations_planned, report.evaluations_attempted, report.evaluations_attempted - initial,
        report.base_samples, report.evaluations_per_row, report.completed_rows,
        next, report.retained_out_of_order, next, next / report.evaluations_per_row,
        slot / 2, slot % 2, quote(&path.to_string_lossy()),
    );
    Ok(ExecutionOutput { stdout, exit_code: exit::BUDGET })
}
