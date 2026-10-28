//! Parallel physical QMC using the original cooling workers and result renderer.
use super::*;
use fs_uq::product_qmc::QmcBatchExecution;
use fs_uq::QmcConfig;
use std::io::Read as _;

pub(super) fn execute(
    worker: &Path, base_text: &str, base: &J, config: &Config, options: &Options,
    workers: usize, replicates: usize,
) -> Result<ExecutionOutput> {
    if replicates == 0 || config.samples % replicates != 0 {
        return Err(bad("samples must divide exactly into the declared QMC replicates"));
    }
    let layout = QmcConfig { replicates, samples_per_replicate: config.samples / replicates };
    let mut plan = config.plan();
    plan.method = PropagationMethod::QuasiMonteCarlo;
    let mut run = QmcBatchExecution::new(&plan, layout).map_err(bad)?;
    let identity = if options.checkpoint.is_some() || options.resume.is_some() {
        Some(checkpoint::model_identity(worker, base_text, &config.render_parameters()?)?)
    } else { None };
    if let Some(path) = &options.resume {
        let cap = 4096 + 8 * MAX_PRODUCT_SAMPLES as u64;
        let file = std::fs::File::open(path).map_err(|e| recovery_error(format!("{}: {e}", path.display())))?;
        let mut bytes = Vec::new();
        file.take(cap + 1).read_to_end(&mut bytes).map_err(|e| recovery_error(e.to_string()))?;
        if bytes.len() as u64 > cap { return Err(recovery_error("parallel QMC checkpoint exceeds its observation cap")); }
        run = QmcBatchExecution::restore(&plan, layout, identity.expect("resume identity"), &bytes)
            .map_err(|e| recovery_error(format!("{}: {e}", path.display())))?;
    }
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    let save = |run: &QmcBatchExecution| -> Result<()> {
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
            super::drain_wave(jobs, slots, &|job| super::evaluate_job(worker, base, config, job, deadline));
        }).map_err(bad)?;
        if report.status == UqStatus::Refused {
            let refusal = model_failure(report.rejection_reason.clone()
                .unwrap_or_else(|| "parallel QMC cooling evaluation refused".into()));
            if let Some(output) = &output {
                output.invalidate(&refusal.to_string()).map_err(|mut error| {
                    error.message.push_str(&format!("; original failure: {refusal}")); error
                })?;
            }
            return Err(refusal);
        }
        save(&run)?;
        if report.status == UqStatus::Cancelled { break; }
    }
    if let Some(result) = &report.result {
        // This is the SERIAL QMC renderer and between-replicate estimator.
        // Terminal restores launch no physics, including under a tiny deadline.
        return Ok(ExecutionOutput {
            stdout: super::super::qmc::render(config, base, result)?, exit_code: exit::SUCCESS,
        });
    }
    let Some(path) = &options.checkpoint else {
        return Err(budget(format!(
            "parallel QMC stopped after {} accepted evaluations; no shortened-net estimate; use --checkpoint to retain completed points",
            report.evaluations_accepted,
        )));
    };
    let termination = if report.status == UqStatus::Cancelled { "wall-time-budget" } else { "sample-chunk" };
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.qmc.parallel.progress.v1\",\"status\":\"budget-truncated\",\"termination\":{},\"qoi\":{},\"workers\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"samples_evaluated_this_run\":{},\"replicates\":{},\"samples_per_replicate\":{},\"contiguous_observations\":{},\"retained_out_of_order\":{},\"next_unpaid_ordinal\":{},\"next_replicate_ordinal\":{},\"next_point_ordinal\":{},\"checkpoint\":{},\"no_claim\":\"retained original QMC points, including successful siblings behind interrupted points; no selected-worker, shortened-net, iid-error or sequential-confidence estimate; resume the original model and layout, optionally changing worker count\"}}\n",
        quote(termination), config.qoi.render(objective_kind(base)?), workers,
        report.samples_planned, report.evaluations_attempted, report.evaluations_attempted - initial,
        layout.replicates, layout.samples_per_replicate,
        report.contiguous_observations, report.retained_out_of_order, report.contiguous_observations,
        report.contiguous_observations / layout.samples_per_replicate,
        report.contiguous_observations % layout.samples_per_replicate, quote(&path.to_string_lossy()),
    );
    Ok(ExecutionOutput { stdout, exit_code: exit::BUDGET })
}
