//! Concurrent original cooling children with ordered, recoverable UQ observations.
mod qmc;
mod sobol;
mod shapley;

use super::*;
use fs_uq::{UqBatchExecution, UqBatchInput};
use std::io::Read as _;

pub(super) fn validate_options(options: &Options) -> Result<()> {
    if options.workers.is_some() && (options.compliance.is_some() || options.design.is_some()
        || options.adjoint_mean_control || options.polynomial.is_some() || options.importance.is_some())
    {
        return Err(bad("--workers supports fixed-count Monte Carlo, replicated QMC and Sobol/Shapley sensitivity; other sampling, control and decision modes are incompatible"));
    }
    Ok(())
}

fn recovery_error(message: impl Into<String>) -> Failure {
    Failure { code: "cooling-network-uq-checkpoint", message: message.into() }
}

pub(super) fn execute(
    worker: &Path, base_text: &str, base: &J, config: &Config, options: &Options, workers: usize,
) -> Result<ExecutionOutput> {
    if options.shapley_sensitivity {
        return shapley::execute(worker, base_text, base, config, options, workers);
    }
    if options.sobol_sensitivity {
        return sobol::execute(worker, base_text, base, config, options, workers);
    }
    if let Some(replicates) = options.qmc_replicates {
        return qmc::execute(worker, base_text, base, config, options, workers, replicates);
    }
    let plan = config.plan();
    let mut run = UqBatchExecution::new(&plan).map_err(bad)?;
    let identity = if options.checkpoint.is_some() || options.resume.is_some() {
        Some(checkpoint::model_identity(worker, base_text, &config.render_parameters()?)?)
    } else { None };
    if let Some(path) = &options.resume {
        let cap = 4096 + 8 * MAX_PRODUCT_SAMPLES as u64;
        let file = std::fs::File::open(path).map_err(|e| recovery_error(format!("{}: {e}", path.display())))?;
        let mut bytes = Vec::new();
        file.take(cap + 1).read_to_end(&mut bytes).map_err(|e| recovery_error(e.to_string()))?;
        if bytes.len() as u64 > cap { return Err(recovery_error("parallel checkpoint exceeds its observation cap")); }
        run = UqBatchExecution::restore(&plan, identity.expect("resume identity"), &bytes)
            .map_err(|e| recovery_error(format!("{}: {e}", path.display())))?;
    }
    // Admission and exact model/law recovery precede output reservation.
    let output = options.checkpoint.as_deref().map(checkpoint::Output::reserve).transpose()?;
    let save = |run: &UqBatchExecution| -> Result<()> {
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
        // All workers and child watchdogs have drained here. Genuine failures
        // win over timeout; never save a resumable prefix past a failed sample.
        if report.status == UqStatus::Refused {
            let refusal = model_failure(report.rejection_reason.clone()
                .unwrap_or_else(|| "parallel cooling evaluation refused".into()));
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
        // Keep exactly the original complete MC result schema and reductions.
        return Ok(ExecutionOutput { stdout: render_result(config, base, result)?, exit_code: exit::SUCCESS });
    }
    let Some(path) = &options.checkpoint else {
        return Err(budget(format!("parallel cooling stopped after {} completed evaluations; no partial distribution; use --checkpoint to retain successful siblings", report.evaluations_accepted)));
    };
    let termination = if report.status == UqStatus::Cancelled { "wall-time-budget" } else { "sample-chunk" };
    let stdout = format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.parallel.progress.v1\",\"status\":\"budget-truncated\",\"termination\":{},\"qoi\":{},\"workers\":{},\"samples_planned\":{},\"samples_evaluated\":{},\"samples_evaluated_this_run\":{},\"contiguous_observations\":{},\"retained_out_of_order\":{},\"next_unpaid_ordinal\":{},\"checkpoint\":{},\"no_claim\":\"retained original-ordinal physical responses, including later successful siblings; no selected-worker or shortened-design statistics; resume with the same model, law and budget, optionally changing worker count\"}}\n",
        quote(termination), config.qoi.render(objective_kind(base)?), workers,
        report.samples_planned, report.evaluations_attempted, report.evaluations_attempted - initial,
        report.contiguous_observations, report.retained_out_of_order, report.contiguous_observations,
        quote(&path.to_string_lossy()),
    );
    Ok(ExecutionOutput { stdout, exit_code: exit::BUDGET })
}

fn evaluate_job(
    worker: &Path, base: &J, config: &Config, job: &UqBatchInput, deadline: Instant,
) -> std::result::Result<Option<f64>, String> {
    if Instant::now() >= deadline { return Ok(None); }
    let request = config.sample_request(base, &job.parameters).map_err(|e| e.to_string())?;
    let result = match config.qoi {
        Qoi::Steady => evaluate_sample(worker, &request, deadline),
        qoi => evaluate_sample_for(worker, &request, deadline, qoi),
    };
    match result {
        Ok(value) => Ok(Some(value)),
        Err(EvaluationError::Budget) => Ok(None),
        Err(EvaluationError::Child(message)) => Err(message),
    }
}

/// Scoped producer threads contain worker panics and join every launched child
/// driver. The child owner itself kills/reaps timed-out processes and I/O workers.
fn drain_wave(
    jobs: &[UqBatchInput], slots: &mut [Option<std::result::Result<f64, String>>],
    evaluate: &(impl Fn(&UqBatchInput) -> std::result::Result<Option<f64>, String> + Sync),
) {
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(jobs.len());
        for (i, job) in jobs.iter().enumerate() {
            match std::thread::Builder::new().spawn_scoped(scope, move || evaluate(job)) {
                Ok(handle) => handles.push((i, handle)),
                Err(error) => {
                    slots[i] = Some(Err(format!("cannot start cooling worker: {error}")));
                    break;
                }
            }
        }
        for (i, handle) in handles {
            slots[i] = match handle.join() {
                Ok(Ok(value)) => value.map(Ok),
                Ok(Err(error)) => Some(Err(error)),
                Err(_) => Some(Err("cooling worker panicked; experiment refused".into())),
            };
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(a: &[&str]) -> Vec<OsString> { a.iter().map(OsString::from).collect() }

    #[test]
    fn worker_options_are_bounded_and_never_silently_change_another_method() {
        for n in ["1", "2", "64"] {
            assert!(Options::parse(&args(&["--workers", n, "--checkpoint", "new", "--max-new-samples", "0"])).is_ok());
        }
        for n in ["0", "65", "-1", "1.5", "NaN"] {
            assert!(Options::parse(&args(&["--workers", n])).is_err());
        }
        for extra in [vec!["--workers", "2"],
            vec!["--pce", "0,2"], vec!["--pce-mean", "0,2,4"],
            vec!["--importance", "0.2,1"], vec!["--importance-pilot", "1,2,0.2,0.5,1,6"],
            vec!["--mean-control", "adjoint"],
            vec!["--compliance-probability", "0.9", "--confidence-alpha", "0.05", "--min-decision-samples", "2"]]
        {
            let own = args(&["--workers", "2"]); let other = args(&extra);
            for input in [own.iter().chain(&other).cloned().collect::<Vec<_>>(), other.iter().chain(&own).cloned().collect()] {
                assert!(Options::parse(&input).is_err());
            }
        }
    }

    #[test]
    fn qmc_workers_preserve_explicit_layout_and_recovery_options() {
        for input in [
            vec!["--workers", "3", "--qmc-replicates", "2"],
            vec!["--qmc-replicates", "2", "--workers", "3"],
            vec!["--workers", "1", "--qmc-replicates", "2", "--resume", "old", "--checkpoint", "new", "--max-new-samples", "0"],
        ] {
            let parsed = Options::parse(&args(&input)).unwrap();
            assert_eq!(parsed.qmc_replicates, Some(2));
            assert!(parsed.workers.is_some());
        }
    }

    #[test]
    fn sobol_workers_retain_pairing_recovery_and_conflict_admission() {
        for input in [
            vec!["--workers", "3", "--sensitivity", "sobol"],
            vec!["--sensitivity", "sobol", "--workers", "3"],
            vec!["--workers", "1", "--sensitivity", "sobol", "--resume", "old", "--checkpoint", "new", "--max-new-samples", "0"],
        ] {
            let parsed = Options::parse(&args(&input)).unwrap();
            assert!(parsed.sobol_sensitivity);
            assert!(parsed.workers.is_some());
        }
        for extra in [vec!["--qmc-replicates", "2"], vec!["--sensitivity", "shapley"],
            vec!["--mean-control", "adjoint"], vec!["--pce", "1,6"],
            vec!["--compliance-probability", "0.9", "--confidence-alpha", "0.05", "--min-decision-samples", "2"]]
        {
            let mut input = args(&["--workers", "3", "--sensitivity", "sobol"]);
            input.extend(args(&extra));
            assert!(Options::parse(&input).is_err());
        }
    }

    #[test]
    fn shapley_workers_preserve_correlated_design_and_recovery_options() {
        for input in [
            vec!["--workers", "3", "--sensitivity", "shapley"],
            vec!["--sensitivity", "shapley", "--workers", "3"],
            vec!["--workers", "1", "--sensitivity", "shapley", "--resume", "old", "--checkpoint", "new", "--max-new-samples", "0"],
        ] {
            let parsed = Options::parse(&args(&input)).unwrap();
            assert!(parsed.shapley_sensitivity && !parsed.sobol_sensitivity);
            assert!(parsed.workers.is_some());
        }
        for extra in [vec!["--sensitivity", "sobol"], vec!["--qmc-replicates", "2"],
            vec!["--importance", "0.2,1"], vec!["--mean-control", "adjoint"],
            vec!["--compliance-probability", "0.9", "--confidence-alpha", "0.05", "--min-decision-samples", "2"]]
        {
            let mut input = args(&["--workers", "3", "--sensitivity", "shapley"]);
            input.extend(args(&extra));
            assert!(Options::parse(&input).is_err());
        }
    }

    #[test]
    fn waves_overlap_drain_and_keep_panics_separate_from_unfinished_calls() {
        let jobs: Vec<_> = (0..2).map(|ordinal| UqBatchInput { ordinal, parameters: vec![ordinal as f64] }).collect();
        let mut slots = vec![None, None];
        let barrier = std::sync::Barrier::new(2);
        drain_wave(&jobs, &mut slots, &|job| {
            barrier.wait(); Ok(Some(job.parameters[0]))
        });
        assert_eq!(slots, [Some(Ok(0.0)), Some(Ok(1.0))]);
        drain_wave(&jobs, &mut slots, &|job| {
            if job.ordinal == 0 { panic!("intentional worker failure"); }
            Ok(None)
        });
        assert!(slots[0].as_ref().unwrap().is_err());
        assert!(slots[1].is_none());
    }
}
