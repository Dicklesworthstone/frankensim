//! Weighted rare-event probabilities from the original cooling producer.
//! Proposal observations are not samples from the physical target distribution.
mod pilot;
mod recovery;
mod decision;

use super::*;
use fs_uq::{ImportanceConfig, ImportanceExecution, ImportanceReport};

#[derive(Debug, Clone)]
pub(super) enum Settings {
    Fixed(ImportanceConfig),
    Pilot(pilot::Settings),
}

impl Settings {
    pub(super) fn parse_pilot(value: &OsString) -> Result<Self> {
        pilot::Settings::parse(value).map(Self::Pilot)
    }

    pub(super) fn parse(value: &OsString) -> Result<Self> {
        let text = value.to_str().filter(|s| s.len() <= 16_384)
            .ok_or_else(|| bad("--importance requires bounded UTF-8 EPSILON,SHIFT..."))?;
        let parts: Vec<_> = text.split(',').take(258).collect();
        if !(2..=257).contains(&parts.len()) {
            return Err(bad("--importance requires EPSILON followed by 1..=256 germ shifts"));
        }
        let fraction = finite(parts[0], "defensive fraction")?;
        if fraction <= 0.0 || fraction > 1.0 || !fraction.recip().is_finite() {
            return Err(bad("importance defensive fraction must be in (0,1] with finite reciprocal"));
        }
        let germ_shift = parts[1..].iter().map(|s| finite(s, "germ shift"))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self::Fixed(ImportanceConfig { defensive_fraction: fraction, germ_shift }))
    }
}

fn finite(text: &str, name: &str) -> Result<f64> {
    text.parse::<f64>().ok().filter(|v| v.is_finite())
        .ok_or_else(|| bad(format!("importance {name} must be a finite number")))
}

pub(super) fn validate_options(options: &Options) -> Result<()> {
    if options.importance.is_some() && (options.design.is_some()
        || options.qmc_replicates.is_some() || options.sobol_sensitivity || options.shapley_sensitivity
        || options.adjoint_mean_control || options.polynomial.is_some())
    {
        return Err(bad("importance cannot combine other sampling modes, mean controls or candidate selection"));
    }
    Ok(())
}

pub(super) fn execute(
    worker: &Path, base_text: &str, base: &J, config: &Config, options: &Options, settings: &Settings,
) -> Result<ExecutionOutput> {
    let proposal = match settings {
        Settings::Fixed(proposal) => proposal,
        Settings::Pilot(settings) => return pilot::execute(worker, base_text, base, config, options, *settings),
    };
    if options.checkpoint.is_some() || options.resume.is_some() || options.compliance.is_some() {
        return recovery::execute(worker, base_text, base, config, options, proposal);
    }
    let deadline = Instant::now() + Duration::from_secs_f64(config.wall_seconds);
    let mut execution = ImportanceExecution::new(&config.plan(), proposal.clone()).map_err(bad)?;
    let report = execution.advance_interruptible(config.samples, || Instant::now() >= deadline,
        |values| evaluate(worker, base, config, values, deadline));
    complete(report.status, report.rejection_reason.as_deref(), report.evaluations_attempted, deadline)?;
    let stdout = render(base, config, proposal, &report)?;
    complete(report.status, None, report.evaluations_attempted, deadline)?;
    Ok(ExecutionOutput { stdout, exit_code: exit::SUCCESS })
}

fn render(base: &J, config: &Config, proposal: &ImportanceConfig, report: &ImportanceReport) -> Result<String> {
    Ok(format!(
        "{{\"schema\":\"frankensim.cooling-network-uq.importance.v1\",\"authority\":\"estimated-model-importance-sampling\",\"method\":\"defensive-gaussian-importance\",\"status\":\"complete\",\"qoi\":{},\"seed\":{},\"proposal\":{}, {},\"correlation\":{},\"parameters\":[{}],\"no_claim\":{}}}\n",
        config.qoi.render(objective_kind(base)?), quote(&config.seed.to_string()), proposal_json(proposal)?,
        estimate_fields(report)?, quote(config.correlation_label), config.render_parameters()?,
        quote(&no_claim(config)),
    ))
}

fn evaluate(worker: &Path, base: &J, config: &Config, values: &[f64], deadline: Instant) -> Result<Option<f64>> {
    let request = config.sample_request(base, values)?;
    let result = match config.qoi {
        Qoi::Steady => evaluate_sample(worker, &request, deadline),
        qoi => evaluate_sample_for(worker, &request, deadline, qoi),
    };
    match result {
        Ok(value) => Ok(Some(value)),
        Err(EvaluationError::Budget) => Ok(None),
        Err(EvaluationError::Child(message)) => Err(model_failure(message)),
    }
}

fn complete(status: UqStatus, reason: Option<&str>, calls: usize, deadline: Instant) -> Result<()> {
    // A real physical refusal takes precedence over a simultaneous deadline.
    if status == UqStatus::Refused {
        return Err(model_failure(reason.unwrap_or("importance experiment refused")));
    }
    if status != UqStatus::Complete || Instant::now() >= deadline {
        return Err(budget(format!("importance experiment stopped after {calls} completed or rejected physical calls; no partial probability published")));
    }
    Ok(())
}

fn proposal_json(proposal: &ImportanceConfig) -> Result<String> {
    let shifts = proposal.germ_shift.iter().map(|&x| number_json(x))
        .collect::<Result<Vec<_>>>()?.join(",");
    Ok(format!("{{\"defensive_fraction\":{},\"independent_normal_germ_shift\":[{}]}}",
        number_json(proposal.defensive_fraction)?, shifts))
}

fn estimate_fields(report: &ImportanceReport) -> Result<String> {
    if report.status != UqStatus::Complete {
        return Err(model_failure("only a complete importance design may publish a probability"));
    }
    let estimate = report.estimate.as_ref().ok_or_else(|| model_failure(
        report.unavailable_reason.unwrap_or("importance probability unavailable")))?;
    Ok(format!(
        "\"samples_planned\":{},\"samples_evaluated\":{},\"temperature_limit_k\":{},\"failure_event\":\"qoi > temperature_limit_k\",\"estimated_failure_probability\":{},\"failure_standard_error\":{},\"failure_effective_samples\":{},\"mean_likelihood_weight\":{},\"proposal_failure_count\":{}",
        report.evaluations_planned, report.evaluations_attempted, number_json(report.threshold)?,
        number_json(estimate.failure_probability)?, optional_number(estimate.standard_error)?,
        optional_number(estimate.failure_effective_samples)?, number_json(estimate.mean_likelihood_weight)?,
        report.proposal_failures,
    ))
}

fn no_claim(config: &Config) -> String {
    format!("{} Target/proposal weighted strict-failure probability, not the proposal failure frequency or a self-normalized estimate. Finite-sample estimates are not clipped to [0,1]. Standard error and effective samples are descriptive, not confidence bounds or valid stopping rules. Zero observed failures or zero empirical error do not prove zero risk. The fixed Gaussian proposal can miss failure regions; no rare-event discovery, variance reduction or physical qualification is guaranteed. Finite sampling-grid, likelihood-rounding and physical-model errors are not bounded. Proposal responses do not supply a target-temperature distribution.", config.qoi.no_claim())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<OsString> { values.iter().map(OsString::from).collect() }

    #[test]
    fn explicit_proposal_is_bounded_and_keeps_signed_germ_coordinates() {
        let settings = Settings::parse(&OsString::from("0.2,3,-2,0")).unwrap();
        let Settings::Fixed(proposal) = settings else { panic!("expected fixed proposal"); };
        assert_eq!(proposal.defensive_fraction, 0.2);
        assert_eq!(proposal.germ_shift, [3.0, -2.0, 0.0]);
        for value in ["", "0.2", "0.2,", "0,1", "1.1,1", "1e-320,1", "NaN,1", "0.2,inf", "0.2,x"] {
            assert!(Settings::parse(&OsString::from(value)).is_err(), "accepted {value}");
        }
        assert!(Settings::parse(&OsString::from(format!("0.2,{}", vec!["0"; 257].join(",")))).is_err());
    }

    #[test]
    fn fixed_and_pilot_recovery_select_their_own_codecs() {
        for extra in [vec!["--checkpoint", "new.bin", "--max-new-samples", "0"],
            vec!["--resume", "old.bin"], vec!["--resume", "old.bin", "--checkpoint", "new.bin"]]
        {
            let mut fixed = args(&["--importance", "0.2,3"]); fixed.extend(args(&extra));
            assert!(Options::parse(&fixed).is_ok());
            let mut pilot = args(&["--importance-pilot", "1,2,0.2,0.5,1,8"]); pilot.extend(args(&extra));
            assert!(Options::parse(&pilot).is_ok());
        }
    }

    #[test]
    fn conflicting_modes_are_rejected_in_either_option_order() {
        assert!(Options::parse(&args(&["--importance", "0.2,3"])).unwrap().importance.is_some());
        for other in [vec!["--importance", "0.5,0"], vec!["--importance-pilot", "1,2,0.2,0.5,1,8"], vec!["--pce", "1,6"],
            vec!["--pce-mean", "1,6,8"], vec!["--sensitivity", "sobol"],
            vec!["--sensitivity", "shapley"], vec!["--qmc-replicates", "2"],
            vec!["--mean-control", "adjoint"], vec!["--max-new-samples", "0"]]
        {
            let fixed = args(&["--importance", "0.2,3"]);
            for reverse in [false, true] {
                let tail = args(&other);
                let joined: Vec<_> = if reverse { tail.into_iter().chain(fixed.clone()).collect() }
                    else { fixed.clone().into_iter().chain(tail).collect() };
                assert!(Options::parse(&joined).is_err(), "accepted {joined:?}");
            }
        }
    }
}
