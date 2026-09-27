//! Native QMC projections use complete independent scrambles as the sampling
//! units. An unfinished net remains paid work, never an iid observation set.

use fs_blake3::ContentHash;
use fs_package::Claim;
use fs_uq::{QmcEstimate, QmcReport};

use super::optional;

pub(super) const SCOPE: &str = "Estimated randomized Sobol quadrature of the declared native numerical cooling model under independent uniform input laws. Independent Owen scrambles, not dependent points within a net, are the units for descriptive standard errors. Only complete equal-sized nets enter mean and compliance estimates; unfinished net points remain retained paid work. These standard errors are not confidence intervals or optional-stopping bounds and do not bound finite-grid quadrature bias, numerical error or physical-model error. Zero between-replicate variation does not prove exactness. Child engineering uncertainty budgets and verdicts remain unchanged. Refused executions publish no estimates; samples are never replaced, clipped or skipped.";

pub(super) fn json(report: &QmcReport) -> String {
    let means = report.replicate_means.iter().map(ToString::to_string)
        .collect::<Vec<_>>().join(",");
    let estimate = report.estimate.as_ref().map_or_else(|| "null".into(), |estimate| {
        format!("{{\"mean_k\":{},\"sampling_standard_error_k\":{}}}",
            estimate.mean, optional(estimate.standard_error))
    });
    let compliance = report.compliance.as_ref().map_or_else(|| "null".into(), |estimate| {
        format!("{{\"probability_of_compliance\":{},\"sampling_standard_error\":{}}}",
            estimate.mean, optional(estimate.standard_error))
    });
    format!(
        "{{\"sampler\":\"owen-scrambled-sobol\",\"replicates_planned\":{},\"samples_per_replicate\":{},\"completed_replicates\":{},\"partial_replicate_samples\":{},\"samples_in_estimate\":{},\"replicate_means_k\":[{means}],\"estimate\":{estimate},\"compliance_estimate\":{compliance}}}",
        report.config.replicates,
        report.config.samples_per_replicate,
        report.completed_replicates,
        report.samples_accepted % report.config.samples_per_replicate,
        report.replicate_means.len() * report.config.samples_per_replicate,
    )
}

fn estimate_text(estimate: Option<&QmcEstimate>, unit: &str) -> String {
    estimate.map_or_else(|| "unavailable".into(), |estimate| {
        format!("{} {unit}; between-replicate standard error {} {unit}",
            estimate.mean,
            estimate.standard_error.map_or_else(|| "unavailable".into(), |value| value.to_string()))
    })
}

pub(super) fn html(report: &QmcReport) -> String {
    format!(
        "<h2>Replicated randomized Sobol quadrature</h2><p>{}/{} complete scrambles, {} points per scramble; {} retained points in an unfinished net are excluded from estimates.</p><p>Mean temperature: {}. Probability of the numerical pass event: {}.</p>",
        report.completed_replicates, report.config.replicates,
        report.config.samples_per_replicate,
        report.samples_accepted % report.config.samples_per_replicate,
        estimate_text(report.estimate.as_ref(), "K"),
        estimate_text(report.compliance.as_ref(), ""),
    )
}

pub(super) fn claim(report: &QmcReport, result: ContentHash) -> Option<Claim> {
    let estimate = report.estimate.as_ref()?;
    let error = estimate.standard_error?;
    Some(Claim::estimated(
        "cooling.uncertainty.sample-mean",
        format!(
            "{} K across {} complete scrambles of {} native solves; descriptive between-replicate standard error {} K. Result {}. {SCOPE}",
            estimate.mean, report.completed_replicates, report.config.samples_per_replicate,
            error, result.to_hex(),
        ),
        "fixed-count-native-randomized-sobol-between-replicate-standard-error",
        error,
    ))
}
