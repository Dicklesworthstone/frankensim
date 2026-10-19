//! Native Latin-hypercube projections: complete independent replicates, never
//! an iid-point Monte Carlo report or a Sobol-scramble claim.

use fs_blake3::ContentHash;
use fs_package::Claim;
use fs_project::uncertainty::LhsLayout;
use fs_uq::{LhsEstimate, LhsReport};

use super::{Model, input_law, optional};

const SCOPE: &str = "Independent randomized Latin hypercubes, not dependent points within one hypercube, are the units for descriptive standard errors. Only complete equal-sized hypercubes enter mean and compliance estimates; unfinished points remain retained paid work. These standard errors are not confidence intervals or optional-stopping bounds and do not bound finite-grid sampling bias, numerical error or physical-model error. Zero between-replicate variation does not prove exactness. Child engineering uncertainty budgets and verdicts remain unchanged. Refused executions publish no estimates; samples are never replaced, clipped or skipped.";

pub(super) fn scope(model: &Model) -> String {
    format!("Estimated replicated Latin-hypercube sampling of the declared native numerical cooling model under {}. {SCOPE}", input_law(model))
}

pub(super) fn layout(model: &Model) -> LhsLayout {
    model.bound.study().lhs().expect("LHS execution has an admitted native layout")
}

pub(super) fn json(report: &LhsReport, layout: LhsLayout) -> String {
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
        "{{\"sampler\":\"philox-latin-hypercube\",\"replicates_planned\":{},\"samples_per_replicate\":{},\"completed_replicates\":{},\"partial_replicate_samples\":{},\"samples_in_estimate\":{},\"replicate_means_k\":[{means}],\"estimate\":{estimate},\"compliance_estimate\":{compliance}}}",
        layout.replicates,
        layout.samples_per_replicate,
        report.completed_replicates,
        report.samples_accepted % layout.samples_per_replicate,
        report.replicate_means.len() * layout.samples_per_replicate,
    )
}

fn estimate_text(estimate: Option<&LhsEstimate>, unit: &str) -> String {
    estimate.map_or_else(|| "unavailable".into(), |estimate| {
        format!("{} {unit}; between-replicate standard error {} {unit}",
            estimate.mean,
            estimate.standard_error.map_or_else(|| "unavailable".into(), |value| value.to_string()))
    })
}

pub(super) fn html(report: &LhsReport, layout: LhsLayout) -> String {
    format!(
        "<h2>Replicated Latin-hypercube sampling</h2><p>{}/{} complete hypercubes, {} strata per coordinate; {} retained points in an unfinished hypercube are excluded from estimates.</p><p>Mean temperature: {}. Probability of the numerical pass event: {}.</p>",
        report.completed_replicates, layout.replicates,
        layout.samples_per_replicate,
        report.samples_accepted % layout.samples_per_replicate,
        estimate_text(report.estimate.as_ref(), "K"),
        estimate_text(report.compliance.as_ref(), ""),
    )
}

pub(super) fn claim(report: &LhsReport, layout: LhsLayout, result: ContentHash, scope: &str) -> Option<Claim> {
    let estimate = report.estimate.as_ref()?;
    let error = estimate.standard_error?;
    Some(Claim::estimated(
        "cooling.uncertainty.sample-mean",
        format!(
            "{} K across {} complete Latin hypercubes of {} native solves; descriptive between-replicate standard error {} K. Result {}. {scope}",
            estimate.mean, report.completed_replicates, layout.samples_per_replicate,
            error, result.to_hex(),
        ),
        "fixed-count-native-latin-hypercube-between-replicate-standard-error",
        error,
    ))
}
