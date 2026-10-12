//! Native weighted DWR policy and report projection; numerical work stays in fs-topols.
use super::*;
use fs_topols::goal_resolution::DwrMeasurement;
use fs_topols::robust_descent::MultiLoadProjectedState;
use fs_topols::robust_resolution::goal::{MultiLoadDwrDecision, MultiLoadDwrMeasurement,
    MultiLoadDwrPolicy, MultiLoadDwrProgress, compare_weighted_dwr_gains};

#[path = "goal/read.rs"]
mod recovery;
pub(super) use recovery::{read, validate_terminal};

pub(super) fn parse(fields: &[Node]) -> Result<Option<MultiLoadDwrPolicy>> {
    let Some(pair) = fields.windows(2).find(|pair|
        matches!(&pair[0].kind, NodeKind::Keyword(key) if key == "goal-check"))
    else { return Ok(None) };
    let items = list(&pair[1], "load-family goal-check")?;
    if items.len() != 3 || !matches!(&items[0].kind, NodeKind::Symbol(s) if s == "elasticity-dwr")
        || !matches!(&items[1].kind, NodeKind::Keyword(s) if s == "max-error-fraction") {
        return Err(malformed("load-family goal-check requires (elasticity-dwr :max-error-fraction FRACTION); work uses max-solves"));
    }
    let policy = MultiLoadDwrPolicy { max_error_fraction: node_number(&items[2], "max-error-fraction")? };
    policy.validate(1, 1).map_err(|error| malformed(&error.to_string()))?;
    Ok(Some(policy))
}
pub(super) fn canonical(policy: MultiLoadDwrPolicy, out: &mut String) {
    let _ = writeln!(out, "      :goal-check (elasticity-dwr :max-error-fraction {})",
        canonical_float(policy.max_error_fraction));
}
pub(super) fn validate(spec: &ElasticitySpec, controls: &Controls, family: &LoadFamily) -> Result<()> {
    if let Some(policy) = family.goal {
        policy.validate(settings(spec, spec.steps).level, family.additional.len() + 1)
            .map_err(|error| malformed(&error.to_string()))?;
        if family.aggregate != RobustAggregate::WeightedSum || family.restoration_reduction.is_some()
            || controls.resolution.is_some() || spec.memory_bytes < 256 * 1024 * 1024 {
            return Err(malformed("weighted DWR needs weighted-sum, at least 256 MiB memory, and no mesh-check or stress restoration"));
        }
    }
    Ok(())
}

fn measurement_json(value: &MultiLoadDwrMeasurement) -> String {
    format!("[{}]", value.cases.iter().map(|(load, m)| format!(concat!(
        "{{\"load\":{},\"snapshot\":\"{:#018x}\",\"level\":{},\"area_m2\":{:.17e},",
        "\"coarse_compliance_j\":{:.17e},\"enriched_compliance_j\":{:.17e},\"eta_signed_j\":{:.17e},",
        "\"coarse_relative_residual\":{:.17e},\"enriched_relative_residual\":{:.17e}}}"),
        case_json(load), m.snapshot, m.level, m.area, m.coarse_compliance, m.enriched_compliance,
        m.eta_signed, m.coarse_relative_residual, m.enriched_relative_residual)).collect::<Vec<_>>().join(","))
}
fn decision_json(value: MultiLoadDwrDecision) -> String {
    let a = value.aggregate;
    format!(concat!("{{\"coarse_gain_j\":{:.17e},\"enriched_gain_j\":{:.17e},",
        "\"signed_correction_change_j\":{:.17e},\"observed_gain_change_j\":{:.17e},",
        "\"weighted_correction_variation_j\":{:.17e},\"weighted_gain_variation_j\":{:.17e},",
        "\"allowance_j\":{:.17e},\"refusal\":{}}}"),
        a.coarse_gain, a.enriched_gain, a.signed_correction_change, a.observed_gain_change,
        value.weighted_correction_variation, value.weighted_gain_variation, a.allowance,
        value.refusal.map_or_else(|| "null".into(), quoted))
}

pub(super) struct LastCheck {
    outcome: &'static str,
    baseline: MultiLoadDwrMeasurement,
    candidate: Option<MultiLoadDwrMeasurement>,
    decision: Option<MultiLoadDwrDecision>,
    checked: usize,
}
impl LastCheck {
    pub(super) fn capture(value: MultiLoadDwrProgress) -> Result<(MultiLoadProjectedProgress, Option<Self>)> {
        let (baseline, candidates, progress) = match value {
            MultiLoadDwrProgress::IterationLimit => return Ok((MultiLoadProjectedProgress::IterationLimit, None)),
            MultiLoadDwrProgress::SolveBudget => return Ok((MultiLoadProjectedProgress::SolveBudget(Vec::new()), None)),
            MultiLoadDwrProgress::Searched { baseline, candidates, progress } => (baseline, candidates, progress),
        };
        let (outcome, candidate, decision) = match &progress {
            MultiLoadProjectedProgress::Accepted(step) => {
                let check = candidates.last().filter(|check| check.refusal.is_none()
                    && step.attempts.last().is_some_and(|attempt| attempt.index == check.index))
                    .ok_or_else(|| malformed("accepted weighted DWR update has no matching candidate assessment"))?;
                let measurement = check.measurement.clone()
                    .ok_or_else(|| malformed("missing accepted DWR measurement"))?;
                let decision = check.decision.filter(|d| d.refusal.is_none())
                    .ok_or_else(|| malformed("accepted weighted DWR update failed its gain test"))?;
                ("accepted", Some(measurement), Some(decision))
            }
            MultiLoadProjectedProgress::NoDescent(_) => ("no-descent", None, None),
            MultiLoadProjectedProgress::SolveBudget(_) => ("solve-budget", None, None),
            MultiLoadProjectedProgress::IterationLimit => return Err(malformed("DWR search exceeded update count")),
        };
        let check = Self { outcome, baseline, candidate, decision, checked: candidates.len() };
        Ok((progress, Some(check)))
    }
    fn json(&self) -> String {
        format!("{{\"outcome\":{},\"baseline\":{},\"accepted_candidate\":{},\"decision\":{},\"candidates_checked\":{}}}",
            quoted(self.outcome), measurement_json(&self.baseline),
            self.candidate.as_ref().map_or_else(|| "null".into(), measurement_json),
            self.decision.map_or_else(|| "null".into(), decision_json), self.checked)
    }
}

pub(super) fn field(policy: Option<MultiLoadDwrPolicy>, check: Option<&LastCheck>) -> String {
    policy.map_or_else(String::new, |policy| format!(concat!(",\"goal_resolution\":{{",
        "\"schema\":\"weighted-dwr-gain-v1\",\"authority\":\"Estimated\",",
        "\"solve_accounting\":\"lifetime; two reserved slots per DWR case entered\",",
        "\"max_error_fraction\":{:.17e},\"last_check\":{}}}"), policy.max_error_fraction,
        check.map_or_else(|| "null".into(), LastCheck::json)))
}
pub(super) fn html(policy: Option<MultiLoadDwrPolicy>, check: Option<&LastCheck>) -> String {
    let Some(policy) = policy else { return String::new() };
    let mut out = format!("<h2>Independent-load DWR gain check</h2><p>Estimated screening at {:.6e} of the proposed weighted compliance gain. Each operating condition is solved independently; signed error changes cannot cancel between cases. Stress remains a separate unweighted sampled constraint. No continuum error bound or optimum is certified. Study solve accounting reserves two slots per DWR case entered, including failed cases.</p>", policy.max_error_fraction);
    if let Some(check) = check {
        let _ = write!(out, "<p>Last complete search: {}. Candidates checked: {}.</p><table><tr><th>Case</th><th>Weight</th><th>Baseline coarse J</th><th>Baseline enriched J</th><th>Candidate coarse J</th><th>Candidate enriched J</th></tr>", check.outcome, check.checked);
        for (i, (load, m)) in check.baseline.cases.iter().enumerate() {
            let candidate = check.candidate.as_ref().map(|c| {
                let c = c.cases[i].1;
                format!("<td>{:.8e}</td><td>{:.8e}</td>", c.coarse_compliance, c.enriched_compliance)
            }).unwrap_or_else(|| "<td>not accepted</td><td>not accepted</td>".into());
            let _ = write!(out, "<tr><td>{i}</td><td>{:.8e}</td><td>{:.8e}</td><td>{:.8e}</td>{candidate}</tr>",
                load.weight(), m.coarse_compliance, m.enriched_compliance);
        }
        out.push_str("</table>");
    } else { out.push_str("<p>No complete DWR search has been retained.</p>"); }
    out
}

#[cfg(test)]
#[path = "goal_tests.rs"]
mod tests;
