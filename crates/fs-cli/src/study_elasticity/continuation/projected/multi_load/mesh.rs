//! Native adapter for the existing all-grid/all-load acceptance gate.
//! Reports retain independent cases; aggregate agreement cannot hide a bad case.
use super::*;
use fs_topols::resolution::ResolutionPolicy;
use fs_topols::robust_resolution::{
    MultiLoadMeshProgress, MultiLoadResolutionReport, MultiLoadResolutionRung,
};

pub(super) struct LastCheck {
    pub(super) outcome: &'static str,
    pub(super) baseline: MultiLoadResolutionReport,
    candidate: Option<MultiLoadResolutionReport>,
    checked: usize,
    pub(super) reason: Option<String>,
}

impl LastCheck {
    pub(super) fn capture(progress: MultiLoadMeshProgress)
        -> Result<(Option<MultiLoadProjectedProgress>, Option<Self>)>
    {
        match progress {
            MultiLoadMeshProgress::IterationLimit => Err(malformed("load-family mesh check exceeded the update budget")),
            MultiLoadMeshProgress::SolveBudget => Ok((Some(MultiLoadProjectedProgress::SolveBudget(Vec::new())), None)),
            MultiLoadMeshProgress::UnresolvedBaseline { baseline, reason } => Ok((None, Some(Self {
                outcome: "baseline-unresolved", baseline, candidate: None, checked: 0, reason: Some(reason),
            }))),
            MultiLoadMeshProgress::Searched { baseline, candidates, progress } => {
                let (outcome, candidate) = match &progress {
                    MultiLoadProjectedProgress::Accepted(step) => {
                        let index = step.attempts.last().ok_or_else(|| malformed("accepted load family has no attempt"))?.index;
                        let report = candidates.iter().find(|c| c.index == index && c.refusal.is_none())
                            .and_then(|c| c.report.clone()).ok_or_else(|| malformed("accepted load family lacks its mesh assessment"))?;
                        ("accepted", Some(report))
                    }
                    MultiLoadProjectedProgress::NoDescent(_) => ("no-descent", None),
                    MultiLoadProjectedProgress::SolveBudget(_) => ("solve-budget", None),
                    MultiLoadProjectedProgress::IterationLimit => return Err(malformed("unexpected load-family mesh iteration limit")),
                };
                let checked = candidates.len();
                Ok((Some(progress), Some(Self { outcome, baseline, candidate, checked, reason: None })))
            }
        }
    }

    fn json(&self) -> String {
        format!("{{\"outcome\":{},\"baseline\":{},\"accepted_candidate\":{},\"candidates_checked\":{},\"reason\":{}}}",
            quoted(self.outcome), report_json(&self.baseline),
            self.candidate.as_ref().map_or_else(|| "null".into(), report_json), self.checked,
            self.reason.as_ref().map_or_else(|| "null".into(), |reason| quoted(reason)))
    }
}

fn policy_json(p: ResolutionPolicy) -> String {
    format!("{{\"extra_levels\":{},\"absolute_compliance_tolerance_j\":{:.17e},\"relative_compliance_tolerance\":{:.17e},\"area_tolerance_m2\":{:.17e}}}",
        p.extra_levels, p.absolute_compliance_tolerance, p.relative_compliance_tolerance, p.area_tolerance)
}

fn report_json(report: &MultiLoadResolutionReport) -> String {
    format!("{{\"rungs\":[{}]}}", report.rungs.iter().map(|rung|
        format!("{{\"level\":{},\"cases\":{}}}", rung.level, states_json(&rung.cases)))
        .collect::<Vec<_>>().join(","))
}

pub(super) fn field(policy: Option<ResolutionPolicy>, check: Option<&LastCheck>) -> String {
    policy.map_or_else(String::new, |policy| format!(
        ",\"mesh_resolution\":{{\"authority\":\"Estimated\",\"scope\":\"every independent case on every grid; observed resolution and sampled stress, not a continuum certificate\",\"policy\":{},\"last_check\":{}}}",
        policy_json(policy), check.map_or_else(|| "null".into(), LastCheck::json)))
}

pub(super) fn html(policy: Option<ResolutionPolicy>, check: Option<&LastCheck>) -> String {
    let Some(policy) = policy else { return String::new() };
    let mut out = format!("<h2>Mesh-checked independent loads</h2><p>{} additional grids. Every case, including zero-weight cases, must pass the unchanged sampled stress limit and its own compliance-resolution tolerance. The declared aggregate must improve at every grid. These are observed numerical checks, not continuous stress or error certificates.</p>", policy.extra_levels);
    let Some(check) = check else { out.push_str("<p>No complete mesh assessment retained yet.</p>"); return out };
    let _ = write!(out, "<p>Last complete check: {}; {} candidates assessed.</p><table><tr><th>Level</th><th>Case</th><th>Baseline J</th><th>Baseline sampled Pa</th><th>Accepted J</th><th>Accepted sampled Pa</th></tr>", check.outcome, check.checked);
    for (i, rung) in check.baseline.rungs.iter().enumerate() {
        for (case, state) in rung.cases.iter().enumerate() {
            let candidate = check.candidate.as_ref().map_or_else(|| "<td>not accepted</td><td>not accepted</td>".into(), |report| {
                let state = &report.rungs[i].cases[case];
                format!("<td>{:.8e}</td><td>{:.8e}</td>", state.compliance, state.sampled_max_von_mises)
            });
            let _ = write!(out, "<tr><td>{}</td><td>{case}</td><td>{:.8e}</td><td>{:.8e}</td>{candidate}</tr>",
                rung.level, state.compliance, state.sampled_max_von_mises);
        }
    }
    out.push_str("</table>");
    out
}

fn read_report(value: &JsonValue, p: ResolutionPolicy, policy: &Controls,
    cases: &[RobustLoadCase]) -> Result<MultiLoadResolutionReport>
{
    let rows = value.get("rungs").and_then(JsonValue::as_array)
        .filter(|rows| rows.len() == p.extra_levels as usize + 1)
        .ok_or_else(|| malformed("incomplete load-family mesh report"))?;
    let rungs = rows.iter().map(|row| {
        let states = row.get("cases").and_then(JsonValue::as_array)
            .filter(|states| states.len() == cases.len())
            .ok_or_else(|| malformed("mesh assessment omitted an independent case"))?;
        Ok(MultiLoadResolutionRung {
            level: u32::try_from(integer(row, "level")?).map_err(|_| malformed("invalid family mesh level"))?,
            cases: states.iter().map(read_stress_measurement).collect::<Result<Vec<_>>>()?,
        })
    }).collect::<Result<Vec<_>>>()?;
    let report = MultiLoadResolutionReport { rungs };
    // The numerical owner validates complete finite measurements even when
    // they genuinely fail the requested stress/resolution gates.
    let _ = report.refusal(p, policy.area.target, policy.stress, cases).map_err(|e| malformed(&e.to_string()))?;
    Ok(report)
}

pub(super) fn read(value: &JsonValue, policy: &Controls, cases: &[RobustLoadCase],
    baseline: &[SampledStressEvaluation], accepted: &[Vec<SampledStressEvaluation>]) -> Result<Option<LastCheck>>
{
    let node = value.get("mesh_resolution");
    let Some(p) = policy.resolution else {
        return if node.is_none() { Ok(None) } else { Err(malformed("undeclared load-family mesh evidence")) };
    };
    let node = node.ok_or_else(|| malformed("missing load-family mesh evidence"))?;
    if node.str_field("authority") != Some("Estimated")
        || node.get("policy") != Some(&document(policy_json(p).as_bytes())?) {
        return Err(malformed("retained load-family mesh policy changed"));
    }
    let last = node.get("last_check").ok_or_else(|| malformed("missing family mesh outcome"))?;
    if matches!(last, JsonValue::Null) {
        return if accepted.is_empty() { Ok(None) } else { Err(malformed("accepted load families lack mesh checks")) };
    }
    let outcome = match last.str_field("outcome") {
        Some("accepted") => "accepted", Some("baseline-unresolved") => "baseline-unresolved",
        Some("no-descent") => "no-descent", Some("solve-budget") => "solve-budget",
        _ => return Err(malformed("unknown load-family mesh outcome")),
    };
    let original = baseline;
    let baseline = read_report(last.get("baseline").ok_or_else(|| malformed("missing mesh baseline"))?, p, policy, cases)?;
    let candidate = match last.get("accepted_candidate") {
        Some(JsonValue::Null) => None, Some(value) => Some(read_report(value, p, policy, cases)?),
        None => return Err(malformed("missing family mesh candidate")),
    };
    let checked = integer(last, "candidates_checked")?;
    let reason = match last.get("reason") {
        Some(JsonValue::Null) => None, Some(JsonValue::Str(s)) => Some(s.clone()),
        _ => return Err(malformed("invalid family mesh refusal")),
    };
    let current = accepted.last().map_or(original, Vec::as_slice);
    let previous = if accepted.len() >= 2 { accepted[accepted.len()-2].as_slice() } else { original };
    let expected = if outcome == "accepted" { previous } else { current };
    if checked > policy.search.max_candidates || (outcome == "accepted") != candidate.is_some()
        || baseline.rungs[0].cases.iter().zip(expected).any(|(a, b)| !same(a, b)) {
        return Err(malformed("load-family mesh baseline disagrees with accepted history"));
    }
    let refused = baseline.refusal(p, policy.area.target, policy.stress, cases).map_err(|e| malformed(&e.to_string()))?;
    if outcome == "baseline-unresolved" {
        if checked != 0 || refused.is_none() || refused != reason {
            return Err(malformed("unresolved load family lacks its measured refusal"));
        }
    } else if refused.is_some() || reason.is_some() {
        return Err(malformed("load-family search began with an unresolved baseline"));
    }
    if let Some(candidate) = &candidate {
        let aggregate = policy.family.as_ref().ok_or_else(|| malformed("missing load policy"))?.aggregate;
        if accepted.is_empty() || checked == 0
            || candidate.rungs[0].cases.iter().zip(current).any(|(a, b)| !same(a, b))
            || candidate.comparison_refusal(&baseline, p, policy.area.target, policy.stress, cases,
                aggregate, policy.search.min_relative_improvement).map_err(|e| malformed(&e.to_string()))?.is_some() {
            return Err(malformed("accepted load family violates its finer-grid gates"));
        }
    }
    Ok(Some(LastCheck { outcome, baseline, candidate, checked, reason }))
}

/// Bind read-only terminals and assessment-only retries to the last complete
/// mesh gate, without spending a recovery solve to inspect retained results.
pub(super) fn validate_terminal(history: &History, policy: &Controls, status: &str, level: u32) -> Result<()> {
    if policy.resolution.is_none() { return Ok(()) }
    let expected = match history.optimizer_terminal.unwrap_or(status) {
        "completed" => Some("accepted"),
        "no-feasible-descent" => Some("no-descent"),
        "mesh-unresolved" => Some("baseline-unresolved"),
        _ => None,
    };
    if history.mesh.as_ref().is_some_and(|check| check.baseline.rungs[0].level != level)
        || expected.is_some_and(|expected| history.mesh.as_ref().is_none_or(|check| check.outcome != expected)) {
        return Err(malformed("load-family mesh evidence disagrees with its terminal or declared grid"));
    }
    Ok(())
}

#[cfg(test)]
mod tests;
