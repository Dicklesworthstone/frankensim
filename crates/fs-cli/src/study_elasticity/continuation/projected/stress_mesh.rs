//! Native declaration/report adapter for the original mesh-checked stress owner.
//! No alternate stress sampler, mechanics solver or candidate search.
use super::*;
pub(super) use fs_topols::resolution::{MeshCheckStage, ResolutionPolicy};
use fs_topols::projected_stress::resolution::{MeshCheckedStressProgress, StressResolutionReport, StressResolutionRung};
use fs_topols::ProjectedStressUpdate;

pub(super) fn parse(fields: &[Node]) -> Result<Option<ResolutionPolicy>> {
    let Some(pair) = fields.windows(2).find(|pair|
        matches!(&pair[0].kind, NodeKind::Keyword(k) if k == "mesh-check"))
    else { return Ok(None) };
    let items = list(&pair[1], "mesh-check")?;
    if items.len() != 9 || !matches!(&items[0].kind, NodeKind::Symbol(k) if k == "refinement") {
        return Err(malformed("stress mesh-check requires refinement and four explicit policy fields"));
    }
    for (i, key) in ["extra-levels", "absolute-compliance-tolerance-j", "relative-compliance-tolerance", "area-tolerance-m2"].iter().enumerate() {
        if !matches!(&items[2*i+1].kind, NodeKind::Keyword(k) if k.as_str() == *key) {
            return Err(malformed("stress mesh-check fields must be complete, unique and in canonical order"));
        }
    }
    let real = |i| super::super::super::number(&items[i], "mesh-check tolerance");
    let policy = ResolutionPolicy {
        extra_levels: u32::try_from(integer_node(&items[2], "extra-levels")?)
            .map_err(|_| malformed("mesh-check levels exceed u32"))?,
        absolute_compliance_tolerance: real(4)?, relative_compliance_tolerance: real(6)?, area_tolerance: real(8)?,
    };
    policy.validate(1).map_err(|e| malformed(&e.to_string()))?;
    Ok(Some(policy))
}
pub(super) fn canonical(p: ResolutionPolicy, out: &mut String) {
    let _ = writeln!(out, "    :mesh-check (refinement");
    let _ = writeln!(out, "      :extra-levels {}", p.extra_levels);
    let _ = writeln!(out, "      :absolute-compliance-tolerance-j {}", canonical_float(p.absolute_compliance_tolerance));
    let _ = writeln!(out, "      :relative-compliance-tolerance {}", canonical_float(p.relative_compliance_tolerance));
    let _ = writeln!(out, "      :area-tolerance-m2 {})", canonical_float(p.area_tolerance));
}
fn policy_json(p: ResolutionPolicy) -> String {
    format!("{{\"extra_levels\":{},\"absolute_compliance_tolerance_j\":{:.17e},\"relative_compliance_tolerance\":{:.17e},\"area_tolerance_m2\":{:.17e}}}",
        p.extra_levels, p.absolute_compliance_tolerance, p.relative_compliance_tolerance, p.area_tolerance)
}
fn report_json(report: &StressResolutionReport) -> String {
    format!("{{\"rungs\":[{}]}}", report.rungs.iter().map(|r|
        format!("{{\"level\":{},\"evaluation\":{}}}", r.level, stress_json(&r.evaluation))).collect::<Vec<_>>().join(","))
}

/// Only a COMPLETE assessment is retained. All earlier ones remain in the
/// ordinary predecessor chain; cancelled work never changes this record.
pub(super) struct LastCheck {
    pub(super) outcome: &'static str,
    pub(super) baseline: StressResolutionReport,
    candidate: Option<StressResolutionReport>,
    checked: usize,
    pub(super) reason: Option<String>,
}
impl LastCheck {
    pub(super) fn capture(progress: MeshCheckedStressProgress) -> Result<(Option<ProjectedStressUpdate>, Self)> {
        match progress {
            MeshCheckedStressProgress::IterationLimit => Err(malformed("stress mesh check ran beyond the update budget")),
            MeshCheckedStressProgress::UnresolvedBaseline { baseline, reason } => Ok((None, Self {
                outcome: "baseline-unresolved", baseline, candidate: None, checked: 0, reason: Some(reason),
            })),
            MeshCheckedStressProgress::Searched { baseline, candidates, update } => {
                let (outcome, candidate) = match &update.progress {
                    ProjectedProgress::Accepted(_) => ("accepted", Some(candidates.last()
                        .filter(|c| c.refusal.is_none()).and_then(|c| c.report.clone())
                        .ok_or_else(|| malformed("accepted stress candidate lacks its complete mesh check"))?)),
                    ProjectedProgress::NoDescent(_) => ("no-descent", None),
                    ProjectedProgress::IterationLimit => return Err(malformed("unexpected stress mesh iteration limit")),
                };
                let checked = candidates.len();
                Ok((Some(update), Self { outcome, baseline, candidate, checked, reason: None }))
            }
        }
    }
    fn json(&self) -> String {
        format!("{{\"outcome\":{},\"baseline\":{},\"accepted_candidate\":{},\"candidates_checked\":{},\"reason\":{}}}",
            quoted(self.outcome), report_json(&self.baseline),
            self.candidate.as_ref().map_or_else(|| "null".into(), report_json), self.checked,
            self.reason.as_ref().map_or_else(|| "null".into(), |r| quoted(r)))
    }
}
pub(super) fn field(policy: Option<ResolutionPolicy>, check: Option<&LastCheck>) -> String {
    policy.map_or_else(String::new, |p| format!(concat!(",\"mesh_resolution\":{{",
        "\"schema\":\"observed-stress-mesh-v1\",\"authority\":\"Estimated\",",
        "\"scope\":\"unchanged stress limit on every grid; sampled, not a continuous stress certificate\",",
        "\"policy\":{},\"last_check\":{}}}"), policy_json(p),
        check.map_or_else(|| "null".into(), LastCheck::json)))
}
pub(super) fn html(policy: Option<ResolutionPolicy>, check: Option<&LastCheck>) -> String {
    let Some(p) = policy else { return String::new() };
    let mut out = format!("<h2>Mesh-checked sampled stress</h2><p>{} additional grids; unchanged material, load and stress allowance at every level. Finer fields are not re-projected. This is observed resolution and sampled stress, not a continuous stress certificate.</p>", p.extra_levels);
    let Some(check) = check else { out.push_str("<p>Not yet assessed; no finer-grid pass is claimed.</p>"); return out };
    let _ = write!(out, "<p>Last complete search: {}. Candidates checked: {}.</p><table><tr><th>Level</th><th>Baseline J</th><th>Baseline sampled Pa</th><th>Accepted sampled Pa</th><th>Baseline probes</th></tr>", check.outcome, check.checked);
    for (i, r) in check.baseline.rungs.iter().enumerate() {
        let candidate = check.candidate.as_ref().map_or_else(|| "not accepted".into(),
            |c| format!("{:.8e}", c.rungs[i].evaluation.sampled_max_von_mises));
        let s = &r.evaluation;
        let _ = write!(out, "<tr><td>{}</td><td>{:.8e}</td><td>{:.8e}</td><td>{candidate}</td><td>{}</td></tr>",
            r.level, s.compliance, s.sampled_max_von_mises, s.sample_count);
    }
    out.push_str("</table>"); out
}
fn read_report(value: &JsonValue, p: ResolutionPolicy, policy: &Controls) -> Result<StressResolutionReport> {
    let rows = value.get("rungs").and_then(JsonValue::as_array)
        .filter(|r| r.len() == p.extra_levels as usize + 1)
        .ok_or_else(|| malformed("incomplete stress mesh report"))?;
    let rungs = rows.iter().map(|row| Ok(StressResolutionRung {
        level: u32::try_from(integer(row, "level")?).map_err(|_| malformed("invalid stress mesh level"))?,
        evaluation: read_stress_measurement(row.get("evaluation")
            .ok_or_else(|| malformed("missing stress mesh evaluation"))?)?,
    })).collect::<Result<Vec<_>>>()?;
    let report = StressResolutionReport { rungs };
    // This validates measurements even for a genuinely unresolved report.
    let _ = report.refusal(p, policy.area.target, policy.stress).map_err(|e| malformed(&e.to_string()))?;
    Ok(report)
}
pub(super) fn read(value: &JsonValue, policy: &Controls, baseline: &SampledStressEvaluation,
    accepted: &[SampledStressEvaluation]) -> Result<Option<LastCheck>>
{
    let node = value.get("mesh_resolution");
    let Some(p) = policy.resolution else {
        return if node.is_none() { Ok(None) } else { Err(malformed("undeclared stress mesh evidence")) };
    };
    let node = node.ok_or_else(|| malformed("missing declared stress mesh evidence"))?;
    if node.str_field("schema") != Some("observed-stress-mesh-v1")
        || node.str_field("authority") != Some("Estimated")
        || node.get("policy") != Some(&document(policy_json(p).as_bytes())?) {
        return Err(malformed("retained stress mesh policy changed"));
    }
    let last = node.get("last_check").ok_or_else(|| malformed("missing stress mesh outcome"))?;
    if matches!(last, JsonValue::Null) {
        return if accepted.is_empty() { Ok(None) } else { Err(malformed("accepted stress updates lack mesh checks")) };
    }
    let outcome = match last.str_field("outcome") {
        Some("accepted") => "accepted", Some("baseline-unresolved") => "baseline-unresolved",
        Some("no-descent") => "no-descent", _ => return Err(malformed("unknown stress mesh outcome")),
    };
    let original = baseline;
    let baseline = read_report(last.get("baseline").ok_or_else(|| malformed("missing stress mesh baseline"))?, p, policy)?;
    let candidate = match last.get("accepted_candidate") {
        Some(JsonValue::Null) => None, Some(v) => Some(read_report(v, p, policy)?),
        None => return Err(malformed("missing stress mesh candidate")),
    };
    let checked = integer(last, "candidates_checked")?;
    let reason = match last.get("reason") {
        Some(JsonValue::Null) => None, Some(JsonValue::Str(s)) => Some(s.clone()),
        _ => return Err(malformed("invalid stress mesh refusal")),
    };
    let current = accepted.last().unwrap_or(original);
    let previous = if accepted.len() >= 2 { &accepted[accepted.len()-2] } else { original };
    if checked > policy.search.max_candidates || (outcome == "accepted") != candidate.is_some()
        || !same(&baseline.rungs[0].evaluation, if outcome == "accepted" { previous } else { current }) {
        return Err(malformed("stress mesh evidence disagrees with accepted history"));
    }
    let refused = baseline.refusal(p, policy.area.target, policy.stress).map_err(|e| malformed(&e.to_string()))?;
    if outcome == "baseline-unresolved" {
        if checked != 0 || refused.is_none() || refused != reason { return Err(malformed("unresolved stress baseline lacks its measured refusal")); }
    } else if refused.is_some() || reason.is_some() { return Err(malformed("stress search started from an unresolved baseline")); }
    if let Some(candidate) = &candidate {
        if accepted.is_empty() || checked == 0 || !same(&candidate.rungs[0].evaluation, current)
            || candidate.comparison_refusal(&baseline, p, policy.area.target, policy.stress,
                policy.search.min_relative_improvement).map_err(|e| malformed(&e.to_string()))?.is_some() {
            return Err(malformed("retained candidate violates finer-grid stress or compliance admission"));
        }
    }
    Ok(Some(LastCheck { outcome, baseline, candidate, checked, reason }))
}

#[cfg(test)]
#[path = "stress_mesh_tests.rs"]
mod tests;
