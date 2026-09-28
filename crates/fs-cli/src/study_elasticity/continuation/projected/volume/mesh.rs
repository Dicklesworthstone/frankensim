//! Optional native mesh-check declaration and projection of actual solver output.
//! This module owns no physics or candidate search. Checks run before acceptance
//! through fs-topols; prior complete checks survive a cancelled later attempt.
use super::*;
pub(super) use fs_topols::resolution::{MeshCheckStage, MeshCheckedProgress, ResolutionPolicy};
use fs_topols::resolution::{ResolutionReport, ResolutionRung};

pub(super) fn parse(fields: &[Node]) -> Result<Option<ResolutionPolicy>> {
    let Some(pair) = fields.windows(2).find(|pair|
        matches!(&pair[0].kind, NodeKind::Keyword(k) if k == "mesh-check"))
    else { return Ok(None) };
    let items = list(&pair[1], "mesh-check")?;
    if items.len() != 9 || !matches!(&items[0].kind, NodeKind::Symbol(k) if k == "refinement") {
        return Err(malformed("mesh-check requires refinement and four explicit policy fields"));
    }
    let keys = ["extra-levels", "absolute-compliance-tolerance-j", "relative-compliance-tolerance", "area-tolerance-m2"];
    for (i, key) in keys.iter().enumerate() {
        if !matches!(&items[2*i+1].kind, NodeKind::Keyword(k) if k.as_str() == *key) {
            return Err(malformed("mesh-check policy fields must be complete, unique and in canonical order"));
        }
    }
    let real = |i| super::super::super::super::number(&items[i], "mesh-check tolerance");
    let policy = ResolutionPolicy {
        extra_levels: u32::try_from(integer_node(&items[2], "extra-levels")?)
            .map_err(|_| malformed("mesh-check level count exceeds u32"))?,
        absolute_compliance_tolerance: real(4)?, relative_compliance_tolerance: real(6)?,
        area_tolerance: real(8)?,
    };
    policy.validate(1).map_err(|e| malformed(&e.to_string()))?;
    Ok(Some(policy))
}
pub(super) fn canonical(policy: ResolutionPolicy, out: &mut String) {
    let _ = writeln!(out, "    :mesh-check (refinement");
    let _ = writeln!(out, "      :extra-levels {}", policy.extra_levels);
    let _ = writeln!(out, "      :absolute-compliance-tolerance-j {}", canonical_float(policy.absolute_compliance_tolerance));
    let _ = writeln!(out, "      :relative-compliance-tolerance {}", canonical_float(policy.relative_compliance_tolerance));
    let _ = writeln!(out, "      :area-tolerance-m2 {})", canonical_float(policy.area_tolerance));
}
fn policy_json(p: ResolutionPolicy) -> String {
    format!("{{\"extra_levels\":{},\"absolute_compliance_tolerance_j\":{:.17e},\"relative_compliance_tolerance\":{:.17e},\"area_tolerance_m2\":{:.17e}}}",
        p.extra_levels, p.absolute_compliance_tolerance, p.relative_compliance_tolerance, p.area_tolerance)
}
fn report_json(report: &ResolutionReport) -> String {
    let rows = report.rungs.iter().map(|r| format!("{{\"level\":{},\"compliance_j\":{:.17e},\"area_m2\":{:.17e},\"snapshot\":\"{:#018x}\"}}",
        r.level, r.compliance, r.area, r.snapshot)).collect::<Vec<_>>().join(",");
    format!("{{\"rungs\":[{rows}],\"max_compliance_change_j\":{:.17e},\"max_area_change_m2\":{:.17e},\"agrees\":{}}}",
        report.max_compliance_change, report.max_area_change, report.agrees)
}

/// Only the last complete search is copied here. Earlier complete assessments
/// remain in the normal predecessor receipt chain. A candidate report is kept
/// only when that actual candidate passed all gates and was accepted.
pub(super) struct LastCheck {
    pub(super) outcome: &'static str,
    pub(super) baseline: ResolutionReport,
    candidate: Option<ResolutionReport>,
    checked: usize,
    pub(super) reason: Option<String>,
}
impl LastCheck {
    pub(super) fn capture(progress: MeshCheckedProgress) -> Result<(Option<ProjectedProgress>, Self)> {
        match progress {
            MeshCheckedProgress::UnresolvedBaseline { baseline, reason } => Ok((None, Self {
                outcome: "baseline-unresolved", baseline, candidate: None, checked: 0, reason: Some(reason),
            })),
            MeshCheckedProgress::Searched { baseline, candidates, progress } => {
                let (outcome, candidate) = match &progress {
                    ProjectedProgress::Accepted(_) => ("accepted", Some(candidates.last()
                        .filter(|c| c.refusal.is_none()).and_then(|c| c.report.clone())
                        .ok_or_else(|| malformed("accepted mesh-checked candidate has no completed assessment"))?)),
                    ProjectedProgress::NoDescent(_) => ("no-descent", None),
                    ProjectedProgress::IterationLimit => return Err(malformed("mesh check started beyond the update budget")),
                };
                let checked = candidates.len();
                Ok((Some(progress), Self { outcome, baseline, candidate, checked, reason: None }))
            }
        }
    }
    fn json(&self) -> String {
        format!("{{\"outcome\":{},\"baseline\":{},\"accepted_candidate\":{},\"candidates_checked\":{},\"reason\":{}}}",
            quoted(self.outcome), report_json(&self.baseline),
            self.candidate.as_ref().map_or_else(|| "null".into(), report_json), self.checked,
            self.reason.as_ref().map_or_else(|| "null".into(), |r| quoted(r)))
    }
    fn endpoint(&self) -> &ResolutionReport { self.candidate.as_ref().unwrap_or(&self.baseline) }
}

pub(super) fn field(policy: Option<ResolutionPolicy>, check: Option<&LastCheck>) -> String {
    policy.map_or_else(String::new, |policy| format!(concat!(",\"mesh_resolution\":{{",
        "\"schema\":\"observed-mesh-check-v1\",\"authority\":\"Estimated\",",
        "\"scope\":\"last-complete-search; measured grid sensitivity, not a continuum error bound\",",
        "\"policy\":{},\"last_check\":{}}}"), policy_json(policy),
        check.map_or_else(|| "null".into(), LastCheck::json)))
}
pub(super) fn html(policy: Option<ResolutionPolicy>, check: Option<&LastCheck>) -> String {
    let Some(policy) = policy else { return String::new() };
    let mut html = format!("<h2>Observed mesh-resolution check</h2><p>{} additional grid levels. Compliance allowance: {:.6e} J + {:.6e} times finer compliance; area allowance: {:.6e} m². Finer fields are not re-projected. Agreement is measured grid sensitivity, not a continuum error bound.</p>",
        policy.extra_levels, policy.absolute_compliance_tolerance, policy.relative_compliance_tolerance, policy.area_tolerance);
    let Some(check) = check else { html.push_str("<p>Not yet evaluated; no mesh-resolution pass is claimed.</p>"); return html };
    let _ = write!(html, "<p>Last complete search: {}. Candidates assessed: {}.</p><table><tr><th>Level</th><th>Baseline J</th><th>Accepted candidate J</th><th>Endpoint area m²</th></tr>", check.outcome, check.checked);
    for (i, r) in check.baseline.rungs.iter().enumerate() {
        let candidate = check.candidate.as_ref().map(|c| c.rungs[i].compliance)
            .map_or_else(|| "not accepted".into(), |v| format!("{v:.8e}"));
        let _ = write!(html, "<tr><td>{}</td><td>{:.8e}</td><td>{candidate}</td><td>{:.8e}</td></tr>",
            r.level, r.compliance, check.endpoint().rungs[i].area);
    }
    html.push_str("</table>"); html
}

fn read_report(value: &JsonValue, policy: ResolutionPolicy) -> Result<ResolutionReport> {
    let rows = value.get("rungs").and_then(JsonValue::as_array)
        .filter(|r| r.len() == policy.extra_levels as usize + 1)
        .ok_or_else(|| malformed("incomplete mesh-resolution levels"))?;
    let mut rungs = Vec::with_capacity(rows.len());
    for row in rows {
        let level = u32::try_from(integer(row, "level")?).map_err(|_| malformed("invalid mesh-check level"))?;
        rungs.push(ResolutionRung { level, compliance: number(row, "compliance_j")?.0,
            area: number(row, "area_m2")?.0, snapshot: hexadecimal(row.str_field("snapshot")
                .ok_or_else(|| malformed("missing mesh-check snapshot"))?, true)? });
    }
    policy.validate(rungs[0].level).map_err(|e| malformed(&e.to_string()))?;
    if rungs.iter().enumerate().any(|(i, r)| r.level != rungs[0].level + i as u32
        || r.compliance < 0.0 || r.area <= 0.0) { return Err(malformed("invalid mesh-resolution rung")); }
    let mut report = ResolutionReport { rungs, max_compliance_change: 0.0, max_area_change: 0.0, agrees: true };
    for pair in report.rungs.windows(2) {
        let dc = (pair[1].compliance-pair[0].compliance).abs();
        let da = (pair[1].area-pair[0].area).abs();
        let allowance = policy.absolute_compliance_tolerance + policy.relative_compliance_tolerance*pair[1].compliance.abs();
        if !dc.is_finite() || !da.is_finite() || !allowance.is_finite() { return Err(malformed("mesh-check arithmetic overflow")); }
        report.max_compliance_change = report.max_compliance_change.max(dc);
        report.max_area_change = report.max_area_change.max(da);
        report.agrees &= dc <= allowance && da <= policy.area_tolerance;
    }
    if &document(report_json(&report).as_bytes())? != value {
        return Err(malformed("mesh-check summary does not reproduce its measurements"));
    }
    Ok(report)
}
fn matches_state(report: &ResolutionReport, expected: Measured) -> bool {
    let r = report.rungs[0];
    r.snapshot == expected.snapshot && r.compliance.to_bits() == expected.compliance.to_bits()
        && r.area.to_bits() == expected.volume.to_bits()
}
pub(super) fn read(value: &JsonValue, policy: &Controls, previous: Measured, current: Measured,
    accepted_count: usize) -> Result<Option<LastCheck>> {
    let node = value.get("mesh_resolution");
    let Some(p) = policy.resolution else {
        return if node.is_none() { Ok(None) } else { Err(malformed("undeclared mesh-resolution evidence")) };
    };
    let node = node.ok_or_else(|| malformed("missing declared mesh-check evidence"))?;
    if node.get("policy") != Some(&document(policy_json(p).as_bytes())?)
        || node.str_field("schema") != Some("observed-mesh-check-v1") || node.str_field("authority") != Some("Estimated") {
        return Err(malformed("retained mesh-check policy changed"));
    }
    let last = node.get("last_check").ok_or_else(|| malformed("missing mesh-check outcome"))?;
    if matches!(last, JsonValue::Null) {
        return if accepted_count == 0 { Ok(None) } else { Err(malformed("accepted updates lack their declared mesh check")) };
    }
    let outcome = match last.str_field("outcome") {
        Some("accepted") => "accepted", Some("baseline-unresolved") => "baseline-unresolved",
        Some("no-descent") => "no-descent", _ => return Err(malformed("unknown mesh-check outcome")),
    };
    let baseline = read_report(last.get("baseline").ok_or_else(|| malformed("missing mesh-check baseline"))?, p)?;
    let candidate = match last.get("accepted_candidate") {
        Some(JsonValue::Null) => None, Some(v) => Some(read_report(v, p)?),
        None => return Err(malformed("missing mesh-check candidate state")),
    };
    let checked = integer(last, "candidates_checked")?;
    let reason = match last.get("reason") {
        Some(JsonValue::Null) => None, Some(JsonValue::Str(s)) => Some(s.clone()),
        _ => return Err(malformed("invalid mesh-check reason")),
    };
    if checked > policy.search.max_candidates || (outcome == "accepted") != candidate.is_some()
        || (outcome == "baseline-unresolved") != reason.is_some()
        || !matches_state(&baseline, if outcome == "accepted" { previous } else { current }) {
        return Err(malformed("mesh check disagrees with accepted study history"));
    }
    if let Some(candidate) = &candidate {
        if accepted_count == 0 || checked == 0 || !baseline.agrees || !candidate.agrees || !matches_state(candidate, current)
            || baseline.rungs.iter().zip(&candidate.rungs).enumerate().any(|(i, (old, new))| old.level != new.level
                || (i != 0 && ((old.area-policy.area.target).abs() > p.area_tolerance
                    || (new.area-policy.area.target).abs() > p.area_tolerance))
                || (old.area-new.area).abs() > p.area_tolerance
                || !(new.compliance < old.compliance*(1.0-policy.search.min_relative_improvement))) {
            return Err(malformed("retained mesh check does not admit the claimed improvement"));
        }
    }
    Ok(Some(LastCheck { outcome, baseline, candidate, checked, reason }))
}

#[cfg(test)]
#[path = "mesh_tests.rs"]
mod tests;
