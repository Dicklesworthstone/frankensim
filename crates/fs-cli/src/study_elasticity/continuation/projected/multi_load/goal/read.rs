//! Read the numerical DWR fields retained by a native load-family study.
use super::*;

fn envelope(policy: MultiLoadDwrPolicy, check: Option<&LastCheck>) -> Result<JsonValue> {
    document(format!("{{{}}}", field(Some(policy), check).trim_start_matches(',')).as_bytes())
}
fn measurement(value: &JsonValue, cases: &[RobustLoadCase]) -> Result<MultiLoadDwrMeasurement> {
    let rows = value.as_array().filter(|rows| rows.len() == cases.len())
        .ok_or_else(|| malformed("incomplete weighted DWR case family"))?;
    let mut measured = Vec::with_capacity(rows.len());
    for (row, &case) in rows.iter().zip(cases) {
        if row.get("load") != Some(&document(case_json(&case).as_bytes())?) {
            return Err(malformed("DWR load differs from the study"));
        }
        measured.push((case, DwrMeasurement {
            snapshot: hexadecimal(row.str_field("snapshot").ok_or_else(|| malformed("missing DWR snapshot"))?, true)?,
            level: u32::try_from(integer(row, "level")?).map_err(|_| malformed("invalid DWR level"))?,
            area: number(row, "area_m2")?.0,
            coarse_compliance: number(row, "coarse_compliance_j")?.0,
            enriched_compliance: number(row, "enriched_compliance_j")?.0,
            eta_signed: number(row, "eta_signed_j")?.0,
            coarse_relative_residual: number(row, "coarse_relative_residual")?.0,
            enriched_relative_residual: number(row, "enriched_relative_residual")?.0,
        }));
    }
    let measured = MultiLoadDwrMeasurement { cases: measured };
    measured.aggregate(cases).map_err(|e| malformed(&e.to_string()))?;
    Ok(measured)
}
fn state(values: &[SampledStressEvaluation], cases: &[RobustLoadCase]) -> Result<MultiLoadProjectedState> {
    let total = summary(values, cases, RobustAggregate::WeightedSum)?;
    Ok(MultiLoadProjectedState { objective: total.compliance, volume: total.volume, snapshot: total.snapshot,
        active_case: None, case_compliances: values.iter().map(|v| v.compliance).collect() })
}

pub(in super::super) fn read(value: &JsonValue, policy: &Controls, cases: &[RobustLoadCase],
    original: &[SampledStressEvaluation], accepted: &[Vec<SampledStressEvaluation>]) -> Result<Option<LastCheck>> {
    let node = value.get("goal_resolution");
    let Some(p) = policy.family.as_ref().and_then(|f| f.goal) else {
        return if node.is_none() { Ok(None) } else { Err(malformed("undeclared DWR evidence")) };
    };
    let node = node.ok_or_else(|| malformed("missing DWR evidence"))?;
    let last = node.get("last_check").ok_or_else(|| malformed("missing DWR check"))?;
    if matches!(last, JsonValue::Null) {
        if !accepted.is_empty() || envelope(p, None)?.get("goal_resolution") != Some(node) {
            return Err(malformed("accepted updates have no DWR check"));
        }
        return Ok(None);
    }
    let outcome = match last.str_field("outcome") {
        Some("accepted") => "accepted", Some("no-descent") => "no-descent", Some("solve-budget") => "solve-budget",
        _ => return Err(malformed("invalid DWR outcome")),
    };
    let baseline = measurement(last.get("baseline").ok_or_else(|| malformed("missing DWR baseline"))?, cases)?;
    let candidate = match last.get("accepted_candidate") {
        Some(JsonValue::Null) => None, Some(value) => Some(measurement(value, cases)?),
        None => return Err(malformed("missing DWR candidate")),
    };
    let checked = integer(last, "candidates_checked")?;
    if checked > policy.search.max_candidates || (outcome == "accepted") != candidate.is_some()
        || (outcome == "accepted" && (checked == 0 || accepted.is_empty())) {
        return Err(malformed("inconsistent DWR candidate history"));
    }
    let current = accepted.last().map(Vec::as_slice).unwrap_or(original);
    let previous = if accepted.len() >= 2 { accepted[accepted.len() - 2].as_slice() } else { original };
    let old_state = state(if outcome == "accepted" { previous } else { current }, cases)?;
    let level = baseline.aggregate(cases).map_err(|e| malformed(&e.to_string()))?.level;
    baseline.matches_state(cases, &old_state, level).map_err(|e| malformed(&e.to_string()))?;
    let decision = if let Some(candidate) = &candidate {
        candidate.matches_state(cases, &state(current, cases)?, level).map_err(|e| malformed(&e.to_string()))?;
        let d = compare_weighted_dwr_gains(&baseline, candidate, cases, p, policy.search.min_relative_improvement)
            .map_err(|e| malformed(&e.to_string()))?;
        if d.refusal.is_some() { return Err(malformed("DWR gain does not meet the declared policy")); }
        Some(d)
    } else { None };
    let check = LastCheck { outcome, baseline, candidate, decision, checked };
    if envelope(p, Some(&check))?.get("goal_resolution") != Some(node) {
        return Err(malformed("DWR summary differs from its measured values"));
    }
    Ok(Some(check))
}

pub(in super::super) fn validate_terminal(history: &History, policy: &Controls, status: &str, level: u32) -> Result<()> {
    if policy.family.as_ref().and_then(|family| family.goal).is_none() { return Ok(()) }
    let Some(check) = &history.goal else {
        return if history.accepted.is_empty() && !matches!(status, "completed" | "no-feasible-descent") {
            Ok(())
        } else { Err(malformed("finished updates have no DWR measurement")) };
    };
    let total = check.baseline.aggregate(&history.cases).map_err(|e| malformed(&e.to_string()))?;
    if total.level != level || (status == "completed" && check.outcome != "accepted")
        || (status == "no-feasible-descent" && check.outcome != "no-descent") {
        return Err(malformed("DWR level or terminal differs from the study"));
    }
    let minimum = history.cases.len() * (1 + 5 * history.accepted.len());
    if history.solves < minimum.max(3 * history.cases.len()) {
        return Err(malformed("DWR work count omits complete assessed families"));
    }
    Ok(())
}
