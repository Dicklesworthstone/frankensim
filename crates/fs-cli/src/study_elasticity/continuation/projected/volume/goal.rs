//! Native signed-DWR admission and reporting on the existing volume-only path.
//! Physics and acceptance belong to fs-topols; this adapter binds controls,
//! completed evidence and conservative work charges to the existing receipts.
use super::*;
pub(super) use fs_topols::goal_resolution::{DwrCheckStage, DwrGatePolicy};
use fs_topols::goal_resolution::{DwrCheckedProgress, DwrGainDecision, DwrMeasurement, compare_dwr_gains};

pub(super) fn parse(fields: &[Node]) -> Result<Option<DwrGatePolicy>> {
    let Some(pair) = fields.windows(2).find(|pair|
        matches!(&pair[0].kind, NodeKind::Keyword(k) if k == "goal-check"))
    else { return Ok(None) };
    let items = list(&pair[1], "goal-check")?;
    if items.len() != 5
        || !matches!(&items[0].kind, NodeKind::Symbol(k) if k == "elasticity-dwr")
        || !matches!(&items[1].kind, NodeKind::Keyword(k) if k == "max-error-fraction")
        || !matches!(&items[3].kind, NodeKind::Keyword(k) if k == "max-solves-per-attempt")
    {
        return Err(malformed("goal-check requires elasticity-dwr, max-error-fraction and max-solves-per-attempt in that order"));
    }
    let policy = DwrGatePolicy {
        max_error_fraction: super::super::super::super::number(&items[2], "max-error-fraction")?,
        max_solves_per_attempt: integer_node(&items[4], "max-solves-per-attempt")?,
    };
    policy.validate(1, 1).map_err(|error| malformed(&error.to_string()))?;
    Ok(Some(policy))
}

pub(super) fn canonical(policy: DwrGatePolicy, out: &mut String) {
    let _ = writeln!(out, "    :goal-check (elasticity-dwr");
    let _ = writeln!(out, "      :max-error-fraction {}", canonical_float(policy.max_error_fraction));
    let _ = writeln!(out, "      :max-solves-per-attempt {})", policy.max_solves_per_attempt);
}

fn measurement_json(m: DwrMeasurement) -> String {
    format!(concat!("{{\"snapshot\":\"{:#018x}\",\"level\":{},\"area_m2\":{:.17e},",
        "\"coarse_compliance_j\":{:.17e},\"enriched_compliance_j\":{:.17e},\"eta_signed_j\":{:.17e},",
        "\"coarse_relative_residual\":{:.17e},\"enriched_relative_residual\":{:.17e}}}"),
        m.snapshot, m.level, m.area, m.coarse_compliance, m.enriched_compliance,
        m.eta_signed, m.coarse_relative_residual, m.enriched_relative_residual)
}
fn decision_json(d: DwrGainDecision) -> String {
    format!(concat!("{{\"coarse_gain_j\":{:.17e},\"enriched_gain_j\":{:.17e},",
        "\"signed_correction_change_j\":{:.17e},\"observed_gain_change_j\":{:.17e},",
        "\"allowance_j\":{:.17e},\"refusal\":{}}}"),
        d.coarse_gain, d.enriched_gain, d.signed_correction_change, d.observed_gain_change,
        d.allowance, d.refusal.map_or_else(|| "null".into(), quoted))
}
fn read_measurement(value: &JsonValue) -> Result<DwrMeasurement> {
    let m = DwrMeasurement {
        snapshot: hexadecimal(value.str_field("snapshot").ok_or_else(|| malformed("missing DWR snapshot"))?, true)?,
        level: u32::try_from(integer(value, "level")?).map_err(|_| malformed("DWR level exceeds u32"))?,
        area: number(value, "area_m2")?.0,
        coarse_compliance: number(value, "coarse_compliance_j")?.0,
        enriched_compliance: number(value, "enriched_compliance_j")?.0,
        eta_signed: number(value, "eta_signed_j")?.0,
        coarse_relative_residual: number(value, "coarse_relative_residual")?.0,
        enriched_relative_residual: number(value, "enriched_relative_residual")?.0,
    };
    m.validate().map_err(|error| malformed(&error.to_string()))?;
    if &document(measurement_json(m).as_bytes())? != value {
        return Err(malformed("DWR measurement does not reproduce its canonical representation"));
    }
    Ok(m)
}
fn matches_state(m: DwrMeasurement, state: Measured) -> bool {
    m.snapshot == state.snapshot && m.coarse_compliance.to_bits() == state.compliance.to_bits()
        && m.area.to_bits() == state.volume.to_bits()
}

/// Only a completed search replaces this evidence. Earlier complete searches
/// remain in the normal predecessor chain. Rejections remain in terminal_refusals.
pub(super) struct LastCheck {
    outcome: &'static str,
    pub(super) baseline: DwrMeasurement,
    candidate: Option<DwrMeasurement>,
    decision: Option<DwrGainDecision>,
    checked: usize,
    reserved: usize,
}
impl LastCheck {
    pub(super) fn capture(checked: DwrCheckedProgress) -> Result<(ProjectedProgress, Self)> {
        let DwrCheckedProgress::Searched { baseline, candidates, reserved_solves, progress } = checked
            else { return Err(malformed("DWR check started beyond the update budget")); };
        if reserved_solves != 2 * (1 + candidates.len()) {
            return Err(malformed("DWR work debit disagrees with the completed assessment family"));
        }
        let (outcome, candidate, decision) = match &progress {
            ProjectedProgress::Accepted(step) => {
                let last = candidates.last().filter(|c| c.refusal.is_none())
                    .ok_or_else(|| malformed("accepted DWR candidate lacks a passed assessment"))?;
                let m = last.measurement.ok_or_else(|| malformed("accepted DWR measurement missing"))?;
                let d = last.decision.filter(|d| d.refusal.is_none())
                    .ok_or_else(|| malformed("accepted DWR decision missing"))?;
                if !matches_state(baseline, step.previous.into()) || !matches_state(m, step.state.into()) {
                    return Err(malformed("DWR check does not bind the actual accepted update"));
                }
                ("accepted", Some(m), Some(d))
            }
            ProjectedProgress::NoDescent(_) => ("no-descent", None, None),
            ProjectedProgress::IterationLimit => return Err(malformed("DWR search exceeded its update budget")),
        };
        Ok((progress, Self { outcome, baseline, candidate, decision,
            checked: candidates.len(), reserved: reserved_solves }))
    }
    fn json(&self) -> String {
        format!("{{\"outcome\":{},\"baseline\":{},\"accepted_candidate\":{},\"decision\":{},\"candidates_checked\":{},\"reserved_solves\":{}}}",
            quoted(self.outcome), measurement_json(self.baseline),
            self.candidate.map_or_else(|| "null".into(), measurement_json),
            self.decision.map_or_else(|| "null".into(), decision_json), self.checked, self.reserved)
    }
}

#[derive(Default)]
pub(super) struct GoalEvidence {
    pub(super) last: Option<LastCheck>,
    attempts: usize,
    pub(super) failure: Option<String>,
}
impl GoalEvidence {
    /// Charge the declared worst case BEFORE physics. The caller persists this
    /// charge before starting the attempt; even an interrupted attempt is not free.
    pub(super) fn begin(&mut self, policy: DwrGatePolicy) -> Result<()> {
        let attempts = self.attempts.checked_add(1).ok_or_else(|| malformed("DWR attempt counter overflow"))?;
        attempts.checked_mul(policy.max_solves_per_attempt)
            .ok_or_else(|| malformed("DWR reserved-work counter overflow"))?;
        self.attempts = attempts;
        self.failure = None;
        Ok(())
    }
    fn json(&self, policy: DwrGatePolicy) -> String {
        format!(concat!("{{\"authority\":\"Estimated\",\"policy\":{{\"max_error_fraction\":{:.17e},",
            "\"max_solves_per_attempt\":{}}},\"attempts_admitted\":{},\"total_reserved_solves\":{},",
            "\"last_check\":{},\"failure\":{}}}"),
            policy.max_error_fraction, policy.max_solves_per_attempt, self.attempts,
            self.attempts * policy.max_solves_per_attempt,
            self.last.as_ref().map_or_else(|| "null".into(), LastCheck::json),
            self.failure.as_ref().map_or_else(|| "null".into(), |s| quoted(s)))
    }
}

pub(super) fn field(policy: Option<DwrGatePolicy>, evidence: Option<&GoalEvidence>) -> String {
    policy.map_or_else(String::new, |p| format!(",\"goal_resolution\":{}",
        evidence.expect("declared DWR evidence").json(p)))
}
pub(super) fn html(policy: Option<DwrGatePolicy>, evidence: Option<&GoalEvidence>) -> String {
    let Some(p) = policy else { return String::new() };
    let e = evidence.expect("declared DWR evidence");
    let mut out = format!("<h2>Signed-DWR improvement check</h2><p>Coarse and enriched compliance must both decrease. Signed-correction and observed-gain changes must each fit {:.6e} times the coarse gain. {} attempts admitted, {} extra solves conservatively reserved; this is not a count of completed solves. Two DWR solves and residual integration are indivisible between cancellation checks. The last complete check survives interruption.</p><p>Estimated numerical evidence, not a continuum error bound, optimality test or physical validation. Common bias may survive. The material constraint is coarse-grid cut quadrature; enrichment does not re-project the geometry. Absolute indicator mass is not an error radius.</p>",
        p.max_error_fraction, e.attempts, e.attempts * p.max_solves_per_attempt);
    if let Some(c) = &e.last {
        let _ = write!(out, "<p>Last complete search: {}. Baseline signed DWR: {:.8e} J; coarse compliance: {:.8e} J; enriched compliance: {:.8e} J. Candidates assessed: {}.</p>",
            c.outcome, c.baseline.eta_signed, c.baseline.coarse_compliance, c.baseline.enriched_compliance, c.checked);
        if let Some(d) = c.decision {
            let _ = write!(out, "<p>Accepted coarse gain: {:.8e} J; enriched gain: {:.8e} J; signed-correction change: {:.8e} J; observed gain change: {:.8e} J; allowance: {:.8e} J.</p>",
                d.coarse_gain, d.enriched_gain, d.signed_correction_change, d.observed_gain_change, d.allowance);
        }
    } else { out.push_str("<p>No completed DWR search is available; no passed check is claimed.</p>"); }
    if e.failure.is_some() { out.push_str("<p>Numerical assessment failed. The retained JSON contains the refusal; accepted geometry is unchanged.</p>"); }
    out
}

pub(super) fn read(value: &JsonValue, policy: &Controls, previous: Measured, current: Measured,
    accepted_count: usize) -> Result<Option<GoalEvidence>> {
    let node = value.get("goal_resolution");
    let Some(p) = policy.goal else {
        return if node.is_none() { Ok(None) } else { Err(malformed("undeclared DWR gate evidence")) };
    };
    let node = node.ok_or_else(|| malformed("missing declared DWR gate evidence"))?;
    let attempts = integer(node, "attempts_admitted")?;
    attempts.checked_mul(p.max_solves_per_attempt).ok_or_else(|| malformed("DWR work counter overflow"))?;
    let failure = match node.get("failure") {
        Some(JsonValue::Null) => None,
        Some(JsonValue::Str(s)) if !s.is_empty() && attempts > 0 => Some(s.clone()),
        _ => return Err(malformed("invalid DWR failure record")),
    };
    let last = match node.get("last_check") {
        Some(JsonValue::Null) if accepted_count == 0 => None,
        Some(last) if attempts > 0 => {
            let baseline = read_measurement(last.get("baseline").ok_or_else(|| malformed("missing DWR baseline"))?)?;
            let checked = integer(last, "candidates_checked")?;
            if checked > policy.search.max_candidates { return Err(malformed("excess DWR candidate count")); }
            let (outcome, candidate, decision) = match last.str_field("outcome") {
                Some("accepted") if accepted_count > 0 && checked > 0 => {
                    let candidate = read_measurement(last.get("accepted_candidate")
                        .ok_or_else(|| malformed("missing accepted DWR candidate"))?)?;
                    let decision = compare_dwr_gains(baseline, candidate, p, policy.search.min_relative_improvement)
                        .map_err(|error| malformed(&error.to_string()))?;
                    if decision.refusal.is_some() || !matches_state(baseline, previous)
                        || !matches_state(candidate, current) {
                        return Err(malformed("retained DWR evidence does not admit the accepted update"));
                    }
                    ("accepted", Some(candidate), Some(decision))
                }
                Some("no-descent") if matches_state(baseline, current) => ("no-descent", None, None),
                _ => return Err(malformed("DWR check disagrees with the accepted trajectory")),
            };
            Some(LastCheck { outcome, baseline, candidate, decision, checked, reserved: 2 * (1 + checked) })
        }
        _ => return Err(malformed("accepted updates lack their declared DWR check")),
    };
    if attempts < accepted_count { return Err(malformed("DWR work charge omits accepted updates")); }
    let evidence = GoalEvidence { last, attempts, failure };
    if &document(evidence.json(p).as_bytes())? != node {
        return Err(malformed("DWR policy, decision or work summary does not reproduce its measurements"));
    }
    Ok(Some(evidence))
}

#[cfg(test)]
#[path = "goal_tests.rs"]
mod tests;
