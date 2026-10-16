//! Native phase accounting for the existing sampled-stress restoration owner.
//! This is receipt admission and presentation, not a second search algorithm.
use super::*;

pub(in super::super) fn reduction(policy: &Controls) -> Option<f64> {
    policy.family.as_ref().and_then(|family| family.restoration_reduction)
}

pub(in super::super) fn feasible(state: &SampledStressEvaluation, policy: &Controls) -> bool {
    state.sampled_max_von_mises <= policy.stress.admitted_max()
}

pub(in super::super) fn baseline_scope(policy: &Controls) -> &'static str {
    if reduction(policy).is_some() { "area_feasible_study_start" }
    else { "feasible_study_start" }
}

/// Check a retained transition under the declared phase's actual measured gate.
/// In restoration, reducing compliance alone is not progress. After feasibility,
/// an improvement in stress cannot excuse a worse objective or renewed violation.
pub(in super::super) fn transition(
    previous: &SampledStressEvaluation, candidate: &SampledStressEvaluation, policy: &Controls,
) -> Result<bool> {
    let bound = policy.stress.admitted_max();
    if let Some(required) = reduction(policy).filter(|_| !feasible(previous, policy)) {
        if !feasible(candidate, policy)
            && !(candidate.sampled_max_von_mises - bound
                < (previous.sampled_max_von_mises - bound) * (1.0 - required))
        {
            return Err(malformed("accepted restoration did not reduce worst sampled stress excess"));
        }
        Ok(true)
    } else {
        if !feasible(previous, policy) || !feasible(candidate, policy)
            || !(candidate.compliance < previous.compliance * (1.0 - policy.search.min_relative_improvement))
        {
            return Err(malformed("accepted feasible descent violates stress or compliance admission"));
        }
        Ok(false)
    }
}

/// The crossing update belongs to restoration; it establishes the baseline for
/// future compliance improvement and must not itself be counted as that gain.
fn first_feasible<'a>(
    baseline: &'a SampledStressEvaluation, accepted: &'a [SampledStressEvaluation], policy: &Controls,
) -> Option<(usize, &'a SampledStressEvaluation)> {
    std::iter::once(baseline).chain(accepted).enumerate()
        .find(|(_, state)| feasible(state, policy))
}

pub(in super::super) fn updates(
    baseline: &SampledStressEvaluation, accepted: &[SampledStressEvaluation], policy: &Controls,
) -> usize {
    if reduction(policy).is_none() { return 0; }
    first_feasible(baseline, accepted, policy).map_or(accepted.len(), |(ordinal, _)| ordinal)
}

pub(in super::super) fn relative_reduction(
    baseline: &SampledStressEvaluation, accepted: &[SampledStressEvaluation], policy: &Controls,
) -> String {
    let Some((_, reference)) = first_feasible(baseline, accepted, policy) else {
        return "null".into();
    };
    let current = accepted.last().unwrap_or(baseline);
    let gain = if reference.compliance > 0.0 {
        (reference.compliance - current.compliance) / reference.compliance
    } else { 0.0 };
    format!("{gain:.17e}")
}

fn metadata(
    baseline: &SampledStressEvaluation, accepted: &[SampledStressEvaluation], policy: &Controls,
) -> Option<String> {
    let required = reduction(policy)?;
    let current = accepted.last().unwrap_or(baseline);
    let first = first_feasible(baseline, accepted, policy);
    let (ordinal, reference) = first.map_or_else(|| ("null".into(), "null".into()),
        |(ordinal, reference)| (ordinal.to_string(), stress_json(reference)));
    let feasible = feasible(current, policy);
    let phase = if feasible { "compliance-descent" } else { "restoring-stress" };
    let repairs = updates(baseline, accepted, policy);
    let excess = (current.sampled_max_von_mises - policy.stress.admitted_max()).max(0.0);
    Some(format!(concat!("{{\"minimum_relative_excess_reduction\":{required:.17e},",
        "\"phase\":\"{phase}\",\"stress_feasible\":{feasible},",
        "\"restoration_updates\":{repairs},\"compliance_updates\":{},",
        "\"first_feasible_update\":{ordinal},\"feasible_baseline\":{reference},",
        "\"remaining_excess_pa\":{excess:.17e}}}"), accepted.len() - repairs,
        required = required, phase = phase, feasible = feasible, repairs = repairs,
        ordinal = ordinal, reference = reference, excess = excess))
}

pub(in super::super) fn json_field(
    baseline: &SampledStressEvaluation, accepted: &[SampledStressEvaluation], policy: &Controls,
) -> String {
    metadata(baseline, accepted, policy).map_or_else(String::new,
        |value| format!(",\"stress_restoration\":{value}"))
}

pub(in super::super) fn check_retained(
    value: &JsonValue, baseline: &SampledStressEvaluation,
    accepted: &[SampledStressEvaluation], policy: &Controls,
) -> Result<()> {
    let expected = metadata(baseline, accepted, policy)
        .map(|text| document(text.as_bytes())).transpose()?;
    if value.get("stress_restoration") != expected.as_ref() {
        return Err(malformed("retained restoration phase or feasible baseline disagrees with measured history"));
    }
    Ok(())
}

pub(in super::super) fn html(
    baseline: &SampledStressEvaluation, accepted: &[SampledStressEvaluation], policy: &Controls,
) -> String {
    if reduction(policy).is_none() { return String::new(); }
    let status = first_feasible(baseline, accepted, policy).map_or_else(
        || "Stress remains infeasible; no compliance-improvement baseline exists.".to_string(),
        |(ordinal, _)| format!("First stress-feasible design: accepted update {ordinal}. Compliance improvement is measured only from that design."));
    format!("<p>Explicit stress restoration: {} accepted repair updates. {status} Restoration uses the governing case's unweighted compliance direction and admits measured worst-excess reduction; compliance may increase. This is not a stress adjoint or a guarantee of finding a feasible design. All phases share the original update and solve allowances.</p>",
        updates(baseline, accepted, policy))
}

#[cfg(test)]
#[path = "restoration_tests.rs"]
mod tests;
