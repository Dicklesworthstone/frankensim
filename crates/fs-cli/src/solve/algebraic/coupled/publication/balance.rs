//! Reuse the physical cooling publisher when adaptive mesh comparisons request
//! a tighter algebraic allowance. Primal work belongs to the mesh rung, not to
//! each request. The pair driver must re-probe every changed field/boundary.

use crate::json_read::JsonValue as Json;
use crate::solve::algebraic::balance::{LinearWork, MAX_RETARGET_CALLS, maximum_bound};
use super::{
    ConductionProblem, Cx, LinearGoalAnalysisConfig, MaximumEvidence, PhysicalAirMaximumPolish,
    RungSolved, SolveRefusal, cancelled, error, prepare_with_budget, receipt,
};

/// All attempts survive physical rejection; preparation has at most five
/// separately bounded calls, while every primal correction shares ONE cap.
#[derive(Clone)]
pub(in crate::solve) struct Work {
    initial_bound_k: Option<f64>,
    primal_limit: usize,
    primal_iterations: usize,
    retarget_calls: usize,
    response_iterations: usize,
    defect_corrections: usize,
    goal_checks: usize,
    physical_checks: usize,
    physical_rejections: usize,
    spectral_work_entries: usize,
    spectral_work_limit_per_call: usize,
    peak_spectral_storage_entries: usize,
    last_target_k: f64,
}

impl Work {
    pub(super) fn new(
        result: &PhysicalAirMaximumPolish, target: f64, limit: usize, spectral_limit: usize,
    ) -> Self {
        let solve = &result.correction.solution.solid;
        Self {
            initial_bound_k: result.correction.initial_analysis.algebraic_half_width_k(),
            primal_limit: limit,
            primal_iterations: solve.primal_iterations,
            retarget_calls: 0,
            response_iterations: solve.analysis.response_iterations(),
            defect_corrections: solve.defect_corrections,
            goal_checks: solve.goal_checks,
            physical_checks: result.physical_checks,
            physical_rejections: result.physical_rejections,
            spectral_work_entries: result.correction.preparation.as_ref().map_or(0, |p| p.work_entries),
            spectral_work_limit_per_call: spectral_limit,
            peak_spectral_storage_entries: result.correction.preparation.as_ref().map_or(0, |p| p.peak_storage_entries),
            last_target_k: target,
        }
    }

    pub(super) fn remaining(&self) -> Result<usize, SolveRefusal> {
        let remaining = self.primal_limit.checked_sub(self.primal_iterations)
            .ok_or_else(|| error("coupled rung exceeded its original primal limit"))?;
        Ok(if self.retarget_calls >= MAX_RETARGET_CALLS { 0 } else { remaining })
    }

    /// Merge only after the numerical call, physical gates and receipt building
    /// finish. Counter overflow never publishes a partially accumulated state.
    pub(super) fn accumulate(&self, next: &Self) -> Result<Self, SolveRefusal> {
        if self.retarget_calls >= MAX_RETARGET_CALLS || next.retarget_calls != 0
            || self.primal_limit != next.primal_limit
            || self.spectral_work_limit_per_call != next.spectral_work_limit_per_call
            || next.primal_iterations > self.remaining()?
        { return Err(error("adaptive cooling attempt changed or exceeded its retained allowance")); }
        let add = |a: usize, b: usize| a.checked_add(b)
            .ok_or_else(|| error("adaptive cooling work count overflow"));
        Ok(Self {
            initial_bound_k: self.initial_bound_k,
            primal_limit: self.primal_limit,
            primal_iterations: add(self.primal_iterations, next.primal_iterations)?,
            retarget_calls: self.retarget_calls + 1,
            response_iterations: add(self.response_iterations, next.response_iterations)?,
            defect_corrections: add(self.defect_corrections, next.defect_corrections)?,
            goal_checks: add(self.goal_checks, next.goal_checks)?,
            physical_checks: add(self.physical_checks, next.physical_checks)?,
            physical_rejections: add(self.physical_rejections, next.physical_rejections)?,
            spectral_work_entries: add(self.spectral_work_entries, next.spectral_work_entries)?,
            spectral_work_limit_per_call: self.spectral_work_limit_per_call,
            peak_spectral_storage_entries: self.peak_spectral_storage_entries.max(next.peak_spectral_storage_entries),
            last_target_k: next.last_target_k,
        })
    }

    pub(super) fn control(
        &self, current: &str, bound: Option<f64>, final_pair: Option<(f64, f64)>,
        last_attempt: Option<&str>, mut checkpoint: impl FnMut() -> Result<(), SolveRefusal>,
    ) -> Result<String, SolveRefusal> {
        checkpoint()?;
        let mut value = parse(current)?;
        let target = final_pair.map_or(self.last_target_k, |(_, allowance)| allowance);
        let met = bound.is_some_and(|bound| bound <= target);
        let preparation_limit = self.primal_limit.checked_mul(MAX_RETARGET_CALLS + 1)
            .ok_or_else(|| error("adaptive cooling preparation allowance overflow"))?;
        let spectral_limit = self.spectral_work_limit_per_call.checked_mul(MAX_RETARGET_CALLS + 1)
            .ok_or_else(|| error("adaptive cooling verification allowance overflow"))?;
        for (key, count) in [
            ("primal_iterations", self.primal_iterations), ("max_primal_iterations", self.primal_limit),
            ("response_iterations", self.response_iterations), ("max_response_iterations", preparation_limit),
            ("max_response_iterations_per_call", self.primal_limit),
            ("max_stability_iterations_per_call", self.primal_limit),
            ("max_stability_iterations", preparation_limit),
            ("defect_corrections", self.defect_corrections),
            ("goal_checks", self.goal_checks), ("preliminary_goal_checks", self.retarget_calls + 1),
            ("physical_checks", self.physical_checks), ("physical_rejections", self.physical_rejections),
            ("spectral_work_entries", self.spectral_work_entries), ("max_spectral_work_entries", spectral_limit),
            ("peak_spectral_storage_entries", self.peak_spectral_storage_entries),
            ("retarget_calls", self.retarget_calls), ("max_retarget_calls", MAX_RETARGET_CALLS),
        ] { receipt::set(&mut value, key, integer(count))?; }
        receipt::set(&mut value, "initial_bound_k", optional(self.initial_bound_k)?)?;
        receipt::set(&mut value, "final_bound_k", optional(bound)?)?;
        receipt::set(&mut value, "requested_tolerance_k", receipt::numeric(target)?)?;
        receipt::set(&mut value, "last_correction_target_k", receipt::numeric(self.last_target_k)?)?;
        receipt::set(&mut value, "goal_met", Json::Bool(met))?;
        receipt::set(&mut value, "work_scope", Json::Str(
            "cumulative per-rung work; candidate, physical refusal, inverse and reference-origin evidence describe the selected published field; spectral_preparation describes its preparation attempt".into()))?;
        if let Some(attempt) = last_attempt {
            receipt::set(&mut value, "last_retarget_attempt", parse(attempt)?)?;
        }
        if let Some((discretization, allowance)) = final_pair {
            receipt::set(&mut value, "status", Json::Str(if met {
                "discretization-allowance-met"
            } else { "discretization-allowance-unresolved" }.into()))?;
            receipt::set(&mut value, "adaptive_balance", Json::Object(vec![
                ("tolerance_basis".into(), Json::Str("measured-adaptive-discretization".into())),
                ("discretization_half_width_k".into(), receipt::numeric(discretization)?),
                ("pair_allowance_k".into(), receipt::numeric(allowance)?),
                ("scope".into(), Json::Str("coarse-plus-fine bound must fit one tenth of the final re-probed Estimated mesh term; each adopted field passes physical cooling gates; not a continuum certificate".into())),
            ]))?;
        }
        let mut out = String::new();
        receipt::encode(&value, &mut out, &mut checkpoint)?;
        checkpoint()?;
        Ok(out)
    }
}

fn integer(value: usize) -> Json { Json::Number { value: value as f64, raw: value.to_string() } }
fn optional(value: Option<f64>) -> Result<Json, SolveRefusal> {
    value.map(receipt::numeric).transpose().map(|value| value.unwrap_or(Json::Null))
}
fn parse(text: &str) -> Result<Json, SolveRefusal> {
    if text.len() > 1024 * 1024 { return Err(error("coupled control receipt exceeds one MiB")); }
    Json::parse(text).map_err(|e| error(e.to_string()))
}
fn control(evidence: &MaximumEvidence) -> Result<&str, SolveRefusal> {
    evidence.control_json.as_deref().ok_or_else(|| error("coupled rung lost its control receipt"))
}

pub(in crate::solve::algebraic) fn retarget(
    cx: &Cx<'_>, solved: &mut RungSolved, vertices: &[usize], target: f64, memory_bytes: u64,
) -> Result<bool, SolveRefusal> {
    cx.checkpoint().map_err(|_| cancelled())?;
    let Some(LinearWork::Coupled(old)) = &solved.algebraic.linear_work else {
        return Err(error("adaptive cooling has no retained coupled work"));
    };
    if old.primal_iterations != solved.algebraic.primal_iterations {
        return Err(error("adaptive cooling primal counters disagree"));
    }
    let remaining = old.remaining()?;
    if remaining == 0 { return Ok(false); }
    let data = solved.adjoint_data.as_ref().ok_or_else(|| error("missing retained cooling operator"))?;
    let problem = ConductionProblem { mesh: &solved.mesh, boundary: &data.boundary,
        material: &data.fallback, element_materials: Some(&data.materials), source: &data.source };
    let config = LinearGoalAnalysisConfig {
        residual_limits: fs_solver::goal::GoalResidualLimits {
            max_rows: usize::try_from(memory_bytes / 256).unwrap_or(usize::MAX),
            max_nonzeros: usize::try_from(memory_bytes / 64).unwrap_or(usize::MAX),
        }, max_stability_iterations: old.primal_limit,
    };
    let history = solved.conjugate_fragment.as_deref().ok_or_else(|| error("missing retained air receipt"))?;
    let publication = prepare_with_budget(cx, problem, data.interfaces.as_ref(), &data.air_paths,
        data.linear, &solved.solution, history, vertices, memory_bytes, config, target, remaining)?;
    let Some(LinearWork::Coupled(next)) = &publication.evidence.linear_work else {
        return Err(error("adaptive cooling attempt lost its work counters"));
    };
    let totals = old.accumulate(next)?;
    let changed = publication.replacement.as_ref().is_some_and(|replacement| {
        replacement.boundary != data.boundary
            || replacement.solution.temperature.iter().zip(&solved.solution.temperature)
                .any(|(a, b)| a.to_bits() != b.to_bits())
    });
    let old_bound = maximum_bound(&solved.algebraic);
    let new_bound = maximum_bound(&publication.evidence);
    // An unchanged field may already have stronger evidence. Preserve the
    // ENTIRE selected analysis/control, not just a bound from another origin.
    let keep_old = !changed && old_bound.is_some_and(|old| new_bound.is_none_or(|new| old <= new));
    let selected = if keep_old { &solved.algebraic } else { &publication.evidence };
    let next_control = totals.control(control(selected)?, if keep_old { old_bound } else { new_bound },
        None, keep_old.then_some(control(&publication.evidence)?),
        || cx.checkpoint().map_err(|_| cancelled()))?;
    cx.checkpoint().map_err(|_| cancelled())?;
    let data = solved.adjoint_data.as_mut().ok_or_else(|| error("retained cooling operator disappeared"))?;
    if !keep_old { solved.algebraic.term = publication.evidence.term; }
    if let Some(replacement) = publication.replacement {
        solved.solution = replacement.solution;
        data.boundary = replacement.boundary;
        solved.conjugate_fragment = Some(replacement.conjugate);
    }
    solved.algebraic.control_json = Some(next_control);
    solved.algebraic.primal_iterations = totals.primal_iterations;
    solved.algebraic.linear_work = Some(LinearWork::Coupled(totals));
    Ok(changed)
}

pub(in crate::solve::algebraic) fn finish(
    solved: &mut RungSolved, discretization_k: f64, allowance_k: f64,
) -> Result<(), SolveRefusal> {
    // The common entry point has validated the downward-rounded pair allowance.
    let Some(LinearWork::Coupled(totals)) = &solved.algebraic.linear_work else {
        return Err(error("adaptive cooling lost its final work accounting"));
    };
    let next_control = totals.control(control(&solved.algebraic)?, maximum_bound(&solved.algebraic),
        Some((discretization_k, allowance_k)), None, || Ok(()))?;
    solved.algebraic.control_json = Some(next_control);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn work(primal_iterations: usize) -> Work {
        Work { initial_bound_k: Some(2.0), primal_limit: 10, primal_iterations,
            retarget_calls: 0, response_iterations: 2, defect_corrections: 0,
            goal_checks: 1, physical_checks: 1, physical_rejections: 1,
            spectral_work_entries: 3, spectral_work_limit_per_call: 100,
            peak_spectral_storage_entries: 4, last_target_k: 0.1 }
    }

    #[test]
    fn zero_progress_calls_and_rejected_iterations_share_hard_rung_limits() {
        let exhausted = work(6).accumulate(&work(4)).unwrap();
        assert_eq!(exhausted.remaining().unwrap(), 0);
        assert_eq!(exhausted.physical_rejections, 2);
        assert!(exhausted.accumulate(&work(1)).is_err());
        let mut calls = work(0);
        for _ in 0..MAX_RETARGET_CALLS { calls = calls.accumulate(&work(0)).unwrap(); }
        assert_eq!(calls.remaining().unwrap(), 0);
        assert!(calls.accumulate(&work(0)).is_err());
        assert_eq!(calls.spectral_work_entries, 3 * (MAX_RETARGET_CALLS + 1));
    }

    #[test]
    fn changed_preparation_policy_and_overflow_refuse_without_altering_totals() {
        let original = work(1);
        let mut changed = work(1);
        changed.spectral_work_limit_per_call += 1;
        assert!(original.accumulate(&changed).is_err());
        let mut overflowing = work(1);
        overflowing.response_iterations = usize::MAX;
        assert!(original.accumulate(&overflowing).is_err());
        assert_eq!(original.primal_iterations, 1);
        assert_eq!(original.retarget_calls, 0);
    }
}
