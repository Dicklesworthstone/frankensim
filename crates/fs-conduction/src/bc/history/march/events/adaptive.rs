//! Event-limited step doubling of the existing implicit thermal endpoints.
//!
//! Compare one backward-Euler interval with two half intervals from the SAME
//! accepted state. Publish the two-half-step trajectory, never a Richardson
//! extrapolation that could violate positivity, phase or the energy balance.
//! The discrepancy estimates local time error; it is not a continuum/global
//! error bound. Solver residual/energy tolerances remain independent and should
//! be tighter than the requested time accuracy. Any solver refusal is returned
//! unchanged, not treated as an error estimate or an excuse to widen its gate.

use super::{boundary_before, next_endpoint, source_before};
use super::super::{ThermalHistoryMethod, ThermalHistoryRun, ThermalHistoryStep, invalid, poll};
use crate::{ConductionError, ConductionProblem};
use crate::transient::backward_euler::StepConfig;
use fs_exec::Cx;

/// Explicit local-temperature accuracy and bounded time-step search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AdaptiveHistoryConfig {
    /// Positive absolute allowance for the estimated local temperature error, K.
    pub absolute_temperature_k: f64,
    /// Relative allowance against the larger coarse/fine temperature CHANGE
    /// from history, not the absolute Kelvin offset. Must be in [0, 1).
    pub relative_temperature: f64,
    /// Smallest controller-selected interval, s. A closer schedule knot or
    /// requested final time may shorten the initial interval below this floor;
    /// that interval is tested once and is never accepted merely for being small.
    pub min_step_s: f64,
    /// Largest proposed interval, s, also capped by every schedule event.
    pub max_step_s: f64,
    /// Maximum full-versus-two-half attempts. Each uses at most three original
    /// endpoint solves with their unchanged individual work budgets.
    pub max_attempts: usize,
}

impl AdaptiveHistoryConfig {
    /// Validate scalar policies and their checked total endpoint-solve count.
    pub fn validate(self) -> Result<(), ConductionError> {
        if !self.absolute_temperature_k.is_finite() || self.absolute_temperature_k <= 0.0
            || !self.relative_temperature.is_finite() || !(0.0..1.0).contains(&self.relative_temperature)
            || !self.min_step_s.is_finite() || self.min_step_s <= 0.0
            || !self.max_step_s.is_finite() || self.max_step_s < self.min_step_s
            || self.max_attempts == 0 || self.max_attempts.checked_mul(3).is_none()
        { return Err(invalid("adaptive history requires positive finite kelvin/step floors, relative tolerance in [0,1), ordered step limits and a bounded positive attempt count")); }
        Ok(())
    }
}

/// One accepted pair of physical intervals and their original solver reports.
#[derive(Debug, Clone)]
pub struct AdaptiveHistoryStep<T> {
    /// Original accepted physical time, s.
    pub start_time_s: f64,
    /// End of the first accepted half step, s.
    pub midpoint_time_s: f64,
    /// End of the second half step; may precede the requested target, s.
    pub end_time_s: f64,
    /// Original report for the first half step, including its own energy audit.
    pub first_half: T,
    /// Original report for the second half step, starting from the first result.
    pub second_half: T,
    /// Maximum scaled step-doubling discrepancy, at most one on success.
    /// Includes the first-order correction for unequal representable half steps.
    pub normalized_error: f64,
    /// Rejected temporal-error trials; no rejected trial advances physical time.
    pub rejected_attempts: usize,
    /// Actual full/half endpoint solves, including discarded coarse/rejected work.
    pub endpoint_solves: usize,
    /// Suggested next proposal within the original step limits, s. The next
    /// call still caps this at its target and actual schedule knots.
    pub suggested_step_s: f64,
}

impl ThermalHistoryRun<'_, '_> {
    /// Advance an event-limited interval with local step-doubling error control.
    ///
    /// Start with `proposed_step_s` (at least the configured floor, capped at
    /// the configured ceiling). A rejected discrepancy halves the interval,
    /// stopping at the floor rather than accepting an inaccurate result. An
    /// accepted discrepancy <= 1/8 suggests doubling the next proposal; all
    /// other accepted intervals suggest retaining their size. This simple
    /// deterministic controller has no hidden memory that can be lost on resume.
    ///
    /// Both accepted half steps retain their energy/flux reports and consume
    /// TWO lifetime-budget slots. All coarse/rejected work is reported but
    /// consumes no accepted slots. The pair publishes atomically after a final
    /// cancellation check. There must be a representable interior midpoint;
    /// direct event stepping remains available for adjacent-float intervals.
    ///
    /// # Errors
    /// Invalid policy/time, fewer than two remaining slots, cancellation,
    /// original endpoint-solver refusal, or unmet time accuracy at the floor or
    /// attempt cap. Every error leaves accepted temperature, time and count intact.
    pub fn advance_adaptive_step(&mut self, cx: &Cx<'_>, target_s: f64,
        proposed_step_s: f64, step: StepConfig, adaptive: AdaptiveHistoryConfig)
        -> Result<AdaptiveHistoryStep<ThermalHistoryStep>, ConductionError>
    {
        let start = self.time_s;
        let cap = admit_proposal(cx, self.accepted_steps, self.max_steps, proposed_step_s, adaptive)?;
        let end = next_endpoint(cx, self.problem.boundary, &self.problem.source,
            std::iter::empty(), start, target_s, cap)?;
        let solve = |old: &[f64], from: f64, to: f64| -> Result<ThermalHistoryStep, ConductionError> {
            let boundary = boundary_before(cx, self.problem.boundary, to)?;
            let source = source_before(&self.problem.source, to, old.len())?;
            let problem = ConductionProblem { mesh: self.problem.boundary.mesh(), boundary: &boundary,
                material: self.problem.material, source: &source, element_materials: self.problem.element_materials };
            match self.method {
                ThermalHistoryMethod::Linear => self.stepper.advance_prescribed(cx, problem,
                    self.problem.interfaces, old, to - from, step).map(ThermalHistoryStep::Linear),
                ThermalHistoryMethod::Nonlinear(policy) => self.stepper.advance_nonlinear_prescribed(cx,
                    problem, self.problem.interfaces, old, to - from, step, policy).map(ThermalHistoryStep::Nonlinear),
            }
        };
        let selected = select_step(cx, start, end, adaptive, |e| e, |middle, end| {
            let coarse = solve(&self.temperature, start, end)?;
            let first = solve(&self.temperature, start, middle)?;
            let second = solve(&first.solution().temperature, middle, end)?;
            let error = temperature_error(cx, &self.temperature, &coarse.solution().temperature,
                &second.solution().temperature, adaptive)?;
            Ok((first, second, error))
        })?;
        self.publish(cx, selected.end_time_s, selected.second_half.solution())?;
        // Admission reserved two slots. The original publisher counted one;
        // the accepted first half is the other, without publishing a partial pair.
        self.accepted_steps += 1;
        Ok(selected)
    }
}

pub(in crate::bc::history::march) fn admit_proposal(cx: &Cx<'_>, accepted: usize, maximum: usize,
    proposed: f64, config: AdaptiveHistoryConfig) -> Result<f64, ConductionError>
{
    poll(cx, accepted)?;
    config.validate()?;
    if maximum.checked_sub(accepted).is_none_or(|left| left < 2) {
        return Err(invalid("adaptive history requires two remaining accepted-step slots"));
    }
    if !proposed.is_finite() || proposed < config.min_step_s {
        return Err(invalid("adaptive proposal must be finite and at least the minimum step"));
    }
    Ok(proposed.min(config.max_step_s))
}

pub(in crate::bc::history::march) fn temperature_error(cx: &Cx<'_>, old: &[f64], coarse: &[f64],
    fine: &[f64], config: AdaptiveHistoryConfig) -> Result<f64, ConductionError>
{
    if old.len() != coarse.len() || old.len() != fine.len() || old.is_empty() {
        return Err(invalid("adaptive temperature fields must have matching nonempty shapes"));
    }
    let mut error = 0.0_f64;
    for (i, ((&old, &coarse), &fine)) in old.iter().zip(coarse).zip(fine).enumerate() {
        if i % 512 == 0 { poll(cx, i)?; }
        error = error.max(component_error(old, coarse, fine,
            config.absolute_temperature_k, config.relative_temperature)?);
    }
    Ok(error)
}

pub(in crate::bc::history::march) fn component_error(old: f64, coarse: f64, fine: f64,
    absolute: f64, relative: f64) -> Result<f64, ConductionError>
{
    let change = (coarse - old).abs().max((fine - old).abs());
    let allowance = absolute + relative * change;
    let error = (fine - coarse).abs() / allowance;
    if ![old, coarse, fine, change, allowance, error].iter().all(|v| v.is_finite()) || allowance <= 0.0 {
        return Err(invalid("adaptive error scale or field difference is not representable"));
    }
    Ok(error)
}

// Only the time controller is shared. Each physical owner supplies its original
// endpoint solves and error coordinates; no alternative physics is defined here.
pub(in crate::bc::history::march) fn select_step<T, E>(cx: &Cx<'_>, start: f64, mut end: f64,
    config: AdaptiveHistoryConfig, convert: impl Fn(ConductionError) -> E,
    mut trial: impl FnMut(f64, f64) -> Result<(T, T, f64), E>)
    -> Result<AdaptiveHistoryStep<T>, E>
{
    for attempt in 0..config.max_attempts {
        poll(cx, attempt).map_err(&convert)?;
        let middle = f64::midpoint(start, end);
        if !(start < middle && middle < end) {
            return Err(convert(invalid("adaptive interval has no representable interior midpoint")));
        }
        let (first, second, discrepancy) = trial(middle, end)?;
        // For first order, fine error / (coarse-fine) is
        // (h1^2+h2^2)/(2*h1*h2), equal to one for exact halves.
        let (h1, h2) = (middle - start, end - middle);
        let error = discrepancy * (0.5 * (h1 / h2 + h2 / h1));
        if !error.is_finite() || error < 0.0 {
            return Err(convert(invalid("adaptive discrepancy is not finite and nonnegative")));
        }
        let dt = end - start;
        if error <= 1.0 {
            let next = if error <= 0.125 { dt.min(config.max_step_s / 2.0) * 2.0 } else { dt };
            poll(cx, attempt).map_err(&convert)?;
            return Ok(AdaptiveHistoryStep { start_time_s: start, midpoint_time_s: middle,
                end_time_s: end, first_half: first, second_half: second, normalized_error: error,
                rejected_attempts: attempt, endpoint_solves: (attempt + 1) * 3,
                suggested_step_s: next.clamp(config.min_step_s, config.max_step_s) });
        }
        if attempt + 1 == config.max_attempts || dt <= config.min_step_s {
            return Err(convert(invalid(&format!("adaptive local error ratio {error} exceeds one at the attempt limit or step floor"))));
        }
        let wanted = (0.5 * dt).max(config.min_step_s);
        let mut smaller = start + wanted;
        // Keep the REAL represented interval at or above the refusal floor.
        if smaller - start < config.min_step_s { smaller = smaller.next_up(); }
        if !(start < smaller && smaller < end) {
            return Err(convert(invalid("adaptive interval cannot be reduced within its represented step floor")));
        }
        end = smaller;
    }
    Err(convert(invalid("adaptive history attempt budget exhausted")))
}
