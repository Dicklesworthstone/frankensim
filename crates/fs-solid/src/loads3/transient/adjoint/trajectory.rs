//! Bounded recording and reverse traversal of the ORIGINAL inertial trajectory.
//! For fixed linear operators, step Jacobians do not depend on primal motion.
//! Retain only accepted durations/jump kinds, not an O(DOFs*steps) state tape.
//! The optional full force-gradient output is separately bounded before forward
//! advancement. This is a fixed-grid derivative: changing knot times, adaptive
//! decisions, material operators, supports or topology requires other pullbacks.
use super::*;
use crate::loads3::transient::events::{LoadJump3, TransientEventStep3};

#[derive(Debug, Clone, Copy)]
pub struct ControlRecordingBudget3 {
    /// Smooth intervals and zero-duration jumps each count as one action.
    pub maximum_actions: usize,
    /// Maximum force-gradient scalars: (actions+1)*free_dofs, including f0.
    /// Other vectors are bounded by the adjoint's independent DOF allowance.
    pub maximum_force_gradient_values: usize,
}
impl Default for ControlRecordingBudget3 {
    fn default() -> Self { Self { maximum_actions: 10_000, maximum_force_gradient_values: 1_000_000 } }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ControlForceKind3 { Initial, IntervalEndpoint, RightLimit }

#[derive(Debug, Clone, PartialEq)]
pub struct ForceCotangent3 {
    pub time_s: f64,
    pub kind: ControlForceKind3,
    /// Objective/N, in the recording's free-DOF order.
    pub values: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TransientControlGradient3 {
    /// Derivatives with respect to independent initial q/v. Initial acceleration
    /// was solved, so its returned independent cotangent is exactly zero.
    pub initial: MotionCotangent3,
    /// Initial load then one endpoint per action, in original forward order.
    /// Equal times intentionally distinguish left and right event loads.
    pub forces: Vec<ForceCotangent3>,
    /// Derivative of the ONE frozen equivalent thermal vector, summed across
    /// initialization, every interval and every acceleration reset.
    pub thermal_equivalent_force: Vec<f64>,
    pub total_adjoint_iterations: usize,
    pub maximum_adjoint_relative_residual: f64,
}

#[derive(Debug, Clone, Copy)]
struct Action3 { time_s: f64, step_size_s: Option<f64> }

/// Owning wrapper prevents unrecorded mutation of the accepted primal state.
/// Clone forks both the original run and its record; assembled operators are
/// shared. Read-only access still permits dynamic stress/reaction observations.
#[derive(Clone)]
pub struct RecordedTransient3<'g, 'm, 'p, 's> {
    run: TetTransientRun3<'g, 'm, 'p, 's>,
    adjoint: TransientAdjoint3<'g, 'm>,
    actions: Vec<Action3>,
    initial_time_s: f64,
    budget: ControlRecordingBudget3,
}

impl<'g, 'm, 'p, 's> TetTransientRun3<'g, 'm, 'p, 's> {
    /// Begin control differentiation before any finite-time advancement. The
    /// initial acceleration is the original mass solve, not an independent
    /// user parameter. No primal matrices, loads or numerical policies change.
    pub fn record_controls(self, budget: ControlRecordingBudget3,
        adjoint_config: TransientAdjointConfig3, cx: &Cx<'_>)
        -> Result<RecordedTransient3<'g, 'm, 'p, 's>, TransientAdjointError3> {
        check(cx)?;
        if self.state.step != 0 {
            return Err(TransientAdjointError3::Invalid("start control recording before the first finite step"));
        }
        if self.free_dofs().len() > budget.maximum_force_gradient_values {
            return Err(TransientAdjointError3::Budget);
        }
        let adjoint = self.adjoint(adjoint_config, cx)?;
        let initial_time_s = self.state.time_s;
        check(cx)?;
        Ok(RecordedTransient3 { run: self, adjoint, actions: Vec::new(), initial_time_s, budget })
    }
}

impl<'g, 'm, 'p, 's> RecordedTransient3<'g, 'm, 'p, 's> {
    pub fn run(&self) -> &TetTransientRun3<'g, 'm, 'p, 's> { &self.run }
    pub fn adjoint(&self) -> &TransientAdjoint3<'g, 'm> { &self.adjoint }
    pub fn force_sample_count(&self) -> usize { self.actions.len()+1 }

    fn reserve(&mut self, extra: usize, cx: &Cx<'_>) -> Result<(), TransientAdjointError3> {
        check(cx)?;
        let actions = self.actions.len().checked_add(extra).ok_or(TransientAdjointError3::Budget)?;
        let samples = actions.checked_add(1).ok_or(TransientAdjointError3::Budget)?;
        let scalars = samples.checked_mul(self.run.free_dofs().len()).ok_or(TransientAdjointError3::Budget)?;
        if actions > self.budget.maximum_actions || scalars > self.budget.maximum_force_gradient_values {
            return Err(TransientAdjointError3::Budget);
        }
        self.actions.try_reserve(extra).map_err(|_| TransientAdjointError3::Budget)?;
        Ok(())
    }

    /// Execute and record one original Newmark interval. Only a fully accepted
    /// candidate changes either the primal state or the action sequence.
    pub fn advance(&mut self, load: &TetLoadCase3<'_, '_>, h_s: f64, cx: &Cx<'_>)
        -> Result<TransientStep3, TransientAdjointError3> {
        self.reserve(1, cx)?;
        let mut candidate = self.run.clone();
        let report = candidate.advance(load, h_s, cx)?;
        check(cx)?;
        self.actions.push(Action3 { time_s: report.time_s, step_size_s: Some(h_s) });
        self.run = candidate;
        Ok(report)
    }

    /// Record an explicit finite force change. q/v are continuous; the reverse
    /// map correctly discards the replaced acceleration's dependence on old load.
    pub fn apply_load_jump(&mut self, load: &TetLoadCase3<'_, '_>, cx: &Cx<'_>)
        -> Result<LoadJump3, TransientAdjointError3> {
        self.reserve(1, cx)?;
        let mut candidate = self.run.clone();
        let report = candidate.apply_load_jump(load, cx)?;
        check(cx)?;
        self.actions.push(Action3 { time_s: report.time_s, step_size_s: None });
        self.run = candidate;
        Ok(report)
    }

    /// One atomic interval plus its right-limit reset. The force-gradient output
    /// has two separate rows at the event time, not a merged or smeared control.
    pub fn advance_to_load_jump(&mut self, before: &TetLoadCase3<'_, '_>,
        after: &TetLoadCase3<'_, '_>, h_s: f64, cx: &Cx<'_>)
        -> Result<TransientEventStep3, TransientAdjointError3> {
        self.reserve(2, cx)?;
        let mut candidate = self.run.clone();
        let report = candidate.advance_to_load_jump(before, after, h_s, cx)?;
        check(cx)?;
        self.actions.push(Action3 { time_s: report.interval.time_s, step_size_s: Some(h_s) });
        self.actions.push(Action3 { time_s: report.jump.time_s, step_size_s: None });
        self.run = candidate;
        Ok(report)
    }

    /// One reverse sweep for ALL endpoint force components and independent
    /// initial q/v. Seed the terminal reduced motion with the derivative of a
    /// scalar objective; direct objective dependence on controls must be added
    /// by that objective's owner. No primal replay or finite-difference solve is
    /// used. Each action performs one bounded transposed solve; failure or
    /// cancellation returns no partial gradient and leaves the record reusable.
    pub fn reverse(&self, terminal: &MotionCotangent3, cx: &Cx<'_>)
        -> Result<TransientControlGradient3, TransientAdjointError3> {
        let n = self.run.free_dofs().len();
        terminal.validate(n, cx)?;
        let mut seed = terminal.clone();
        let mut forces = Vec::with_capacity(self.force_sample_count());
        let mut thermal = vec![0.0; n];
        let mut iterations = 0usize;
        let mut maximum_relative = 0.0_f64;
        for action in self.actions.iter().rev() {
            check(cx)?;
            let pb = match action.step_size_s {
                Some(h) => self.adjoint.pullback_step(h, &seed, cx)?,
                None => self.adjoint.pullback_acceleration(&seed, cx)?,
            };
            for (i, &value) in pb.force.iter().enumerate() {
                if i % 256 == 0 { check(cx)?; }
                thermal[i] = number(thermal[i]+value)?;
            }
            iterations = iterations.checked_add(pb.iterations).ok_or(TransientAdjointError3::Budget)?;
            maximum_relative = maximum_relative.max(pb.true_relative_residual);
            seed = pb.previous;
            forces.push(ForceCotangent3 { time_s: action.time_s,
                kind: if action.step_size_s.is_some() { ControlForceKind3::IntervalEndpoint }
                    else { ControlForceKind3::RightLimit }, values: pb.force });
        }
        let initial = self.adjoint.pullback_acceleration(&seed, cx)?;
        for (i, &value) in initial.force.iter().enumerate() {
            if i % 256 == 0 { check(cx)?; }
            thermal[i] = number(thermal[i]+value)?;
        }
        iterations = iterations.checked_add(initial.iterations).ok_or(TransientAdjointError3::Budget)?;
        maximum_relative = maximum_relative.max(initial.true_relative_residual);
        forces.push(ForceCotangent3 { time_s: self.initial_time_s,
            kind: ControlForceKind3::Initial, values: initial.force });
        forces.reverse();
        check(cx)?;
        Ok(TransientControlGradient3 { initial: initial.previous, forces,
            thermal_equivalent_force: thermal, total_adjoint_iterations: iterations,
            maximum_adjoint_relative_residual: maximum_relative })
    }
}
