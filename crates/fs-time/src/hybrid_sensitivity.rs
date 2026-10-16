//! Forward sensitivities of smooth hybrid ODE trajectories, including event
//! time and reset derivatives. This adapter runs the production hybrid driver;
//! it does not finite-difference whole trajectories or differentiate RK45's
//! adaptive controller. Its tangent equations approximate continuous-model
//! derivatives, not a discrete adjoint of the numerical event locator.
//!
//! For a transverse guard g(t, x, p) = 0, dt/dp = -(g_x S + g_p)/(g_t + g_x f-).
//! At a continuing reset, S+ = R_x(S- + f- dt/dp) + R_p + R_t dt/dp - f+ dt/dp.
//! Terminal resets omit the last term: their result is the derivative of the
//! state AT the parameter-dependent terminal time, not at a fixed clock time.
//! This is the saltation/event-time correction (Kong et al., arXiv:2306.06862).
//!
//! Requires correct supplied derivatives and a locally unchanged, isolated,
//! transverse event sequence. Initial/same-time events and competing brackets
//! are refused, not assigned a fictitious smooth gradient. Finite guard scans
//! do not prove isolation; numerical gradients are not physical certificates.

use crate::AdaptiveState;
use crate::adaptive::events::multiple::{EventSetHit, EventSpec};
use crate::adaptive::events::multiple::run::{
    HybridConfig, HybridError, HybridReport, HybridReset, HybridState, HybridSystem,
    ResetAction, run_hybrid,
};

/// Derivative actions for model-defined parameter directions. A direction can
/// seed an initial condition, a physical parameter, or a linear combination.
/// Callbacks must overwrite every output component, including exact zeros.
/// Parameters, their direction seeds, and the model must remain fixed on resume.
/// Long callback kernels must handle their own cancellation; the hybrid driver
/// polls around complete RHS/reset calls, not inside these derivative actions.
pub trait HybridDerivatives: HybridSystem {
    /// f_x * tangent + f_p * seed[direction], at fixed time.
    #[allow(clippy::too_many_arguments)]
    fn rhs_tangent(
        &self, mode: &Self::Mode, time: f64, state: &[f64], direction: usize,
        tangent: &[f64], out: &mut [f64],
    );

    /// Write g_x and return g_t, holding state and parameters fixed.
    fn guard_gradient(
        &self, mode: &Self::Mode, id: u64, time: f64, state: &[f64], gradient: &mut [f64],
    ) -> Result<f64, String>;

    /// Direct g_p * seed[direction], holding time and state fixed.
    fn guard_parameter(
        &self, mode: &Self::Mode, id: u64, time: f64, state: &[f64], direction: usize,
    ) -> Result<f64, String>;

    /// R_x * tangent + R_p * seed[direction] + R_t * time_tangent.
    /// `tangent` ALREADY includes the pre-event f- * dt/dp shift. Differentiate
    /// the actual reset, including any projection, not an unrelated ideal map.
    #[allow(clippy::too_many_arguments)]
    fn reset_tangent(
        &self, mode: &Self::Mode, event: &EventSetHit, state: &[f64],
        reset: &HybridReset<Self::Mode>, direction: usize, tangent: &[f64],
        time_tangent: f64, out: &mut [f64],
    ) -> Result<(), String>;

    /// Explicit nonnegative finite floor on |g_t + g_x f-|, in this guard's
    /// units per time. Guard rescaling requires rescaling this floor as well.
    fn crossing_speed_floor(&self, mode: &Self::Mode, id: u64) -> f64;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Layout { n: usize, directions: usize, tangent_end: usize, total: usize }

/// Whether returned state derivatives hold clock time fixed or follow the
/// parameter-dependent terminal event. Pending resets still use FixedTime.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TangentConvention { FixedTime, TerminalEvent }

/// Clone checkpoints retain the augmented ODE state AND pending-reset/guard
/// activation state. The underlying hybrid state is exposed read-only.
#[derive(Debug, Clone)]
pub struct SensitivityState<Mode> {
    trajectory: HybridState<Mode>,
    layout: Layout,
}

impl<Mode> SensitivityState<Mode> {
    pub fn trajectory(&self) -> &HybridState<Mode> { &self.trajectory }
    pub fn values(&self) -> &[f64] { &self.trajectory.integration().u[..self.layout.n] }
    pub fn tangent(&self, direction: usize) -> Option<&[f64]> {
        if direction >= self.layout.directions { return None; }
        let start = self.layout.n * (direction + 1);
        Some(&self.trajectory.integration().u[start..start + self.layout.n])
    }
    /// Derivatives of the MOST RECENTLY COMMITTED reset's event time. During
    /// a pending reset these still refer to the preceding committed event.
    pub fn last_event_time_tangents(&self) -> Option<&[f64]> {
        (self.trajectory.transitions() > 0)
            .then(|| &self.trajectory.integration().u[self.layout.tangent_end..self.layout.total - 1])
    }
    pub fn convention(&self) -> TangentConvention {
        if self.trajectory.is_terminated() { TangentConvention::TerminalEvent }
        else { TangentConvention::FixedTime }
    }
}

/// Augment a model with forward tangent columns without changing its guards,
/// reset laws, or the hybrid driver's checkpoint and cancellation machinery.
/// All augmented components participate in RK45 error control: choose units
/// and parameter seed scales compatible with the scalar absolute tolerance.
/// The adaptive mesh need not match a primal-only solve bit-for-bit.
pub struct ForwardSensitivity<'a, S> { system: &'a S, layout: Layout }

impl<'a, S: HybridDerivatives> ForwardSensitivity<'a, S> {
    /// `max_components` caps augmented STATE length, not total RK workspace
    /// memory. The layout is [x, tangent columns, last event-time derivatives, last reset time].
    pub fn new(system: &'a S, n: usize, directions: usize, max_components: usize)
        -> Result<Self, HybridError>
    {
        let tangent_end = directions.checked_add(1).and_then(|m| n.checked_mul(m));
        let total = tangent_end.and_then(|v| v.checked_add(directions)).and_then(|v| v.checked_add(1));
        let Some((tangent_end, total)) = tangent_end.zip(total) else {
            return Err(HybridError::Model("sensitivity dimension overflow".into()));
        };
        if n == 0 || directions == 0 || total > max_components {
            return Err(HybridError::Model("nonempty sensitivities must fit the component budget".into()));
        }
        Ok(Self { system, layout: Layout { n, directions, tangent_end, total } })
    }

    /// Initial fixed-time tangent columns, each of length n, must be supplied
    /// explicitly. Initial time is held fixed in all parameter directions.
    pub fn initial(
        &self, mode: S::Mode, time: f64, values: &[f64], step: f64, tangents: &[f64],
    ) -> Result<SensitivityState<S::Mode>, HybridError> {
        if values.len() != self.layout.n || tangents.len() != self.layout.tangent_end - self.layout.n
            || !time.is_finite() || !step.is_finite() || step <= 0.0
            || values.iter().chain(tangents).any(|v| !v.is_finite())
        {
            return Err(HybridError::Model("invalid initial sensitivity state".into()));
        }
        self.system.validate_state(&mode, values).map_err(HybridError::Model)?;
        let mut augmented = Vec::with_capacity(self.layout.total);
        augmented.extend_from_slice(values);
        augmented.extend_from_slice(tangents);
        augmented.resize(self.layout.total, 0.0);
        augmented[self.layout.total - 1] = time;
        Ok(SensitivityState {
            trajectory: HybridState::new(mode, AdaptiveState::new(time, &augmented, step)),
            layout: self.layout,
        })
    }

    pub fn run<Cancel: FnMut() -> bool>(
        &self, state: &mut SensitivityState<S::Mode>, config: &HybridConfig, cancelled: &mut Cancel,
    ) -> Result<HybridReport, HybridError> {
        if state.layout != self.layout {
            return Err(HybridError::Model("sensitivity checkpoint layout mismatch".into()));
        }
        run_hybrid(&mut state.trajectory, self, config, cancelled)
    }
}

fn finite(values: &[f64], context: &str) -> Result<(), String> {
    if values.iter().all(|v| v.is_finite()) { Ok(()) }
    else { Err(format!("non-finite or unwritten {context}")) }
}

impl<S: HybridDerivatives> HybridSystem for ForwardSensitivity<'_, S> {
    type Mode = S::Mode;

    fn rhs(&self, mode: &S::Mode, time: f64, state: &[f64], out: &mut [f64]) {
        out.fill(f64::NAN);
        if state.len() != self.layout.total || out.len() != self.layout.total { return; }
        let n = self.layout.n;
        self.system.rhs(mode, time, &state[..n], &mut out[..n]);
        for direction in 0..self.layout.directions {
            let start = n * (direction + 1);
            self.system.rhs_tangent(mode, time, &state[..n], direction,
                &state[start..start + n], &mut out[start..start + n]);
        }
        out[self.layout.tangent_end..].fill(0.0);
    }
    fn events(&self, mode: &S::Mode) -> &[EventSpec] { self.system.events(mode) }
    fn guard(&self, mode: &S::Mode, id: u64, time: f64, state: &[f64]) -> f64 {
        if state.len() != self.layout.total { return f64::NAN; }
        self.system.guard(mode, id, time, &state[..self.layout.n])
    }
    fn validate_state(&self, mode: &S::Mode, state: &[f64]) -> Result<(), String> {
        if state.len() != self.layout.total { return Err("sensitivity state dimension mismatch".into()); }
        finite(state, "sensitivity state")?;
        self.system.validate_state(mode, &state[..self.layout.n])
    }
    fn reset(&self, mode: &S::Mode, event: &EventSetHit, state: &[f64])
        -> Result<HybridReset<S::Mode>, String>
    {
        self.validate_state(mode, state)?;
        if event.contenders.len() != 1 || event.contenders[0] != event.selected || event.accepted == 0 {
            return Err("sensitivity requires one isolated event reached by a positive-time step".into());
        }
        let n = self.layout.n;
        let values = &state[..n];
        let time = event.selected.occurrence.time;
        let id = event.selected.id;
        if !time.is_finite() || time <= state[self.layout.total - 1] {
            return Err("sensitivity requires positive time between resets".into());
        }
        let mut normal = vec![f64::NAN; n];
        let gt = self.system.guard_gradient(mode, id, time, values, &mut normal)?;
        let floor = self.system.crossing_speed_floor(mode, id);
        finite(&normal, "guard gradient")?;
        if !gt.is_finite() || !floor.is_finite() || floor < 0.0 {
            return Err("invalid guard time derivative or crossing-speed floor".into());
        }
        let mut before = vec![f64::NAN; n];
        self.system.rhs(mode, time, values, &mut before);
        finite(&before, "pre-event dynamics")?;
        let denominator = normal.iter().zip(&before).fold(gt, |sum, (g, f)| g.mul_add(*f, sum));
        if !denominator.is_finite() || denominator.abs() <= floor {
            return Err("grazing or unresolved guard crossing speed".into());
        }
        let reset = self.system.reset(mode, event, values)?;
        if reset.state.len() != n || reset.consumed_ids.as_slice() != [id] {
            return Err("sensitivity reset must retain dimension and consume its single guard".into());
        }
        finite(&reset.state, "reset state")?;
        self.system.validate_state(&reset.mode, &reset.state)?;
        let continuing = reset.action == ResetAction::Continue;
        let mut after = vec![0.0; n];
        if continuing {
            after.fill(f64::NAN);
            self.system.rhs(&reset.mode, time, &reset.state, &mut after);
            finite(&after, "post-event dynamics")?;
        }
        let mut next = vec![f64::NAN; self.layout.total];
        next[..n].copy_from_slice(&reset.state);
        let mut at_event = vec![0.0; n];
        for direction in 0..self.layout.directions {
            let start = n * (direction + 1);
            let tangent = &state[start..start + n];
            let gp = self.system.guard_parameter(mode, id, time, values, direction)?;
            let numerator = normal.iter().zip(tangent).fold(gp, |sum, (g, s)| g.mul_add(*s, sum));
            let dt = -numerator / denominator;
            if !gp.is_finite() || !dt.is_finite() { return Err("non-finite event-time derivative".into()); }
            for i in 0..n { at_event[i] = before[i].mul_add(dt, tangent[i]); }
            finite(&at_event, "event-state tangent")?;
            self.system.reset_tangent(mode, event, values, &reset, direction,
                &at_event, dt, &mut next[start..start + n])?;
            finite(&next[start..start + n], "reset tangent")?;
            if continuing {
                for i in 0..n { next[start + i] = (-after[i]).mul_add(dt, next[start + i]); }
            }
            next[self.layout.tangent_end + direction] = dt;
        }
        next[self.layout.total - 1] = time;
        finite(&next, "saltation update")?;
        Ok(HybridReset { mode: reset.mode, state: next, action: reset.action, consumed_ids: reset.consumed_ids })
    }
}

#[cfg(test)]
mod tests;
