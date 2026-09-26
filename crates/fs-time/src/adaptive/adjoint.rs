//! Discrete reverse derivatives of the Dormand--Prince fifth-order update.
//!
//! This differentiates the numerical step with its time and step size HELD
//! FIXED, not a separately integrated continuous adjoint. The model supplies
//! matrix-free transposed derivative actions; no Jacobian is assembled and no
//! state-sized tangent column is allocated per parameter. As usual in numeric
//! AD, rounding decisions themselves are not differentiated.
//!
//! Smooth ODEs only: PI decisions, rejected trials, event localization and reset
//! laws are outside this map. Hybrid continuous sensitivities remain available
//! separately through `crate::hybrid_sensitivity`. Matched, deterministic RHS
//! and derivative callbacks are caller obligations, not verified by a trait.

use super::{A, B5, AdaptiveError, Workspace, stage_time};

/// A fixed model/parameter point and its RHS vector-Jacobian product (VJP).
/// Callbacks must overwrite EVERY output component, including exact zeros.
/// A VJP computes f_x^T * seed and f_p^T * seed at fixed time; it must not
/// include initial-condition or objective derivatives. Long kernels must
/// provide their own cancellation; the driver polls around each callback.
pub trait OdeVjp {
    fn dimension(&self) -> usize;
    fn parameter_count(&self) -> usize;
    fn rhs(&self, time: f64, state: &[f64], out: &mut [f64]);
    fn rhs_vjp(
        &self, time: f64, state: &[f64], seed: &[f64],
        state_bar: &mut [f64], parameter_bar: &mut [f64],
    ) -> Result<(), String>;
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum AdjointError {
    InvalidInput(&'static str),
    WorkspaceLimit { required: usize, limit: usize },
    Integration(AdaptiveError),
    Derivative(String),
    NonFiniteDerivative,
    NonFiniteAccumulation,
    Cancelled,
}

impl From<AdaptiveError> for AdjointError {
    fn from(error: AdaptiveError) -> Self { Self::Integration(error) }
}
impl std::fmt::Display for AdjointError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RK45 discrete adjoint failed: {self:?}")
    }
}
impl std::error::Error for AdjointError {}

/// The fifth-order endpoint and the pullback of its supplied cotangent.
#[derive(Debug, Clone, PartialEq)]
pub struct StepVjp {
    pub value: Vec<f64>,
    pub initial: Vec<f64>,
    pub parameters: Vec<f64>,
}

fn finite(values: &[f64]) -> bool { values.iter().all(|x| x.is_finite()) }

// Scalar scratch ceiling, including the returned vectors while being built:
// Workspace 9n, stage cotangents 6n, initial/VJP state 2n, parameter/VJP 2p.
// Excludes caller inputs, model-owned memory and Vec allocation metadata.
fn workspace_size(n: usize, p: usize, limit: usize) -> Result<usize, AdjointError> {
    let required = n.checked_mul(17).and_then(|v| p.checked_mul(2).and_then(|w| v.checked_add(w)))
        .ok_or(AdjointError::InvalidInput("workspace dimension overflow"))?;
    if required > limit { return Err(AdjointError::WorkspaceLimit { required, limit }); }
    Ok(required)
}

fn check<M: OdeVjp>(model: &M, state: &[f64], limit: usize) -> Result<(), AdjointError> {
    if state.is_empty() || state.len() != model.dimension() || !finite(state) {
        return Err(AdjointError::InvalidInput("state must be finite, nonempty and dimension-matched"));
    }
    workspace_size(state.len(), model.parameter_count(), limit)?;
    Ok(())
}

fn poll<Cancel: FnMut() -> bool>(cancelled: &mut Cancel) -> Result<(), AdjointError> {
    if cancelled() { Err(AdjointError::Cancelled) } else { Ok(()) }
}

/// Evaluate a frozen step using exactly the production trial's stages, their
/// summation order, and the production fifth-order endpoint expression. `end`
/// is retained separately from h: subtracting floating-point time endpoints
/// does not in general recover the original proposed h bit-for-bit.
fn replay<M: OdeVjp, Cancel: FnMut() -> bool>(
    model: &M, state: &[f64], time: f64, end: f64, h: f64, cancelled: &mut Cancel,
) -> Result<Workspace, AdjointError> {
    if !time.is_finite() || !end.is_finite() || end <= time || !h.is_finite() || h <= 0.0 {
        return Err(AdjointError::InvalidInput("frozen step must have finite increasing times and positive h"));
    }
    poll(cancelled)?;
    let mut work = Workspace::new(state.len());
    if !work.stages(state, time, end, h, &|t, u, out| model.rhs(t, u, out), cancelled)? {
        return Err(AdjointError::Cancelled);
    }
    for (i, value) in work.next.iter_mut().enumerate() {
        let mut du = 0.0f64;
        for (j, kj) in work.k.iter().enumerate() { du = B5[j].mul_add(kj[i], du); }
        *value = h.mul_add(du, state[i]);
        if !value.is_finite() {
            return Err(AdaptiveError::NonFiniteState { stage: 7, component: i }.into());
        }
    }
    poll(cancelled)?;
    Ok(work)
}

/// Pull a cotangent through ONE fifth-order step, with all time coordinates
/// frozen. `max_workspace_components` bounds owned scalar scratch (17n + 2p),
/// not allocator metadata or model-owned memory. Zero model parameters are
/// allowed. Invalid, cancelled or non-finite evaluations return no gradient.
/// This is a step map, not an accuracy-controlled advance: no tolerance or
/// acceptance claim is attached to its endpoint.
pub fn step_vjp<M: OdeVjp, Cancel: FnMut() -> bool>(
    model: &M, time: f64, state: &[f64], h: f64, seed: &[f64],
    max_workspace_components: usize, cancelled: &mut Cancel,
) -> Result<StepVjp, AdjointError> {
    check(model, state, max_workspace_components)?;
    if seed.len() != state.len() || !finite(seed) {
        return Err(AdjointError::InvalidInput("cotangent must be finite and dimension-matched"));
    }
    let end = time + h;
    let work = replay(model, state, time, end, h, cancelled)?;
    reverse(model, state, time, end, h, seed, work, cancelled)
}

#[allow(clippy::too_many_arguments)]
fn reverse<M: OdeVjp, Cancel: FnMut() -> bool>(
    model: &M, state: &[f64], time: f64, end: f64, h: f64, seed: &[f64],
    mut work: Workspace, cancelled: &mut Cancel,
) -> Result<StepVjp, AdjointError> {
    let n = state.len();
    let p = model.parameter_count();
    let mut initial = seed.to_vec();
    let mut parameters = vec![0.0; p];
    let mut bars = vec![vec![0.0; n]; 6];
    let mut state_bar = vec![0.0; n];
    let mut parameter_bar = vec![0.0; p];
    // k6 contributes only to error control, which is frozen. k1 has b1=0,
    // but must still be reversed because later stages depend on it.
    for (j, row) in bars.iter_mut().enumerate() {
        for (value, seed) in row.iter_mut().zip(seed) { *value = (h * B5[j]) * seed; }
        if !finite(row) { return Err(AdjointError::NonFiniteAccumulation); }
    }
    for stage in (0..6).rev() {
        poll(cancelled)?;
        work.stage_values(state, h, stage);
        state_bar.fill(f64::NAN);
        parameter_bar.fill(f64::NAN);
        model.rhs_vjp(stage_time(time, end, h, stage), &work.stage, &bars[stage],
            &mut state_bar, &mut parameter_bar).map_err(AdjointError::Derivative)?;
        poll(cancelled)?;
        if !finite(&state_bar) || !finite(&parameter_bar) {
            return Err(AdjointError::NonFiniteDerivative);
        }
        for (value, update) in initial.iter_mut().zip(&state_bar) { *value += update; }
        for (value, update) in parameters.iter_mut().zip(&parameter_bar) { *value += update; }
        if !finite(&initial) || !finite(&parameters) { return Err(AdjointError::NonFiniteAccumulation); }
        for (j, row) in bars.iter_mut().enumerate().take(stage) {
            let a = A[stage - 1][j];
            if a != 0.0 {
                for (value, update) in row.iter_mut().zip(&state_bar) {
                    *value = (h * a).mul_add(*update, *value);
                }
                if !finite(row) { return Err(AdjointError::NonFiniteAccumulation); }
            }
        }
    }
    poll(cancelled)?;
    Ok(StepVjp { value: work.next, initial, parameters })
}

#[cfg(test)]
mod tests;
