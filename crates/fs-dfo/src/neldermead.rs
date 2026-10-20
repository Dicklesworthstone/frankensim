//! Nelder-Mead simplex (the deterministic local derivative-free baseline).
//! Standard reflection/expansion/contraction/shrink coefficients, total-order
//! objective ranking and lowest-index tie breaks. `state` provides hard budgets,
//! fallible callbacks, optional boxes and cancellation-correct in-process resume.

/// Persistent simplex state with explicit stopping and callback-work allowances.
pub mod state;

/// Minimize `f` from `x0` with signed, nonzero initial axis step `h`.
/// Returns the best actually evaluated (point, value, callback count). The hard
/// limit includes initialization: a small allowance returns a partial-simplex
/// best rather than silently exceeding `max_evals`. At least one call is needed.
/// A target of negative infinity disables target stopping. Other targets must
/// be finite. Fully collapsed simplexes stop; no stationarity is inferred.
///
/// This callback-only compatibility surface shares the stateful worker, but has
/// no execution context. Use `state::NelderMeadState::advance` for cancellable,
/// resumable physical work and original typed callback failures.
///
/// # Panics
/// Empty/nonfinite data, zero/absorbed initial step, zero budget, invalid target,
/// nonfinite generated arithmetic/objective, allocation refusal or callback panic.
#[must_use]
pub fn nelder_mead<F: FnMut(&[f64]) -> f64>(
    f: &mut F, x0: &[f64], h: f64, max_evals: usize, f_target: f64,
) -> (Vec<f64>, f64, usize) {
    state::run_legacy(f, x0, h, max_evals, f_target)
}
