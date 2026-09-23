//! Prepared version of the canonical Gonzalez discrete-gradient step.
//!
//! Storage, operators, discrete-gradient formula and energy accounting are the
//! same as the reference path. Only scratch ownership and linear factorization
//! traversal change. No polynomial truncation, fast-math authority or linearized
//! replacement of nonlinear storage is introduced here.
use crate::{PhsError, PortHamiltonian, NEWTON_MAX, NEWTON_TOL, discrete_gradient_into_unchecked};
use fs_la::LuWorkspace;
mod structure;
mod analytic;
mod dissipation;
mod relaxation;
use dissipation::Dissipation;
use structure::FlowPattern;

/// Scalar ledger from a prepared step; state and port output use caller buffers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PreparedStepRecord {
    /// Actual change in Hamiltonian [J].
    pub delta_h: f64,
    /// `dt * (dg^T R dg + dg^T D(midpoint,dg))` [J].
    /// The additional nonlinear port is zero on the original stepping APIs.
    pub dissipated: f64,
    /// `dt * u^T G^T dg` [J].
    pub supplied: f64,
    /// Number of Newton updates performed.
    pub newton_iters: usize,
    /// Infinity norm of the accepted equation residual.
    pub solver_residual: f64,
}

impl PreparedStepRecord {
    /// Equation-balance diagnostic, not an independent structural passivity proof.
    #[must_use]
    pub fn balance_residual(&self) -> f64 { self.delta_h + self.dissipated - self.supplied }

    /// Independent supply-rate defect, subject to the disclosed solver residual.
    #[must_use]
    pub fn supply_defect(&self) -> f64 { self.delta_h - self.supplied }
}

/// Refusal from cooperative prepared stepping. Neither output is published.
#[derive(Debug, Clone, PartialEq)]
pub enum PreparedStepError {
    /// The caller requested cancellation at a solver boundary.
    Cancelled,
    /// The numerical owner refused the candidate.
    Solver(PhsError),
}

impl From<PhsError> for PreparedStepError {
    fn from(error: PhsError) -> Self { Self::Solver(error) }
}

impl core::fmt::Display for PreparedStepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => f.write_str("prepared step cancelled"),
            Self::Solver(error) => core::fmt::Display::fmt(error, f),
        }
    }
}

impl std::error::Error for PreparedStepError {}

/// Reusable implicit-step scratch with transactional, caller-owned outputs.
///
/// Construction allocates. Successful or failed [`Self::step_into`] calls do not
/// allocate in the solver itself. That guarantee is conditional on the supplied
/// [`crate::Storage`] callbacks not allocating. Storage must be read-only during
/// a trial; material history belongs at the caller's accepted-step boundary.
///
/// Newton uses central differences by default or explicit analytic tangents,
/// with the same dense LU. Neither is a hard-real-time
/// performance certificate. The LU traversal is unblocked and may differ from
/// [`crate::step`] by floating-point roundoff. Both paths use the same Gonzalez
/// kernel. Exact nonzero J-R dependency rows are retained across Newton probes;
/// they are refreshed from the current operators before every step. No numerical
/// threshold drops couplings. The Jacobian and LU remain dense. No existing
/// reference-path bit semantics are changed.
#[derive(Debug)]
pub struct StepWorkspace {
    n: usize,
    m: usize,
    max_iterations: usize,
    x: Vec<f64>,
    best: Vec<f64>,
    trial: Vec<f64>,
    midpoint: Vec<f64>,
    effort: Vec<f64>,
    nonlinear_loss: Vec<f64>,
    nonlinear_effort: Vec<f64>,
    residual: Vec<f64>,
    plus: Vec<f64>,
    minus: Vec<f64>,
    forcing: Vec<f64>,
    delta: Vec<f64>,
    jacobian: Vec<f64>,
    output: Vec<f64>,
    lu: LuWorkspace,
    flow: FlowPattern,
}

fn dimensions(what: &'static str) -> PhsError { PhsError::Dimension { what } }

fn norm(values: &[f64]) -> Result<f64, PhsError> {
    let mut result = 0.0_f64;
    for &value in values {
        if !value.is_finite() { return Err(dimensions("nonfinite prepared-step value")); }
        result = result.max(value.abs());
    }
    Ok(result)
}

fn zeroed(len: usize) -> Result<Vec<f64>, PhsError> {
    let mut values = Vec::new();
    values.try_reserve_exact(len).map_err(|_| dimensions("prepared workspace capacity"))?;
    values.resize(len, 0.0);
    Ok(values)
}

#[allow(clippy::too_many_arguments)] // one residual with explicitly borrowed scratch
fn residual_into(sys: &PortHamiltonian, x0: &[f64], x: &[f64], gu: &[f64],
    dt: f64, pattern: &FlowPattern, midpoint: &mut [f64], effort: &mut [f64], out: &mut [f64],
    dissipation: Option<Dissipation<'_>>, nonlinear_loss: &mut [f64]) -> Result<f64, PhsError>
{
    norm(x)?;
    discrete_gradient_into_unchecked(sys.storage.as_ref(), x0, x, midpoint, effort);
    norm(effort)?;
    if let Some(port) = dissipation { port.evaluate(midpoint, effort, nonlinear_loss)?; }
    for (row, value) in out.iter_mut().enumerate() {
        let mut flow = 0.0;
        for &col in pattern.row(row) {
            flow += (sys.j[row * sys.n + col] - sys.r[row * sys.n + col]) * effort[col];
        }
        if dissipation.is_some() { flow -= nonlinear_loss[row]; }
        *value = x[row] - x0[row] - dt * (flow + gu[row]);
    }
    norm(out)
}

impl StepWorkspace {
    /// Prepare for systems of this fixed state and port dimension.
    /// Operators are read from `sys` at each call, never cached across models.
    ///
    /// # Errors
    /// Refuses unrepresentable dimensions, malformed raw structure or capacity failure.
    pub fn new(sys: &PortHamiltonian) -> Result<Self, PhsError> {
        let n = sys.n;
        let m = sys.m;
        let square = n.checked_mul(n).ok_or_else(|| dimensions("workspace state extent"))?;
        let ports = n.checked_mul(m).ok_or_else(|| dimensions("workspace port extent"))?;
        if sys.j.len() != square || sys.r.len() != square || sys.g.len() != ports {
            return Err(dimensions("prepared J/R/G sizes"));
        }
        Ok(Self {
            n, m, max_iterations: NEWTON_MAX,
            x: zeroed(n)?, best: zeroed(n)?, trial: zeroed(n)?,
            midpoint: zeroed(n)?, effort: zeroed(n)?, residual: zeroed(n)?,
            nonlinear_loss: zeroed(n)?, nonlinear_effort: zeroed(n)?,
            plus: zeroed(n)?, minus: zeroed(n)?, forcing: zeroed(n)?,
            delta: zeroed(n)?, jacobian: zeroed(square)?, output: zeroed(m)?,
            lu: LuWorkspace::new(n).map_err(|_| dimensions("prepared LU capacity"))?,
            flow: FlowPattern::new(n)?,
        })
    }

    /// Set a bounded Newton-update budget (at most the reference path's budget).
    /// Zero permits only a state already satisfying the step equation.
    ///
    /// # Errors
    /// Refuses budgets above the canonical reference maximum.
    pub fn set_iteration_limit(&mut self, limit: usize) -> Result<(), PhsError> {
        if limit > NEWTON_MAX { return Err(dimensions("prepared Newton budget")); }
        self.max_iterations = limit;
        Ok(())
    }

    /// Solve one held-input step, publishing `x_next` and `y` only on success.
    ///
    /// `y` is the discrete-gradient port output, NOT the endpoint gradient.
    /// Physical callers must still apply their independent energy and validity
    /// gates before making the candidate state authoritative. This routine has
    /// no simulation clock; `dt` is a physical duration, not a time accumulator.
    ///
    /// # Errors
    /// Rejects mismatched lengths, invalid duration, nonfinite data/intermediates,
    /// singular Newton systems and exhausted Newton convergence budgets.
    /// Both output slices remain unchanged on every error.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_into(&mut self, sys: &PortHamiltonian, x0: &[f64], u: &[f64],
        dt: f64, x_next: &mut [f64], y: &mut [f64]) -> Result<PreparedStepRecord, PhsError>
    {
        // This convenience path never requests cancellation.
        self.step_into_controlled(sys, x0, u, dt, x_next, y, || false)
            .map_err(|error| match error {
                PreparedStepError::Solver(error) => error,
                PreparedStepError::Cancelled => dimensions("unexpected cancellation"),
            })
    }

    /// Step with cooperative cancellation, including inside Newton/Jacobian work.
    ///
    /// `cancelled` must be read-only with respect to the physical model. It is
    /// polled before numerical work, once per Newton iteration and Jacobian
    /// column, around the dense solve, and immediately before publication.
    /// Individual storage callbacks and the LU kernel are not preemptible.
    /// This is bounded cooperative polling, not a hard wall-time guarantee.
    ///
    /// # Errors
    /// Returns [`PreparedStepError::Cancelled`] or the original numerical refusal.
    /// Both caller output buffers remain unchanged on every refusal.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_into_controlled<F: FnMut() -> bool>(
        &mut self, sys: &PortHamiltonian, x0: &[f64], u: &[f64], dt: f64,
        x_next: &mut [f64], y: &mut [f64], cancelled: F,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        self.step_with_hessian(sys, x0, u, dt, x_next, y, None, None, cancelled)
    }

    /// Solve the same Gonzalez equation using the caller's exact storage Hessian
    /// action instead of finite-difference residual probes. No heap allocation is
    /// added. The action must fill `out` with `H''(x) * direction` for THIS system,
    /// including every contact, internal state and frozen material history.
    /// Return false for an unsupported action. Missing/nonfinite actions refuse;
    /// there is no silent fallback. Piecewise laws use their branch tangent.
    /// The dense LU and all original equation/energy acceptance limits remain.
    ///
    /// # Errors
    /// The original step refusals plus an unsupported/nonfinite Hessian action.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_into_analytic(
        &mut self, sys: &PortHamiltonian, x0: &[f64], u: &[f64], dt: f64,
        x_next: &mut [f64], y: &mut [f64],
        hessian: &dyn Fn(&[f64], &[f64], &mut [f64]) -> bool,
    ) -> Result<PreparedStepRecord, PhsError> {
        self.step_into_analytic_controlled(sys, x0, u, dt, x_next, y, hessian, || false)
            .map_err(|error| match error {
                PreparedStepError::Solver(error) => error,
                PreparedStepError::Cancelled => dimensions("unexpected cancellation"),
            })
    }

    /// Analytic Newton with the same cooperative polling and transactional output
    /// contract as `step_into_controlled`. Callbacks must be read-only and must
    /// not allocate for allocation-free execution. Neither this action nor LU is
    /// preemptible; this is not a hard-real-time qualification.
    ///
    /// # Errors
    /// Same refusals as `step_into_analytic`, plus caller cancellation.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_into_analytic_controlled<F: FnMut() -> bool>(
        &mut self, sys: &PortHamiltonian, x0: &[f64], u: &[f64], dt: f64,
        x_next: &mut [f64], y: &mut [f64],
        hessian: &dyn Fn(&[f64], &[f64], &mut [f64]) -> bool, cancelled: F,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        self.step_with_hessian(sys, x0, u, dt, x_next, y, Some(hessian), None, cancelled)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines, clippy::type_complexity)]
    fn step_with_hessian<F: FnMut() -> bool>(
        &mut self, sys: &PortHamiltonian, x0: &[f64], u: &[f64], dt: f64,
        x_next: &mut [f64], y: &mut [f64],
        hessian: Option<&dyn Fn(&[f64], &[f64], &mut [f64]) -> bool>,
        dissipation: Option<Dissipation<'_>>, mut cancelled: F,
    ) -> Result<PreparedStepRecord, PreparedStepError> {
        let mut poll = || {
            if cancelled() { Err(PreparedStepError::Cancelled) } else { Ok(()) }
        };
        poll()?;
        if hessian.is_some() && dissipation.is_some_and(|d| d.tangent.is_none()) {
            return Err(dimensions("analytic nonlinear dissipation requires its exact tangent").into());
        }
        let n = self.n;
        if sys.n != n || sys.m != self.m || x0.len() != n || x_next.len() != n
            || u.len() != self.m || y.len() != self.m
            || sys.j.len() != self.jacobian.len() || sys.r.len() != self.jacobian.len()
            || sys.g.len() != n * self.m
        { return Err(dimensions("prepared state/input/output/structure dimensions").into()); }
        if !dt.is_finite() || dt < 0.0 { return Err(dimensions("dt must be finite and non-negative").into()); }
        let mut scale = norm(x0)?.max(1.0e-30);
        norm(u)?;
        norm(&sys.j)?; norm(&sys.r)?; norm(&sys.g)?;
        // Inspect structure once per step, not once per Newton residual probe.
        // Refresh even on workspace reuse with another same-dimension system.
        self.flow.refresh(sys)?;
        let h0 = sys.hamiltonian(x0);
        if !h0.is_finite() { return Err(dimensions("nonfinite initial Hamiltonian").into()); }
        self.forcing.fill(0.0);
        for row in 0..n {
            for (port, &input) in u.iter().enumerate() {
                self.forcing[row] += sys.g[row * self.m + port] * input;
            }
            let increment = dt * self.forcing[row];
            if !increment.is_finite() { return Err(dimensions("nonfinite input increment").into()); }
            scale = scale.max(increment.abs());
        }
        self.x.copy_from_slice(x0);
        let mut best_norm = f64::INFINITY;
        let mut previous_norm = f64::INFINITY;
        let mut initial_norm = 0.0;
        let mut stagnant = 0;
        let mut iterations = 0;
        loop {
            poll()?;
            let rnorm = residual_into(sys, x0, &self.x, &self.forcing, dt, &self.flow,
                &mut self.midpoint, &mut self.effort, &mut self.residual,
                dissipation, &mut self.nonlinear_loss)?;
            if iterations == 0 { initial_norm = rnorm; }
            let improved = rnorm < best_norm;
            if improved { self.best.copy_from_slice(&self.x); best_norm = rnorm; }
            if rnorm <= NEWTON_TOL * scale { break; }
            // Contact onset may first increase the residual. Falling residuals
            // after that are progress even before they beat the entry iterate.
            // Keep the global best for admission, not for the stagnation clock.
            if rnorm < previous_norm { stagnant = 0; } else { stagnant += 1; }
            previous_norm = rnorm;
            if stagnant >= 3 || iterations >= self.max_iterations {
                // Preserve the reference's disclosed FD-noise-floor acceptance.
                // A zero budget must not silently admit a nonzero initial residual.
                let acceptance_scale = norm(&self.best)?.max(scale).max(initial_norm);
                if iterations > 0 && best_norm <= 1.0e-6 * acceptance_scale {
                    self.x.copy_from_slice(&self.best);
                    break;
                }
                return Err(PhsError::NewtonStalled { residual: best_norm }.into());
            }
            if let Some(action) = hessian {
                self.analytic_jacobian_into(sys, x0, dt, action, dissipation, &mut poll)?;
            } else {
                for col in 0..n {
                    poll()?;
                    let h = 1.0e-6 * (scale + self.x[col].abs());
                    if !h.is_finite() || h <= 0.0 { return Err(dimensions("invalid finite-difference increment").into()); }
                    self.trial.copy_from_slice(&self.x);
                    self.trial[col] += h;
                    residual_into(sys, x0, &self.trial, &self.forcing, dt, &self.flow,
                        &mut self.midpoint, &mut self.effort, &mut self.plus,
                        dissipation, &mut self.nonlinear_loss)?;
                    self.trial[col] = self.x[col] - h;
                    residual_into(sys, x0, &self.trial, &self.forcing, dt, &self.flow,
                        &mut self.midpoint, &mut self.effort, &mut self.minus,
                        dissipation, &mut self.nonlinear_loss)?;
                    for row in 0..n {
                        self.jacobian[row * n + col] = (self.plus[row] - self.minus[row]) / (2.0 * h);
                    }
                }
            }
            poll()?;
            self.lu.solve_into(&self.jacobian, &self.residual, &mut self.delta)
                .map_err(|_| PhsError::NewtonStalled { residual: rnorm })?;
            poll()?;
            for row in 0..n { self.x[row] -= self.delta[row]; }
            iterations += 1;
        }
        // Re-evaluate at the accepted candidate: FD probes and best-iterate
        // restoration must never leave a stale effort in the energy ledger.
        let solver_residual = residual_into(sys, x0, &self.x, &self.forcing, dt, &self.flow,
            &mut self.midpoint, &mut self.effort, &mut self.residual,
                dissipation, &mut self.nonlinear_loss)?;
        let mut dissipated = 0.0;
        for row in 0..n {
            let mut value = 0.0;
            for col in 0..n { value += sys.r[row * n + col] * self.effort[col]; }
            dissipated += self.effort[row] * value;
        }
        if dissipation.is_some() {
            dissipated += dissipation::power(&self.effort, &self.nonlinear_loss)?;
        }
        dissipated *= dt;
        self.output.fill(0.0);
        for port in 0..self.m {
            for row in 0..n { self.output[port] += sys.g[row * self.m + port] * self.effort[row]; }
        }
        let supplied = dt * u.iter().zip(&self.output).map(|(a, b)| a * b).sum::<f64>();
        let delta_h = sys.hamiltonian(&self.x) - h0;
        norm(&self.output)?;
        norm(&[delta_h, dissipated, supplied, solver_residual])?;
        let record = PreparedStepRecord { delta_h, dissipated, supplied, newton_iters: iterations, solver_residual };
        poll()?;
        x_next.copy_from_slice(&self.x);
        y.copy_from_slice(&self.output);
        Ok(record)
    }
}
